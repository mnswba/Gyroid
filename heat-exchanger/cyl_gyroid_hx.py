"""
cyl_gyroid_hx.py
------------------------------------------------------------------------
Cylindrical two-fluid GYROID heat exchanger with a deep groin-vault arch
and four tapered pillars ending in 3/8"-18 NPT MALE threaded tips.

Builds the geometry as implicit fields (value < 0 = inside) and streams
binary STL files you import into nTop:
    Import Mesh  ->  Implicit Body from Mesh   (or run run_ntop_automate.py)

Install once:   pip install numpy scikit-image
                (optional: pip install trimesh matplotlib  for check/preview)
Run (VS Code):  python cyl_gyroid_hx.py            full export
                python cyl_gyroid_hx.py --preview  quick PNG sections only
                python cyl_gyroid_hx.py --voxel 0.4 --check   coarse test

Outputs (./hx_output):
    HX_solid.stl     printable part (walls, shell, pillars, NPT threads)
    Hot_fluid.stl    hot fluid domain  (for CFD)
    Cold_fluid.stl   cold fluid domain (for CFD)

Layout (mm, origin on the axis at the body datum z = 0):
    Cylinder: radius R, from z = -SKIRT to z = HB, top edge rounded.
    Groin-vault arch (two crossing elliptical arches) cut up to z = Z_APEX
    leaves four pillars at 45/135/225/315 deg. Pillars taper to R_NECK at
    z = -SKIRT, then a 3/8-18 NPT male thread runs to the tip at z = -LEG_H.
    HOT : 45 and 225 deg pillars   COLD: 135 and 315 deg pillars.
    Each pillar holds a header that follows the arch; the four headers are
    quadrants separated by the arch crown and divider walls. Gyroid core
    from z = ZP to the top. A seal slab at the core base lets each fluid's
    gyroid labyrinth open only into its own two quadrants.
------------------------------------------------------------------------
"""

import argparse
import os
import time
import numpy as np
from skimage.measure import marching_cubes

# ======================================================================
# 1. PARAMETERS (mm)
# ======================================================================
VOXEL = 0.2           # 0.4 quick test, 0.2 default, 0.15 crisp threads
CHUNK = 40            # x-slices per chunk (lower it if you run out of RAM)
EXPORT_FLUIDS = True
OUT_DIR = "hx_output"

# Body
R = 70.0              # cylinder radius
HB = 90.0             # top of body above datum z = 0
R_TOP = 10.0          # top edge rounding
SHELL = 2.0           # outer wall thickness

# Arch and pillars
LEG_H = 40.0          # tip depth below the datum
SKIRT = 28.0          # pillar shoulder depth (thread below it)
R_NECK = 11.0         # pillar radius at the shoulder
LEG_POS = 48.0        # radius of the pillar axes
Z_APEX = 20.0         # arch apex height, cut up into the body
ARCH_A = 23.5         # arch half-width at the ellipse base
Z_ARCH0 = -40.0       # ellipse base (below the tips: no vertical walls)
K_ARCH = 6.0          # fillet where arch meets wall / pillar
TAPER_K = 55.0        # pillar flare: r = R_NECK + TAPER_K * u^2
BORE_R = 5.0          # 10 mm flow bore through the thread
ZP_LOW = -20.0        # headers reach down into the pillars to here

# Headers and core
ZP = 30.0             # gyroid core starts here (above the arch apex)
DV = 2.0              # divider wall between quadrant headers
S_SEAL = 1.5          # seal slab at the core base
CELL = (8.0, 8.0, 8.0)  # gyroid cell size x, y, z
T_WALL = 0.6          # gyroid wall thickness (check your printer minimum)

# 3/8"-18 NPT external thread (ASME B1.20.1 basic values)
IN = 25.4
NPT_P = IN / 18.0             # 1.411 pitch
NPT_E0 = 0.61201 * IN         # 15.545 pitch dia at the small end
NPT_H = 0.8 * NPT_P           # truncated thread height
NPT_TAPER = 1.0 / 16.0        # on diameter
# threaded length = LEG_H - SKIRT = 12 mm (standard L2 = 10.2 mm)

# Pillars: (angle deg, fluid)
LEGS = [(45, "hot"), (135, "cold"), (225, "hot"), (315, "cold")]

# ======================================================================
# 2. FIELD HELPERS
# ======================================================================
def smin(a, b, k):
    """Smooth union (blend radius ~k)."""
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0.0, 1.0)
    return b * (1 - h) + a * h - k * h * (1 - h)


def smax(a, b, k):
    """Smooth intersection."""
    return -smin(-a, -b, k)


def leg_xy(a):
    t = np.radians(a)
    return LEG_POS * np.cos(t), LEG_POS * np.sin(t)


def cylinder(x, y, z):
    """Vertical cylinder z = -SKIRT .. HB with rounded edges."""
    r = np.hypot(x, y)
    hh = (HB + SKIRT) / 2
    qx = r - (R - R_TOP)
    qz = np.abs(z - (HB - SKIRT) / 2) - (hh - R_TOP)
    return (np.hypot(np.maximum(qx, 0), np.maximum(qz, 0))
            + np.minimum(np.maximum(qx, qz), 0) - R_TOP)


def arches(x, y, z):
    """Groin vault: two crossing elliptical tunnels (< 0 inside the cut)."""
    hz = Z_APEX - Z_ARCH0
    e = lambda w: np.sqrt((w / ARCH_A) ** 2 + ((z - Z_ARCH0) / hz) ** 2) - 1
    return np.minimum(e(y), e(x)) * min(ARCH_A, hz)


def taper(x, y, z):
    """Envelope that narrows each pillar to R_NECK at the shoulder."""
    u = np.clip((z + SKIRT) / (SKIRT + 10.0), 0, None)
    r_env = R_NECK + TAPER_K * u ** 2
    slope = 2 * TAPER_K * u / (SKIRT + 10.0)
    d = [(np.hypot(x - lx, y - ly) - r_env) / np.sqrt(1 + slope ** 2)
         for lx, ly in (leg_xy(a) for a, _ in LEGS)]
    return np.minimum.reduce(d)


def npt_tip(x, y, z, a):
    """3/8-18 NPT male thread, right-handed about the outward (-z) axis."""
    lx, ly = leg_xy(a)
    dx, dy = x - lx, y - ly
    rho, t = np.hypot(dx, dy), -z
    s_tip = LEG_H - t                               # distance from the tip
    th = np.arctan2(dx, dy)
    rp = (NPT_E0 + np.maximum(s_tip, 0) * NPT_TAPER) / 2
    w = np.mod(t - NPT_P * th / (2 * np.pi), NPT_P) - NPT_P / 2
    r_thr = rp + np.tan(np.radians(60)) * (NPT_P / 4 - np.abs(w))
    r_thr = np.clip(r_thr, rp - NPT_H / 2, rp + NPT_H / 2)
    f = 0.5 * (rho - r_thr)                         # ~ distance to flank
    f = np.maximum(f, -s_tip)                       # cut at the tip
    f = np.maximum(f, (rho - rp) - (s_tip - 0.8))   # tip chamfer
    return np.maximum(f, s_tip - (LEG_H - SKIRT + 1.0))  # into the shoulder


def sector(x, y, a):
    """Quarter sector centred on angle a, shrunk by DV/2 (< 0 inside)."""
    t = np.radians(a)
    u = x * np.cos(t) + y * np.sin(t)
    v = -x * np.sin(t) + y * np.cos(t)
    return (np.abs(v) - u) / np.sqrt(2) + DV / 2


def gyroid(x, y, z):
    """Gyroid divided by its gradient -> approx. distance in mm."""
    kx, ky, kz = (2 * np.pi / c for c in CELL)
    sx, cx = np.sin(kx * x), np.cos(kx * x)
    sy, cy = np.sin(ky * y), np.cos(ky * y)
    sz, cz = np.sin(kz * z), np.cos(kz * z)
    g = sx * cy + sy * cz + sz * cx
    gx = kx * (cx * cy - sz * sx)
    gy = ky * (-sx * sy + cy * cz)
    gz = kz * (-sy * sz + cz * cx)
    grad = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)
    return g / np.maximum(grad, 0.25 * 2 * np.pi / np.mean(CELL))


# ======================================================================
# 3. FIELDS
# ======================================================================
def fields(x, y, z):
    outer = smax(cylinder(x, y, z), -arches(x, y, z), K_ARCH)
    outer = smax(outer, taper(x, y, z), K_ARCH)
    for a, _ in LEGS:
        outer = np.minimum(outer, npt_tip(x, y, z, a))
    inner = outer + SHELL

    G = gyroid(x, y, z)
    core = np.maximum(ZP - z, inner)
    slab_free = (ZP + S_SEAL) - z
    sec = {f: np.minimum.reduce([sector(x, y, a) for a, ff in LEGS if ff == f])
           for f in ("hot", "cold")}
    hot = np.maximum(np.maximum(T_WALL / 2 - G, core),
                     np.minimum(slab_free, sec["hot"]))
    cold = np.maximum(np.maximum(G + T_WALL / 2, core),
                      np.minimum(slab_free, sec["cold"]))

    for a, f in LEGS:
        lx, ly = leg_xy(a)
        plen = np.maximum(np.maximum(ZP_LOW - z, z - (ZP + 0.5)),
                          np.maximum(inner, sector(x, y, a)))
        bore = np.maximum(np.hypot(x - lx, y - ly) - BORE_R, z - (ZP_LOW + 1.0))
        v = np.minimum(plen, bore)
        if f == "hot":
            hot = np.minimum(hot, v)
        else:
            cold = np.minimum(cold, v)

    hot, cold = np.maximum(hot, outer), np.maximum(cold, outer)
    solid = np.maximum(outer, -np.minimum(hot, cold))
    return {"HX_solid": solid, "Hot_fluid": hot, "Cold_fluid": cold}


# ======================================================================
# 4. MESH CHUNK BY CHUNK -> BINARY STL (low memory)
# ======================================================================
class StlWriter:
    def __init__(self, path):
        self.path = path
        self.f = open(path, "wb")
        self.f.write(b"cylindrical gyroid HX".ljust(80, b" "))
        self.f.write(np.uint32(0).tobytes())        # patched in close()
        self.n = 0

    def add(self, verts, faces):
        tri = verts[faces].astype(np.float32)
        nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-12)
        rec = np.zeros(len(tri), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)),
                                         ("a", "<u2")])
        rec["n"], rec["v"] = nrm, tri
        self.f.write(rec.tobytes())
        self.n += len(tri)

    def close(self):
        self.f.seek(80)
        self.f.write(np.uint32(self.n).tobytes())
        self.f.close()


def grid(voxel):
    pad = 1.0 + 0.137 * voxel        # odd offset keeps nodes off flat faces
    xs = np.arange(-R - pad, R + pad + voxel, voxel)
    zs = np.arange(-LEG_H - pad, HB + pad + voxel, voxel)
    return xs, xs.copy(), zs


def export(voxel=VOXEL, out_dir=OUT_DIR, fluids=EXPORT_FLUIDS):
    os.makedirs(out_dir, exist_ok=True)
    xs, ys, zs = grid(voxel)
    print(f"Grid {len(xs)} x {len(ys)} x {len(zs)} "
          f"= {len(xs) * len(ys) * len(zs) / 1e6:.0f} M voxels")
    Y, Z = np.meshgrid(ys, zs, indexing="ij")
    names = ["HX_solid"] + (["Hot_fluid", "Cold_fluid"] if fluids else [])
    writers = {n: StlWriter(os.path.join(out_dir, f"{n}.stl")) for n in names}
    t0 = time.time()
    for i0 in range(0, len(xs) - 1, CHUNK):
        i1 = min(i0 + CHUNK, len(xs) - 1)          # chunks share one slice
        X = xs[i0:i1 + 1][:, None, None]
        f = fields(X, Y[None], Z[None])
        for n in names:
            vol = f[n].astype(np.float32)
            # nodes lying (almost) exactly on a surface make zero-area
            # slivers: nudge them just outside
            eps = 1e-3 * voxel
            vol[np.abs(vol) < eps] = eps
            if vol.min() > 0 or vol.max() < 0:
                continue
            v, faces, _, _ = marching_cubes(vol, level=0.0)
            v = v.astype(np.float64)
            v[:, 0] += i0                          # identical seam vertices
            v = v * voxel + (xs[0], ys[0], zs[0])
            writers[n].add(v, faces)
        del f
        print(f"  x {xs[i0]:7.1f} .. {xs[i1]:7.1f} mm  ({time.time() - t0:5.0f} s)")
    for w in writers.values():
        w.close()
        print(f"{os.path.basename(w.path)}: {w.n:,} triangles -> {w.path}")
    return [w.path for w in writers.values()]


def check(path):
    """Watertightness and connectivity check (pip install trimesh).
    Run on a coarse export (voxel >= 0.4) - full-size meshes are huge."""
    import trimesh
    m = trimesh.load(path, process=False)
    m.merge_vertices(digits_vertex=4)
    _, c = np.unique(m.edges_sorted, axis=0, return_counts=True)
    parts = m.split(only_watertight=False)
    print(f"{os.path.basename(path)}: {len(m.faces):,} tris, bad edges "
          f"{(c != 2).sum()}, watertight {m.is_watertight}, "
          f"bodies {len(parts)}, volume {abs(m.volume) / 1e3:.1f} cm3")


def preview(out_png="hx_preview.png", res=0.25):
    """Colour sections (grey solid, red hot, blue cold) - fast sanity check."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["white", "#555b66", "#E4572E", "#2E86AB"])

    def lab(f):
        m = np.zeros(f["HX_solid"].shape, int)
        m[f["HX_solid"] < 0] = 1
        m[f["Hot_fluid"] < 0] = 2
        m[f["Cold_fluid"] < 0] = 3
        return m

    g = np.arange(-R - 3, R + 3, res)
    gz = np.arange(-LEG_H - 2, HB + 2, res)
    fig, axs = plt.subplots(1, 4, figsize=(20, 5.6))
    for ax, zc in [(axs[0], (Z_APEX + ZP) / 2), (axs[1], (ZP + HB) / 2)]:
        Xs, Ys = np.meshgrid(g, g, indexing="ij")
        ax.imshow(lab(fields(Xs, Ys, np.full_like(Xs, zc))).T, origin="lower",
                  extent=(g[0], g[-1], g[0], g[-1]), cmap=cmap, vmin=0, vmax=3)
        ax.set_title(f"z = {zc:.0f} mm")
    for ax, ang in [(axs[2], 45), (axs[3], 135)]:
        S, Zg = np.meshgrid(g, gz, indexing="ij")
        t = np.radians(ang)
        ax.imshow(lab(fields(S * np.cos(t), S * np.sin(t), Zg)).T, origin="lower",
                  extent=(g[0], g[-1], gz[0], gz[-1]), cmap=cmap, vmin=0, vmax=3)
        ax.set_title(f"section through {ang} / {ang + 180} deg pillars")
    plt.tight_layout()
    plt.savefig(out_png, dpi=110)
    print(f"preview -> {out_png}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[2])
    ap.add_argument("--voxel", type=float, default=VOXEL)
    ap.add_argument("--out", default=OUT_DIR)
    ap.add_argument("--no-fluids", action="store_true")
    ap.add_argument("--preview", action="store_true", help="PNG sections only")
    ap.add_argument("--check", action="store_true", help="check STLs after export")
    a = ap.parse_args()
    if a.preview:
        preview()
    else:
        paths = export(a.voxel, a.out, not a.no_fluids)
        if a.check:
            for p in paths:
                check(p)
