"""
gyroid_hx_model.py
------------------------------------------------------------------------
Two-fluid Gyroid COUNTER-FLOW heat exchanger with ALL FOUR PORTS ON THE
SAME SIDE FACE (+y) as 3/8"-18 NPT MALE threaded nipples.
Builds geometry as implicit fields and exports STL (import into nTop with
Import Mesh -> Implicit Body from Mesh, SolidWorks, or a slicer).

Install:  pip install numpy scikit-image   (optional check(): trimesh)
Run:      python gyroid_hx_model.py

Layout (origin at body centre, mm):
    x = length / counter-flow axis, y = width, z = height
    All nipples point +y and sit on the y = +W/2 face.
    HOT : in  x- end, upper half -> out x+ end, upper half
    COLD: in  x+ end, lower half -> out x- end, lower half
    Each end has a split header: hot plenum on top, cold plenum below,
    separated by a divider plate. The gyroid core itself runs end to end,
    so the two fluids flow in opposite directions (counter-flow).
------------------------------------------------------------------------
"""

import os
import time
import numpy as np
from skimage.measure import marching_cubes

# ======================================================================
# 1. PARAMETERS (mm)
# ======================================================================
VOXEL = 0.2           # 0.5 quick check, 0.2 default, 0.12 crisp threads
CHUNK = 80            # x-slices per chunk (lower if RAM is short)
EXPORT_FLUIDS = True
OUT_DIR = "hx_output"

# Envelope (body only; nipples add NIPPLE_LEN on the +y side)
LB = 130.0            # body length (x)
W = 100.0             # body width (y)
H = 50.0              # body height (z); two 21 mm collars must fit stacked
SHELL_T = 1.5
R_EDGE = 5.0

# Gyroid core
CELL = (9.0, 6.0, 6.0)
T_MIN = 0.5
T_MAX = 0.8           # thicker walls along the port (+y) side, where flow is
                      # fastest. Set T_MAX = T_MIN to turn grading off.
G_W = 30.0            # distance over which the thickening fades
LC = 84.0             # core length

# Headers
H_MIN = 5.0           # plenum depth on the far (-y) side
DV = 1.5              # divider plate between hot and cold plenums
S_SEAL = 1.5          # seal slab capping the wrong fluid at each core face

# 3/8"-18 NPT external thread (ASME B1.20.1 basic values)
IN = 25.4
NPT_P = IN / 18.0
NPT_E0 = 0.61201 * IN
NPT_H = 0.8 * NPT_P
NPT_TAPER = 1.0 / 16.0
NIPPLE_LEN = 14.0
BORE_D = 10.0
COLLAR_D = 21.0

# Derived
H_MAX = LB / 2 - SHELL_T - LC / 2                # plenum depth at port side
X_P = LC / 2 + H_MAX / 2                         # |x| of port centres
Z_P = (H / 2 - SHELL_T + DV / 2) / 2             # |z| of port centres
CM = np.mean(CELL)
Y_TIP = W / 2 + NIPPLE_LEN

# Ports: (name, sx = end, sz = upper/lower, fluid). All on the +y face.
PORTS = [
    ("hot_in",   -1, +1, "hot"),
    ("hot_out",  +1, +1, "hot"),
    ("cold_in",  +1, -1, "cold"),
    ("cold_out", -1, -1, "cold"),
]

# ======================================================================
# 2. FIELD HELPERS  (value < 0 means inside the body)
# ======================================================================
def smin(a, b, k):
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0.0, 1.0)
    return b * (1 - h) + a * h - k * h * (1 - h)


def round_box(x, y, z, hx, hy, hz, r):
    qx = np.abs(x) - (hx - r)
    qy = np.abs(y) - (hy - r)
    qz = np.abs(z) - (hz - r)
    outside = np.sqrt(np.maximum(qx, 0) ** 2 + np.maximum(qy, 0) ** 2
                      + np.maximum(qz, 0) ** 2)
    inside = np.minimum(np.maximum(qx, np.maximum(qy, qz)), 0)
    return outside + inside - r


def gyroid(x, y, z):
    kx, ky, kz = (2 * np.pi / c for c in CELL)
    sx, cx = np.sin(kx * x), np.cos(kx * x)
    sy, cy = np.sin(ky * y), np.cos(ky * y)
    sz, cz = np.sin(kz * z), np.cos(kz * z)
    g = sx * cy + sy * cz + sz * cx
    gx = kx * (cx * cy - sz * sx)
    gy = ky * (-sx * sy + cy * cz)
    gz = kz * (-sy * sz + cz * cx)
    grad = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)
    return g / np.maximum(grad, 0.25 * 2 * np.pi / CM)


def wall_thickness(y):
    """Thicker walls toward the port side (+y)."""
    return T_MIN + (T_MAX - T_MIN) * np.exp(-((W / 2 - y) / G_W) ** 2)


def port_local(x, y, z, sx, sz):
    """Nipple axis is +y. t = outward axial, r = radius, theta in a
    right-handed frame about +y (z x x = y)."""
    dx = x - sx * X_P
    dz = z - sz * Z_P
    return y, np.hypot(dx, dz), np.arctan2(dx, dz)


def npt_male(x, y, z, sx, sz):
    """Right-hand 3/8-18 NPT external thread on a tapered nipple."""
    t, r, th = port_local(x, y, z, sx, sz)
    s = Y_TIP - t                                  # distance from tip
    rp = (NPT_E0 + np.maximum(s, 0) * NPT_TAPER) / 2
    w = np.mod(t - NPT_P * th / (2 * np.pi), NPT_P) - NPT_P / 2
    r_prof = rp + np.tan(np.radians(60)) * (NPT_P / 4 - np.abs(w))
    r_prof = np.clip(r_prof, rp - NPT_H / 2, rp + NPT_H / 2)
    f = 0.5 * (r - r_prof)
    f = np.maximum(f, -s)
    f = np.maximum(f, (r - rp) - (s - 0.8))
    return np.maximum(f, s - (NIPPLE_LEN + 4.0))


def collar(x, y, z, sx, sz):
    t, r, _ = port_local(x, y, z, sx, sz)
    a0, a1 = W / 2 - 9.0, W / 2 + 1.5
    return np.maximum(r - COLLAR_D / 2, np.maximum(a0 - t, t - a1))


def bore(x, y, z, sx, sz):
    t, r, _ = port_local(x, y, z, sx, sz)
    a0 = W / 2 - SHELL_T - 1.0
    return np.maximum(r - BORE_D / 2, a0 - t)


def plenum(x, y, z, sx, sz, inner):
    """Tapered header: deepest at the port side (+y), H_MIN at -y."""
    depth = H_MIN + (H_MAX - H_MIN) * (W / 2 + y) / W
    u = sx * x
    p = np.maximum(LC / 2 - 0.5 - u, u - (LC / 2 + depth))
    p = np.maximum(p, DV / 2 - sz * z)
    return np.maximum(p, inner)


# ======================================================================
# 3. FIELDS
# ======================================================================
def fields(x, y, z):
    body = round_box(x, y, z, LB / 2, W / 2, H / 2, R_EDGE)
    inner = body + SHELL_T

    G = gyroid(x, y, z)
    T = wall_thickness(y)
    ax = np.abs(x)
    in_core = ax - LC / 2
    slab_free = ax - (LC / 2 - S_SEAL)
    hot_core = np.maximum(np.maximum(T / 2 - G, in_core),
                          np.maximum(np.minimum(slab_free, DV / 2 - z), inner))
    cold_core = np.maximum(np.maximum(G + T / 2, in_core),
                           np.maximum(np.minimum(slab_free, DV / 2 + z), inner))

    hot, cold = hot_core, cold_core
    outer = body
    for _, sx, sz, fluid in PORTS:
        outer = smin(outer, collar(x, y, z, sx, sz), 3.0)
        outer = np.minimum(outer, npt_male(x, y, z, sx, sz))
        v = np.minimum(plenum(x, y, z, sx, sz, inner), bore(x, y, z, sx, sz))
        if fluid == "hot":
            hot = np.minimum(hot, v)
        else:
            cold = np.minimum(cold, v)

    solid = np.maximum(outer, -np.minimum(hot, cold))
    return {"HX_solid": solid,
            "Hot_fluid": np.maximum(hot, outer),
            "Cold_fluid": np.maximum(cold, outer)}


# ======================================================================
# 4. MESHING -> BINARY STL (chunked)
# ======================================================================
class StlWriter:
    def __init__(self, path):
        self.f = open(path, "wb")
        self.f.write(b"gyroid counter-flow HX".ljust(80, b" "))
        self.f.write(np.uint32(0).tobytes())
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


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    pad = 1.0 + 0.137 * VOXEL
    xs = np.arange(-LB / 2 - pad, LB / 2 + pad + VOXEL, VOXEL)
    ys = np.arange(-W / 2 - pad, Y_TIP + pad + VOXEL, VOXEL)
    zs = np.arange(-H / 2 - pad, H / 2 + pad + VOXEL, VOXEL)
    print(f"Grid {len(xs)} x {len(ys)} x {len(zs)} "
          f"= {len(xs) * len(ys) * len(zs) / 1e6:.0f} M voxels")
    Y, Z = np.meshgrid(ys, zs, indexing="ij")

    names = ["HX_solid"] + (["Hot_fluid", "Cold_fluid"] if EXPORT_FLUIDS else [])
    writers = {n: StlWriter(os.path.join(OUT_DIR, f"{n}.stl")) for n in names}
    t0 = time.time()
    for i0 in range(0, len(xs) - 1, CHUNK):
        i1 = min(i0 + CHUNK, len(xs) - 1)
        X = xs[i0:i1 + 1][:, None, None]
        f = fields(X, Y[None], Z[None])
        for n in names:
            vol = f[n].astype(np.float32)
            vol[vol == 0] = 1e-6
            if vol.min() > 0 or vol.max() < 0:
                continue
            v, faces, _, _ = marching_cubes(vol, level=0.0)
            v = v.astype(np.float64)
            v[:, 0] += i0
            v = v * VOXEL + (xs[0], ys[0], zs[0])
            writers[n].add(v, faces)
        del f
        print(f"  x {xs[i0]:7.1f} .. {xs[i1]:7.1f} mm  ({time.time() - t0:5.0f} s)")

    for n, w in writers.items():
        w.close()
        print(f"{n}: {w.n:,} triangles -> {w.f.name}")


def check(path):
    """Optional mesh check (pip install trimesh)."""
    import trimesh
    m = trimesh.load(path, process=False)
    m.merge_vertices(digits_vertex=5)
    _, c = np.unique(m.edges_sorted, axis=0, return_counts=True)
    print(f"{path}: {len(m.faces):,} triangles, "
          f"bad edges: {(c != 2).sum()}, watertight: {m.is_watertight}")


if __name__ == "__main__":
    main()
