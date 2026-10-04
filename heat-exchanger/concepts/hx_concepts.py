"""
hx_concepts.py
------------------------------------------------------------------------
Concept geometries for compact two-fluid gyroid heat exchangers with
different port / header layouts. Each concept is an implicit model
(value < 0 = inside) rendered as an exterior view plus colour sections
(grey solid, red hot, blue cold). Walls are drawn thicker than a
printable design so they show at concept resolution.

    A  Inline cartridge   : axial hot ports on the pipe line, tangential
                            cold ports, concentric split headers.
                            Pure counter-flow.
    B  Radial disc        : hot enters at the hub and flows outwards,
                            cold enters at the rim and flows inwards.
                            Radial counter-flow, very low height.
    C  U-turn block       : all four ports on one face; both fluids make
                            two passes around a central baffle.
                            Counter-flow in each pass.
    D  Coaxial nested     : ONE coaxial connection per end (hot in the
                            centre tube, cold in the outer annulus);
                            two stacked header levels route the cold
                            annulus over the hot header. Pure counter-flow.

Run:  python hx_concepts.py      -> figures/concept_*.png, metrics printed
------------------------------------------------------------------------
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes

CELL = 10.0
T = 1.6                 # display wall thickness
SHELL = 2.0
S_SEAL = 1.5
DV = 2.0
os.makedirs("figures", exist_ok=True)


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def smin(a, b, k):
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0, 1)
    return b * (1 - h) + a * h - k * h * (1 - h)


def gyroid(x, y, z):
    k = 2 * np.pi / CELL
    sx, cx, sy, cy, sz, cz = (np.sin(k * x), np.cos(k * x), np.sin(k * y),
                              np.cos(k * y), np.sin(k * z), np.cos(k * z))
    g = sx * cy + sy * cz + sz * cx
    grad = k * np.sqrt((cx * cy - sz * sx) ** 2 + (cy * cz - sx * sy) ** 2
                       + (cz * cx - sy * sz) ** 2)
    return g / np.maximum(grad, 0.25 * k)


def rbox(x, y, z, hx, hy, hz, r):
    q = [np.abs(x) - (hx - r), np.abs(y) - (hy - r), np.abs(z) - (hz - r)]
    out = np.sqrt(sum(np.maximum(c, 0) ** 2 for c in q))
    return out + np.minimum(np.maximum(q[0], np.maximum(q[1], q[2])), 0) - r


def rcyl(r, z, R, z0, z1, rr):
    """Cylinder (axis z) radius R from z0 to z1 with rounded edges."""
    hh, zc = (z1 - z0) / 2, z - (z0 + z1) / 2
    qx, qz = r - (R - rr), np.abs(zc) - (hh - rr)
    return (np.hypot(np.maximum(qx, 0), np.maximum(qz, 0))
            + np.minimum(np.maximum(qx, qz), 0) - rr)


def tube(p, axis, c, r, t0, t1):
    """Solid rod along unit axis through point c, radius r, t in [t0, t1]."""
    d = [p[i] - c[i] for i in range(3)]
    t = sum(d[i] * axis[i] for i in range(3))
    rad = np.sqrt(np.maximum(sum(d[i] ** 2 for i in range(3)) - t ** 2, 0))
    return np.maximum(rad - r, np.maximum(t0 - t, t - t1))


def assemble(outer, core, G, hot_open, cold_open, hot_extra, cold_extra):
    """Two gyroid labyrinths inside `core`. In the seal slabs a labyrinth
    is kept only where its *_open field < 0. *_extra: plenums / bores."""
    inner = outer + SHELL
    hot = np.maximum(np.maximum(T / 2 - G, np.maximum(core, inner)), hot_open)
    cold = np.maximum(np.maximum(G + T / 2, np.maximum(core, inner)), cold_open)
    for v in hot_extra:
        hot = np.minimum(hot, v)
    for v in cold_extra:
        cold = np.minimum(cold, v)
    hot, cold = np.maximum(hot, outer), np.maximum(cold, outer)
    solid = np.maximum(outer, -np.minimum(hot, cold))
    return dict(outer=outer, solid=solid, hot=hot, cold=cold,
                core=np.maximum(core, inner))


# ----------------------------------------------------------------------
# A  Inline cartridge (axis z): axial hot, tangential cold
# ----------------------------------------------------------------------
A = dict(R=35.0, ZC=55.0, ZE=70.0, NOZ=84.0, RS=23.5)


def concept_A(x, y, z):
    R, ZC, ZE, NOZ, RS = A["R"], A["ZC"], A["ZE"], A["NOZ"], A["RS"]
    r = np.hypot(x, y)
    outer = rcyl(r, z, R, -ZE - 2, ZE + 2, 8.0)
    p = (x, y, z)
    for s in (1, -1):                                  # axial hot nozzles
        outer = smin(outer, tube(p, (0, 0, s), (0, 0, 0), 9.0, ZE - 6, NOZ), 4)
    for s, yy in ((1, 15.0), (-1, -15.0)):             # tangential cold nozzles
        outer = smin(outer, tube(p, (s, 0, 0), (0, yy, s * (ZC + 7.5)),
                                 8.0, 0, R + 12), 4)
    inner = outer + SHELL
    G = gyroid(x, y, z)
    core = np.abs(z) - ZC
    slab = (ZC - S_SEAL) - np.abs(z)                   # < 0 inside seal slabs
    hot_open = np.minimum(-slab, r - (RS - DV / 2))
    cold_open = np.minimum(-slab, (RS + DV / 2) - r)
    hot_x, cold_x = [], []
    for s in (1, -1):
        zz = s * z
        plen = np.maximum(np.maximum(ZC - 0.5 - zz, zz - (ZE - SHELL)), inner)
        hot_x.append(np.maximum(plen, r - (RS - DV / 2)))
        cold_x.append(np.maximum(plen, (RS + DV / 2) - r))
        hot_x.append(tube(p, (0, 0, s), (0, 0, 0), 6.0, ZE - 8, NOZ + 1))
    for s, yy in ((1, 15.0), (-1, -15.0)):
        cold_x.append(tube(p, (s, 0, 0), (0, yy, s * (ZC + 7.5)), 5.0, RS + 2, R + 13))
    return assemble(outer, core, G, hot_open, cold_open, hot_x, cold_x)


# ----------------------------------------------------------------------
# B  Radial counter-flow disc (axis z)
# ----------------------------------------------------------------------
B = dict(R=62.0, H=13.0, RH=10.0, RC=50.0, RR=57.0)


def concept_B(x, y, z):
    R, H, RH, RC, RR = B["R"], B["H"], B["RH"], B["RC"], B["RR"]
    r = np.hypot(x, y)
    p = (x, y, z)
    outer = rcyl(r, z, R, -H, H, 6.0)
    outer = smin(outer, tube(p, (0, 0, 1), (0, 0, 0), 9.0, 0, H + 14), 5)
    outer = smin(outer, tube(p, (0, 0, -1), (0, 0, 0), 9.0, 0, H + 14), 5)
    outer = smin(outer, tube(p, (0, 1, 0), (R - 5, 0, 6.0), 6.0, 0, 26), 3)
    outer = smin(outer, tube(p, (0, -1, 0), (-(R - 5), 0, -6.0), 6.0, 0, 26), 3)
    inner = outer + SHELL
    G = gyroid(x, y, z)
    core = np.maximum(RH - r, r - RC)
    slab = np.minimum(r - RH - S_SEAL, RC - S_SEAL - r)   # < 0 near hub / rim
    hot_open = np.minimum(-slab, DV / 2 - z)               # hot: upper half
    cold_open = np.minimum(-slab, DV / 2 + z)              # cold: lower half
    hub = np.maximum(r - (RH + 0.5), inner)
    ring = np.maximum(np.maximum(RC - 0.5 - r, r - RR), inner)
    hot_x = [np.maximum(hub, DV / 2 - z), np.maximum(ring, DV / 2 - z),
             tube(p, (0, 0, 1), (0, 0, 0), 6.0, 0, H + 15),
             tube(p, (0, 1, 0), (R - 5, 0, 6.0), 4.0, 0, 27)]
    cold_x = [np.maximum(hub, DV / 2 + z), np.maximum(ring, DV / 2 + z),
              tube(p, (0, 0, -1), (0, 0, 0), 6.0, 0, H + 15),
              tube(p, (0, -1, 0), (-(R - 5), 0, -6.0), 4.0, 0, 27)]
    return assemble(outer, core, G, hot_open, cold_open, hot_x, cold_x)


# ----------------------------------------------------------------------
# C  U-turn block: all ports on the top face
# ----------------------------------------------------------------------
C = dict(HX=42.0, HY=27.0, HZ=62.0, ZT=46.0, ZB=-44.0)


def concept_C(x, y, z):
    HX, HY, HZ, ZT, ZB = C["HX"], C["HY"], C["HZ"], C["ZT"], C["ZB"]
    p = (x, y, z)
    ports = [(-20, -13, "hot"), (-20, 13, "cold"), (20, -13, "hot"), (20, 13, "cold")]
    outer = rbox(x, y, z, HX, HY, HZ, 6.0)
    for px, py, _ in ports:
        outer = smin(outer, tube(p, (0, 0, 1), (px, py, 0), 7.0, 0, HZ + 14), 4)
    inner = outer + SHELL
    G = gyroid(x, y, z)
    baffle = np.maximum(np.abs(x) - 1.0, ZB - z)          # solid, open below ZB
    core = np.maximum(z - ZT, -baffle)
    slab = (ZT - S_SEAL) - z
    hot_open = np.minimum(slab * -1, y + DV / 2)           # hot: y < 0 at top
    cold_open = np.minimum(slab * -1, DV / 2 - y)          # cold: y > 0 at top
    hot_x, cold_x = [], []
    for px, py, f in ports:
        plen = np.maximum(np.maximum(ZT - 0.5 - z, inner),
                          np.maximum(DV / 2 - np.sign(px) * x, DV / 2 - np.sign(py) * y))
        bore = tube(p, (0, 0, 1), (px, py, 0), 4.5, ZT, HZ + 15)
        (hot_x if f == "hot" else cold_x).extend([plen, bore])
    out = assemble(outer, core, G, hot_open, cold_open, hot_x, cold_x)
    out["solid"] = np.minimum(out["solid"], np.maximum(baffle, outer))
    out["hot"] = np.maximum(out["hot"], -baffle)
    out["cold"] = np.maximum(out["cold"], -baffle)
    return out


# ----------------------------------------------------------------------
# D  Coaxial nested: one coaxial connection per end
# ----------------------------------------------------------------------
D = dict(R=35.0, ZC=50.0, ZE=68.0, NOZ=86.0, RS=23.5, ZS=60.0)


def concept_D(x, y, z):
    """Each end: centre tube (hot) inside an annulus (cold). The hot tube
    drops through a sleeve into a central hot header next to the core; the
    cold annulus feeds an upper distribution disk that sits over the hot
    header (divider plate between them) and overflows into the outer ring
    header next to the core."""
    R, ZC, ZE, NOZ, RS, ZS = D["R"], D["ZC"], D["ZE"], D["NOZ"], D["RS"], D["ZS"]
    r = np.hypot(x, y)
    p = (x, y, z)
    outer = rcyl(r, z, R, -ZE - 2, ZE + 2, 8.0)
    for s in (1, -1):
        outer = smin(outer, tube(p, (0, 0, s), (0, 0, 0), 14.0, ZE - 6, NOZ), 6)
    inner = outer + SHELL
    G = gyroid(x, y, z)
    core = np.abs(z) - ZC
    slab = (ZC - S_SEAL) - np.abs(z)
    hot_open = np.minimum(-slab, r - (RS - DV / 2))
    cold_open = np.minimum(-slab, (RS + DV / 2) - r)
    hot_x, cold_x = [], []
    for s in (1, -1):
        zz = s * z
        # lower level (next to the core): hot centre disk, cold outer ring
        lvl = np.maximum(np.maximum(ZC - 0.5 - zz, zz - (ZS - DV / 2)), inner)
        hot_x.append(np.maximum(lvl, r - (RS - DV / 2)))
        cold_x.append(np.maximum(lvl, (RS + DV / 2) - r))
        # upper level: cold distribution disk fed by the annulus, overlapping
        # the outer ring; the hot tube passes through it in a sleeve
        up = np.maximum(np.maximum(ZS + DV / 2 - zz, zz - (ZE - SHELL)), inner)
        sleeve = 8.5 - r                                   # keep out of hot tube
        cold_x.append(np.maximum(np.maximum(up, sleeve), r - (RS + 6)))
        hot_x.append(tube(p, (0, 0, s), (0, 0, 0), 6.5, ZC, NOZ + 1))
        ann = np.maximum(np.maximum(r - 11.5, 8.5 - r), np.maximum(ZS + 2 - zz, zz - (NOZ + 1)))
        cold_x.append(ann)
    return assemble(outer, core, G, hot_open, cold_open, hot_x, cold_x)


CONCEPTS = {
    "A": ("Inline cartridge: axial hot, tangential cold (counter-flow)",
          concept_A, (-48, 48, -48, 48, -90, 90)),
    "B": ("Radial disc: hub-to-rim counter-flow, 26 mm tall",
          concept_B, (-90, 90, -90, 90, -30, 30)),
    "C": ("U-turn block: all four ports on the top face",
          concept_C, (-45, 45, -30, 30, -65, 80)),
    "D": ("Coaxial nested: one coaxial connection per end",
          concept_D, (-40, 40, -40, 40, -90, 90)),
}


# ----------------------------------------------------------------------
# rendering
# ----------------------------------------------------------------------
GREY = (0.62, 0.66, 0.72)
CMAP = ListedColormap(["white", "#555b66", "#E4572E", "#2E86AB"])


def shade(v, f):
    tri = v[f]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
    L = np.array([-0.4, -0.6, 0.7]); L /= np.linalg.norm(L)
    s = 0.35 + 0.65 * np.abs(n @ L)
    return tri, np.c_[np.outer(s, GREY), np.ones(len(s))]


def labels(f):
    m = np.zeros(f["solid"].shape, int)
    m[f["solid"] < 0] = 1; m[f["hot"] < 0] = 2; m[f["cold"] < 0] = 3
    return m


def section(ax, fn, plane, lim, res=0.3, title=""):
    x0, x1, y0, y1, z0, z1 = lim
    if plane == "xz":
        a, b = np.arange(x0, x1, res), np.arange(z0, z1, res)
        A_, B_ = np.meshgrid(a, b, indexing="ij")
        f = fn(A_, np.full_like(A_, 0.0), B_)
    elif plane == "yz":
        a, b = np.arange(y0, y1, res), np.arange(z0, z1, res)
        A_, B_ = np.meshgrid(a, b, indexing="ij")
        f = fn(np.full_like(A_, 0.0), A_, B_)
    else:                                              # ("xy", zc)
        a, b = np.arange(x0, x1, res), np.arange(y0, y1, res)
        A_, B_ = np.meshgrid(a, b, indexing="ij")
        f = fn(A_, B_, np.full_like(A_, plane[1]))
    ax.imshow(labels(f).T, origin="lower", extent=(a[0], a[-1], b[0], b[-1]),
              cmap=CMAP, vmin=0, vmax=3, interpolation="nearest")
    ax.set_title(title, fontsize=10)
    ax.tick_params(labelsize=8)


def metrics(fn, lim, step=0.8):
    x0, x1, y0, y1, z0, z1 = lim
    xs, ys, zs = (np.arange(x0, x1, step), np.arange(y0, y1, step),
                  np.arange(z0, z1, step))
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    f = fn(X, Y, Z)
    dv = step ** 3 / 1e3                               # cm3
    env = (f["outer"] < 0).sum() * dv
    core = (f["core"] < 0).sum() * dv
    solid = (f["solid"] < 0).sum() * dv
    leak = ((f["hot"] < 0) & (f["cold"] < 0)).sum()
    assert leak == 0, f"hot and cold domains overlap in {leak} voxels"
    return env, core, solid, f, (xs, ys, zs)


def render(key, step=0.8):
    title, fn, lim = CONCEPTS[key]
    env, core, solid, f, (xs, ys, zs) = metrics(fn, lim, step)
    o = np.array([xs[0], ys[0], zs[0]])
    v, fc, _, _ = marching_cubes(f["outer"].astype(np.float32), 0.0, spacing=(step,) * 3)
    ext = shade(v + o, fc)
    X = np.meshgrid(xs, ys, zs, indexing="ij")
    cut = np.maximum(f["solid"], -np.maximum(-X[0], X[1]))
    del X, f
    v, fc, _, _ = marching_cubes(cut.astype(np.float32), 0.0, spacing=(step,) * 3)
    cw = shade(v + o, fc)

    fig = plt.figure(figsize=(20, 5.4), dpi=105)
    for i, (tri, view, t) in enumerate([(ext, (20, -55), "Exterior"),
                                        (cw, (22, -45), "Cut-away")]):
        ax = fig.add_subplot(1, 4, i + 1, projection="3d")
        ax.add_collection3d(Poly3DCollection(tri[0], facecolors=tri[1], edgecolor="none"))
        ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3]); ax.set_zlim(lim[4], lim[5])
        ax.set_box_aspect((lim[1] - lim[0], lim[3] - lim[2], lim[5] - lim[4]))
        ax.view_init(*view); ax.set_axis_off(); ax.set_title(t, fontsize=10)
    secs = {"A": [("xz", "Section y = 0"), ("yz", "Section x = 0")],
            "B": [("xz", "Section y = 0 (height exaggerated x2)"), (("xy", 6.0), "z = +6 (hot level)")],
            "C": [("xz", "Section y = 0 (baffle, U-turn)"), ("yz", "Section x = -20 (left pass)")],
            "D": [("xz", "Section y = 0"), (("xy", D["ZS"] + 4), "Upper header level: cold disk around the hot tube")]}[key]
    for j, (pl, t) in enumerate(secs):
        ax = fig.add_subplot(1, 4, 3 + j)
        if key == "C" and pl == "yz":
            fn2 = lambda a, b, c: concept_C(np.full_like(a, -20.0) + a * 0, b, c)
            section(ax, fn2, "yz", lim, title=t)
        else:
            section(ax, fn, pl, lim, title=t)
        if key == "B" and pl == "xz":
            ax.set_aspect(2.0)
    fig.suptitle(f"Concept {key}: {title}", fontsize=13)
    plt.tight_layout(rect=(0, 0, 1, 0.94))
    plt.savefig(f"figures/concept_{key}.png"); plt.close()
    return env, core, solid


def effectiveness_chart():
    """epsilon-NTU for Cr = 1 (balanced streams), standard closed forms."""
    ntu = np.linspace(0.01, 6, 200)
    cf = ntu / (1 + ntu)                                   # counter-flow, Cr = 1
    pf = 0.5 * (1 - np.exp(-2 * ntu))                      # parallel, Cr = 1
    xf = 1 - np.exp(ntu ** 0.22 * (np.exp(-ntu ** 0.78) - 1))   # cross, both unmixed
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=110)
    for yv, lab, c in [(cf, "Counter-flow (A, B, D; C per pass)", "#2E86AB"),
                       (xf, "Cross-flow, unmixed (approx.)", "#8E6C8A"),
                       (pf, "Parallel flow", "#E4572E")]:
        ax.plot(ntu, yv, lw=2.2, label=lab, color=c)
    ax.set_xlabel("NTU = UA / C_min"); ax.set_ylabel("effectiveness  ε")
    ax.set_title("ε-NTU, balanced streams (C_r = 1)")
    ax.grid(alpha=0.3); ax.legend(frameon=False); ax.set_ylim(0, 1)
    plt.tight_layout(); plt.savefig("figures/eps_ntu.png"); plt.close()
    for n in (1, 2, 3, 5):
        print(f"NTU {n}: counter {n / (1 + n):.3f}  "
              f"cross {1 - np.exp(n ** 0.22 * (np.exp(-n ** 0.78) - 1)):.3f}  "
              f"parallel {0.5 * (1 - np.exp(-2 * n)):.3f}")


if __name__ == "__main__":
    import sys
    keys = sys.argv[1:] or list(CONCEPTS)
    for k in keys:
        env, core, solid = render(k)
        print(f"{k}: envelope {env:7.1f} cm3  core {core:7.1f} cm3  "
              f"core/envelope {core / env:5.1%}  solid {solid:6.1f} cm3")
    effectiveness_chart()
