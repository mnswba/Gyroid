"""
cyl_hx_concept.py
------------------------------------------------------------------------
CONCEPT render of a cylindrical two-fluid gyroid heat exchanger with four
arched legs (ports) blended smoothly into the cylinder underside.
Geometry is built as implicit fields (value < 0 = inside) and rendered to
PNG. This is a visual concept, not the print/CFD export: walls are drawn
thicker than a printable design would use so they show at this resolution.

Layout (mm, origin at centre of cylinder base):
    Cylinder: radius R, from z = 0 to z = HB, top edge rounded.
    Legs at 45, 135, 225, 315 deg, LEG_H = 40 mm long (z = -40 to 0).
    Each leg is a slow trumpet taper out of the cylinder underside that
    ends in a 3/8"-18 NPT MALE threaded tip (ASME B1.20.1 basic sizes).
    HOT  legs: 45 and 225 deg (in / out)   COLD legs: 135 and 315 deg.
    Bottom header split into four quadrant plenums, one per leg.
    Gyroid core above the headers; a seal slab lets each fluid's
    labyrinth open only into its own quadrants.
------------------------------------------------------------------------
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes

R, HB, R_TOP = 70.0, 90.0, 10.0       # cylinder radius, height, top rounding
SHELL = 2.0
LEG_H = 40.0                          # leg length below the cylinder
R_BASE = 22.0                         # leg radius where it leaves the body
R_NECK = 11.0                         # leg radius just above the thread
LEG_POS = R - R_BASE                  # leg base flush with the outer wall
BLEND = 15.0                          # extra fillet at the leg/body joint
BORE_R = 5.0                          # 10 mm flow bore

# 3/8"-18 NPT external thread (ASME B1.20.1 basic values)
IN = 25.4
NPT_P = IN / 18.0                     # 1.411 mm pitch
NPT_E0 = 0.61201 * IN                 # 15.545 mm pitch dia at the tip
NPT_H = 0.8 * NPT_P                   # truncated thread height
NPT_TAPER = 1.0 / 16.0                # on diameter
THREAD_L = 12.0                       # threaded length (L2 = 10.2 mm)
ZP = 14.0                             # header (plenum) height
DV = 2.0                              # divider between quadrant plenums
S_SEAL = 1.5
CELL = 12.0
T_WALL = 1.6                          # display thickness (print ~0.5-0.8)
LEGS = [(45, "hot"), (135, "cold"), (225, "hot"), (315, "cold")]


def smin(a, b, k):
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0, 1)
    return b * (1 - h) + a * h - k * h * (1 - h)


def cylinder(x, y, z):
    """Vertical cylinder with rounded top edge (sharp bottom, legs blend there)."""
    r = np.hypot(x, y)
    zc = z - HB / 2
    qx, qz = r - (R - R_TOP), np.abs(zc) - (HB / 2 - R_TOP)
    top = np.hypot(np.maximum(qx, 0), np.maximum(qz, 0)) \
        + np.minimum(np.maximum(qx, qz), 0) - R_TOP
    return np.maximum(top, np.maximum(r - R, -z))


def leg_xy(a):
    t = np.radians(a)
    return LEG_POS * np.cos(t), LEG_POS * np.sin(t)


def leg(x, y, z, a):
    """Tapered leg: trumpet flare from R_BASE at the body to R_NECK at the
    thread shoulder, then a 3/8-18 NPT male thread down to the tip."""
    lx, ly = leg_xy(a)
    dx, dy = x - lx, y - ly
    rho, t = np.hypot(dx, dy), -z                  # t = depth below body
    t1 = LEG_H - THREAD_L                          # end of the flare
    u = np.clip(t / t1, 0, 1)
    r_prof = R_NECK + (R_BASE - R_NECK) * (1 - u) ** 2
    slope = 2 * (R_BASE - R_NECK) * (1 - u) / t1
    flare = (rho - r_prof) / np.sqrt(1 + slope ** 2)
    flare = np.maximum(flare, np.maximum(z - 1.0, t - t1))

    # NPT thread, right-handed about the outward (-z) axis
    s_tip = LEG_H - t                              # distance from the tip
    th = np.arctan2(dx, dy)
    rp = (NPT_E0 + np.maximum(s_tip, 0) * NPT_TAPER) / 2
    w = np.mod(t - NPT_P * th / (2 * np.pi), NPT_P) - NPT_P / 2
    r_thr = rp + np.tan(np.radians(60)) * (NPT_P / 4 - np.abs(w))
    r_thr = np.clip(r_thr, rp - NPT_H / 2, rp + NPT_H / 2)
    thr = 0.5 * (rho - r_thr)
    thr = np.maximum(thr, -s_tip)                  # cut at the tip
    thr = np.maximum(thr, (rho - rp) - (s_tip - 0.8))   # tip chamfer
    thr = np.maximum(thr, s_tip - (THREAD_L + 1.0))     # runs into shoulder
    return np.minimum(flare, thr)


def sector(x, y, a):
    t = np.radians(a)
    u = x * np.cos(t) + y * np.sin(t)
    v = -x * np.sin(t) + y * np.cos(t)
    return (np.abs(v) - u) / np.sqrt(2) + DV / 2


def gyroid(x, y, z):
    k = 2 * np.pi / CELL
    g = (np.sin(k * x) * np.cos(k * y) + np.sin(k * y) * np.cos(k * z)
         + np.sin(k * z) * np.cos(k * x))
    return g / (k * 1.0)              # rough distance scaling (|grad| ~ k)


def fields(x, y, z):
    outer = cylinder(x, y, z)
    for a, _ in LEGS:
        outer = smin(outer, leg(x, y, z, a), BLEND)
    # trim the blend to the cylinder wall: legs stay flush, the blend only
    # forms the arches underneath instead of bulging outwards
    outer = np.maximum(outer, np.hypot(x, y) - R)
    inner = outer + SHELL

    G = gyroid(x, y, z)
    core = np.maximum(ZP - z, inner)
    slab_free = (ZP + S_SEAL) - z
    sec = {f: np.minimum(*[sector(x, y, a) for a, ff in LEGS if ff == f])
           for f in ("hot", "cold")}
    hot = np.maximum(np.maximum(T_WALL / 2 - G, core), np.minimum(slab_free, sec["hot"]))
    cold = np.maximum(np.maximum(G + T_WALL / 2, core), np.minimum(slab_free, sec["cold"]))

    for a, f in LEGS:
        lx, ly = leg_xy(a)
        plen = np.maximum(np.maximum(SHELL - z, z - (ZP + 0.5)),
                          np.maximum(inner, sector(x, y, a)))
        bore = np.maximum(np.hypot(x - lx, y - ly) - BORE_R, z - (SHELL + 1.0))
        v = np.minimum(plen, bore)
        if f == "hot":
            hot = np.minimum(hot, v)
        else:
            cold = np.minimum(cold, v)
    hot, cold = np.maximum(hot, outer), np.maximum(cold, outer)
    solid = np.maximum(outer, -np.minimum(hot, cold))
    return outer, solid, hot, cold


def mesh(field, step, origin):
    v, f, _, _ = marching_cubes(field.astype(np.float32), 0.0, spacing=(step,) * 3)
    return v + origin, f


def shade(v, f, rgb):
    tri = v[f]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
    L = np.array([-0.4, -0.6, 0.7]); L /= np.linalg.norm(L)
    s = 0.35 + 0.65 * np.abs(n @ L)
    return tri, np.c_[np.outer(s, rgb), np.ones(len(s))]


def draw3d(ax, parts, view, title):
    for tri, col in parts:
        ax.add_collection3d(Poly3DCollection(tri, facecolors=col, edgecolor="none"))
    ax.set_xlim(-75, 75); ax.set_ylim(-75, 75); ax.set_zlim(-42, 92)
    ax.set_box_aspect((150, 150, 134))
    ax.view_init(*view); ax.set_axis_off(); ax.set_title(title, fontsize=11)


def main():
    step = 0.8
    xs = np.arange(-73, 73.01, step); zs = np.arange(-43, 93.01, step)
    X, Y, Z = np.meshgrid(xs, xs, zs, indexing="ij")
    outer, solid, hot, cold = fields(X, Y, Z)
    o = np.array([xs[0], xs[0], zs[0]])
    grey = (0.62, 0.66, 0.72)

    # 1) exterior
    v, f = mesh(outer, step, o); ext = shade(v, f, grey)
    # 2) cutaway: remove the quadrant x>0, y<0 (through the 315 deg cold leg)
    cut = np.maximum(solid, -np.maximum(-X, Y))
    v, f = mesh(cut, step, o); cw = shade(v, f, grey)
    hot_cut = np.maximum(hot, -np.maximum(-X, Y))
    del X, Y, Z

    fig = plt.figure(figsize=(18, 7), dpi=110)
    draw3d(fig.add_subplot(131, projection="3d"), [ext], (18, -60), "Exterior")
    draw3d(fig.add_subplot(132, projection="3d"), [ext], (-35, -60),
           "Underside: four arched legs = four ports")
    draw3d(fig.add_subplot(133, projection="3d"), [cw], (22, -45),
           "Cut-away: gyroid core, quadrant headers, leg bores")
    plt.tight_layout(); plt.savefig("figures/cyl_hx_3d.png"); plt.close()

    # close-up of one leg: exterior (outer field, so no hidden internal
    # faces) and a section through the leg axis
    lx, ly = leg_xy(45)
    fs = 0.15
    gx = np.arange(lx - 26, lx + 26, fs); gy = np.arange(ly - 26, ly + 26, fs)
    gz = np.arange(-41.5, 2, fs)
    Xl, Yl, Zl = np.meshgrid(gx, gy, gz, indexing="ij")
    ol, _, _, _ = fields(Xl, Yl, Zl)
    ol = np.maximum(ol, Zl - 0.5)                  # leg only, below the body
    ol = np.maximum(ol, np.hypot(Xl - lx, Yl - ly) - 25)
    del Xl, Yl, Zl
    v, f = mesh(ol, fs, np.array([gx[0], gy[0], gz[0]]))
    tri, col = shade(v, f, grey)
    fig = plt.figure(figsize=(13, 6), dpi=110)
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax.add_collection3d(Poly3DCollection(tri, facecolors=col, edgecolor="none"))
    ax.set_xlim(lx - 26, lx + 26); ax.set_ylim(ly - 26, ly + 26); ax.set_zlim(-42, 2)
    ax.set_box_aspect((52, 52, 44)); ax.view_init(-12, -60); ax.set_axis_off()
    ax.set_title("Leg: slow taper from the body to a 3/8\"-18 NPT male tip", fontsize=11)
    ax = fig.add_subplot(1, 2, 2)
    q = np.arange(-28, 28, 0.08); qz = np.arange(-42, 18, 0.08)
    Q, QZ = np.meshgrid(q, qz, indexing="ij")
    tt = np.radians(45)
    _, s2, h2, c2 = fields(lx + Q * np.cos(tt), ly + Q * np.sin(tt), QZ)
    m = np.zeros(s2.shape, int); m[s2 < 0] = 1; m[h2 < 0] = 2; m[c2 < 0] = 3
    ax.imshow(m.T, origin="lower", extent=(q[0], q[-1], qz[0], qz[-1]),
              cmap=ListedColormap(["white", "#555b66", "#E4572E", "#2E86AB"]),
              vmin=0, vmax=3, interpolation="nearest")
    ax.set_title("Section through the leg axis (hot leg, 45 deg)", fontsize=11)
    ax.set_xlabel("mm"); ax.set_ylabel("z (mm)")
    plt.tight_layout(); plt.savefig("figures/cyl_hx_leg.png"); plt.close()

    # 3) colour sections: solid grey, hot red, cold blue
    cmap = ListedColormap(["white", "#555b66", "#E4572E", "#2E86AB"])
    def lab(s, h, c):
        m = np.zeros(s.shape, int)
        m[s < 0] = 1; m[h < 0] = 2; m[c < 0] = 3
        return m
    g = np.arange(-75, 75.01, 0.25)
    fig, axs = plt.subplots(1, 4, figsize=(20, 5.6), dpi=110)
    for ax, zc, t in [(axs[0], 8.0, "Header level z = 8 mm\n(4 quadrant plenums)"),
                      (axs[1], 45.0, "Core z = 45 mm\n(interleaved gyroid channels)")]:
        Xs, Ys = np.meshgrid(g, g, indexing="ij")
        _, s, h, c = fields(Xs, Ys, np.full_like(Xs, zc))
        ax.imshow(lab(s, h, c).T, origin="lower", extent=(g[0], g[-1], g[0], g[-1]),
                  cmap=cmap, vmin=0, vmax=3, interpolation="nearest")
        ax.set_title(t); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
    gz = np.arange(-43, 93.01, 0.25)
    for ax, ang, t in [(axs[2], 45, "Vertical section through HOT legs\n(45 / 225 deg)"),
                       (axs[3], 135, "Vertical section through COLD legs\n(135 / 315 deg)")]:
        Sg, Zg = np.meshgrid(g, gz, indexing="ij")
        tt = np.radians(ang)
        _, s, h, c = fields(Sg * np.cos(tt), Sg * np.sin(tt), Zg)
        ax.imshow(lab(s, h, c).T, origin="lower", extent=(g[0], g[-1], gz[0], gz[-1]),
                  cmap=cmap, vmin=0, vmax=3, interpolation="nearest")
        ax.set_title(t); ax.set_xlabel("distance along section (mm)"); ax.set_ylabel("z (mm)")
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color="#555b66", label="solid wall"),
                        Patch(color="#E4572E", label="hot fluid"),
                        Patch(color="#2E86AB", label="cold fluid")],
               loc="lower center", ncol=3, fontsize=11)
    plt.tight_layout(rect=(0, 0.06, 1, 1)); plt.savefig("figures/cyl_hx_sections.png")
    print("saved")


if __name__ == "__main__":
    main()
