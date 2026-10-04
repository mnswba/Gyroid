"""
radial_disc_hx.py
------------------------------------------------------------------------
Radial counter-flow GYROID disc heat exchanger ("pancake").

This script does NOT make an STL. It:
  1. holds every design parameter in one place,
  2. checks the design (sizes, flow areas, hot/cold separation),
  3. writes  ntop_recipe.md : a step-by-step, block-by-block process for
     building the model NATIVELY in nTop with these exact numbers,
  4. renders preview.png (colour sections) so you can see the result
     before you build it.

Run in VS Code:   pip install numpy matplotlib
                  python radial_disc_hx.py

Flow (z = disc axis, origin at disc centre, mm):
  HOT : in  at the TOP centre port   -> hub (upper half) -> flows OUTWARD
        through the gyroid           -> rim ring (upper) -> out tangential
        port on the rim at +x.
  COLD: in  tangential rim port at -x -> rim ring (lower) -> flows INWARD
        through the gyroid            -> hub (lower half) -> out at the
        BOTTOM centre port.
  -> radial counter-flow. A divider plate at z = 0 splits the hub and the
     rim ring; thin seal rings at the core edges let each fluid's gyroid
     labyrinth open only into its own half.
------------------------------------------------------------------------
"""
import numpy as np

# ======================================================================
# 1. PARAMETERS (mm)
# ======================================================================
P = dict(
    R_OUT=62.0,        # disc outer radius
    H=36.0,            # disc height (z = -H/2 .. +H/2)
    R_EDGE=6.0,        # rounding of the disc edges
    SHELL=2.0,         # outer wall thickness
    R_HUB=10.0,        # hub (centre header) radius
    R_CORE=50.0,       # outer radius of the gyroid core
    R_RING=57.0,       # outer radius of the rim ring header
    DIV=2.0,           # divider plate thickness (z = 0, hub and ring)
    SEAL=1.5,          # radial thickness of the seal rings at core edges
    CELL=6.0,          # gyroid cell size (x = y = z)
    T_WALL=0.6,        # gyroid wall thickness
    # ports: 3/8"-18 NPT male (OD 17.145 mm at the large end)
    PORT_OD=17.145,
    PORT_BORE=10.0,
    PORT_LEN=14.0,     # length beyond the body
    COLLAR_D=21.0,     # boss diameter where the axial ports meet the disc
    Z_RIM_PORT=9.0,    # |z| of the rim port axes (hot +z, cold -z)
)


def derived(p=P):
    d = {}
    d["k"] = 2 * np.pi / p["CELL"]
    d["z_top"] = p["H"] / 2
    d["z_bore"] = p["H"] / 2 - p["SHELL"] - 1     # axial bores start here
    d["half_in"] = p["H"] / 2 - p["SHELL"] - p["DIV"] / 2     # inner height per half
    d["x_rim_port"] = (p["R_CORE"] + p["R_RING"]) / 2         # rim port axis x
    d["rim_port_exit_y"] = np.sqrt(p["R_OUT"] ** 2 - d["x_rim_port"] ** 2)
    # flow area each fluid sees (~50 % of the open gyroid section)
    open_frac = 1 - p["T_WALL"] * 3.1 / p["CELL"]              # ~ void share
    h_core = p["H"] - 2 * p["SHELL"]
    for name, r in (("hub", p["R_HUB"]), ("rim", p["R_CORE"])):
        d[f"A_{name}_cm2"] = 2 * np.pi * r * h_core * open_frac / 2 / 100
    d["A_port_cm2"] = np.pi * (p["PORT_BORE"] / 2) ** 2 / 100
    d["v_ratio_hub_to_rim"] = p["R_CORE"] / p["R_HUB"]
    # gyroid level that gives ~T_WALL: |F| <= g0, with |grad F| ~ 1.2 k
    d["g0"] = 0.5 * p["T_WALL"] * d["k"] * 1.2
    d["core_volume_cm3"] = np.pi * (p["R_CORE"] ** 2 - p["R_HUB"] ** 2) * h_core / 1e3
    d["area_density_m2_m3"] = 3.1 / (p["CELL"] / 1000)      # ~3.1/a for a gyroid
    d["HT_area_cm2"] = d["area_density_m2_m3"] * d["core_volume_cm3"] / 100
    return d


def checks(p=P, d=None):
    d = d or derived(p)
    msgs = []
    top_of_rim_port = p["Z_RIM_PORT"] + p["PORT_OD"] / 2
    if top_of_rim_port > p["H"] / 2:
        msgs.append(f"rim port ({top_of_rim_port:.1f}) sticks out above the disc "
                    f"top ({p['H'] / 2:.1f}): raise H or lower Z_RIM_PORT")
    if p["Z_RIM_PORT"] - p["PORT_BORE"] / 2 < p["DIV"] / 2:
        msgs.append("rim port bore cuts the divider plate")
    if p["Z_RIM_PORT"] + p["PORT_BORE"] / 2 > p["H"] / 2 - p["SHELL"]:
        msgs.append("rim port bore cuts the top/bottom wall")
    if p["R_HUB"] < p["PORT_BORE"] / 2 + 1:
        msgs.append("hub smaller than the axial port bore")
    if p["R_RING"] > p["R_OUT"] - p["SHELL"] - p["R_EDGE"] / 2:
        msgs.append("rim ring too close to the outer wall / edge rounding")
    if p["T_WALL"] < 0.4:
        msgs.append("gyroid wall below ~0.4 mm: check your printer minimum")
    return msgs


# ======================================================================
# 2. IMPLICIT MODEL (for the preview and the leak check only)
#    value < 0 = inside. Mirrors the nTop recipe step by step.
# ======================================================================
def _smin(a, b, k):
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0, 1)
    return b * (1 - h) + a * h - k * h * (1 - h)


def _rcyl(r, z, R, z0, z1, rr):
    hh, zc = (z1 - z0) / 2, z - (z0 + z1) / 2
    qx, qz = r - (R - rr), np.abs(zc) - (hh - rr)
    return (np.hypot(np.maximum(qx, 0), np.maximum(qz, 0))
            + np.minimum(np.maximum(qx, qz), 0) - rr)


def _rod(x, y, z, axis, c, r, t0, t1):
    d = (x - c[0], y - c[1], z - c[2])
    t = d[0] * axis[0] + d[1] * axis[1] + d[2] * axis[2]
    rad = np.sqrt(np.maximum(d[0] ** 2 + d[1] ** 2 + d[2] ** 2 - t ** 2, 0))
    return np.maximum(rad - r, np.maximum(t0 - t, t - t1))


def model(x, y, z, p=P):
    d = derived(p)
    r = np.hypot(x, y)
    H2 = p["H"] / 2
    xr = d["x_rim_port"]
    zr = p["Z_RIM_PORT"]
    ro = p["PORT_OD"] / 2
    rb = p["PORT_BORE"] / 2

    # Step 1-2: disc + port bosses
    outer = _rcyl(r, z, p["R_OUT"], -H2, H2, p["R_EDGE"])
    for s in (1, -1):
        outer = _smin(outer, _rod(x, y, z, (0, 0, s), (0, 0, 0), ro, 0,
                                  H2 + p["PORT_LEN"]), 4)
        outer = _smin(outer, _rod(x, y, z, (0, 0, s), (0, 0, 0), p["COLLAR_D"] / 2,
                                  0, H2 + 1.5), 3)
    outer = _smin(outer, _rod(x, y, z, (0, 1, 0), (xr, 0, zr), ro, -2,
                              d["rim_port_exit_y"] + p["PORT_LEN"]), 3)
    outer = _smin(outer, _rod(x, y, z, (0, -1, 0), (-xr, 0, -zr), ro, -2,
                              d["rim_port_exit_y"] + p["PORT_LEN"]), 3)
    inner = outer + p["SHELL"]

    # Step 3: gyroid field and core
    k = d["k"]
    F = (np.sin(k * x) * np.cos(k * y) + np.sin(k * y) * np.cos(k * z)
         + np.sin(k * z) * np.cos(k * x))
    g0 = d["g0"]
    core = np.maximum(np.maximum(p["R_HUB"] - r, r - p["R_CORE"]), inner)
    seal = np.minimum(r - p["R_HUB"] - p["SEAL"], p["R_CORE"] - p["SEAL"] - r)
    upper, lower = p["DIV"] / 2 - z, p["DIV"] / 2 + z           # < 0 in that half
    hot = np.maximum(np.maximum((g0 - F) / k, core), np.minimum(-seal, upper))
    cold = np.maximum(np.maximum((F + g0) / k, core), np.minimum(-seal, lower))

    # Step 4: hub, rim ring and port bores
    hub = np.maximum(r - (p["R_HUB"] + 0.5), inner)
    ring = np.maximum(np.maximum(p["R_CORE"] - 0.5 - r, r - p["R_RING"]), inner)
    # rim bores are wider than the ring: keep them outside R_RING - 1 so they
    # open only into their own ring and never reach the seal / core
    rim_clip = (p["R_RING"] - 1.0) - r
    for v in (np.maximum(hub, upper), np.maximum(ring, upper),
              _rod(x, y, z, (0, 0, 1), (0, 0, 0), rb, H2 - p["SHELL"] - 1,
                   H2 + p["PORT_LEN"] + 1),
              np.maximum(_rod(x, y, z, (0, 1, 0), (xr, 0, zr), rb, -2,
                              d["rim_port_exit_y"] + p["PORT_LEN"] + 1),
                         rim_clip)):
        hot = np.minimum(hot, v)
    for v in (np.maximum(hub, lower), np.maximum(ring, lower),
              _rod(x, y, z, (0, 0, -1), (0, 0, 0), rb, H2 - p["SHELL"] - 1,
                   H2 + p["PORT_LEN"] + 1),
              np.maximum(_rod(x, y, z, (0, -1, 0), (-xr, 0, -zr), rb, -2,
                              d["rim_port_exit_y"] + p["PORT_LEN"] + 1),
                         rim_clip)):
        cold = np.minimum(cold, v)
    hot, cold = np.maximum(hot, outer), np.maximum(cold, outer)
    solid = np.maximum(outer, -np.minimum(hot, cold))
    return solid, hot, cold


def preview(path="preview.png", res=0.2):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    cmap = ListedColormap(["white", "#555b66", "#E4572E", "#2E86AB"])
    H2, L = P["H"] / 2, P["R_OUT"] + P["PORT_LEN"] + 6
    g = np.arange(-L, L, res)
    gz = np.arange(-H2 - P["PORT_LEN"] - 3, H2 + P["PORT_LEN"] + 3, res)

    def lab(s, h, c):
        m = np.zeros(s.shape, int)
        m[s < 0], m[h < 0], m[c < 0] = 1, 2, 3
        assert not np.any((h < 0) & (c < 0)), "hot and cold overlap!"
        # touching without a wall also counts as a leak: no cell may be
        # within half a pixel of both fluids
        assert not np.any((h < res / 2) & (c < res / 2)), "hot and cold touch!"
        return m

    fig = plt.figure(figsize=(18, 10), dpi=100)
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.25])
    X, Z = np.meshgrid(g, gz, indexing="ij")
    for i, (yc, t) in enumerate([(0.0, "Section y = 0 (hub, divider, rim ring)")]):
        ax = fig.add_subplot(gs[0, :])
        ax.imshow(lab(*model(X, np.full_like(X, yc), Z)).T, origin="lower",
                  extent=(g[0], g[-1], gz[0], gz[-1]), cmap=cmap, vmin=0, vmax=3,
                  interpolation="nearest")
        ax.set_title(t); ax.set_xlabel("x (mm)"); ax.set_ylabel("z (mm)")
    Xp, Yp = np.meshgrid(g, g, indexing="ij")
    for j, (zc, t) in enumerate([(P["Z_RIM_PORT"], "z = +%.0f: HOT half (flows outward)"),
                                 (-P["Z_RIM_PORT"], "z = %.0f: COLD half (flows inward)"),
                                 (0.0, "z = 0: divider plate")]):
        ax = fig.add_subplot(gs[1, j])
        ax.imshow(lab(*model(Xp, Yp, np.full_like(Xp, zc))).T, origin="lower",
                  extent=(g[0], g[-1], g[0], g[-1]), cmap=cmap, vmin=0, vmax=3,
                  interpolation="nearest")
        ax.set_title(t % zc if "%" in t else t); ax.set_xlabel("x (mm)")
    fig.legend(handles=[Patch(color="#555b66", label="solid"),
                        Patch(color="#E4572E", label="hot"),
                        Patch(color="#2E86AB", label="cold")],
               loc="upper right", ncol=3)
    plt.tight_layout()
    plt.savefig(path)
    print(f"preview -> {path}")


# ======================================================================
# 3. nTop RECIPE (plain English, filled with the current numbers)
# ======================================================================
RECIPE = """# Radial disc gyroid heat exchanger: native nTop build recipe

Generated by `radial_disc_hx.py`. All numbers are in mm. If you change a
parameter in the script, re-run it and this file updates.

Block names are given as they appear in nTop's block search. If your nTop
version names a block differently, search for the operation (for example
"union" or "cylinder"). Create every number below as a **Scalar variable
(Length)** first. To drive the design from nTop Automate later, right-click
a variable and choose **Make Input**.

## 0. Variables
| Name | Value | Name | Value |
|---|---|---|---|
| R_OUT | {R_OUT} | H | {H} |
| R_EDGE | {R_EDGE} | SHELL | {SHELL} |
| R_HUB | {R_HUB} | R_CORE | {R_CORE} |
| R_RING | {R_RING} | DIV | {DIV} |
| SEAL | {SEAL} | CELL | {CELL} |
| T_WALL | {T_WALL} | PORT_OD | {PORT_OD} |
| PORT_BORE | {PORT_BORE} | PORT_LEN | {PORT_LEN} |
| COLLAR_D | {COLLAR_D} | Z_RIM_PORT | {Z_RIM_PORT} |
| X_RIM_PORT | {x_rim_port:.2f} | Y_RIM_EXIT | {rim_port_exit_y:.2f} |

The disc axis is Z and the origin is the disc centre. The disc spans
z = -{z_top:.1f} to +{z_top:.1f}.

## 1. Outer body
1. **Cylinder** (implicit). Centre (0, 0, -{z_top:.1f}), axis +Z,
   radius R_OUT, height H. Name it `Disc`.
2. Round its edges: **Fillet Body** (or smooth the edge with an offset
   pair: **Offset Body** −{R_EDGE} then **Offset Body** +{R_EDGE}).
   Keep the result named `Disc`.

## 2. Port bosses (3/8"-18 NPT male, OD {PORT_OD})
3. **Cylinder** `Boss_HotIn`: centre (0, 0, +{z_top:.1f}), axis +Z,
   radius PORT_OD/2, length PORT_LEN.
4. **Cylinder** `Boss_ColdOut`: centre (0, 0, -{z_top:.1f}), axis -Z,
   same size.
5. **Cylinder** `Boss_HotOut` (tangential rim port): centre
   (X_RIM_PORT, -2, +Z_RIM_PORT), axis +Y, radius PORT_OD/2, length
   Y_RIM_EXIT + PORT_LEN + 2. It starts inside the rim ring so the
   bore opens into it.
6. **Cylinder** `Boss_ColdIn`: centre (-X_RIM_PORT, +2, -Z_RIM_PORT),
   axis -Y, same size.
7. Two collars, **Cylinder** radius COLLAR_D/2, length 3, at the top and
   bottom centre.
8. **Boolean Union** of `Disc` with all bosses and collars, smoothing
   (blend) radius 3-4. This gives `Outer`.

## 3. Inner volume
9. **Offset Body** `Outer` by -SHELL. This gives `Inner`, the region
   that can hold fluid.

## 4. Gyroid field (the two labyrinths)
10. Build the scalar field
    F = sin(kx)·cos(ky) + sin(ky)·cos(kz) + sin(kz)·cos(kx),
    with k = 2π / CELL = {k:.4f} 1/mm. Use either:
    - nTop's TPMS gyroid field (Periodic Lattice / TPMS field with unit
      cell Gyroid, cell size CELL), or
    - the math blocks (X, Y, Z coordinate fields → Multiply → Sin, Cos
      → Multiply → Add).
11. Wall level g0 = {g0:.4f}, which gives a wall of about {T_WALL} mm.
    Make two bodies:
    - `HotSide` = **Body from Field** (F > +g0), i.e. field (g0 − F).
    - `ColdSide` = **Body from Field** (F < −g0), i.e. field (F + g0).

    Alternative: make the wall directly with **Walled TPMS** (Gyroid,
    thickness T_WALL, cell CELL). If you do, use the same field for both
    sides so the phase matches.

## 5. Core region and seals
12. `CoreRing` = **Boolean Intersect** of `Inner` with an annulus
    (**Cylinder** R_CORE minus **Cylinder** R_HUB, full height).
13. `UpperHalf` = half-space z > DIV/2. `LowerHalf` = half-space
    z < −DIV/2. Make each as a large box or a **Plane** body.
14. `SealRings` = annulus R_HUB … R_HUB+SEAL **Boolean Union**
    annulus R_CORE−SEAL … R_CORE, full height.
15. `HotCore` = (`CoreRing` ∩ `HotSide`) − (`SealRings` − `UpperHalf`).
    Inside the seal rings, hot is kept only in the upper half.
16. `ColdCore` = (`CoreRing` ∩ `ColdSide`) − (`SealRings` − `LowerHalf`).

## 6. Headers and port bores
17. `Hub` = **Cylinder** radius R_HUB+0.5, full height, ∩ `Inner`.
18. `Ring` = annulus R_CORE−0.5 … R_RING, full height, ∩ `Inner`.
19. `HotVoid` = `HotCore` ∪ (`Hub` ∩ `UpperHalf`) ∪ (`Ring` ∩ `UpperHalf`)
    ∪ bore (**Cylinder** radius PORT_BORE/2 through `Boss_HotIn`,
    from z = +{z_bore:.1f} up to the tip; it must NOT reach the divider
    plate) ∪ (bore through `Boss_HotOut` along +Y **minus** a cylinder of
    radius R_RING−1). The rim bore (Ø{PORT_BORE}) is wider than the ring,
    so this clip stops it cutting into the seal ring and the core. Without
    it, hot and cold would touch.
20. `ColdVoid` = `ColdCore` ∪ (`Hub` ∩ `LowerHalf`) ∪ (`Ring` ∩ `LowerHalf`)
    ∪ bore through `Boss_ColdOut` (from z = -{z_bore:.1f} down to the tip)
    ∪ (bore through `Boss_ColdIn` minus the same radius-(R_RING−1) cylinder).

## 7. Final part
21. `HX` = **Boolean Subtract**: `Outer` − (`HotVoid` ∪ `ColdVoid`).
22. Check it with **Section** views at y = 0 and z = ±Z_RIM_PORT. They
    should match `preview.png`.
23. For CFD, `HotVoid` and `ColdVoid` (each ∩ `Outer`) are the two fluid
    domains.
24. To print, run **Mesh from Implicit Body** (or your print-prep blocks)
    on `HX`.

## 8. Threads
nTop has no simple NPT-thread block that I can confirm. Two options:
- Leave the bosses as plain Ø{PORT_OD} and cut the 3/8"-18 NPT thread
  after printing (die), or
- Import a thread body from CAD and **Boolean Union** it onto each boss.

## Design numbers (from the script)
- Core volume: {core_volume_cm3:.0f} cm³. Estimated wall area: about
  {HT_area_cm2:.0f} cm² (gyroid area density about 3.1/CELL; estimate).
- Flow area per fluid: about {A_hub_cm2:.1f} cm² at the hub and
  {A_rim_cm2:.1f} cm² at the rim. The port bore is {A_port_cm2:.2f} cm².
- Velocity falls about {v_ratio_hub_to_rim:.0f}× from hub to rim (flow
  area grows with radius). Expect higher heat-transfer coefficients near
  the hub. A graded cell, finer toward the rim, can even this out.
"""


def write_recipe(path="ntop_recipe.md"):
    d = derived()
    with open(path, "w") as f:
        f.write(RECIPE.format(**P, **d))
    print(f"recipe  -> {path}")


if __name__ == "__main__":
    d = derived()
    print("Radial disc HX")
    for k_ in ("core_volume_cm3", "HT_area_cm2", "A_hub_cm2", "A_rim_cm2",
               "A_port_cm2", "g0", "x_rim_port"):
        print(f"  {k_:16s} {d[k_]:.3f}")
    for m in checks(P, d) or ["all geometry checks passed"]:
        print("  check:", m)
    write_recipe()
    preview()
