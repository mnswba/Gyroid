#!/usr/bin/env python3
"""Geometry toolkit for a graded, stretched, offset sheet-gyroid two-fluid heat exchanger core.

The core is a sheet (double-network) gyroid: the wall is the region
    |G(X, Y, Z) - c| <= t,   G = sin X cos Y + sin Y cos Z + sin Z cos X
which splits space into two interpenetrating, non-intersecting channels:
    side A (hot stream):  G - c >  t
    side B (cold stream): G - c < -t

Three design knobs are graded along the flow axis x:
    s(x)  streamwise stretch = streamwise cell length / transverse cell length
    t(x)  level-set half-thickness (set so the printed wall never drops below a minimum)
    c(x)  level-set offset (shifts fluid volume between side A and side B)

Grading uses a phase-continuous streamwise coordinate
    X(x) = integral_0^x 2*pi / (s(x') * L_t) dx'
so the local cell length is exactly s(x) * L_t. The naive form X = 2*pi*x / L(x) does not do
this: its local period is L / (1 - x L'/L), which distorts cells far from x = 0.

Everything here is geometry only (volume fractions, wetted area, hydraulic diameter, wall
thickness). It does not predict heat transfer or pressure drop; those need CFD or tests.

Usage:
    python gyroid_hx_geometry.py                 # print tables and write figures
    python gyroid_hx_geometry.py --stl core.stl  # also export a watertight STL test coupon
"""

import argparse
import os
import struct

import numpy as np
from skimage.measure import marching_cubes, mesh_surface_area

TWO_PI = 2.0 * np.pi
FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")

# Colours from the validated categorical palette (slot 2 orange, slot 1 blue) plus a neutral wall.
COLOR_HOT = "#eb6834"
COLOR_COLD = "#2a78d6"
COLOR_WALL = "#52514e"


def gyroid(X, Y, Z):
    return np.sin(X) * np.cos(Y) + np.sin(Y) * np.cos(Z) + np.sin(Z) * np.cos(X)


def gyroid_grad(X, Y, Z):
    """Gradient of G with respect to the phase coordinates (X, Y, Z)."""
    gx = np.cos(X) * np.cos(Y) - np.sin(Z) * np.sin(X)
    gy = -np.sin(X) * np.sin(Y) + np.cos(Y) * np.cos(Z)
    gz = -np.sin(Y) * np.sin(Z) + np.cos(Z) * np.cos(X)
    return gx, gy, gz


def _phase_grid(n, endpoint):
    if endpoint:
        return np.linspace(0.0, TWO_PI, n + 1)
    return np.linspace(0.0, TWO_PI, n, endpoint=False) + np.pi / n


def volume_fractions(t, c, n=96):
    """Fluid/wall volume fractions of one periodic cell. Independent of stretch (affine map)."""
    ph = _phase_grid(n, endpoint=False)
    G = gyroid(*np.meshgrid(ph, ph, ph, indexing="ij"))
    eps_a = np.mean(G - c > t)
    eps_b = np.mean(G - c < -t)
    return eps_a, eps_b, 1.0 - eps_a - eps_b


def unit_cell_metrics(t, c, stretch, n=96):
    """Geometry of one cell with transverse size 1 and streamwise size `stretch`.

    Lengths are returned in units of the transverse cell size L_t; multiply by L_t to get mm.
    """
    eps_a, eps_b, wall = volume_fractions(t, c, n)

    ph = _phase_grid(n, endpoint=True)
    X, Y, Z = np.meshgrid(ph, ph, ph, indexing="ij")
    G = gyroid(X, Y, Z)
    spacing = (stretch / n, 1.0 / n, 1.0 / n)
    cell_volume = stretch

    def area(level):
        verts, faces, _, _ = marching_cubes(G, level=level, spacing=spacing)
        return mesh_surface_area(verts, faces), verts

    area_a, _ = area(c + t)
    area_b, _ = area(c - t)
    area_mid, verts_mid = area(c)

    # Local wall thickness on the mid-surface, first order in t: delta = 2 t / |grad_phys G|.
    # Phase -> physical scaling: dX/dx = 2 pi / stretch, dY/dy = dZ/dz = 2 pi.
    Xm = verts_mid[:, 0] * TWO_PI / stretch
    Ym = verts_mid[:, 1] * TWO_PI
    Zm = verts_mid[:, 2] * TWO_PI
    gx, gy, gz = gyroid_grad(Xm, Ym, Zm)
    grad_phys = np.sqrt((gx * TWO_PI / stretch) ** 2 + (gy * TWO_PI) ** 2 + (gz * TWO_PI) ** 2)
    delta_local = 2.0 * t / grad_phys

    a_a = area_a / cell_volume  # wetted area per unit core volume, side A
    a_b = area_b / cell_volume
    return {
        "eps_a": eps_a,
        "eps_b": eps_b,
        "wall": wall,
        "a_a": a_a,
        "a_b": a_b,
        "a_mid": area_mid / cell_volume,
        "dh_a": 4.0 * eps_a / a_a,
        "dh_b": 4.0 * eps_b / a_b,
        "delta_mean": wall / (area_mid / cell_volume),
        "delta_p05": float(np.percentile(delta_local, 5)),
    }


def _bisect(f, lo, hi, iters=40):
    flo = f(lo)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if (fmid > 0) == (flo > 0):
            lo, flo = mid, fmid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def offset_for_split(t, split_a, n=96):
    """Level-set offset c giving side A the fraction `split_a` of all fluid volume."""
    def f(c):
        eps_a, eps_b, _ = volume_fractions(t, c, n)
        return eps_a / (eps_a + eps_b) - split_a
    return _bisect(f, -0.6, 0.6)


def thickness_for_min_wall(delta_min_rel, c, stretch, n=96):
    """Half-thickness t whose 5th-percentile local wall thickness equals delta_min_rel (units of L_t).

    delta_p05 is linear in t to first order, so one evaluation at a probe t suffices.
    """
    t_probe = 0.1
    m = unit_cell_metrics(t_probe, c, stretch, n)
    return t_probe * delta_min_rel / m["delta_p05"]


def solve_cell(stretch, delta_min_mm, L_t_mm, split_a, n=96):
    """Find (t, c) for one station: wall >= delta_min everywhere (5th pct), fluid split = split_a."""
    c, t = 0.0, 0.1
    for _ in range(3):
        t = thickness_for_min_wall(delta_min_mm / L_t_mm, c, stretch, n)
        c = offset_for_split(t, split_a, n)
    return t, c, unit_cell_metrics(t, c, stretch, n)


def stretch_profile(xi, s_hot, s_cold):
    """Smooth monotone stretch from the hot end (xi = 0) to the cold end (xi = 1)."""
    w = 0.5 - 0.5 * np.cos(np.pi * np.clip(xi, 0.0, 1.0))
    return s_hot + (s_cold - s_hot) * w


def streamwise_phase(x, cell_len):
    """Phase-continuous X(x) = integral 2 pi / cell_len dx (trapezoid rule)."""
    k = TWO_PI / cell_len
    return np.concatenate([[0.0], np.cumsum(0.5 * (k[1:] + k[:-1]) * np.diff(x))])


def graded_field(x, y, z, L_t, s_of_x, t_of_x, c_of_x):
    """Signed wall field W = t - |G - c| (W >= 0 is solid) and G - c on an (x, y, z) grid."""
    X = streamwise_phase(x, s_of_x * L_t)
    Xg, Yg, Zg = np.meshgrid(X, TWO_PI * y / L_t, TWO_PI * z / L_t, indexing="ij")
    Gc = gyroid(Xg, Yg, Zg) - c_of_x[:, None, None]
    return t_of_x[:, None, None] - np.abs(Gc), Gc


def write_binary_stl(path, verts, faces, weld_tol=1e-5):
    # Marching cubes emits a few sliver triangles whose vertices coincide once stored as float32.
    # Weld vertices on a weld_tol grid and drop the collapsed triangles so the mesh stays manifold.
    key = np.round(verts / weld_tol).astype(np.int64)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    verts = uniq * weld_tol
    faces = inv.reshape(-1)[faces]
    keep = (faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2]) & (faces[:, 0] != faces[:, 2])
    faces = faces[keep]
    tri = verts[faces]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    norm = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = np.divide(normals, norm, out=np.zeros_like(normals), where=norm > 0)
    with open(path, "wb") as fh:
        fh.write(b"graded sheet-gyroid heat exchanger coupon".ljust(80, b" "))
        fh.write(struct.pack("<I", len(faces)))
        rec = np.zeros(len(faces), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
        rec["n"] = normals
        rec["v"] = tri
        rec.tofile(fh)
    return len(faces)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Lt", type=float, default=4.0, help="transverse cell size, mm (default 4)")
    ap.add_argument("--length", type=float, default=100.0, help="core length along flow, mm (default 100)")
    ap.add_argument("--s-hot", type=float, default=2.0, help="streamwise stretch at the hot end (default 2.0)")
    ap.add_argument("--s-cold", type=float, default=1.2, help="streamwise stretch at the cold end (default 1.2)")
    ap.add_argument("--wall-min", type=float, default=0.30, help="minimum printed wall, mm (default 0.30)")
    ap.add_argument("--split-a", type=float, default=0.60, help="side-A share of fluid volume (default 0.60)")
    ap.add_argument("--stations", type=int, default=6, help="stations in the printed schedule (default 6)")
    ap.add_argument("--n", type=int, default=96, help="voxels per cell edge for metrics (default 96)")
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--stl", help="write a watertight STL coupon of the graded wall to this path")
    ap.add_argument("--stl-res", type=float, default=0.08, help="STL voxel size, mm (default 0.08)")
    args = ap.parse_args()

    # 1. Sanity checks against known values.
    m0 = unit_cell_metrics(0.0, 0.0, 1.0, args.n)
    print("Check: G = 0 surface area per cell volume (unit cube) = %.4f  (literature ~3.09)" % m0["a_mid"])
    print("Check: side split at t = 0, c = 0 = %.4f / %.4f  (symmetry: 0.5 / 0.5)\n" % (m0["eps_a"], m0["eps_b"]))

    # 2. Effect of streamwise stretch at a fixed minimum wall and a 50/50 split.
    print("Stretch sweep: L_t = %.1f mm, min wall %.2f mm, 50/50 split" % (args.Lt, args.wall_min))
    print("  s     t      porosity  a [1/m]   Dh [mm]  wall mean [mm]")
    sweep = []
    for s in (1.0, 1.25, 1.5, 2.0, 2.5, 3.0):
        t, c, m = solve_cell(s, args.wall_min, args.Lt, 0.5, args.n)
        a_si = m["a_a"] / (args.Lt * 1e-3)
        sweep.append((s, t, m["eps_a"] + m["eps_b"], a_si, m["dh_a"] * args.Lt, m["delta_mean"] * args.Lt))
        print("  %.2f  %.4f  %.3f     %7.0f   %.3f    %.3f" % sweep[-1])

    # 3. Graded schedule along the core (hot end at x = 0).
    print("\nGraded schedule: length %.0f mm, s %.2f (hot end) -> %.2f (cold end), side-A share %.2f"
          % (args.length, args.s_hot, args.s_cold, args.split_a))
    print("  x [mm]  s     t       c       eps_A  eps_B  wall   a_A [1/m]  Dh_A [mm]  Dh_B [mm]")
    xi_st = np.linspace(0.0, 1.0, args.stations)
    sched = []
    for xi in xi_st:
        s = float(stretch_profile(xi, args.s_hot, args.s_cold))
        t, c, m = solve_cell(s, args.wall_min, args.Lt, args.split_a, args.n)
        sched.append((xi * args.length, s, t, c, m))
        print("  %6.1f  %.3f  %.4f  %+.4f  %.3f  %.3f  %.3f  %7.0f    %.3f      %.3f"
              % (xi * args.length, s, t, c, m["eps_a"], m["eps_b"], m["wall"],
                 m["a_a"] / (args.Lt * 1e-3), m["dh_a"] * args.Lt, m["dh_b"] * args.Lt))

    # 4. Why the phase-continuous coordinate matters.
    L0, L1 = args.s_hot * args.Lt, args.s_cold * args.Lt
    dL = (L1 - L0) / args.length
    naive_period_end = L1 / (1.0 - args.length * dL / L1)
    print("\nNaive X = 2 pi x / L(x) with linear L: local period at x = L is %.2f mm, intended %.2f mm (%+.0f%%)."
          % (naive_period_end, L1, 100.0 * (naive_period_end / L1 - 1.0)))
    print("The phase-continuous X(x) used here gives exactly the intended period at every x.")

    if not args.no_figures:
        make_figures(args, sweep, sched)

    if args.stl:
        export_stl(args, xi_st, sched)


def _schedule_interp(args, sched, x):
    xs = np.array([row[0] for row in sched])
    s = stretch_profile(x / args.length, args.s_hot, args.s_cold)
    t = np.interp(x, xs, [row[2] for row in sched])
    c = np.interp(x, xs, [row[3] for row in sched])
    return s, t, c


def make_figures(args, sweep, sched):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    os.makedirs(FIG_DIR, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#8a8984", "axes.labelcolor": "#0b0b0b",
                         "xtick.color": "#52514e", "ytick.color": "#52514e"})

    # Figure 1: stretch sweep, one measure per panel (no dual axes).
    sw = np.array(sweep)
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.3), constrained_layout=True)
    panels = [(3, "Wetted area per side [1/m]"), (4, "Hydraulic diameter [mm]"), (2, "Porosity [-]")]
    for ax, (col, title) in zip(axs, panels):
        ax.plot(sw[:, 0], sw[:, col], color=COLOR_COLD, lw=2, marker="o", ms=5)
        ax.set_title(title, loc="left", fontsize=10)
        ax.set_xlabel("Streamwise stretch s = L_x / L_t")
        ax.grid(axis="y", color="#e6e5e1", lw=0.8)
    axs[2].set_ylim(0, 1)
    fig.suptitle("Sheet gyroid, L_t = %.0f mm, min wall %.2f mm, 50/50 split" % (args.Lt, args.wall_min),
                 x=0.01, ha="left", fontsize=11)
    fig.savefig(os.path.join(FIG_DIR, "stretch_sweep.png"), dpi=150)
    plt.close(fig)

    # Figure 2: mid-plane (z = 0) section of the graded core vs. a naive grading.
    res = 0.05
    x = np.arange(0.0, args.length + res, res)
    y = np.arange(0.0, 3 * args.Lt + res, res)
    z = np.array([0.0])
    s_x, t_x, c_x = _schedule_interp(args, sched, x)
    W, Gc = graded_field(x, y, z, args.Lt, s_x, t_x, c_x)
    phase = np.where(W[:, :, 0] >= 0, 1, np.where(Gc[:, :, 0] > 0, 0, 2))

    L_x = s_x * args.Lt
    Xn = TWO_PI * x / L_x
    Xg, Yg = np.meshgrid(Xn, TWO_PI * y / args.Lt, indexing="ij")
    Gn = gyroid(Xg, Yg, np.zeros_like(Xg)) - c_x[:, None]
    phase_naive = np.where(t_x[:, None] - np.abs(Gn) >= 0, 1, np.where(Gn > 0, 0, 2))

    cmap = ListedColormap([COLOR_HOT, COLOR_WALL, COLOR_COLD])
    fig, axs = plt.subplots(2, 1, figsize=(11, 3.6), constrained_layout=True, sharex=True)
    for ax, ph, title in ((axs[0], phase, "Phase-continuous grading (this design)"),
                          (axs[1], phase_naive, "Naive grading X = 2πx / L(x): cells distort toward the cold end")):
        ax.imshow(ph.T, origin="lower", cmap=cmap, vmin=0, vmax=2, interpolation="nearest",
                  extent=(x[0], x[-1], y[0], y[-1]), aspect="auto")
        ax.set_title(title, loc="left", fontsize=10)
        ax.set_ylabel("y [mm]")
    axs[1].set_xlabel("x [mm]  (hot end at x = 0, stretch %.1f → %.1f)" % (args.s_hot, args.s_cold))
    handles = [plt.Rectangle((0, 0), 1, 1, color=col) for col in (COLOR_HOT, COLOR_WALL, COLOR_COLD)]
    fig.legend(handles, ["Side A (hot stream)", "Wall", "Side B (cold stream)"], loc="outside upper right",
               ncol=3, frameon=False, fontsize=9)
    fig.savefig(os.path.join(FIG_DIR, "graded_core_section.png"), dpi=150)
    plt.close(fig)
    print("\nFigures written to %s" % FIG_DIR)


def export_stl(args, xi_st, sched):
    res = args.stl_res
    length = min(args.length, 30.0)
    x = np.arange(0.0, length + res, res)
    y = np.arange(0.0, 2 * args.Lt + res, res)
    z = np.arange(0.0, 2 * args.Lt + res, res)
    s_x, t_x, c_x = _schedule_interp(args, sched, x)
    W, _ = graded_field(x, y, z, args.Lt, s_x, t_x, c_x)
    W = np.pad(W, 1, constant_values=-1.0)  # close the solid at the coupon boundary
    verts, faces, _, _ = marching_cubes(W, level=0.0, spacing=(res, res, res))
    verts -= res
    n_tri = write_binary_stl(args.stl, verts, faces)
    print("STL coupon (%.0f x %.0f x %.0f mm, %d triangles) written to %s"
          % (length, 2 * args.Lt, 2 * args.Lt, n_tri, args.stl))


if __name__ == "__main__":
    main()
