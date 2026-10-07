"""Micro-level (unit-cell) analysis of gyroid structures for thermal management.

Everything here is computed from the level-set gyroid

    f(x, y, z) = sin X cos Y + sin Y cos Z + sin Z cos X,   X = 2*pi*x/a

using only NumPy, so the numbers in docs/micro-level-thermal-management.md can be
reproduced and checked.

Sub-commands
------------
geometry   porosity, surface-area density, hydraulic diameter and wall thickness
           of sheet and network gyroids versus the level-set parameter
keff       effective (homogenised) conductivity of the sheet-gyroid solid frame
           by a periodic finite-volume conduction solve
coldplate  illustrative porous-fin estimate of a sheet-gyroid cold plate at several
           cell sizes; Nu and f*Re are ASSUMED inputs, not gyroid correlations
all        run everything

Usage: python3 scripts/gyroid_micro.py [geometry|keff|coldplate|all]
"""

import sys

import numpy as np

TAU = 2.0 * np.pi
F_MAX = 1.5  # max |f| of the level-set gyroid


def gyroid(X, Y, Z):
    return np.sin(X) * np.cos(Y) + np.sin(Y) * np.cos(Z) + np.sin(Z) * np.cos(X)


def gyroid_grad_norm(X, Y, Z):
    """|grad_X f| with X the angular coordinate (period 2*pi)."""
    gx = np.cos(X) * np.cos(Y) - np.sin(Z) * np.sin(X)
    gy = np.cos(Y) * np.cos(Z) - np.sin(X) * np.sin(Y)
    gz = np.cos(Z) * np.cos(X) - np.sin(Y) * np.sin(Z)
    return np.sqrt(gx * gx + gy * gy + gz * gz)


def level_stats(levels, N=256, eps=0.04):
    """For each level c return (P(f < c), S(c)).

    P(f < c) is the volume fraction below the level; S(c) is the area of the
    iso-surface f = c per unit volume, made dimensionless with the cell size a
    (physical area density = S / a).  Area uses the smoothed-delta form of the
    co-area formula: A/V = < delta_eps(f - c) |grad_x f| >, |grad_x f| = 2*pi |grad_X f| / a.
    """
    levels = np.asarray(levels, dtype=float)
    s = (np.arange(N) + 0.5) / N * TAU
    X, Y = np.meshgrid(s, s, indexing="ij")
    below = np.zeros_like(levels)
    area = np.zeros_like(levels)
    for z in s:
        Z = np.full_like(X, z)
        f = gyroid(X, Y, Z).ravel()
        g = gyroid_grad_norm(X, Y, Z).ravel()
        for i, c in enumerate(levels):
            below[i] += np.count_nonzero(f < c)
            d = f - c
            m = np.abs(d) < eps
            area[i] += np.sum((1.0 + np.cos(np.pi * d[m] / eps)) / (2.0 * eps) * g[m])
    n = float(N) ** 3
    return below / n, TAU * area / n


def sheet_table(N=256):
    """Sheet gyroid: solid where |f| <= t, two congruent fluid labyrinths."""
    ts = np.array([0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80])
    below, area = level_stats(np.concatenate([[0.0], ts, -ts]), N=N)
    s_mid = area[0]
    k = len(ts)
    phi_s = below[1:1 + k] - below[1 + k:]
    s_face = area[1:1 + k]  # wetted area density of ONE labyrinth (f = +t surface)
    rows = []
    for t, ph, sf in zip(ts, phi_s, s_face):
        eps1 = (1.0 - ph) / 2.0          # fluid fraction of one labyrinth
        dh = 4.0 * eps1 / sf              # hydraulic diameter of one labyrinth / a
        wall = ph / s_mid                 # mean wall thickness / a (solid volume / mid-surface area)
        rows.append((t, ph, 1.0 - ph, 2.0 * sf, dh, wall))
    return s_mid, rows


def network_table(N=256):
    """Network (skeletal) gyroid: solid where f > c, one fluid labyrinth f < c."""
    cs = np.array([-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9])
    below, area = level_stats(cs, N=N)
    rows = []
    for c, b, s in zip(cs, below, area):
        phi_s = 1.0 - b
        dh = 4.0 * b / s
        rows.append((c, phi_s, b, s, dh))
    return rows


def solve_sheet_t(phi_target, N=192):
    """Bisection for the sheet half-level t giving solid fraction phi_target."""
    lo, hi = 0.0, F_MAX
    for _ in range(30):
        t = 0.5 * (lo + hi)
        b, _ = level_stats([t, -t], N=N, eps=0.05)
        if b[0] - b[1] < phi_target:
            lo = t
        else:
            hi = t
    return 0.5 * (lo + hi)


def sheet_props(phi_s, N=192):
    """(t, S_face, S_mid, Dh/a, wall/a) for a sheet gyroid of solid fraction phi_s."""
    t = solve_sheet_t(phi_s, N=N)
    _, area = level_stats([t, 0.0], N=N, eps=0.05)
    s_face, s_mid = area
    dh = 4.0 * ((1.0 - phi_s) / 2.0) / s_face
    return t, s_face, s_mid, dh, phi_s / s_mid


def keff_sheet(t, N=96, contrast=1e-4, tol=1e-7, maxit=20000):
    """Homogenised conductivity k_eff/k_s of a sheet gyroid (|f| <= t solid).

    Periodic cell, unit mean gradient imposed along z, face conductances are
    harmonic means of voxel conductivities, Jacobi-preconditioned CG.
    Fluid voxels get conductivity `contrast` * k_s.  By cubic symmetry the
    effective conductivity tensor is isotropic, so one direction suffices.
    """
    s = (np.arange(N) + 0.5) / N * TAU
    X, Y, Z = np.meshgrid(s, s, s, indexing="ij")
    k = np.where(np.abs(gyroid(X, Y, Z)) <= t, 1.0, contrast)
    phi = float(np.mean(k == 1.0))
    del X, Y, Z
    kf = [2.0 * k * np.roll(k, -1, axis=ax) / (k + np.roll(k, -1, axis=ax)) for ax in range(3)]
    diag = sum(kf[ax] + np.roll(kf[ax], 1, axis=ax) for ax in range(3))

    def apply(T):
        out = diag * T
        for ax in range(3):
            out -= kf[ax] * np.roll(T, -1, axis=ax) + np.roll(kf[ax] * T, 1, axis=ax)
        return out

    # grid spacing h = 1/N, mean gradient G = 1 along z (axis 2): flux balance gives
    # A T' = h * (k_{i+1/2} - k_{i-1/2}) for the periodic fluctuation T'
    b = (kf[2] - np.roll(kf[2], 1, axis=2)) / N
    T = np.zeros_like(k)
    r = b - apply(T)
    zv = r / diag
    p = zv.copy()
    rz = np.vdot(r, zv)
    bn = np.linalg.norm(b)
    for it in range(maxit):
        Ap = apply(p)
        alpha = rz / np.vdot(p, Ap)
        T += alpha * p
        r -= alpha * Ap
        if np.linalg.norm(r) < tol * bn:
            break
        zv = r / diag
        rz_new = np.vdot(r, zv)
        p = zv + (rz_new / rz) * p
        rz = rz_new
    # mean z-flux: q = k_face * ((T[k+1]-T[k])*N + 1)
    q = kf[2] * ((np.roll(T, -1, axis=2) - T) * N + 1.0)
    return phi, float(np.mean(q)), it + 1


# ---------------------------------------------------------------- reporting

def report_geometry():
    s_mid, rows = sheet_table()
    print("SHEET GYROID  (solid |f| <= t; two separate fluid labyrinths)")
    print(f"  mid-surface (f = 0) area density S0 = {s_mid:.4f} / a")
    print("  t     phi_s   porosity  S_wet,total*a  Dh/a (each)  wall/a")
    for t, ph, por, st, dh, w in rows:
        print(f"  {t:4.2f}  {ph:6.3f}  {por:7.3f}   {st:10.3f}     {dh:8.3f}     {w:6.3f}")
    print()
    print("NETWORK GYROID  (solid f > c; one fluid labyrinth f < c)")
    print("  c      phi_s   porosity  S*a     Dh/a")
    for c, ph, por, s, dh in network_table():
        print(f"  {c:5.2f}  {ph:6.3f}  {por:7.3f}  {s:6.3f}  {dh:6.3f}")
    print()


def report_keff():
    print("SHEET GYROID effective conductivity (fluid/solid conductivity ratio 1e-4)")
    print("  t     N    phi_s   k_eff/k_s   k_eff/(phi_s k_s)   CG iters")
    for t in (0.2, 0.4, 0.6):
        for N in (64, 96, 128):
            phi, ke, its = keff_sheet(t, N=N)
            print(f"  {t:3.1f}  {N:4d}  {phi:6.3f}   {ke:7.4f}       {ke / phi:6.3f}           {its}")
    print()


# Water near 30 C (rounded handbook values)
WATER = dict(k=0.615, rho=996.0, cp=4178.0, mu=7.97e-4)


def report_coldplate(phi_s=0.30, k_s=380.0, H=2.0e-3, q_flux=1.0e6,
                     u=0.5, L=10.0e-3, keff_factor=None, nus=(5.0, 10.0, 20.0), f_re=100.0):
    """Porous-fin (volume-averaged) model of a base-heated sheet-gyroid layer.

    Both labyrinths carry the same coolant.  Fluid is taken at inlet temperature
    (no caloric rise) and solid conduction is lumped into k_eff; this is a
    first-order sizing estimate, not a substitute for conjugate CFD.
    """
    t, s_face, s_mid, dh_a, wall_a = sheet_props(phi_s)
    if keff_factor is None:
        _, ke, _ = keff_sheet(t, N=96)
        keff_factor = ke  # k_eff / k_s
    k_eff = keff_factor * k_s
    print(f"SHEET-GYROID COLD PLATE, phi_s={phi_s}, t={t:.3f}, k_s={k_s} W/mK, "
          f"k_eff={k_eff:.1f} W/mK, H={H * 1e3:.1f} mm, q''={q_flux / 1e4:.0f} W/cm^2")
    print(f"  water, interstitial velocity u={u} m/s, flow length L={L * 1e3:.0f} mm, ASSUMED f*Re={f_re}")
    print("  a[mm]  wall[um]  Dh[um]  Sv[m2/m3]  Re   |  Nu  h[W/m2K]  mH    eta   R''[K cm2/W]  dT[K]  |  dp[kPa]")
    for a_mm in (2.0, 1.0, 0.5, 0.25):
        a = a_mm * 1e-3
        dh = dh_a * a
        sv = 2.0 * s_face / a           # both faces of the sheet are wetted
        re = WATER["rho"] * u * dh / WATER["mu"]
        dp = f_re * WATER["mu"] * u / (2.0 * dh ** 2) * L
        for j, nu in enumerate(nus):
            h = nu * WATER["k"] / dh
            m = np.sqrt(h * sv / k_eff)
            eta = np.tanh(m * H) / (m * H)
            r = 1.0 / (k_eff * m * np.tanh(m * H))
            lead = (f"  {a_mm:4.2f}   {wall_a * a * 1e6:6.0f}   {dh * 1e6:6.0f}  {sv:8.0f}  {re:5.0f} |"
                    if j == 0 else " " * 45 + "|")
            tail = f"  |  {dp / 1e3:6.2f}" if j == 0 else ""
            print(f"{lead} {nu:4.0f}  {h:8.0f}  {m * H:5.2f}  {eta:5.2f}   {r * 1e4:8.4f}   {q_flux * r:6.1f}{tail}")
    print()


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd in ("geometry", "all"):
        report_geometry()
    if cmd in ("keff", "all"):
        report_keff()
    if cmd in ("coldplate", "all"):
        report_coldplate()
