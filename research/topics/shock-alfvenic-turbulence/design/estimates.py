"""Design estimates for the topic (answers to the review of 2026-11-13, Majors 2, 3 and 9).

Run: .venv/bin/python research/topics/shock-alfvenic-turbulence/design/estimates.py > .../design/estimates.md

Units: lengths in d_i, times in 1/Omega_ci, speeds in v_A, fields in B0. Shock frame inflow speed U1 = M_A v_A.
Geometry: x shock normal, B0 along y (in the plane, theta_Bn = 90 deg), slab fluctuations delta B_x(y) along B0.

Assumptions (stated, not derived here):
- Shock thickness L_sh := ion foot width 0.68 M_A d_i (reflected-ion turnaround of a perpendicular shock; the
  standard estimate attributed to Woods 1971 / Livesey et al. 1982, not read for this topic).
- Spectrum of the slab: P(k_n) ~ 1 / (1 + (k_n L_c)^(5/3)) on k_n = 2 pi n / L_y, n = 1..N (guo2010effect).
- Thermal speeds: v_t = sqrt(T/m) (one-dimensional rms). beta_s = 2 mu0 n T_s / B0^2, so v_ti = sqrt(beta_i / 2) v_A,
  v_te = sqrt(beta_e m_i / (2 m_e)) v_A.
- Throughput 1.1e10 particle-steps/s on 8 GPUs (blast_wave/full, notebook 2026-11-13).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import erfc, wofz


# ---------------------------------------------------------------- dispersion of parallel waves
def zeta_fn(z):
    """Plasma dispersion function Z(z) = i sqrt(pi) w(z)."""
    return 1j * math.sqrt(math.pi) * wofz(z)


def disp_parallel(omega, k, sign, mr, beta_e, beta_i, wpe_wce):
    """D(omega, k) for parallel circularly polarized waves in a Maxwellian e-i plasma (nonrelativistic).
    sign=+1: Stix R (electron-resonant, whistler branch), -1: L (ion-resonant). Units Omega_ci, d_i, v_A."""
    c_va = wpe_wce * math.sqrt(mr)  # c / v_A = omega_pi / Omega_ci
    tot = k * k * c_va**2 - omega**2
    # species: (omega_p^2 / omega_pi^2, Omega_s / Omega_ci signed, sqrt(2 T/m) / v_A)
    for wp2, om, vt in ((1.0, 1.0, math.sqrt(beta_i)), (mr, -mr, math.sqrt(beta_e * mr))):
        z0 = omega / (k * vt)
        zs = (omega + sign * om) / (k * vt)
        tot -= c_va**2 * wp2 * z0 * zeta_fn(zs)
    return tot


def solve(k, sign, mr, beta_e, beta_i, wpe_wce):
    """Complex root by secant iteration from the cold Hall-MHD guess."""
    h = 0.5 * k
    w0 = k * (math.sqrt(1 + h * h) + sign * h) + 0j
    w1 = w0 * (1 + 1e-3) - 1e-4j
    f0 = disp_parallel(w0, k, sign, mr, beta_e, beta_i, wpe_wce)
    f1 = disp_parallel(w1, k, sign, mr, beta_e, beta_i, wpe_wce)
    for _ in range(200):
        if f1 == f0:
            break
        w2 = w1 - f1 * (w1 - w0) / (f1 - f0)
        w0, f0 = w1, f1
        w1, f1 = w2, disp_parallel(w2, k, sign, mr, beta_e, beta_i, wpe_wce)
        if abs(w1 - w0) < 1e-12 * max(1.0, abs(w1)):
            break
    return w1


def modes(L_y, L_c, n_max):
    n = np.arange(1, n_max + 1)
    k = 2 * math.pi * n / L_y
    p = 1.0 / (1.0 + (k * L_c) ** (5.0 / 3.0))
    return k, p / p.sum()  # weights of <delta B^2>


def polarization_fraction(t, k, w, rot_rate):
    """Energy fraction still in the x-polarization after time t (linear start along x, each mode rotating its
    plane at rot_rate(k))."""
    return float(np.sum(w * np.cos(rot_rate * t) ** 2))


def time_to(frac, k, w, rot_rate, t_max=2000.0):
    ts = np.linspace(0, t_max, 200001)
    for t in ts:
        if polarization_fraction(t, k, w, rot_rate) < frac:
            return t
    return float("inf")


# ---------------------------------------------------------------- field-line excursion and access
def excursion(db, L_y, L_c, n_max, seeds=200, ny=4096):
    """rms and median peak-to-peak displacement xi_x(y) = int delta B_x / B0 dy for random phases;
    db = rms of delta B_x / B0."""
    k, w = modes(L_y, L_c, n_max)
    a = np.sqrt(2.0 * w) * db  # amplitudes: sum a^2 / 2 = db^2
    xi_rms = math.sqrt(float(np.sum(a * a / (2 * k * k))))
    rng = np.random.default_rng(1)
    y = np.linspace(0, L_y, ny, endpoint=False)
    pp = []
    for _ in range(seeds):
        ph = rng.uniform(0, 2 * math.pi, size=k.size)
        xi = np.sum((a / k)[:, None] * np.sin(k[:, None] * y[None, :] + ph[:, None]), axis=0)
        pp.append(xi.max() - xi.min())
    return xi_rms, float(np.median(pp))


def subluminal_fraction(db, M_A, mr, beta_e, nsamp=200001):
    """Mean over the front (delta B_x Gaussian, rms db) of the fraction of thermal upstream electrons with
    v_par cos(theta_loc) > U1, i.e. able to stream upstream along the local field (HTF-subluminal access)."""
    rng = np.random.default_rng(2)
    bx = rng.normal(0.0, db, nsamp)
    cos_t = np.abs(bx) / np.sqrt(1.0 + bx * bx)
    vte = math.sqrt(beta_e * mr / 2.0)
    with np.errstate(divide="ignore"):
        thr = np.where(cos_t > 0, M_A / cos_t, np.inf)
    return float(np.mean(0.5 * erfc(thr / (math.sqrt(2.0) * vte)))), float(np.median(cos_t))


# ---------------------------------------------------------------- cost
def cost(mr, wpe_wce, L_x, L_y, cells_per_de, ppc, t_run, rate=1.1e10, cfl=0.5):
    di_cells = math.sqrt(mr) * cells_per_de
    cells = (L_x * di_cells) * (L_y * di_cells)
    particles = cells * ppc
    dt_wpe = cfl / cells_per_de  # dt in 1/omega_pe with dx = d_e / cells_per_de and c = 1
    steps = t_run * mr * wpe_wce / dt_wpe
    return particles, steps, particles * steps / rate / 3600.0


def main():
    beta_e = beta_i = 1.0
    wpe_wce = 10.0
    M_A = 4.5
    db = 0.6
    print("# Design estimates (generated by `design/estimates.py`)\n")
    print(f"Fixed unless stated: beta_e = beta_i = {beta_e}, omega_pe/Omega_ce = {wpe_wce}, M_A = {M_A}, "
          f"rms delta B_x / B0 = {db}, theta_Bn = 90 deg with B0 in the plane. "
          "Assumptions: see the module docstring.\n")

    # 1. Hall / kinetic rotation of a linear polarization
    print("## 1. Rotation of a linear polarization (review, Major 2)\n")
    print("Rotation rate of the plane of a linear polarization = (Re omega_R - Re omega_L) / 2. "
          "Cold Hall MHD: k^2 d_i v_A / 2. Kinetic: roots of the parallel Maxwellian dispersion relation "
          "(R and L), m_i/m_e = 100.\n")
    print("| k d_i | Hall rate [rad Ω_ci] | kinetic rate [rad Ω_ci] | Im ω_L / Ω_ci | Im ω_R / Ω_ci |")
    print("|---|---|---|---|---|")
    for kk in (0.05, 0.1, 0.2, 0.3, 0.5, 0.63):
        wr = solve(kk, +1, 100, beta_e, beta_i, wpe_wce)
        wl = solve(kk, -1, 100, beta_e, beta_i, wpe_wce)
        print(f"| {kk:.2f} | {kk * kk / 2:.4f} | {(wr.real - wl.real) / 2:.4f} | {wl.imag:.2e} | {wr.imag:.2e} |")
    print()
    print("Energy-weighted depolarization of a linear δB_x slab spectrum (kinetic rates, m_i/m_e = 100): "
          "time until the x-polarized share of δB² falls below 0.75 and 0.5.\n")
    print("| L_y [d_i] | L_c [d_i] | t(0.75) [Ω_ci⁻¹] | t(0.5) [Ω_ci⁻¹] |")
    print("|---|---|---|---|")
    for L_y, L_c in ((50, 10), (50, 20), (100, 20), (200, 40), (400, 80)):
        n_max = int(L_y / 1.0)  # modes down to k d_i = 2 pi
        k, w = modes(L_y, L_c, n_max)
        sel = k < 2.0
        k, w = k[sel], w[sel] / w[sel].sum()
        rate = np.array([(solve(kk, +1, 100, beta_e, beta_i, wpe_wce).real
                          - solve(kk, -1, 100, beta_e, beta_i, wpe_wce).real) / 2 for kk in k])
        print(f"| {L_y} | {L_c} | {time_to(0.75, k, w, rate):.1f} | {time_to(0.5, k, w, rate):.1f} |")
    print()

    # 2. Xi
    print("## 2. Field-line excursion Ξ = ξ_x / L_sh (review, Major 3)\n")
    L_sh = 0.68 * M_A
    print(f"L_sh = 0.68 M_A d_i = {L_sh:.2f} d_i at M_A = {M_A}. ξ_pp: median peak-to-peak over 200 random-phase "
          "realizations. t_conn = ξ_pp / U1: time during which one field line meets the front more than once. "
          "t_along = (L_y / 2) / v_te: thermal transit between two intersections (v_te at m_i/m_e = 100 and 1836).\n")
    print("| L_y [d_i] | L_c [d_i] | δB_x/B0 | ξ_rms [d_i] | ξ_pp [d_i] | Ξ = ξ_pp/L_sh | t_conn [Ω_ci⁻¹] | "
          "t_along 100 | t_along 1836 |")
    print("|---|---|---|---|---|---|---|---|---|")
    rows = [(50, 10, db), (50, 20, db), (100, 25, db), (200, 50, db), (400, 100, db), (200, 50, 0.3), (200, 50, 1.0)]
    for L_y, L_c, d in rows:
        xr, xpp = excursion(d, L_y, L_c, n_max=int(L_y))
        ta = [(L_y / 2) / math.sqrt(beta_e * mr / 2) for mr in (100, 1836)]
        print(f"| {L_y} | {L_c} | {d} | {xr:.1f} | {xpp:.1f} | {xpp / L_sh:.1f} | {xpp / M_A:.1f} | "
              f"{ta[0]:.1f} | {ta[1]:.1f} |")
    # Guo & Giacalone 2010 for reference: delta B^2/B0^2 = 0.3 in two components (circular), L_c = L_z = 1024
    xr, xpp = excursion(math.sqrt(0.15), 1024, 1024, n_max=1024)
    print(f"| 1024 (guo2010effect) | 1024 | {math.sqrt(0.15):.2f} | {xr:.1f} | {xpp:.1f} | "
          f"{xpp / (0.68 * 5.6):.1f} (M_A 5.6) | {xpp / 5.6:.1f} | – | – |")
    print()

    # 3. Access
    print("## 3. Access of thermal electrons along the local field (review, Major 3)\n")
    print("Mean over the front of the fraction of thermal upstream electrons with v_∥ cos θ_loc > U1 "
          "(δB_x Gaussian with the given rms; δB_z ignored). Median cos θ_loc in brackets.\n")
    print("| δB_x/B0 | M_A | m_i/m_e = 25 | 100 | 400 | 1836 |")
    print("|---|---|---|---|---|---|")
    for d in (0.3, 0.6, 1.0):
        for M in (3.0, 4.5, 8.0):
            vals = []
            for mr in (25, 100, 400, 1836):
                f, med = subluminal_fraction(d, M, mr, beta_e)
                vals.append(f"{f:.2e}")
            print(f"| {d} (cos θ ≈ {med:.2f}) | {M} | " + " | ".join(vals) + " |")
    print()

    # 4. thresholds
    print("## 4. Thresholds at simulated and real mass ratios (review, Major 3)\n")
    print("| Threshold | Formula (source) | m_i/m_e = 100 | 400 | 1836 |")
    print("|---|---|---|---|---|")
    rows = [
        ("whistler critical, θ = 90 / 70 / 45°", "0.5 √(m_i/m_e) cos θ (lindberg2024electron form)",
         lambda mr: " / ".join(f"{0.5 * math.sqrt(mr) * math.cos(math.radians(t)):.1f}" for t in (90, 70, 45))),
        ("whistler injection, θ = 70°, β_e = 1", "(cos θ/2) √(m_i/m_e) √β_e (amano2010critical Eq. 2)",
         lambda mr: f"{0.5 * math.cos(math.radians(70)) * math.sqrt(mr):.1f}"),
        ("surfing trapping (α = 0.2)", "(1+α)(m_i/m_e)^(2/3) (bohdan2017electron Eq. 2)",
         lambda mr: f"{1.2 * mr ** (2 / 3):.0f}"),
        ("SSDA injection, θ = 70°", "M_A ≳ 55 cos θ (amano2022theory Eq. 14, factors 1)",
         lambda mr: f"{55 * math.cos(math.radians(70)):.0f} (no m_i/m_e)"),
        ("ion-scale ripples", "M_A ≈ 3.5 (trotta2019electron, hybrid)", lambda mr: "3.5"),
    ]
    for name, form, fn in rows:
        print(f"| {name} | {form} | {fn(100)} | {fn(400)} | {fn(1836)} |")
    print()

    # 5. cost
    print("## 5. Cost of one 2-D run (review, Major 9 and matrix)\n")
    print("Piston box preloaded over its full length L_x; 20 particles per cell in total; CFL 0.5 dx/c; "
          "1.1e10 particle-steps/s. Debye check: λ_De / dx with v_te from β_e.\n")
    print("| m_i/m_e | cells per d_e | L_x × L_y [d_i] | t_run [Ω_ci⁻¹] | particles | steps | hours | λ_De/dx |")
    print("|---|---|---|---|---|---|---|---|")
    for mr, cpd, Lx, Ly, T in ((100, 5, 200, 100, 40), (100, 5, 200, 200, 40), (100, 5, 200, 400, 40),
                               (100, 10, 200, 100, 40), (400, 5, 200, 100, 40), (400, 5, 100, 50, 20)):
        p, s, h = cost(mr, wpe_wce, Lx, Ly, cpd, 20, T)
        vte_c = math.sqrt(beta_e * mr / 2) / (wpe_wce * math.sqrt(mr))  # v_te / c
        print(f"| {mr} | {cpd} | {Lx} × {Ly} | {T} | {p:.2e} | {s:.2e} | {h:.0f} | {vte_c * cpd:.2f} |")


if __name__ == "__main__":
    main()
