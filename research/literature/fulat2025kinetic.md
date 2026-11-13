# Fulat, Moore, Al-Awashra, Bohdan, Tsirou, Pohl & Zweibel (2025): Kinetic simulations of particle acceleration at shock waves in plasma with pre-existing turbulence

- Ref: Proc. 39th Int. Cosmic Ray Conf. (ICRC2025), PoS(ICRC2025)042 (Crossref); DOI 10.22323/1.501.0042;
  conference proceedings, not peer-reviewed in the journal sense
- Read: 2026-11-13, full (4 pages; the relativistic part skimmed); topics: shock-alfvenic-turbulence

## Claim
In 2D3V PIC, a nonrelativistic high-Mach oblique shock (theta_Bn = 60 deg) with pre-existing compressive
turbulence (dn/n0 about 15 %) has stronger foreshock whistlers and more non-thermal electrons than the same shock
in a homogeneous upstream. (Preliminary version of the oblique-shock results of Fulat et al. 2026.)

## Method
- Turbulent slabs simulated separately and inserted into the shock run with field/current matching (Fulat et al.
  2023 framework).
- Table 1: 20 particles per cell, M_s = 36, M_A = 32, theta_Bn = 60 deg, beta = 1, T = 1e7 K, m_i/m_e = 100,
  lambda_se = 20 cells, Omega_i^-1 = 4.8e4 steps; runs H (homogeneous) and T (turbulent); end about 20 Omega_i^-1.

## Numbers to compare against
- Whistler (right-hand) growth rate from |B_R|^2 ~ exp(-2 gamma x / v0): gamma_H = (6.13 +- 0.25)e-3 Omega_e,
  gamma_T = (4.36 +- 0.43)e-3 Omega_e (k = 2 uncertainties); growth starts at x - x_sh about 90 lambda_si.
- Whistler level at x - x_sh about 20 lambda_si: |B_R|^2/B0^2 about 0.04 (H) and 0.1 (T).
- Nonlinear foreshock structures appear at about 9 (H) and 7 (T) Omega_i^-1, size up to 7 lambda_si (T) vs
  about 2 lambda_si (H).
- Downstream electrons (20-90 lambda_si behind, t about 20 Omega_i^-1): fraction with gamma - 1 >= 3 is 0.14 % (H)
  vs 0.23 % (T); energy share 1.9 % vs 3.6 %; max Lorentz factor 40 % higher in T; tail index about -2.1 (T)
  (Fig. 2).

## Assessment
- Compressive, kinetic-scale turbulence (not Alfvenic, not L_c >= 10 d_i) and a 20 Omega_i^-1 run: the effect is
  attributed to stronger whistlers, i.e. a microphysical channel, not field-line meandering.
- One run per case; no seed or ensemble spread given for the 0.14 vs 0.23 % difference.
- Superseded for the oblique case by the refereed/preprint Fulat et al. (2026); kept as the dated first report.
