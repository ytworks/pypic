# Jikei, Grošelj & Sironi (2026): Dependence of particle acceleration efficiency on shock velocity in weakly magnetized electron-ion shocks

- Ref: journal reference not checked (the arXiv metadata gives none, and no DOI); arXiv:2607.05778v1,
  7 July 2026. Preprint only, as far as the metadata says.
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
At fixed upstream magnetization, quasi-parallel shocks change from Bell-dominated at low speed to
Weibel-dominated at high speed, with the transition near M_A ~ 100. Ions are accelerated with similar
efficiency in both regimes; electrons are not: ~15 % of the shock energy goes to non-thermal electrons
in the Weibel regime, less than ~2 % in the Bell regime.

## Method
- PIC, OSIRIS, 2-D (x-y), downstream (wall) frame, B0 in the simulation plane at theta_B = 20 deg.
- sigma = B0^2/[4 pi n0 (m_i + m_e) c^2] = 10^-4.5 = 3.16e-5 (upstream frame, Eq. 1); upstream
  four-velocity u0/c = 1/6, 1/3, 2/3, 1, 4/3; m_i/m_e = 100. M_A = v0/V_A with V_A/c ~ sqrt(sigma)
  (Sect. 3.4); values per run are not tabulated. Upstream temperature, M_s and beta: not stated.
- Width 163.84 d_i, length not stated; dx = 0.4 d_e, c dt = 0.5 dx, 16 particles per cell per species,
  cubic splines, 4 filter passes; run to omega_pi t = 10000 (not given in units of 1/Omega_ci).
- Checks: narrow boxes (10.24 d_i) with 64 particles per cell for u0/c = 1/3 and 2/3; the u0/c = 4/3
  run continued to omega_pi t = 25000 in a narrow box (App. B). 1-D periodic SMILEI runs of the Bell
  instability driven by an electron-ion beam at sigma = 1e-4 (App. A).

## Numbers to compare against
- Low-current Bell regime needs eta = J_cr/(n0 e c) < eta_crit = 2 sqrt(sigma) (Eqs. 5-6), i.e. a
  returning-ion density ratio alpha < 2/M_A (Eq. 10). Bell-dominated runs settle at
  eta/eta_crit <~ 0.5 and alpha << 0.1 (after omega_pi t ~ 7000 for u0/c = 1/6, ~9000 for 1/3);
  Weibel-dominated runs keep alpha >~ 0.1 (Fig. 2, Sect. 3.2).
- Transition at u0/c ~ 2/3, M_A ~ 100; M_A = 100 (v/0.1 c) (sigma/1e-6)^(-1/2) (Eq. 13).
- Non-thermal electrons = p > 3 p_peak, p_peak from the peak of (gamma - 1) f_e(gamma); energy in units
  of the incoming ion energy, 50-100 d_i behind the shock, omega_pi t = 7000-10000 (Eq. 7, Fig. 3c):
  ~0.15 (~0.3 of the electron energy) for u0/c = 1 and 4/3; ~0.015 (~0.06 of the electron energy)
  for u0/c = 1/3. Electrons hold 0.25-0.45 of the downstream energy.
- Electron thermal peak at p/(m_i u0) ~ 0.07 (Bell runs) and ~0.2 (Weibel runs) (Fig. 3b).
- E_max,i ~ t in Bell-dominated shocks, ~ t^(1/2) in Weibel-dominated ones (Fig. 4).
- No number fraction, spectral index or K_ep is given.

## Assessment
- Shown by simulation: the regime change, the self-regulation of the returning-ion current and the
  drop of the electron non-thermal energy in the Bell regime, at m_i/m_e = 100 and v0 >~ 0.16 c.
- Stated without a separate test: Bell fields of large amplitude and wavelength above the thermal
  electron gyroradius make the shock "effectively perpendicular for electrons" and suppress their
  injection (Sect. 3.3, taken from their earlier paper, ApJ 998, 149).
- Possible inconsistency: Sect. 3.4 puts the electron tail of the u0/c <~ 1/3 runs "about five orders
  of magnitude below the peak", next to an energy fraction of 0.015 for u0/c = 1/3 (not checked
  against the figure). The u0/c = 1/6 electron tail may not be numerically converged (Sect. 3.3);
  the 2/3 run is still evolving at the end.
- Stated limits and open points (Sects. 3.4, 4.4): whether M_A ~ 100 holds for non-relativistic
  shocks needs runs with v0 << 0.1 c and sigma << 1e-6 (Weibel growth may be weaker there); the
  SNR-like regime M_A >~ 100, v0 < 0.1 c is not covered; multi-D PIC of SNR-like shocks has mostly
  been quasi-perpendicular and 1-D parallel runs cannot hold Weibel modes, so a multi-D quasi-parallel
  survey is missing; the electron efficiency may have a minimum near v0 ~ 0.1 c.
- Disagreement: Sect. 4.4 reads gupta2024electron as "efficiency is higher for slower shocks" for
  v0 << 0.1 c. That survey covers v_sh = 0.067-0.267 c, finds eta_nt,e nearly independent of speed
  and ties the drop at its fastest shock to a superluminal effect of the 1-D set-up.
- Agreement in kind with zekovic2025slams: strong transverse fields ahead of a parallel shock block
  the injection of magnetized electrons.
