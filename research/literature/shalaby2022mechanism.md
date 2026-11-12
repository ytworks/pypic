# Shalaby, Lemmerz, Thomas & Pfrommer (2022): The mechanism of efficient electron acceleration at parallel non-relativistic shocks

- Ref: ApJ 932, 86 (2022) (volume and page as printed in the reference lists of gupta2024electron and
  gupta2025speed; the arXiv metadata gives only "Accepted to ApJ"); DOI 10.3847/1538-4357/ac6ce7;
  arXiv:2202.05288v2. Peer-reviewed publication.
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
At a strictly parallel shock, electrons are pre-accelerated efficiently only when the
intermediate-scale instability (Shalaby et al. 2021) can grow at the shock front, i.e. when the
Alfvén Mach number there is below sqrt(m_i/m_e)/2. A higher M_A or a reduced mass ratio switches the
instability off and with it the electron acceleration; a reduced mass ratio also gives wrong electron
and ion heating. Hence shock PIC runs "require realistic mass ratios".

## Method
- 1D3V PIC, code SHARP; shock by reflection at x = 0 (rest frame of the contact discontinuity).
- Three runs (Table 1): Ma5Mr1836 (M_A = 5.3, m_i/m_e = 1836), Ma5Mr100 (5.3, 100), Ma21Mr1836
  (21.3, 1836). M_A = v_sh/v_A, with v_sh the shock speed in the upstream rest frame.
- theta_Bn = 0. Upstream speed 0.1 c in the simulation frame; v_sh = (4/3) 0.1 c = 0.133 c.
- M_s ~ 365 (sound speed from Te + Tp, adiabatic index 5/3); Te = Tp = 4e-8 m_i c^2/k_B; upstream
  beta_e = beta_i = 1.28e-4 (M_A = 5.3), ~2e-3 (M_A = 21.3) (Sect. 7.1). omega_pe/Omega_ce not stated.
- 10 cells per upstream skin depth (species not specified), c dt = 0.45 dx, 200 particles per cell per
  species. Right boundary expands; box length not stated. Spectra are shown at t Omega_i = 260
  (Fig. 4) and 325 (Fig. 8); the total run time is not stated in the text.
- Linear theory: cold-beam dispersion relation of parallel modes at the shock front (Eq. 1, Fig. 2).

## Numbers to compare against
- Instability condition at the front: M_Af < sqrt(m_r)/2 (Eq. 2). With a front density of 4 n0,
  M_Af = 2 M_A: 10.6 < 21.4 (satisfied); 40.6 > 21.4 and, for m_r = 100, 10.3 > 5 (violated).
- "Non-thermal" = u > 5 u_m, u_m the peak of u^4 f(u) (Sect. 6, App. A.4); efficiency
  eps_sh = E_nonth / [0.5 (m_e + m_i) v_sh^2] (Eq. A10). No injection momentum is defined.
- Shock front (Sect. 6, Fig. 7): non-thermal electrons carry 0.05-0.1 % of the upstream kinetic energy
  in Ma5Mr1836 and 0.01-0.05 % in Ma21Mr1836; m_r = 100 lowers it by 2-3 orders of magnitude.
  Non-thermal ions: 20-30 %, roughly independent of m_r.
- Downstream (> 200 c/omega_i behind the shock): the non-thermal electron-to-ion energy ratio K_ei is
  2-3 orders of magnitude larger in Ma5Mr1836 than in the other two runs (Fig. 7); the conclusions
  quote an electron efficiency lower by "more than 2 orders of magnitude" at M_A = 21.3.
- Te/Ti downstream = 0.4 (both m_r = 1836 runs) and 1.4 (m_r = 100) (Sect. 5.2, Fig. 6);
  k_B T_i ~ 2.5 m_e c^2 against the MHD value 3.06 m_e c^2 (Eqs. 4-5). Density jump n/n0 ~ 6 (Sect. 3).
- SN 1006 (Sect. 7.3): v_s = 3000 km/s, B ~ 35 muG, n = 0.12 cm^-3 give a post-shock M_A = 15.2 <
  sqrt(m_r)/2 (Eq. 6); B_min = 22.7 muG (v_s/3000 km/s) (n/0.1 cm^-3)^(1/2) A^(-1/2) (Eqs. 7-8).

## Assessment
- Shown by simulation: three 1-D runs. Small-scale (k c/omega_i > 2) transverse field power survives
  into the downstream only in Ma5Mr1836 (Figs. 4-5), the only run with a strong electron tail.
  Shown by linear theory: only there are the sub-ion-scale modes comoving with the upstream (Fig. 2).
- Stated, not shown: that scattering on these waves is the acceleration mechanism ("possibly in a
  first-order Fermi type process", Sect. 7.2); there are no particle trajectories. The scattering
  also involves waves reflected at the contact discontinuity (Sect. 7.2), a feature of the set-up.
- Limits: 1-D, theta_Bn = 0 only, very cold upstream (beta ~ 1e-4), one run per case; the m_r = 100
  run also differs from the reference in its electron thermal momentum. Left to future work:
  low-M_s cluster shocks (Sect. 2.1) and beta_e >~ 1 (Sect. 7.1). Obliquity and 2-D are not discussed.
- Disagreements: with whistler-based injection, which needs high M_A (Riquelme & Spitkovsky 2011;
  Niemiec et al. 2012 found no significant acceleration); with Park et al. (2015), where SDA plus
  Bell-mode scattering relies on accelerated ions, whereas here electrons accelerate before ions do.
  gupta2024electron later finds efficient injection at M_A = 20 and no intrinsic m_r dependence.
- Geometry: the introduction rests on the polar caps of SN 1006 and on the global MHD model of
  Winner et al. (2020), which needs an electron efficiency at least 10x higher for theta_B < 45 deg.
