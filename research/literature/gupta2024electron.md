# Gupta, Caprioli & Spitkovsky (2024): Electron acceleration at quasi-parallel non-relativistic shocks: a 1D kinetic survey

- Ref: ApJ 976, 10 (2024), DOI 10.3847/1538-4357/ad7c4c (as printed in the reference lists of
  gupta2025speed and jikei2026dependence; the arXiv metadata says "Submitted to ApJ");
  arXiv:2408.16071v1 (the version read). Peer-reviewed publication.
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
Quasi-parallel shocks inject both species into DSA for all parameters surveyed: about 0.5 % of the
electrons and 1 % of the protons, carrying ~2 % and ~10 % of the post-shock energy. The tail of each
species starts at p ~ 3 sqrt(m_alpha/m_i) m_i v_pt. Whether the electron tail keeps growing depends
on M_A (it stalls for M_A <~ 10), not on the mass ratio.

## Method
- PIC, Tristan-MP, quasi-1-D (5 cells in y), upstream rest frame, moving reflecting piston.
- theta_Bn = 30 deg, B0 in the x-y plane. 16 runs (Table 1): piston speed v_pt/c = 0.05, 0.1, 0.2
  (v_sh = 4 v_pt/3 = 0.067-0.267 c); m_i/m_e = 16, 100 (default), 400, 1836; M_A = 5-40;
  M_s = v_sh/v_th,i = 5-160 with v_th,i = sqrt(k T/m_i); beta = 4 (M_A/M_s)^2 = 0.06-64; T_e = T_i.
  Mach numbers use v_sh in the upstream frame for a compression of 4 (+-10 % in practice).
  Fiducial run B1: v_pt = 0.1 c, m_i/m_e = 100, M_A = 20, M_s = 40, beta = 1.
- dx = d_e/10, dt = 0.045 /omega_pe, 200 particles per cell per species, 32 current-filter passes;
  receding right boundary, box up to ~1e4 d_i for 300 /omega_ci. omega_pe/Omega_ce not stated.
- Periodic-box runs of the non-resonant streaming instability for three mass ratios (App. A).

## Numbers to compare against
- Thermal/non-thermal boundary: p_0,alpha = 5 p_th,alpha ~ 3 sqrt(2 mu_alpha m_alpha/m_i) m_i v_pt,
  mu_alpha = T_alpha/(T_i + T_e) ~ 1/2 (Eq. 3; Fig. 1 d1). eta_nt = number fraction above p_0 (Eq. 6);
  eps_nt = share of the post-shock energy density above p_0 (Eq. 8); both over 1000 d_i downstream.
- eta_nt,i ~ 1e-2 and eta_nt,e <~ 0.5e-2, nearly independent of v_sh and m_i/m_e; eta_nt,e drops ~5x
  at low M_A or M_s (Fig. 7, Sect. 6.1). eps_nt,i ~ 10 %; eps_nt,e ~ 2 %, ~1 % at m_i/m_e = 1836
  (Sect. 6.2.1, Fig. 8).
- K_ei = f_e(p)/f_i(p) ~ (eta_nt,e/eta_nt,i) (p/p_0,i)^(q_i-q_e) m_R^(-(q_e-3)/2) (Eq. 7):
  <~ 1e-2 for m_R = 1836 and q ~ 4.
- Slope ~4 only for M_A, M_s >~ 10; for M_A <~ 10 little growth beyond p ~ m_i v_pt (Fig. 5, Sect. 8).
  v_pt = 0.2 c: tail saturates at ~10 m_i v_pt once patches exceed theta* = acos(v_sh/c) (Sect. 5.1).
- Te/Ti ~ 0.3-0.5 downstream (Sect. 6.3, Fig. 9); eps_B ~ 0.03-0.06 (Eq. 9, Fig. 8, Sect. 8).
- Non-resonant (Bell) regime expected for 15 <~ M_A <~ 200; seen for M_A = 20-40, resonant modes
  for M_A = 5 and 10 (App. B, Fig. 11).

## Assessment
- Shown by simulation: the scalings above, at one obliquity. M_A = 5 runs with m_R = 100 and 1836 have
  comparable cut-off energies (Fig. 6); the steeper tail at m_R = 100 is attributed to the thermal
  peak lying close to the cut-off, not to a different injection efficiency (Sects. 5.3, 7).
- Stated without a model: electrons are reflected at the shock and scattered by proton-driven
  turbulence; the mechanism is deferred to later work (see gupta2025speed).
- Extrapolation: K_ei at m_R = 1836 comes from Eq. 7; the m_R = 1836 runs are short (fields averaged
  over 100-150 /omega_ci, where the streaming instability has not yet saturated; Sect. 4, App. A).
- Stated limits: 1-D exaggerates the superluminal stall because B_x cannot change and the front
  cannot ripple; whether the precursor shortens with speed in multi-D needs checking (Sect. 5.1).
  Open (Sect. 7): long-term regulation of injection once the tail spans many decades; weakly
  magnetized non-relativistic shocks. Only theta_Bn = 30 deg.
- Disagreements: with shalaby2022mechanism - at M_A = 20 > sqrt(m_R)/4 the efficiency is >> 10x
  theirs, so the intermediate-scale threshold "may not be a threshold for electron injection";
  possible causes named are theta_Bn = 0 and M_s ~ 600 (beta <~ 0.004) there (Sect. 7). eta_nt,e is
  <~ 2x above Park et al. (2015) (frame choice) and ~2x below the theta_Bn = 63 deg runs of Xu et al.
  (2020). Relativistic subluminal shocks give ~4x higher eta_nt,e and eps_nt,e (Sect. 7).
- eta_nt,i ~ 1e-2 (above 5 p_th) is 10-100x the eta ~ 1e-4-1e-3 (tail from xi ~ 3.2-3.7) of
  arbutina2021nonlinear; Sect. 7 nevertheless calls the two in qualitative agreement.
