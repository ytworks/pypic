# Gupta, Caprioli & Spitkovsky (2025): Speed-dependent threshold for electron injection into diffusive shock acceleration

- Ref: journal reference not checked (the arXiv metadata gives none, and no DOI); arXiv:2506.09134v1,
  10 June 2025. Preprint only, as far as the metadata says.
- Read: 2026-11-12, full (with appendices A-C); topics: spherical-shock-injection

## Claim
Electrons enter DSA once their speed, not their momentum or Larmor radius, lets them outrun the shock
along the field: u > v_sh. From there on they gain momentum per cycle as DSA predicts. Heating, mirror
reflection, a speed threshold and DSA losses reproduce the PIC electron spectra and give K_ep << 1.

## Method
- PIC, Tristan-MP, set-up of gupta2024electron (quasi-1-D, upstream frame, piston v_pt = 0.1 c,
  v_sh = 0.133 c, theta_Bn = 30 deg, Te = Ti): M_A = 20 and 5, M_s = 40, m_i/m_e = 100, run to
  ~300 /omega_ci. Checks with other mass ratios (1836 in App. B), resolution and particle number are
  mentioned. beta not stated here (1 and 0.06 for these runs in gupta2024electron).
- ~2e5 tracked electrons: momentum change per upstream-downstream-upstream cycle (Fig. 2, App. A).
- Test electrons in field snapshots give the diffusion coefficient D_e(p) (App. B); test electrons at
  a prescribed tanh ramp give the reflected fraction (App. C).
- Monte Carlo "minimal model": thermal electrons with M_s,e = 2.5; M_A = 20 with theta_Bn ~ 45 deg
  (mean inclination in strong turbulence) and M_A = 5 with 30 deg (Sect. III, Fig. 4).

## Numbers to compare against
- Threshold: the mean dp/p follows Eq. 1, G ~ (2V/3u)[1 + sqrt(1 + (V/u)^2) + 3V/2u], for all
  electrons with u > v_sh (Fig. 2), slightly better with V = v_sh than with V = v_pt. Same at
  M_A = 5 (Fig. 3). Protons: v_inj >~ 3 v_sh, p_inj ~ 3 m_i v_sh (Sect. I, after Caprioli et al. 2015).
- Reflection: u_min = v_sh sec(theta_Bn) / sqrt{[1 - (v_sh sec(theta_Bn)/c)^2](R_B - 1) + 1} (Eq. 4);
  needs sqrt(3) v_th,e >~ u_min, i.e. an effective M_s,e <~ 3 just upstream. Then ~10 % of the electrons
  are reflected at quasi-parallel shocks and <~ 1/2 of those are injected (App. C). The reflected
  fraction falls steeply for M_s,e > 3 and with obliquity (Fig. 11).
- Fraction measurable above the downstream thermal peak (p >~ 3 sqrt(m_e/m_i) m_i v_pt):
  eta_e ~ 0.002 (eta_inj/0.05) / sqrt(m_R/100) (Sect. III).
- K_ep = f_e/f_i ~ (eta_e/eta_i) (p/p_inj)^(q_i-q_e) m_R^(-(q_e-3)/2) ~ 1e-3 for eta_e ~ 1e-3,
  eta_i ~ 1e-2, q = 4, m_i/m_e = 1836 (Eq. 5).
- P_rem ~ 1 - 4 v_pt / ((R - 1) u); q = 3 - ln(P_rem)/ln(1 + <dp/p>) (App. A.2).
- p_max,e/(m_i v_pt) ~ 9 at M_A = 20 (D_1,e ~ D_Bohm(B0), D_2,e = 0.1 D_Bohm(B0)) and ~0.9 at
  M_A = 5 (~0.8 at m_R = 1836) (Eq. 3, Sect. II D, App. B). Electrons reach the proton p_max only for
  t >> (c/v_pt)^2 /omega_ci.

## Assessment
- Shown by simulation: the energy gain per cycle matches DSA above u = v_sh, for two Mach numbers at
  m_R = 100 and one obliquity; trajectories show mirroring and return to the shock (Figs. 1-3).
- Model, not first principles: the Monte Carlo spectra start from M_s,e = 2.5 at the shock and use an
  imposed upper cut-off; the heating that brings M_s,e below 3 is not modelled (referred to Gupta et
  al. 2024, ApJ 968, 17), and the thermal peak is not reproduced (Fig. 4).
- Extrapolation: K_ep ~ 1e-3 at the real mass ratio follows from Eq. 5, not from a run. Strong shocks
  with M_s >> 3 sqrt(m_i/m_e) need non-adiabatic pre-heating before any electron reflects (Sect. III).
- Stated limits and open points (Sect. IV): weak shocks; more oblique shocks, where transport is more
  complicated; non-linear magnetic fluctuations upstream and downstream (zekovic2025slams is cited
  here); weak self-generated turbulence may stop electrons from returning to the shock.
  Superluminal structures or drifting scattering centres would steepen the slope (App. A.2).
- Disagreements: against the usual momentum or Larmor-radius threshold (Blasi et al. 2005; Malkov &
  Drury 2001; Kang et al. 2002). The eta_e quoted here (0.002 at m_R = 100, falling as m_R^-1/2)
  differs from the eta_nt,e <~ 0.5e-2 "independent of m_R" of gupta2024electron, although that paper
  is cited for it.
