# Arbutina & Zeković (2021): Non-linear diffusive shock acceleration: A recipe for injection of electrons

- Ref: Astroparticle Physics 127, 102546 (2021); DOI 10.1016/j.astropartphys.2020.102546;
  arXiv:2012.15117v1. Peer-reviewed publication.
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
In long PIC runs of quasi-parallel shocks the ion and electron power laws become parallel and both
leave their Maxwellian at about the same multiple xi of the thermal momentum, so the injected number
fractions are similar (eta_e ~ eta_i). Used as the injection rule in Blasi's semi-analytic non-linear
DSA model with constant electron heating, this gives K_ep as a function of Mach number, and ~1:100 at
Mach ~ 100 in the test-particle regime.

## Method
- PIC: TRISTAN-MP, 2-D boxes in the x-y plane, all narrower than an ion gyroradius; "initially
  parallel" field (along the flow, hence in-plane). Four runs (Table I):
  run 1: m_i/m_e = 16, M_A = 13, M_s = 1800, u0 = 0.33 c, width 6.4 ion inertial lengths, 250 /omega_ci;
  run 2: 50, 11, 40, 0.33 c, 7.2, 370; run 3: 100, 13, 1800, 0.33 c, 2.6, 130;
  run 4: 16, 16, 35, 0.4 c, < 1 (about 2.5 c/omega_pe, "almost 1D"), 1930 /omega_ci.
  sigma = 0.6e-3 (1.0e-3 in run 2). u0 is called the shock velocity in the lab frame; the frame of the
  Mach numbers is not stated. Te ~ Ti upstream; beta not stated.
- Fixed box 25600 x 25.6 c/omega_pe (runs 1, 3), expanding box (runs 2, 4); 10 cells per c/omega_pe
  (5 in run 3); 8 particles per cell in total; 32 current-filter passes.
- Model: Blasi (2002, 2004) non-linear DSA for protons and electrons, thermal-leakage injection
  p_inj = xi p_th,2 (Eq. 11), Alfvén-heating parameter zeta, and a constant energy dE ~ 0.3 keV moved
  from protons to electrons ahead of the subshock (Eqs. 32-33).

## Numbers to compare against
- Injection: p_inj is where the power law meets the Maxwellian; xi = p_inj/p_th, p_th = sqrt(2 m k T);
  eta = 4 (R - 1) xi^3 exp(-xi^2) / (3 sqrt(pi)) (Eqs. 4, 12). Recipe: xi_e ~ xi_p, Gamma_e ~ Gamma_p,
  eta_e ~ eta_p (Sect. II, non-relativistic case).
- Run 4 at 1930 /omega_ci: far downstream xi ~ 3.7, eta_i ~ 1e-4, f ~ p^-4.2; near downstream
  xi ~ 3.2, eta ~ 1e-3, slope ~5 (Sects. II and IV, Figs. 10-11).
- K_ep = (eta_e/eta_p) (p_inj,e/p_inj,p)^(3/(R-1)) (Eq. 6) = (m_e T_2e / m_p T_2p)^(3/(2(R_sub-1)))
  for eta_e = eta_p (Eq. 40); temperature ratio from Eq. 41.
- Model at u0 = 5000 km/s, n_H ~ 0.1 cm^-3, T0 = 1e5 K, B0 = 5.3775 muG, M_s = M_A = 135, zeta = 0.5
  (Sect. IV, Fig. 8): xi = 3.3 gives R_sub = 3.093, R_tot = 10.210, K_ep = 0.0007; xi = 4.3 gives
  R_sub = 3.999, R_tot = 4.018, K_ep = 0.0019.
- Test-particle K_ep ~ 1:100 at Mach ~ 100, u0 ~ 1000 km/s, for a sound speed of 10 km/s (Fig. 12).
- Matching the model to run 4 needs electron heating of ~150 eV (xi = 3.2) and ~208 eV (xi = 3.7).

## Assessment
- Shown by simulation: parallel tails and similar eta profiles of ions and electrons late in one long
  run with m_i/m_e = 16 (Figs. 7, 10, 11). Electrons are relativistic there, so xi_e does not match
  xi_i in that run; xi_e ~ xi_p is asserted for the non-relativistic case (Sect. II).
- Assumed: eta_e = eta_p at the real mass ratio and at SNR speeds. The electron heating dE is taken
  from observations (Ghavamian et al. 2007, 2013) and sets K_ep through Eq. 41.
- Conjectured: the pre-acceleration is "Fermi-like" SDA between the shock and upstream waves, almost
  equivalent to SSDA; its microphysics "is still to be shown" (Sect. II). No trajectories.
- Run 2 shows a temporary switch-off of electron acceleration when the precursor wave amplitude
  changes (Sect. II, Fig. 6): injection is not steady.
- Stated limits: boxes narrower than an ion gyroradius, so no shock rippling (Sect. II);
  u0 = 0.33-0.4 c; shocks only weakly modified. Open (Sect. V): whether eta_e ~ eta_p holds at other
  Mach numbers and whether the shock becomes modified; synchrotron losses, Galactic transport.
- Disagreements: against models that inject electrons and protons at the same momentum or energy with
  K_ep set by hand. Mass ratio 16 is taken as sufficient (citing Sironi & Spitkovsky 2011), the
  opposite of shalaby2022mechanism. gupta2024electron finds eta_e <~ 0.5 eta_i, not equality;
  gupta2025speed places the electron threshold at a speed, not at a multiple of p_th.
