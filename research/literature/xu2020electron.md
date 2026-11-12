# Xu, Spitkovsky & Caprioli (2020): Electron acceleration in one-dimensional non-relativistic quasi-perpendicular collisionless shocks

- Ref: Astrophys. J. Lett. 897, L41 (2020) (volume and page as printed in the reference list of Kumar & Reville
  2021; the arXiv metadata gives only the DOI); DOI 10.3847/2041-8213/aba11e; arXiv:1908.07890v2; peer-reviewed
  (comment: accepted by ApJL). The title above is that of v2; the arXiv metadata title omits "one-dimensional".
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
At a 63 deg shock with high M_A and M_s, electrons reflected by magnetic mirroring drive their own upstream waves
(a firehose-type heat-flux instability), are trapped between these waves and the shock for repeated shock drift
cycles, and enter DSA with f(p) ~ p^-4, while protons are hardly accelerated.

## Method
- PIC, TRISTAN-MP, 1-D along x, all field and velocity components; upstream rest frame, conducting wall moved at
  v_0 = 0.15 c into the plasma.
- m_p/m_e = 100 (limited runs at 400). Reference run: M_A = 63, M_s = 55 with M_s = v_sh/sqrt(5/3 (T_p+T_e)/m_p);
  T_p = T_e = 4e-4 m_e c^2. Scan: M_s = 3, 27, 55 times M_A = 11, 28, 63 (Fig. 4), by changing T and B_0.
- beta and omega_pe/Omega_ce not stated (my arithmetic for the reference run: beta ~ 1.6, omega_pe/Omega_ce ~ 32).
- theta_Bn = 63 deg in every run. Final box ~3.5e4 c/omega_pe, 10 cells per c/omega_pe, 256 particles per cell
  per species. Run to ~45/Omega_cp.
- Checks: linear kinetic dispersion relation of a beam-plasma system; 2-D periodic beam-plasma PIC run (Fig. 5);
  one short 3-D shock run, 50 x 50 (c/omega_pe)^2 across, 8 cells per c/omega_pe, to ~3/Omega_ci (Fig. 6).

## Numbers to compare against
- Upstream waves: peak at k ~ 0.1-0.2 omega_pe/c, left-hand circularly polarised, hence non-resonant with
  electrons (Fig. 2b-c). Linear theory with beam fraction n_b/n_0 = 0.03, drift 0.35 c along B, beam 100 times
  hotter than the background, propagation at 63 deg to B: same fastest-growing k (Fig. 2d).
- Proton current upstream below 10 % of the electron current (Sect. 3).
- Electron track: pre-heated in the foot and reflected at t ~ 0.5-0.6e5/omega_pe; trapped SDA cycles for
  ~0.6-0.9e5/omega_pe; standard DSA from p_inj ~ 30-80 m_e c (Sect. 3, Fig. 3).
- Reference run, downstream (200-4000 c/omega_pe behind the ramp): index 4 in f(p); non-thermal number fraction
  eta_e ~ 7 %, energy fraction eps_e ~ 20 % at the end; "non-thermal" = energy above 5 times that of the
  downstream thermal peak (footnote 1); E_max grows about linearly in time (Sect. 3, Fig. 4i).
- eta_e ~ 0.7 % at Mach numbers ~30 (Fig. 4e); delta_B/B_0 > 1 (peak upstream value) only when both Mach numbers
  are large; at low M_A in hot plasma only an SDA tail and no DSA by the end of the run (Fig. 4a, Sect. 3).
- Fig. 4 panel labels as extracted (symbols lost; e and i agree with the text): eta_e = <0.1, 0.1, 0.1 % (M_s = 3),
  0.6, 0.7, 6.3 % (M_s = 27), 0.7, 1.0, 7.1 % (M_s = 55) for M_A = 11, 28, 63; delta_B/B_0 = 0.4, 0.9, 1.6;
  0.2, 0.7, 7.8; 0.4, 0.8, 8.9.
- Instability threshold v_dr/v_A = lambda sqrt(m_i/m_e) beta^kappa, 0 < kappa < 0.1 (Shaaban et al. 2018; Sect. 3).
- 3-D check: electron reflection lower than in 1-D "by a factor of few" at ~3/Omega_ci (Fig. 6).

## Assessment
- Shown by simulation (1-D): mirror reflection; electron-driven oblique waves; SDA cycles between waves and shock;
  a p^-4 electron tail without proton DSA; strong dependence on M_A and weaker dependence on M_s.
- The protons form a steep tail only: the electron-driven waves change the obliquity for a time and let a few
  protons escape upstream. So here the electrons help the protons, not the reverse.
- Stated or conjectured: the best angle for electron acceleration lies between 60 and 70 deg, and electron DSA
  must vanish towards superluminal angles; neither is tested (one angle only). Low-M_s shocks would inject on
  much longer time scales. At the real mass ratio the M_A needed is higher, "several hundred".
- The comparison with quasi-parallel shocks (eta_e < ~0.1 % at Mach ~20, read from Fig. 4 of Park et al. 2015)
  uses different angle, Mach numbers and speed; it is not a controlled comparison.
- Remnant picture offered (Sect. 4): quasi-parallel regions inject ions, amplify the field strongly and are
  X-ray bright; quasi-perpendicular regions inject electrons (eta_e ~ 5-7 %) without ions, to lower energies,
  and are radio bright. This is an interpretation, with SN 1006 named as the example.
- Limits stated: 1-D admits only waves along x (the dominant mode is oblique, which the 2-D beam test supports);
  multi-dimensional variation of the field angle (shock corrugation) can lower the mean injection; m_i/m_e = 100,
  and runs at 400 need higher M_A to drive the waves; enough particles per cell are needed against numerical drag.
- Related work it discusses: Guo et al. (2014) found firehose-mediated acceleration at low Mach number and high
  beta; earlier quasi-perpendicular PIC runs had not shown DSA spectra.
