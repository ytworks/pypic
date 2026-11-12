# Park, Caprioli & Spitkovsky (2015): Simultaneous Acceleration of Protons and Electrons at Nonrelativistic Quasiparallel Collisionless Shocks

- Ref: Phys. Rev. Lett. 114, 085003 (2015); DOI 10.1103/PhysRevLett.114.085003; arXiv:1412.0672v2;
  peer-reviewed Letter (journal reference from the arXiv metadata)
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
A 1-D PIC run of a quasi-parallel (30 deg), M_A = 20 shock carries protons and electrons together into diffusive
shock acceleration (DSA) with f(p) ~ p^-4. Electrons get there by shock drift acceleration (SDA) and then a "hybrid"
stage, SDA cycles plus scattering on proton-driven Bell waves. The electron/proton ratio is K_ep ~ 1e-3 to 1e-2.

## Method
- PIC, TRISTAN-MP, 1-D along x, all field and velocity components kept; plasma sent against a reflecting wall, so
  the simulation frame is the downstream frame.
- m_p/m_e = 100 (one run with 400). v_u = 0.1 c (one run with 0.05 c), flow speed in the simulation frame.
- M_A = 20, M_s = 40, both with v_sh = v_u r/(r-1), the upstream speed in the shock frame; M_s = v_sh/sqrt(T_p/m_p).
- T_e = T_p = 1.12e-3 m_e c^2. beta is defined as 2 M_A^2/M_s^2 (Conclusions); the value 0.5 is my arithmetic.
- theta_Bn = 30 deg only. omega_pe/Omega_ce not stated (about 15 from 4.6e5/omega_pe = 310/Omega_cp).
- Box grows to ~4e5 cells at 10 cells per c/omega_pe; dt = 0.045/omega_pe; 200 particles per cell per species.
- Run time 4.6e5/omega_pe ~ 310/Omega_cp (reference run); not stated for the other two runs.

## Numbers to compare against
- Precursor: deltaB_perp/B_1 ~ 2 over ~5e3 c/omega_pe ahead of the shock (Fig. 1d).
- Wave spectrum peaks at k_max = 4e-3 omega_pe/c, above k_res = 1/rho_CR ~ 1.1e-3 omega_pe/c; right-handed at
  k_max (non-resonant Bell), left-handed resonant mode at 1/rho_CR (Fig. 2). Electrons resonate: k_max rho_e ~ 1.
- Protons: injected at p_inj ~ 3 m_p v_u = 30 m_e c after a few gyro-cycles of SDA (Fig. 3a-d).
- Electron track: mirror reflection at t ~ 3.65e5/omega_pe, trapped between shock and nearest upstream wave until
  ~4.25e5/omega_pe; energy from 5e-3 to 60 m_e c^2 (factor > 1e4) in < 8e4/omega_pe ~ 50/Omega_cp (Fig. 3e-f).
- Electron entry into standard DSA at p_inj ~ 50 m_e c (Fig. 3g-h); 30-100 m_e c over other tracks (text).
- Downstream spectra (region 2000 c/omega_pe behind the shock): both species ~ p^-4 after ~4e5/omega_pe ~
  268/Omega_cp; electrons already p^-4 for 2.5 < p/m_e c < 30; T_e ~ T_p ~ 0.12 m_e c^2 = 0.12 m_p v_u^2 (Fig. 4).
- K_ep = f_e(p)/f_p(p) for p > p_inj, independent of p: 3.8e-3 (m_p/m_e = 100, v_u = 0.1 c); 5.5e-3 (400, 0.1 c);
  1.2e-3 (100, 0.05 c) (section "Electron/Proton Acceleration Efficiency", from Fig. 4).
- Extrapolations: K_ep < ~0.01 at the real mass ratio and 0.1 c; 1-3e-3 for young remnants (same section).
- Injected number fractions, energy fractions: not given.

## Assessment
- Shown by simulation: p^-4 tails of both species; Bell-type (right-handed, k > 1/rho_CR) precursor waves; the
  electron sequence reflection -> trapped SDA cycles -> DSA at a momentum close to the proton injection momentum.
  The waves help in two ways: they make the shock locally oblique (larger pitch angles for mirror reflection) and
  they scatter electrons back to the shock.
- Stated from few runs: K_ep roughly proportional to v_sh/c (two speeds); weak rise with mass ratio (two values),
  read as thermal electrons becoming trans-relativistic; K_ep did not vary at later times.
- Conjectured, not simulated: shocks of all inclinations reflect electrons; quasi-perpendicular shocks should make
  electron-driven upstream turbulence and reach electron DSA with virtually no proton acceleration; the absence
  of electron DSA in quasi-perpendicular PIC runs so far is put down to box size.
- Limits it states: one theta, M_A and beta; a scan over M_A, beta, theta, mass ratio and speed is still needed.
  Not discussed: the 1-D geometry itself, which allows only wave vectors along x.
- Relation to others: the hybrid stage is likened to Guo, Sironi & Narayan (2014), where electron-driven firehose
  waves confine electrons at quasi-perpendicular shocks; here the confining waves are proton-driven.
- K_ep is compared with ~1.6e-3 inferred for Tycho (Morlino & Caprioli 2012) and ~0.01 at Earth near 10 GeV.
