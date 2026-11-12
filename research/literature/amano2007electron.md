# Amano & Hoshino (2007): Electron injection at high Mach number quasi-perpendicular shocks: surfing and drift acceleration

- Ref: Astrophys. J. 661, 190-202 (2007); DOI 10.1086/513599; arXiv:astro-ph/0612204 (v2); peer-reviewed (ApJ)
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
At a high-M_A quasi-perpendicular shock electrons are injected in two steps. Shock surfing acceleration (SSA) on
Buneman waves at the leading edge of the foot makes a suprathermal tail with mostly perpendicular energy; the part
of the tail outside the loss cone is mirror-reflected at the overshoot (shock drift acceleration, SDA = fast Fermi)
and streams upstream. A semi-analytic model of this chain gives an injected fraction that peaks near theta_Bn = 80 deg.

## Method
- 1-D electromagnetic PIC (code not named), injection / reflecting wall, downstream rest frame.
- m_i/m_e = 100; omega_pe/Omega_ce = 20; beta_e = beta_i = 0.08; v_A = 5e-3 c; injection four-velocity U_0 = 0.05 c;
  M_A ~ 15 in the shock rest frame (upstream speed there U_1 ~ 1.5 U_0); sonic Mach number not stated.
- theta_Bn: "several runs"; Fig. 6-7 show points at 70, 75, 80 deg and the text reports 85 deg (Sect. 3.3). Main run 80 deg.
- 51200 cells of about one Debye length (length in c/omega_pe not stated); 100 particles per cell per species.
- Fig. 1 is at omega_pe t = 12000 = 5.5 Omega_ci^-1; averages are over the last ~1 Omega_ci^-1 of each run (Sect. 3.3).
- Model (Sect. 3): bi-kappa distribution in the foot (Eq. 6); adiabatic mirror reflection including the cross-shock
  potential in the de Hoffmann-Teller frame (Eq. 8, 12-13); corrections for an SSA cut-off energy E_max and an
  escape probability P_esc during the reflection (Eq. 23-24).

## Numbers to compare against
- Buneman threshold M_A >~ sqrt(beta_e m_i/m_e) (Eq. 2). "Unlimited" SSA: 2 M_A sqrt(alpha m_e/m_i) >= 1 (Eq. 25,
  taken from Hoshino & Shimada 2002; alpha = wave energy conversion factor of order one, Eq. 3).
- 80 deg run: SSA raises the energy ~100 times within ~5 Omega_ce^-1, nearly all perpendicular, mu not conserved;
  the front reforms every ~2 Omega_ci^-1 (Sect. 2, Fig. 4). Maximum gamma_e ~ 2 (Sect. 3.1).
- Spectrum in the transition region: power law of index 3-4 for 10 <~ K_e/K_e0 <~ 100; a hump above ~100 only
  upstream and at the leading edge, i.e. reflected, not leaked (Fig. 3). Reflected beam u_par/U_0 ~ 15 (Fig. 2).
- Mirror speed U_s = U_1/cos(theta_Bn) ~ 8.6 U_0; reflection needs u_perp >~ 2.7 U_0 for B_max/B_1 ~ 10 while the
  upstream thermal speed is 0.2 U_0 (Eq. 1, Sect. 2): thermal electrons cannot be reflected without SSA.
- Potential matters when phi~ (m_i/m_e) cos^2(theta_Bn) >~ 1 (Eq. 11): theta_Bn <~ 81 deg at m_i/m_e = 100, <~ 88 deg at
  1836, for phi~ = e phi/(m_i V_1^2/2) = 0.4. A reduced mass ratio therefore overestimates reflection (Sect. 3.2).
- Reflected density over upstream density (Fig. 7, read off a clean plot, approximate): 5e-2 (70 deg), 1.5e-2
  (75 deg), 2e-3 (80 deg); none at 85 deg, where the uncorrected model gives ~4e-3 (Sect. 3.3).
- Model parameters: kappa = 2.5 (index 3.5), phi~ = 0.4, B_f/B_1 = 2.5, B_max/B_1 = 10; fitted E_max = 100 (units of
  upstream electron bulk energy) and P_esc = 4 Omega_ci (Sect. 3.3-3.4).
- SNR extrapolation (Sect. 3.5): needs E_max >~ 750; with E_max = 1000 (10-100 keV at V_1/c ~ 1e-2) and P_esc = 4 the
  injected fraction peaks at ~2e-4 near 80 deg, energy density ~10 % (paper's normalisation); cut-off angles ~85 deg
  (from P_esc) and ~88 deg (from E_max). Set against ~1e-3 and ~30 % of the thermal energy density inferred for
  SN 1006 (Bamba et al. 2003), and against theta_Bn >~ 80 deg required there by the filament widths.
- Reflected beam can drive Alfven waves if v_r/V_1 >~ m_i/(2 M_A m_e) (Eq. 26): >~ 1-10 for M_A = 100-1000 (Fig. 9).

## Assessment
- Shown by simulation: the SSA-then-SDA sequence on single trajectories and by back-tracing reflected electrons;
  low-energy electrons are never reflected; reflected electrons at 70-80 deg and none at 85 deg; the beam drives
  Langmuir waves upstream (absent at 90 deg).
- Model, not simulation: P_esc is a free parameter fitted to the runs and attributed to self-reformation; E_max is
  assumed to grow with mass ratio (qualitative argument only); all SNR numbers are extrapolations.
- Stated limits: 1-D (no oblique Buneman modes, no whistlers, no ripples); reduced mass ratio; weakly relativistic;
  scattering of the reflected beam and DSA are not followed; the index may depend on omega_pe/Omega_ce, theta_Bn and
  m_i/m_e and varies with reformation phase; unknown whether reformation survives at SNR shocks (Sect. 4).
- Works discussed: Levinson (1992) whistler route at parallel shocks for M_A > 43/sqrt(beta_e); Vandas (2001):
  adiabatic SDA alone does not match bow-shock electrons; SDA alone needs unrealistic seed energies at high M_A.
