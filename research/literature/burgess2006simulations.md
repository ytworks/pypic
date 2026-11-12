# Burgess (2006): Simulations of electron acceleration at collisionless shocks: the effects of surface fluctuations

- Ref: Astrophys. J. 653, 316-324 (2006); DOI 10.1086/508805; arXiv:astro-ph/0610714; peer-reviewed
  (the arXiv header says accepted by ApJ).
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
At a supercritical quasi-perpendicular shock the front is rippled on ion scales, and these moving ripples
scatter suprathermal electrons inside the shock layer. Mirror reflection (fast Fermi) then becomes a
multi-step process: higher maximum energies, power-law spectra, and acceleration over a wider range of
theta_Bn than adiabatic reflection at a planar steady shock allows.

## Method
- Two steps: 2-D hybrid run (kinetic ions, massless adiabatic electron fluid; method of Matthews 1994),
  all fields stored; then non-relativistic test electrons integrated in those fields (4th order, adaptive
  step 0.1 to 2e-5 Omega_ce^-1, bicubic interpolation in space, linear in time).
- Shock by injection against a reflecting wall; periodic in y. Box 100 x 20 c/omega_pi, cells 0.2, time
  step 0.01 Omega_ci^-1, 200 ions per cell, small ad hoc resistivity.
- beta_i = beta_e = 0.5. Cases: A in-plane B, theta_Bn = 87 deg, inflow 3.5 v_A, shock M_A = 5.1;
  B out-of-plane B, 87 deg, M_A = 5.1; C in-plane, 87 deg, M_A = 2.1; D in-plane, 80 deg, M_A = 5.1;
  further in-plane runs at 85 and 88 deg (Sect. 3.1, 3.2). Sonic Mach number: not stated.
- c/v_A = 5000 (v_A = 60 km/s; printed as "v_A/c of 5000"). m_i/m_e for the test electrons: not stated.
- Electrons: mono-energetic shells, E_i = 100 eV to 10 keV, released 5 c/omega_pi upstream at
  T = 10 Omega_ci^-1, collected at +7 (upstream) and -4 c/omega_pi (downstream). Total run time: not stated.

## Numbers to compare against
- Planar-like cases (B out of plane, or M_A = 2.1): energy gain only a factor 2-3, sharp cut-off, as in
  adiabatic reflection theory (abstract; Fig. 5, 7).
- Case A, E_i = 1 keV: maximum energy almost 10 times that of cases B and C; nearly flat differential
  energy flux from about 2 to 6 keV; above 6 keV a power law of slope about -4.6 (Fig. 5).
- Plateau reaches 2-5 times the injection energy (Sect. 4); clearest for E_i = 500 eV to 2 keV (Fig. 6).
- Slope against obliquity at E_i = 1 keV: -3.1 at 88 deg, -7.0 at 80 deg (Fig. 10). Intermediate-energy
  flux about equal for 85-88 deg and lower by a factor 0.1-0.5 at 80 deg (Fig. 10).
- theta_Bn = 88 deg, E_i = 5 and 10 keV: one power law over 5 decades in flux up to about 100 keV (Fig. 8).
- Upstream and downstream spectra coincide above the plateau when the front is rippled (Fig. 11); without
  ripples the downstream electrons are hardly energized (Fig. 12).
- Ripple pattern changes on 0.1 to 1.0 Omega_ci^-1; ripples travel along the front at about the Alfven speed
  of the overshoot (Sect. 3.1, after Lowe & Burgess 2003). Ripple wavelength: not quantified.

## Assessment
- Shown by simulation: the contrast between rippled (in-plane field, M_A = 5.1) and unrippled shocks at
  equal Mach number isolates the ripples as the cause. Shown: the local normal direction changes along
  the front and in time, so the local theta_Bn fluctuates.
- Interpretation, not measured: first mirror energization, then scattering that keeps electrons in the
  layer and gives stochastic gains (Sect. 4). No trajectories or scattering rates are shown.
- Conjecture: a low-Mach shock moving through turbulence could behave like a rippled one (Sect. 4).
- Stated limits: no electron-scale fluctuations (hybrid fields), so not valid for shocks with strong
  electron heating; mono-energetic injection only; 2-D; needs checking with 2-D full-particle runs.
- On full PIC it notes that reduced mass ratios can change the behaviour (Scholer et al. 2003).
- The front is planar on average. Curved shocks appear only in the introduction: Krauss-Varban & Burgess
  (1991), Vandas (1995, 2001), Krauss-Varban (1994), Savoini & Lembege (2001).
- Trotta & Burgess (2019) later find that 2-D ripples trap electrons too well; the plateau is partly a
  2-D artefact.
- Doubts: Fig. 11 text names 88 and 85 deg, its caption 87 and 85 deg; the summary prints the energy
  range as "100 - 10 keV".
