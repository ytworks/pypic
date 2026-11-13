# Behar, Pucci, Simon Wedlund, Henri, Ballerini, Preisser & Califano (2024): Impact of solar wind turbulence on a planetary bow shock

- Ref: Astron. Astrophys. 692, A53 (2024); DOI 10.1051/0004-6361/202451520; arXiv:2409.17942;
  peer-reviewed.
- Read: 2026-11-13, skimmed (ar5iv HTML, set-up, figures, limits); topics: shock-alfvenic-turbulence

## Claim
First global 3-D hybrid runs of turbulent solar wind hitting a magnetized planet: turbulence breaks the
coherence of foreshock waves, makes large ripples on the quasi-perpendicular bow shock, and creates
bubble-like structures at the nose; none of these appear with a laminar wind.

## Method
- Hybrid PIC code Menura. Two stages: (1) periodic decaying turbulence (2000 d_i box, 400^3 cells,
  dx = 5 d_i, 600 ppc, dt = 0.5 Omega_ci^-1, hyper-resistivity 0.01), oblique B0, beta_i = 0.5,
  T_i/T_e = 1, isothermal electrons; delta B0/B0 = delta v0/v_A = 0.54, incompressible modes with
  polarization perpendicular to B0 and k, k = 1-5 (2 pi/L_box), random phases. Taken as developed when
  the rms current saturates (t ~ 650 Omega_ci^-1): delta B/B0 = 0.45, delta v/v_A = 0.33, ~ -5/3 spectrum
  (Fig. 1).
- (2) That state (E, B, ion distribution) initializes the planet run. The solver stays in the solar-wind
  frame: the planet moves through the plasma and the box is shifted to follow it (dipole recomputed,
  extra Lorentz term in Faraday's law). No inflow boundary for turbulence is needed; method details
  deferred to Behar et al. (2022) and Behar & Henri (2023) (not read).
- U_sw = 10 v_A, M_A = 10, beta = 1; laminar reference run with identical numerics.

## Numbers to compare against
- Shock-position spread 0.8 d_i (laminar) vs 2.0 d_i (turbulent); quasi-perpendicular oscillations
  +- 10 d_i vs +- 2 d_i (Fig. 6). Comparison at t = 250 Omega_ci^-1.

## Assessment
- Design idea relevant to us: run in the plasma frame and move the obstacle (here a planet; for us a
  piston or a driven wall) through a pre-relaxed periodic turbulent state. This avoids a turbulent
  inflow boundary and the matching of slabs, at the cost of a moving boundary.
- Limits: hybrid, ion scales poorly resolved (dx = 5 d_i); decay of turbulence between start and shock
  not quantified; no electron physics; separation of turbulence from laminar effects left to future
  work (stated).
