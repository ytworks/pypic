# Demidem, Nattila & Veledina (2023): Relativistic collisionless shocks in inhomogeneous magnetized plasmas

- Ref: Astrophys. J. Lett. 947, L10 (2023) (Crossref); DOI 10.3847/2041-8213/acc84a; arXiv:2212.06053;
  peer-reviewed
- Read: 2026-11-13, full (arXiv v2); topics: shock-alfvenic-turbulence

## Claim
A relativistic magnetized perpendicular pair shock (sigma = 0.1), which does not accelerate particles in a uniform
upstream, does so when the upstream carries kinetic-scale density waves with drho/rho >~ 0.5: the waves corrugate
the front and drive downstream shear vortices that accelerate particles stochastically.

## Method
- PIC (Runko), pair plasma, Gamma0 = 10, sigma0 = 0.1, B0 out of plane (perpendicular, superluminal); reflecting
  wall; 2D production runs, some 3D checks (skin depth 10 cells in 2D, 3 in 3D; >= 16 particles per cell per
  species).
- Upstream density wave amplitude delta = 0.25, 0.5, 0.75; lambda_y = 32 d_e (also 16 d_e), k_x = k_y/100;
  L_y = min(256 d_e, 8 lambda_y) in 2D.

## Numbers to compare against
- Downstream rms dB/B0 about 0.1 for delta = 0.5 (stated as the lowest amplitude probed in kinetic relativistic
  turbulence).
- delta = 0.25: only a small excess and a cut-off about 30 % higher than the planar case; delta = 0.75: power law
  p about 2.5 with gamma_max proportional to t (Bohm-like) (Figs. 3-5; spectra in 400 <= x/d_e <= 900).
- First-crossing energy kick larger by about 2 in the corrugated case; acceleration most efficient when
  lambda_y is comparable to the downstream Larmor radius (about 10 d_e here).
- Requirement stated: front displacement delta X_sh / lambda_y about 1.

## Assessment
- Relativistic pair plasma and compressive (density) perturbations at kinetic scales: not our regime (electron-ion,
  nonrelativistic, Alfvenic, L_c >= 10 d_i). Relevant as a fully kinetic demonstration that upstream
  inhomogeneity can unlock acceleration that is absent for a laminar shock.
- The mechanism is downstream turbulent (shear) acceleration, not field-line meandering at the front.
- Out-of-plane B0 in 2D: particles tied to field lines in the plane-normal sense (Jokipii et al. 1993).
