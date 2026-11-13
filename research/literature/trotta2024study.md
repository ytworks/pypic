# Trotta, Valentini, Burgess & Servidio (2024): A study of the transition to a turbulent shock using a coarse-graining approach to ion phase-space transport

- Ref: Mon. Not. R. Astron. Soc. 536, 2825-2838 (Crossref issue date 2024-12-19; volume 536 is the 2025 volume);
  DOI 10.1093/mnras/stae2750; arXiv:2412.08404; peer-reviewed
- Read: 2026-11-13, full for set-up and results (coarse-graining section skimmed); topics: shock-alfvenic-turbulence

## Claim
With increasing pre-existing turbulence the shock goes from laminar to turbulent: the front becomes irregular, the
y-averaged transition widens, the shock speeds up, downstream structures carry transmitted turbulence, and ions
reach higher maximum energies. The coarse-grained Vlasov flux (Trotta et al. 2021) gives consistent results when the
filter scale is in the inertial range.

## Method
- Same simulation campaign as Trotta et al. (2021): 2.5D compressible MHD decaying turbulence, windowed via the
  vector potential and a Helmholtz-decomposed velocity field, imposed on a 2D hybrid shock (HYPSI, CAM-CL, adiabatic
  fluid electrons, gamma = 5/3).
- B0 in the x-y plane, theta_Bn = 45 deg; dB/B0 = 0, 0.4, 0.8, 2.1 (the text says "three simulations" but lists four
  levels plus the laminar case in figures).
- V_in = 3.5 v_A, M_A about 6 (Trotta et al. 2021 gave about 5.5 for the same set-up); beta_i = beta_e = 1;
  box 256 x 256 d_i, dx = 0.5 d_i, dt = 0.01 Omega_ci^-1, resistivity 0.01, 500 ions per cell, end 90 Omega_ci^-1.

## Numbers to compare against
- Shock speeds in the simulation frame: 1.4, 1.55, 1.7, 2.0 v_A for dB/B0 = 0, 0.4, 0.8, 2.1 (Sect. 3).
- Ion energy spectra in 10 x 256 d_i^2 boxes up- and downstream: higher maximum energy and a low-energy upstream
  tail with turbulence (Fig. 4); compatible with Nakanotani et al. (2022) downstream.
- Velocity-space circle w = 5 v_A in this paper (10 v_A in 2021).

## Assessment
- Ions only; electrons fluid. Same 2D in-plane-field restriction as the 2021 paper.
- The turbulent patch is finite: once it has passed, unperturbed plasma arrives (stated limit), so the shock sees
  turbulence only for part of the run.
- States that pre-existing turbulence may inhibit self-generated upstream fluctuations (foreshock) - shown only
  qualitatively from snapshots (Fig. 3).
