# Trotta, Pezzi, Burgess, Preisser, Blanco-Cano, Kajdic, Hietala, Horbury, Vainio, Dresing, Retino, Marcucci, Sorriso-Valvo, Servidio & Valentini (2023): Three-dimensional modelling of the shock-turbulence interaction

- Ref: Mon. Not. R. Astron. Soc. 525, 1856-1866 (2023) (Crossref); DOI 10.1093/mnras/stad2384; arXiv:2305.15168;
  peer-reviewed
- Read: 2026-11-13, full (arXiv v1; Sect. 3.3-3.4 skimmed); topics: shock-alfvenic-turbulence

## Claim
First fully 3D hybrid simulation of a supercritical perpendicular shock propagating into pre-existing MHD
turbulence. Turbulence adds front corrugation at the turbulent scales on top of the self-generated ion-scale
ripples; ripples survive at dB/B0 = 0.5 but become hard to separate at 1.0. Across the shock the fluctuation level
rises and the spectral anisotropy of the upstream turbulence is reduced.

## Method
- 3D compressible pseudo-spectral MHD, 256^3, decaying turbulence with dB/B0 = 0.5 and 1.0, B0 along z; snapshot at
  maximum <|j|^2> imposed (B and ion flow) on the upstream of a hybrid shock (HYPSI, CAM-CL; adiabatic fluid
  electrons). Perturbation windowed to zero at the boundaries.
- theta_Bn = 90 deg (B0 along z, normal along x); V_in = 4.5 v_A, M_A about 6; beta_i = 1.
- Box 128^3 d_i^3; dx = 0.5 d_i; dt = 0.01 Omega_ci^-1 (B substeps dt/10); > 50 ions per cell; end time
  20 Omega_ci^-1.

## Numbers to compare against
- Upstream spectra: short k^-5/3 range perpendicular to B0, about k^-2.8 at sub-ion scales; downstream: k^-1
  range then steepening; anisotropy reduced, most for dB/B0 = 0.5 (Fig. 2; at T Omega_ci = 14, upstream
  75-90 d_i, downstream 105-120 d_i, shock at 95 d_i).
- Local theta_Bx along the front (front = where B > 3 B0): departures up to about 40 deg already without
  turbulence (ripples); with turbulence departures spread over a wider range of scales and short-wavelength
  ripple departures weaken (Fig. 4, PDF).
- Unperturbed downstream parallel spectrum: bump at k_z d_i about 1 from ripples.

## Assessment
- Ions only; no energetic-particle or electron analysis. Useful as the reference for front corrugation in 3D with
  an imposed turbulent upstream, and for the local-theta_Bn spread an electron sees.
- Correlation length of the turbulence in d_i is not stated explicitly (MHD box 2 pi L_A mapped to the hybrid box);
  with 128 d_i and an MHD spectrum of limited range, large-scale field-line meandering over >= 10 d_i is only
  partly represented. Run time 20 Omega_ci^-1 is short.
- Ripple patches with different orientations are claimed to matter for acceleration; not tested.
