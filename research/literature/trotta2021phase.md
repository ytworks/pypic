# Trotta, Valentini, Burgess & Servidio (2021): Phase space transport in the interaction between shocks and plasma turbulence

- Ref: Proc. Natl. Acad. Sci. USA 118, e2026764118 (2021) (Crossref); DOI 10.1073/pnas.2026764118; arXiv:2011.10417
  (v1, titled "The interaction between shocks and plasma turbulence: phase space transport", read here);
  peer-reviewed
- Read: 2026-11-13, full (arXiv v1; the PNAS version may differ); topics: shock-alfvenic-turbulence

## Claim
Pre-existing, fully developed MHD turbulence carried into an oblique supercritical shock strongly enhances proton
phase-space transport upstream: the field-aligned beam (FAB) broadens and shifts, the maximum energy grows, and
core and beam mix. A coarse-grained Vlasov flux through a velocity circle maps where energization and deceleration
occur; it anti-correlates with E_parallel (correlation about 0.7 in magnitude).

## Method
- Two stages: 2.5D compressible MHD decaying turbulence (Perri et al. 2017), windowed and imposed as the upstream
  of a 2D hybrid PIC shock (HYPSI; kinetic ions, fluid electrons); injection method, reflecting wall.
- Mean field B0 in the x-y simulation plane at theta_Bn = 45 deg; dB/B0 = 0, 0.4, 0.8, 2.1.
- V_in = 3.5 v_A, M_A about 5.5; beta_i = 1; box 256 x 256 d_i; dx = 0.5 d_i; dt = 0.01 Omega_ci^-1;
  > 500 ions per cell upstream.
- Coarse-graining: box filter l = 5 d_i (inertial range), velocity radius w = 10 v_A (between core and beam);
  2D-2V reduction (integrated over v_z).

## Numbers to compare against
- Upstream ion spectra: maximum energy rises with dB/B0; FAB peak shifts and spreads; a low-energy (decelerated)
  population appears (Fig. 1d). No numbers given in the text.
- Correlation coefficient between velocity-space flux and E_parallel: C about 0.7 (anti-correlated; Fig. 4).
- PDFs of spatial and velocity-space fluxes widen with dB/B0; even small turbulence is an efficient mixer (Fig. 5).

## Assessment
- Ions only; electrons are a massless fluid, so nothing about electron injection.
- 2D with in-plane mean field at 45 deg: field lines lie in the plane, so the ignorable-coordinate restriction of
  Jokipii et al. (1993) applies to perpendicular transport across the out-of-plane direction.
- The turbulence patch is finite in space and time (windowed); no statement on how long the shock sees it here
  (Trotta et al. 2024 discusses this limit).
- The coarse-grained flux is a diagnostic of where energy changes, not a separation of mechanisms (meandering vs
  ripple vs wave-particle); E_parallel correlation is the only mechanism hint.
- Stated future work: full phase space and kinetic electrons.
