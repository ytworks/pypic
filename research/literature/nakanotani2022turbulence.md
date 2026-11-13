# Nakanotani, Zank & Zhao (2022): Turbulence-dominated shock waves: 2D hybrid kinetic simulations

- Ref: Astrophys. J. 926, 109 (2022) (Crossref); DOI 10.3847/1538-4357/ac4781; no arXiv version found;
  peer-reviewed
- Read: 2026-11-13, full (author manuscript PDF from par.nsf.gov, text extraction); topics: shock-alfvenic-turbulence

## Claim
A low-Mach shock moving into fully developed turbulence with no mean field (tearing-unstable force-free fields
injected upstream) amplifies magnetic, velocity and density fluctuations, the Zank et al. (2021) transmission model
reproduces the spectral change, and protons form a DSA-like power law downstream, scattered by the upstream
turbulence.

## Method
- 2D hybrid (kinetic protons, massless electron fluid; CAM-CL with a Boris pusher); injection method, rigid wall,
  simulation frame = downstream frame.
- Upstream: B = (0, B0 cos(2 pi (x+V0 t)/lambda), B0 sin(...)), lambda = 8 c/omega_pp, injected continuously; it
  becomes turbulent by the tearing instability from x about 300; no mean field (dB/B0 effectively >> 1).
- V0 = 2 V_A0; beta_p = 0.4, beta_e = 0.1; box 1000 x 102.4 c/omega_pp; dx = 0.2 c/omega_pp;
  dt = 0.005 Omega_cp^-1; resistivity 2e-3.

## Numbers to compare against
- Shock speed about 1.1 V_A0 (simulation frame); sonic Mach number 2.4 (c_s = 1.3 V_A0 just upstream);
  compression about 2.65 (Fig. 1).
- Across the shock: total turbulence energy x10; spectra of B, u, rho amplified by factors about 10, 2, 10;
  spectral index unchanged except density; sigma_c -> 0 and sigma_r -> 1 downstream; decay within 50 c/omega_pp.
- Proton spectrum (X = 600-900, T = 350): non-thermal from E = 10 to about 100 m_p V_A0^2, index close to the DSA
  value 2.82 for r = 2.65, slightly harder (Fig. 5); one trajectory shows repeated reflection and upstream
  scattering (Fig. 6).

## Assessment
- No mean field: neither quasi-perpendicular nor field-line meandering around a mean direction; distinct from our
  regime (dB/B0 0.3-1 with a defined theta_Bn).
- Electrons are fluid; no electron result. 2D only, narrow box (102.4 c/omega_pp) compared with the 8 c/omega_pp
  current-sheet scale.
- Stated open: whether the power law extends beyond E about 200 in longer runs.
