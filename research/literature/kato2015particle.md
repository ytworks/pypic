# Kato (2015): Particle acceleration and wave excitation in quasi-parallel high-Mach-number collisionless shocks: Particle-in-cell simulation

- Ref: Astrophys. J. 802, 115 (2015) (volume and page as printed in the reference lists of Xu et al. 2020 and
  Kumar & Reville 2021; the arXiv metadata gives only the DOI); DOI 10.1088/0004-637X/802/2/115;
  arXiv:1407.1971v2; peer-reviewed (comment: accepted by ApJ)
- Read: 2026-11-12, introduction, set-up, electron results (Sect. 3.3.2), discussion and conclusions;
  topics: spherical-shock-injection

## Claim
In a long 1-D PIC run of a quasi-parallel (30 deg), M_A ~ 28 shock, accelerated protons drive large Alfvenic
upstream waves that make the shock locally quasi-perpendicular. Protons are injected steadily into a fast,
non-diffusive Fermi-like process; electrons join the same process but only in one burst of injection.

## Method
- PIC, own 1D3V code after Birdsall & Langdon; injection method with two reflecting walls; downstream rest frame.
- m_p/m_e = 30. Flow V = 0.37 c (Lorentz factor 1.08). Shock speed ~0.1 c (downstream frame), ~0.46 c (upstream
  frame). M_A ~ 28 (matches 0.46 c / v_A, v_A = 1.59e-2 c). Sonic Mach number not stated.
- beta_e = beta_p = 0.5 (total 1); T_e = T_p = 1.96e-3 m_e c^2; sigma = 1.9e-3; omega_ce/omega_pe = 8.847e-2.
- Theta_0 = 30 deg, measured in the downstream rest frame; field in the x-y plane.
- Box 1.5e5 c/omega_pe = 2.74e4 c/omega_pi (2.5e6 cells, dx = 0.06 c/omega_pe); dt = 0.05/omega_pe; 160 particles
  per cell per species. Run to 2e5/omega_pe (= 590/omega_cp with the stated omega_cp = 2.95e-3 omega_pe).

## Numbers to compare against
- Upstream waves: appear at ~1.5e4/omega_pe; B_t >> B_0 near the shock; wavelength 200-400 c/omega_pe, about the
  gyroradius of reflected protons; local field angle ~ perpendicular, local M_A ~ 3.6 (Sect. 3.1-3.2, Figs. 1-3).
- Downstream energy spectra at 2e5/omega_pe: power-law index ~2.4 (protons, E_kin > 3 m_e c^2) and ~3.0
  (electrons, E_kin > 5 m_e c^2); maximum E_kin ~1000 (protons) and ~400 m_e c^2 (electrons) (Sect. 3.3, Fig. 5a).
- Non-thermal electrons ~ two orders of magnitude fewer than non-thermal protons (Sect. 3.3).
- Energy shares downstream: non-thermal protons ~10 %, non-thermal electrons ~0.1 % of the total particle
  energy; total proton/electron energy ratio ~2.5 (Fig. 6). T_p = 0.46, T_e = 0.24 m_e c^2 (Fig. 5a).
- Acceleration rates dE/dt ~ 7.5e-3 (protons, Fig. 7) and ~1.1e-2 m_e c^2 omega_pe (electrons, Fig. 13).
- Electron injection: nearly all of the 321 electrons ending above 60 m_e c^2 were injected near 6e4/omega_pe,
  when the wave amplitude dropped after a maximum (Figs. 15-16). The amount of non-thermal electrons falls with
  time (Fig. 5b).
- Pre-acceleration in the foreshock: for t > 1.5e5/omega_pe about half of the injected electrons of the 50-60
  m_e c^2 sample exceed 10 m_e c^2 before reaching the shock (Fig. 17); hot upstream component ~0.35 m_e c^2 (Fig. 18).
- Amano & Hoshino criterion: M_A^inj = 0.8-1.5 for local angle 30-60 deg; subluminal only for Theta < 63 deg
  at 0.46 c (Sect. 4).
- Numerical drag: dE/dt ~ -(1/2N_e) m_e c^2 omega_pe (Eq. 4); acceleration almost stops for N_e = 50 and 133
  particles per skin depth, still affected at 667; the run uses 2667 (Fig. 19).

## Assessment
- Shown by simulation: proton-driven waves; steady proton injection with phase trapping; a Fermi-like cycle whose
  return times are of order a gyration time, not the Bohm diffusion time (Fig. 10); bursty electron injection.
- No reflected electron beam and no sign of whistler scattering were found, so the author concludes the electron
  injection differs from the Amano & Hoshino (2010) model although its Mach-number criterion is met (Sect. 4).
- Conjectured: pre-accelerated foreshock electrons may later give steady injection; the single burst may be an
  artefact of the initial condition (Sect. 4).
- Limits stated: 1-D (the regular, near-monochromatic wave field may be a 1-D effect; 2-D/3-D ion-Weibel modes are
  missing); fast shock, so the growing waves make it superluminal and electron reflection impossible, runs at
  lower speed needed; M_A still below remnant values (~100); small mass ratio.
- Contrast with Park et al. (2015): there a p^-4 tail and standard DSA; here spectra in energy with other indices
  and a non-diffusive process. Neither paper compares the two in detail.
