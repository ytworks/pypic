# Hoshino (2024): A Hard Energy Spectrum in 3D Guide-Field Magnetic Reconnection

- Ref: Phys. Plasmas (2024), volume and pages not in the arXiv record; arXiv:2404.15662
- Read: 2026-11-12, full (arXiv v1); topics: (none yet)

## Claim
In relativistic pair-plasma reconnection (T0/mc^2 = 1, sigma = 2), the nonthermal energy share falls
with the guide field in both 2D and 3D, but less in 3D. In 3D the kappa index stays about 3-4 up to
BG/B0 = 5/2, while in 2D it softens to kappa ~ 8. The paper attributes the 3D difference to patchy,
turbulent reconnection driven by oblique tearing modes across the sheet.

## Method
- 3D PIC (semi-implicit time integration, the author's MRI-study code) on Fugaku. Periodic in x and
  z, conducting walls at y = +-80.7 lambda; box 322.5 x 161.4 x 161.4 lambda on
  4032 x 2016 x 2016 cells; 1.6e11 particles, 81 per cell at the neutral sheet. The 2D companion
  runs use the same x-y box with 972 per cell.
- Harris pair sheet plus uniform BG; Nb = 0.05 N(0), T0/mc^2 = 1 (sigma = 2), u_d/c = 0.52,
  rg/lambda = 0.45, c/omega_p = 6.5 cells. BG/B0 from 0 to 5/2 (0, 1/4, 1/2, 1, 2 in Fig. 5).
- Diagnostics: density isosurfaces and slices (Fig. 1); x-averaged density profiles (Fig. 2);
  B_y(kx, y, kz) spectra and 3D/2D ratios (Figs. 3, 4, 9); total kinetic energy (Fig. 5).
- Spectra cover the whole domain, since 3D has no separatrix surface. They are fitted to
  N_M + N_kappa with cutoff exp(-((gamma - gamma_cut)/gamma_cut)^s) by Levenberg-Marquardt (MPFIT)
  on log bins (Figs. 6-8, App. A). eps_den and eps_ene as in the 2D papers (Eqs. 4-5).

## Numbers to compare against
- Kinetic energy at t VA/lambda = 1000: about 4x the initial value for BG = 0, and the gain
  "remains at a level of 50%" for BG/B0 = 2; 2D and 3D similar (Fig. 5).
- Onset ~ 100 lambda/VA in 2D, ~ 70 lambda/VA in 3D for BG/B0 < 1/4 (drift-kink) (Sec. IV).
- Residual < 2e-1 over wide ranges (Fig. 6, Eq. 3).
- kappa (3D) ~ 3-4 for every BG (Fig. 7b; "kappa ~ 4" in Sec. VII). kappa (2D) is hard for
  BG/B0 < 1 and ~ 8 at BG/B0 = 5/2 (Fig. 7b, Sec. VII). Fits are averaged over
  t VA/lambda = 914-1045; the 1-sigma bars are mostly smaller than the symbols.
- BG/B0 = 2: (T_kappa/mc^2, gamma_cut) = (0.79, 14.1) in 3D and (1.65, 22.9) in 2D (Fig. 8).
- eps_ene ~ 0.9 and eps_den ~ 0.5 for BG/B0 < 1/4 (Fig. 7c). eps_ene(2D) falls from 0.9 at
  BG = 0 to 0.27 at BG/B0 = 5/2; eps_ene(3D) stays ~ 0.6 at strong BG (Fig. 7c, Sec. VI).
- T_M/mc^2 ~ 1 (unheated) in both. T_kappa(2D) rises with BG, while T_kappa(3D) stays flat or falls
  slightly for BG/B0 > 3/2 (Fig. 7a).
- 3D/2D ratio of B_y(kx): < 1 at y = 0 but > 1 at y = 20 lambda for BG/B0 = 1; for BG/B0 = 2 the 3D
  waves are stronger from y = 0 to 20 lambda (Fig. 9; also Fig. 4).
- Peak electric-field energy is 0.1-0.25 of the initial kinetic energy (Sec. VI).

## Assessment
- Shown by simulation: the 3D-2D difference in kappa and eps_ene against BG. Oblique tearing
  (k.B = 0) appears off the neutral sheet, the 3D density profile is broader, and 3D has more wave
  power off the sheet (Figs. 2-4, 9).
- Attributed, not demonstrated: that the turbulence causes the extra acceleration (second-order
  Fermi, or island first-order Fermi as in hoshino2012stochastic). There is no particle tracing and
  no energisation diagnostic, only correlation with wave power.
- Definitions change the answer (acknowledged). At BG/B0 = 2 the thermal and nonthermal parts are
  poorly separated (Fig. 8, Sec. VII), so the 3D kappa rests on the fit model. Whole-domain spectra
  include upstream plasma, so the fractions are not directly comparable with the exhaust-only
  hoshino2022efficiency and hoshino2023energy.
- Assumed without new runs: weak sigma dependence for 1 < sigma < 10, taken from the earlier 2D
  work (Sec. VII).
- Limits (stated): one temperature (sigma), pair plasma, and a narrow dynamic range (a 3D box
  ~10^3x larger would be needed). The final state is taken as saturated because one large island has
  formed by t VA/lambda = 1045. The paper does not discuss the 12x fewer particles per cell in 3D.
