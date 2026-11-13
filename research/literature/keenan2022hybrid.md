# Keenan, Le, Winske, Stanier, Wetherton, Cowee & Guo (2022): Hybrid Particle-in-Cell Simulations of Electromagnetic Coupling and Waves From Streaming Burst Debris

- Ref: Phys. Plasmas 29, 012107 (2022) (volume, article number and year from the Crossref record of the
  DOI); DOI 10.1063/5.0075482; arXiv:2112.12114 (v1, 22 Dec 2021, "accepted by Physics of Plasmas");
  peer-reviewed. The arXiv title reads "Paricle-in-Cell" (sic).
- Read: 2026-11-13, skimmed (abstract, set-up, conclusions); topics: spherical-shock-injection
- Source: the arXiv full text (HTML of v1). Read: abstract, Sects. I-II (3-D set-up and result), Sect. III
  (linear analysis, by figure), Sect. IV (2-D finite beams) and Sect. V (conclusions).

## Claim
A point-like explosion of debris into a magnetized background couples in two ways. Across the field most of
the debris is stopped by its Larmor interaction; a diamagnetic cavity forms and a shock leaves it across the
field. Along the field the coupling is weak: debris beams stream out and drive ion-ion beam instabilities
(right-hand resonant, RHI; non-resonant, NRI) with a transverse structure that the 1-D problem does not
have. Whether a quasi-parallel shock forms depends on the length, width, density and speed of the beam; in
the runs shown, none is identified without doubt.

## Method
- Hybrid-VPIC: kinetic ions, fluid electrons. No kinetic electrons.
- 3-D run (Sect. II): 1296 x 540 x 540 cells of 1 d_i; dt omega_ci = 0.005; ~150 particles per cell;
  background beta = 0.1; B along x. Debris: m_d/m_i = 3, a Gaussian cloud, Alfvén Mach number
  M_A = u_d/v_A = 5, equal-mass radius R_m = 40 d_i, directed gyroradius rho_d = u_d/Omega_cd = 15 d_i.
- 2-D runs of finite debris beams along the field (Sect. IV): 2700 x 720 cells of 1 d_i; m_d/m_b = 3,
  n_d0/n_0 = 0.05, u_dx/v_A = 7.5 in the text (the caption of Fig. 7 gives 0.03 and M_A = 7).

## Numbers to compare against
Conditions for coupling and for a shock:
- Coupling across the field: rho_d/R_m = 0.375 with M_A = 5, "chosen to clearly reproduce strong debris ion
  to background coupling" (Sect. II). No threshold is derived here; for that see `le2021astrophysical`.
- Cavity: "reaches its peak size of ~R_m = 40 d_i at t omega_ci ~ 25"; at t omega_ci = 62 it "has already
  begun collapsing" (Sect. II, Fig. 1).
- Perpendicular shock: "a shock wave propagates out across the magnetic field" (Fig. 1); at t omega_ci = 90
  "A shock front propagating primarily in the z direction perpendicular to the background magnetic field is
  also clearly visible" (Fig. 2a). Its Mach number and its formation time are NOT stated. The debris Mach
  number is 5; the paper does not compare the shock's Mach number with it.
- Along the field: ~9 % of the debris density has u_dx > 0.9 u_d = 4.5 v_A (Sect. II).
- Linear analysis in a window of the 3-D run (Sect. III): mean n_d/n_e = 0.023, M_A = 5.66.
- 2-D beams (Sect. IV): local linear growth rates gamma_RHI ~ 0.12 omega_ci and gamma_NRI ~ 0.15 omega_ci.
  Beam length 50 d_i: interaction time tau ~ Lambda_x/v_dx ~ 7/omega_ci, no saturation, only "a weak wave
  train" (Fig. 7). Beam length 500 d_i (tau ~ 70/omega_ci by the same estimate): the instability saturates
  and the beam slows (Fig. 8). "This type of saturated coupling is required to generate a parallel shock
  front"; this "places lower bounds on the total amount of ion debris required to generate a non-linear
  quasi-parallel shock".
- Even then: "an unambiguous identification of a shock (i.e., by checking whether or not the disturbance
  obeys the Rankine-Hugoniot jump conditions) is not possible", and "the apparent lack of a reflected ion
  population" is "evidence against shock formation" (end of Sect. IV). "the coupling threshold for parallel
  shock formation may be rather steep" (Sect. V).

## Assessment
- For this topic: direct prior work on the drive, not on the electrons. It is the same set-up as ours (a
  debris cloud expanding into a uniform field, in 3-D) at ion scales. It gives one working point for strong
  coupling across the field: M_A = 5, rho_d/R_m = 0.375, R_m = 40 d_i, cavity maximum at t omega_ci ~ 25.
- It warns about the polar regions: with a debris mass of this size the parallel part did not become a
  shock. A front of ours may then be a shock near the equator only, and theta_Bn "from 90 deg to 0" would
  not all be shock. A run has to show a shock at each latitude before it reads electron spectra there.
  (This reading is ours.)
- Limits: fluid electrons, so nothing on electron acceleration; the perpendicular shock is shown but not
  measured (no Mach number, no jump conditions); grid of 1 d_i.
