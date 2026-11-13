# Le, Winske, Stanier, Daughton, Cowee, Wetherton & Guo (2021): Astrophysical explosions revisited: collisionless coupling of debris to magnetized plasma

- Ref: J. Geophys. Res. Space Physics 126, e2021JA029125 (2021) (volume and article number from the Crossref
  record of the DOI); DOI 10.1029/2021JA029125; arXiv:2109.00583 (v1, "To appear in JGR: Space Physics");
  peer-reviewed.
- Read: 2026-11-13, skimmed (abstract, set-up, conclusions); topics: spherical-shock-injection
- Source: the arXiv full text (HTML of v1). Read: abstract, Sect. 2 (energy balance), Sect. 3 (set-up),
  Sects. 4-5 (cavity size, free-streaming fraction, by figure), Sects. 6-8 (flute modes, 3-D remarks,
  summary).

## Claim
For a debris cloud that expands into a magnetized plasma, simple scaling laws in the total debris mass and
the explosion speed give the maximal size of the debris cloud and the fraction of debris that streams freely
along the field. The laws hold while the debris gyroradius rho_d is smaller than the equal-mass radius R_m.
For rho_d > R_m the debris is "decoupled": single-particle orbits, not a collective cavity, limit its
motion. The fraction that escapes with little coupling grows with the debris gyroradius.

## Method
- 2-D Hybrid-VPIC (kinetic ions, fluid electrons), field in the plane; some runs with the field out of the
  plane for flute modes. No kinetic electrons.
- Scan (Sects. 3-4): M_A = V_d/v_A = 0.5, 1, 2, 5, 10; R_m = 10, 40, 150 d_i; m_d = 3 m_b, equal charge;
  beta_i = beta_e = 0.1. Test cases with a heavy minority species, m_h = 15 m_b.
- R_m is "the radius that contains a mass of ambient plasma equal to the total debris mass" (Sect. 2).

## Numbers to compare against
Conditions for coupling and for a shock:
- Stopping (Sect. 2, Eq. 1): for M_A > 1 the debris is stopped mainly by the swept-up background mass and
  reaches about R_m. For M_A < 1 the field pressure stops it at the magnetic radius R_B, with
  R_B/R_m = M_A^(2/D) in D dimensions.
- Fit of the maximal cloud size across the field (Eq. 4, Fig. 3):
  Delta_max/R_m ~ (1 + c0 M_A) / (1 + c1/M_A^2)^(1/D), with D = 2, c0 = 0.03, c1 = 2 for m_d/m_b = 3.
  "Eq. 4 only applies for systems with rho_d < R_m".
- Decoupled regime: "roughly when rho_d > R_m" (Sect. 4). In the scan these are the runs with R_m = 10 d_i
  at M_A = 5 and 10 (Fig. 3). Example with R_m = 40 d_i and M_A = 10 (Fig. 4): the main species, rho_d = 30
  d_i < R_m, couples; the heavy species, rho_h = 150 d_i > R_m, partly gyrates freely out to about rho_h.
  Here the gyroradius is taken with the gyrofrequency of the debris ion.
- Free streaming along the field (Sect. 5): ions within an angle theta* ~ asin(1.4 rho_d/R_m) of the field
  escape (after Bashurin et al. 1983). Fit for M_A > 2 (Eq. 6, Fig. 6): F_90 ~ C_90 sigma(x), with
  sigma = (x / sqrt(1 + x^2))^1.1, x = rho_d/Delta_max, C_90 ~ 0.29 in 2-D and 0.1 in 3-D (Sect. 7).
  In all runs some ions stream "nearly freely" for "at least several hundred ion inertial lengths".
- Flute modes move debris across the field by "at most a few debris gyroradii"; weaker in 3-D (Sect. 7).
- Shock: only mentioned. "The magnetic field even away from the debris is perturbed by a magnetosonic shock
  that is launched by the expansion" (Sect. 3, on Fig. 1); the cavity collapses "as energy is radiated away,
  possibly in the form of a shock" (Sect. 2). "Details of the streaming ions and expanding shock waves will
  be treated elsewhere" (Sect. 4). NOT stated: the Mach number of the shock, its formation time in ion
  gyro-times, any criterion for its formation, the shock's Mach number against that of the debris.

## Assessment
- For this topic: direct prior work on the drive. It gives the design rule for a debris-driven front: keep
  rho_d < R_m for collective coupling; then the piston stops near R_m (for M_A > 1), with the correction of
  Eq. 4. A run of ours can be placed in its (M_A, R_m/d_i) plane; our 3-D geometry has D = 3 in the energy
  balance, while the constants c0, c1 were fitted in 2-D only.
- It does not answer whether a shock forms, or how strong it is. The scaling laws are for the debris cloud
  and the cavity. A claim that our front is a shock cannot lean on this paper; `keenan2022hybrid` shows a
  perpendicular shock in 3-D for one parameter set, also without measuring it.
- The free-streaming fraction grows with rho_d/R_m. Debris that streams along the field ahead of the front
  is a source of energetic ions, and of beam-driven waves, in the polar regions that is not due to a shock.
  An analysis of electrons there has to separate the two. (This remark is ours.)
