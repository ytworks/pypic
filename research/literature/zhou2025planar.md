# Zhou, Chen, Dong, Wang, Zou, Walsh & Tóth (2025): Planar Collisionless Shock Simulations with the Semi-implicit Particle-in-cell Model FLEKS

- Ref: ApJS 282, 12; DOI 10.3847/1538-4365/ae1ef6; arXiv:2506.08384. Peer reviewed. Read in the arXiv
  version.
- Read: 2026-11-14, arXiv PDF text: abstract, set-up, the 1-D and 2-D quasi-parallel results, conclusions;
  topics: qpar-electron-conditions

## Claim
A semi-implicit, energy-conserving PIC code (FLEKS, the GL-ECSIM scheme) reproduces the known structure of
planar collisionless shocks with cells of 0.5-4 electron skin depths. Quasi-parallel shocks need at least
two spatial dimensions.

## Method
- Semi-implicit PIC; the shock is set up in its rest frame by relaxing the MHD jump states.
- m_i/m_e = 25; 100 and 400 in 1-D. Δx = 0.5-4 d_e. 400 particles per cell, noise stated as ~5 %.
- Each run lasts 20 ω_ci⁻¹.

## Numbers to compare against
- Whistler precursors appear at m_i/m_e = 100 and 400, not at 25.
- Run length 20 ω_ci⁻¹; our runs last 50-300 ω_ci⁻¹.

## Assessment
- Stated on dimension: "Unlike the 1D simulations, variations in B_x are now present"; "at least two
  spatial dimensions are critical for … capturing the complex dynamics of quasi-parallel shocks, including
  surface rippling, shocklets, SLAMS, magnetic reconnection and jets".
- No nonthermal electrons are studied; the runs are too short for a tail.
- For this topic: (1) a second launching method, relaxation in the shock rest frame, against our wall
  piston; (2) evidence that a reduced mass ratio of 25 removes the whistler precursor; (3) a statement from
  a method paper that a 1-D quasi-parallel run misses rippling, SLAMS in their 2-D form, reconnection and
  jets. Our quasi-1-D results carry that limit; run X is the only check of it.
- The implicit scheme allows coarse cells; whether a tail of electrons survives cells of several d_e is not
  shown there.
