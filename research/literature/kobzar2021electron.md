# Kobzar, Niemiec, Amano, Hoshino, Matsukiyo, Matsumoto & Pohl (2021): Electron Acceleration at Rippled Low-Mach-number Shocks in High-beta Collisionless Cosmic Plasmas

- Ref: ApJ 919, 97 (2021) (volume and article number from the Crossref record of the DOI);
  DOI 10.3847/1538-4357/ac1107; arXiv:2107.00508 (v1, "to be published in ApJ"); peer-reviewed.
- Read: 2026-11-13, abstract; topics: spherical-shock-injection
- ABSTRACT ONLY (arXiv abstract). Nothing below goes beyond the abstract.

## Claim
At low Mach number shocks in high-beta plasmas, "the electron acceleration rate increases considerably when
the rippling modes appear". "The main acceleration mechanism is stochastic shock-drift acceleration, in
which electrons are confined at the shock by pitch-angle scattering off turbulence and gain energy from the
motional electric field." Multi-scale turbulence "at the shock transition and the region immediately behind
the main shock overshoot is essential for electron energization".

## Method
- "large-scale fully-kinetic two-dimensional particle-in-cell simulations", "in application to merger shocks
  in galaxy clusters". Mach number, beta, obliquity and mass ratio: not in the abstract.

## Numbers to compare against
- Downstream electron spectrum: a power law with index p ~ 2.5, "in agreement with observations"; the
  abstract says this is shown "for the first time".
- "Wide-energy non-thermal electron distributions are formed both upstream and downstream of the shock."
- "The maximum energy of the electrons is sufficient for their injection into diffusive shock acceleration";
  the value is not in the abstract.

## Assessment
- For this topic: stochastic shock drift acceleration (theory in `katou2019theory`; see also
  `amano2022nonthermal`) is a main candidate at the quasi-perpendicular part of a front. This work adds that
  the rate depends on rippling of the shock surface, an ion-scale corrugation of the front.
- A curved front is not a rippled front. But the result shows that the shape of the surface at ion scales
  changes the electron acceleration rate. On a small sphere the curvature and the rippling scales are not
  far apart, so a run of this topic has to tell them apart. (This remark is ours.)
- Regime: low Mach number and high beta. Values are not in the abstract, so the distance to our parameters
  cannot be stated from it.
