# Gupta, Caprioli & Spitkovsky (2023): Return Currents in Collisionless Shocks

- Ref: arXiv:2312.13365 (v1, 20 Dec 2023). The arXiv record has no journal reference and no DOI; its comment
  field says "Submitted to ApJ". A published version was not verified here; treat it as a preprint until the
  record is checked.
- Read: 2026-11-13, abstract; topics: spherical-shock-injection
- ABSTRACT ONLY (arXiv abstract). Nothing below goes beyond the abstract.

## Claim
The current of particles that a quasi-parallel shock sends upstream is neutralized by the background plasma
through "distinct processes in different regions: the far upstream, the shock precursor, and the shock
foot". Far upstream, nonthermal protons carry the current, "drive electrostatic modes and produce
supra-thermal electrons that move towards upstream infinity". In the precursor, "electromagnetic streaming
instabilities" contribute to the thermalization of these electrons. At the foot, the electrons meet
shock-reflected protons, "resulting in a two-stream type instability".

## Method
- "kinetic particle-in-cell simulations" of "quasi-parallel non-relativistic electron-proton shocks", with
  "particle tracking and controlled simulations". Dimensionality, Mach number and mass ratio: not in the
  abstract.

## Numbers to compare against
- No numbers in the abstract.
- "the instability at the shock foot can make the effective thermal speed of electrons comparable to the
  drift speed of the reflected protons".

## Assessment
- For this topic: it concerns the quasi-parallel part of a front (the polar regions of our sphere). It gives
  a chain from the proton current to supra-thermal electrons that needs accelerated or reflected protons
  first. The abstract presents it as important for "the injection of supra-thermal electrons into diffusive
  shock acceleration".
- The foot result (electron thermal speed near the drift speed of the reflected protons) is a heating level
  that can be checked in a run, once reflected protons are present there.
- Later works of the same authors are `gupta2024electron` and `gupta2025speed`.
