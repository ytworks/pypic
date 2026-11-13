# Torralba Paz, Bohdan & Niemiec (2023): Prediction and Anomaly Detection of accelerated particles in PIC simulations using neural networks

- Ref: Proceedings of Science, PoS(ICRC2023)341 (from the comment field of the arXiv record; the record has
  no DOI); arXiv:2308.15835 (v1). A conference proceeding. The arXiv record gives the first author as
  "Paz, Gabriel Torralba"; the citekey uses the family name "Torralba Paz".
- Read: 2026-11-13, abstract; topics: spherical-shock-injection
- ABSTRACT ONLY (arXiv abstract). Nothing below goes beyond the abstract.

## Claim
Neural networks can replace the "by-eye inspection of particle trajectories" in PIC simulations, which
"includes a high level of bias and uncertainty in pinpointing specific acceleration mechanisms". Regression
and classification "predict the final energies of particles with high accuracy", and anomaly detection "is
able to discern between energetic and non-energetic particles", "regardless of how noisy and imbalanced the
datasets are".

## Method
- A Convolutional Neural Network for classification, regression and anomaly detection.
- Test data: "252,000 electrons which have been traced in a PIC simulation of a non-relativistic high Mach
  number perpendicular shock", in which the "two-stream electrostatic Buneman instability" pre-accelerates
  "a portion of electrons to nonthermal energies".

## Numbers to compare against
- 252,000 traced electrons. No accuracy figure in the abstract.

## Assessment
- For this topic: a method only; no physics result on injection. It is relevant to the analysis: this topic
  has to assign a mechanism to traced electrons, and the abstract names the weakness of doing that by eye.
- The abstract says the networks predict the final energy and separate energetic from non-energetic
  particles. It does not say that they identify the mechanism.
