# Haggerty & Caprioli (2019): dHybridR: A Hybrid Particle-in-cell Code Including Relativistic Ion Dynamics

- Ref: ApJ 887, 165; DOI 10.3847/1538-4357/ab58c8; arXiv:1909.05255. Peer reviewed. Read in the arXiv
  version.
- Read: 2026-11-14, arXiv PDF text: abstract, the model and its range of validity, the electron closure,
  the heliospheric example, conclusions; topics: qpar-electron-conditions

## Claim
A hybrid code with relativistic ions (kinetic ions, massless fluid electrons) follows ion acceleration at
shocks into the relativistic range. Results obtained at shock speeds of ~0.1 c apply to slower shocks as
long as the scales are separated.

## Method
- Hybrid; electrons are a massless fluid with a prescribed polytropic equation of state. The effective
  index is chosen to give equipartition downstream, because adiabatic electrons give too low an electron
  pressure at high Mach number.
- Valid for non-relativistic bulk flows and n_cr ≪ n_i. The time step scales as 1/(c/v_A), so runs use
  V_sh ~0.1 c.
- A 3-D low-Mach-number run with "the trajectory of a synthetic probe through the simulation box, mimicking
  in-situ spacecraft observations".

## Numbers to compare against
- None taken.

## Assessment
- On the speed of light in the runs: "There are numerous astrophysical systems where c/v_A and c/v_sh are
  considerably larger than the simulations presented in this work; however, as long as there is a clear
  separation of scales between the thermal/Alfvénic speed, the speed of the shock, and the speed of light,
  the underlying physics can be studied fruitfully." This is stated, not tested, and it is stated for ions.
- On electrons: "the electrons physical description and evolution is more ad hoc".
- Link to observations: the synthetic probe "reproduces both qualitatively and quantitatively the shock
  rippling that has been found with in-situ satellite observations (Johlander et al. 2018)".
- For this topic: (1) the argument of scale separation is the usual licence for running at v_sh ~0.1 c. For
  electrons it is weaker: an electron of given p̃ moves at a different u/v_sh in the run than in nature, so
  the dependence on v_sh/c has to be measured (idea 2, null test N4) and cannot rest on this argument. (2) A
  published way of confronting a run with a spacecraft: fly a synthetic probe through the box. (3) Hybrid
  runs carry no electron kinetics; what they give us is the ion-scale turbulence.
