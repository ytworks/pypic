# Caprioli, Pop & Spitkovsky (2015): Simulations and theory of ion injection at non-relativistic collisionless shocks

- Ref: ApJL 798, L28; DOI 10.1088/2041-8205/798/2/L28; arXiv:1409.8291. Peer reviewed. Read in the arXiv
  version.
- Read: 2026-11-14, arXiv PDF text, whole letter; topics: qpar-electron-conditions

## Claim
Ion injection at quasi-parallel shocks is set by the periodic reformation of the shock: ions that arrive
when the barrier is steepest are reflected, gain energy by shock drift, and are injected if they survive
enough cycles. A minimal model with one loss probability per cycle reproduces the injected fraction and the
spectrum of the simulations.

## Method
- 2-D hybrid (dHybrid), 4 particles per cell; ion trajectories followed through the shock.
- Electron pressure polytropic (γ_e ≈ 4.3 at M = 20): "this parametrization of the electron physics has
  little effect on the overall shock structure and on ion injection".

## Numbers to compare against
- ~25 % of the ions are reflected at each reformation; reflection "is a necessary condition for being
  injected, but it is not sufficient".
- Injection velocity v_inj ≳ 2.5-3.5 V_sh; E_inj ≳ 5-10 E_sh. Supra-thermal ions reach E ≲ 6 E_sh by drift.
- Minimal model: N ≈ 2.4 cycles, injected fraction ≈ 3.6 % (≲ 4 %).
- Injection is suppressed for θ ≳ 45-50° (at 50°, N ≈ 3.8 and a fraction ~5 × 10⁻³).
- Quasi-parallel shocks amplify and isotropize the field ahead of them, so the effective inclination is
  ~45°. For M ≲ 4 the overshoot vanishes and the efficiency is lower.

## Assessment
- "In this work we did not account for injection of electrons, and of heavy ions".
- On method: "A characterization of particle injection without free parameters requires a self-consistent
  calculation of the shock structure on microphysical scales, which can be achieved only with kinetic plasma
  simulations". Monte Carlo models reproduce the supra-thermal bridge at Earth's bow shock "but need an a
  priori parametrization of the ion mean free path".
- On real shocks: "At real quasi-parallel shocks (e.g., Burgess et al. 2005), reformation is unlikely
  coherent along the shock surface"; pre-existing turbulence may make injection patchy; "our results are
  still expected to be valid in a time/space-averaged sense, as injection is a local process".
- "the dependence of ion injection on the shock obliquity is not an artifact of 2D simulations".
- For this topic: the model of what a finished injection result looks like: a definition of "injected"
  (a threshold in energy, a number of cycles, a loss probability), a fraction, and a minimal model that
  reproduces the spectrum. No such result exists for electrons at event conditions in what we read. Three
  points bear on our set-up: reformation is coherent across a quasi-1-D box by construction; an effective
  inclination of ~45° means the PSP case (θ 8°) and the θ 30° runs may differ less at the shock than far
  upstream; the loss of the overshoot at M ≲ 4 lies in the range of the 1-au case (run W, M_A 3).
