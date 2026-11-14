# Orusa & Valenzuela-Villaseca (2025): Criteria for ion acceleration in laboratory magnetized quasi-perpendicular collisionless shocks: When are 2D simulations enough?

- Ref: Phys. Plasmas 32 (5), 052901; DOI 10.1063/5.0269035; arXiv:2503.00163. Peer reviewed. Read in the
  arXiv version.
- Read: 2026-11-14, arXiv PDF text: abstract, set-up, results of runs A-D, the criteria, the survey of
  experiments, conclusions; topics: qpar-electron-conditions

## Claim
At perpendicular shocks, ions are accelerated in 3-D hybrid runs where the matching 2-D runs show none. The
difference depends on the regime: for M_A < 19 and M_s < 7 "the particle spectra from 2D and 3D simulations
are indistinguishable". Present laboratory shocks lie in that regime, "implying that 1D and 2D simulations
should be enough".

## Method
- dHybridR, non-relativistic regime, θ_Bn = 90° only, 2-D and 3-D, 10 cells per d_i, 8 (3-D) or 4 (2-D)
  particles per cell, adiabatic electrons, t = 10 ω_c⁻¹.
- 2-D runs repeated with 64 and 121 particles per cell.

## Numbers to compare against
- Efficiency at 10 ω_c⁻¹ in 3-D: run A (M_A 25, β 2, M_s 19) 0.3 %; B (25, 5, 13) 0.2 %; C (28, 18, 7)
  0.05 %; D (19, 2, 15) 0.04 %. "No accelerated particles are found in the corresponding 2D simulations".
- Criteria for ion acceleration: M_A ≳ 25, β ≲ 5, duration ≳ 5-10 ω_c⁻¹.
- Experiments surveyed: OMEGA (Schaeffer 2019: ~750 km/s, M_A 15, M_s 15, ~4 proton gyrations); NIF (Fiuza
  2020: non-thermal electrons to ~500 keV, more than 100 times thermal).

## Assessment
- The absence of a tail in 2-D is "not an artifact of particle statistics" (more particles per cell change
  nothing); the stated reason: "in 2D, magnetic field lines effectively act as walls, significantly reducing
  the probability of return".
- Quasi-parallel shocks (θ ≤ 60°) "develop more gradually … leading to a slower shock formation mediated by
  multi-scale plasma processes"; full modelling of the experiments is "computationally prohibitively
  expensive".
- Ions only; perpendicular only. Quasi-perpendicular shocks in 2-D PIC and hybrid runs show "a preference
  for electron rather than ion injection".
- For this topic: a worked example of the question "when does a lower dimension give the same answer", with
  a control for statistics. The result is specific to perpendicular shocks, where field lines block return
  in 2-D; it does not transfer to quasi-parallel shocks. What transfers is the form of the test: the same
  run in two dimensions with the particle number varied (our run X against P0, with null test N2), and a
  statement of the regime in which the lower dimension suffices. It is also the entry to the laboratory
  shocks (OMEGA, NIF) as a boundary area.
