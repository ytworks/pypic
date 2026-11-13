# Shankar, Li, Dorfman & Fu (2026): Particle-in-Cell Simulation of the Parametric Decay Instability of Alfven Waves with Absorbing Boundary Conditions

- Ref: preprint, arXiv:2606.04932 (submitted 2026-06-03); no journal reference on arXiv; not
  peer-reviewed.
- Read: 2026-11-13, abstract; topics: shock-alfvenic-turbulence

## Claim
Fully kinetic 1-D (k_perp = 0) PIC of the parametric decay of an Alfven wave at very low beta with
absorbing wave boundaries (open system, SMILEI, parameters tied to the LAPD device per the search
snippet): about 92 % of the pump energy goes to the backward daughter wave, 6-7 % to ions and 1-2 % to
electrons. Heating starts only once the instability is well developed, at a rate about twice the linear
growth rate.

## Method
- Fully kinetic PIC, 1-D, k_perp = 0, absorbing boundaries for waves (first such study, per authors).
- beta = 5e-4, pump amplitude delta B/B0 = 0.01 (main case); scans in plasma and wave parameters.
- Left-hand circular polarization (search snippet). Box, mass ratio, resolution: not read.

## Numbers to compare against
- Energy partition at beta = 5e-4, delta B/B0 = 0.01: daughter wave ~92 %, ions ~6-7 %, electrons ~1-2 %.
- Heating rate ~ 2 x linear growth rate (consistent with energy ~ amplitude^2).
- Theory vs simulation growth rates: "qualitative" agreement over the scanned parameters (abstract).

## Assessment
- Far from our regime (beta and amplitude two orders of magnitude lower), but gives the open-boundary
  counterpart of gonzalez2023particle: with absorbing boundaries the electron share is smaller
  (1-2 % vs ~10 % at delta B/B0 = 0.5, beta 0.25 periodic). The two differ in amplitude, beta and
  boundaries at once, so the cause of the difference is not separable from the abstracts.
- Authors expect stronger electron heating with finite k_perp: not shown.
- Preprint; check for a journal version before citing numbers.
