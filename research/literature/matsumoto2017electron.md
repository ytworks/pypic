# Matsumoto, Amano, Kato & Hoshino (2017): Electron surfing and drift accelerations in a Weibel-dominated high-Mach-number shock

- Ref: Phys. Rev. Lett. 119, 105101 (2017); DOI 10.1103/PhysRevLett.119.105101; arXiv:1709.03673 (v1);
  peer-reviewed (PRL)
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
In a 3-D quasi-perpendicular high-M_A shock, coherent Buneman waves at the leading edge and ion-Weibel turbulence in
the foot coexist in separate regions. Electrons are first surfing-accelerated (SSA); the Weibel turbulence then
scatters them in pitch angle and confines them in the transition region while they drift against the motional
electric field ("stochastic SDA"). The result is a power-law tail that neither 2-D geometry produces: the
in-plane run lacks SSA, the out-of-plane run lacks the turbulence.

## Method
- 3-D PIC (SIMD-optimised, hybrid parallel, quadratic particle weighting; code not named), injection method,
  downstream rest frame; energies quoted in the shock rest frame.
- M/m = 64; beta = 1 (beta_i = beta_e = 0.5); theta_Bn = 74.3 deg (B_0 in the x-z plane, E_0 along y); M_s ~ 22.8,
  M_A ~ 20.8 and V_up/c ~ 0.26 in the shock rest frame; subluminal, arctan(c/V_up) > theta_Bn.
  omega_pe/Omega_ce not stated (~10 follows from V_up, M_A and M/m; my arithmetic).
- c/omega_pe = 20 cells; dt = 0.025 omega_pe^-1; L_y = L_z = 4.8 c/omega_pi, L_x grows; 8800 x 768 x 768 cells at the
  end; 20 particles per cell per species; ~1e12 particles (K computer).
- Structure shown at Omega_gi T = 7.62; ~1e7 tracer electrons followed over 6.8 < T < 8.8.
- Companion 2-D runs with in-plane and out-of-plane B_0 "under the same upstream conditions"; superluminal runs
  are mentioned without parameters.

## Numbers to compare against
- Buneman waves with |E| > B_0 at the leading edge (8 < x < 10 c/omega_pi); fronts oblique in x-y; no structure
  along z (Fig. 1a, 2a).
- Foot (0 < x < 6): Weibel rib structures; B components amplified up to ~20 times upstream after compression
  (Fig. 1b). delta B/B ~ 1 near the front in 3-D, ~0.3 in the 2-D out-of-plane run (Fig. 4a, c).
- Stage 1, SSA at 7 < x < 8.5: d(gamma) ~ 1, gain perpendicular to B, mu x ~15 (Fig. 2c-d, 3c), during 6.9 < T < 7.3.
- Stage 2 in 0 <= x <= 4: mean energy of the top-1000 electrons follows d(gamma) = e E_0y dy/mc^2 (Fig. 3b), as in
  SDA, but mu keeps rising and parallel and perpendicular momentum both grow: not adiabatic (Fig. 2d, 3c).
- Tracers: thermal peak settles at gamma-1 ~ 0.2 by T ~ 8.1; tail cut-off gamma-1 ~ 5 at T ~ 8.8 (Fig. 3d).
- Electrons with gamma > 5 diffuse upstream, beyond x > 10 (Fig. 4a).
- Downstream spectrum (-5.0 < x < -4.5): power law of index 3.5 in the energy distribution in 3-D; Maxwellian in
  2-D in-plane; small nonthermal part in 2-D out-of-plane (Fig. 4d).

## Assessment
- Shown by simulation: one 3-D run and two 2-D runs; one traced orbit and statistics of the 1000 most energetic
  against 1000 randomly chosen tracers.
- Not given: injected fraction, energy fraction, maximum energy of the full (non-tracer) spectrum as a number.
- Stated without data: superluminal cases accelerate little (as for magnetised relativistic shocks, Sironi et
  al. 2013); the mechanism should work over a "wide range" of obliquity if c/V_up >> 1; reconnection in Weibel
  turbulence (Matsumoto et al. 2015) was absent because M_A was too low.
- Limits I note: M/m = 64; V_up/c = 0.26, so the subluminal limit is arctan(1/0.26) = 75.4 deg (my arithmetic)
  and theta_Bn = 74.3 deg lies just inside it; one obliquity; transverse box 4.8 c/omega_pi; tracers followed
  for 2 Omega_gi^-1.
- Left open by the authors: longer evolution, self-excited upstream waves and entry into the DSA cycle.
- Relation to other works: the weak Buneman mode with in-plane field agrees with Riquelme & Spitkovsky (2011)
  and Wieland et al. (2016); the limited Weibel growth with out-of-plane field is attributed to Wieland et al.
  (2016). The result contradicts the "naive" expectation that Weibel turbulence would destroy SSA. My reading:
  the second stage is diffusive here, where Amano & Hoshino (2007) and Matsumoto et al. (2012) had it adiabatic.
