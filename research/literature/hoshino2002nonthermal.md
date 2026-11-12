# Hoshino & Shimada (2002): Nonthermal Electrons at High Mach Number Shocks: Electron Shock Surfing Acceleration

- Ref: ApJ 572, 880–887 (2002); DOI 10.1086/340454; arXiv:astro-ph/0203073
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
In a high-Mach-number perpendicular shock, the Buneman instability between reflected ions and incoming
electrons builds large-amplitude electrostatic solitary waves (ESWs, electron holes) in the foot. Electrons
trapped in an ESW drift along the shock surface and gain energy from the motional electric field
("electron shock surfing"). Above a critical Alfvén Mach number the ESW force always beats the Lorentz
force, so trapping, and acceleration, would be unlimited in a stable, planar shock.

## Method
- 1D (x; three velocity components, y also integrated for orbits) relativistic EM PIC (Hoshino et al. 1992
  code). Plasma injected at x = 0 against a reflecting wall; shock runs in −x, downstream at rest. Perpendicular
  shock, B along z (Sect. 2).
- One run: βe = βi = 0.01, ωpe/ωce = 19, mi/me = 20, L = 300 c/ωpe, v0 = 0.25c (≈ 3c/8 in shock frame),
  VA/c = 1.2e-2, MA ≈ 32; cell ≈ downstream Debye length, ≈ 1800 particles/species/cell downstream;
  run to ωpe t = 1815 (Sect. 2).
- Diagnostics: phase space and fields (Fig. 1); downstream spectrum summed over 538 snapshots with Maxwellian
  fit and one-count level (Fig. 2); two tracked electrons: Ex stack plot and |P| history (Fig. 3), first
  adiabatic invariant (Fig. 4), x–y orbit (Fig. 5).
- Analytic model (Sect. 4): ESW amplitude from a fraction α of the electron drift energy (Eqs. 2–4,
  α ~ 1/4 or (me/mi)^(1/3)), force balance eE_esw = e vy B0/c (Eqs. 5–7), unlimited-trapping threshold (Eq. 8),
  reformation-limited maximum energy (Eq. 9).

## Numbers to compare against
- Nonthermal excess above the Maxwellian for γ/γ0 > 1.8, γ0 ≈ 1.03, MA ≈ 32, mi/me = 20 (Fig. 2, Sect. 2);
  spectrum hits the one-count level near γ/γ0 ≈ 4.5 (Sect. 5).
- Crossing one ESW alone: Δγ/γ0 ≈ 0.26 (E_esw ≈ 20 E0, width ≈ c/ωpe from Fig. 1; Sect. 2).
- Surfing gain of the trapped orbit: ΔY ≈ 60 c/ωpe, Δε ≈ 1.6 me c², γ ≈ 2.6, P/P0 ≈ 9 (Sect. 3, Figs. 3, 5).
- E_esw/E0: 10–30 in the run vs 19–23 from Eq. 4; Earth bow shock observed 10–100 vs 13–23 from Eq. 4
  with VA ≈ 300 km/s (Sect. 4).
- Unlimited trapping for MA > (1/2)(mi/(me α))^(1/2) (Eq. 7), i.e. MA > 43–75 at mi/me = 1830 (Eq. 8).
- Reformation limit ε_max ~ mi c² (v0/c) (Eq. 9); γ_max/γ0 ≈ 7.5 for the run (Sect. 5).
- Cited earlier PIC (Shimada & Hoshino 2000): suprathermal tail at MA ≈ 10, none at MA ≈ 3 (Sect. 2).

## Assessment
- Shown (simulation): ESW chain in the foot and a weak suprathermal tail in one 1D run. For two
  representative electrons, the energy gain coincides with ESW trapping and transverse drift, and the
  magnetic moment is broken only while trapped (Figs. 3–5). The attribution of the whole tail to surfing
  rests on these orbits plus an order-of-magnitude estimate; there is no statistics over tail particles and
  no parameter scan.
- Shown (theory vs data): the Eq. 4 ESW amplitude matches the run and bow-shock wave data within about a factor of 2.
- Stated/conjectured: unlimited acceleration above MA ≈ 43–75 (assumes a stable ESW; not simulated); SNR
  relevance; surfing faster than DSA. Whether surfing gives a power law is left open (Sect. 5).
- Note (my evaluation of Eq. 7): at mi/me = 20 the threshold is MA ≈ 4–5, so the MA ≈ 32 run is nominally in
  the unlimited regime, yet the tracked electron detraps. The paper attributes this to time-variable ESWs
  that slow down to the plasma frame (Sects. 3, 5) and does not discuss the tension explicitly.
- Stated limits (Sect. 5): 1D only (finite transverse ESW size in 2D/3D may limit acceleration, or hopping
  between ESWs may extend it); reduced mi/me = 20; too few particles to see the cutoff (≈ 100× needed);
  reformation at about the ion gyroperiod may detrap electrons; only weakly magnetized ωce < ωpe; result
  depends on ESW speed relative to the plasma.
- Open: 2D shocks, dependence on ωce/ωpe, power-law formation, in-situ tests at high-MA interplanetary shocks.
