# Caprioli & Spitkovsky (2014): Simulations of ion acceleration at non-relativistic shocks. I. Acceleration efficiency

- Ref: Astrophys. J. 783, 91 (2014) (as printed in Pohl et al. 2020, ref. 58, and Marcowith et al. 2016,
  ref. 28); DOI 10.1088/0004-637X/783/2/91; arXiv:1310.2943; peer-reviewed (arXiv comment: accepted).
- Read: 2026-11-12, introduction, set-up, Sect. 3-5, 7-9 (efficiency against obliquity and Mach number,
  field amplification, SN 1006, drift acceleration, 3-D runs, conclusions); topics: spherical-shock-injection

## Claim
With a laminar, ordered upstream field, ions are injected into diffusive shock acceleration, and drive
upstream magnetic turbulence, only at quasi-parallel shocks. Above about 45 deg obliquity they gain a
factor of a few by shock drift, the maximum energy stops growing, and no upstream turbulence appears.

## Method
- Hybrid code dHybrid (kinetic ions, fluid electrons with a polytropic law), non-relativistic, 2-D and
  3-D; initial field in the simulation plane. No electron physics; m_i/m_e and v_sh/c do not enter.
- Shock by a flow against a reflecting wall. M = v_sh/v_A with v_sh the upstream speed in the downstream
  (wall) frame; in the shock frame M is 5/4 larger for compression 4 (Eq. 1). M_A ~ M_s (ion thermal
  speed = v_A, T_e = T_i; beta is not quoted).
- Reference run: M = 20, parallel, box 1e5 x 1e2 (c/omega_p)^2, 2 cells per ion skin depth, 4 particles
  per cell, step 5e-3 omega_c^-1, run to 2500 omega_c^-1 (Sect. 2).
- Survey: obliquity 0, 20, 30, 45, 50, 60, 80 deg; M = 5, 10, 30, 50; box 40000 x 500 (c/omega_p)^2;
  step 0.01/M omega_c^-1; run to 200 omega_c^-1 (Sect. 4).
- 3-D: M = 6; 0, 45, 80 deg; box 2000 x 200 x 200 (c/omega_p)^3; 8 particles per cell (Sect. 8).

## Numbers to compare against
- Efficiency = share of post-shock energy density in ions above 10 E_sh, E_sh = m v_sh^2 / 2, at
  200 omega_c^-1 (Fig. 3): 5-15 % for low obliquity; nearly flat below about 40 deg, most clearly for
  M = 5 and 10; drops above about 45 deg; tends to zero above about 60 deg at every M; highest for
  M = 50, parallel. Abstract and conclusions quote 10-20 % for strong quasi-parallel shocks.
- 3-D, M = 6: about 12 %, 3 %, 1 % at 0, 45, 80 deg (Sect. 8).
- Parallel shocks: injected fraction 1e-3 to 1e-4 (Sect. 3), about 1e-3 (Sect. 9); injection momentum
  3-4 times the downstream thermal momentum (Eq. 5: 3-3.5); spectrum p^-4, i.e. E^-1.5 (Fig. 1).
- Above 45 deg: ions stay below about 10 E_sh (Sect. 7, 9); drift-heated component about 2.2 T_d
  for compression 4 (Sect. 7, after Park et al. 2012).
- Upstream field, M = 30, quasi-parallel: dB/B_0 = 3-4 ahead of the shock, post-shock 6-10 B_0 on
  average; quasi-perpendicular: only compression, about 4 B_0 (Fig. 4). M = 50: up to dB/B_0 ~ 6 (Sect. 9).
- 3-D out-of-plane field upstream: about B_0 at 0 deg, below 0.1 B_0 at 80 deg (Fig. 14, Sect. 8).
- Self-generated transverse field at a parallel shock is not expected to raise the upstream obliquity
  above about 45 deg (Sect. 7; an argument, not a measurement).

## Assessment
- Shown by simulation: the drop of ion efficiency and of self-generated turbulence with obliquity, for
  M = 5-50 in 2-D and M = 6 in 3-D.
- Reasoning, not measured: at high obliquity returning ions reach only about one gyroradius upstream, so
  waves have about one cyclotron time to grow (Sect. 7). The injection mechanism is deferred.
- Stated limit: valid for laminar inflow and ordered field; pre-existing turbulence of dB/B ~ 1 on scales
  above the gyroradius of ions with v ~ v_sh may allow injection at oblique shocks (Sect. 7, 9).
  Efficiencies at 200 omega_c^-1 may not be saturated (Sect. 4).
- Whole-remnant geometry (Sect. 5.1): in SN 1006 the bright synchrotron caps are where the shock is
  quasi-parallel (Reynoso et al. 2013); the radio polarization is low there and up to 70 % elsewhere,
  as expected from turbulent against merely compressed fields. The authors add that brightness also
  depends on how well electrons are accelerated, which they do not model.
- Electrons (Sect. 5.1, statements only): observed to be accelerated at oblique shocks; whistler waves
  are excited preferentially at low-Mach-number oblique shocks; whether quasi-parallel shocks accelerate
  electrons is "under investigation".
- Disagreements it discusses: Gargate & Spitkovsky (2012) found the maximum at 15 deg and M = 6, here
  blamed on short boxes and fixed time step (Sect. 4); efficient ion acceleration at perpendicular shocks
  (Baring et al. 1995; Giacalone & Ellison 2000) relied on prescribed or seeded scattering (Sect. 8).
- Doubt: injection energy is 4-5 E_sh in the text of Sect. 3 but 3-4 E_sh in the caption of Fig. 2.
