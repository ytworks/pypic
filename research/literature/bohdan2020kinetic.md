# Bohdan, Pohl, Niemiec, Vafin, Matsumoto, Amano & Hoshino (2020): Kinetic simulations of nonrelativistic perpendicular shocks of young supernova remnants. III. Magnetic reconnection

- Ref: ApJ 893, 6 (2020) (volume/article as cited in Bohdan et al. 2021; meta.txt has no journal_ref);
  DOI 10.3847/1538-4357/ab7cd6; arXiv:2003.01879 (v1, "Submitted to ApJ", read here); peer-reviewed
- Read: 2026-11-12, full (the dispersion-relation appendix only skimmed); topics: spherical-shock-injection

## Claim
In 2-D perpendicular shocks with an in-plane field component, the ion Weibel filaments of the foot turn into
current sheets that break up by reconnection in the ramp. The number of reconnection sites rises with mi/me and
MA. Electrons are energised at these sites, and at the highest mass ratio (400) about 80 % of all electrons pass
through one, so reconnection supplies most of the (small, ~0.5 %) downstream nonthermal population.

## Method
- PIC, 2D3V, code THATMPI (from TRISTAN; triangular-shaped-cloud particles, Friedman filter, Vay pusher); flow-flow
  set-up, beams at ±0.2c; θBn = 90°; βe = 5·10⁻⁴ (left, runs *1) and 0.5 (right, runs *2).
- Table 1, runs (φ; mi/me; ωpe/Ωe; MA; Ms cold/warm): A (0°; 50; 12; 22.6; 1104/35), B (0°; 100; 12; 31.8; 1550/49),
  C (0°; 100; 17.3; 46; 2242/71), D (0°; 200; 8.5; 32; 1550/49), E (0°; 200; 12; 44.9; 2191/69),
  F (0°; 400; 13; 68.7; 3353/106), G (45°; 100; 12; 31.8; 1550/49). Mach numbers in the upstream frame.
- Shock speed 0.263c (upstream frame), 0.067c (simulation frame); r = 4; λse = 20Δ; 20 particles per cell per
  species; time step 1/40 ωpe⁻¹ (Sect. 2). Box sizes and run times are not restated (see Paper I).
- Diagnostics: "magnetic vortex" = local maximum of Az; vortex number VN = vortices per λsi of transverse box
  size; electron density inside vortices (Sect. 3.1). Linear ion-Weibel growth rate from the measured foot
  distributions (App. A, Fig. 4). Hand-drawn reconnection regions, 6-11 per warm shock, for the bookkeeping of
  Table 3 (Sect. 3.4, Fig. 6). Traced electrons split into energy parallel and perpendicular to B (Fig. 5).

## Numbers to compare against
- VN, time-averaged (Table 2), cold/warm: A 0.1/0.11; B 0.61/0.68; C 1.53/1.36; D 0.89/0.87; E 1.97/2.34;
  F 5.03/6.4; G (45°) 1.15/1.17. No dependence on βe. VN swings 0.2-1.5 in B2 with the reformation phase (Fig. 2).
- Weibel mode: wavelength ≈ λsi, growth rate ≈ 0.1 ωpi in all in-plane runs (Fig. 4). Field |BW| ≈ 0.1 MA |B0|
  (Eq. 1); guide field in the sheets ≲ 8 |B0| (Sect. 3.1).
- Reconnection condition: foot length rgi exceeds LW ≈ 100 v0/ωpi, i.e. 10⁻² sqrt(mi/me) ωpe/Ωe > 1 (Eqs. 2-3);
  left-hand side 0.8 (A), 1.2 (B, D, G), 1.7 (C, E), 2.6 (F) (Table 2).
- Prediction VN ≈ 12.5 for the real mass ratio at MA ≈ 150, from VN ∝ MA at high mass ratio (Sect. 3.2).
- Time available for reconnection ≈ 0.5 Ωi⁻¹ ≤ 3.5 local ion gyro-times (Eq. 9), cut off by reformation.
- Table 3, warm shocks (kBT/me c²; NTEF; from surfing; from reconnection; fMR): A2 0.09; 0.28 %; 0.15 %; 0.1 %;
  10 %. B2 0.18; 0.55 %; 0.05 %; 0.4 %; 26 %. C2 0.22; 0.36 %; 0.05 %; 0.13 %; 38 %. D2 0.3; 0.7 %; 0.025 %;
  0.55 %; 38 %. E2 0.37; 0.56 %; 0.025 %; 0.5 %; 43 %. F2 0.73; 0.57 %; 0.0015 %; 0.4 %; 79 %.
  fMR = share of all electrons crossing the shock that pass a reconnection region (Eq. 10).
- Reconnection share of the tail saturates near 0.8 at high mass ratio (Sect. 3.4).
- Guiding-centre ordering needs 0.1 (mi/me) vsh/(γe ve) ≫ 1 (Eq. 8): marginal in the runs, violated for energetic
  electrons at vsh ≲ 0.02c whatever the mass ratio (Sect. 3.3).

## Assessment
- Shown by simulation: chains of X-points and islands at φ = 0° and 45° (Fig. 1); VN rises with mi/me and MA and is
  almost twice as large at 45° as at 0° (there many vortices come from turbulence, not filament decay); four kinds
  of trajectory among electrons reaching γ > 3: Ez at an X-point (Speiser orbit), capture in an island ("island
  surfing" plus compression), bouncing between merging islands, and stochastic scattering off vortices (Fig. 5).
- The attribution in Table 3 rests on manually chosen regions scaled up through VN; no error bars are given.
  Surfing plus reconnection always sum to less than NTEF; the rest is put down to chaotic scattering in the ramp
  and overshoot turbulence (Sect. 3.4).
- Stated or conjectured: reconnection dominates at realistic parameters (extrapolated from mi/me ≤ 400 at 0.263c);
  the 3-D scaling with mi/me and ωpe/Ωe "should be similar"; reconnection is more efficient in 3-D (Dahlin et al.
  2015); it may take part in electron reflection at quasi-perpendicular shocks (to be clarified, Sect. 4).
- Limits stated: only warm shocks analysed in detail (cold ones "similar"); reconnection reaches only early stages;
  at mi/me ≲ 100 only some filaments reconnect; no ion heating by reconnection is seen.
- Against other work: the 3-D run of Matsumoto et al. (2017) shows no reconnection, explained by Eq. 3 failing at
  its small mass ratio (Sects. 1, 4). Own cross-reference, not in the paper: Wieland et al. (2016), 45° and
  mi/me = 50, also saw none.
- Doubts: the VN extrapolation is printed with the MA ratio inverted, but the number 12.5 follows from VN ∝ MA;
  ωpe/Ωe of run F is 13 here and 12 in Paper I; MA of E2 is 44.6 in Table 3 and 44.9 in Table 1.
