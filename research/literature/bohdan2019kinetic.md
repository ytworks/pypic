# Bohdan, Niemiec, Pohl, Matsumoto, Amano & Hoshino (2019): Kinetic simulations of nonrelativistic perpendicular shocks of young supernova remnants. I. Electron shock-surfing acceleration

- Ref: ApJ 878, 5 (2019) (volume/article as cited in Papers II and III; meta.txt has no journal_ref);
  DOI 10.3847/1538-4357/ab1b6d; arXiv:1904.13153 (v1, "Submitted to ApJ", read here); peer-reviewed
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
The energy an electron gains while surfing at the leading edge of the foot always comes from the electrostatic
field of a Buneman wave it travels with. The wave strength is set by how far MA lies above the trapping threshold
(1+α)(mi/me)^(2/3), not by MA or mi/me alone. In-plane 2-D runs cannot reach the surfing efficiency of out-of-plane
(or 3-D) runs even above the φ-corrected threshold, because the in-plane relative speed of reflected ions and
incoming electrons, and with it the wave phase speed and the number of electrons picked up, is lower.

## Method
- PIC, 2D3V, modified TRISTAN with particle tracing; flow-flow set-up of Bohdan et al. (2017): beams at ±0.2c,
  equal density, temperatures differing by 1000 (βe = 5·10⁻⁴ left, runs *1; βe = 0.5 right, runs *2).
- θBn = 90°. Table 1, runs (φ; mi/me; ωpe/Ωe; MA; Ly/λsi): A (0°; 50; 12; 22.6; 10.9), B (0°; 100; 12; 31.8; 24),
  C (0°; 100; 17.3; 46; 12), D (0°; 200; 8.5; 32; 11.9), E (0°; 200; 12; 44.9; 11.9), F (0°; 400; 12; 68.7; 8.2),
  G (90°; 100; 12; 35.5; 12). Ms (*1/*2): A 949/30, B 1342/42.4, C 1941/61.4, D 1342/42, E 1898/60, F 2904/91.8,
  G 1369/43.3. Mach numbers in the upstream frame.
- Shock speed in the upstream frame 0.263c (A-F, Γ = 5/3) and 0.294c (G, Γ = 2) (Sect. 2).
- λse = 20Δ; 20 particles per cell per species; Lx up to ≈ 280 λsi; run time (6-8) Ωi⁻¹; time step 1/40 ωpe⁻¹.
- Purpose of the runs: A, B, E, F exceed the trapping threshold by a similar margin at four mass ratios; C also
  exceeds the φ-corrected threshold; D lies below the threshold (Sect. 2, Fig. 3).
- Appendix A: non-relativistic equations of motion of an electron in a Buneman wave plus B0, for φ = 90° and 0°.

## Numbers to compare against
- Trapping threshold MA ≥ (1+α)(mi/me)^(2/3) (Eq. 1): 16, 26, 41, 65 for mi/me = 50, 100, 200, 400 at α = 0.2;
  φ-corrected (Eq. 2): 22.4, 36, 58, 92 (α = 0.2) and 28, 46, 72, 115 (α = 0.5) (Table 1).
- Proposed extra factor 1.5 for the in-plane case (Eq. 3): MA ≃ 60.3 needed at mi/me = 100; not tested (Sect. 3.1).
- Table 2, βe = 0.5 shocks, peak |EES|/(cB0) and fraction pre-accelerated to γ−1 > 0.1: A2 1.1, 0.43 %; B2 1.3,
  0.46 %; C2 2.3, 0.6 %; D2 0.4, 0.34 %; E2 1.3, 0.49 %; F2 1.1, 0.44 %; G2 2.7, 6.8 %; G2* (another reformation
  phase, same field as C2) 2.3, 2.7 %.
- A2, B2, E2, F2 give statistically identical foot spectra, fraction ~0.45 % (Fig. 5). At fixed mi/me a higher MA
  gives 30-40 % more pre-accelerated electrons, not the factor 10 to the out-of-plane case (Sect. 3.2).
- Largest Lorentz factor after surfing ≈ 3-4, comparable in C2 (0°) and G2* (90°) (Sect. 3.2, Fig. 6).
- Buneman wavelength ≈ 1.9 λse (in plane) and 3.3 λse (out of plane) (Sect. 3.1); phase speed to match 0.1c and
  0.4c (Sect. 3.2); relative speed of electrons and reflected ions up to ~0.6c out of plane (Sect. 3.1).
- Trapping times ≈ 13 ωpe⁻¹ (0°) and 11 ωpe⁻¹ (90°) (Eq. 4); mean gains Δγ ≈ 0.18 and 0.42 (Eq. 5); estimate
  Δε ≈ Δv² ttr ωpe me (me/mi)^(1/6) (Eq. 6). Potential-well depth Δγ ≈ 0.05 and 0.17 (Sect. 3.2).
- At βe = 0.5 the Buneman area is 20-30 % larger and the intensity 20-50 % lower than at βe = 5·10⁻⁴ (Sect. 3.1).

## Assessment
- Shown by simulation: wave amplitude and pre-accelerated fraction depend on MA relative to the threshold and not
  on mi/me over 50-400; run D (below threshold) has three times weaker waves than B at the same MA; run C meets
  Eq. 2 but stays 20 % below G in field and a factor ~10 below in fraction; at equal field (G2*) the out-of-plane
  run still pre-accelerates four times more electrons.
- Shown with traced electrons and the analytic model: in the upstream flow frame the energy gain is from the wave
  field alone (Eqs. A5, A8); the interaction time is limited by the intermittency of the waves (Fig. 8, Eq. 4).
- Stated or conjectured: the factor 1.5 of Eq. 3; that out-of-plane 2-D reproduces the 3-D surfing, taken from
  Matsumoto et al. (2017), which is a 3-D run of an oblique subluminal shock, not a perpendicular one.
- Limits stated: only the first stage in the Buneman zone is analysed (ramp and overshoot in Papers II, III);
  strictly perpendicular; 2-D; φ = 45° not rerun because it behaved like 0° in Bohdan et al. (2017).
- Not addressed: dependence on shock speed (all runs 0.263c or 0.294c); numbers are for βe = 0.5 shocks.
- Doubts: Ms in Table 1 differs from the values for the same runs in Papers II and III and in Bohdan et al. 2017
  (B: 1342/42.4 here, 1550/49 there). ωpe/Ωe of run F is 12 here and 13 in Paper III; 13 is the value that
  reproduces MA = 68.7 at 0.263c (my check). The summary says the efficiency is "proportional to MA" at fixed mass
  ratio, which rests on two pairs of runs.
- Introduction: observations do not fix the field geometry at SNR shocks; models of the same remnant favour
  quasi-perpendicular or quasi-parallel fields (Sect. 1).
