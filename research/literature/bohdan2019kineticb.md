# Bohdan, Niemiec, Pohl, Matsumoto, Amano & Hoshino (2019): Kinetic simulations of nonrelativistic perpendicular shocks of young supernova remnants. II. Influence of shock-surfing acceleration on downstream electron spectra

- Ref: ApJ 885, 10 (2019) (volume/article as cited in Paper III; meta.txt has no journal_ref);
  DOI 10.3847/1538-4357/ab43cf; arXiv:1909.05294 (v1, header "Accepted September 13, 2019", read here); peer-reviewed
- Read: 2026-11-12, introduction, set-up and conclusions (results read for the numbers); topics: spherical-shock-injection

## Claim
For in-plane 2-D perpendicular shocks the downstream nonthermal electron fraction does not follow the surfing
efficiency in the foot. The tail is made mainly in the turbulent ramp and overshoot (stochastic Fermi-like
scattering, reconnection sites). Surfing gives a fixed energy (<γ−1> ≈ 0.15) to ~0.5 % of the electrons, while the
downstream temperature grows almost linearly with mi/me, so the share of the tail owed to surfing falls with mass
ratio and should be minor at the real mass ratio for a strictly perpendicular shock.

## Method
- Same runs as Paper I, in-plane only (φ = 0°): PIC 2D3V, modified TRISTAN, flow-flow, beams at ±0.2c, θBn = 90°.
- Table 1, runs (mi/me; MA; Ms cold/warm): A (50; 22.6; 1104/35), B (100; 31.8; 1550/49), C (100; 46; 2242/71),
  D (200; 32; 1550/49), E (200; 44.9; 2191/69), F (400; 68.7; 3353/106); βe = 5·10⁻⁴ (*1) and 0.5 (*2).
  Mach numbers in the upstream frame; shock speed 0.263c there; ωpe/Ωe = 8.5-17.3; 20 particles per cell per
  species; λse = 20Δ; run time (6-8) Ωi⁻¹ (Sect. 2). Box sizes: see Paper I.
- Spectra are taken behind the overshoot over ≈ 25 λsi (two reformation cycles). Nonthermal fraction (NTEF) =
  excess over a relativistic Maxwellian fit; with a 5 kBT cut, NTEF = N(E > 5 kBT) − 1.85 % (Sects. 3.1, 3.2).
- About 5·10⁵ traced electrons in run E2; two subsets: "BI" (γ−1 > 0.1 right after the Buneman zone) and "HE"
  (E ≥ 5 kBT downstream) (Sect. 3.2, Fig. 4).

## Numbers to compare against
- Table 2, NTEF (%) cold/warm shock: A 0.1±0.03 / 0.28±0.08; B 0.17±0.03 / 0.55±0.08; C 0.23±0.03 / 0.36±0.05;
  D 0.4±0.2 / 0.7±0.07; E 0.17±0.03 / 0.56±0.05; F 0.46±0.05 / 0.57±0.07.
- Table 2, downstream kBT/me c² (warm): 0.091, 0.183, 0.217, 0.28, 0.368, 0.732 for A2-F2, a factor 3.5±1 below
  the equipartition value 6.5·10⁻³ mi/me (Eq. 3); kBT ≃ 0.05 mi vsh²/2 (Sect. 3.3).
- Table 2, mean kinetic energy of the 1000 most energetic electrons, in me c²: ≈ 2 (A), 5.3-6.1 (B, C), 8.3-11.6
  (D, E), 22.7 (F).
- Only 6 % of the BI electrons are in the HE sample in run E2 (mi/me = 200), about three times the chance rate
  (Sect. 3.2). HE electrons gain non-adiabatically at the ramp and overshoot (magnetic moment rises, Fig. 4f).
- Mean surfing gain <γ−1> ≈ 0.155 (0.146-0.163 over all runs), independent of mi/me (Sect. 3.3); available beam
  energy ≈ 2 me vsh² ≃ 0.15 me c² (Eq. 5); 0.6 % of it goes to electrons in plane, 13 % out of plane (Sect. 3.3).
- Fraction of surfing electrons still in the tail downstream: ~30 % (A), 10 % (B, C), 5 % (D, E), 0.3 % (F);
  NTEF owed to surfing: 0.15 % (A), 0.05 % (B), 0.025 % (E), 0.0015 % (F) (Sect. 3.3).
- Probability to end above 5 kBT against energy after the Buneman zone is one curve for all runs once the energy
  is scaled by the downstream temperature; ≥ 10 % above kBT, ≈ 30 % above 3 kBT (Fig. 5b).
- Estimate for a 3-D perpendicular shock at mi/me = 100, combining run B with the out-of-plane run G of Paper I:
  NTEF ≈ 2 %, of which ≈ 1.5 % from surfing; the 3-D oblique run of Matsumoto et al. (2017) has ≈ 1.5 % (Sect. 3.3).
- DSA injection needs pe ≫ mi vsh ≈ 25 MeV/c for vsh = 8000 km/s (Sect. 4).

## Assessment
- Shown by simulation: no monotonic relation between pre-accelerated fraction and NTEF (C against B, E against D);
  NTEF is lower and the temperature slightly higher at βe = 5·10⁻⁴ than at 0.5; the tail forms at the ramp and
  overshoot; a simple additive model of surfing plus Weibel-zone heating (Eq. 4) underpredicts the data, so
  surfing helps later energisation but does not dominate it.
- Stated or conjectured: the 3-D estimate of ≈ 2 % (it mixes two different 2-D geometries); the extrapolation to
  the real mass ratio (largest run mi/me = 400); surfing "may remain important" at realistic shock speed because
  Te ∝ mi vsh² is then at most mildly relativistic, but the scaling with shock speed is called unclear (Sect. 4).
- Limits stated: vsh = 0.263c, ten times a realistic 8000 km/s; in-plane 2-D; strictly perpendicular. The
  conclusions "apply only to strictly perpendicular shocks": at oblique shocks surfing helps reflect electrons
  upstream and feeds stochastic shock drift acceleration (SSDA), which gives power laws and higher energies.
- Shock drift acceleration as in Matsumoto et al. (2017) is not seen; 2-D or the 90° geometry is blamed (Sect. 3.2).
- The 5 kBT cut is called ad hoc; being suprathermal by a few kBT is far from DSA injection (Sect. 4).
- Disagreement with the group's 2017 summary, which named the Buneman zone the main acceleration site.
- Doubts: Ms in Table 1 differs from Paper I for the same runs; the out-of-plane mean gain is quoted as 0.31
  here and 0.42 in Paper I (Eq. 5); one sentence on run E2 (0.56 % against 0.35 % for HE electrons) is garbled.
