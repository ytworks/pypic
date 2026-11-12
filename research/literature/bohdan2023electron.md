# Bohdan (2023): Electron acceleration at supernova remnants

- Ref: Plasma Physics and Controlled Fusion 65 (2023) (volume as cited in Fulat et al. 2026; article number not
  checked; meta.txt comment "Accepted for publication in Plasma Physics and Controlled Fusion");
  DOI 10.1088/1361-6587/aca5b2; arXiv:2211.13992 (v1, read here); peer-reviewed review
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
A short review of what PIC simulations say about electron injection at SNR shocks, ordered by obliquity. It holds
that perpendicular shocks are understood but cannot start diffusive shock acceleration (DSA) by themselves, and
that oblique and parallel shocks do produce power-law electron tails in simulations, but with parameters whose
scaling to real remnants is unknown.

## Method
- Review of fully kinetic simulations, mostly by the author's group. Example runs: perpendicular 2-D (vsh = 0.263c,
  MA = 68.9, Ms = 106, mi/me = 400; Fig. 3) and (MA = 22.6, Ms = 35, mi/me = 50; Fig. 6); oblique 2-D, θBn = 60°
  (vsh = 0.263c, MA = 30, Ms = 32.5, mi/me = 50; Figs. 5, 6; Bohdan et al. 2022); parallel 2-D narrow box
  (vsh = 0.4c, MA = 16, Ms = 35, mi/me = 16; Figs. 7, 8; Arbutina & Zeković 2021).
- Classes (Sect. 2): perpendicular 75° ≲ θBn ≤ 90°, oblique 50° ≲ θBn ≲ 75°, parallel 0° ≤ θBn ≲ 50°; the
  boundaries are "not very strict" and may depend on the shock parameters. SNR values assumed:
  vsh = 1000-10,000 km/s = (0.003-0.03)c; Ms, MA ≈ 20-1000.

## Numbers to compare against
- Injection gap: εinj,e/εth,e ≈ 2ξ mpc/(mevsh) ≈ 10⁵-10⁶ with ξ ≈ 3; for protons εinj/εth = ξ² ≈ 10 (Eq. 1, Sect. 1).
  In simulations: ≈ 10³, (γ − 1)inj ≈ 35 for the oblique run (Sect. 4.3); ≈ 50, ≈ 10 for the parallel run (Sect. 5.3).
- Perpendicular (Sect. 3): width ≈ 0.5 ri (ri: upstream ion gyroradius); Bsh/B0 ∝ √MA; Weibel growth rate in
  units of Ωi independent of vsh/c and mi/me; Te/Ti ≈ 0.09-0.25 predicted for real parameters; suprathermal
  fraction about 0.5 %; maximum energy 10-100 times the downstream thermal energy; almost all electrons take part
  in reconnection for MA ≥ 70, yet its share of acceleration falls.
- Oblique (Sect. 4): transition ≈ 12.5 ri wide in the run shown; steady state of the foreshock needs slightly more
  than 300 Ωi⁻¹ at MA = 30, width then ≈ 200 ri; whistler turbulence δB/B ~ 1; downstream temperature within 20 %
  of the perpendicular run; tail ε^−2.7; expected Te/Ti ≈ 0.1-0.3.
- Electron-acoustic waves raise energies by a factor ~10-50 but affect less than 0.1 % of incoming electrons;
  whistler interaction is said to accelerate "up to 7 % of electrons" at high Mach numbers (Sect. 4.3).
- Parallel (Sect. 5): size with foreshock ≈ 160 ri (up to 10³ ri in hybrid runs); Te/Ti ≈ 0.5-1; tail ε^−2.
- Table 1 (status): perpendicular: structure, heating, acceleration all "well studied"; oblique and parallel:
  structure and acceleration "partially studied", heating "not studied".

## Assessment
- Mechanisms recognised: perpendicular: shock surfing on Buneman waves (needs the trapping condition; spoiled by
  heating at the real mass ratio), reconnection of Weibel filaments, stochastic Fermi acceleration (the most
  efficient; maximum energy ∝ mi/me) (Sect. 3.3). Oblique: these three plus electron-acoustic waves, whistlers and
  stochastic shock-drift acceleration (Sect. 4.3). Parallel: mirror reflection, then repeated shock-drift cycles
  between the front and Bell-type upstream waves (Sect. 5.3).
- Simulations by obliquity: perpendicular 2-D (refs. 25-36) and the 3-D run of Matsumoto et al. (2017), which has
  θBn ≈ 74.3° but is called too short to form a foreshock; oblique 1-D (Xu et al. 2020; Kumar & Reville 2021) and
  2-D (Bohdan et al. 2022; Morris et al. 2022); parallel 1-D (Park et al. 2015; Kato 2015), 2-D narrow box
  (Arbutina & Zeković 2021) and hybrid (Caprioli & Spitkovsky 2014). No curved or global shock simulation is
  mentioned; all results are for planar shocks in a homogeneous medium, and rippling appears only as something
  the narrow parallel run cannot capture (Sect. 5.1).
- Open questions listed: heating at oblique and parallel shocks; scaling of every mechanism with vsh, mi/me and
  Mach number; multi-dimensional effects at parallel shocks; the long run times oblique shocks need; how to
  compare Te/Ti with remote observations (Sects. 3.2, 4-6).
- Numerical limits named: reduced injection gap; 1-D and out-of-plane 2-D runs lack the Weibel instability; small
  mass ratios miss the intermediate-scale instability (Shalaby et al. 2021); the parallel 2-D run is effectively 1-D.
- States, without new evidence, that perpendicular shocks are "not so important in astrophysical sense" because
  nothing escapes upstream; DSA there would need seed cosmic rays, neutrals or pre-existing turbulence.
- "7 % of electrons" disagrees with Morris et al. (2023): 0.8 % of electrons holding 7 % of the energy. The
  factor 10-50 for electron-acoustic waves is not in Morris et al. (2022), which quotes about 1 %, not 0.1 %.
  The threshold used (gyroradius above the shock width) is what Amano et al. (2022, Sect. 2.5) call pessimistic.
