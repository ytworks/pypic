# Bohdan, Niemiec, Kobzar & Pohl (2017): Electron Pre-Acceleration at Nonrelativistic High-Mach-Number Perpendicular Shocks

- Ref: ApJ 847, 71 (2017) (volume/article as cited in Bohdan et al. 2019-2021; meta.txt has no journal_ref);
  DOI 10.3847/1538-4357/aa872a; arXiv:1708.05528 (v1, read here); peer-reviewed
- Read: 2026-11-12, introduction, set-up and conclusions (results skimmed for numbers); topics: spherical-shock-injection

## Claim
In 2-D PIC runs of perpendicular shocks at MA ≈ 32-36 the orientation of B0 relative to the simulation plane decides
the electron injection. Shock surfing on Buneman waves in the foot is the first stage for every orientation. Out of
plane it is strong and followed by adiabatic compression; in plane and at 45° it is weak and followed by
second-order Fermi-like scattering off Weibel turbulence. The nonthermal fraction is about ten times larger out of
plane, and it is modulated by the reformation cycle.

## Method
- PIC, 2D3V, modified TRISTAN; flow-flow set-up of Wieland et al. (2016) with two equal-density beams at ±0.2c
  (relative 0.38c). The simulation frame is the downstream frame of both shocks (Sect. 2).
- θBn = 90°. φ (angle of B0 to the plane) = 0° (runs A), 45° (B), 90° (C). Left shock cold (βp = 5·10⁻⁴, runs *1),
  right shock warm (βp = 0.5, runs *2).
- mi/me = 100; ωpe/Ωe = 12; vA = 8.29·10⁻³ c; 20 particles per cell per species; λse = 20Δ, λsi = 200Δ.
- Table 1: A1/A2 and B1/B2: MA 31.7, Ms 1550/49; C1/C2: MA 35.5, Ms 1581/55 (upstream frame). Shock speed in the
  upstream frame 0.263c (A, B; Γ = 5/3, r = 3.97) and 0.294c (C; Γ = 2, r = 2.98) (Sect. 2).
- Box: Ly = 24 λsi (A, B), 12 λsi (C), periodic; Lx up to 250 λsi; run time 8 Ωi⁻¹ (Ωi⁻¹ = 1200 ωpe⁻¹).
- Convergence checked with 10-40 particles per cell and mi/me = 50-400 (Sect. 2, no figures shown).

## Numbers to compare against
- Buneman driving condition Ms ≳ (1+α)/2 · sqrt(mi/me) · sqrt(Te/Ti) with α ≃ 0.2, Ti far upstream, Te local
  (Eq. 1); here Ms ≥ 6. Trapping condition MA ≥ (1+α)(mi/me)^(2/3) ≃ 26 (Eq. 2, from Matsumoto et al. 2012).
- Revised for a field at angle φ to the plane: MA ≥ sqrt(2/(1+sin²φ)) (1+α)(mi/me)^(2/3) (Eq. 6), giving MA ≃ 39
  (φ = 0°) and 35 (45°) with the measured reflection rates; driving condition revised likewise (Eq. 7).
- Electrostatic field |Eest|/B0: ~2 (C1, coherent wave trains), peaks ~1.7 but average ≲ 1 (A1, B1, patchy);
  20-50 % weaker but over a 20-30 % larger area at βp = 0.5 (Sect. 3.3, Fig. 6).
- Electrons with γ > 1.5 in the Buneman zone: about 2200 (A1), 7600 (B1), 3.4·10⁶ (C1) (Sect. 3.4.1, Fig. 9).
  Their number follows the area of strong field: Spearman 0.7-0.9 in all six runs (Fig. 10).
- Traced electrons: out of plane γ ≃ 2.5 after surfing, ≃ 6 at the overshoot, falling behind it (Fig. 11); in
  plane γ ≃ 1.5 after surfing, ≃ 9 after many scatterings (Fig. 13).
- Downstream nonthermal fraction (Table 2): A1 0.2±0.1 %, A2 0.7±0.1 %, B1 0.2±0.1 %, B2 0.5±0.1 %, C1 4±1 %,
  C2 7±1 %. A traced C2 population has 4.8±0.2 % carrying ~26 % of the electron energy (Sect. 3.4.2).
- Largest Lorentz factor γ ≈ 9 in all runs (Fig. 16). Downstream kT/mc² = 0.03-0.054 (Table 2), called 2-4 times
  the Rankine-Hugoniot prediction (Sect. 3.4.5; Paper II instead compares with the equipartition value, Eq. 3 there).
- Reformation period 1.55 Ωi⁻¹, independent of φ and βp; shock speed swings between 0.03c and 0.15c (Sect. 3.2).

## Assessment
- Shown by simulation: the φ dependence of Buneman coherence and amplitude, traced to how much of the reflected-ion
  motion lies in the plane (Figs. 6-8); two different second stages (adiabatic against stochastic); patchy
  downstream nonthermal fraction that follows the reformation phase (Fig. 15); in A and B the fraction at a fixed
  place drops with time, in C it does not (Sect. 3.4.4); reconnection islands at 0° and 45° (Fig. 4).
- Stated or conjectured: Eq. 6 explains the weak surfing of in-plane runs here and in the literature (Table 3);
  the 90° fraction is "probably a better indicator" of the 3-D efficiency in the foot; Γ = 2 may distort the 3-D
  shock structure at 90°; the role of reconnection is deferred.
- Limits stated: 2-D; mi/me = 100; ripples are not captured (wavelength exceeds Ly, Sect. 3.2); downstream is not
  homogenised within 8 Ωi⁻¹; the nonthermal fraction depends on the Maxwellian fit and on the region chosen.
- No clear trend of the number of γ > 3 electrons with βp, but a higher nonthermal fraction at βp = 0.5 for every φ.
- Later revisions by the same group: Paper I shows that meeting Eq. 6 does not restore the out-of-plane efficiency;
  Paper II finds that surfing contributes little to the in-plane downstream tail; Paper III assigns most of that
  tail to reconnection. The summary here still calls the Buneman zone the main acceleration site.
- Against Wieland et al. (2016): 45° behaves like the in-plane case, not like 3-D; that run meets Eqs. 2 and 6
  yet shows almost no tail, which is attributed to spectra taken far downstream (Table 3).
