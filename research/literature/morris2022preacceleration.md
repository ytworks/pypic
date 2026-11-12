# Morris, Bohdan, Weidl & Pohl (2022): Pre-acceleration in the Electron Foreshock I: Electron Acoustic Waves

- Ref: ApJ 931, 129 (2022) (volume/article as cited in Morris et al. 2023 and Fulat et al. 2026; meta.txt has no
  journal_ref); DOI 10.3847/1538-4357/ac69c7; arXiv:2204.11569 (v1, read here); peer-reviewed
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
At oblique high-Mach-number shocks the share of incoming electrons that is reflected upstream grows steeply as θBn
decreases (none at 74.3°, about 5 % at 30°). The reflected beam drives electrostatic electron-acoustic waves (EAWs) in
the outer foreshock. These waves change the field-aligned momentum of incoming electrons; they turn a small share
back before it reaches the shock and send some reflected electrons back to the shock with more than thermal energy.

## Method
- PIC, 2D3V, code derived from TRISTAN (Vay pusher); one shock from a reflecting wall, box grows with time (Sect. 2).
- mi/me = 64; MA ≈ 21 (defined in the upstream frame); Ms not stated (cs = 0.00814c is given); β = 1; Te = Ti;
  vth,e = 0.0713c; inflow 0.2c, shock speed 0.264c (upstream frame), 0.0673c (simulation frame); ωpe/Ωe = 9.9.
- B0 = B0 (cos θBn, sin θBn cos φ, sin θBn sin φ); φ is the angle to the simulation plane. Table 1: θBn = 30°, 45°,
  63° with φ = 0°; 63° also with φ = 45°, 90°; 74.3° with φ = 0°, 90°. A run at 55° (φ = 0°) is used in Sect. 5 and
  listed in Sect. 6 but is missing from Table 1.
- Transverse box 768Δ = 4.8 λsi (about 55 thermal electron gyroradii); λse = 20Δ; 20 particles per cell per
  species; δt = 1/40 ωpe⁻¹; run time 7.8 Ωci⁻¹ (chosen to match the 3-D run of Matsumoto et al. 2017).
- Reflection rate R = Nref/(Nref + Nups), from the two-peaked electron spectrum in the region 2-5 λsi ahead of the
  foot (simulation frame), split at the minimum between the peaks (Sect. 4.1). 10,000 electrons traced per run.

## Numbers to compare against
- Superluminal limit for these runs: θcrit = arctan(c/vsh) ≈ 75.2° (Sect. 2).
- Reflection needs ve cos θBn ≥ vsh in the shock frame (Sect. 4.1.1, after Amano & Hoshino 2007).
- In plane (φ = 0°): R ≈ 5 % at 30°, ≈ 0.05 % at 63°, none at 74.3°; later onset at larger θBn (Sect. 4.1.1, Fig. 4).
  The value at 45° is only in Fig. 4. Reflected electrons are fewer but more energetic at larger θBn (Fig. 3).
- θBn = 63°: R ≈ 0.03 % for φ = 0° and 45°, ≈ 4 % for φ = 90°; no reflection at 74.3° for any φ (Sect. 4.1.2, Fig. 5).
- EAW check at 30°: beam p0 = 0.92 mec, pth = 0.309 mec, nref = 0.01 nups give a peak growth rate 0.024 ωpe at
  k∥λse = 1.76 (Eq. 2); measured k∥λse ≈ 1.79 (Sect. 4.2, Fig. 7). Noise floor UES/U0 ≈ 7·10⁻⁴ (Fig. 8).
- 70-80 % of reflected electrons meet the shock once. Two or more encounters: ≈ 30 % at 63°, ≈ 9 % at 45° (Sect. 5).
- Traced foreshock electrons turned back by EAWs without reaching the shock: ≈ 15 % at 45°, 70 of 10,000 at 55°,
  negligible above (Sect. 5). Conclusion: injection is compromised for ≈ 1 % of upstream electrons at θBn ≤ 45°.
- Quoted from Bohdan et al. (2022; θBn = 60°, φ = 90°, 50.4 Ωci⁻¹): the EAW region saturates at about 270 Ωci⁻¹;
  ion reflection and whistlers ahead of the ramp start near 10 Ωci⁻¹ (Sect. 5).

## Assessment
- Shown by simulation: R against θBn and φ; agreement of the wavenumber with linear theory of the electron-acoustic
  instability; correlation of wave energy with the energy of reflected electrons (Fig. 8); single-particle tracks
  with a two-step reflection (surfing on Buneman waves, then mirroring at constant magnetic moment; Figs. 9-11);
  a parallel electric field of the EAWs turning an electron far ahead of the shock (Fig. 14); the mean first
  reflection point moving away from the foot with time (Fig. 15).
- Stated, not shown: that electrons returning with higher energy are injected more easily into stochastic
  shock-drift or stochastic Fermi acceleration. No downstream spectrum and no injected fraction are measured.
- The authors' verdict on EAWs is mixed: they accelerate and reflect some upstream electrons, and they keep about
  1 % of them from reaching the shock. Later papers of the group quote this as "less than 0.1 %" (Bohdan 2023;
  Fulat et al. 2026).
- Limits: 7.8 Ωci⁻¹ is before whistlers and ion reflection appear; narrow box; mi/me = 64; vsh = 0.264c puts the
  superluminal limit at 75°, so the absence of reflection at 74.3° belongs to this shock speed (my reading; the
  paper does not discuss realistic speeds). In-plane runs capture EAWs but weaken Buneman surfing, out-of-plane
  runs do the opposite, so R differs by a factor of about 100 at 63°. The 30° shock is still forming and is
  dominated by an ion foreshock (Sect. 5).
- Inconsistencies: 0.05 % (Sect. 4.1.1) against 0.03 % (Sect. 4.1.2) for 63°, φ = 0°; Fig. 9 is for 55°. The
  abstract gives wave energy ∝ R^0.6 and a wavelength ≈ 2 λse; the body states neither, and k∥λse ≈ 1.79
  corresponds to 2π/k∥ ≈ 3.5 λse.
- Credits: mirror reflection Honda & Honda (2005); two-step reflection Amano & Hoshino (2007); EAW theory
  Gary (1987) and Bohdan et al. (2022).
