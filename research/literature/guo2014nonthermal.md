# Guo, Sironi & Narayan (2014): Non-Thermal Electron Acceleration in Low Mach Number Collisionless Shocks. I. Particle Energy Spectra and Acceleration Mechanism

- Ref: ApJ 794, 153 (2014); DOI 10.1088/0004-637X/794/2/153; arXiv:1406.5190 (v2, "final version accepted to ApJ",
  read here); peer-reviewed publication
- Read: 2026-11-12, introduction, set-up, Sect. 4, conclusions and appendices; topics: spherical-shock-injection

## Claim
At a low Mach number, high-beta, quasi-perpendicular shock, thermal electrons are mirror-reflected and energized
by shock drift acceleration (SDA). The reflected electrons stream upstream, make the electron pressure
anisotropic (T∥ > T⊥) and drive magnetic waves; the waves scatter later reflected electrons back to the shock for
further SDA cycles. This Fermi-like loop needs no downstream scattering and produces a power-law tail.

## Method
- PIC, TRISTAN-MP. Reference run 2-D with the upstream field in the simulation plane (ϕB = 0°); comparison runs
  in 1-D, 2-D out-of-plane (ϕB = 90°) and 3-D (App. A). Shock made by reflection off a wall; the simulation frame
  is the downstream frame.
- mi/me = 100 (one run with 400 and θB = 77°, App. B). Ti = Te = 1e9 K = 86 keV, upstream flow u0 = 0.15c.
- Mach numbers: M = u0/cs = 2 in the simulation frame (Eq. 1), Ms = 3 in the shock rest frame (Eq. 3).
  Magnetization σ = 0.03, with MA ≡ u0/vA = (2/σ)^(1/2) (Sect. 2), i.e. ≈ 8.2 in the simulation frame (my
  arithmetic; the value is not printed). Plasma beta βp = 20 (total, Te = Ti).
- θBn = 63°. Shock speed 0.076c in the simulation frame (Sect. 3), ≈ 0.22c in the upstream frame (my arithmetic
  from Eq. 9). ωpe/Ωce not stated; ≈ 5.4 from the times in the Fig. 2 caption (my arithmetic).
- Box: transverse 76 c/ωpe = 7.6 c/ωpi ≈ one ion Larmor radius; expanding along x with a moving injector
  (length not stated; shock at x ≃ 1115 c/ωpe at Ωci t = 26.9). 10 cells per c/ωpe, dt = 0.045/ωpe, 32 particles
  per cell (16 per species). Run time up to Ωci t = 29.8 (Fig. 4). 3-D run: 76 c/ωpe in y and z, 5 cells per
  c/ωpe, 2 particles per cell, to Ωci t = 7.4.
- Theory: relativistic SDA in the de Hoffmann–Teller (HT) frame with a cross-shock potential; 12800 traced electrons.

## Numbers to compare against
- HT frame speed relative to the upstream: u_t = u_sh^up secθB = (u_sh + u0)/(1 + u_sh u0/c²) secθB (Eq. 9). It
  exists only for u_t ≤ c; the reference parameters are superluminal for θB ≳ 77° (Sect. 4.2.1).
- Reflection in the HT frame: v∥ < 0 (Eq. 13) and v⊥ ≥ c {(B0/B)[(γ + Δφ)² − 1]/γ²}^(1/2) (Eq. 14), equivalently
  a minimum pitch angle (Eq. 15); Δφ = e(Φ − Φ0)/me c². Upstream frame, Δφ → 0: v∥ < u_t (Eq. 19) and
  v⊥ ≥ γ_t (u_t − v∥) tanα0, sin²α0 = B0/B (Eq. 20); least speed u_t [tan²α0/(tan²α0 + 1/γ_t²)]^(1/2) (Eq. 22).
- Energy after reflection: γ_r = γ_i [1 + 2u_t(u_t − v_i∥)/(c² − u_t²)] (Eq. 24).
- Drift picture: Δγ = (−e/me c²)∫Ez dz ≃ −(σ/2)^(1/2) (mi/me)^(1/2) (u0/c)² sinθB Δz/(c/ωpe) (Eqs. 7–8). One
  cycle lasts ~1/Ωci; the example electron drifts 200 c/ωpe and gains Δγ ≃ 6 (Fig. 5).
- Fit of the reflected population: field compression b = 4, Δφ = 0.5 (Fig. 6).
- About 15 % of electrons in the non-thermal tail, slope p ≃ 2.4 in dN/dγ ∝ (γ − 1)^-p (Fig. 4c); they carry
  ≃ 10 % of the incoming ion bulk kinetic energy density (Sect. 5).
- γ_max ≳ 20 at Ωci t = 28.9 against ~6 for a single SDA cycle (Fig. 7b). A downstream tail appears after ~20/Ωci.
- Ions: ~20 % reflected, no acceleration at this obliquity (Fig. 4a).
- Mass-ratio rule: results unchanged if secθB (me/mi)^(1/2) is held fixed (Eq. B1): 63° at 100, 77° at 400,
  84° at 1836 (App. B).

## Assessment
- Shown (simulation): reflected electrons obey the SDA reflection region and energy mapping; a one-cycle
  synthetic spectrum matches the early spectrum (Fig. 7b); later growth comes from repeated SDA after upstream
  scattering (Fig. 8); returning electrons, not ions, occupy the wave region (Sect. 3). 2-D in-plane agrees with
  3-D up to Ωci t = 7.4, while 1-D and 2-D out-of-plane give a higher normalization and slower growth (Fig. 9).
- Stated, deferred to paper II: the nature of the waves; validity across obliquity, temperature (1e7–1e9 K) and
  σ = 0.01–0.03; waves need the parallel electron pressure above the magnetic pressure and vanish at low beta.
- Stated without demonstration: SDA "cannot operate efficiently" at high Mach number because the reflecting
  region of velocity space shrinks with u_t (Sect. 5, from Eq. 22); the reflected fraction also falls with Δφ.
- Limits: electrons start mildly relativistic (kTe/me c² = 0.17), far from SNR upstream conditions; reduced
  mass ratio with rescaled angle; transverse box of one ion Larmor radius; planar shock; one obliquity here.
- Disagreements it states: thermal-leakage injection (Kang & Ryu 2010) needs far larger momenta; a
  κ-distribution upstream (Park et al. 2013) inflates injection; whistler injection of Riquelme & Spitkovsky
  (2011) is sensitive to mass ratio and belongs to low beta; Matsukiyo et al. (2011, 1-D) and Park et al. (2012,
  2013) saw SDA but no sustained Fermi process.
