# Amano & Hoshino (2022): Theory of Electron Injection at Oblique Shock of Finite Thickness

- Ref: ApJ (2022), volume/article number not checked; DOI 10.3847/1538-4357/ac4f49; arXiv:2201.11416 (v1, submitted version, read here)
- Read: 2026-11-12, full; topics: (none yet)

## Claim
Stochastic shock drift acceleration (SSDA) and DSA are two limits of one diffusion-convection equation at an
oblique shock of finite thickness Ls. The ratio of diffusion length l_diff = κ cos²θ/Us to Ls separates them.
SSDA (l_diff ≲ Ls) acts inside the transition layer and is fast, but its spectra are steep (q ≳ 6). Electrons
enter DSA when the SSDA maximum energy reaches a threshold set by the turbulence dissipation scale k* outside
the layer. That threshold is about 0.1–1 MeV for typical conditions, and injection is favoured at high
MA/cosθBn (the HTF Alfvén Mach number).

## Method
- Theory: diffusion-convection equation in the de Hoffmann–Teller frame (Eq. 1), derived from relativistic
  focused transport with weak anisotropy (App. A); integral spectral-index formula from known profiles (Eq. 6);
  constant-gradient estimate (Eq. 8); SSDA maximum energy (Eq. 10); resonance threshold (Eq. 11); injection
  condition (Eqs. 12–14).
- Numerics: steady state of Eq. 1 with a tanh Bz profile (Eq. 15), Vx Bz = const (Eq. 16), and a free-escape
  boundary at −L_FEB. Analytic external solutions (Eq. 19) are matched to a Chebyshev pseudospectral interior
  (App. B). Fixed r = 4, X2 = −X1 = 4Ls. Diffusion models: constant κ (Figs. 2–3); Bohm-like κ ∝ p (Fig. 4);
  a "realistic" model with pitch-angle rate Dμμ (Eq. 20) whose layer part follows MMS bow-shock data and whose
  upstream part is cut below p*^u (Figs. 5–9).
- No PIC and no new observations. The paper says Sect. 3 does not aim at comparison with specific data or simulations.

## Numbers to compare against
- DSA limit q = 3r/(r − 1) (Eq. 5); SSDA q ≈ 3[1 + (l_diff ∂lnB/∂x)^(-1)] (Eq. 8), so q ≳ 6 for l_diff ≲ Ls (Sect. 2.4).
- Ls = η Us/Ωci, η ~ 1 (Sect. 2.5); E_max,SSDA/E0 ≃ 6η Dμμ/Ωci, E0 = me (Us/cosθBn)²/2, non-relativistic (Eq. 10).
- Injection: MA/cosθBn ≳ 55 η^(-1/2) [(Dμμ/Ωce)/0.1]^(-1/2) (k*c/ωpi)^(-1) (Eq. 14). This gives MA ≳ 17 for a
  quasi-parallel shock even with Bohm scattering, and MA ~ 5 at θBn ~ 85° with Dμμ/Ωce ~ 0.1 (Sect. 2.7).
- Threshold energy 0.1–1 MeV for VA ~ 100 km/s, k*c/ωpi ~ 0.1–1 (Fig. 1, Sect. 2.6).
- Fiducial model: ν = 2, D̃ = 3e-2, p*/me c = 0.1 (≈ 2.5 keV); upstream ν = 6, D̃ = 3e-3; MA = 10, θBn = 85°,
  VA/c = 2e-4, η = 0.5 (Sect. 3.4). SSDA dominates below ≈ 25 keV, pure power law below ≈ 5 keV (Sect. 3.4, Fig. 5).
- Spectra: q ~ 6 at low energy, q = 4 at high energy for k*c/ωpi = 0.32 and 1.0, cutoff for 0.1 (Fig. 6);
  hardening above E/me c² ~ 3e-2 only for θBn = 85° (Fig. 7); Bohm-like case reaches l_diff/Ls ≳ 10 at p/p1 ≳ 30 (Fig. 4).
- p*/p1 ~ (k*c/ωpi)^(-1)(mi/me)(MA/cosθBn)^(-1). Bow shock MA = 8, θBn = 85° gives ≈ 20, i.e. an energy gain
  of ≈ 400. The real mass ratio needs MA/cosθBn ≳ 1e3 for p*/p1 ~ 1 (Eq. 21, Sect. 4.4).
- Bow-shock curvature limit E ≲ 31 keV for Dμμ/Ωce = 0.1, obliquity change δ = 1°, Us = 400 km/s, B = 10 nT, Rc = 20 RE (Eq. 22, Sect. 4.5).

## Assessment
- Shown (theory and numerical solutions of the model): DSA and SSDA emerge from the same equation. Eq. 6 agrees
  with numerical slopes (Fig. 3). The spectrum turns from steep to the DSA slope where l_diff crosses Ls,
  provided l_diff ≪ L_FEB there (Figs. 4, 6, 7). Spatial variation of spectra (Figs. 8–9) is proposed as a test.
- Inputs, not results: the scattering rates (inside the layer from a single MMS crossing; upstream form and
  ν = 6 chosen); results depend on L_FEB partly through that choice (Sect. 3.4).
- Stated/conjectured: SSDA called the most promising injection candidate (Sect. 1); self-sustained whistler
  generation keeping l_diff/Ls ≈ const (speculation, Sect. 4.1); non-stationarity changes efficiency only
  quantitatively (Sect. 4.3); small-mass-ratio PIC artificially eases electron injection (Sect. 4.4); bow shock
  not an efficient injector because of curvature (Sect. 4.5); a laser-driven Weibel shock spectrum fits SSDA (Sect. 4.6).
- Stated limits: steady, planar, locally subluminal shock (HTF); weak anisotropy fails for relativistic HTF
  speed (Sects. 2.6, 4.4) and for electrons slower than Us/cosθBn (pitch-angle equation needed, Sect. 4.2); O(1) factors dropped;
  k* and Dμμ are not known from first principles; reduced PIC parameters (fast shocks, small mi/me) distort comparisons (Sect. 4.4).
- Open: time-dependent solutions for PIC comparison; wave intensity and direction vs MA/cosθBn; the
  magnetized-to-Weibel shock transition; heating as the seed normalization.
