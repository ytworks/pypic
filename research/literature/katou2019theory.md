# Katou & Amano (2019): Theory of Stochastic Shock Drift Acceleration for Electrons in the Shock Transition Region

- Ref: ApJ 874, 119 (2019) (volume/article as cited by Amano et al. 2020, ref. 21; the arXiv metadata gives only the
  DOI); DOI 10.3847/1538-4357/ab0d8a; arXiv:1903.02277 (v1, submitted version, read here); peer-reviewed (accepted by ApJ)
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
Pitch-angle scattering added to adiabatic shock drift acceleration (SDA) keeps electrons inside the transition layer
of a quasi-perpendicular shock, where the gradient-B drift energizes them continuously (stochastic SDA, SSDA). In
the strong-scattering limit the layer-averaged steady spectrum is a power law whose index depends only on the field
gradient scale; it cuts off where the diffusion length reaches the layer thickness, at an energy linear in the
scattering rate and in (u0/cosθBn)². The authors propose SSDA as the injection step at fast SNR shocks.

## Method
- Theory in the de Hoffmann–Teller frame (HTF): focused transport equation with a phenomenological, elastic
  (plasma-frame) pitch-angle diffusion term (Eqs. 3–6); flow along B with u∥ = u_sh = u0/cosθBn held constant,
  justified by Eq. 5 (velocity-gradient term is O(cos²θBn) of the field-gradient term; App. A, MA ≫ 1 assumed).
- Box model: average over the layer 0 < s < L_sh, L_sh ≈ u0/(Ωci cosθBn) along B (Eq. 9); constant gradient scale
  1/L = ∂lnB/∂s (Eq. 7); isotropic limit with downstream escape rate u_sh/L_sh (Eq. 16) and an upstream escape rate
  that switches on when l_diff = v²/(6 Dμμ u_sh) approaches L_sh (Eqs. 17–19). Dμμ constant.
- Monte Carlo test of the box model with anisotropy: stochastic differential equations (Eqs. 22–23), Milstein
  scheme, 1e5 particles injected isotropically at v = u_sh, Dμμ Δt = 5e-4, particle splitting; escape rate taken
  from the isotropic estimate (Sect. 3.1). No PIC, no new data.
- Comparison with the bow shock through quasi-linear theory for parallel whistlers (Eqs. 28–30) and published wave
  amplitudes; B0 = 10 nT, ne = 10 cm⁻³, u0 = 400 km/s (Sect. 4). Cross-shock potential and second-order Fermi neglected.

## Numbers to compare against
- Adiabatic SDA: Δv = 2u0/cosθBn = 2u_sh in the upstream frame (Eq. 1), valid while the shock is subluminal;
  transmitted electrons gain at most a factor ≈ 4 at a strong shock (Sect. 2.1). Reflection: HTF pitch angle outside
  the loss cone set by the field compression (Sect. 2.1, no equation).
- Spectrum: N(v) ∝ v^-p, p = 1 + 3L/L_sh (Eq. 20); p ~ 4, i.e. f(ε) ∝ ε^-3 for L ~ L_sh; independent of Dμμ.
  Monte Carlo agrees for L/L_sh = 0.5–2.0 at Dμμ/Ωce = 0.01 (Fig. 3a).
- Cutoff: ε_max ~ ε_sh (mi/me)(Dμμ/Ωce) = (1/2) mi u_sh² Dμμ/Ωce, ε_sh = me u_sh²/2 (Eq. 21); linear in Dμμ for
  Dμμ/Ωce = 0.01–0.1 (Fig. 3b; cutoff defined as a 10 % deviation from the power-law fit).
- Sustaining condition at energy ε: Dμμ/Ωce ≳ (me/mi)(ε/ε_sh) (Eq. 31); wave power I(k_res)k_res/B0² ≳
  (4/π)(me/mi)(ε/ε_sh) (Eq. 32). Controlling parameter MA* = u_sh/VA = MA/cosθBn at fixed VA/c (Sect. 4, Fig. 5,
  VA/c = 2e-4, MA* = 1e2–1e4).
- Bow shock: electrons of 0.1 keV / 10 keV resonate at ≈ 40 Hz / 1 Hz in the upstream field (Fig. 6a); required
  I(f)f/B0² ~ 1e-4–1e-2 for θBn ≳ 80° (Fig. 6b) against observed δB/B0 ~ 0.01–0.1 at ω/Ωce ~ 0.1–0.4 (Eq. 33):
  enough for ε ≲ 1 keV only at θBn ≳ 80°.
- ε_max ≈ 11 keV (u0/400 km/s)² (cosθBn/cos85°)^-2 (Dμμ/0.1Ωce) (Eq. 34); Dμμ/Ωce ~ 0.1 taken as realistic, the
  Bohm value ~1 called too optimistic. Shock speed ten times the bow shock's gives ~1 MeV (Sect. 4).
- Bow-shock curvature limits SDA and SSDA alike to ~10 keV (Sect. 4, citing Krauss-Varban & Burgess 1991).
- Resonance with MHD-scale waves (kλi ≪ 1) needs v ~ c (Fig. 4), the usual injection threshold.

## Assessment
- Shown (theory + Monte Carlo of the same box model): power law with index set by L/L_sh, cutoff ∝ Dμμ; weak,
  nearly linear anisotropy that grows toward the cutoff and does not change the index (Fig. 2b).
- Not shown: the scattering itself. Dμμ is an input; its size at the bow shock is an order-of-magnitude
  quasi-linear estimate from literature wave amplitudes, and the wave source (modified two-stream, temperature
  anisotropy, loss-cone beam) is left open (Sect. 4).
- Stated/conjectured: SNR shocks reach mildly relativistic energies by SSDA and then DSA takes over (Sect. 4, 5).
- Stated limits: spatial structure of the layer not solved; escape rate isotropic although anisotropy grows near
  the cutoff; plane, steady, subluminal, quasi-perpendicular (cosθBn ≪ 1) shock; v ≫ u_sh; energies well above
  the cross-shock potential; quasi-linear assumptions doubtful inside the layer.
- Geometry: the curvature limit is said not to be intrinsic; larger shocks reach higher energies. The escape is
  attributed to the change of local θBn along the field line as the electron drifts (Sect. 4).
- Later correction: Amano et al. (2020) note a missing factor 1/(6η) in the threshold of Eq. 31 (η: fraction of
  the layer with strong waves). Eq. 21 follows from l_diff ~ L_sh with the factor 6 dropped (my check).
- Open: self-consistent origin of ≳ 0.1 Ωce whistlers; PIC and in-situ tests (Sect. 5).
