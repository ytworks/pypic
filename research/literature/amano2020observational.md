# Amano et al. (2020): Observational Evidence for Stochastic Shock Drift Acceleration of Electrons at the Earth's Bow Shock

- Ref: Phys. Rev. Lett. 124, 065101 (2020); DOI 10.1103/PhysRevLett.124.065101; arXiv:2002.06787 (v1, read here);
  peer-reviewed publication. The Supplemental Material is not part of the arXiv text and was not read.
- Read: 2026-11-12, full (main text); topics: spherical-shock-injection

## Claim
One MMS crossing of a nearly perpendicular bow shock shows electrons diffusing inside the shock transition layer:
intensities rise exponentially toward downstream, the distribution is nearly isotropic, and whistler power is
high. The scattering rate inferred from the profiles and the measured wave power both exceed the SSDA thresholds
exactly over the energy range of the power law, and the predicted cutoff matches the fitted one. The authors take
this as direct evidence for SSDA and extrapolate to mildly relativistic electrons at fast SNR shocks.

## Method
- Observation: MMS (four spacecraft, 1 s averages), 2016-12-09 around 10:29 UT; electrons from FPI and FEEPS,
  fields from FGM/SCM/EDP. Shock parameters: MA ≃ 8.9, θBn ≃ 85°, normal solar-wind speed u0 ≃ 590 km/s
  (spacecraft frame), shock speed relative to the spacecraft U0 ≃ 30 km/s from the foot thickness (p. 1).
  Plasma beta: not stated. fce ≃ 700 Hz (p. 4). MA/cosθBn ≈ 102 (my arithmetic).
- Transport model: diffusion–convection equation for the isotropic part in the HTF (Eq. 1), from the focused
  transport equation of Katou & Amano (2019) under weak anisotropy; κ = v²/(6Dμμ).
- Test 1: fit f(E,t) ∝ exp[t/τ(E)] over 10:29:10–10:29:20, convert to Dμμ/Ωce =
  (1/6)(me/mi)(E/E_sh)(u0/Ωci)/(U0 τ) assuming a stationary structure (p. 4), compare with the threshold.
- Test 2: quasi-linear threshold wave power, with finite-frequency and dispersion corrections, against the
  measured magnetic spectrum (Fig. 4C–D).

## Numbers to compare against
- Threshold: D*μμ/Ωce = (1/6η)(me/mi)(E/E_sh), E_sh = me u_sh²/2, u_sh = u0/cosθBn (unnumbered, p. 3);
  η ~ 0.4 from ≈ 10 s of strong waves within a ≈ 25 s layer; the acceleration region is η u0/Ωci thick (p. 4).
- Threshold wave power: P*(f)f/B0² = (2/3πη)(me/mi)(B/B0)(E/E_sh)[(1 + |ω/kvμ|) dlnω/dlnk]^-1 (unnumbered, p. 4).
- Cutoff: predicted ~20 keV from Fig. 4B; fitted E_cutoff ≃ 22 ± 1 keV with f ∝ E^-p exp(−E/E_cutoff) (p. 4).
  The legend of Fig. 4A gives p = 3.5 (read from the rendered figure; not in the text).
- Scaling: E_cutoff ≃ 20 keV (u0/600 km/s)² (cosθBn/cos85°)^-2 (Dμμ/0.03Ωce) (unnumbered, p. 5).
- Waves: the tightest condition is at f ≳ 0.1 fce, where whistlers resonate with 0.1–1 keV electrons; the inferred
  scattering rate saturates at a few keV where δB ~ B (p. 4–5, Fig. 4B).
- Anisotropy: 0.2–1 keV electrons stream away from the shock before 10:29:10 and become nearly isotropic after
  it; weak perpendicular anisotropy near 1 keV inside the acceleration region, variable at 30 ms (p. 2, 5).
- Extrapolation: heliospheric shocks (~400 km/s) do not normally give relativistic electrons; young SNR shocks
  (≳ 3000 km/s) with relatively large obliquity give more than a few hundred keV inside the layer (p. 5).

## Assessment
- Shown (observation, one event): the three diffusion signatures coincide in time; two independent estimates
  (Fig. 4B, 4D) lie above the SSDA thresholds over the power-law range; the cutoff agrees within the
  order-of-magnitude accuracy the authors claim.
- The authors themselves say the power-law index is not convincing evidence, because power laws are common; the
  cutoff–scattering relation is the test.
- Stated/conjectured: whistlers come from SDA-produced perpendicular anisotropy, so SSDA may sustain itself; SNR
  injection; turbulence ahead of a quasi-parallel shock can make locally near-perpendicular patches, which is
  offered as the link to the radial fields seen in young SNRs (p. 5). None of these is tested here.
- Stated limits: stationarity during the crossing; order-of-magnitude threshold; quasi-linear estimate for
  parallel waves while lower-frequency waves are oblique and nonlinear.
- Not addressed: bow-shock curvature. Katou & Amano (2019) quote a curvature limit of ~10 keV for the bow shock,
  the same order as the 22 keV cutoff attributed here to scattering alone.
- Corrects Katou & Amano (2019): the factor 1/(6η) was missing there (p. 4).
- Doubts: the printed normal vector (+0.979, −0.827, −0.525) is not a unit vector; the text says the threshold
  wave power scales as (MA/cosθBn)², while the formula (∝ E/E_sh) gives the power −2 at fixed energy; the caption
  of Fig. 4D gives an averaging interval starting at 10:19:15, presumably 10:29:15; Eq. 1 as printed has no
  factor u_sh in its energy-gain term, which is dimensionally inconsistent (compare Katou & Amano 2019, Eq. 14).
