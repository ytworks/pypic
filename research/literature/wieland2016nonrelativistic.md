# Wieland, Pohl, Niemiec, Rafighi & Nishikawa (2016): Non-relativistic perpendicular shocks modeling young supernova remnants: nonstationary dynamics and particle acceleration at forward and reverse shocks

- Ref: ApJ 820, 62 (2016) (volume/article from the DOI and from the reference lists of Bohdan et al. 2017-2021;
  meta.txt has no journal_ref); DOI 10.3847/0004-637X/820/1/62; arXiv:1602.05064 (v1, "ApJ, in press", read here);
  peer-reviewed
- Read: 2026-11-12, introduction, set-up and conclusions (results skimmed for numbers); topics: spherical-shock-injection

## Claim
Two colliding plasma slabs (density ratio 10) with a strictly perpendicular field tilted 45° out of the 2-D plane
form a forward and a reverse shock of MA ≈ 28. Both are mediated by Weibel-type filaments, reform every ~1.5 Ωi⁻¹
and ripple on reflected-ion gyro-scales. Ions get a suprathermal tail, but electrons are only heated: the Buneman
waves in the foot are too weak to hold relativistic electrons, and neither reconnection nor gradient drift is seen.

## Method
- PIC, 2D3V, modified TRISTAN (4th-order FDTD, weak Friedman filter, Vay pusher). Flow-flow set-up with a
  tapered-field transition zone between the slabs that removes the start-up electromagnetic transient (Sect. 2.1, App. A).
- θBn = 90°; B0 at φ = 45° to the simulation plane, chosen to approximate 3-D physics (Sects. 1, 2.2).
- mi/me = 50; 10 particles per cell per species (weights make the density ratio); Ωe/ωpe,L = 0.032 (Sect. 2.2).
- Table 1, dense slab (reverse shock) | dilute slab (forward shock): λse 7.9Δ | 25Δ; electron thermal speed 0.002c
  both; streaming 0.0354c | −0.354c in the centre-of-momentum frame (relative 0.38c); vA 0.00447c | 0.0142c;
  shock speed in the upstream frame 0.127c | 0.39c; MA 28.5 | 27.6; Ms 252 | 755 (both in the upstream frames).
- β is not tabulated: "βe ≪ 1" (Sect. 1); electron β "of order 1/25" after far-upstream heating (Sect. 3.2.1).
- Box: Ly = 324.1 λse = 42.9 λsi,L, periodic; Lx grows to 9493.7 λse (moving injector, radiating x boundaries).
  Run time 20 Ωi⁻¹ = 31645.5 ωpe,L⁻¹ (Sect. 2.2). Lengths are in the skin depth of the dense plasma.

## Numbers to compare against
- Buneman threshold Ms ≥ (1+α)/2 · sqrt(mi/me) · sqrt(Te/Ti), α = reflected/incoming ion density (Eq. 4).
- Buneman wavenumber k∥ ≈ c/(vrel λse,local) (Eq. 3); seen at k∥λse ≃ 0.8 (early) and 0.7 (late, Fig. 11);
  k⊥/k∥ ≃ 0.15 here against ≃ 0.5 in Matsumoto et al. 2012 (Sect. 3.2.1).
- Trapping threshold after Matsumoto et al. 2012: MA ≈ 180 for mi/me = 1836, MA ≥ 16 for mi/me = 50 (Sect. 2.2).
- Expected momentum gain 0.9 me c from the saturated Buneman field over a Larmor time (Eqs. 5-6); beam-energy
  estimate 0.1 me c (Eq. 7); observed px spread ~0.1 me c (Fig. 5, Sect. 3.1).
- Compression 3.86 (forward) and 4.02 (reverse), as the Γ = 5/3 jump conditions (Fig. 6); overshoot nov/nup ≈ 8,
  rms field ≈ 8.5 B0 (Sect. 3.1).
- Reformation period ≈ 1.5 Ωi⁻¹ (forward, Fig. 18) and ≈ 1.4 Ωi⁻¹ (reverse) (Sect. 3.4).
- Ripples at the reverse shock: wavelength ~150 λse ≈ 20 λsi (100-300 λse over time), amplitude 2-5 λsi, speed
  0.08c ≈ 18 vA,L along the surface, changes within 0.1 Ωi⁻¹ (Sect. 3.2.2).
- Gradient-drift estimate |vgd| ≲ 0.075c · m/mi (Eq. 9); the measured bulk drift equals the E×B drift (Figs. 10, 15).
- Downstream Te/Ti ≃ 0.28 from momentum spreads (Eq. 10), Ti/Te ≈ 2 from Maxwellian fits (Sect. 3.3). Electron
  spectra: Maxwellian (reverse), "marginal" excess (forward) (Figs. 16-17); steady for t ≈ 5-20 Ωi⁻¹.

## Assessment
- Shown by simulation: Weibel-mediated perpendicular shocks at MA ≈ 28; cyclic reformation driven by unsteady ion
  reflection and filament merging (Fig. 19); ripples matching the Burgess & Scholer (2007) reflected-ion mechanism;
  bulk drift is E×B, with a small electric field in the drift frame (summary item 4); ion tail but no electron
  tail up to 20 Ωi⁻¹; no chains of magnetic islands (Sect. 3.2.1).
- Stated or conjectured: the 45° field approximates 3-D; β and mi/me are "not decisive" for electron acceleration
  (argued from other papers' runs, not from a scan); the missing electron tail points to lower efficiency in 3-D
  than in out-of-plane 2-D; ripples should enhance local Buneman heating (analysis deferred, Sect. 3.2.2).
- Strictly perpendicular shocks are superluminal, so no de Hoffmann-Teller frame exists and injection into shock
  drift acceleration should be suppressed (Sect. 3.2.1, citing Ball & Melrose 2001). Oblique shocks are left to a
  future paper (Sect. 1).
- Limits stated: one run; mi/me = 50; shock speeds a few times above SNR values (Sect. 4); reverse-shock foot poorly
  resolved (Buneman wavelength a few cells); forward-shock ripples (380-630 λse) do not fit in Ly; some upstream
  electron heating may come from the initial collision (Sect. 3.2.1).
- Disagreements: Matsumoto et al. 2012 (out-of-plane, similar MA and mi/me) find efficient surfing, this run none;
  their trapping estimate is called too optimistic (Sect. 3.3). Matsumoto et al. 2015 report turbulent reconnection
  (in-plane), absent here. Bohdan et al. 2017 later find that 45° behaves like the in-plane case, see reconnection
  at 45° for mi/me = 100, and attribute the missing tail partly to spectra taken far downstream.
