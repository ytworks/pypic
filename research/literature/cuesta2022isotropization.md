# Cuesta, Chhiber, Roy, Goodwill, Pecora, Jarosik, Matthaeus, Parashar & Bandyopadhyay (2022): Isotropization and Evolution of Energy-containing Eddies in Solar Wind Turbulence: Parker Solar Probe, Helios 1, ACE, WIND, and Voyager 1

- Ref: Astrophys. J. Lett. 932, L11 (2022); DOI 10.3847/2041-8213/ac73fd; arXiv:2205.00526; peer-reviewed.
- Read: 2026-11-13, full (arXiv v1 PDF text; figures not viewed); topics: shock-alfvenic-turbulence

## Claim
Magnetic correlation lengths grow by two orders of magnitude from ~5×10⁴ km at 0.1 au to ~3×10⁶ km at 10 au.
Near the Sun λ∥ < λ⊥ (≈ 0.75 at 0.1 au); the two roughly isotropise within ~0.4–1 au and reverse to
λ∥/λ⊥ ≈ 1.29 at 1 au (ACE/WIND).

## Method
- PSP orbits 1–8, Helios 1, ACE, WIND, Voyager 1; 0.08–10 au.
- λ_C = V_sw τ_e with τ_e the e-folding time of R(τ) = ⟨b(t)·b(t+τ)⟩/⟨b·b⟩ (Eq. 1); ACE values from an
  exponential fit (Roy et al. 2021). Taylor hypothesis; ~1 % of PSP intervals (R < 0.2 au) have V_A/V_sw > 0.66.
- Parallel channel 0° < Θ_BV < 40°, perpendicular 50° < Θ_BV < 90°.
- Interval lengths: 1 h (R < 0.3 au), 3 h (0.3–1 au, PSP), 1 day (ACE/WIND, Voyager).

## Numbers to compare against
- Table 1 (ACE+WIND, 1 au, mean ± s.e.): all wind λ∥ = 1.98 ± 0.09, λ⊥ = 1.55 ± 0.02 (10⁶ km), ratio 1.28;
  slow (< 450 km/s) 1.98 / 1.67, ratio 1.19; fast (> 600 km/s) 1.78 / 1.05, ratio 1.70.
- Table 2 (λ ∝ R^α): λ_C α = 0.97 (R < 0.3 au), 0.29 (0.3–1 au), 0.27 (> 1 au); λ∥ 1.03 / 0.28 / 0.64;
  λ⊥ 0.61 / 0.27 / 0.23.
- Quoted: Dasso et al. (2005) λ∥/λ⊥ = 0.71 (fast) and 1.18 (slow) at 1 au.
- Derived (my arithmetic, d_i ≈ 100 km at 1 au): λ_C ≈ 1.5–2×10⁶ km ≈ 1.5–2×10⁴ d_i at 1 au (1-day
  intervals). At 0.1 au, ~5×10⁴ km; with d_i ≈ 33 km at 0.24 au (Jebaraj et al. 2024) the ratio is also ~10³–10⁴.

## Assessment
- Shown: the outer scale is ≳ 10⁴ d_i at 1 au and ≳ 10³ d_i near the Sun; the correlation scale is never within
  reach of a fully kinetic box, so PIC "correlation lengths" of 10–100 d_i represent only the inertial-range
  part of the turbulence.
- The values depend on the interval length (jumps at 0.3 au and 1 au in Fig. 2; Isaacs et al. 2015), so the
  1 au numbers with 1-day intervals are upper-side estimates.
- Fast-wind results disagree with Dasso et al. (2005); method differences are offered as the reason.
- Single-spacecraft, Taylor-hypothesis estimates; not the two-point method of Matthaeus et al. (2005).
