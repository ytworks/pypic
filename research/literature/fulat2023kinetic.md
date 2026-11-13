# Fulat, Bohdan, Torralba Paz & Pohl (2023): Kinetic simulations of non-relativistic high-Mach-number perpendicular shocks propagating in a turbulent medium

- Ref: ApJ 959, 119 (2023) (volume/article as cited in Fulat et al. 2026; meta.txt has no journal_ref, comment
  "Accepted for publication in ApJ"); DOI 10.3847/1538-4357/ad04dc (same source); arXiv:2310.14388 (v1, read
  here); peer-reviewed
- Read: 2026-11-12, introduction, set-up and conclusions (results skimmed for numbers); topics: spherical-shock-injection

## Claim
Compressive upstream turbulence on kinetic scales, with density fluctuations of up to 10 % and clumps of a few ion
skin depths, changes neither the structure of a perpendicular high-Mach-number shock nor the heating and
pre-acceleration of electrons there. The paper also introduces the method that makes such runs possible.

## Method
- PIC, 2D3V, code THATMPI; reflecting wall, simulation in the downstream frame (Sect. 2).
- θBn = 90°, B0 in the simulation plane (B0 along y); mi/me = 100; MA ≈ 32, Ms ≈ 36 (from the shock speed 0.264c in
  the upstream frame); βp ≈ 1; Te ≈ Ti ≈ 1.6·10⁻³ mec²/kB; inflow 0.2c. ωpe/Ωce is not stated; Ωi⁻¹ ≈ 48,000 δt
  with δt = 1/40 ωpe⁻¹ implies 12 (my arithmetic).
- Box width Ly ≈ 6 λsi (the reference run B2 of Bohdan et al. 2019a had 24 λsi); λse = 20Δ; 20 particles per cell
  per species; run time 10 Ωi⁻¹ (Sect. 2.3).
- Runs: H (homogeneous), T1 (δn/n = 3.5 %), T2 (δn/n = 10 %); δn is the rms density in tiles of λsi/4 (Table 1).
- Turbulence: slabs of decaying compressive turbulence are prepared in periodic square boxes from 50 velocity
  modes with 2π/L < kx, ky < 6·2π/L and no imposed field perturbation, then appended upstream. Fields are
  interpolated in a matching zone, particles are re-drawn from interpolated moments, Gauss's law is restored with
  an added electric field and div B is cleaned by projection (Sect. 2.1, 2.2, App. A, B).

## Numbers to compare against
- Table 1 (H, T1, T2): vsh/c = 0.263±0.002, 0.265±0.003, 0.263±0.003; downstream kBTe/mec² = 0.168±0.004,
  0.165±0.006, 0.170±0.006.
- Reformation period 1.64 Ωi⁻¹ in all runs; compression 3.97; reflected ions reach at most 12 λsi ahead; field
  amplification 3.3 at the overshoot in all runs (Sect. 3.1, Figs. 2, 3).
- Largest density clumps ≈ 3 λsi against a shock transition of ≈ 15 λsi (Sect. 3.1).
- Buneman waves: electrostatic energy density 8.0, 5.5, 6.5 ·10⁻⁴ nemec² (H, T1, T2); peak |k|λse ≈ 4.1, 4.3, 4.0,
  as expected from kλse ≈ c/Δv with Δv ≈ 0.25c, 0.24c, 0.24c (Sect. 3.2, Eq. 1, Fig. 4).
- Weibel field energy and the number of magnetic islands: no significant difference (Sect. 3.3, Fig. 5).
- rms deviation of the local field angle ahead of the foot: 1.1° (T1), 2.3° (T2); the same in all runs inside
  the foot, ramp and overshoot (Sect. 3.5, Fig. 7).
- Upstream velocity fluctuations in T2 ≈ 0.2 % of vsh; density fluctuations reach ≈ 250 % at the overshoot and
  ≈ 10 % downstream in all runs (Sect. 3.6, Fig. 8).
- Downstream electron spectra of H and T2 are the same kappa distribution, κ = 13 (Sect. 3.7, Fig. 9).
- To keep Ms ≳ 30 the density fluctuations must stay near 10 % (Sect. 4, item 2).

## Assessment
- Shown by simulation: a null result for one geometry (perpendicular, in plane) and one kind of turbulence
  (compressive, kinetic scale, weak magnetic part), within the cycle-to-cycle scatter of three reformation cycles.
- Stated limits: the box is too narrow for front corrugation or ripples; clumps are smaller than the shock
  width; the amplitude is capped by heating and the cap "should be even lower" at realistic shock speeds; Alfvénic
  turbulence is left for later; ion temperature is not discussed; 2-D only.
- Conjectured: turbulence on scales above the shock width may matter (hybrid results of Trotta et al. 2021 and
  Nakanotani et al. 2022 are cited); the effect should be larger at lower MA, where Weibel modes are weaker, and at
  oblique shocks, where a foreshock leaves room for interaction (taken up in Fulat et al. 2026).
- The weak tail (κ = 13) is "not sufficient for injection" in all three runs; no nonthermal fraction is tabulated.
- Contrast drawn by the authors: test-particle work (Guo & Giacalone 2010, 2012) found that magnetic upstream
  turbulence raises electron acceleration; pair-plasma PIC (Demidem et al. 2023; Bresci et al. 2023) found the same
  for relativistic shocks.
- The introduction divides shocks into perpendicular (90°), oblique (about 50°-75°) and quasi-parallel (0°).
