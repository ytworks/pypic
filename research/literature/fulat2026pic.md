# Fulat, Moore, Alawashra, Tsirou, Bohdan, Amano & Pohl (2026): PIC simulations of nonrelativistic high-Mach-number oblique shocks propagating in a turbulent medium

- Ref: accepted for publication in ApJ (arXiv comment; volume and page not in meta.txt, journal reference not
  checked); arXiv:2605.21986 (v1, read here); no DOI in the arXiv metadata; peer-reviewed according to the comment
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
At an oblique (θBn = 60°) high-Mach-number shock, compressive upstream turbulence with 15 % density fluctuations
changes the electron foreshock: whistler-driven cavities appear earlier and grow larger, reflected electrons are
hotter and the foreshock is shorter. After 20 Ωi⁻¹ the downstream nonthermal electrons are more numerous, carry
more energy and reach higher energy than with a homogeneous upstream medium; the slope of the tail is unchanged.

## Method
- PIC, 2D3V, code THATMPI; reflecting wall; turbulent slabs prepared and matched as in Fulat et al. (2023).
- θBn = 60°, B0 out of the simulation plane: B0 = B0 (cos θ0, 0, sin θ0).
- Table 1: mi/me = 100; MA = 32; Ms = 36; βp = 1; vsh = 0.264c (upstream frame); vth,e = 0.057c; Te ≈ Ti ≈ 10⁷ K;
  20 particles per cell per species; λse = 20Δ. ωpe/Ωce not stated (12 follows from Ωi⁻¹ ≈ 48,000 δt, my arithmetic).
- Box width Ly ≈ 12 λsi; run time about 20 Ωi⁻¹; shock at x ≈ 177 λsi at the end.
- Runs: H (homogeneous) and T (δn/n0 = 15 % rms in tiles of λsi/4, with δB/B0 ≈ 10 % and bulk-velocity fluctuations).
- Whistlers are identified by circular polarisation in the upstream frame (Eqs. 1, 2; App. A). Reflected
  electrons are separated from the inflow by removing a 3σ sphere around a fitted drifting Maxwellian (Sect. 3.3).

## Numbers to compare against
- Cavities appear at ≈ 9 Ωi⁻¹ (H) and ≈ 7 Ωi⁻¹ (T). Largest size 7 λsi (T) against "roughly 2 λsi" (H) in Sect. 3.1,
  but "3 λsi" for H in the summary (item 6).
- Example cavity in run T: R ≈ 2 λsi, field about twice B0, electrons inside with p⊥ ≈ 4 mec and γ ≈ 4 (Fig. 2).
- Field amplification starts about 30 λsi ahead of the front in both runs and is about 25 % stronger in T
  (Sect. 3.1, Fig. 3). Reformation period ≈ 1.6 Ωi⁻¹ after a start-up phase of about 7 Ωi⁻¹ (Fig. 4).
- Whistler growth rate from |BR|² ∝ exp[−2γ(x − xsh)/vsh] (Eq. 3) between 20 and 90 λsi ahead of the shock:
  γH = (8.16 ± 0.33)·10⁻³ Ωe, γT = (5.8 ± 0.58)·10⁻³ Ωe (2σ). Amplitude at 20 λsi: |BR|²/B0² ≈ 0.04 and 0.1, the
  runs not named, presumably H and T (Sect. 3.2, Fig. 5). Background turbulence: BR² ≈ BL² ≈ 0.01 B0² (App. B).
- Foreshock length at 20 Ωi⁻¹: 670 λsi (H), 540 λsi (T), described as "25 %" or about five upstream ion gyroradii
  closer (Sect. 3.3, Fig. 6). The return-current drift levels off at x ≲ 400 λsi (H) and ≲ 300 λsi (T) (Fig. 7).
- Reflected electrons: similar density gradient, larger perpendicular rms momentum in T; equal energy density
  within about 100 λsi of the front, more of it further out in H (Fig. 8).
- Downstream spectrum 20-90 λsi behind the shock at 20 Ωi⁻¹: electrons with γ − 1 ≥ 3 are 0.14 % (H) and 0.23 % (T)
  of all electrons and hold 1.9 % and 3.6 % of the electron energy; the largest Lorentz factor is 40 % higher
  in T; tail N(γ) ∝ (γ − 1)^−3.1 in both runs (Sect. 3.4, Fig. 9; summary item 8).
- Scaling argument: the exponent of Eq. (3) is −2(γ/Ωe)(mi/me)/MA · (x − xsh)/λsi, so the growth region shrinks in
  ion units as mi/me rises and becomes comparable to the shock thickness for the real mass ratio (Sect. 4).

## Assessment
- Shown by simulation: one pair of runs, compared at 20 Ωi⁻¹ (earlier runs of the group reached 50 Ωi⁻¹). The
  nonthermal number rises by a factor 1.6 and its energy share by 1.9. No error bars, seeds or time series of
  these fractions are given, and the foreshock is not in a steady state.
- The paper cannot separate two causes of the hotter reflected beam: scattering on the added magnetic
  fluctuations, or stronger front corrugation by the density clumps; it suggests both act (Sect. 3.3).
- Conjectured: that higher turbulence levels help stochastic shock-drift acceleration (summary item 8); that
  the final cavity size may be the same in both runs or be set by the box (Sect. 4).
- Stated limits: 2-D with B0 out of the plane suppresses the Weibel instability, so its interplay with whistlers
  and turbulence needs 3-D; electron-acoustic waves are not captured; the turbulence covers about one decade of
  wavenumbers in the dissipation range; its amplitude is capped by heating (Ms ≳ 30); only compressive turbulence.
- Differences from the perpendicular study (Fulat et al. 2023), which found no effect: obliquity, but also B0 out
  of plane instead of in plane, a box twice as wide, 15 % instead of 10 %, and twice the run time.
- The introduction quotes Morris et al. (2022) as "less than 0.1 %" of incoming electrons affected by
  electron-acoustic waves; that paper's conclusion says about 1 % for θBn ≤ 45°.
- Sect. 3.2 names "the pre-existing turbulence in run H"; run T must be meant. 670 → 540 λsi is 19 % of the
  longer foreshock; with an ion gyroradius MA λsi = 32 λsi it is four gyroradii (my arithmetic).
