# Amano, Matsumoto, Bohdan, Kobzar, Matsukiyo, Oka, Niemiec, Pohl & Hoshino (2022): Nonthermal Electron Acceleration at Collisionless Quasi-perpendicular Shocks

- Ref: Reviews of Modern Plasma Physics 6, 29 (2022) (as cited in Morris et al. 2023; meta.txt has no journal_ref,
  comment "To appear in Reviews of Modern Plasma Physics as an invited review"); DOI 10.1007/s41614-022-00093-1
  (as cited in Fulat et al. 2026); arXiv:2209.03521 (v1, read here); invited review, peer-reviewed
- Read: 2026-11-12, sections 1-4, 5.2-5.4 and 6 (5.1 skimmed); topics: spherical-shock-injection

## Claim
For quasi-perpendicular shocks (θBn ≳ 45°) the review argues that stochastic shock-drift acceleration (SSDA) is the
most probable route of electron injection into diffusive shock acceleration (DSA), that shock surfing (SSA) works
in 3-D as a fast pre-acceleration ahead of it, and that shocks with Alfvén and sound Mach numbers above about
20-40 are dominated by the ion Weibel instability, with reconnection as a further channel. Which mechanism
dominates for given shock parameters is left open.

## Method
- Review of theory, in-situ data and PIC runs, mainly the authors' own. Mach numbers are in the normal incidence
  frame (shock at rest); the theory uses the de Hoffmann-Teller value MA^HTF = MA/cos θBn (Sect. 2.1).
- Runs quoted with parameters: 3-D, Matsumoto et al. (2017): MA = 20.8, θBn = 74°, β = 0.5 per species,
  ωpe/Ωce = 10, mi/me = 64 (Sect. 3.5); 2-D, Kobzar et al. (2021): MA = 6.1, θBn = 75°, β = 2.5 per species,
  mi/me = 100; 2-D perpendicular, Matsumoto et al. (2012): MA = 30, mi/me = 100 (Sect. 4.4); 1-D perpendicular,
  Hoshino & Shimada (2002): MA = 32, mi/me = 20 (Sect. 4.3). Model solutions: MA = 10, θBn = 85° (Fig. 8).

## Numbers to compare against
- Mach numbers: young SNR shocks ~10-100, heliospheric shocks ~2-10 (Sect. 1). Shock speed/c: ≳ 0.1 in typical PIC
  runs, ~0.01 for young SNRs, ~0.001 for planetary bow shocks (Sect. 5.2). Superluminal only for θBn > 89° at
  young SNR speeds (Sect. 2.1). Whistler critical Mach number MA^w* ≈ (mi/me)^1/2 cos θBn / 2 (Sect. 2.2, Fig. 1).
- Injection thresholds: p/(meVA) ≈ (mi/me) MA, ~10 MeV at SNRs (Eq. 2, gyroradius = shock width), against
  pinj/(meVA) ≈ (mi/me)(k*c/ωpi)^−1, ~0.1-1 MeV (Eq. 3, gyroradius = dissipation scale of the turbulence) (Sect. 2.5).
- Shock-drift acceleration: reflection if the pitch angle exceeds asin √(B0/Bm); gain Δp ~ 2meUs/cos θBn; only
  electrons with v ≳ Us/(2 cos θBn) are reflected (for Bm/B0 ~ 4); time scale ~Ωci⁻¹ (Sect. 3.2, Eq. 5).
- SSDA: needs diffusion length ≲ shock thickness (Eq. 7); works between pmin ≈ meVA MA^HTF (Eq. 12) and pmax
  (Eq. 13); spectral index q ≳ 6 against 4 for DSA (Eq. 8); bow-shock spectrum cut off near 20 keV (Fig. 4).
  Injection into DSA needs MA^HTF above the bound of Eq. (14); acceleration is expected for
  MA^HTF ≳ (mi/me)^1/2 βe^χ, −1/2 ≤ χ ≤ 1/2 (Eq. 15); one reflection suffices for MA^HTF ≳ (mi/me)(k*c/ωpi)^−1 (Eq. 16).
- SSA: Buneman instability for MA ≳ (1+αr)/2 (mi/me)^1/2 βe^1/2, αr ~ 0.2 (Eq. 21); efficient for
  MA ≳ (1+αr)(mi/me)^2/3 (Eq. 22); time scale ~Ωce⁻¹ (Eq. 19); gain per reflection ~40 times the upstream electron
  flow energy, relativistic only if Us/c ≳ 0.1 (Sect. 4.5).
- Weibel-dominated regime and reconnection: MA ≳ 10 (nr/2n0)^−1/2 ≈ 32 for nr/n0 = 0.2 (Eq. 26); growth rate
  γmax/Ωci ≈ (nr/2n0)^1/2 MA (Eq. 25), insensitive to mi/me and ωpe/Ωce (Sect. 5.2, Fig. 17).

## Assessment
- Supported by observation: SSDA at Earth's bow shock (Amano et al. 2020; Sect. 3.4). Supported by simulation:
  SSA followed by SSDA in the one 3-D run, scattering by Weibel turbulence (Sects. 3.5, 4.4); reconnection only in
  2-D runs with an in-plane field, not found in 3-D (Sect. 5.4).
- Theory with unknown inputs: the scattering rate inside and outside the shock sets every SSDA threshold; the
  exponent χ and the wave source (Levinson 1992; Amano & Hoshino 2010; modified two-stream) are undecided (Sect. 3.7).
- Heliosphere against SNRs: heliospheric shocks accelerate electrons poorly (a few hundred keV at most; Sect. 2.5);
  the contrast is expected to follow from the Mach-number thresholds; they are rarely Weibel-dominated (Sect. 5.2).
- Geometry: θBn enters only through MA^HTF, which favours nearly perpendicular shocks; no PIC run of a curved
  shock is discussed. The curved bow shock appears only in observations: the most energetic upstream electrons
  come from where the field is tangent to it (Anderson et al. 1979; Sect. 3.1). Rippling helps electrons in
  test-particle runs on hybrid fields (Burgess 2006; Sect. 3.5) and is the low-Mach form of the Weibel mode (5.2).
- Numerical limits named: reduced mi/me lowers the SSA threshold (Eq. 22); high Us/c makes one SSA reflection
  relativistic; 2-D favours SSA (field out of plane) or reconnection (in plane), never both; a low mass ratio
  may favour drift-kink over tearing (Sect. 5.4); 1-D SSA is unrealistic (Sect. 4.4).
- Open questions listed: relative roles of SSDA, SSA and reconnection; wave generation and scattering rates; the
  exact condition for reconnection and the fate of current sheets in 3-D; whether SSDA survives in a
  Weibel-dominated shock; scaling of SSA with shock speed and mass ratio; multiple SSA cycles; βe dependence.
- Not covered: quasi-parallel shocks (said to contain locally quasi-perpendicular portions, Sect. 1) and the
  foreshock waves of Bohdan et al. (2022). Levinson's threshold (Sect. 3.7) is partly illegible in the text file.
