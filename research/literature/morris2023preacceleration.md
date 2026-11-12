# Morris, Bohdan, Weidl, Tsirou, Fulat & Pohl (2023): Pre-acceleration in the Electron Foreshock II: Oblique Whistler Waves

- Ref: ApJ 944, 13 (2023) (volume/article as cited in Fulat et al. 2026; meta.txt has no journal_ref);
  DOI 10.3847/1538-4357/acaec8; arXiv:2301.00872 (v1, read here); peer-reviewed
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
In the inner electron foreshock of an oblique high-Mach-number shock, reflected electrons drive oblique whistler
waves that grow into non-linear magnetic cavities of ion scale. The cavities leave incoming thermal electrons almost
untouched but scatter and trap reflected electrons and carry them back to the shock, where they gain energy by a
process like stochastic shock-drift acceleration (SSDA), about three times faster than the motional field explains.

## Method
- PIC, 2D3V, code derived from TRISTAN (Vay pusher), reflecting wall; the same run as in Bohdan et al. (2022).
- θBn = 60°, B0 out of the simulation plane (φ = 90°): B0 = B0 (0.5, 0, √3/2).
- mi/me = 50; MA = 30, Ms = 32.5, both with the shock speed 0.263c in the upstream frame (0.067c in the simulation
  frame); β = 1; kTe = kTi = 9.86·10⁻⁴ mec²; Ωce = 0.06 ωpe; inflow 0.2c.
- 40 particles per cell per species; λse = 8Δ; δt = 1/16 ωpe⁻¹; transverse box 41.6 λsi; run time 51.5 Ωci⁻¹; the
  box reaches about 2000 λsi at 49.9 Ωci⁻¹ (Sects. 2, 3).
- 10,000 upstream electrons chosen at random at 30 Ωci⁻¹ and traced to the end (Sect. 4.2).

## Numbers to compare against
- Whistlers appear at about 10 Ωci⁻¹. The reflection rate (measured 2-10 λsi ahead of the shock) drops after that
  and settles near 4 % for t > 30 Ωci⁻¹; the mean energy and the energy density of the reflected beam grow about
  linearly in time (Sect. 3.1, Fig. 3). Reflected ion density about 10⁻⁴ n0.
- Linear stage (from Bohdan et al. 2022): fastest growth at k⊥λse ≈ 0.2; phase and group speeds of order 10⁻³ c,
  so the waves move with the upstream flow (Sect. 4.1).
- Cavities: radii 1-5 λsi, density down to 0.15 n0, field up to about 4 B0 in a ring at the edge (Sect. 4.2,
  Figs. 4, 5). Gyroradii of reflected electrons 1-5 λsi, up to about 10 λsi (Sect. 4.4).
- Interaction probability with a cavity, Eq. (7); it needs v∥' cos θBn' ≥ vsh, hence v⊥' ≤ 0.848c. Trapping is most
  likely for v⊥' ≈ 0.40c, v∥' ≈ 0.54c, about 18 times the upstream thermal energy (Sect. 4.5).
- Extrapolation from Bohdan et al. (2022): the whistler region becomes steady at about 125 Ωci⁻¹ with a length of
  about 2000 λsi (Sect. 4.5).
- Energy gain against drift distance: the most energetic electrons gain 3-4 times qE0yΔy (Sect. 4.6, Fig. 9). For
  the electron of Fig. 7, Wz ≈ 2.1 Wy and Wx ≈ 0.17 Wy; the local field is tilted by about 26° from B0.
- Of the 70 traced electrons with γ ≥ 5, 20 % have a pitch angle above π/3 and are caught by the shock again; with
  R ≈ 4 % this gives 0.8 % of the incoming electrons, called an approximate upper limit (Sect. 4.7, Fig. 10).
- Downstream: Te/Ti ≈ 0.17; power law with index up to −2.8; 0.8 % of the electrons lie outside the Maxwellian fit
  and hold about 7 % of the electron energy; the most energetic reach γ ≈ (mi/me)(vsh/c) ≈ 13, the threshold for
  diffusive shock acceleration (DSA) being "a factor of a few" higher (Sect. 4.8, Fig. 11).

## Assessment
- Shown by simulation: the cavity structure; the sign of the Lorentz force for a test electron moving with the
  flow (outward) and for a reflected one (inward); one fully analysed trapped electron (Figs. 7, 8); the excess of
  the energy gain over the SSDA estimate for the traced sample (Fig. 9).
- The authors' verdict: the cavities do not energise electrons themselves; their "primary role" is to keep
  reflected electrons near the shock. The gain occurs at the shock front, where the cavities also perturb the fields.
- Stated or estimated: the 0.8 % (a pitch-angle argument, not a count of returned electrons); the trapping-probability
  model Ptrap ∝ Pinteract exp(−rge/rnls); the origin of the factor 3-4 "requires further investigation".
- Limits stated: DSA and a nonthermal ion tail are not reached; 2-D out of plane, so the electron-acoustic waves
  of Paper I are poorly captured; one parameter set; mi/me = 50 and vsh = 0.263c.
- Sect. 4.4 gives the escape condition as vx > vsh cos θBn = vsh/2, Sect. 4.5 as v∥' cos θBn' ≥ vsh (which matches
  v∥,min' = 0.53c); the first looks like a slip.
- Bohdan (2023, Sect. 4.3) quotes "up to 7 % of electrons" for this mechanism; here 7 % is the energy share.
- Fulat et al. (2026), with mi/me = 100, a 12 λsi box and 20 Ωi⁻¹, find 0.14 % of electrons above γ − 1 = 3 and
  an index −3.1 for a homogeneous upstream medium; the definitions of "nonthermal" differ.
- Credits: SDA Wu (1984), Leroy & Mangeney (1984); SSDA Katou & Amano (2019) and Matsumoto et al. (2017);
  limits of adiabatic SDA Vandas (2001); whistler scattering at the bow shock Oka et al. (2017, 2019).
