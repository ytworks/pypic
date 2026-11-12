# Electron injection at supernova remnant shocks: survey of prior work

Topic `spherical-shock-injection`. Written 2026-11-12. Status: first survey, not yet reviewed.

**Scope.** The electron injection problem at the shocks of supernova remnants (SNRs): what observations
require, which mechanisms have been proposed, what kinetic simulations have shown, and what is open. Heliospheric
shocks appear only where they serve as a benchmark for SNR shocks.

**How to read it.**
- A name such as `park2015simultaneous` is a citekey. Its note is `research/literature/<citekey>.md`, and its full
  reference is in Sect. 10.1. Locations (Sect., Fig., Eq., Table) refer to the arXiv version that was read.
- "Not read" means the work was not read for this survey. Such a work is cited through a work that was read, and
  its reference is in Sect. 10.2 as that work prints it.
- "Our arithmetic" marks a value derived here from numbers in the papers. It is not a statement of the paper.
- Abbreviations: θ_Bn, angle between the shock normal and the upstream field; M_A and M_s, Alfvén and sonic Mach
  numbers; SSA, shock surfing acceleration; SDA, shock drift acceleration; SSDA, stochastic shock drift
  acceleration; DSA, diffusive shock acceleration; K_ep, electron-to-proton ratio of the accelerated particles;
  PIC, particle-in-cell.

## 1. Summary

1. **The problem.** DSA accelerates only particles that can cross the shock. For thermal electrons the gap to that
   threshold is 10⁵–10⁶ in energy at a real SNR shock (`bohdan2023electron` Eq. 1). How electrons cross the gap is
   the injection problem. The works do not agree on what the threshold is (Sect. 2).
2. **Observations** give the result of injection, not the process. Young SNRs hold 10–100 TeV electrons
   (`vink2012supernova` Sect. 6.3.1). Models of Tycho and SN 1006 need K_ep ≈ 10⁻³ (`morlino2012strong` Sect. 3;
   `mcginness2026sn` Table 1). No observation measures the injected fraction of electrons.
3. **Obliquity in remnants is contested.** In SN 1006 most works place the efficient acceleration where the shock
   is quasi-parallel (the polar caps). The radio sample of bilateral remnants favours quasi-perpendicular
   acceleration (`west2017connection` Sect. 5). Every global model takes its efficiency law ε(θ_Bn) as an
   assumption (Sect. 3.3).
4. **Mechanisms depend on obliquity** (Sect. 4):
   - quasi-perpendicular: surfing on Buneman waves, shock drift with scattering (SSDA), whistlers, reconnection
     and stochastic scattering in Weibel turbulence;
   - oblique: an electron foreshock with electron-acoustic and whistler waves;
   - quasi-parallel: mirror reflection, then scattering between the shock and proton-driven waves.

   Electron spectra consistent with DSA have been reached in 1-D (`park2015simultaneous`; `xu2020electron`;
   `kumar2021nonthermal`; `gupta2024electron`), in a thin 2-D box (`arbutina2021nonlinear`) and in 2-D at a mildly
   relativistic speed (`crumley2019kinetic`). None of these is a quasi-perpendicular shock in more than one
   dimension.
5. **No kinetic simulation covers the range of θ_Bn with one set of parameters**, and none found has a front that
   is curved on the large scale at SNR-like parameters (Sect. 6). Curved fronts were treated for the Earth's bow
   shock only, in works that could not be read here.
6. **Numerical limits shape the answers** (Sect. 7): reduced mass ratio, shock speeds of 0.03–0.4 c, and 1-D or
   2-D geometry. The field orientation in 2-D decides which instability exists. The thresholds of the mechanisms
   scale with different powers of the mass ratio.
7. **Open** (Sects. 8.2, 8.3): the form of ε(θ_Bn); which mechanism dominates at the real mass ratio and shock
   speed; whether quasi-perpendicular or quasi-parallel regions inject more electrons; the transition from
   pre-acceleration to DSA in more than one dimension.

## 2. The injection problem in numbers

**What "injected" means.** The works use different thresholds:

| Threshold | Statement | Source |
|---|---|---|
| Gyroradius ≥ shock width | gap of 10⁵–10⁶ in energy for real shocks; 10³ in an oblique PIC run | `bohdan2023electron` Eq. 1, Sect. 4.3 |
| Momentum | p_e ≫ m_i V_sh: about 10 MeV/c at V_sh = 0.01 c | `pohl2020pic` Eq. 34; `bohdan2019kineticb` Sect. 4 |
| Diffusion length ≲ shock thickness | the gyroradius condition is "the most pessimistic"; injection at 0.1–1 MeV | `amano2022nonthermal` Sect. 2.5, Eq. 3 |
| Proton injection momentum | electrons join DSA near p_inj ≈ 3 m_p v_u of the protons | `park2015simultaneous` Fig. 3 |
| Speed | an electron is injected once its speed exceeds v_sh, whatever its gyroradius | `gupta2025speed` Sect. II B, Fig. 2 |
| Scale of upstream structures | gyroradius ≈ λ/2π of the steepened waves (SLAMS) at high M_A | `zekovic2025slams` Eq. 1 |
| Multiple of the thermal momentum | p_inj = ξ p_th with ξ ≈ 3.2–3.7, the same for both species | `arbutina2021nonlinear` Eqs. 4, 12, Sect. IV |

**Parameters of SNR shocks** that a simulation would have to match:
- Shock speed 2000–6000 km/s in the remnants with X-ray synchrotron rims (`vink2012supernova` Sect. 11.3), that
  is 0.007–0.02 c (our arithmetic). SN 1006: 5000 ± 400 km/s, 3000 ± 400 km/s in the north-west
  (`mcginness2026sn` Sect. 2.1). Tycho: 4990 km/s in the model (`morlino2012strong` Sect. 2.1).
- Field: 3–5 μG in the interstellar medium (`reynoso2013radio` Sect. 4; `vink2012supernova` Sect. 11.3) and
  100–600 μG at the rims (`vink2012supernova` Sect. 11.2).
- Density 0.02–0.3 cm⁻³ (`mcginness2026sn` Sect. 2.1; `morlino2012strong` Sect. 2.1).
- Mach numbers: only `morlino2012strong` gives them, M_s ≃ 443 and M_A ≃ 15 in the amplified field (Sect. 2.3).
  With the unamplified 5 μG and 0.3 cm⁻³ the far-upstream values are M_A ≈ 250, ω_pe/Ω_ce ≈ 350 and β ≈ 0.8
  (our arithmetic).
- A shock is superluminal, so that no electron can outrun it along the field, only for θ_Bn > 89°
  (`amano2022nonthermal` Sect. 2.1).
- Radius over shock thickness is about 10¹² (our arithmetic from `vink2012supernova` Sect. 5.2).

**Parameters of the simulations** (Sect. 5): m_i/m_e = 16–400 in most; shock speed 0.03–0.4 c; M_A = 5–70 in
most, up to 150 in 2-D (`bohdan2021magnetic`) and 300 in 1-D (`zekovic2025slams`). At 0.26 c a shock is superluminal above θ_Bn ≈ 75° (`morris2022preacceleration`
Sect. 2). The superluminal part of a sphere in a uniform field has the area fraction v_sh/c: 26 % at 0.26 c
against under 2 % at 5000 km/s (our arithmetic).

## 3. Observational constraints

### 3.1 Electrons are accelerated at SNR shocks
- Thin X-ray synchrotron rims show 10–100 TeV electrons at the forward shocks of young remnants
  (`vink2012supernova` Sect. 6.3.1). The first synchrotron shell was found in SN 1006 (Koyama et al. 1995, not
  read; cited by `vink2012supernova` and `reynolds2004microphysics`).
- The cut-off photon energy is hν = 1.4 η⁻¹ (V_s / 5000 km/s)² keV, with η the diffusion coefficient in units of
  the Bohm value (`vink2012supernova` Eq. 52).
- Thermal electrons stay cool: kT_e ≤ 5 keV in all young remnants; in SN 1006 kT_e ≈ 1.5 keV against
  528 ± 150 keV for O VII (`vink2012supernova` Sect. 11.1).
- `vink2012supernova` (Sect. 5.5) calls injection an "outstanding uncertainty".

### 3.2 How many electrons
- Tycho: K_ep = 1.6 × 10⁻³ from the radio and X-ray synchrotron emission, a free normalisation of the model
  (`morlino2012strong` Sect. 3). The paper states that a theory of electron injection is "still lacking".
- SN 1006: K_ep = 9.3 × 10⁻⁴ in all four quadrants of a quadrant-by-quadrant model (`mcginness2026sn` Table 1).
  Electrons are tied to protons by construction there. The same paper quotes a range of 10⁻²–10⁻⁴ for K_ep in
  the literature (Sect. 2.2).
- At Earth K_ep ≈ 0.01, and young SNRs suggest 10⁻³ or less (`park2015simultaneous`, Introduction).
- K_ep is degenerate with the magnetic field in these fits (reading of `morlino2012strong` by this survey).
- No paper measures the fraction of thermal electrons that is injected.

### 3.3 Where on the shell: obliquity
A spherical remnant in a uniform field has every θ_Bn. SN 1006 has two bright limbs, and the question is whether
they are polar caps (quasi-parallel) or an equatorial belt (quasi-perpendicular).

**How models put obliquity in.** None of these laws is derived from kinetic theory:

| Law | Form | Source |
|---|---|---|
| quasi-parallel / quasi-perpendicular / isotropic | ε ∝ cos²Θ, sin²Θ, constant; Θ is taken to the post-shock field | `bocchino2011constraints` Sect. 2; `west2017connection` Sect. 2; origin: Fulbright & Reynolds 1990, Leckband et al. 1989 (not read) |
| obliquity in the acceleration rate | κ = κ_∥ [cos²θ_Bn + sin²θ_Bn / (1 + η²)], faster at perpendicular shocks | `reynolds2004microphysics` Eq. 3; Jokipii 1987 (not read) |
| step for protons | ζ_p = (0.15/2) [tanh((θ_crit − θ)/δ) + 1], θ_crit = 45°, δ = 10° | `winner2020evolution` Eq. 1, after `caprioli2014simulations` |
| step for electrons, three trials | the same tanh between ζ_∥ and ζ_⊥ with ζ_∥/ζ_⊥ = 10, 1, 0.1 | `winner2020evolution` Eq. 4, Fig. 1 |
| pressure of accelerated particles | tanh steps at 45° (thermal protons, 12 %) and 70° (re-accelerated seeds, 6 %) | `giuffrida2022supernova` Methods, Eq. 1 |
| per quadrant, fitted | efficiency 20 % in the caps, 0.4–0.8 % elsewhere, K_ep constant | `mcginness2026sn` Table 1 |

**Evidence.**

| Work | Data or method | Finding | Its own caveat |
|---|---|---|---|
| Rothenflug et al. 2004 (not read) | X-ray centre-to-rim ratio of SN 1006 | excludes an equatorial belt | via `bocchino2011constraints` Sect. 1 |
| `reynolds2004microphysics` | synthetic images with obliquity-dependent diffusion | bright rims "almost certainly perpendicular" (Conclusions) | mean field 10 μG assumed |
| `bocchino2011constraints` | 3-D MHD and radio maps of SN 1006 | quasi-parallel fits best with the field at 38° ± 4° to the line of sight and a gradient of the field strength (Sect. 3.4.1) | in a uniform field quasi-parallel cannot fit (Sect. 3.3); quasi-perpendicular is "worse", not excluded (Sect. 3.4.2) |
| `reynoso2013radio` | radio polarimetry of SN 1006 | bright lobes: radial, disordered field, polarisation 0.17; south-east rim: tangential, ordered field, 0.6 ± 0.2 (Sects. 3.1, 3.2) | calls the gradient model "somewhat contrived" (Sect. 1) |
| `zhou2023magnetic` | X-ray polarimetry of the north-east limb | downstream field radial; polarisation degree 22.4 ± 3.5 % (Table 1) | "consistent with" quasi-parallel; radial and uniform fields cannot be separated (Sect. 3.2) |
| `giuffrida2022supernova` | post-shock density against azimuth | compression rises from 4 to about 7 towards θ ≈ 0° (Fig. 3b) | about hadrons; no electron number |
| `winner2020evolution` | 3-D MHD with electron spectra | maps need electron efficiency at least 10 times higher at quasi-parallel shocks (Sects. 3.1, 7) | three laws tried; the field model is imposed |
| `west2017connection` | 32 bilateral radio remnants in a Galactic field model | the quasi-perpendicular weight matches about 75 % (Sect. 5); quasi-parallel fails in at least 10 (Sect. 3) | comparison by eye; SN 1006 fits better with quasi-parallel (Sect. 4.2) |

**Disagreements.**
- SN 1006 is taken as quasi-parallel by most works, and as quasi-perpendicular by `reynolds2004microphysics` and
  the works it follows.
- The population of bilateral remnants points the other way (`west2017connection`). Most bilateral remnants have
  their axis along the Galactic plane, and SN 1006 is "rather peculiar" (`bocchino2011constraints` footnote 5).
- `winner2020evolution` (Sect. 7) states that its result "precludes extrapolation" of 1-D PIC runs that favour
  quasi-perpendicular shocks for electrons. Its quasi-perpendicular trial was motivated by exactly those runs
  (`xu2020electron`), against `park2015simultaneous` for quasi-parallel (Sect. 2.2).
- The field of SN 1006 differs between models: 140 μG (`winner2020evolution`) against 6 μG
  (`mcginness2026sn`) at the quasi-perpendicular rim.

**What stays unknown, by the papers' own words.** The form of the electron law is not known from first
principles (`winner2020evolution` Sect. 2.2). Electron injection "is not completely understood", and the
mechanism at quasi-perpendicular regions is open (`mcginness2026sn` Sects. 2.2, 3.2).

## 4. Mechanisms and what kinetic simulations show

`bohdan2023electron` (Sects. 2–5) sorts the field by obliquity: perpendicular (75–90°), oblique (50–75°),
parallel (below 50°). This section follows that order.

### 4.1 Quasi-perpendicular shocks at high Mach number

**Surfing on Buneman waves (SSA).** Reflected ions drift against incoming electrons in the foot and drive the
Buneman instability; electrons trapped in the waves gain energy from the motional field.
- Origin: Shimada & Hoshino 2000 and McClements et al. 2001 (not read); `hoshino2002nonthermal`.
- Conditions:
  - instability: M_A ≥ [(1+α)/2] (β_e m_i/m_e)^(1/2), with α ≈ 0.2 the reflected-ion fraction
    (`matsumoto2012electron` Eq. 5);
  - trapping: M_A ≥ (1+α) (m_i/m_e)^(2/3) (`matsumoto2012electron` Eq. 8). This is about 26 at m_i/m_e = 100
    and about 180 at 1836 (`wieland2016nonrelativistic` Sect. 2.2);
  - the earlier 1-D condition was M_A > 43–75 at the real mass ratio (`hoshino2002nonthermal` Eq. 8).
- Dimensionality decides the outcome:
  - 1-D: planar waves trap electrons (`amano2007electron`).
  - 2-D with the field out of the plane: surfing works, as a stochastic process (`amano2009electron`;
    `matsumoto2012electron`; `matsumoto2013electron`).
  - 2-D with the field in the plane: Buneman waves are suppressed; the pre-accelerated fraction is about 0.45 %
    against 6.8 % out of the plane (`bohdan2019kinetic` Table 2).
  - 3-D: surfing exists and feeds a second stage (`matsumoto2017electron`).
- Against it:
  - a 2-D electrostatic run at the real mass ratio finds heating only, because the wave potential is coherent
    over one wavelength (`ohira2007absence` Sects. 2.4, 3);
  - its maximum energy depends on the shock speed and not on the mass ratio, so it is negligible at 1836
    (`riquelme2011electron` Appendix);
  - the mean gain per electron is independent of the mass ratio while the downstream temperature grows with it,
    so the share of the tail owed to surfing falls from 0.15 % to 0.0015 % of electrons for m_i/m_e from 50
    to 400 (`bohdan2019kineticb` Sect. 3.3).

**Shock drift (SDA) and stochastic shock drift (SSDA).** An electron mirrored by the shock gains energy from
the motional field; with pitch-angle scattering inside the shock layer the process repeats.
- Origin: Wu 1984; Leroy & Mangeney 1984 (not read).
- One reflection gives Δv = 2 u₀ / cos θ_Bn (`katou2019theory` Eq. 1). The gain grows with θ_Bn while the number
  of reflected electrons falls (`amano2007electron` Sect. 1).
- With scattering: a power law N(v) ∝ v^(−p), p = 1 + 3L/L_sh ≈ 4 (`katou2019theory` Eq. 20), and a maximum
  energy ε_max ~ ½ m_i u_sh² D_μμ/Ω_ce with u_sh = u₀/cos θ_Bn (Eq. 21). That is about 11 keV at the Earth's
  bow shock (Eq. 34) and about 1 MeV at ten times its speed (Sect. 4).
- In-situ support: one bow-shock crossing with M_A ≃ 8.9 and θ_Bn ≃ 85° shows a cut-off at 22 ± 1 keV against
  about 20 keV predicted (`amano2020observational` p. 4, Fig. 4).
- Thresholds:
  - whistlers driven by the reflected electrons: M_A ≳ (cos θ_Bn / 2) (β_e m_i/m_e)^(1/2)
    (`amano2010critical` Eq. 2);
  - injection through SSDA: M_A / cos θ_Bn ≳ 55, for the scattering rate and thickness assumed there
    (`amano2022theory` Eq. 14);
  - wave generation by reflected electrons: M_A / cos θ_Bn ≳ (α/4) (m_i/m_e)^(1/2) ≈ 50
    (`wang2026generation` Eq. 22);
  - reflection weakens again at very high M_A / cos θ_Bn: the reflected density falls as
    exp[−(M_A/cos θ_Bn)² (m_e/m_i) / ((B₂/B₁ − 1) β_e)] (`wang2026generation` Eq. 24).
- 3-D PIC at θ_Bn = 74.3°: surfing first, then SSDA with scattering by Weibel turbulence; a power law of index
  3.5; electrons with γ > 5 escape upstream (`matsumoto2017electron` Figs. 3, 4).
- A reduced mass ratio makes this injection artificially easy (`amano2022theory` Eq. 21, Sect. 4.4).

**Whistler waves.** Oblique whistlers in the foot accelerate electrons along the field; the electrons escape
upstream and return (`riquelme2011electron`).
- Condition: M_A ≲ (m_i/m_e)^(1/2), tested up to m_i/m_e = 1600 (Sect. 3.1, Fig. 12).
- Outcome: 2–10 % of electrons by number and 10–30 % by energy in a tail of index 2.7–4 (Sect. 4.4).
- For SNRs this needs M_A ≲ 20, that is an amplified upstream field (Sect. 5).
- This contradicts the surfing condition, which needs M_A above (m_i/m_e)^(2/3).

**Weibel turbulence: reconnection and stochastic scattering.** At M_A of a few tens and above, the ion Weibel
instability shapes the shock.
- A perpendicular shock at M_A ~ 130 with the field in the plane is Weibel-mediated and shows no electron
  acceleration (`kato2010nonrelativistic`).
- Reconnection of the Weibel current sheets (Matsumoto et al. 2015, not read) supplies most of the tail at
  m_i/m_e ≥ 100: 0.4–0.55 % of electrons out of a nonthermal fraction of 0.55–0.7 % (`bohdan2020kinetic`
  Table 3). The rest is put down to stochastic Fermi-like scattering (`bohdan2019kineticb`).
- The nonthermal fraction at a perpendicular shock is about 0.5 %, reaching at most 10–100 kT
  (`bohdan2023electron` Sect. 3.3): this is pre-acceleration, far from the injection threshold.
- The Weibel field grows as B_sh ≈ 2 √M_A B₀, which raises the injection momentum as √M_A
  (`bohdan2021magnetic` Eqs. 2, 6).

**Electron-driven waves and DSA in 1-D.**
- θ_Bn = 63°, M_A = 63, m_i/m_e = 100: reflected electrons drive waves far upstream and reach DSA with a p⁻⁴
  spectrum; about 7 % of electrons by number and 20 % by energy; protons are not injected (`xu2020electron`
  Sect. 3, Fig. 4).
- θ_Bn = 75°, M_A = 63, m_i/m_e = 256: both species start DSA; reflection depends on the phase of the shock
  reformation (`kumar2021nonthermal` Sects. 3, 4).
- Both are 1-D. A 3-D check of 3 Ω_ci⁻¹ reflects fewer electrons "by a factor of few" (`xu2020electron` Fig. 6).

### 4.2 Oblique shocks: the electron foreshock
- Reflection needs v_e cos θ_Bn ≥ v_sh in the shock frame (`morris2022preacceleration` Sect. 4.1.1).
- The reflected fraction falls with obliquity: about 5 % at 30°, 0.05 % at 63°, none at 74.3° (m_i/m_e = 64,
  M_A = 21, v_sh = 0.264 c; `morris2022preacceleration` Sect. 4.1, Figs. 4, 5). At 63° it is 4 % with the field
  out of the plane.
- Outer foreshock: electron-acoustic waves. They keep about 1 % of the upstream electrons from the shock for
  θ_Bn ≤ 45° (Sect. 6): a small effect.
- Inner foreshock: oblique whistlers that grow into cavities of 1–5 ion skin depths. They confine the reflected
  electrons; the energy gain happens at the shock. Downstream, 0.8 % of the electrons hold 7 % of the energy,
  with an index up to −2.8 (`morris2023preacceleration` Sects. 4.6, 4.8). DSA is not reached.
- A foreshock needs 125–300 Ω_ci⁻¹ to become steady (`morris2023preacceleration` Sect. 4.5;
  `morris2022preacceleration` Sect. 5). The 2-D runs of Sect. 5 last 6–52 Ω_ci⁻¹.
- Pre-existing upstream turbulence: no change at a perpendicular shock (`fulat2023kinetic` Table 1, Fig. 9).
  At 60° it raises the nonthermal fraction from 0.14 % to 0.23 % by number, with the same slope
  (`fulat2026pic` Sect. 3.4). This is one pair of runs.

### 4.3 Quasi-parallel shocks
- **Reference picture** (`park2015simultaneous`): 1-D, θ_Bn = 30°, M_A = 20, m_i/m_e = 100.
  - Electrons are mirror-reflected and gain energy by shock drift.
  - They are then trapped between the shock and the Bell waves that the protons drive, until they reach the
    proton injection momentum.
  - Both species then follow p⁻⁴.
  - K_ep = 3.8 × 10⁻³, and 5.5 × 10⁻³ at m_i/m_e = 400.
- **A different result at the same angle** (`kato2015particle`): electrons are injected in one burst only, with
  an energy index of about 3.0 (Sect. 3.3, Figs. 15, 16).
- **A 1-D survey** at 30° (`gupta2024electron`):
  - the nonthermal electron fraction is at most 0.5 % by number;
  - the spectrum reaches p⁻⁴ only for M_A, M_s ≳ 10 (Sects. 5, 6);
  - the fraction is independent of the mass ratio (Fig. 7d);
  - the electron-to-ion ratio is at most 10⁻² at the real mass ratio, from a formula (Eq. 7).
- **A speed threshold** (`gupta2025speed`): electrons enter DSA once their speed exceeds v_sh. Reflection needs
  an effective electron sonic Mach number below about 3 just upstream (Eq. 4); about 10 % are reflected and at
  most half of those injected (App. C). K_ep ~ 10⁻³ (Eq. 5).
- **Parallel shocks, a competing mechanism** (`shalaby2022mechanism`): ion-cyclotron waves of the
  intermediate-scale instability scatter electrons. It needs M_A at the front below (m_i/m_e)^(1/2)/2 (Eq. 2),
  so a reduced mass ratio switches it off: m_i/m_e = 100 gives two to three orders of magnitude less (Sect. 6).
- **High M_A** (`zekovic2025slams`): proton-driven waves steepen into large-amplitude structures (SLAMS) with
  B_⊥/B₀ from 3–5 at M_A = 20 to at least 30 at M_A = 300 (Table 1). They trap electrons, and a faster
  Fermi-type process gives a steep spectrum, slope 5–5.7, up to the DSA threshold (Eq. 7, Fig. 4). No run
  reaches the transition to DSA (Sect. 5).
- **Regime change in 2-D** (`jikei2026dependence`): at θ_Bn = 20°, shocks with M_A below about 100 are
  dominated by Bell-type fields and put about 1.5 % of the shock energy into nonthermal electrons. Above it
  they are Weibel-dominated, with about 15 % (Sect. 3.3). The runs are at trans-relativistic speeds.
- **A recipe for models** (`arbutina2021nonlinear`): the same ξ for electrons and protons gives
  K_ep = 0.0007–0.0019 at 5000 km/s (Sect. IV).
- **Mildly relativistic, 2-D** (`crumley2019kinetic`): at 10° a p^(−4.2) electron tail with 5 × 10⁻⁴ of the
  energy; at 55°, which is superluminal at that speed, no tail (Sects. 4, 6).

### 4.4 Sequences of mechanisms
Every work that reaches high energy has a sequence: a first stage that lifts electrons out of the thermal pool
and a second that takes over.

| Work | First stage | Second stage | Reached |
|---|---|---|---|
| `amano2007electron` | surfing | mirror reflection and shock drift | γ ~ 2; upstream beam |
| `matsumoto2012electron` | surfing, several encounters | adiabatic gain at the overshoot | γ ~ 7–9 |
| `matsumoto2017electron` | surfing | SSDA with Weibel scattering | γ − 1 ~ 5; upstream escape |
| `bohdan2020kinetic` | (surfing minor) | reconnection, stochastic scattering | 0.5 % nonthermal |
| `riquelme2011electron` | whistler parallel field | escape, return, drift-like gain | Γ ≈ 50 |
| `morris2023preacceleration` | reflection | confinement in whistler cavities, gain at the shock | γ ≈ 13 |
| `xu2020electron` | foot heating, reflection | electron-driven waves, shock-drift cycles | DSA |
| `park2015simultaneous` | reflection, shock drift | trapping by proton-driven waves | DSA |
| `zekovic2025slams` | reflection | Fermi cycles between shock and SLAMS | below the DSA threshold |

### 4.5 Injection without shock microphysics
Two proposals supply seed electrons from atomic or nuclear processes: ionisation of accelerated nuclei
(Morlino 2009) and radioactive decay in the ejecta (Zirakashvili & Aharonian 2011). They were found by the
search and not read (Sect. 10.3). They lie outside this topic.

## 5. Kinetic simulations at a glance

All are PIC unless noted. "in" and "out" give the field orientation to the plane in 2-D. Mach numbers are as the
papers state them; their frames differ (see the notes). Outcomes are for electrons.

**Perpendicular and quasi-perpendicular (θ_Bn ≥ 60°)**

| Work | Dim. | m_i/m_e | M_A | β | θ_Bn | v_sh/c | Outcome |
|---|---|---|---|---|---|---|---|
| `amano2007electron` | 1-D | 100 | ≈15 | 0.08 each | 70–85° | ≈0.075 | reflected fraction ≈5e-2, 1.5e-2, 2e-3 at 70, 75, 80°; none at 85° (Fig. 7, Sect. 3.3) |
| `amano2009electron` | 2-D out | 25 | ≈14 | 0.5 each | 90° | 0.2 (box) | tail of index 2.0–2.5 (Fig. 4) |
| `ohira2007absence` | 2-D electrostatic, periodic foot | 1836 | not stated | not stated | 90° | 0.02 | heating only |
| `matsumoto2012electron` | 2-D out | 25, 100 | 14–30 | β_e 0.5, 4.5 | 90° | 0.2 (box) | tails if the trapping condition and low β_e hold (Table 1) |
| `matsumoto2013electron` | 2-D out | 225 | 44.8 | β_e 0.5 | 90° | not stated | γ_max ~ 12 (Fig. 2b) |
| `matsumoto2017electron` | 3-D | 64 | 20.8 | 0.5 each | 74.3° | 0.26 | power law, index 3.5 (Fig. 4d) |
| `kato2010nonrelativistic` | 2-D in | 30 | ~130 | ~26 | 90° | 0.33 | Maxwellian; T_e/T_i ≈ 0.33 (Fig. 14) |
| `riquelme2011electron` | 2-D in, out; thin 3-D | 25–1600 | 3.5–14 | 0.005–0.5 each | 45–90° | 0.042, 0.14 | 2–10 % by number; index 2.7–4; nothing at 90° (Sect. 4.4) |
| `wieland2016nonrelativistic` | 2-D, field at 45° to the plane | 50 | 28 | β_e ≪ 1 | 90° | 0.13, 0.39 | no tail (Figs. 16, 17) |
| `bohdan2017electron` | 2-D in, 45°, out | 100 | 32–36 | 5e-4, 0.5 | 90° | 0.26–0.29 | nonthermal 0.2–0.7 % (in, 45°); 4–7 % (out) (Table 2) |
| `bohdan2019kinetic`, `bohdan2019kineticb`, `bohdan2020kinetic` | 2-D in (one out) | 50–400 | 23–69 | 5e-4, 0.5 | 90° | 0.26 | nonthermal 0.1–0.7 %; most from reconnection for m_i/m_e ≥ 100 |
| `bohdan2021magnetic` | 2-D in | 50–400 | 23–150 | 5e-4, 0.5 | 90° | 0.26 | field amplification only |
| `fulat2023kinetic` | 2-D in, turbulent upstream | 100 | 32 | 1 | 90° | 0.264 | no change by turbulence |
| `fulat2026pic` | 2-D out, turbulent upstream | 100 | 32 | 1 | 60° | 0.264 | 0.14 % → 0.23 % with turbulence (Sect. 3.4) |
| `morris2023preacceleration` | 2-D out | 50 | 30 | 1 | 60° | 0.263 | 0.8 % of electrons, 7 % of energy (Sect. 4.8) |
| `xu2020electron` | 1-D | 100 | 63 (scan 11–63) | not stated | 63° | 0.15 (wall) | DSA; 7 % by number (Sect. 3) |
| `kumar2021nonthermal` | 1-D | 256 | 63 | ≈1.8 each | 75° | 0.15 | onset of DSA for both species (Sect. 3) |

**Oblique and quasi-parallel (θ_Bn < 60°)**

| Work | Dim. | m_i/m_e | M_A | β | θ_Bn | v_sh/c | Outcome |
|---|---|---|---|---|---|---|---|
| `morris2022preacceleration` | 2-D in (out at 63°, 74.3°) | 64 | 21 | 1 | 30–74.3° | 0.264 | reflected 5 % (30°) → 0.05 % (63°) → 0 (74.3°) |
| `park2015simultaneous` | 1-D | 100, 400 | 20 | not stated | 30° | 0.1 (flow) | DSA; K_ep = 3.8e-3 |
| `kato2015particle` | 1-D | 30 | ≈28 | 1 | 30° | 0.46 | one injection burst; energy index ≈3.0 |
| `fang2019early` | 1-D | 30 | 16.5 | not stated | 15, 30, 45° | 0.2 | tails at all three; DSA "not well developed" (Sect. 4) |
| `crumley2019kinetic` | 2-D in | 64 | 15 | 1 | 10°, 55° | 0.83 | 10°: p^(−4.2); 55°: no tail |
| `shalaby2022mechanism` | 1-D | 1836, 100 | 5.3, 21.3 | ~1e-3 | 0° | 0.133 | efficient only at low M_A and real mass ratio (Sect. 6) |
| `arbutina2021nonlinear` | 2-D thin | 16–100 | 11–16 | not stated | 0° | 0.33–0.4 | p^(−4.2); η ≈ 1e-4 (Sect. IV) |
| `gupta2024electron` | quasi-1-D, 16 runs | 16–1836 | 5–40 | 0.06–64 | 30° | 0.067–0.267 | nonthermal ≲ 0.5 % by number; p⁻⁴ for M_A, M_s ≳ 10 |
| `gupta2025speed` | quasi-1-D | 100 | 20, 5 | not stated | 30° | 0.133 | threshold: speed above v_sh |
| `zekovic2025slams` | 1-D, 2-D in | 32 | 20–300 | not stated | 0° | 0.033–0.33 | slope 5–5.7 below the DSA threshold |
| `jikei2026dependence` | 2-D in | 100 | regime change near 100 | not stated | 20° | trans-relativistic | 1.5 % (Bell) against 15 % (Weibel) of the shock energy |

**Test electrons in hybrid or prescribed fields**

| Work | Method | M_A | θ_Bn | Outcome |
|---|---|---|---|---|
| `burgess2006simulations` | 2-D hybrid, test electrons | 5.1, 2.1 | 80–88° | ripples raise the maximum energy almost tenfold (Fig. 5, Sect. 4) |
| `trotta2019electron` | 2-D and 3-D hybrid, test electrons | 2.9–6.6 | 80–87° | threshold near M = 3.5 (Sect. 5) |
| `guo2010effect` | 2-D hybrid with imposed turbulence, test electrons | 4.0 | mean 60–90° | up to 200–300 times the initial 100 eV at 90° (Fig. 4) |
| `guo2015acceleration` | 3-D test electrons, prescribed shock and turbulence | assumed ≫ 1 | mean 0–90° | fraction above 10 keV from 0–9 % (weak turbulence) to 7–11 % (strong) across 0–90° (Table 1) |

**Low-Mach reference** (not SNR-like): `guo2014nonthermal`, 2-D, M_s = 3, β_p = 20, θ_Bn = 63°. Shock drift
plus waves driven by the returning electrons gives repeated cycles and a tail with about 15 % of the electrons
(Sect. 5, Fig. 4c).

**The definitions of "nonthermal" differ** between the works, so the fractions above cannot be compared
directly:
- excess over a Maxwellian fit (`bohdan2017electron`);
- above five times the thermal momentum (`gupta2024electron` Eq. 3);
- above five times the energy of the thermal peak (`xu2020electron`);
- above 4 m_e c (`kumar2021nonthermal`);
- above three times the peak momentum (`jikei2026dependence`).

## 6. Dependence on θ_Bn, non-planar fronts and global geometry

### 6.1 Which angles have been simulated
- 90° only: most of the multi-dimensional work (Sect. 5).
- Scans within one study:
  - 45–90° at M_A = 7 and m_i/m_e = 400, 2-D, in the whistler regime (`riquelme2011electron` Sect. 4.1);
  - 30–74.3° at M_A = 21, 2-D, reflection and foreshock only, short runs (`morris2022preacceleration`);
  - 70–85° in 1-D (`amano2007electron`);
  - 15–45° in 1-D (`fang2019early`).
- Single angles elsewhere: 0°, 10°, 20°, 30°, 55°, 60°, 63°, 74.3°, 75°.
- The one 3-D run is at 74.3° and is "too short" for a foreshock (`bohdan2023electron` Sect. 4).
- No study follows electrons from the thermal pool to DSA through the range of θ_Bn at fixed parameters.
  No controlled scan crosses the oblique range between 45° and 63° to DSA.

### 6.2 What is stated about the dependence
- Shock drift gains more and reflects fewer electrons as θ_Bn grows (`amano2007electron` Sect. 1). The
  controlling quantity is M_A / cos θ_Bn (`katou2019theory` Eq. 21; `amano2022theory` Eq. 14;
  `wang2026generation` Eq. 22).
- The dependence is not monotonic at M_A = 10: 75° and 80° are marginal, 85° is unstable to the whistlers that
  SSDA needs, 87° is weaker again (`wang2026generation` Sect. IV C).
- In the whistler regime the hardest spectrum is at 70–75° for M_A = 7 and at 60° for M_A = 14; there is no
  tail at 90° (`riquelme2011electron` Sects. 4.1, 4.4). Electrons escape upstream if
  100 cos²θ_Bn ≳ M_A / (m_i/m_e)^(1/2) (Sect. 5).
- Conjectures, not tested:
  - quasi-parallel shocks accelerate both species and quasi-perpendicular ones electrons only
    (`park2015simultaneous`, Conclusions);
  - the optimum for electrons lies between 60° and 70° (`xu2020electron` Sect. 3);
  - a growing precursor makes the transport independent of the initial angle (`kumar2021nonthermal` Sect. 4).
- Cut-offs in a model for SNRs: near 85° electrons no longer escape, and near 88° the maximum energy is limited
  (`amano2007electron` Sect. 3.5).

### 6.3 Protons decide where proton-driven waves exist
- Hybrid simulations give ion acceleration efficiencies of 5–15 % below about 40°, a drop above 45°, and none
  above 60° (`caprioli2014simulations` Fig. 3, Sect. 7). Upstream δB/B₀ reaches 3–6 at quasi-parallel shocks
  and stays below 0.1 at 80° in 3-D (Figs. 4, 14).
- So on a sphere in a uniform field, proton-driven turbulence exists on the quasi-parallel caps only
  (reading of this survey).
- For a uniform field the caps with θ_Bn < 45° are 29 % of the sphere, the band from 45° to 75° is 45 %, and
  the band above 75° is 26 % (our arithmetic).
- This is disputed: other methods find ion acceleration at quasi-perpendicular shocks (van Marle et al. 2018;
  Reville & Bell 2013; not read; via `pohl2020pic` Sect. 5.3). Orusa & Caprioli 2023 find it in 3-D hybrid runs
  (not read; cited by `mcginness2026sn`).
- Electron injection without injected protons was seen at 63° (`xu2020electron` Sect. 3).

### 6.4 Local obliquity: ripples, turbulence and self-generated fields
The local θ_Bn along one front differs from the mean.
- **Ripples.** Ion-scale ripples of a quasi-perpendicular front scatter electrons within the layer. The local
  angle spreads by 1.4° at M = 2.9 and by 15.2° at M = 6.6 (`trotta2019electron` Fig. 2). An electron can enter
  where the front is locally perpendicular and leave where it is not (Sect. 5). 2-D ripples over-trap; 3-D
  gives a shallower fall-off (Figs. 4, 9, 10). These are hybrid simulations with test electrons.
- **Pre-existing turbulence.** Meandering field lines carry electrons back to the front several times
  (`guo2010effect`). With strong turbulence, a shock with a parallel mean field accelerates electrons where the
  local field is quasi-perpendicular (`guo2015acceleration` Table 1).
- **Self-generated fields at quasi-parallel shocks.**
  - Proton-driven waves raise the local angle to about 80° (`kato2015particle` Sects. 3.1, 4).
  - Patches beyond arccos(v_sh/c) are superluminal (`gupta2024electron` Sect. 4).
  - For magnetised electrons the parallel shock then resembles a quasi-perpendicular one
    (`zekovic2025slams` Sects. 4.4, 5; `jikei2026dependence` Sect. 3.3).
  - The same idea is offered as the link to the radial fields of young remnants (`amano2020observational` p. 5).

### 6.5 Fronts curved on the large scale
- No kinetic simulation of electron injection at a spherical or otherwise curved shock with SNR-like parameters
  was found.
- The curved-shock literature concerns the Earth's bow shock and was not read (Sect. 10.2):
  - Krauss-Varban & Burgess 1991: a scatter-free curved shock; electrons that enter near 90° leave at a slightly
    lower θ_Bn (via `burgess2006simulations`);
  - Vandas 1995, 2001: a curved shock of finite width; adiabatic theory fails on the observed spectra;
  - Krauss-Varban 1994: scattering at a curved shock;
  - Savoini & Lembège 2001: a full-particle simulation of a curved shock and its electron foreshock.
- Curvature limits the energy at the bow shock to about 10 keV for shock drift with or without scattering; this
  is called "not an intrinsic limit" (`katou2019theory` Sect. 4). `amano2022theory` gives the limit as a formula
  (Eq. 22, Sect. 4.5).
- Accelerated electrons travel along the field line and escape "due to the change in local θ_Bn"
  (`katou2019theory` Sect. 4). The most energetic upstream electrons at the bow shock come from the tangent
  point of the field line (Anderson et al. 1979, not read; via `amano2022nonthermal` Sect. 3.1).
- Whole remnants are modelled with MHD or semi-analytic methods and an assumed ε(θ_Bn) (Sect. 3.3). The
  brightness then peaks about 35° from the field axis, because field compression rises with θ faster than the
  assumed efficiency falls (`winner2020evolution` Sect. 3.2, Fig. 5).

### 6.6 Non-stationarity
- Shock reformation has a period of 1.5–2 Ω_ci⁻¹ (`bohdan2017electron` Sect. 3.2; `amano2007electron`
  Sect. 3.4; `fulat2026pic` Sect. 3.1).
- Buneman intensity, surfing, reconnection and the reflected fraction follow that cycle
  (`bohdan2017electron` Fig. 15; `kumar2021nonthermal` Sect. 4; `gupta2025speed` App. C).
- Injection can be intermittent: one burst in `kato2015particle` (Figs. 15, 16).

## 7. Numerical limits that the literature acknowledges

- **Mass ratio.** The thresholds scale differently: surfing as (m_i/m_e)^(2/3), whistlers and SSDA as
  (m_i/m_e)^(1/2), the Buneman instability as (β_e m_i/m_e)^(1/2). A run can therefore sit in a regime that the
  real shock is not in (reading of this survey).
  - Electron energy fractions are more than ten times too large at m_i/m_e = 16 (`crumley2019kinetic` App. A).
  - The parallel-shock mechanism of `shalaby2022mechanism` vanishes at 100.
  - `gupta2024electron` finds the fraction independent of the mass ratio.
- **Shock speed.** Most PIC runs use more than 0.1 c, ten times the real speed (`pohl2020pic` Sects. 2.2, 5).
  A few 1-D runs go down to 0.03–0.07 c (`zekovic2025slams`; `gupta2024electron`). Consequences:
  - a larger superluminal range of θ_Bn (Sect. 2);
  - a single surfing reflection is already relativistic, which fails at 0.01 c (`amano2022nonthermal` Sect. 4.5);
  - an injection gap of 10³ in place of 10⁵–10⁶ (`bohdan2023electron` Sect. 4.3).
  - The scaling with the shock speed is "not entirely clear" (`bohdan2019kineticb` Sect. 4).
- **Dimensionality.**
  - 1-D holds planar Buneman waves and no Weibel modes.
  - 2-D with the field out of the plane holds Buneman modes; whistlers and ripples are excluded.
  - 2-D with the field in the plane holds whistlers, Weibel modes and reconnection; Buneman modes are suppressed.
  - Only 3-D holds all of them (`riquelme2011electron` Appendix; `matsumoto2017electron`;
    `bohdan2017electron` Sect. 4).
  - "Only 3D simulations can provide a complete picture" (`pohl2020pic` Sect. 5.3).
- **Run time and box.** A few to a few tens of Ω_ci⁻¹ in 2-D; foreshocks need 125–300 (Sect. 4.2). Several
  boxes are too narrow to hold ripples (`bohdan2017electron` Sect. 3.2; `arbutina2021nonlinear` Sect. II).
- **Numerical drag** on fast particles in PIC (Kato 2013, arXiv:1312.5507, not read; cited by `kato2015particle`,
  `xu2020electron`, `kumar2021nonthermal`).

## 8. Disagreements, open questions and gaps

### 8.1 Disagreements between the works
1. **Does surfing matter at real parameters?** Yes above a Mach threshold (`matsumoto2012electron`;
   `matsumoto2017electron`). No (`ohira2007absence`; `riquelme2011electron`; `bohdan2019kineticb`).
2. **Which Mach number do SNR shocks need?** M_A ≲ 20 for whistlers (`riquelme2011electron`) against
   M_A ≳ (m_i/m_e)^(2/3) for surfing (`matsumoto2012electron`).
3. **Shock drift at high Mach number.** It "cannot operate efficiently" (`guo2014nonthermal` Sect. 5) against
   favoured at high M_A / cos θ_Bn (`katou2019theory`; `amano2020observational`).
4. **Parallel shocks.** Proton-driven turbulence and reflection (`park2015simultaneous`; `gupta2024electron`)
   against ion-cyclotron waves that need the real mass ratio and low M_A (`shalaby2022mechanism`) against
   steepened structures at high M_A (`zekovic2025slams`).
5. **Reduced mass ratio.** Trustworthy (`gupta2024electron`; `arbutina2021nonlinear`) or not
   (`shalaby2022mechanism`).
6. **Ions above 60°.** Only a steep tail at 63° (`xu2020electron`) against an ion tail at 75°
   (`kumar2021nonthermal`).
7. **Quasi-parallel or quasi-perpendicular for electrons.** PIC efficiencies are higher at 63° than at 30°
   (`xu2020electron` Sect. 3; `gupta2024electron` Sect. 7). The SN 1006 models need the opposite
   (`winner2020evolution`). The radio population favours quasi-perpendicular (`west2017connection`).
8. **The injection threshold** itself (Sect. 2).
9. **What limits the energy at the bow shock:** curvature (`katou2019theory`) or scattering
   (`amano2020observational`).

### 8.2 Open questions that the works state
- The functional form of the electron efficiency against obliquity (`winner2020evolution` Sect. 2.2).
- The relative roles of surfing, SSDA and reconnection; the scattering rates; whether SSDA survives Weibel
  turbulence (`amano2022nonthermal` Sects. 5.4, 6).
- The scaling with shock speed, mass ratio and Mach number; multi-dimensional parallel shocks
  (`bohdan2023electron` Sect. 6).
- How the injection efficiency depends on shock conditions; the passage from whistler-accelerated electrons to
  p⁻⁴ (`marcowith2016microphysics` Sect. 3.4.1).
- More oblique shocks, where transport is more complicated (`gupta2025speed` Sect. IV).
- Multi-dimensional confirmation of electron DSA at quasi-perpendicular shocks (`xu2020electron` Sect. 4).
- Fully kinetic oblique runs; 3-D (`pohl2020pic` Sect. 5.3).
- Which part of upstream turbulence acts, magnetic scattering or corrugation of the front
  (`fulat2026pic` Sect. 3.3).

### 8.3 Gaps relative to this topic (reading of this survey)
- **No curved front.** Electron injection on a front that is curved on the large scale, with θ_Bn varying along
  it, has no kinetic simulation at SNR-like parameters among the works found. The bow-shock works of
  Sect. 6.5 are the direct predecessors and are not read. Until they are, a novelty claim about curvature is
  not certain.
- **No single-parameter scan.** The θ_Bn dependence of injection is known piecewise, from studies with
  different mass ratios, Mach numbers, dimensions and definitions.
- **The global models need exactly that function.** ε(θ_Bn) for electrons is assumed in every model of
  Sect. 3.3, and observations of SN 1006 are used to choose between assumed forms.
- **Transport along the front.** Statements that electrons leave their region through the change of the local
  θ_Bn exist for the bow shock only (Sect. 6.5).
- **Regime.** SNR shocks have M_A of tens to hundreds and are Weibel-mediated in the simulations above. A shock
  with M_A of a few is in the regime of `trotta2019electron`, `burgess2006simulations` and
  `guo2014nonthermal`. Which regime a spherical run reaches decides which of these works it can be compared
  with.
- **Superluminal band.** At simulation speeds a large part of a sphere is superluminal (Sect. 2); at SNR speeds
  it is under 2 %.

These gaps do not decide the research question of the topic; that decision is the owner's.

## 9. Search log and limits of this survey

**Searches** (2026-11-12, arXiv API, abstracts and titles; 60 newest results each; 238 results not already
known). Query terms:
- "electron injection" and shock;
- "electron acceleration", "supernova remnant", "particle-in-cell";
- "quasi-parallel", electron, shock, "particle-in-cell";
- "shock drift acceleration" and electron;
- "shock surfing" and electron;
- "curved shock", "shock curvature" or "spherical shock", with electron and acceleration;
- obliquity, electron, shock, acceleration;
- "electron-to-proton", shock, acceleration;
- "SN 1006" with "quasi-parallel", obliquity or "polar cap";
- "electron foreshock";
- "pre-acceleration", electron, shock;
- Weibel, electron, "supernova remnant", shock.

A second round searched by author and title (the groups of the works above; X-ray polarimetry of remnants;
bilateral remnants; laboratory shocks; blast shells in PIC). Works were also taken from the reference lists
of the works read.

**Selection.** 54 works were read for this survey (Sect. 10.1), chosen for direct relevance to electron
injection at SNR-like shocks and to obliquity. Five further notes existed already.

**Reading depth.** Each note states it. Reviews were read in the named sections only. For papers not marked
as direct prior work, the introduction, set-up and conclusions were read and the results skimmed for numbers.

**Checks.** The headline numbers of 20 works were checked against the paper texts after the notes were
written. Numbers were not read off figures, except where a note says so (`amano2007electron` Fig. 7).

**Limits.**
- Only arXiv could be reached. Works that are not there were not read: the theory of the 1980s and 1990s,
  the bow-shock literature in J. Geophys. Res., Matsumoto et al. 2015 (Science), and several observational
  papers. They are in Sect. 10.2.
- Journal volumes and pages come from the arXiv metadata, or from reference lists of other works read here.
  Where neither gave them, the entry says "not checked".
- The curved-shock query found no kinetic work on SNR-like shocks. One search is not proof that none exists.
- Not covered: low-Mach shocks in galaxy clusters beyond `guo2014nonthermal`; relativistic shocks; laboratory
  experiments; studies of electron heating alone; in-situ studies of interplanetary shocks.
- Forward citations (who cites the key works) were searched by keyword only, not by a citation index.

## 10. References

### 10.1 Works read (notes in `research/literature/`)

Observations and global models:
- `vink2012supernova` — J. Vink 2012, Astron. Astrophys. Rev. 20, 49 (article number as printed in another
  work's reference list). arXiv:1112.0576. DOI 10.1007/s00159-011-0049-1.
- `reynolds2004microphysics` — S. P. Reynolds 2004, Adv. Space Res. (COSPAR proceedings; volume not checked).
  arXiv:astro-ph/0308217.
- `morlino2012strong` — G. Morlino & D. Caprioli 2012, A&A 538, A81. arXiv:1105.6342.
  DOI 10.1051/0004-6361/201117855.
- `bocchino2011constraints` — F. Bocchino, S. Orlando, M. Miceli & O. Petruk 2011, A&A 531, A129.
  arXiv:1105.2689. DOI 10.1051/0004-6361/201016341.
- `reynoso2013radio` — E. M. Reynoso, J. P. Hughes & D. A. Moffett 2013, AJ 145, 104. arXiv:1302.4678.
  DOI 10.1088/0004-6256/145/4/104.
- `west2017connection` — J. L. West, S. Safi-Harb & G. Ferrand 2017, A&A 597, A121. arXiv:1609.02887.
  DOI 10.1051/0004-6361/201628079.
- `zhou2023magnetic` — P. Zhou et al. 2023, ApJ (accepted; volume not checked). arXiv:2309.01879.
- `giuffrida2022supernova` — R. Giuffrida et al. 2022, Nat. Commun. 13, 5098. arXiv:2208.14491.
- `mcginness2026sn` — E. McGinness, R. Diesing, D. Caprioli & F. Acero 2026, ApJ (volume not checked).
  arXiv:2606.12111. DOI 10.3847/1538-4357/aea3f3.
- `winner2020evolution` — G. Winner, C. Pfrommer, P. Girichidis, M. Werhahn & M. Pais 2020, MNRAS 499, 2785.
  arXiv:2006.06683. DOI 10.1093/mnras/staa2989.

Quasi-perpendicular shocks, PIC:
- `amano2007electron` — T. Amano & M. Hoshino 2007, ApJ 661, 190. arXiv:astro-ph/0612204. DOI 10.1086/513599.
- `amano2009electron` — T. Amano & M. Hoshino 2009, ApJ 690, 244. arXiv:0805.1098.
  DOI 10.1088/0004-637X/690/1/244.
- `ohira2007absence` — Y. Ohira & F. Takahara 2007, ApJ 661, L171. arXiv:0705.2061.
- `matsumoto2012electron` — Y. Matsumoto, T. Amano & M. Hoshino 2012, ApJ 755, 109. arXiv:1204.6312.
  DOI 10.1088/0004-637X/755/2/109.
- `matsumoto2013electron` — Y. Matsumoto, T. Amano & M. Hoshino 2013, Phys. Rev. Lett. 111, 215003.
  arXiv:1310.8367. DOI 10.1103/PhysRevLett.111.215003.
- `matsumoto2017electron` — Y. Matsumoto, T. Amano, T. N. Kato & M. Hoshino 2017, Phys. Rev. Lett. 119, 105101.
  arXiv:1709.03673. DOI 10.1103/PhysRevLett.119.105101.
- `kato2010nonrelativistic` — T. N. Kato & H. Takabe 2010, ApJ 721, 828. arXiv:1008.0265.
  DOI 10.1088/0004-637X/721/1/828.
- `riquelme2011electron` — M. A. Riquelme & A. Spitkovsky 2011, ApJ 733, 63. arXiv:1009.3319.
  DOI 10.1088/0004-637X/733/1/63.
- `wieland2016nonrelativistic` — V. Wieland, M. Pohl, J. Niemiec, I. Rafighi & K.-I. Nishikawa 2016, ApJ 820, 62.
  arXiv:1602.05064. DOI 10.3847/0004-637X/820/1/62.
- `bohdan2017electron` — A. Bohdan, J. Niemiec, O. Kobzar & M. Pohl 2017, ApJ 847, 71. arXiv:1708.05528.
  DOI 10.3847/1538-4357/aa872a.
- `bohdan2019kinetic` — A. Bohdan et al. 2019, ApJ 878, 5. arXiv:1904.13153. DOI 10.3847/1538-4357/ab1b6d.
- `bohdan2019kineticb` — A. Bohdan et al. 2019, ApJ 885, 10. arXiv:1909.05294. DOI 10.3847/1538-4357/ab43cf.
- `bohdan2020kinetic` — A. Bohdan et al. 2020, ApJ 893, 6. arXiv:2003.01879. DOI 10.3847/1538-4357/ab7cd6.
- `bohdan2021magnetic` — A. Bohdan et al. 2021, Phys. Rev. Lett. 126, 095101. arXiv:2102.04328.
  DOI 10.1103/PhysRevLett.126.095101.

Oblique shocks, foreshock, turbulent upstream:
- `morris2022preacceleration` — P. J. Morris, A. Bohdan, M. S. Weidl & M. Pohl 2022, ApJ 931, 129.
  arXiv:2204.11569. DOI 10.3847/1538-4357/ac69c7.
- `morris2023preacceleration` — P. J. Morris et al. 2023, ApJ 944, 13. arXiv:2301.00872.
  DOI 10.3847/1538-4357/acaec8.
- `fulat2023kinetic` — K. Fulat, A. Bohdan, G. Torralba Paz & M. Pohl 2023, ApJ 959, 119. arXiv:2310.14388.
- `fulat2026pic` — K. Fulat et al. 2026, accepted by ApJ (volume not checked). arXiv:2605.21986.

Quasi-parallel shocks and injection thresholds:
- `park2015simultaneous` — J. Park, D. Caprioli & A. Spitkovsky 2015, Phys. Rev. Lett. 114, 085003.
  arXiv:1412.0672. DOI 10.1103/PhysRevLett.114.085003.
- `kato2015particle` — T. N. Kato 2015, ApJ 802, 115. arXiv:1407.1971.
- `crumley2019kinetic` — P. Crumley, D. Caprioli, S. Markoff & A. Spitkovsky 2019, MNRAS 485, 5105.
  arXiv:1809.10809. DOI 10.1093/mnras/stz232.
- `xu2020electron` — R. Xu, A. Spitkovsky & D. Caprioli 2020, ApJ Lett. 897, L41. arXiv:1908.07890.
  DOI 10.3847/2041-8213/aba11e.
- `kumar2021nonthermal` — N. Kumar & B. Reville 2021, ApJ Lett. 921, L14. arXiv:2110.09939.
  DOI 10.3847/2041-8213/ac30e0.
- `fang2019early` — J. Fang, C.-Y. Lu, J.-W. Yan & H. Yu 2019, Res. Astron. Astrophys. 19, 182. arXiv:1908.08170.
- `shalaby2022mechanism` — M. Shalaby, R. Lemmerz, T. Thomas & C. Pfrommer 2022, ApJ 932, 86. arXiv:2202.05288.
  DOI 10.3847/1538-4357/ac6ce7.
- `arbutina2021nonlinear` — B. Arbutina & V. Zeković 2021, Astropart. Phys. 127, 102546. arXiv:2012.15117.
  DOI 10.1016/j.astropartphys.2020.102546.
- `gupta2024electron` — S. Gupta, D. Caprioli & A. Spitkovsky 2024, ApJ 976, 10. arXiv:2408.16071.
  DOI 10.3847/1538-4357/ad7c4c.
- `gupta2025speed` — S. Gupta, D. Caprioli & A. Spitkovsky 2025, preprint arXiv:2506.09134. Probably published
  as ApJ Lett. 994, L34 (2025): `mcginness2026sn` cites that entry for the K_ep formula; the identity was not
  verified.
- `zekovic2025slams` — V. Zeković, A. Spitkovsky & Z. Hemler 2025, ApJ 988, 40. arXiv:2408.02084.
  DOI 10.3847/1538-4357/ade05f.
- `jikei2026dependence` — T. Jikei, D. Grošelj & L. Sironi 2026, preprint arXiv:2607.05778.

Shock drift theory, heliospheric benchmark, geometry:
- `katou2019theory` — T. Katou & T. Amano 2019, ApJ 874, 119. arXiv:1903.02277. DOI 10.3847/1538-4357/ab0d8a.
- `amano2020observational` — T. Amano et al. 2020, Phys. Rev. Lett. 124, 065101. arXiv:2002.06787.
  DOI 10.1103/PhysRevLett.124.065101.
- `wang2026generation` — R. Wang & T. Amano 2026, accepted by Phys. Plasmas (volume not checked).
  arXiv:2603.19018.
- `raptis2025revealing` — S. Raptis et al. 2025, Nat. Commun. 16, 488. arXiv:2502.10643.
  DOI 10.1038/s41467-024-55641-9.
- `guo2014nonthermal` — X. Guo, L. Sironi & R. Narayan 2014, ApJ 794, 153. arXiv:1406.5190.
  DOI 10.1088/0004-637X/794/2/153.
- `burgess2006simulations` — D. Burgess 2006, ApJ 653, 316. arXiv:astro-ph/0610714. DOI 10.1086/508805.
- `guo2010effect` — F. Guo & J. Giacalone 2010, ApJ 715, 406. arXiv:1003.5946.
  DOI 10.1088/0004-637X/715/1/406.
- `guo2015acceleration` — F. Guo & J. Giacalone 2015, ApJ 802, 97. arXiv:1409.5854.
  DOI 10.1088/0004-637X/802/2/97.
- `trotta2019electron` — D. Trotta & D. Burgess 2019, MNRAS 482, 1154. arXiv:1808.00812.
  DOI 10.1093/mnras/sty2756.
- `caprioli2014simulations` — D. Caprioli & A. Spitkovsky 2014, ApJ 783, 91. arXiv:1310.2943.
  DOI 10.1088/0004-637X/783/2/91.

Reviews:
- `bohdan2023electron` — A. Bohdan 2023, Plasma Phys. Control. Fusion 65 (article number not checked).
  arXiv:2211.13992. DOI 10.1088/1361-6587/aca5b2.
- `amano2022nonthermal` — T. Amano et al. 2022, Rev. Mod. Plasma Phys. 6, 29. arXiv:2209.03521.
- `marcowith2016microphysics` — A. Marcowith et al. 2016, Rep. Prog. Phys. 79, 046901. arXiv:1604.00318.
  DOI 10.1088/0034-4885/79/4/046901.
- `pohl2020pic` — M. Pohl, M. Hoshino & J. Niemiec 2020, Prog. Part. Nucl. Phys., article 103751 (volume not
  checked). arXiv:1912.02673. DOI 10.1016/j.ppnp.2019.103751.

Notes that existed before this survey:
- `hoshino2002nonthermal` — M. Hoshino & N. Shimada 2002, ApJ 572, 880. arXiv:astro-ph/0203073.
- `amano2010critical` — T. Amano & M. Hoshino 2010, Phys. Rev. Lett. 104, 181102. arXiv:1004.3838.
- `amano2022theory` — T. Amano & M. Hoshino 2022, ApJ 927, 132. arXiv:2201.11416.

### 10.2 Works cited above and not read
Entries are as printed in the reference lists of the works named in brackets.

Observations and models of remnants:
- Fulbright, M. S. & Reynolds, S. P. 1990, ApJ, 357, 591 [`bocchino2011constraints`, `west2017connection`].
- Leckband, J. A., Spangler, S. R., & Cairns, I. H. 1989, ApJ, 338, 963 [`west2017connection`].
- Rothenflug, R., Ballet, J., Dubner, G., et al. 2004, A&A, 425, 121 [`bocchino2011constraints`].
- Koyama, K., et al. 1995, Nature, 378, 255 [`vink2012supernova`].
- Jokipii, J. R. 1987, ApJ, 313, 842 [`reynolds2004microphysics`].
- Cassam-Chenaï, G., Hughes, J. P., Reynoso, E. M., Badenes, C., & Moffett, D. 2008, ApJ, 680, 1180
  [`reynoso2013radio`].
- Gaensler, B. M. 1998, ApJ, 493, 781 [`bocchino2011constraints`]. On arXiv: astro-ph/9707149.
- Völk, H. J., Berezhko, E. G., & Ksenofontov, L. T. 2003, A&A, 409, 563 [`bocchino2011constraints`,
  `winner2020evolution`]. On arXiv: astro-ph/0306016.
- Ellison, D. C., Baring, M. G., & Jones, F. C. 1995, ApJ, 453, 873 [`bocchino2011constraints`].

Surfing, Buneman instability, Weibel turbulence:
- Buneman, O. 1958, Phys. Rev. Lett., 1, 8 [`bohdan2017electron`].
- Papadopoulos, K. 1988, Ap&SS, 144, 535 [`amano2007electron`].
- Shimada, N., & Hoshino, M. 2000, ApJ, 543, L67 [`amano2007electron`, `wieland2016nonrelativistic`].
- McClements, K. G., Dieckmann, M. E., Ynnerman, A., Chapman, S. C., & Dendy, R. O. 2001, Phys. Rev. Lett.,
  87, 255002 [`amano2007electron`].
- Matsumoto, Y., Amano, T., Kato, T. N., & Hoshino, M. 2015, Science, 347, 974 [`bohdan2020kinetic`].
- Bohdan, A., Weidl, M. S., Morris, P. J., & Pohl, M. 2022, Physics of Plasmas, 29, 052301
  [`morris2023preacceleration`, `fulat2026pic`].

Shock drift, whistlers, heliospheric observations:
- Wu, C. S. 1984, J. Geophys. Res., 89, 8857 [`katou2019theory`, `amano2007electron`].
- Leroy, M. M., & Mangeney, A. 1984, Annales Geophysicae, 2, 449 [`katou2019theory`].
- Krauss-Varban, D., & Wu, C. S. 1989, J. Geophys. Res., 94, 15367 [`katou2019theory`, `guo2014nonthermal`].
- Ball, L., & Melrose, D. B. 2001, PASA, 18, 361 [`park2015simultaneous`, `xu2020electron`].
- Levinson, A. 1992, ApJ, 401, 73 [`katou2019theory`, `riquelme2011electron`].
- Krasnoselskikh, V. V., Lembège, B., Savoini, P., & Lobzin, V. V. 2002, Physics of Plasma, 9, 1192
  [`riquelme2011electron`].
- Oka, M., Terasawa, T., Seki, Y., et al. 2006, Geophys. Res. Lett., 33, L24104 [`katou2019theory`].
- Oka, M., Wilson, III, L. B., Phan, T. D., et al. 2017, ApJL, 842, L11 [`katou2019theory`].
- Masters, A., Stawarz, L., Fujimoto, M., et al. 2013, Nature Physics, 9, 164 [`katou2019theory`,
  `guo2015acceleration`].
- Anderson, K. A., Lin, R. P., Martel, F., Lin, C. S., Parks, G. K., & Rème, H. 1979, Geophys. Res. Lett., 6,
  401 [`burgess2006simulations`, `amano2022nonthermal`].

Curved and rippled shocks:
- Krauss-Varban, D., & Burgess, D. 1991, J. Geophys. Res., 96, 143 [`burgess2006simulations`,
  `katou2019theory`].
- Krauss-Varban, D. 1994, J. Geophys. Res., 99, 2537 [`burgess2006simulations`].
- Vandas, M. 1995, J. Geophys. Res., 100, 23499; Vandas, M. 2001, J. Geophys. Res., 106, 1859
  [`burgess2006simulations`, `katou2019theory`].
- Savoini, P., & Lembège, B. 2001, J. Geophys. Res., 106, 12975 [`burgess2006simulations`].
- Lembège, B., & Savoini, P. 2002, J. Geophys. Res., 107, 1037 [`burgess2006simulations`,
  `trotta2019electron`].
- Lowe, R. E., & Burgess, D. 2003, Ann. Geophys., 21, 671 [`burgess2006simulations`, `trotta2019electron`].
- Jokipii, J. R., & Giacalone, J. 2007, ApJ, 660, 336 [`guo2010effect`].

Protons and quasi-parallel shocks:
- Bell, A. R. 2004, MNRAS, 353, 550 [`park2015simultaneous`].
- Caprioli, D., Pop, A.-R., & Spitkovsky, A. 2015, ApJ, 798, L28 [`gupta2024electron`]. On arXiv: 1409.8291.
- Blasi, P., Gabici, S., & Vannoni, G. 2005, MNRAS, 361, 907 [`arbutina2021nonlinear`, `morlino2012strong`].
- Orusa, L., & Caprioli, D. 2023, PhRvL, 131, 095201 [`mcginness2026sn`]. On arXiv: 2305.10511.
- van Marle et al. 2018 and Reville & Bell 2013: named in `pohl2020pic` Sect. 5.3; entries not copied.

### 10.3 Found by the search and not read (for follow-up)
Titles and identifiers are from the arXiv metadata.
- 2606.19979, R. Ferrazzoli et al. 2026, "Supernova Remnants in the IXPE era: a review".
- 2505.09712, P. Zhou et al. 2025, "X-ray polarization in SN 1006 southwest shows spatial variations and
  differences with the radio band".
- 1702.02054, S. Katsuda 2017, "Supernova of 1006 (G327.6+14.6)".
- 1502.03133, E. M. Schneiter et al. 2015, "3D MHD simulation of polarized emission in SN 1006".
- 0811.2319, O. Petruk et al. 2009, "Aspect angle for interstellar magnetic field in SN 1006".
- 0704.0890, S. Orlando et al. 2007, "On the Origin of Asymmetries in Bilateral Supernova Remnants".
- 1510.08536, J. L. West et al. 2016, "The Connection between Supernova Remnants and the Galactic Magnetic
  Field: A Global Radio Study of the Axisymmetric Sample".
- 1805.00128, M. Pais et al. 2018, "The effect of cosmic-ray acceleration on supernova blast wave dynamics".
- 1905.07414, R. Diesing & D. Caprioli 2019, "On the Spectrum of Electrons Accelerated in Supernova Remnants".
- 2103.02375, P. Cristofari et al. 2021, "Cosmic ray protons and electrons from supernova remnants".
- 1209.1466, K. Maeda 2012, "Injection and Acceleration of Electrons at A Strong Shock: Radio and X-ray Study
  of Young Supernova 2011dh".
- 0908.3217, G. Morlino 2009, "Cosmic-ray electron injection from the ionization of nuclei".
- 1011.4775, V. N. Zirakashvili & F. A. Aharonian 2010, "Radioactivity and Electron Acceleration in Supernova
  Remnants".
- astro-ph/0002347, M. E. Dieckmann et al. 2000, "Electron acceleration due to high frequency instabilities at
  supernova remnant shocks".
- 0808.3195, Y. Ohira & F. Takahara 2008, "Oblique Ion Two-Stream Instability in the Foot Region of a
  Collisionless Shock".
- 1004.2771, N. Shimada, M. Hoshino & T. Amano 2010, "Structure of a strong supernova shock wave and rapid
  electron acceleration confined in its transition region".
- 2008.05920, A. Bohdan et al. 2020, "Kinetic simulation of nonrelativistic perpendicular shocks of young
  supernova remnants. IV. Electron heating".
- 2308.15835, G. Torralba Paz et al. 2023, "Prediction and Anomaly Detection of accelerated particles in PIC
  simulations using neural networks".
- 1209.3938, J. Niemiec et al. 2012, "Nonrelativistic parallel shocks in unmagnetized and weakly magnetized
  plasmas".
- 1409.7393, X. Guo, L. Sironi & R. Narayan 2014, "Non-Thermal Electron Acceleration in Low Mach Number
  Collisionless Shocks. II.".
- 2202.02689, J. Fang et al. 2022, "Kinetic Simulation on Electron, Proton and Helium Acceleration in a
  Nonrelativistic Quasiparallel shock".
- 2312.13365, S. Gupta, D. Caprioli & A. Spitkovsky 2023, "Return Currents in Collisionless Shocks".
- 2312.07933, T. Jikei et al. 2023, "Enhanced Magnetic Field Amplification by Ion-Beam Weibel Instability in
  Weakly Magnetized Astrophysical Shocks".
- 2512.03169, T. Jikei et al. 2025, "Magnetic Field Amplification and Particle Acceleration in Weakly Magnetized
  Trans-relativistic Electron-ion Shocks".
- 2507.13436, L. Orusa et al. 2025, "The role of three-dimensional effects on ion injection and acceleration in
  perpendicular shocks".
- 2205.01715, J. Ng et al. 2022, "Electron-scale reconnection in three-dimensional shock turbulence".
- 2001.04945, K. V. Lezhnin et al. 2020, "Kinetic simulations of electron pre-energization by magnetized
  collisionless shocks in expanding laboratory plasmas".
- 2107.00508, O. Kobzar et al. 2021, "Electron Acceleration at Rippled Low-Mach-number Shocks in High-beta
  Collisionless Cosmic Plasmas".
- 1602.08170, X. Kong et al. 2016, "Electron Acceleration at a Coronal Shock Propagating Through a Large-scale
  Streamer-like Magnetic Field".
- 2112.12114, B. D. Keenan et al. 2021, "Hybrid Particle-in-Cell Simulations of Electromagnetic Coupling and
  Waves From Streaming Burst Debris" (prior work of the engine's blast-wave example).
