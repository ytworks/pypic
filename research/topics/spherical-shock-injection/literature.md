# Prior work — spherical-shock-injection

This file places the topic in its field. Each work has its own file in `research/literature/<citekey>.md`; here they are related to the question.

**Survey:** `survey-snr-electron-injection.md` (2026-11-12) covers the electron injection problem at supernova
remnant (SNR) shocks, with 54 works read and full references. This file is the short map and the reading list.
A one-page summary is in `survey-snr-electron-injection-summary.md`.
Japanese versions: `survey-snr-electron-injection.ja.md` and `survey-snr-electron-injection-summary.ja.md`; the
English files are the originals.
Section numbers below refer to the survey.

## Map
- **What observations require** (Sect. 3): `vink2012supernova`, `reynolds2004microphysics`, `morlino2012strong`.
  - Obliquity in whole remnants: `bocchino2011constraints`, `reynoso2013radio`, `west2017connection`,
    `zhou2023magnetic`, `giuffrida2022supernova`, `mcginness2026sn`, `winner2020evolution`.
  - Every global model assumes its efficiency law ε(θ_Bn); SN 1006 and the radio population point in opposite
    directions.
- **Quasi-perpendicular shocks** (Sect. 4.1):
  - surfing on Buneman waves: `hoshino2002nonthermal`, `amano2007electron`, `amano2009electron`,
    `matsumoto2012electron`, `matsumoto2013electron`, `matsumoto2017electron`; against it: `ohira2007absence`;
  - shock drift with scattering: `katou2019theory`, `amano2020observational`, `amano2010critical`,
    `amano2022theory`, `wang2026generation`;
  - whistlers: `riquelme2011electron`;
  - Weibel turbulence, reconnection, stochastic scattering: `kato2010nonrelativistic`,
    `wieland2016nonrelativistic`, `bohdan2017electron`, `bohdan2019kinetic`, `bohdan2019kineticb`,
    `bohdan2020kinetic`, `bohdan2021magnetic`;
  - electron-driven waves and diffusive acceleration in 1-D: `xu2020electron`, `kumar2021nonthermal`.
- **Oblique shocks and the electron foreshock** (Sect. 4.2): `morris2022preacceleration`,
  `morris2023preacceleration`, `fulat2023kinetic`, `fulat2026pic`.
- **Quasi-parallel shocks and injection thresholds** (Sect. 4.3): `park2015simultaneous`, `kato2015particle`,
  `fang2019early`, `crumley2019kinetic`, `shalaby2022mechanism`, `arbutina2021nonlinear`, `gupta2024electron`,
  `gupta2025speed`, `zekovic2025slams`, `jikei2026dependence`.
- **Geometry beyond a planar uniform shock** (Sect. 6):
  - where protons drive waves: `caprioli2014simulations`;
  - ripples and upstream turbulence that change the local θ_Bn: `burgess2006simulations`, `trotta2019electron`,
    `guo2010effect`, `guo2015acceleration`;
  - fronts curved on the large scale (read 2026-11-13; see "Curved fronts" below; Sect. 6.5 of the survey was
    written before this reading).
- **Reviews:** `bohdan2023electron`, `amano2022nonthermal`, `marcowith2016microphysics`, `pohl2020pic`.
- **Benchmarks outside the SNR regime:** `raptis2025revealing`, `guo2014nonthermal`.
- **Second-order Fermi acceleration in islands**, for reconnection: `hoshino2012stochastic`.

## Still to read
- **Curved shocks: full texts not obtained** (publisher hosts refused the requests; no open copy found). Known by
  abstract only, so every statement about them here is provisional:
  - `kraussvarban1991electron` (the first treatment of electrons at a curved shock) and
    `giacalone2017acceleration` (protons at a spherical blast wave): the two that matter most for novelty;
  - `kraussvarban1989electron`, `kraussvarban1994electron`, `vandas1995acceleration`, `vandas1995accelerationb`,
    `vandas2001shock`, `anderson1979thin`, `wu1984fast`, `holman1983solar`, `knock2001theory`,
    `vandas2004acceleration`, `zhou2026crossscale`;
  - not obtained at all: Leroy & Mangeney 1984 (Ann. Geophys. 2, 449); Savoini & Lembège 1999 (Adv. Space Res.
    24, 13); the time-of-flight papers (Filbert & Kellogg 1979; Cairns 1986, 1987), known through citations.
- **Wakefield acceleration in relativistic shocks:** arXiv:0709.1783 (Hoshino 2008), arXiv:1909.03337. The
  survey did not cover relativistic shocks.
- **Spherical or expanding shocks (the driver)**, prior work of `examples/blast_wave`:
  - Le et al. 2021 (arXiv:2109.00583), hybrid;
  - Keenan et al. 2022 (arXiv:2112.12114), hybrid;
  - Dieckmann et al. 2018, 2-D PIC;
  - Lezhnin et al. 2020 (arXiv:2001.04945), electron pre-energization in expanding laboratory plasmas.

  To check: what these show about electrons.
- **Found by the search and not read:** the list in the survey, Sect. 10.3.
- **Not on arXiv:** Matsumoto et al. 2015 (Science), Rothenflug et al. 2004, Fulbright & Reynolds 1990.

## Open
From the survey, Sects. 8.2 and 8.3:
- The dependence of electron injection on θ_Bn is known piecewise only, from studies with different mass ratios,
  Mach numbers, dimensions and definitions.
- The global models of remnants need that dependence and assume it.
- No kinetic simulation of electron injection on a front curved on the large scale was found for SNR-like
  parameters. Two full-PIC studies exist at bow-shock parameters, in 2-D (`savoini2001twodimensional`,
  `savoini2010origin`); see "Curved fronts" below. Two key works are known by abstract only.
- Which mechanism dominates at the real mass ratio and shock speed is disputed.

## Curved fronts (read 2026-11-13)
The direct predecessors of the exploratory runs. Reading level: full for `savoini2001twodimensional`,
`savoini2010origin`, `vandas1994acceleration`, `vandas2001shockb`, `kong2016electron`; skimmed for
`pulupa2008structure`, `pulupa2012stereo`, `kong2015possible`, `lezhnin2021kinetic`, `dieckmann2018expansion`,
`liu2026efficient`, `savoini2013origin`, `savoini2015production`, `savoini2020deep`; abstract only for the works
listed under "Still to read".

- **Already published (no novelty in these statements):**
  - the source of reflected electrons lies just behind the tangent field line, with an energy set by
    V_sh/cos θ_Bn: `anderson1979thin`, `wu1984fast`, `holman1983solar`, `vandas2001shock`, `kong2016electron`
    (injection at 74°–87° for an expanding circular front, none above 86.5°);
  - θ_Bn at the foot of one field line changes in time as a curved front expands:
    `savoini2001twodimensional` (Sect. 3.4), `savoini2010origin`, `kong2016electron`, and the time-of-flight
    literature;
  - a curved front lowers the gain of one reflection and focuses the reflected flux:
    `kraussvarban1991electron` (abstract), `vandas1994acceleration` (mean gain about 1/4 of the plane shock for
    R_c/d = 10²–10³);
  - an electron that enters near 88° and leaves near 48°, in full PIC: `savoini2010origin` ("leaked"
    electrons, about 10 % of the backstreaming ones; the mechanism given is bouncing in the downstream
    fluctuations);
  - without scattering the gain is confined to 5°–10° from the tangent point; with scattering, higher phase-space
    densities at lower θ_Bn than a local mapping gives: `kraussvarban1994electron` (abstract);
  - several encounters and collapsing traps, where a bent field line or a non-convex front closes the path:
    `kong2015possible`, `kong2016electron`, `liu2026efficient`, `pulupa2008structure`, `burgess2006simulations`;
  - protons at a spherical blast wave gather at the poles after a first encounter at the equator:
    `giacalone2017acceleration` (abstract).
- **Full-PIC studies of electrons at a globally curved front:** `savoini2001twodimensional` and
  `savoini2010origin` only (2.5-D, cylindrical piston, m_i/m_e = 42, M_A = 3, θ_Bn from 90° to 65° and to 45°,
  runs of 1.1 and 2.2 ion gyro-periods); `zhou2026crossscale` (global 2-D, abstract only).
  `dieckmann2018expansion` (2-D radial blast) and `lezhnin2021kinetic` (planar) do not analyse electron
  acceleration at a curved front.
- **Not found in what could be reached:**
  - electrons in a full-PIC run of a spherical (3-D) expanding and slowing front;
  - a comparison of the energy of an escaping electron with the one-reflection energy at its exit θ_Bn;
  - a path across latitudes by drift inside the ramp, as opposed to leakage from downstream;
  - a run above m_i/m_e = 84 or M_A ≈ 5.
- **How firm this is:** weak until `kraussvarban1991electron` and `giacalone2017acceleration` are read in full.
  A forward-citation search was possible through OpenAlex only (ADS refused requests).
- **Corrections to the reference list of the survey (Sect. 10.2), kept here because the survey is not
  rewritten:** "Vandas 1995" is two papers (`vandas1995acceleration`, zero thickness; `vandas1995accelerationb`,
  non-zero thickness); "Vandas 2001" in JGR is `vandas2001shock`, and a second paper of 2001 is
  `vandas2001shockb`; `vandas1994acceleration` is an earlier, open paper; Schmidt & Cairns 2014 (JGR 119, 69) is
  marked retracted in OpenAlex and is not to be cited without a check.

## This topic
To be written once the research question is decided (the owner's decision, open since the review of
2026-11-12). The gaps relative to the topic are in the survey, Sect. 8.3.
