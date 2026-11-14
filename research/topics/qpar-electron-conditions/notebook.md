# Notebook — qpar-electron-conditions

Dated lab notes, newest first. Write down what was tried, what was seen, and what was decided, failed attempts included.
Never rewrite a past entry. To correct one, add a new entry that refers to it.

## 2026-11-14 — Step 3, second part: what a simulation can state, and the boundary areas

Later on the same day as the entry below. The first part of step 3 collected the observed side only. The research
is done by simulation, so the question "what can a simulation state here, and what not" has to be surveyed with the
same care, together with the boundary areas. This entry records that survey. Result file: `simulation-claims.md`.

### What was done
- **Extraction from what was already read.** The notes and stored texts of the fully kinetic works on quasi-parallel
  shocks (`gupta2024electron`, `gupta2025speed`, `park2015simultaneous`, `kato2015particle`, `shalaby2022mechanism`,
  `zekovic2025slams`, `fang2019early`, all on `research/spherical-shock-injection`; `fang2022kinetic`,
  `yu2024nonthermal`, `bessho2023electron` here) were gone through for four things: what each states, on what basis,
  which limits its authors name, and how it is tied to an observation.
- **A second search, 26 keyword queries on OpenAlex** (10-12 results each, with abstracts):
  - methods and numerics: "particle-in-cell collisionless shock simulation mass ratio dependence electron
    acceleration"; "numerical heating particle-in-cell shock simulation number of particles per cell electron
    acceleration"; "shock simulation injection method reflecting wall piston comparison particle-in-cell";
    "three-dimensional particle-in-cell quasi-parallel shock electron"; "finite grid instability numerical heating
    drifting cold plasma particle-in-cell"; "simulation box transverse size cosmic ray current filamentation shock
    convergence"; "semi-implicit particle-in-cell collisionless shock realistic mass ratio";
  - neighbouring methods: "hybrid simulation test particle electrons quasi-parallel shock acceleration"; "MHD-PIC
    shock electron acceleration injection"; "magnetohydrodynamic-particle-in-cell method cosmic rays shock Bai";
    "Vlasiator electron foreshock hybrid-Vlasov electrons"; "eVlasiator electrons hybrid-Vlasov global simulation";
    "Vlasov-Fokker-Planck shock precursor cosmic ray electron injection simulation"; "continuum Vlasov-Maxwell
    collisionless shock electrons simulation field-particle correlation"; "dHybridR hybrid relativistic ions shock
    simulation"; "focused transport electron shock acceleration interplanetary model relativistic electrons"; "free
    escape boundary diffusive shock acceleration maximum energy time-dependent finite size electrons heliospheric
    shock"; "minimal model ion injection quasi-parallel shocks reformation";
  - runs against events: "global hybrid simulation foreshock transient electron acceleration test particle";
    "comparison particle-in-cell simulation MMS observations quasi-parallel shock electrons"; "interplanetary shock
    particle-in-cell simulation observed parameters electron acceleration Solar Orbiter Wind"; "Saturn bow shock
    particle-in-cell simulation Cassini high Mach number";
  - laboratory and border regimes: "laboratory magnetized collisionless shock electron acceleration experiment";
    "shock obliquity dependence electron acceleration particle-in-cell quasi-parallel oblique transition"; "weak
    quasi-parallel shock electron preacceleration intracluster particle-in-cell critical Mach number"; "electron
    injection shock speed dependence scaling non-relativistic shock simulation extrapolation realistic shock
    velocity".
- **Screening and reading.** Abstracts of about 90 works were kept and screened. 29 notes were written:
  - read in the body (arXiv PDF text; the read range is in each note): `zhou2025planar`, `marcowith2020multi`,
    `simon2026faux`, `bai2015magnetohydrodynamic`, `vanmarle2022diffusive`, `haggerty2019dhybridr`,
    `caprioli2018diffusive`, `caprioli2015simulations`, `orusa2025criteria`;
  - read as abstracts: `swisdak2023comparison`, `kropotina2023weibel`, `zhang2024kinetic`, `tsubouchi2004full`,
    `johlander2022quasi`, `liu2017fermi`, `nyberg2026electron`, `john2024effects`, `battarbee2021vlasov`,
    `lezhnin2021kinetic`, `kang2019electron`, `reville2013universal`, `bessho2025electron`,
    `liseykina2019electron`, `kylivnyk2026hybrid`, `gingell2017mms`, `schaeffer2019direct`,
    `vanthieghem2024electron`, `lindberg2025fermi`, `fiuza2020electron`.
- **Sources that did not work:** no open full text for `swisdak2023comparison`, `zhang2024kinetic`,
  `tsubouchi2004full`; our PDF text extractor gave nothing for arXiv:2310.13779 (`kropotina2023weibel`); the journal
  page of `nyberg2026electron` could not be read. All were read as abstracts and are marked so.
- **Passages re-read in the stored texts:** `kato2015particle`, the drag on energetic electrons in 1-D (their Eq. 4
  and the particle numbers at which acceleration fails); `gupta2024electron`, the definition of M_A, the relation of
  piston and shock speed, the frame, the dependence of spectra and efficiency on the shock speed; `gupta2025speed`,
  p_max in time.

### What the survey shows (detail and sources in `simulation-claims.md`)
- No published simulation in what was read measures D̃_e or an electron injection rate at the conditions of our four
  cases, or confronts one of the events with numbers from self-consistent electron kinetics. This repeats the
  result of the first survey with a wider base; it is still not exhaustive.
- The limits that authors themselves name recur in eight kinds: dimension, mass ratio, v_sh/c, run time, box and
  boundaries (with the launching method and the frame), particle noise, extrapolation by a formula, and ties to
  observations that are indirect (a global ratio at a remnant, a slope, a mechanism seen in situ).
- Published ways of tying a run to an observation: a global ratio at remnants; a synthetic probe through the box;
  hybrid runs at event parameters; test electrons in prescribed or hybrid fields; a model in prescribed fields
  fitted to one event; a mechanism compared with in-situ data; transport from data; laboratory shocks; Monte Carlo
  with prescribed inputs.
- The methods next to PIC either have no kinetic electrons or take the electron injection (or the fields at
  electron scales) as an input.
- Facts that bear on statements of the plan:
  - `gupta2024electron` already reports the electron injection efficiency for v_sh ≈ 0.067-0.267 c and mass ratios
    up to 1836: it "does not show much dependence on shock speed and m_R", and drops by a factor ~5 at low Mach
    number. The energetic electrons become fewer at high speed; the authors state that 1-D exaggerates this.
  - `gupta2024electron` runs in the upstream rest frame, chosen against numerical heating; its efficiency is a
    factor ≲ 2 above the downstream-frame run of `park2015simultaneous`. Our set-up is in the downstream frame.
  - `kato2015particle`: in 1-D, energetic electrons lose energy at a rate set by the number of electrons per
    electron skin depth; with too few, no acceleration appears.
  - `caprioli2018diffusive`: seeds drive waves through their current and can switch on the injection of thermal
    ions. Our seeds are a passive species.
  - `vanmarle2022diffusive`: "injection is, fundamentally, a 3D process" (ions, oblique shocks).

### Our calculations (own; not results)
- **Numerical drag against a Bohm-rate gain.** Loss from Kato's Eq. 4, gain from the formula of `gupta2025speed`
  with compression 4 and D_2 = D_1/10. Table by run in `simulation-claims.md` §5: loss/gain is 0.02-0.11 for W, V2,
  P0, P2; 0.3-0.5 for P5a, V1, E; above 1 for S1, P5b, P6, S2, at 5000 electrons per d_e and D̃ = 1.
- **A correction of our own first number.** The first evaluation for Gupta's run B1 used the piston speed (0.1 c) in
  place of the shock speed. `gupta2024electron` defines M_A with v_sh and sets v_sh = v_pt R/(R−1), so v_sh is
  0.133 c. With it the formula gives 2.1e-4 m_e c² ω_pe, against a mean rate of 2.0e-4 from the published p_max at
  300 ω_ci⁻¹. The earlier value (8.9e-5, never written to a file) is void.
- **A tension that is not resolved.** The same formula gives loss/gain ≈ 0.24 for Gupta's own B1, while that work
  reports agreement with its p_max formula and no sensitivity to the particle number. Either the 1-D loss formula
  overstates the drag in a filtered quasi-1-D run, or the effect is inside their agreement. Not checked; a short
  test run with the particle number varied would decide it.

### Changes of the README, each with its fact
| Fact | Change |
|---|---|
| The second search and its 29 notes | Related work: "Second search" added, with a pointer to `simulation-claims.md` |
| `gupta2024electron` reports the injection efficiency against shock speed and mass ratio (quotation above) | "Not new": that item added; what remains open is named |
| `liseykina2019electron`: test electrons in 1-D hybrid quasi-parallel shocks | "Not new": added beside `hanusch2020electron` |
| `kato2015particle` Eq. 4 and our estimate | "Strongest objections": numerical drag |
| `gupta2024electron` on the frame, against `park2015simultaneous` | "Strongest objections": the frame |
| `caprioli2018diffusive` | "Strongest objections": passive seeds cannot drive waves |
| `vanmarle2022diffusive` | "Strongest objections": injection as a 3-D process |
| `inputs.md`: v_sh/c is 0.0013-0.0020 at Earth and 0.0015 at Saturn | "Out of reach": the range corrected (it said 0.0013-0.0015 for both) |
| The points of `simulation-claims.md` §8 | Approach: a line that names them as open, for the owner; set-up, definitions and criteria unchanged |
| The second review is done; this survey is done | Status: "Done" and "Next" updated |

`literature.md`: five Map bullets for the 29 works and a pointer to `simulation-claims.md`.

### Second review by `research-reviewer`
Done on README, `inputs.md` and the plan. It lists nine corrections to `inputs.md` (attributions, a crossing time,
the kind of some values), finds that the published M_A values of the PSP event straddle M_A ≈ 10, that three of the
four "no" classes are limited by detection according to their own authors, and that part of the dependence on shock
speed is published; it proposes thresholds for the criteria. None of it is applied in this commit: the corrections
are to be checked in their sources first, and the thresholds are proposals until the owner decides.

### Not done
- The search is keyword-based and not exhaustive; no forward-citation pass was made for the 29 new works.
- 20 of the 29 works were read as abstracts only.
- The engine's implicit scheme as a route to larger cells was not examined (`zhou2025planar` is the one work read
  on semi-implicit PIC for shocks).
- The drag estimate is not tested by a run.
- The reviewer's corrections to `inputs.md` and the instrument backgrounds as intensities.

### Next
The owner's decisions on the open points (`simulation-claims.md` §8) and on the criteria; the corrections to
`inputs.md`; then the engine pilot (step 4).

## 2026-11-14 — Step 3, first part: the inputs of the four cases from primary sources

Later on the same day as the entry below.

### What was done
- The values the runs and the transport model need were taken from the stored full texts and from five further
  works, and written to `inputs.md` with source and kind (measured, derived, model, assumed, own).
- Passages re-read for this: `jebaraj2024acceleration` (A), `jebaraj2026wave`, `kouloumvakos2026wave`,
  `jebaraj2026diffusive`, the table of `chiappetta2026evolution` (PSP); `masters2013electron`, `masters2017situ`,
  `masters2016suprathermal` (Saturn); `wilson2016relativistic`, `shi2025compound`, `liu2019relativistic`,
  `liu2026efficient` (A), the table of `oka2025maximum` (Earth, Saturn); table 1 of `kartavykh2025statistical`
  (now added to its note), `dresing2016efficiency`, `talebpoursheshvan2023role` (1 au).
- New notes: `lalti2022database` (body, arXiv PDF text); `dresing2022role`, `masters2011electron`,
  `achilleos2006orientation`, `went2011new` (abstracts; full texts not obtained).
- This closes the item "Open from the plan's list of unread works" of the entry below: Dresing et al. 2022 is
  identified (ApJL 925, L21) and read as an abstract.
- Looked at, no note: Sulaiman et al. 2015 (PRL 115, 125001; quasi-perpendicular high-M_A shocks at Saturn; abstract).

### What the sources say about the kind of each input
- **PSP:** B, density and T_i are measured. T_e is not ("≲ 40 eV, unresolved"; 20 ± 5 eV is adopted from earlier
  orbits). M_A is 7-7.5, 9.1 ± 1.35, 4.68 ± 3.45 or ~4. Reasons given by the authors: a coarser estimate (9.1); an
  underestimated upstream density from SPAN-i (4.68; the compression rises from 1.39 to ~2.54 with the Langmuir-wave
  density). The hybrid model with M_A ~4 and ~1350 km/s was read as an abstract only.
- **Saturn:** only B (~0.1 nT), the times and the electron spectra are measured. M_A ~100 comes from a propagated
  solar wind model, θ_Bn ~20° from a model normal, M_s ~25 from typical temperatures, the flow speed (~450-500 km/s)
  and the shock motion (100 km/s) are typical or assumed values. No density, no β_e.
- **Earth:** the M_A range (most crossings at 5-15) is from 2797 MMS crossings. The limit without transients is the
  noise above 12 keV. No source read gives a connection time.
- **1 au:** averages of 51 Solar Orbiter shocks over all geometries (M_A 3.3 ± 0.3, gas compression 2.0 ± 0.1,
  β 0.84 ± 0.12). No values for quasi-parallel shocks alone.

### Our calculations (unverified; not claims)
- ω_ci: 7.66 rad/s at 80 nT; 9.58 × 10⁻³ rad/s at 0.1 nT; 0.479 and 0.766 rad/s at 5 and 8 nT.
- PSP: 4 h at 80 nT is ω_ci t ≈ 1.1 × 10⁵; β_e 0.04-0.08 for T_e 20-40 eV.
- Saturn: the quasi-parallel connection at the spacecraft lasts ≥ 60 min before the crossing (00:05 to 01:05 UT):
  ω_ci t ≥ 34. Size over speed (50-60 R_S, 300-600 km/s) gives 46-112. M_A 100 at 0.1 nT and 450 km/s implies a
  density of 0.24 cm⁻³. M_A 100 with M_s 25 implies a total β of ~19.
- Earth: size over speed (15-25 R_E, 400-600 km/s) at 5-8 nT gives ω_ci t 76-305. β_e ≈ 1 for 12 eV, 5 cm⁻³, 5 nT.
- 1 au: 1-3 days at 5 nT is ω_ci t 4 × 10⁴ - 1.2 × 10⁵.

### Changes of the README, each with its fact
No decision of the plan (question, hypotheses, runs, proposed criteria) was changed.
- No source read gives Mach numbers for the 15 quasi-parallel shocks of `yang2019electron`; the measured values are
  M_ms 1.25-2.85 in other samples and M_A 3.3 ± 0.3 for 51 shocks of all geometries. → The case table no longer
  states "M_A ~1.5-3" for that subset.
- `lalti2022database`: most crossings of Earth's bow shock have M_A 5-15; `oka2025maximum` lists 400-600 km/s.
  → The Earth row carries the ranges; ~8 stays as the representative value of run E.
- Saturn's M_A and speed are model and typical values (above). → Said so in the case table and in the list of
  objections.
- The connection estimates of the entry below (Earth 120-240, Saturn 30-70) are replaced in the README by the values
  above (76-305; ≥ 34 and 46-112). The planned run times (E: 240, S1-S3: 50) lie inside them.

### Not done
- Instrument backgrounds as intensities for each case (needed for H8 and for the confrontation by intensity).
- Still not found: the Gupta et al. 2025 conference contribution; the journal reference of `jebaraj2026diffusive`; the
  exact title of `galishnikova2026self`.
- A caution for later reading: each entry of table 1 of `kartavykh2025statistical` has a second value that is not the
  median of the four listed M_A in our check; only averages and standard errors are used.

### Next
- Re-review of the README, `inputs.md` and the plan by `research-reviewer`; the owner's decision on the criteria and
  thresholds.

## 2026-11-14 — Topic started: decisions, survey, differences from the approved plan

### Decisions (owner)
- Theme: electron acceleration at strong quasi-parallel shocks (Saturn/Cassini, Parker Solar Probe, SN 1006).
- Question: what decides whether a quasi-parallel shock accelerates electrons to relativistic energy.
- The plan (question in operational form, hypotheses H1a-H8, runs, proposed criteria) was approved after one review
  by `research-reviewer`. The criteria and null tests stay proposals until the review of step 3 and the owner's
  approval.
- Earlier candidates of this session (four themes, and selection criteria set without the owner) were rejected by the
  owner as too small and driven by the simulation instead of by observations. They are not reused.

### Survey (plan step 1)
- **Forward citations** (OpenAlex and INSPIRE) of 15 seed works; 337 citing works after merging. Seeds:
  `masters2013electron`, `masters2016suprathermal`, `masters2017situ`, `jebaraj2024acceleration`,
  `jebaraj2024direct`, `jebaraj2026wave`, `kouloumvakos2026wave` (no citing work yet), `liu2019relativistic`,
  `wilson2016relativistic`, `raptis2025revealing`, `gupta2024electron`, `gupta2025speed`, `guo2015acceleration`,
  `shalaby2022mechanism`, `zekovic2025slams`.
- **Venue sweep** (OpenAlex, title and abstract, 2016-2026): ApJ, ApJL, ApJS, JGR Space Physics, GRL, MNRAS, A&A, PRL,
  Physics of Plasmas, Nature Communications. Eight queries: electron acceleration with "quasi-parallel" / with
  "parallel shock"; "electron injection" with shock; electron "diffusion coefficient" with shock and kinetic or
  particle-in-cell; "relativistic electrons" with bow shock or foreshock; energetic electrons with "interplanetary
  shock"; electron "diffusive shock acceleration" with particle-in-cell or hybrid; electron, foreshock and
  acceleration with Saturn, Jupiter, Mars, Mercury or Venus. 84 works not already among the citing works.
- **Reading:** all titles screened; 50 abstracts; the text of about 22 works in part or whole. Each note states its
  range in its `Read:` line. 37 notes were written (the files added to `research/literature/` together with this entry).
- **Result:** no work was found that measures the electron injection rate or D̃_e at the conditions of an observed
  event from self-consistent kinetics, or that confronts an event with a transport model built on such a measurement.
- **Limits of the search:** most of the 84 sweep works were screened by title only. ADS and Semantic Scholar did not
  answer (HTTP 405 and 429). Publisher pages of IOP, A&A (HTML), Elsevier and PubMed Central refused automated
  reading; where no arXiv or repository version existed, the abstract from OpenAlex or Crossref was read. Installing a
  PDF library was not permitted; PDF text was extracted with a small script of our own, which loses digits in some
  files (`grassi2023electron`: only statements checked in context are quoted).
- **Not found:** a 2025 conference contribution of Gupta et al. named in the plan; the journal reference of
  `jebaraj2026diffusive` (arXiv:2610.08034); the exact title of `galishnikova2026self` (arXiv:2611.15879), to be
  checked before citing. The journal abstract of `bebesi2019observations` (a conference abstract was read instead).
- **Looked at, no note** (off the question): arXiv:2303.10969 (Dresing et al. 2023), arXiv:2409.03174 (Guo et al.
  2024), arXiv:2007.09194; Do et al. 2025 (ApJ 979, 50; time-dependent acceleration and escape of protons at
  travelling shocks near the Sun; abstract); Burkholder et al. 2025 (10.3847/1538-4357/ada440; foreshock
  fluctuations against Mach number at planetary bow shocks; abstract); Yang et al. 2018 and Qin et al. 2023 (shock
  drift at interplanetary shocks; abstracts). `winner2020evolution` and `morris2022preacceleration` have notes on
  `research/spherical-shock-injection`.
- **Open from the plan's list of unread works:** the plan names "Dresing et al. 2022" (electron intensities against
  the Mach number of coronal shocks). No note was written and this notebook cannot say that it was read. To do at
  step 3.

### Facts from the survey that bear on the hypotheses
Stated without ranking. Sources and ranges are in the notes.
- H1a/H1b: for protons at the PSP event, a steady balance holds only near 1.3 MeV and confinement fails between 2.3
  and 6 MeV; the age limit (~8.9 MeV) and the size limit (~9.5 MeV) lie above (`jebaraj2026diffusive`, preprint).
- H1b: the observed maxima often reach the Hillas limit (`oka2025maximum`); a Bohm scaling with the size matches
  Earth, Saturn and Jupiter if dB/B_0 ≫ 1 and enough time are assumed (`raptis2026relativistic`).
- H2/H3: mass ratio 100 suppresses electron acceleration at M_A 5.3 and inflates it at M_A 21.3 in 1-D parallel shocks
  (`shalaby2024energy`). Runs at M_A ~17-20 with mass ratio 30 show electron tails (`fang2022kinetic`,
  `yu2024nonthermal`).
- H4: seeds before every interplanetary shock (`tsurutani1985acceleration`); > 500 keV at Earth only with 1-5 keV seeds
  (`raptis2025revealing`); seeds to ~100 keV at Saturn (`masters2013electron`).
- H5: at Earth, a model without transients stops at 10-20 keV and with them reaches ≥ 200 keV, at a
  quasi-perpendicular shock (`shi2025compound`); a planar interplanetary shock lacks mature SLAMS
  (`raptis2026compressive`); persistent SLAMS need M_A ≳ 20 (`hemler2025shock`, abstract).
- H6: the 4 electron responses among 51 shocks have θ_Bn 88°, 42°, 80°, 61° (`kartavykh2025statistical`); reconnection
  in the transition energizes electrons in a 2-D run at Earth's parameters (`bessho2023electron`, abstract).
- H7/H8: the PSP crossing is a flank crossing with a saturated low-energy instrument (`dresing2025reason`,
  `kouloumvakos2026wave`); its M_A is 9.1 ± 1.35, 7-7.5, 4.68 ± 3.45 or ~4 depending on the work
  (`jebaraj2024acceleration`, `jebaraj2026wave`, `chiappetta2026evolution`, `giacalone2025source`).
- On the long-term state: "largely insensitive to the details of CR injection" (`galishnikova2026self`, abstract).

### Differences between the README and the approved plan
Each line names the fact first. No decision of the plan (question, hypotheses, runs, proposed criteria) was changed.
- `wilson2016relativistic` is now read (the plan listed it as unread): outside the disturbances THEMIS shows "only
  noise > 12 keV". → The Earth "no" is written as an instrument limit.
- `shi2025compound` confronts an Earth event with a data-constrained transport model; `raptis2026relativistic` gives a
  size scaling over three planets; `hemler2025shock` states the strategy "extend PIC by test particles";
  `grassi2023electron` (the LLNL report of the plan) and `hanusch2020electron` measure or follow electrons in kinetic
  or hybrid fields. → All five are added to the "not new" list.
- `gupta2024electron` varied the piston speed (0.05, 0.1, 0.2 c) for spectra and efficiency; `gupta2025speed`
  measured D_e at one speed. → Gap 2 reads "D̃_e is measured at one shock speed", not "the dependence on v_sh/c is
  unmeasured".
- The statement "> 500 keV only with 1-5 keV seeds" and the Bohm length of 10-100 R_E are in `raptis2025revealing`
  (Nat. Commun. 16, 488), not in arXiv:2411.12815 (`raptis2025multi`). → Citations corrected; the size limit is cited
  to `oka2025maximum` and `raptis2026relativistic`.
- `kartavykh2025statistical`, body: the table has 58 forward shocks, the electron search covers 51, and 4 respond.
  → "4 of 51" replaces "about 8 % of 58".
- `giacalone2025source` (abstract) and `chiappetta2026evolution` give M_A ~4 and a speed of ~1350 km/s for the PSP
  crossing. → The case table carries the spread; fixing the input is part of step 3.
- `bessho2023electron` (abstract): reconnection in the shock transition. → Named under H6 as a further non-DSA
  process and under the limits of the quasi-1-D set-up.
- `jebaraj2026diffusive`: wave growth limits proton confinement at the PSP event. → Written as an open point under
  the question; no candidate was added.
- DOI 10.3847/2041-8213/aeaf56 belongs to `kouloumvakos2026wave` (arXiv:2610.10067).

### Estimates (ours, unverified; not claims)
The p_max formula of `gupta2025speed` with D_1 = D_Bohm and D_2 = D_1/10, compression 4 (the reviewer reproduced the
values). Conventions (1/3, 1/2, the non-relativistic part of the integral) leave a factor 2-4.

| Case | Input (u_1, ω_ci t) | Reached with Bohm | Observed | Reading |
|---|---|---|---|---|
| Earth | 400 km/s, 240 | ~40 keV | none above the noise without transients | inside the factor 2-4: not decided |
| PSP | 2400 km/s, ~10⁵ (4 h, 80 nT) | 6.7 MeV at ω_ci t ≈ 700 | ≥ 6.7 MeV | H1 suffices if D̃ ≲ ~150; the measured proton κ_eff is ~180 u r_L at 2 MeV (our calculation) |
| Saturn | 450 km/s, 30-70 | ~10 keV | ≥ ~1 MeV | two orders short with H1 alone |
| weak shock at 1 au | 150 km/s, ~10⁵ | a few MeV; ~20 keV at 100 × Bohm | mostly 2-10 keV | no conflict if D̃ is two orders larger; needs compression and intensity (H8) |

- The largest energies of the two "yes" cases both lie at p ≈ 1.0 m_i v_sh (Saturn 2.78 against 2.75 m_e c; PSP 14.1
  against 14.7 m_e c). Two readings: the instrument limits happen to lie there; or the stall scale p ≈ m_i v_pt of
  `gupta2024electron`.
- Connection times can change by a factor 2.5-3 with the geometry of the field-line sweep (the reviewer's judgement).

### Next
- Step 3: the input table from primary sources; re-review of the README and the plan; the owner's decision on the
  criteria and thresholds.
