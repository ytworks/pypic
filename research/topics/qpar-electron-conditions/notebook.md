# Notebook — qpar-electron-conditions

Dated lab notes, newest first. Write down what was tried, what was seen, and what was decided, failed attempts included.
Never rewrite a past entry. To correct one, add a new entry that refers to it.

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
