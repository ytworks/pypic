# What decides whether a quasi-parallel shock accelerates electrons to relativistic energy

Branch: `research/qpar-electron-conditions`. Started: 2026-11-14. Notebook: `notebook.md`. Survey: `literature.md`.
Inputs: `inputs.md`.
Sections follow `.claude/rules/research/stance.md`. Mark a section as open rather than guessing.

Theme and question were fixed by the owner on 2026-11-14. The plan was approved by the owner on the same day and
reviewed once by `research-reviewer` before approval. Citekeys marked (S) have their note on
`research/spherical-shock-injection`, (A) on `research/shock-alfvenic-turbulence`; neither is merged into trunk yet.

## Background
- **Diffusive shock acceleration (DSA) and injection.** Quasi-parallel shocks accelerate ions by DSA in turbulence
  the ions drive themselves (`caprioli2014simulationsb`). Electrons must first be injected: brought from the thermal
  pool to a momentum from which they can cross the shock repeatedly.
- **The view before 2013.** In the solar system, electron acceleration at quasi-parallel shocks was held to be
  "effectively absent", and MeV electrons were rare even at quasi-perpendicular shocks (`masters2013electron`).
- **Observations since then do not fit one rule.**
  - Saturn, Cassini, 2007-02-03: quasi-parallel (θ_Bn ~20°), M_A ~100; electrons to ~1 MeV, the top of the
    detector range (`masters2013electron`). Of 508 usable crossings 31 show a signal and 3 reach ~1 MeV
    (`masters2016suprathermal`). A quasi-perpendicular crossing at M_A ~15 has the same index above 100 keV within
    the errors (`masters2017situ`).
  - Parker Solar Probe (PSP), 2023-03-13, 0.24 au: nearly parallel (θ_Bn ~8°); electrons isotropic to 6.7 MeV, near
    the top of the instrument range (`jebaraj2024acceleration` (A)). M_A is 7-7.5 in the later analysis
    (`jebaraj2026wave`); other works give ~4 for the same crossing (`giacalone2025source`,
    `chiappetta2026evolution`).
  - Earth: electrons above 100 keV appear with foreshock transients (`liu2019relativistic`, `raptis2025multi`,
    `raptis2025revealing` (S)). Outside them THEMIS shows "only noise > 12 keV" (`wilson2016relativistic`). A
    data-constrained model of an Earth event stops near 20 keV without transients (`shi2025compound`).
  - Jupiter, Juno, 2023-10-01: ≥ 1 MeV electrons, inside a foreshock transient (`raptis2026relativistic`).
  - Interplanetary shocks near 1 au: electron effects are rare and belong to oblique or quasi-perpendicular shocks.
    Shock-associated 65-75 keV increases at 5 of 475 crossings, one of them quasi-parallel
    (`dresing2016efficiency`); response at 4 of 51 shocks, θ_Bn 42-88° (`kartavykh2025statistical`); MeV
    increases at 9 of 587, θ_Bn 45-88° (`talebpoursheshvan2023role`); effects mostly limited to 2-10 keV
    (`tsurutani1985acceleration`); at 40 keV ~7 % of 15 quasi-parallel shocks (`yang2019electron`).
  - SN 1006 (remote): the field just behind the north-eastern shock is nearly parallel to the normal; which geometry
    is the efficient one is called a long-standing question (`zhou2023magnetic` (S)).
- **Simulations.** Fully kinetic runs of quasi-parallel shocks give DSA electron tails for M_A ≳ 10 and a stall below
  (`gupta2024electron` (S)); a 1-D parallel study finds the opposite order (`shalaby2022mechanism` (S)). The electron
  diffusion coefficient has been measured in such runs at two conditions (`gupta2025speed` (S)).

## Research question
Fixed by the owner on 2026-11-14:

> What decides whether a quasi-parallel shock accelerates electrons to relativistic energy?

Operational form, the part a simulation can decide:

> With a uniform upstream and only the observed shock parameters given, do the injection rate and the electron
> diffusion coefficient in Bohm units, D̃_e(x, p̃), that the shock produces by itself reproduce, together with the
> observed age and size, the presence or absence and the spectrum of relativistic electrons in each event? Where they
> do not, what has to be added (seed electrons, pre-existing structure, the history of the shock)?

The reviewer's judgement, which this form makes explicit: the plan can answer the narrow question "do the limits of
time and size, with the measured D_e, suffice".

The cases. Every value, its source and its kind (measured, model, assumed, own) are in `inputs.md`:

| Case | M_A | Speed against the upstream | Relativistic electrons | Remarks |
|---|---|---|---|---|
| Earth, no transient | 5-15 (most crossings); ~8 as representative | 400-600 km/s | no (noise above 12 keV) | instrument limit, not a zero; no published connection time found |
| PSP 2023-03-13 | 7-7.5 (other estimates ~4, 9.1) | ~2400 km/s (other estimate 1350) | yes (≥ 6.7 MeV) | shock ~4 h old; flank crossing |
| Saturn 2007-02-03 | ~100 (solar wind model) | ~450 km/s (typical value, not measured) | yes (≥ ~1 MeV) | 3 of 508 crossings; seeds to ~100 keV upstream; density and shock motion not measured |
| 1 au quasi-parallel | not published for the subset; M_A 3.3 ± 0.3 and M_ms 1.3-2.9 in samples of all geometries | ~500 km/s (spacecraft frame, all geometries) | 40 keV increase at ~7 % of 15 | MeV cases at 1 au are oblique |

The outcomes do not follow the order of M_A.

Candidate explanations, all kept side by side:
- **H1a time:** injection occurs in every case; the energy reached is set by D_e and the time of connection.
- **H1b size:** electrons escape when the diffusion length exceeds the size of the shock (`oka2025maximum`,
  `raptis2026relativistic`).
- **H2 lower bound in M_A:** below M_A ~10 the self-generated turbulence is too weak (`gupta2024electron` (S)).
- **H3 upper bound in M_A:** electron-scale instabilities act only for M_A ≲ √(m_i/m_e)/4 (`shalaby2022mechanism` (S)).
- **H4 seeds:** only with keV-100 keV electrons upstream (`raptis2025revealing` (S), `tsurutani1985acceleration`,
  the magnetospheric electrons of `masters2013electron`).
- **H5 upstream structure:** pre-existing turbulence (`guo2015acceleration` (S)) or foreshock transients
  (`liu2019relativistic`, `shi2025compound`, `liu2026efficient` (A)) are needed.
- **H6 not DSA:** shock drift acceleration (`yang2019electron`, `kartavykh2025statistical`), trapping or
  quasi-periodic acceleration by SLAMS (`masters2013electron`, `zekovic2025slams` (S), `hemler2025shock`),
  reconnection in the transition (`bessho2023electron`).
- **H7 history:** the PSP shock travelled ~4 h from near the Sun; M_A, θ_Bn and B changed on the way
  (`dresing2025reason`, `chiappetta2026evolution`).
- **H8 detection:** a small compression gives a steep spectrum; "yes" or "no" is an intensity against an instrument
  background (`masters2016suprathermal`, `kouloumvakos2026wave` for the saturation at PSP).

What the runs can test: H1-H3 and part of H6 directly; H4 with seed electrons as a test species; H8 through the
predicted intensity. H5 and H7 only in the form "is the shock alone enough"; if not, they stay outside this topic.
Nine candidates face four cases of different kinds, so a single answer is not assured: the paper states which
candidates the kinetic measurement excludes and which remain.

Open at the review of step 3: for protons at the PSP event, confinement is limited by wave growth, a factor of two
below the age limit and below the size limit (`jebaraj2026diffusive`, a preprint). Whether this is a separate
candidate for electrons, or is covered by measuring D̃_e in a growing precursor (null test N8), is not decided.

## Prior work
- **Direct predecessors:** `gupta2024electron` (S) and `gupta2025speed` (S).
  - 2024: quasi-1-D PIC at θ_Bn 30°, 300 ω_ci⁻¹, M_A 5-40, mass ratio 16-1836, piston speed 0.05-0.2 c, β 0.06-64.
    DSA slopes and growing tails for M_A ≳ 10; for M_A ≲ 10 the tail stops, "at least" within the run time. The run
    nearest to the PSP case (B2_II: M_A 5, mass ratio 1836, β 0.06) stops near p ≈ m_i v_pt.
  - 2025: a formula for p_max(t) from the upstream and downstream diffusion coefficients; D_e(p) measured with test
    electrons in field snapshots of the PIC runs, at M_A 5 and 20 (θ_Bn 30°, one piston speed). At M_A 20,
    D_e ≈ 0.6-1 Bohm just upstream. "A full characterization of weak and more oblique shocks" is left as future work.
  - Why it was left (our reading): computing time. The paper notes that tails grow slowly at low M_A.
- **Events already confronted with models:**
  - Saturn: test electrons in prescribed turbulence reach ~1 MeV for every θ_Bn when the variance is 10 B_0², and
    7 keV at θ_Bn 0° when it is 0.1 B_0² (`guo2015acceleration` (S)). The turbulence is an input.
  - PSP: for protons, the precursor and the transport are measured (`jebaraj2026wave`, `kouloumvakos2026wave`,
    `jebaraj2026diffusive`), and hybrid simulations are compared with the event (`giacalone2025source`).
  - Earth: a transport model with wave properties taken from data reproduces an event with transients; the shock
    there is quasi-perpendicular (`shi2025compound`).
  - Planets together: a size-limited Bohm scaling is shown consistent with the maxima at Earth, Saturn and Jupiter,
    with Bohm scattering and sufficient time assumed (`raptis2026relativistic`).
- **Where they stop:**
  1. D_e is measured at θ_Bn 30°, M_A 5 and 20 only; not at the conditions of an event (M_A 7.5 with θ_Bn 8° and
     β 0.2; M_A ~100; M_A ~3).
  2. D̃_e is measured at one shock speed. `gupta2024electron` (S) varied the speed (0.05, 0.1, 0.2 c) for spectra and
     efficiency, not for D_e. Real shocks are 1.2-1.8 orders of magnitude slower than the simulated ones.
  3. No work was found that confronts an event with electron transport taken from a self-consistent kinetic
     simulation: Saturn and Earth use prescribed fields, PSP has protons only.

## Idea
1. **Scaling leaves one quantity unmatched.** With M_A, β_e, β_i, θ_Bn and the mass ratio 1836 set to the event, the
   electron sonic Mach number is matched as well. Only v_sh/c is not. Where the event has no measurement (T_e at
   PSP; density, β_e and the speed at Saturn) the run takes a stated value from `inputs.md`.
2. **Measure the dependence on v_sh/c, then transfer.** D̃_e(x, p̃) and the injection rate are measured at three
   values of v_sh/c; the exponent and its error carry the result to the real v_sh/c through a transport model.
3. **The age of the short cases is within reach.** The connection times at Earth and Saturn can be run. Our
   estimates are ω_ci t 76-305 for Earth (size over speed; no published value found) and 46-112 for Saturn, where the
   spacecraft saw the quasi-parallel connection for ≥ 60 min before the crossing (ω_ci t ≥ 34) (`inputs.md`). What a
   run covers directly is the turbulence at that age, not the momentum the electrons reach: at equal ω_ci t and D̃ the
   formula of `gupta2025speed` (S) keeps p̃·u/v_pt, and simulated electrons have a smaller u/v_sh than real ones.
   Every case is therefore confronted through the transport model.
4. **PSP offers a measured precursor.** The proton κ_eff and the wave amplitudes are known (`jebaraj2026wave`); the
   precursor of the run can be compared with them directly.
5. Said plainly: most of the difficulty is computing time. The idea is in 1-4: bring the runs into reach and carry
   what is out of reach with a stated error.

## Related work and novelty
- **Search done (2026-11-14):** forward citations of 15 seed works (337 citing works), a sweep of 10 journals for
  2016-2026 with 8 queries (84 further works), screening by title, 50 abstracts, and the text of about 22 works in
  part or whole. Record and the list of failed sources: `notebook.md`. Limits: titles only for most of the sweep;
  four full texts could not be obtained and were read as abstracts.
- **Second search (2026-11-14), on the simulation side:** 26 keyword queries on OpenAlex for the methods next to
  PIC, for runs confronted with in-situ events, for laboratory shocks and for the border regimes; the abstracts of
  about 90 works screened, 29 further notes (9 read in the body, 20 as abstracts). It is not exhaustive. What published simulations
  state, with which limits, how they were tied to observations, and what the planned runs can and cannot state:
  `simulation-claims.md`.
- **Novelty claimed, in two parts:**
  1. The injection rate and D̃_e(x, p̃) of electrons, measured in self-consistent kinetic simulations at the
     conditions of observed events, with the dependence on v_sh/c. The range that is certainly new is
     m_e/m_i ≲ p̃ ≲ 1; above it an electron has the rigidity of an ion and D is largely known
     (`caprioli2014simulationsb`).
  2. A transport model fed with these measurements, confronted with the events that did and did not accelerate
     electrons, with errors, stating which of H1-H8 it excludes.
- **Not new, stated as such:**
  - D_e(p) at M_A 5 and 20, θ_Bn 30° (`gupta2025speed` (S)); the electron mean-square-displacement method in PIC
    (`grassi2023electron`, unmagnetized shocks); test electrons in hybrid fields (`hanusch2020electron`,
    `liseykina2019electron`).
  - The electron injection efficiency against shock speed and mass ratio at θ_Bn 30°: `gupta2024electron` (S) covers
    v_sh ≈ 0.067-0.267 c and mass ratios up to 1836 and reports that it "does not show much dependence on shock
    speed and m_R", with a drop by a factor ~5 at low Mach number. What remains open is the dependence of D̃_e on
    v_sh/c, the event values of M_A, θ_Bn and β, and speeds below 0.067 c.
  - Saturn reproduced with prescribed turbulence (`guo2015acceleration` (S)).
  - Proton transport at the PSP event (`jebaraj2026wave`, `kouloumvakos2026wave`, `jebaraj2026diffusive`,
    `giacalone2025source`).
  - The size limit as an idea and as a scaling over three planets (`oka2025maximum`, `raptis2026relativistic`).
  - Confronting an event with a data-constrained transport model (`shi2025compound`).
  - Extending PIC results with test particles to longer scales (`hemler2025shock`, prescribed SLAMS).
- **Reviewer's verdict before the survey:** part 1 "needs evidence, likely to hold after the survey"; part 2
  "needs rethinking" (answered by idea 3 and by confronting through the transport model only).
- **Strongest objections an expert would raise:**
  - The transfer in v_sh/c over more than one order of magnitude rests on three points.
  - A quasi-1-D run cannot hold the oblique whistlers seen at the PSP shock, nor reconnection in the transition
    (`bessho2023electron`).
  - The parameters of the PSP crossing are disputed by a factor of two (M_A ~4 to 9.1).
  - At Saturn, M_A (a propagated solar wind model), θ_Bn (a model normal) and the speed (a typical value) are not
    measured at the event; density and β_e are not given.
  - Mass ratio 100, needed at M_A 100, changes electron acceleration with a sign that depends on M_A
    (`shalaby2024energy`).
  - The long-term state may not depend on injection details (`galishnikova2026self`, a preprint).
  - Numerical drag on energetic electrons in 1-D PIC, dE/dt ~ −(1/2N_e) m_e c² ω_pe (`kato2015particle` (S)). Our
    estimate, unverified, puts it at a few per cent of a Bohm-rate gain for P0 and above the gain for runs P5b,
    P6, S1 and S2 at the planned particle number (`simulation-claims.md` §5).
  - The frame: `gupta2024electron` (S) runs in the upstream rest frame to avoid numerical heating and finds an
    efficiency a factor ≲ 2 above the downstream-frame run of `park2015simultaneous` (S). The set-up below is in the
    downstream frame.
  - Passive seed electrons cannot drive waves; in hybrid runs the current of reflected seeds changes the turbulence
    and the injection (`caprioli2018diffusive`, ions).
  - "injection is, fundamentally, a 3D process" (`vanmarle2022diffusive`, on ions at oblique shocks).

## Approach
Criteria and null tests below are proposals. They are fixed, with thresholds, after the review of step 3 and the
owner's approval, before any production run.

- **Set-up:** quasi-1-D (long in x, a few cells across), explicit scheme, a piston against a wall (downstream frame),
  inflow through the open x+ face (PIC-D13), Δx = d_e/10, ~1000 particles per cell column over both species,
  precision `compact` or `mixed`, restart checkpoints (PIC-D25). Outputs under `/mnt/nvme/pypic`.
- **Open points of the set-up after the second search** (the owner's to decide; nothing above is changed by them):
  the frame, the particle number, the definition of the injection rate and its region, how H4 is stated with
  passive seeds, and which border regimes belong to the question (`simulation-claims.md` §8).
- **Limits of 1-D:** no oblique whistlers; B_x cannot change, which overstates the superluminal cut-off
  (`gupta2024electron` (S)). Run X measures the difference to 2-D; the results name the limit.
- **Measured quantities:**
  1. Momentum spectra of both species upstream and downstream in time; p_max(t); compression.
  2. Injection rate: the share of electrons per unit time that pass p_inj, the momentum at which an electron can
     outrun the shock upstream (the condition of `gupta2025speed` (S)).
  3. D̃_e(x, p̃) by the method of `gupta2025speed` (S), appendix B (test electrons in saved fields), as a function of
     the distance from the shock up to the end of the box; independently from the tracked PIC electrons; with and
     without the wave electric field.
  4. Wavenumber spectra of the field, δB/B in space, precursor waves; at the PSP conditions against the measured
     amplitudes (total ~10 %, resonant 5-6 %).
  5. Energy gain of single electrons by field component, to separate diffusive gain from drift and trapping (H6).
  6. Seed electrons: a passive test species (PIC-D9) of keV-100 keV equivalent energy in the same run (H4).
- **Topic code (no engine change):** the D̃_e measurement; a transport model taking D̃(x, p̃) and the injection rate,
  with time, size (free-escape boundary) and compression, giving the spectrum j(E, t) as an intensity against each
  instrument's background; the comparison table.
- **Runs** (cost estimates for 8 GPUs, to be replaced by pilot timings; read them with a factor 1.5-2):

  | Stage | Run | Conditions | Purpose | ω_ci t | Estimate |
  |---|---|---|---|---|---|
  | 1 | V1 | Gupta B1: M_A 20, mass ratio 100, θ 30°, v_pt 0.1 c | reproduce slope, p_max(t), D_e; fix conventions | 150 | 0.9 d |
  | 1 | V2 | Gupta B2_II: M_A 5, mass ratio 1836, β 0.06 | reproduce the stall at low M_A | 150 | 1.0 d |
  | 2 | P0 | PSP: M_A 7.5, β_i 0.2, T_e 20-40 eV, θ 8°, mass ratio 1836, v_sh 0.133 c, seeds | injection, D̃_e, δB/B against data | 200 | 4.1 d |
  | 2 | P0b | P0, another random seed | scatter | 100 | 1.0 d |
  | 2 | P5a, P5b | P0 at v_sh 0.067 c and 0.033 c | exponent in v_sh/c (N4) | 100 | 2.0, 4.1 d |
  | 2 | P6 | V1 at half v_sh/c | the same at M_A 20 | 100 | 0.8 d |
  | 2 | P1, P3, P4 | P0 with one change: M_A 5 / θ 30° / β ~2 | which quantity matters | 100 | 0.5, 1.0, 1.0 d |
  | 2 | W | 1 au quasi-parallel: M_A 3, β ~1, θ 30°, mass ratio 1836 | D̃_e, compression, intensity | 300 | 1.5 d |
  | 3 | P2 | P0 at M_A 20 | H2, H3 | 100 | 7.2 d |
  | 3 | E | Earth: M_A 8, β ~1, θ ~25°, mass ratio 1836, v_sh 0.05-0.07 c, seeds | turbulence and D̃_e over the full age | 240 | 13-18 d |
  | 3 | S1 | Saturn: M_A 100, θ 20°, β ~20, mass ratio 100, v_sh 0.133 c | the regime of `zekovic2025slams` (S), as control | 50 | 2.5 d |
  | 3 | S2 | S1 at v_sh 0.033 c (v_sh B_⊥/B_0 < c), seeds | nearer to Saturn's real regime | 50 | 9.9 d |
  | 3 | S3 | S2 with the real electron sonic Mach number (~1) | separate mass ratio from M_se | 50 | 9.9 d |
  | 3 | X | P0 in 2-D (width and time from the pilot) | injection and D̃_e (p̃ ≪ 1) against 1-D | 30-60 | 7-26 d |

  Cost formula (ours, reproduced by the reviewer): particle-steps ≈ 180 N_x (m_i/m_e) M_A² (ω_ci t)² / (v_sh/c),
  at ~1.6 × 10¹⁰ particle-steps/s. Stage 3 may shrink or change after stage 2; each stage ends with a review.
- **Out of reach, stated:** the real v_sh/c (PSP 0.008; Earth 0.0013-0.0020; Saturn 0.0015); the real length of the
  electron precursor; Saturn with mass ratio ≥ 400 at a subluminal speed (Saturn claims are limited to mass ratio
  100); the age (ω_ci t ~10⁵) and history of the PSP shock.
- **Reproduction (V1, V2), proposed:** slope within ±0.3, p_max(t) and D_e within a factor 2 of the published values;
  conventions (D_B = u r_L, ⟨Δx²⟩/t) fixed in V1. No further stage until a deviation is explained.
- **Null tests, proposed:**
  - N1 numerical heating: a uniform flow without a shock for the same time gives the detection floor.
  - N2 statistics: twice the particles and another seed (P0b) leave the tail fraction and D̃_e unchanged.
  - N3 resolution: half Δx (short).
  - N4 extrapolation: exponents of D̃_e and of the injection rate in v_sh/c with errors, carried to the real
    v_sh/c. A case whose predicted margin is smaller than this error is reported as "not decided".
  - N5 mass ratio: one set at fixed electron sonic Mach number, one at fixed β_e (`shalaby2024energy`).
  - N6 boundary: twice the box length (short).
  - N7 method: the measurement code on synthetic fields with a known D; intervals where ⟨Δx²⟩ is not linear in
    time are not used.
  - N8 stationarity: D̃_e(p̃ ≤ 1) does not change in the second half of the run; otherwise the value is reported as
    "at that age" only.
- **Confrontation, proposed:** all comparisons are between the transport model's intensity and the instrument
  background. H1 survives if the prediction falls on the observed side in every case within the errors of N4 and of
  the inputs. H2, H3: the order of injection rate and D̃_e over P1, P0, P2, S. H4: seeds against thermal electrons in
  one run. H6: the energy-gain budget. H8: predicted compression and intensity.
- **Outcomes that leave the question unanswered, named in advance:**
  - P0 stalls near p ≈ m_i v_pt as B2_II does: D̃ above is not measured and nothing transfers to age 10⁵. The
    statement is then limited to "a uniform upstream does not grow the turbulence within 200 ω_ci⁻¹" and to the
    comparison of δB/B with the measured precursor.
  - The N4 error is so large that every case is "not decided": the claim is limited to novelty part 1.
  - 1-D and 2-D injection rates differ strongly: the 1-D transport results are bounds only.

## Claims and evidence
No records yet. Everything below is a hypothesis.

- (hypothesis) The injection rate and D̃_e(x, p̃) at the PSP conditions differ from the M_A 5 and M_A 20 values of
  `gupta2025speed` (S) by more than the run-to-run scatter. Records: none. Planned: V1, V2, P0, P0b.
- (hypothesis) D̃_e and the injection rate follow a power law in v_sh/c over 0.033-0.133 c. Records: none. Planned:
  P0, P5a, P5b, P6.
- (hypothesis) Time and size limits with the measured D̃_e place each of the four cases on its observed side.
  Records: none. Planned: all runs and the transport model.
- (estimate, not a claim) With Bohm scattering upstream and a tenth of it downstream, the formula of
  `gupta2025speed` (S) gives ~40 keV for Earth, 6.7 MeV for PSP at ω_ci t ≈ 700, ~10 keV for Saturn and a few MeV
  for a weak 1-au shock. Convention factors of 2-4 are open; the inputs are not fixed. Details: `notebook.md`.

## Publication
Open. Proposal in the plan: a peer-reviewed journal where the direct predecessors appeared (ApJ or ApJL) and arXiv.
The owner decides.

## Status
Snapshot of 2026-11-14.

- Done: theme and question (owner); plan approved (owner) after one review; the survey (step 1); the literature notes
  (37 new notes, all dated 2026-11-14); this README.
- Result of the survey: no work was found that measures the electron injection rate or D̃_e at event conditions from
  self-consistent kinetics, or that confronts an event with a transport model built on such a measurement. Points
  where the survey changed a statement of the plan are listed in `notebook.md`, each with its fact.
- Step 3, first part: the inputs of the four cases are collected in `inputs.md` with source and kind; five further
  notes (`dresing2022role`, `lalti2022database`, `masters2011electron`, `achilleos2006orientation`, `went2011new`).
- Step 3, second part: the simulation side. `simulation-claims.md` records what published simulations of
  quasi-parallel shocks state about electrons, the limits that recur, how runs were tied to observations, what each
  planned statement can and cannot rest on, an own estimate of numerical drag, and the boundary areas (neighbouring
  methods, laboratory shocks, border regimes, domain boundaries); 29 further notes.
- Second review by `research-reviewer` (README, `inputs.md`, plan): done. Its corrections to `inputs.md` and its
  proposed thresholds are not applied; they are proposals until the owner decides.
- Missing: the inputs listed as open in `inputs.md` (instrument backgrounds as intensities; a published connection
  time at Earth; Saturn's density, β_e and speed; T_e at PSP; Mach numbers of the quasi-parallel shocks at 1 au); the
  fixed criteria and thresholds; Publication.
- Full texts not obtained (read as abstracts): `qin2026electron`, `bai2023electron`, `yang2019electron`,
  `giacalone2025source`, `wilson2025large`, `hemler2025shock`, `mitchell2021energetic`, `bessho2023electron`,
  `gupta2023regulates`, `galishnikova2026self`, `tsurutani1985acceleration`; `bebesi2019observations` only through a
  conference abstract.
- Next: the owner's decisions on the open points of `simulation-claims.md` §8 and on the criteria and thresholds;
  the reviewer's corrections to `inputs.md`, each checked in its source. Then the engine pilot (step 4), which
  shares the GPUs with the pilot of `shock-alfvenic-turbulence`.
