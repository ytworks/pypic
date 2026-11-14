# What decides whether a quasi-parallel shock accelerates electrons to relativistic energy

Branch: `research/qpar-electron-conditions`. Started: 2026-11-14. Notebook: `notebook.md`. Survey: `literature.md`.
Inputs: `inputs.md`.
Sections follow `.claude/rules/research/stance.md`. Mark a section as open rather than guessing.

Theme and question were fixed by the owner on 2026-11-14. The plan was approved by the owner on the same day and
reviewed once by `research-reviewer` before approval. A second review (README, `inputs.md`, plan) followed step 3;
the owner approved its recommendations on 2026-11-14, and the criteria and thresholds in Approach are fixed since
then. Citekeys marked (S) have their note on
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
  - Parker Solar Probe (PSP), 2023-03-13, 0.24 au: nearly parallel (θ_Bn ~8°); an electron power law (δ ~2.7) from
    50 keV to 2.5 MeV and similar intensities in two viewing directions up to ~6.7 MeV, near the top of the
    instrument range (`jebaraj2024direct`). M_A is 7-7.5 in the later analysis (`jebaraj2026wave`); other works give
    9.1 (`jebaraj2024acceleration` (A)) and ~4 for the same crossing (`giacalone2025source`;
    `chiappetta2026evolution`, from a density its authors call a lower limit).
  - Earth: electrons above 100 keV appear with foreshock transients (`liu2019relativistic`, `raptis2025multi`,
    `raptis2025revealing` (S)). Outside them three example THEMIS distributions show "only noise > 12 keV"; the same
    work finds enhancements at ≥ 30 keV in 10 of 30 foreshock disturbances, 5 of them SLAMS
    (`wilson2016relativistic`). A
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
| Earth, no transient | 5-15 (highest counts of crossings, all geometries); ~8 as representative | 400-600 km/s | no (three example distributions: noise above 12 keV) | not a statistical upper limit; bursts ≥ 30 keV at SLAMS (two readings, `inputs.md`); no published connection time found |
| PSP 2023-03-13 | 7-7.5 (other estimates ~4, 4.68 from a lower-limit density, 9.1) | ~2400 km/s (other estimate 1350) | yes (power law to 2.5 MeV, seen to ~6.7 MeV) | shock ~4 h old; flank crossing |
| Saturn 2007-02-03 | ~100 (solar wind model) | ~450 km/s (typical value, quoted for another crossing; not measured) | yes (≥ ~1 MeV) | 3 of 508 crossings; seeds to ~100 keV upstream; density and shock motion not measured |
| 1 au quasi-parallel | not published for the subset; M_A 3.3 ± 0.3 and M_ms 1.3-2.9 in samples of all geometries | ~500 km/s (spacecraft frame, all geometries) | 40 keV increase at ~7 % of 15 | MeV cases at 1 au are oblique |

The outcomes do not follow the order of M_A.

Electron sonic Mach number, M_se = u_1/√(kT_e/m_e) (own values, `inputs.md`): PSP 0.90-1.28 (T_e 40-20 eV), Earth
0.27-0.41, Saturn 0.34 (T_e 10 eV, its upper bound) to ~1 (T_e = T_i). In every case a large share of the thermal
electrons already moves faster than the shock.

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
  background (`masters2016suprathermal`, `kouloumvakos2026wave` for the saturation at PSP). The measured δ ~2.7 at
  PSP corresponds in DSA to a compression of ~2.8, not 4 (own; the authors: "DSA predicts δ=2").

What the runs can test: H1-H3 and part of H6 directly; H4 with seed electrons as a test species; H8 through the
predicted intensity. H5 and H7 only in the form "is the shock alone enough"; if not, they stay outside this topic.
Nine candidates face four cases of different kinds, so a single answer is not assured: the paper states which
candidates the kinetic measurement excludes and which remain.

Open at the review of step 3: for protons at the PSP event, confinement is limited by wave growth, a factor of two
below the age limit and below the size limit (`jebaraj2026diffusive`, a preprint). Whether this is a separate
candidate for electrons is not decided. What is fixed: null test N8 cannot carry D̃_e from the age of a run to the
age of the PSP shock (Approach), so every transferred value is stated for the run's age.

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
   electron sonic Mach number is matched as well (convention M_se = u_1/√(kT_e/m_e), as in `gupta2024electron` (S)).
   Only v_sh/c is not. Where the event has no measurement (T_e at PSP; density, β_e and the speed at Saturn) the run
   takes a stated value from `inputs.md`. T_e at PSP is 20 ± 5 eV as adopted in one work and "≲ 40 eV, unresolved"
   in another, which gives M_se 1.28-0.90; the single value of P0 is not yet fixed and is fixed before P0 starts.
2. **Measure the dependence on v_sh/c, then transfer.** D̃_e(x, p̃) and the injection rate are measured at three
   values of v_sh/c; the exponent and its error carry the result to the real v_sh/c through a transport model.
3. **The age of the short cases is within reach.** The connection times at Earth and Saturn can be run. Our
   estimates are ω_ci t 76-305 for Earth and 46-112 for Saturn; both are size over speed, a crossing-time proxy and
   not a connection time (no published value found). At Saturn the spacecraft saw the quasi-parallel connection for
   ≥ 60-65 min before the crossing (ω_ci t ≥ 34-37), a lower bound (`inputs.md`). What a
   run covers directly is the turbulence at that age, not the momentum the electrons reach: at equal ω_ci t and D̃ the
   formula of `gupta2025speed` (S) keeps p̃·u/v_pt, and simulated electrons have a smaller u/v_sh than real ones.
   Every case is therefore confronted through the transport model.
4. **PSP offers a measured precursor.** The proton κ_eff and the wave amplitudes are known (`jebaraj2026wave`); the
   precursor of the run can be compared with them, in the same window and band only. The scales differ: the ~10 %
   and 5-6 % are values of the 58-minute precursor in 0.05-1 Hz, while a box of 10³ d_i is the last ~17 s before
   the ramp (own), where `jebaraj2024acceleration` (A) reports δB/B ≥ 0.5.
5. Said plainly: most of the difficulty is computing time. The idea is in 1-4: bring the runs into reach and carry
   what is out of reach with a stated error.

## Related work and novelty
- **Search done (2026-11-14):** forward citations of 15 seed works (337 citing works), a sweep of 10 journals for
  2016-2026 with 8 queries (84 further works), screening by title, 50 abstracts, and the text of about 22 works in
  part or whole. Record and the list of failed sources: `notebook.md`. Limits: titles only for most of the sweep;
  four full texts could not be obtained and were read as abstracts.
- **Second search (2026-11-14), on the simulation side:** 26 keyword queries on OpenAlex for the methods next to
  PIC, for runs confronted with in-situ events, for laboratory shocks and for the border regimes; the abstracts of
  about 90 works screened, 29 further notes (9 read in the body, 20 as abstracts). It is not exhaustive. What
  published simulations state, with which limits, how they were tied to observations, and what the planned runs can
  and cannot state: `simulation-claims.md`.
- **Novelty claimed, in two parts:**
  1. The injection rate and D̃_e(x, p̃) of electrons, measured in self-consistent kinetic simulations at the
     conditions of observed events, with the dependence on v_sh/c. The range that is certainly new is
     m_e/m_i ≲ p̃ ≲ 1; above it an electron has the rigidity of an ion and D is largely known
     (`caprioli2014simulationsb`). The measured proton κ_eff at PSP applies at p̃ ≈ 6-9, outside that range. The
     step is narrower than "event conditions" suggests for PSP: run B2_II of `gupta2024electron` (S) (M_A 5,
     M_se 0.9, β 0.06, mass ratio 1836, θ_Bn 30°) is close, and the step from it is M_A 5 → 7.5, θ_Bn 30° → 8°,
     β_i 0.03 → 0.2.
  2. A transport model fed with these measurements, confronted with each event as an intensity or an upper limit at
     a named energy, with errors, stating which of H1-H8 it excludes. The observed "no" is an upper limit, not an
     absence: at Saturn "electron injection may have been taking place at all these crossings"
     (`masters2016suprathermal`); at Earth bursts ≥ 30 keV occur at SLAMS (`wilson2016relativistic`); near 1 au
     electron channels can be contaminated by protons.
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
- **Reviewer's verdict after step 3 (second review):** part 1 "needs evidence; narrower than written"; part 2
  "needs rethinking in wording, not in aim". Both are applied in the two parts above. The topic keeps one question
  provided H5 and H7 stay outside it.
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
Criteria, thresholds and null tests below are fixed since 2026-11-14: they are the proposals of the second review,
approved by the owner. They hold for every production run. A change needs a new fact, named first in `notebook.md`.
Values marked "(third review, 2026-11-15)" were fixed on that date, before V1; their facts are in that notebook entry.
Where the review required an item to be settled but gave no value, the item is marked "not yet fixed" with the
facts, and is fixed before the first run that needs it.

- **Set-up:** quasi-1-D (long in x, a few cells across), explicit scheme, a piston against a wall (downstream frame),
  inflow through the open x+ face (PIC-D13), Δx = d_e/10, ~1000 particles per cell column over both species,
  precision `compact` or `mixed`, restart checkpoints (PIC-D25). Outputs under `/mnt/nvme/pypic`. V1 alone runs at
  the particle number of `gupta2024electron` (S): 200 per cell per species over 5 cells, 2000 per column.
- **Points of the set-up raised by the second search** (`simulation-claims.md` §8):
  - Changed on the review's recommendation: the definition of the injection rate (measured quantity 2); the
    particle number of V1; N1 at each v_sh/c.
  - Kept as planned, because the review gave no recommendation and no new fact requires a change. The limit of each
    is named with the results:
    - the frame (downstream): published runs differ by a factor ≲ 2 in efficiency between frames; N1 measures the
      numerical heating at each speed;
    - ~1000 particles per column in the runs other than V1, the same at each speed of the v_sh/c scan: numerical
      drag (`simulation-claims.md` §5); N2 doubles the number for P0;
    - H4 with passive seeds: the test is whether seeds in the shock's own fields reach higher energy than thermal
      electrons; the waves that seeds would drive are not included;
    - the run list: the border regimes of `simulation-claims.md` §6.3 outside it are named as limits.
- **Limits of 1-D:** no oblique whistlers; B_x cannot change, which overstates the superluminal cut-off
  (`gupta2024electron` (S)). Run X measures the difference to 2-D; the results name the limit.
- **Measured quantities:**
  1. Momentum spectra of both species upstream and downstream in time; p_max(t); compression.
  2. Injection rate, first definition: the non-thermal fraction η_nt,e of `gupta2024electron` (S) (its Eqs. 3 and
     6), the share of downstream electrons with p ≥ p_0,e = 5 p_th,e.
     - The fact behind the change from the plan's definition (the share passing p_inj, the momentum at which an
       electron outruns the shock, after `gupta2025speed` (S)): at the events' M_se thermal electrons already pass
       that momentum, so the share does not measure injection.
     - Second definition (third review, 2026-11-15): φ_e ≡ 4π p_ref³ f_2(p_ref) / n_2 at the fixed momentum
       p_ref = 6 √(m_e/m_i) m_i v_pt, f averaged over ±0.1 dex. Valid only where its bins hold ≥ 100 counts and
       ≥ 5 times the N1 floor. No published φ_e exists: it compares our runs with each other.
     - Third number, reported beside the two: the fixed-threshold fraction η_e(p ≥ 3 √(m_e/m_i) m_i v_pt) of
       `gupta2025speed` (S).
     - Decision rule: a difference between two runs counts only if η_nt,e and φ_e agree in sign and each exceeds
       2 σ_pt.
     - Downstream region (third review, 2026-11-15): D_fix = [x_sh − 15 M_A d_i, x_sh − 2.5 M_A d_i], used only
       when the wall is ≥ 10 M_A d_i behind its far edge (from t = 100 ω_ci⁻¹ on); the comparison value is the
       last block before t = 150 ω_ci⁻¹, and η(t) from t = 100 is shown. p_th,e = (momentum of the peak of p⁴f in
       D_fix) / √2 in every run.
     - Sensitivity, fixed in advance: the same numbers in [x_sh − 15, x_sh] M_A d_i and in the two halves of
       D_fix; a spread above 30 % is stated as the error of the definition.
     - Limit: D_fix is not the region of `gupta2024electron` (S). Its "1000 d_i" does not exist before t = 200 in
       B1 and never in its M_A 5 runs; the region those runs used is not stated. Where it is stated (spectra "at
       late times", density profile at t = 275 ω_ci⁻¹ with x_sh ≈ 1400 d_i), the far edge of the region is
       ≈ 400 d_i = 20 M_A d_i from its piston (own arithmetic).
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
- **Runs** (cost estimates for 8 GPUs from the cost formula below, read with a factor 1.5-2; V1's is from the
  step-4 timings):

  | Stage | Run | Conditions | Purpose | ω_ci t | Estimate |
  |---|---|---|---|---|---|
  | 1 | V1 | Gupta B1: M_A 20, mass ratio 100, θ 30°, v_pt 0.1 c, 2000 particles per column | reproduce slope, p_max(t), D_e at t = 135; fix conventions | 150 | 6.6 d |
  | 1 | V2 | Gupta B2_II: M_A 5, mass ratio 1836, β 0.06 | reproduce the stall at low M_A | 150 | 1.0 d |
  | 2 | P0 | PSP: M_A 7.5, β_i 0.2, one T_e in 20-40 eV (not yet fixed), θ 8°, mass ratio 1836, v_sh 0.133 c, seeds | injection, D̃_e, δB/B against data | 200 | 4.1 d |
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
  at ~1.6 × 10¹⁰ particle-steps/s, for ~1000 particles per column. V1's 6.6 d is not from this formula: it is from
  the step times measured in step 4 on V1's box (2700 d_i, 2000 particles per column; `notebook.md`, 2026-11-15),
  with a line fitted to a uniform plasma used for the compressed rank, so it is an estimate until a run with a
  shock is timed. The formula gave 1.8 d for V1; the other rows have not been re-estimated from the timings.
  V1's box (third review, 2026-11-15): L_box ≥ x_sh(t_end) + k λ_max(t_end), k ≥ 3,
  λ_max = M_A t̃ / [4 (1 + 4 D_2/D_1)] d_i; that is ≥ 2360 d_i, and V1 keeps 2700 d_i (k = 3.6), 2000 particles per
  column, no seed species. k = 3 is the reviewer's choice, not checked against a source. Stage 3 may shrink or change after stage 2; each stage ends with a review. A run
  of the PSP case at M_A 9-10 is not in the list; whether it is added is decided after P0.
- **Out of reach, stated:** the real v_sh/c (PSP 0.008; Earth 0.0013-0.0020; Saturn 0.0015); the real length of the
  electron precursor; Saturn with mass ratio ≥ 400 at a subluminal speed (Saturn claims are limited to mass ratio
  100); the age (ω_ci t ~10⁵) and history of the PSP shock; an absolute intensity at Saturn as long as the density
  and T_e of the event have no source.
- **Reproduction (V1, V2), fixed.** Conventions (D_B = u r_L, ⟨Δx²⟩/t) are fixed in V1. No further stage until a
  deviation is explained.
  - V1 slope: q within ±0.3 of 4 at t = 150 ω_ci⁻¹ (last block, D_fix), fitted over the fixed window
    [0.6, 2.4] m_i v_pt (third review, 2026-11-15). ±0.3 is a factor 2 in f per decade and can fail; the slope 4
    holds only for M_A ≳ 10. It is a test against the DSA value, not against a published number: the published
    spectra of B1 are at t = 225 and 275 ω_ci⁻¹, which V1 does not reach.
  - V1 early stop: at t = 20 ω_ci⁻¹ the mean shock speed over t = 10-20 is within ±10 % of v_pt/3 and the
    compression within 3.6-4.4; otherwise the run stops and the numbers are reported.
  - V1 D_e: D_1,e in [0.3, 1.2] D_B0 at t = 135 ω_ci⁻¹ and p = m_i v_pt; D_2/D_1 in [0.05, 0.2]. That is a factor
    2 around the published 0.6 and 0.1. D_1 is the value at x_sh + 2.5 M_A d_i and D_2 the value at
    x_sh − 2.5 M_A d_i (Fig. 7 of `gupta2025speed` (S): 0.55 and 0.042 D_B0, read off the figure), each in a domain
    of 200 d_i; the other six positions of that figure are reported against its values.
    - Limit: the published time average at +2.5 M_A d_i is about 1.4 D_B0, 2.5 times the t = 135 value, and the
      published variation over t = 100-275 spans about 0.2-4 (its Fig. 8; the "1.7" written here before came from
      the text's "0.6-1"). A value outside the band at t = 135 is therefore not by itself separable from that
      variation; D_1 and D_2 are also reported over t = 100-150 against the figure's band.
    - Limit: `gupta2025speed` (S) interpolates the fields with a second-order shape, `diffusion.py` linearly.
      Not tested.
  - V2: no slope criterion, because there is no power law at M_A 5. p_max within a factor 2 of 0.8 m_i v_pt at
    100-150 ω_ci⁻¹; η_nt,e within a factor 2 of the published value. The operational definition of p_max is
    not yet fixed; it is fixed before V2.
- **Scatter σ_pt:** the scatter of one measured point. It is taken from the seed pair (P0, P0b) and from ≥ 5 time
  blocks of ≥ 20 ω_ci⁻¹. In V1 (third review, 2026-11-15): five blocks of 10 ω_ci⁻¹ over t = 100-150 with a linear
  trend removed, because D_fix is clear of the wall only in those 50 ω_ci⁻¹. Neither Gupta paper states a seed scatter; with ~10⁵ tail electrons the counting error is
  negligible, so realisation and time variation dominate.
- **Null tests, fixed:**
  - N1 numerical heating: a uniform flow without a shock for the same time gives the detection floor. It is run at
    each v_sh/c, because Δx/λ_D changes along the scan (0.75, 1.5 and 3.0 at 0.133, 0.067 and 0.033 c for
    M_se 1). A signal counts if it is ≥ 5 times that floor and has ≥ 100 counts per bin.
  - N2 statistics: twice the particles, and another seed (P0b). "Unchanged" means a difference ≤ 30 % in η_nt,e and
    in the time-averaged D̃_e, and ≤ 0.15 in the slope. 30 % keeps the N4 factor ≤ 1.8 (PSP) and ≤ 2.8 (Earth,
    Saturn).
  - N3 resolution: half Δx (short). Change ≤ 2σ_pt.
  - N4 extrapolation: exponents of D̃_e and of the injection rate in v_sh/c with errors, carried to the real v_sh/c.
    - The power law is accepted if the residual of the three points is ≤ 2σ_pt, and if the 0.067 and 0.133 c
      points agree with the weak dependence that `gupta2024electron` (S) reports.
    - The extrapolated 1σ is 2.2σ_pt at PSP (0.008 c) and 3.9σ_pt at Earth and Saturn (0.0015 c). A case is
      "decided" only if its margin exceeds twice that; otherwise it is reported as "not decided".
    - Limit: with three points a pure power law is an assumption with one degree of freedom.
  - N5 mass ratio: one set at fixed electron sonic Mach number, one at fixed β_e (`shalaby2024energy`).
  - N6 boundary: twice the box length (short). Change ≤ 2σ_pt. In V1 a box of twice the length cannot be run
    to the end; in its place (third review, 2026-11-15), at t = 100, 125 and 150: (a) the energy of the
    particles with p ≥ p_0 that left through x+ is ≤ 5 % of the non-thermal energy in the box, in the form
    below, and (b) δB_⊥ in the last 200 d_i before the absorbing layer is ≤ 2 times the N1 level.
    - N6 (a) for V1 (form fixed 2026-11-15, before V1; `notebook.md`, "V1's diagnostics"). For each species
      α = e, i separately, at the spectra outputs nearest to t = 100, 125 and 150:
      R_α(t) = E_esc,α(t) / E_nt,α(t) ≤ 5 %. Momenta and energies in the box frame.
      - p_0,i = 5 p_th,i with p_th,i from the peak of p⁴f of the ions in D_fix, as p_0,e.
      - E_nt,α(t): the kinetic energy of all particles of the species with p ≥ p_0,α in the whole box at t,
        from the (x, |u|) histogram.
      - E_esc,α(t): from dumps (x, u, weight, id) of the particles with p ≥ p_d,α in the zone from 200 d_i
        inside the absorbing layer to the face, every 25,000 steps (0.75 Ω_ci⁻¹), p_d,e = 0.1 m_i v_pt,
        p_d,i = 1.5 m_i v_pt. Slab S = 5 to 200 d_i inside the layer, length L_S.
        N⁺(t) = Σ_dumps≤t (Δt/L_S) Σ_{S, p ≥ p_0, u_x > 0} w E_kin v_x; N(t) the same without the condition on
        u_x (signed); U_far(t) = the kinetic energy with p ≥ p_0 between the inner edge of S and the face at t.
        Then N − U_far ≤ E_esc ≤ N ≤ N⁺.
      - Verdict per species and time: pass if N⁺/E_nt ≤ 5 %, or if N/E_nt ≤ 5 %; fail if
        (N − U_far)/E_nt > 5 %; otherwise "not decided". (a) passes only if both species pass. From the first
        time that is not a pass, the results are labelled "outside the N6 tolerance" (or "N6 not decided") with
        that time, the box rule (k) is revisited before P0, and a sharper record (an escape diagnostic of the
        engine, or a longer box) is decided then; the estimator is not changed after the run.
      - Reported beside: R_α with p_0 replaced by p_0/2 and by 2 p_0; the escaped number fraction; the momentum
        spectrum of the outward flux; R against time at every dump; the engine's cumulative energy of removed
        particles (signed, all species, no criterion).
      - Limits: 5 % and the stand-in itself are the reviewer's choice without a source; the estimate holds if
        nothing above p_d enters at x+ and no particle crosses p_0 inside the zone (in the box frame an ion's
        energy changes by gyration in the motional field, own estimate ≲ 17 % in p at p_0,i); the receding
        boundary of `gupta2024electron` and its length at t ≤ 150 are not stated, so R_α has no published
        counterpart; a pass does not show that the turbulence far upstream equals that of a longer box (part
        (b) and Q monitor the field only at the edge).
  - Upstream field in V1 (third review, 2026-11-15): Q = (⟨B_y⟩/B_0y) / (⟨n⟩/n_0) in the zone 5 to 200 d_i
    inside the absorbing layer, blocks of 5 Ω_ci⁻¹: |Q − 1| ≤ 0.05 up to t = 150 (5 % is 1.2° in θ_Bn and
    −1.3 % in M_A). V1 starts only after the flow-through control passes its criteria (`notebook.md`,
    2026-11-15, third review).
  - N7 method: the measurement code on synthetic fields recovers a known D within 10 % over 0.1-10 D_B; the
    log-slope of ⟨Δx²⟩ is within 1 ± 0.1 over a factor 3 in time, and intervals outside that are not used. The
    statistical error is ~√(2/N), about 1.4 % for 10⁴ electrons.
  - N8 stationarity: |d ln D̃_e/d ln t| < 0.3 over the second half of the run gives the label "stationary within
    the run"; otherwise the value is reported as "at that age" only. N8 alone licenses no transfer beyond the age
    of the run: carrying D̃ from t = 200 to 1.1 × 10⁵ within a factor 2 needs a slope ≤ 0.11, which a run cannot
    resolve.
  - Added by the second review:
    - N9 frozen against tracked: D̃_e from frozen field snapshots against D̃_e from tracked PIC electrons in the
      evolving fields, at the same x and p̃. Agreement within a factor 1.5 at p̃ ≤ 1, which is below the 1.7 time
      spread.
    - N10 test-electron domain: its length 2L doubled, because the far-upstream curves of `gupta2025speed` (S),
      appendix B, do not reach an asymptote. The review gives no threshold; the change is reported.
    - N11 mechanism off: the measurement without the wave electric field is a control. The review gives no
      threshold; the difference is reported with σ_pt.
    - N12 wall: electrons that touched the piston region are left out of the spectra and the injection rate.
      Form fixed 2026-11-15, before V1 (`notebook.md`, "V1's diagnostics"). Two exclusions, both evaluated in
      D_fix and its sensitivity regions at every spectra output of t = 100–150, from a dump with ids of the
      electrons with p ≥ p_d,e = 0.1 m_i v_pt (box frame) in x = 0–900 d_i taken at the steps of the spectra
      outputs.
      - (L) origin, the rule: an electron is left out if it was loaded at x_0 ≤ X_c = 20 M_A d_i (V1: 400 d_i),
        i.e. if its fluid element met the shock before the nominal time t_c = x_0/(M_A d_i) Ω_ci⁻¹ = 20 Ω_ci⁻¹.
        Decided from the id alone. Injected electrons are never left out. Set reported: X_c = 10, 20,
        40 M_A d_i.
      - (E) wall contact, a bound: an electron is tagged if its id appears in any dump of the zone
        x ≤ L_p = 2.5 M_A d_i (V1: 50 d_i; same threshold p_d,e) taken every 20,000 steps
        (0.60 Ω_ci⁻¹ ≤ 2 L_p/c) from t = 0. The tagged set contains every electron that was reflected by the
        wall while above p_d,e; it also contains electrons that entered the zone and turned back. Removing it
        gives a lower bound of the number without wall-reflected electrons; the number without exclusion is the
        upper bound.
      - Numbers: η_nt,e, φ_e, η_e(p ≥ 3 √(m_e/m_i) m_i v_pt) and q are given (1) with (L) at X_c = 20 M_A d_i,
        the value used in comparisons between our runs and in the V1 slope criterion; (2) without exclusion,
        the value compared with published numbers of `gupta2024electron`, which states no exclusion; (3) with
        (L) at 10 and 40 M_A d_i; (4) with (E); (5) with (L) and (E) together. p_th,e, p_0,e and the particle
        number of the region are taken from the histogram without exclusion in all five.
      - Verdict, last block before t = 150: "no wall effect" if (1), (4) and (5) each differ from (2) by
        ≤ 2 σ_pt. Otherwise "wall effect detected": the range of (1)–(5) is reported, and where it exceeds 30 %
        (0.15 in q) it is stated as the error of the definition. If (1) and (2) fall on different sides of the
        V1 slope criterion, the slope is "not decided" and counts as a deviation to be explained before the
        next stage.
      - Checks: the dump reproduces the histogram counts above p_d,e in D_fix within 1 %; a dump of every 101st
        electron id without energy threshold in x = 150–900 d_i at t = 102, 126 and 150 gives the share of
        electrons below p_d,e in D_fix with x_0 ≤ X_c (expected zero) and the change of p_th,e under (L); if
        p_th,e changes by more than 2 %, that is reported and p_0,e of (1) is recomputed with the exclusion.
      - If the measured p_0,e is below p_d,e: N12 is evaluated for φ_e, η_e(p ≥ 3 √(m_e/m_i) m_i v_pt) and q
        only, and η_nt,e is reported without exclusion with that limit.
      - Limits: (L) and (E) also remove electrons accelerated by the formed shock, so a wide range does not
        show that the wall caused it; (E) misses an electron that was below p_d,e while in the zone; X_c and
        L_p have no value in the literature (`gupta2024electron` and `gupta2025speed` state no treatment of
        their piston) and are our choices; the wall as a reflecting downstream boundary is common to V1 and the
        source.
    - N13 particle number of the reproduction: V1 at the particle number of `gupta2024electron` (S).
- **Confrontation, fixed:** all comparisons are between the transport model's intensity and the instrument
  background at a named energy.
  - H1, observed side: "yes" if the predicted intensity is ≥ 3 times the background at the named energy; "no" if
    it is ≤ 1/3; otherwise "not decided". A factor 3 is at or above the N4 factor at 30 % scatter. H1 survives if
    the prediction falls on the observed side in every case within the errors of N4 and of the inputs.
    - PSP: the spectrum of `jebaraj2024direct` up to 2.5 MeV. The result is stated "for M_A ≤ 7.5". If P0 stalls
      there are two readings: the shock alone does not reproduce the event, so H5 or H7 is needed; or M_A is
      really 9-10 and a run there grows.
    - Saturn: the LEMMS channel backgrounds; a signal ≥ 10 times the seeds below 100 keV. Until the density and T_e
      of the event have a source, S1-S3 decide orderings (H3), not an intensity. Which B enters H1b is not yet
      fixed: e V B L is 87-210 keV at the event's 0.1 nT (own), below the observed ~1 MeV, and 0.2-2 MeV for
      0.2-1.0 nT (`oka2025maximum`).
    - Earth: the "no" in both readings of `inputs.md`: A, only noise above 12 keV outside the enhancements; B,
      bursts ≥ 30 keV at SLAMS, which a run of the shock alone can form. H1a cannot be tested without a
      connection time.
    - 1 au: only shocks without proton contamination of the electron channels. Run W states a result for an M_A 3
      shock at 30°, not for every quasi-parallel shock near 1 au.
  - H2, H3: the order of injection rate and D̃_e over P1, P0, P2, S. An order counts only if adjacent runs differ
    by more than 3σ_pt.
  - H4: seeds against thermal electrons in one run.
  - H6: the energy-gain budget over all tracked electrons above p_0 (≥ 10³ of them). "DSA" if ≥ 2/3 of the gain
    comes in cycles with Δp/p within a factor 2 of the DSA prediction (`gupta2025speed` (S), its Eq. 1 and Fig. 2);
    "not DSA" if ≤ 1/3; otherwise "mixed".
  - H8: predicted compression and intensity.
  - 1-D against 2-D (run X): "differ strongly" means by more than a factor 2 or by 3σ_pt, whichever is larger.
  - δB/B of P0 against the measured precursor: in the same window and band only (idea 4).
- **Outcomes that leave the question unanswered, named in advance:**
  - P0 stalls near p ≈ m_i v_pt as B2_II does (the reviewer holds this the likeliest of the three): D̃ above is not
    measured and nothing transfers to age 10⁵. The
    statement is then limited to "a uniform upstream does not grow the turbulence within 200 ω_ci⁻¹" and to the
    comparison of δB/B with the measured precursor.
  - The N4 error is so large that every case is "not decided": the claim is limited to novelty part 1.
  - 1-D and 2-D injection rates differ strongly (as defined above): the 1-D transport results are bounds only.

## Claims and evidence
No claim has a record yet. Everything below is a hypothesis. Records so far are of null tests and pilots:
N7 passes on synthetic fields (`runs/2026-11-15-n7-synthetic.md`: nine rows within 0.1-10 D_B deviate by at
most 1.6 %, criterion 10 %; the ballistic control returns no value).

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
Snapshot of 2026-11-15.

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
- Second review by `research-reviewer` (README, `inputs.md`, plan): done. The owner approved its recommendations on
  2026-11-14. Applied: its nine corrections to `inputs.md`, each checked in its source; its thresholds and added
  null tests (Approach); the narrower wording of the novelty.
- Third review (2026-11-15) applied: the second definition of the injection rate, the downstream region, V1's
  box, slope time, early stop and upstream-field tolerance are fixed (Approach; `notebook.md`).
- Not yet fixed, to be fixed before the first run that needs them: the operational definition of p_max (before
  V2); the face rule at M_se ≲ 1 and each plasma's own flow-through control (before V2, P0, E); the single T_e of
  P0 (before P0); which B enters H1b at Saturn (before the confrontation).
- Missing: the inputs listed as open in `inputs.md` (instrument backgrounds as intensities; a published connection
  time at Earth; Saturn's density, β_e and speed; T_e at PSP; Mach numbers of the quasi-parallel shocks at 1 au);
  Publication. The review lets these wait until stage 3 or the write-up.
- Full texts not obtained (read as abstracts): `qin2026electron`, `bai2023electron`, `yang2019electron`,
  `giacalone2025source`, `wilson2025large`, `hemler2025shock`, `mitchell2021energetic`, `bessho2023electron`,
  `gupta2023regulates`, `galishnikova2026self`, `tsurutani1985acceleration`; `bebesi2019observations` only through a
  conference abstract.
- Step 4: the engine change is merged (PR #52) and taken; pilots and N7 are recorded in `runs/`. Open: the
  flow-through control to 40 / Ω_ci and its companion with velocity histograms; V1 starts only if they pass the
  criteria fixed in `notebook.md` (2026-11-15, third review).
- Blocker (2026-11-15): V1 as configured needs 2.40 times the engine's particle ids (2^53) and would be stopped
  at 62.5 / Ω_ci of 150 (`notebook.md`). `plan.py` refuses such a run. None of V1's numbers is changed; the id
  block of the entering particles is engine work on main (step 4 of the plan), to be merged and taken before V1.
- Next: the engine change; close step 4 with the control; then stage 1 (V1, V2). The GPUs are shared with the
  pilot of `shock-alfvenic-turbulence`.
