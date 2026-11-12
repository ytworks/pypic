# Notebook — spherical-shock-injection

Dated lab notes, newest first. Write down what was tried, what was seen, and what was decided, failed attempts included.
Never rewrite a past entry. To correct one, add a new entry that refers to it.

## 2026-11-13 — Null run of L: the numerical floor; front measurement revised to protocol Sect. 4
- **Run L-null.** `configs/explore-L-null.yaml` (no debris), 8 GPUs, commit `0114ef9`, seed 12345, 2960 steps in
  778 s (262 ms per step, 1.16 ns per particle and step). Output
  `output/research/spherical-shock-injection/L-null`, 15 GB, not in git (tracks 11 GB, reduced 4.3 GB).
  Gauss law ≤ 2.2 × 10⁻¹⁵ and div B ≤ 1.6 × 10⁻¹⁵ at the 10 steps that record them; fill ≤ 0.667.
- **Energy.** The budget grows linearly to +12.7 % at t = 740: numerical heating at 2 particles per cell, all
  of it in the electrons (electron kinetic energy +10.1 %, ion +0.2 %).
- **N1 (counts).** No electron with γ ≥ 1.5 in any of the 41 dumps, whatever its origin; largest γ 1.374; at
  most 15 per dump above γ = 1.25. With a floor of 0 the count criterion of a signal in L is the fixed one:
  100 or more.
- **N2 (temperature).** T_e in the corner column goes 0.01000 → 0.01095 at t = 740, nearly linearly (0.01004,
  0.01029, 0.01054, 0.01075 at t = 148, 296, 444, 592). The meridional map is uniform at t = 740: median
  0.0110–0.0111 in every ring, 5–95 % range 0.0104–0.0117. Above E = 0.1 the corner spectrum holds 1.37 times
  the Maxwell–Jüttner count at the fitted T (214 against 156 in weight); almost nothing above 0.2, none above
  0.5. T_i stays 0.0099. The sound speed from the measured temperatures goes 0.0455 → 0.0466.
- **N3 (scattering).** 191 335 tracked electrons, 167 026 of them in the box to the end. Pitch-angle rate ν
  at lags 4 / 16 / 64, by speed: v < 0.1: 2.6e-2 / 1.6e-2 / 8.4e-3; 0.1–0.2 (the bulk, 158 731 orbits):
  1.2e-2 / 7.8e-3 / 5.0e-3; 0.2–0.3: 5.0e-3 / 3.5e-3 / 2.6e-3; v > 0.3 (484 orbits): 2.5e-3 / 1.8e-3 / 1.3e-3.
  ν falls with the lag, so the motion is not a pure diffusion in pitch angle. |γ(end) − γ(0)|: median 0.010,
  99 % 0.058, largest 0.143; 63 orbits change by 0.1 or more.
  - Consequence, before L is seen: for the bulk electrons ν τ exceeds 0.1 once a layer takes 20 or more to
    cross (ν τ = 0.16–0.32 at τ = 20), and the dwell of 66 gives 0.33–0.55. For v > 0.3 it is 0.09–0.17 at 66.
    Unless L's measured layer is thinner than that, B2 will carry the note "loss cone filled numerically" and
    no mechanism that needs scattering can be assigned from these runs.
- **N5 (birth radius and ρ_P).** The birth cell decoded from the id lies within 0.25 d_e of the position at
  t = 0 for every tracked electron: the decoding is right. |ρ_P − ρ_birth| (median / 99 % / 99.9 % / max):
  0.26 / 0.99 / 1.26 / 1.92 d_e at t = 0 (the gyroradius), 0.76 / 2.90 / 3.71 / 5.29 at t = 740. The drift of
  ρ_P itself grows as √t: 99.9 % within 1.07 (t = 74), 2.20 (t = 296), 3.51 d_e (t = 740). A Δρ of one
  d_i = 4 d_e in L is therefore above the 99.9 % level of the null at every time.
- **Front measurement, revised to the protocol's Sect. 4** (`analyze.py front`, `shockkin.py`). Run on L-null
  as its null test: 0 of 738 (time, θ) bins show a front in both planes; nothing valid, nothing usable.
  23 of 23 self-tests pass; the new ones are a clean spherical front (flags all true, speed error against its
  own σ_V 1.3), and three fronts that must fail (planes 15 % apart, a ramp 20 d_e wide, debris 5 d_e behind).
- **Choices made while writing the code, where the protocol leaves the reading open.** None changes a number
  of the protocol.
  1. Ion-flow check (Sect. 4): from `flux:i` / `n:i`, not from `ekin:i`. The energy density cannot separate
     the bulk flow from the thermal motion. The criterion (30 %) is unchanged.
  2. Agreement of V_HT: bins in which every compared value is superluminal (V_HT ≥ 1) count as agreeing.
     A ratio of two infinite speeds carries no information.
  3. "The four quadrants agree within 15 %": the spread among the four single-quadrant fits and the folded
     fit, per plane.
  4. R is the mean of the two planes (`slice_y`, `slice_x`) where both show a front, and missing otherwise.
  5. Outputs in a row for criterion d: ceil(2/Ω_ci/Δt − 0.1) = 6 / 11 / 15 for L / M / H, as the protocol says.
  6. The mid-ramp B surface is searched only up to 3 d_i ahead of the density front, and is missing when it
     lies more than 3 d_i behind it.
  7. The comparison of V_HT between thresholds uses the fits, not the differences.
  8. E_th and v_f use the measured upstream temperature at each output (N2), since T_e rises by 9.5 %.
- **Units confirmed on the null** (needed by the checks above): the ion spectrum of `spec_up` is in m_e c²
  (T_i = 0.0099); `n:i` is in n₀ (noise 0.36 per cell at 2 per cell); `flux:i` is n v in c (noise 0.0096).
- Next: run L from the commit that holds this entry.

## 2026-11-13 — Review R1 of the protocol; cadence check; protocol revised before the first run
- **Review R1** (`reviews/2026-11-13-protocol-r1.md`, the research reviewer, on Sects. 3–8 at `ec6470e`). Its
  verdict: the structure is right and the algebra of S0 correct, but three things would let a false "hypothesis
  obtained" through: S0 with one V_HT per reflection on a curved front, two published classes missing from the
  known side (leakage, trapping), and no control of the particle number. All eight required changes are taken;
  two with a deviation, named below. The revised `exploratory-protocol.md` lists them in its Sect. 11.
- **Cadence check (R1, change 5).** `configs/check-L-cadence.yaml`, 8 GPUs, commit `2d309c8`, 600 steps in
  175 s; output `output/research/spherical-shock-injection/check-L-cadence` (4.9 GB, not in git).
  - Work of the recorded field against the energy change, per orbit with |Δγ| ≥ 0.1: read every step the
    median error is 0.4 % (all 210 orbits within 10 %); every 2 steps 0.6 % (98 % within 10 %); every 3 steps
    2.5 % (88 %); every 4 steps 5.4 % (73 %).
  - Decision: tracks every 2 steps.
- **Deviations from R1, change 1 (track stride).**
  - 947 in L and M, not 463: at every 2 steps 463 is 110 GB per run, and 291 GB are free. R1 names 947 as its
    fallback. With it the nominal prediction puts about 1 / 70 / 420 reflections on tracks (b = 2 / 3 / 4).
  - 4735 = 5 × 947 in the null run of L and in H, not the stride of L: the null holds no reflection and H
    supports no hypothesis. Its tracked electrons are every fifth one of L's.
- **Not done, against R1, change 6:** S0′ (guiding-centre orbits in the fitted front). The range of V_HT over
  the dwell is used, which R1 gives as the minimum. A residual in a band with ε > 0.3 stays "candidate".
- **Correction to the entry "16-step check of L" below:** a run does not take 45–50 minutes. A step costs
  290 ms (measured over 600 steps); the 0.90–1.06 s of the 16-step check was the cost of writing every output
  every 4 steps. A run is about 15 minutes plus the denser tracks; the null run of L will measure it.
- **Seen in the check, to keep in mind for the census:**
  - `hot_e` (γ ≥ 1.25) is dominated by debris electrons: 3.3 × 10⁶ at t = 148 against 3.5 × 10⁴ ambient ones.
    The dump levels off near 105 MB, so a run holds about 5 GB of it;
  - ambient electrons with γ ≥ 1.5 number 300–800 before any shock can be valid (t < 112). They belong to the
    start of the driver, not to a shock;
  - the work along B is about 1 % of the gain of the 210 orbits: the early gain is across the field.
- **Code state.** `analyze.py` has `front` (density front, valid domain with the numbers of the first version)
  and `cadence`; 18 of 18 self-tests pass. The R1 changes to the front measurement, and the commands for the
  floor, the census, the tracks and the records, are not written yet. They read outputs that the configs
  already write, so the null run of L starts now and the code follows while it runs.
- **Still limiting the verdict to "candidate":** the two unread predecessors (entry below), S0′, and the pair of
  runs at 2 and 8 particles per cell (made only if L shows a residual or a K1 signal).

## 2026-11-13 — Curved-shock predecessors read before the first run
- **Done:** 27 notes added to `research/literature/` (full text: 5; skimmed: 9; abstract only: 13).
  `literature.md` has a new part, "Curved fronts"; the protocol's Sect. 10 says what it means for S0, K1, K2.
- **Not obtained in full:** Krauss-Varban & Burgess 1991 (JGR 96, 143) and Giacalone 2017 (ApJ 848, 123), the
  two that matter most. The publishers' hosts refused the requests and no open copy was found. Leroy &
  Mangeney 1984 and Savoini & Lembège 1999 were not obtained at all. ADS refused requests, so the
  forward-citation search used OpenAlex only.
- **Found, against what the survey said** ("no kinetic simulation of a curved front was found"):
  - two full-PIC studies of electrons at a curved front exist, in 2-D at bow-shock parameters
    (`savoini2001twodimensional`, `savoini2010origin`; m_i/m_e = 42, M_A = 3);
  - the second shows electrons that enter near θ_Bn = 88° and leave near 48°. This is signature K1 of the
    protocol, so K1 as stated is not new;
  - the sweep of the foot-point θ_Bn in time, the source just behind the tangent line, and repeated encounters
    in closed geometries are all published.
- **Not found in what could be reached:** electrons in a 3-D spherical, slowing front in full PIC; the energy
  of an escaping electron against the one-reflection energy at its exit angle; a path across latitudes by
  drift in the ramp.
- **Consequence:** with the two key works unread, criterion C5 cannot be met; the best verdict of this
  exercise is "candidate". The owner can lift this limit with a library copy of the two papers.
- **Corrections to the survey's reference list** are kept in `literature.md` ("Vandas 1995" is two papers;
  "Vandas 2001" is two papers; Schmidt & Cairns 2014 is marked retracted in OpenAlex).
- The reading was done by a reader agent; five notes are from full texts. I checked the format of the notes,
  not each statement against the papers.

## 2026-11-13 — 16-step check of L; protocol and baseline written before the first run
- **Check run (run 0).** `configs/check-L.yaml`, 8 GPUs, commit `3e97c6c`, 16 steps in 18.8 s. All nine
  diagnostics wrote files; Gauss law 4 × 10⁻¹⁵, div B ≤ 1.3 × 10⁻¹⁶. Output:
  `output/research/spherical-shock-injection/check-L` (2.9 GB, not in git).
- **Correction to the entry of 2026-11-12 and to the plan:** a run takes about 45–50 minutes (0.90–1.06 s per
  step, measured), not 15–20. The whole run set is about 3–5 GPU-hours.
- **Changed against the plan:** the check uses `gamma_min: 1.05`, not 1.02. At 1.02 the dump would hold a
  quarter of the ambient electrons.
- **Found in the files:**
  - the loaded electrons follow a Maxwell–Jüttner distribution (T = 0.01), not a Gaussian in u: tail fraction
    at γ ≥ 1.05 measured 0.019628, Maxwell–Jüttner 0.019623, Gaussian 0.0166. The baseline uses the first;
  - `id // 2` is the birth cell of an ambient electron (checked on every electron of the dump); injected
    electrons have ids ≥ 20N; debris electrons are the species `de`, with small, varying weights;
  - reading a dump needs `import hdf5plugin` first (zstd).
- **Written:** `shockkin.py` (reflection map, foot-point formulas for an axisymmetric front, count model,
  back-tracing), `analyze.py` (`selftest`: 10 of 10 pass; `predict`), `exploratory-protocol.md`.
- **Prediction before the runs** (nominal sphere at 0.87 V_d, `data/predict-nominal.json`): one-step reflection
  gives γ ≥ 1.5 almost only at foot-point θ_Bn = 50–70°, nothing at γ ≥ 3, and a count that goes from
  5 × 10² (b = 2) to 2 × 10⁵ (b = 4) macro-electrons in the box at the end.
  - Consequence: counts cannot carry a claim. The per-electron test (does the electron trace back to a thermal
    electron reflected once?) is the primary one.
- **Not done yet:** the commands for the front, the census and the tracks; the review R1 of the protocol; no
  production run.

## 2026-11-13 — Review of the exploratory conditions (issue #48): record corrected, run set changed
- Asked by the owner: run the exploratory conditions, analyse them, and test whether a novel and certain
  hypothesis about the acceleration mechanism can be obtained; take the findings of issue #48 into the trials.
- The review is kept in `reviews/2026-11-13-exploratory-conditions.md`. Answers, point by point:
  - **A (record).** Corrected in `exploratory-conditions.md` and its Japanese version, Sects. 2.2 and 2.3:
    - each remnant value names its speed and electron temperature;
    - M_A/√β_e against the threshold of `amano2010critical` is √2 times the flow speed over the electron thermal
      speed. A remnant has about 18 at 5000 km/s and 10⁴ K and 3.2 at 3000 km/s and 10 eV; the runs have 3.5. It
      is no longer listed as matched;
    - "θ_Bn runs from 0° to 90°" now describes the initial geometry only;
    - only β of run L is listed as matched.
  - **B.1–B.5 (design).** Accepted. They go into a protocol that is committed before the first run
    (`exploratory-protocol.md`): the reflection baseline, the definition of a shock and its valid domain, the
    numerical floor from a null run, the split of energetic electrons by origin, the foot-point θ_Bn, the
    measured luminal angle, the closure test of the work, one energy threshold.
  - **C (decisions; the owner approved the plan of 2026-11-13).**
    - The record is corrected.
    - The null run of L is added as a further run; H is not replaced. The order becomes: a 16-step check, the
      null of L, L, M, then H as a contrast case only.
    - The baseline and the shock definition are written before the first run.
- The physical conditions of L, M and H are unchanged.
- Sect. 6 (checks 2 and 6) and Sect. 8 (order of the runs) of the conditions record are replaced by the protocol;
  the record keeps them as designed.

## 2026-11-12 — Three exploratory conditions designed (not run)
- Asked by the owner: two or three conditions as close to a supernova remnant as possible, for trials that
  build hypotheses: a large box, few particles, a short run. Record them in markdown; do not run yet.
- The owner chose: only what `examples/blast_wave` can express; vary the field strength only; 768³ cells with
  2 particles per cell.
- Done:
  - read the example's planner, the engine's limits and the output system;
  - designed L, M, H (nominal M_A = 3.5, 6.0, 8.7; m_i/m_e = 16; V_d = 0.25 c; v_te = 0.1 c; R_B ≈ 148 d_e);
  - wrote `configs/explore-{L,M,H}.yaml`; the planner accepts all three (`--plan-only`, CPU, plan for 8 GPUs);
  - an independent check recomputed the numbers from the planner's formulas and reviewed the output block;
  - wrote `exploratory-conditions.md` and its Japanese version.
- Found by the check, and now in the record:
  - the front runs at about 0.87 V_d and slows, so H is expected at or below the surfing threshold;
  - a track stride that shares a factor with 768 would line the tracked electrons up in planes (4733 is used);
  - the energy selection of the particle dump has no test in the repository.
- Corrected against the plan: the run time is about 15–20 minutes per run, not 7–10. The first figure counted
  the particles only, and the grid costs more here.
- Not done: no run; the electron analysis (`analyze.py`) and the run-recording tools are not written.

## 2026-11-12 — Citation check of the survey: 16 corrections
- An independent check compared about 290 located statements of Sects. 1–8 of the survey with the paper texts.
  About 280 were confirmed. Three were wrong and 13 imprecise or mis-located; all 16 are corrected in the
  English and Japanese surveys, and one wording in the summaries.
- The three wrong statements:
  - "only one paper gives Mach numbers for SNR shocks": four others quote values (M_A ≈ 20–1000);
  - "K_ep is degenerate with the magnetic field": the Tycho paper says the X-ray flux fixes K_ep independently
    of the field (the note `morlino2012strong` is corrected too);
  - "the 2-D runs last 6–52 Ω_ci⁻¹": two 2-D studies of parallel shocks run longer.
- Not checked: statements that rest on `hoshino2002nonthermal`, `amano2010critical` and `amano2022theory`
  (no text at hand), and three values that are garbled or only in figures; the survey marks them.

## 2026-11-12 — One-page summary of the survey, in English and Japanese
- Asked by the owner. `survey-snr-electron-injection-summary.md` and its `.ja.md` give the whole picture on one
  page: the problem, what observations fix, the mechanisms by θ_Bn, what decides the answer in a simulation,
  what is supported, disputed and unknown, and what it means for this topic.
- Every statement is taken from the survey, with its citekeys.

## 2026-11-12 — Japanese version of the survey
- Asked by the owner. `survey-snr-electron-injection.ja.md` has the same sections as the English survey, which
  stays the original.
- Its reference lists (Sect. 10) are taken from the English file by script; only headings and remarks are
  translated.
- Checked by script, section by section: the citekeys are identical, and no number of the English version is
  missing.

## 2026-11-12 — Survey: electron injection at supernova remnant shocks
- Asked by the owner: survey the prior work on the electron injection problem at SNR shocks, in markdown, with
  the citations kept.
- Done:
  - arXiv searches by keyword (12 queries, 238 new results), then by author and title, then reference lists;
  - 54 works selected, full texts from arXiv, read by nine readers in groups (observations; surfing line;
    perpendicular SNR series; foreshock and reviews; quasi-parallel; thresholds 2021–2026; shock drift theory;
    geometry);
  - one note per work in `research/literature/`;
  - the headline numbers of 20 works checked against the paper texts;
  - the survey written: `survey-snr-electron-injection.md`.
- Main findings:
  - mechanisms are sorted by obliquity, and each is known from planar runs at one or a few angles;
  - no study covers the range of θ_Bn with one parameter set;
  - global models of remnants assume the efficiency law ε(θ_Bn);
  - SN 1006 points to quasi-parallel acceleration and the radio population to quasi-perpendicular;
  - no kinetic simulation of a curved front at SNR-like parameters was found.
- Not done:
  - the bow-shock works on curved shocks (not on arXiv) are not read, so the last finding is not certain;
  - relativistic shocks and wakefield acceleration were not covered;
  - the survey has not been reviewed by the research reviewer.
- Corrected along the way: the reference of `amano2022theory` (ApJ 927, 132).
- No decision on the research question; it stays with the owner.

## 2026-11-12 — First review (purpose stage)
- Review: `reviews/2026-11-12-purpose.md`. Verdict:
  - "every mechanism over all θ_Bn" needs rethinking (a programme, not one question);
  - "curvature-only channels" needs the controlling number R/L_s and a target system first;
  - "the switch-on angle θ_c and its scaling with m_i/m_e and β_e" is answerable.
- Fixed in the README: the Background said that known mechanisms are framed only for planar
  shocks. amano2022theory Sect. 4.5 (Eq. 22) already has a curvature limit.
- Added to `literature.md`: the curved-shock works the review found by search (not read).
- Points to keep whatever the question becomes:
  - the mass ratio moves every threshold, with different powers;
  - θ_Bn is to be taken from the measured front normal;
  - the superluminal band;
  - Buneman waves against resolution and the current filter;
  - tail statistics per θ_Bn band;
  - tracked sets fixed at the start.
- Open: the owner's decision on the research question.

## 2026-11-12 — Topic opened
- The owner set the purpose: the electron injection problem, studied with PIC of a spherical
  shock, considering every acceleration mechanism, unknown ones included.
  - Mechanisms named: surfing, drift, wakefield, first- and second-order Fermi, combinations,
    geometric effects of the curved surface, and the shape changes as θ_Bn goes from
    perpendicular to parallel.
- Branch cut from research/trunk at 6f00c11.
- README: purpose, background and a draft research question. Prior work lists the four notes
  read so far; the survey list is in `literature.md`.
- Open: the survey, a narrowed question, the idea, and the driver and its parameters.
  `examples/blast_wave` (M_A = 3) is likely below the Mach numbers that surfing needs.
