# Notebook — spherical-shock-injection

Dated lab notes, newest first. Write down what was tried, what was seen, and what was decided, failed attempts included.
Never rewrite a past entry. To correct one, add a new entry that refers to it.

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
