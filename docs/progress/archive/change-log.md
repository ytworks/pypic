# Archived Change Log

Older `PROGRESS.md` Change Log entries, newest first. Frozen: insert only at the marker
(procedure: `docs/progress/archiving.md`).

<!-- ARCHIVE-INSERT-POINT -->

### 2026-11-10 — [progress] implicit_push with reciprocals instead of FP64 divisions (issue #22)
- What: the particle map of the implicit kernel computes `eps / gbar`, `1/(1+t^2)` and `DT/(g0+g1)` once each, the Picard residual norms multiply by `1/h`, and `eps = q dt / 2m`, `q/m` come from the host. Implicit two-stream 1 GPU 586 → 513.5 ms/step (−12.4 %), 8 GPUs 109.5–110.6 → 99.7–101.0. Picard iterations (145 in 30 steps) and the energy budget error (3.03e-11) are unchanged.
- Why: same reason as #21 — the kernel is FP64-pipe bound and the map runs once per sub-iteration of every particle. Rounding-level change, scheme unchanged.
- Refs: #22, #4; `pypic/kernels/implicit.cu`, `pypic/simulation.py`
### 2026-11-10 — [progress] push_deposit with reciprocals instead of FP64 divisions (issue #21)
- What: `1/gm`, `2/(1+t^2)` and `DT/g1` are one reciprocal each (9 divisions → 3), `h = q dt / 2m` comes from the host, the periodic wrap is a comparison (same value for moves below one period), `g0` only on diagnostics steps. Two-stream 1 GPU 43.0–43.9 → 38.0–38.1 ms/step, 8 GPUs 6.6–6.8 → 6.1–6.3 (mixed 6.1–6.2 → 5.5), shock 8 GPUs 12.0–12.1 → 10.9–11.0. Not adopted (no gain): compiling out the open-boundary loss reductions.
- Why: the kernel is FP64-pipe bound and 75–96 % of the explicit step; divisions are its most expensive operations. Rounding-level change, scheme unchanged.
- Refs: #21, #4; `pypic/kernels/particles.cu`; `test_boris_gyration_conserves_gamma_and_period`
### 2026-11-10 — [progress] Implicit field solve: one ghost exchange per CG matvec (issue #20)
- What: `curl_E` is also taken on the ghost layer that `curl_B` reads, so the matvec exchanges ghosts once (bitwise the same on the owned unknowns); the CG vector updates are three fused kernels (`vec_add_scaled`, `cg_update`, `cg_direction`); the unused `itmax` atomic of `implicit_push` is gone (no measurable effect). Implicit two-stream, 8 GPUs: `_cg` 58–60 → 31 ms, step 136–145 → 108–110 ms/step; 1 GPU 597 → 586. CG iterations unchanged.
- Why: on 8 GPUs the CG was 42 % of the implicit step, host-side sequencing of two exchanges and ~10 small launches per iteration, not grid work.
- Refs: #20, #4; `Simulation._matvec`, `Simulation._cg`; `pypic/kernels/fields.cu`
### 2026-11-10 — [progress] Particle migration by dedicated kernels (issue #19)
- What: `particles.migrate` takes, fills, classifies, packs and appends with five kernels (`pypic/kernels/migrate.cu`) on a particle-major buffer; one host sync fewer and no allreduce; particle order unchanged (old vs new arrays bitwise equal on 1/4/8 ranks). 8 GPUs: migrate 3.8 → 0.6 ms; two-stream 9.7–10.0 → 6.7–6.8 ms/step (mixed 9.2–9.4 → 6.2–6.5), shock 14.7–15.0 → 11.8.
- Why: migration was 23–40 % of the 8-GPU step: host time of ~100 small CuPy launches for a few thousand particles (0.3 ms on 1 GPU).
- Refs: #19, #4; `scripts/bench_step.py`; `tests/gpu/test_engine.py::test_migration_across_faces_and_corners_keeps_every_particle`
### 2026-11-10 — [design] PIC-D11 mixed-precision field gather, Proposed (issue #10)
- What: `precision: mixed` gathers E, B for the explicit push from float32 copies with rank-independent weights; positions, pusher and deposit stay FP64; default `fp64`. Push 18.4 → 14.2 ms; two-stream 1 GPU 42.7 → 38.6 ms/step; shock 8 GPUs 14.4 → 13.4. All 19 explicit gates pass in mixed mode.
- Why: `push_deposit` is FP64-pipe bound; the gather is the part that can leave FP64 without touching Gauss (deposit) or 1 vs N GPU (global weights).
- Refs: #10, #4; PIC-D11 in `docs/design/engine.md`; `test_mixed_gather_is_decomposition_independent`
### 2026-11-10 — [amendment] Esirkepov deposit on the 3-node CIC support (issue #6)
- What: the explicit deposit loops over the 3 nodes per axis that a sub-cell CIC move touches (2 prefix sums per component), and the gather shares its weights. Electron `push_deposit` 37.9 → 30.8 ms; two-stream 1 GPU 51.6 → 42.7 ms/step, 8 GPUs 10.7 → 9.2–10.2; shock 8 GPUs 16.0 → 14.6–14.9.
- Why: `push_deposit` is FP64-pipe bound (85 %) on SM120; same scheme, fewer FP64 operations.
- Refs: #6, #4; PIC-D2 review amendment in `docs/design/engine.md`; `test_esirkepov_current_moment_equals_the_displacement`
### 2026-11-10 — [progress] Ghost exchange with fused slab kernels (issue #8)
- What: ghost copy / add / zero / pack / unpack of up to 3 arrays per launch, one message per face; a self-wrapping axis is copied in place. 8 GPUs: two-stream 11.7 → 10.5–10.9 ms/step (diag every step 13.7 → 11.6–12.4), shock 17.8 → 16.0 ms/step.
- Why: the 8-GPU step was bound by host-side sequencing of ~30 small operations per ghost call, not by bandwidth; no compute/communication overlap (the field kernels take < 0.1 ms per step, nothing to hide behind).
- Refs: #8, #4
### 2026-11-10 — [progress] Batched particle migration (issue #7)
- What: `particles.migrate` moves every species at once: one message group per decomposed axis, sort-based partition, 3 host syncs per step (was ~20 for two species), one allreduce for the diagnostics record. 8-GPU two-stream: NCCL calls 16.2 → 14.2 per step; step time unchanged within noise (11.3–12.2 → 11.7 ms).
- Why: the issue's premise; the profile after it shows the 8-GPU step is bound by host-side sequencing (GPU idle ~5 of 14 ms under nsys), not by the syncs' count.
- Refs: #7, #4; `tests/gpu/test_engine.py::test_migration_across_faces_and_corners_keeps_every_particle`
### 2026-11-10 — [progress] Current filter as fused ping-pong pass kernels (issue #9)
- What: one kernel per pass filters up to 3 arrays, passes alternate with work buffers; bitwise the former result. Shock `pair_2d_long` on 8 GPUs: filter cost 2.1 → ~1.6 ms of a 17.8 → ~17.6 ms step; ~250 → ~35 launches per step.
- Why: the filter was ~250 small launches and two sweeps per pass; it is memory-bound (~1 TB/s per pass), and the halo exchange per pass (15 per step) is now most of what remains (follow-up: two passes per sweep).
- Refs: #9, #4; PIC-D10; `tests/gpu/test_engine.py::test_current_filter_kernels_match_the_array_expressions`
### 2026-11-10 — [progress] CG inner product in two stages (issue #5)
- What: `Simulation._dot` reduces over the inner axes, then sums; implicit two-stream (50.3M particles) 1017 → ~570 ms/step on 1 GPU, 151 → 144 ms/step on 8 GPUs.
- Why: a CuPy ReductionKernel to one scalar runs as one thread block (3.1 ms per call, 49 % of the implicit GPU time).
- Refs: #5, #4; `tests/gpu/test_engine.py::test_cg_inner_product_matches_the_owned_unknowns`
### 2026-11-10 — [progress] Particle row offsets are 64-bit (issue #11)
- What: the particle kernels take `cap` as 64-bit, so `row * cap` no longer overflows above ~2.7e8 particles of one species per GPU (it was an illegal address). `Species` rejects a capacity above 2**31 - 1 (the kernels' int particle index).
- Why: a 96 GB GPU holds ~1e9 particles; the int32 offsets capped a species at ~1/4 of that (R11 needs large runs). The regression test failed on the old code.
- Refs: #11, #4; `tests/gpu/test_engine.py::test_particle_rows_beyond_int32_offsets`

### 2026-11-10 — [progress] GitHub wiki articles of the examples (EN + JA)
- What: one page per example and language: conditions table, GIF animations and comparison charts with theory or the published results, all generated from the runs (`scripts/build_wiki.py run|collect|render`, templates and numbers in `docs/wiki/`). New charts in the examples (growth ratio, branch errors, convergence, jump conditions, p vs σ, reconnection field lines with X/O points); `pypic.analysis.animate` lays out with the first frame's titles and spans the spectrum axis over all frames. README numbers corrected.
- Why: owner request; the numbers come from runs (fresh gate configs at a8caa8f, recorded production runs), so the pages cannot drift from the data.
- Refs: `docs/wiki/`, `tests/cpu/test_wiki.py`, `bench/results/2026-11-10-d03ab601a-gpu-gates.json`
### 2026-11-10 — [amendment] PR #2 fifth review: rows of a failed HDF5 flush
- What: `SeriesWriter._write` keeps the rows buffered until the HDF5 flush succeeds, so a retry repeats the write and the flush.
- Why: the review found that after a failed flush the next flush() returned without flushing, and SWMR readers saw no row. The new CPU test failed on that code.
- Refs: review amendments in `docs/design/output-analysis.md`; `tests/cpu/test_output.py`
### 2026-11-10 — [amendment] PR #2 fourth review: scalar rows of a failed write
- What: scalar history rows count as written once the scalars writer has taken them (rank 0), not before the write. A failed write is retried, and rows already taken stay buffered.
- Why: the review found a regression from the third-review change: after a failed write, a retry wrote nothing. The new test failed on that code.
- Refs: review amendments in `docs/design/output-analysis.md`; `tests/gpu/test_output.py`
### 2026-11-10 — [progress] Second recorded realizations of the production runs
- What: `pair_sigma10_large` (run e879cfd6) and `pair_2d_long` (run 40c0ca85) were run again and recorded, one record per run.
  - Reconnection, two recorded runs: p = 2.10 / 2.10 and rate 0.064 / 0.074 v_A.
  - Shock, two recorded runs: β_sh −5.9 / −7.6 %, n₂/n₁ +6.9 / +8.1 %, tail p 3.91 / 3.85, γ_max ∝ t^0.36 / t^0.41.
- Why: the 2026-11-09 entry "Production runs of the acceleration examples" quoted a first realization of each that was never recorded (second review). Those values are withdrawn: p = 2.15 and rate 0.080; −4.2 %, +7.8 %, p 3.93 and t^0.42.
- Refs: `bench/results/*-run-e879cfd6.json`, `*-run-40c0ca85.json`; `examples/{reconnection,relativistic_shock}/README.md`; ANA-D3
### 2026-11-10 — [amendment] PR #2 third review: rank-0 file errors, untracked code
- What:
  - Every rank-0 file operation of the output shares its outcome with all ranks (`OutputError` on every rank).
  - Provenance has one definition (`pypic.provenance`), used by run.json, gate records and example records. Untracked code files make a commit dirty; changes outside the code paths do not.
- Why: the review reproduced a rank-0-only `PermissionError` that left the other ranks in a collective, and an untracked analysis helper that was missing from `analysis_commit`. Each fix has a test that failed on the old code.
- Refs: OUT-D2, ANA-D3 and review amendments in `docs/design/output-analysis.md`

### 2026-11-09 — [amendment] PR #2 second review: recorded results
- What:
  - `run.json` carries a run id and a start time, both fixed at start-up. Example records are one per run, with the run id in the file name and index entry, so realizations no longer overwrite each other.
  - The analysis commit is taken before the analysis, and nothing is recorded if the code changes during it.
  - Evidence is strict JSON (NaN → null; the five existing files were normalized), and the verifier enforces it.
- Why: the review reproduced 3 issues; each fix has a test that failed on the old scripts.
- Refs: ANA-D3 and review amendments in `docs/design/output-analysis.md`; `scripts/bench_records.py`
### 2026-11-09 — [amendment] PR #2 review fixes
- What:
  - The density-perturbation inverse is a bisection with a residual check reduced over ranks. Newton stalled at α → 1.
  - A test species' virtual move deposits `qdep`.
  - One run per output directory (`output.overwrite` / `--overwrite`).
  - Also fixed: empty tracks, typed `output.*` overrides, `fit_two_wave` (finite weights, coefficients at t = 0), the work split at B = 0, the drift-frame orbit theory up to t_end, and `record_example.py` (repository examples only).
- Why: the review reproduced 10 bugs (1 P1, 9 P2). Each fix has a regression test that failed on the old code.
- Refs: review amendments in `docs/design/engine.md` and `docs/design/output-analysis.md`; PR #2
### 2026-11-09 — [progress] Production runs of the acceleration examples
- What:
  - pair_3d: 3-D jump conditions, β_sh −3.7 %, n₂/n₁ +3.5 %.
  - pair_2d_long: early Fermi tail p ≈ 3.9 at t = 2000, γ_max ∝ t^0.36–0.42.
  - isothermal_hydrogen_long: E_max −12 % up to τ = 30.
  - pair_sigma10_large: p = 2.10–2.15 vs Werner+2016's 2.12; rate 0.06–0.08 v_A (two realizations: the atomic deposit order is not deterministic).
  - All recorded as typed example evidence (`scripts/record_example.py`, validated in CI).
- Why: the reported (non-analytic) results planned for PR #2.
- Fixes: the shock speed is fitted over 0.3–0.7 t_valid (the precursor reaches the plasma edge later); run.json takes the git state once at start-up.
- Refs: `bench/results/*-example-*.json`; `examples/{relativistic_shock,reconnection,plasma_expansion}/README.md`

### 2026-11-09 — [design] Analytic-solution examples and their set-up features (PR #2, phases 2-4)
- What: examples plasma_waves, landau_damping (+ Dawson, wave breaking, relativistic), filamentation, particle_orbits, plasma_expansion, relativistic_shock, reconnection, each gated against theory (ANA-D2); engine PIC-D9 (Gauss-exact perturbations, explicit half kick, test species and particle lists, E0 / exact plane waves / field profiles, density profiles) and PIC-D10 (binomial current filter, Gauss- and implicit-energy-exact); `pypic.theory` waves/streaming.
- Why: owner request (one example per analytic problem); each feature keeps the conservation contracts (GPU tests: Gauss <= 1e-12, implicit energy <= 1e-11, 1 vs 2 GPUs).
- Findings: Bernstein needs theta_e = 1e-3 (relativistic shift at 1e-2); PIC expansion ~15 % below Mora's isothermal E_max (kinetic electrons; gate set above the measured offset); Buneman needs remove_mean_current.
- Refs: PIC-D9, PIC-D10, ANA-D2; `examples/*/README.md`; plan `docs/plans/2026-11-09-output-analysis-examples.md`

### 2026-11-09 — [design] Research output, analysis package, particle ids (PR #2, phase 1)
- What: HDF5 output — openPMD dumps + reduced series, cadence from physics, collective diagnostics written by rank 0, files independent of the GPU count; particle ids (PIC-D8); `pypic.analysis` (reader, omega-k / growth / power-law fits, plots, movies, XDMF, CLI); shared example runner; Buneman config with `remove_mean_current`; `Simulation.save()` / `mode_diag` / `phase_space_*` removed.
- Why: owner request — keep research data in light files at the cadence the physics needs, with general analysis tools and one example per analytic problem (plan `docs/plans/2026-11-09-output-analysis-examples.md`).
- Fixes: implicit error flags reduced over ranks (rank-local raise hung callers); particle on the upper reflecting wall of a decomposed axis no longer dropped.
- Refs: OUT-D1..D4, ANA-D1, PIC-D8; `docs/design/output-analysis.md`; `tests/gpu/test_output.py`

### 2026-11-09 — [amendment] PR #1 second review fixes
- What: the implicit deposit folds a segment about the wall chosen from the unfolded midpoint (was the folded one, wrong when |d| > L/2 → charge not conserved); `SimConfig.from_dict` accepts `SpeciesConfig` instances so typed planner overrides work again.
- Why: two P2 findings from the second review, each reproduced by a regression test that fails on the pre-fix code.
- Refs: PR #1 review on 96ec343; PIC-D3, PIC-D7

### 2026-11-09 — [progress] CPU CI and development harness
- What: GitHub Actions CPU job (ruff, PROGRESS budget, bench evidence index, NVRTC compile check of all kernels for sm_80/sm_120, `tests/cpu`); tests split into `tests/cpu` (no CuPy) and `tests/gpu`; CuPy became the `gpu` extra; GPU gates script recording evidence in `bench/results/`; CLAUDE.md/AGENTS.md, rules, design record, roadmap, runbook.
- Why: keep the conservation guarantees verifiable across sessions and contributors; CI runners have no GPU, so GPU evidence is recorded and validated instead.
- Refs: PR #1; `.github/workflows/ci.yml`, `scripts/{check_kernels,gpu_gates,verify_bench_results,check_progress_size}.py`; first evidence `bench/results/2026-11-09-fcbc4ced5-gpu-gates.json`

### 2026-11-09 — [amendment] PR #1 review fixes
- What: CG 0/0 guard; exact Maxwell–Jüttner sampling (Gamma-mixture rejection for theta < 1); non-converged implicit steps raise and are not committed; ghost width validated against cells per rank; nested planner overrides via `SimConfig.from_dict`; launcher polls all ranks.
- Why: six review findings (2× P1, 4× P2), each reproduced by a regression test that fails on the pre-fix code.
- Refs: PR #1 review; PIC-D3, PIC-D5, PIC-D6, PIC-D7

### 2026-11-09 — [design] Initial engine (PIC-D1..D7)
- What: explicit and implicit schemes, three boundary types, NCCL multi-GPU, neutral MJ loading, two-stream planner/example; verified on 8 GPUs (growth vs theory +2.6 / −3.2 / −3.7 / +2.1 %).
- Why: owner request; implicit scheme chosen in addition to explicit to guarantee total energy.
- Refs: PR #1; `docs/design/engine.md`
