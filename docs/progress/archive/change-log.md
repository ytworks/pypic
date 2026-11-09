# Archived Change Log

Older `PROGRESS.md` Change Log entries, newest first. Frozen: insert only at the marker
(procedure: `docs/progress/archiving.md`).

<!-- ARCHIVE-INSERT-POINT -->

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
