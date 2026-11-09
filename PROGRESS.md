# Progress

Cross-session memory. Rules: `.claude/rules/progress-log.md`.
Older Change Log entries: `docs/progress/archive/change-log.md`.

## Product

**Goal** (owner request 2026-11-09; roadmap `docs/roadmap.md`): a GPU-accelerated, multi-GPU,
3-D special-relativistic electromagnetic PIC engine for plasma experiments.

- Two schemes, switchable: explicit (relativistic Buneman–Boris, Esirkepov, Yee; Gauss and
  div B to round-off) and implicit (Crank–Nicolson + implicit midpoint; total energy and
  charge conserved exactly). Cubic (CIC) macro-particles.
- Boundaries per face: periodic, open, reflecting. Species: electron, ion (mass ratio, Z),
  positron; e-i / e-p / e-i-p mixtures, charge neutral by construction.
- Experiments are set up from physics (linear theory → box, grid, dt, run length); first
  example: 3-D two-stream instability.
- Hardware: 8× RTX PRO 6000 Blackwell (SM120), single node, NCCL; CPU-only CI on GitHub.

## Current Status

Snapshot date: 2026-11-09. Decisions: `docs/design/engine.md` (PIC-D1..D7, Accepted).

| Milestone | Status |
|---|---|
| M1 Explicit core (Yee, Boris, Esirkepov) | Complete |
| M2 Boundaries (periodic, CPML/open, reflecting) | Complete; open-boundary injection not implemented (R1) |
| M3 Multi-GPU (NCCL decomposition, launcher) | Complete; 1 vs 2/8 GPU equal to ~1e-15 |
| M4 Mixtures + planner + two-stream example | Complete; growth within 3.5 % of theory |
| M5 Implicit energy-conserving scheme | Complete; energy error 1.9e-10 over 679 steps |
| M6 Performance | Partial: 1.0 ns/particle/step on 1 GPU; 4.0× on 8 GPUs (R3–R5 open) |
| Harness: CPU CI + GPU gates | Introduced 2026-11-09 (this PR) |

- GPU gates: ALL PASSED on 0aac3a7, 8 GPUs (`bench/results/2026-11-09-0aac3a717-gpu-gates.json`):
  growth error +2.9 / −1.8 / −2.6 / +3.5 % (ei/ep/eip/ei_implicit), Gauss ≤ 3.5e-12, div B ≤ 1.0e-14,
  implicit energy 1.9e-10; 23 GPU tests pass.
- Known gaps: no particle injection at open faces (R1); no field/particle output or restart (R2).
- Blockers: none. GPUs are shared with kairyu vLLM containers (see `docs/gpu-runbook.md`).

## Change Log

Newest first; only the most recent entries are kept here (see the size budget in
`.claude/rules/progress-log.md`).

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
