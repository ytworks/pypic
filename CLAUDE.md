# pypic

pypic is a GPU-accelerated, multi-GPU, 3-D special-relativistic electromagnetic Particle-in-Cell
engine (CuPy + NVRTC CUDA kernels, NCCL domain decomposition), with an explicit scheme
(Buneman–Boris + Esirkepov + Yee) and an exactly energy-conserving implicit scheme. Package: `pypic/`.

## Session start

ALWAYS read `PROGRESS.md` first — it is the cross-session memory of design changes,
milestone status, gate evidence and blockers.

## Progress log rules

@.claude/rules/progress-log.md

## Numerical and parallel invariants

@.claude/rules/numerics-invariants.md

## Test policy

Apply these rules whenever a change adds, removes, moves or reclassifies behavior or tests.

1. Keep only tests that protect current behavior from a concrete regression: necessary and
   sufficient, not maximal. Do not mirror constants or restate an invariant that product code
   already enforces at construction (e.g. the CFL check, the ghost-width check).
2. Test observable results: conservation laws, agreement with analytic/linear theory,
   1-GPU vs N-GPU equivalence, end-to-end behavior of the planner and launcher.
3. `tests/cpu/` must run without CuPy, a GPU or a driver (CI). `tests/gpu/` needs a CUDA
   device and is skipped as a whole otherwise; a GPU change is verified by running it and the
   GPU gates (`scripts/gpu_gates.py`), not by CI.
4. A verdict depends only on the code under test: never on wall-clock time, speed ratios or
   runner load. Use deterministic quantities (fixed seeds, exact identities, counts); use
   timeouts only to detect hangs; record timing as a diagnostic.
5. Statistical checks use fixed seeds and explicit sigma bounds (see the Maxwell–Jüttner test).
6. A regression test for a bug must fail on the code before the fix; check that once.
7. When a feature is deleted, delete its tests and helpers with it; file moves are not coverage.

## Where things live

- Engine: `pypic/` (kernels in `pypic/kernels/*.cu`, compile-time constants in `pypic/kernelsrc.py`)
- Set-up layer (pure Python, CPU-testable): `pypic/config.py`, `pypic/setup/`
- Experiments: `examples/<name>/` (configs, run, analyze)
- Design decisions (D-IDs, `Status:` lines): `docs/design/*.md`
- Roadmap: `docs/roadmap.md`; implementation plans: `docs/plans/YYYY-MM-DD-<topic>.md`
- GPU gate procedure: `docs/gpu-runbook.md`; recorded evidence: `bench/results/` (`index.json`)
- Archived progress history: `docs/progress/archive/`

## Dev commands

- CPU (any machine): `pip install -e ".[dev,kernelcheck]"`, then `ruff check .`,
  `python scripts/check_kernels.py`, `pytest tests/cpu` (what CI runs).
- GPU: `pip install -e ".[gpu,multigpu,dev,kernelcheck]"`, then `pytest tests/gpu` and
  `python scripts/gpu_gates.py --gpus 8` (see `docs/gpu-runbook.md`).
