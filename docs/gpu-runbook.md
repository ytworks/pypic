# GPU runbook

How to verify a change on a GPU machine. CI is CPU-only (lint, kernel compile check, CPU
tests); everything that executes kernels runs here. Criteria: `.claude/rules/numerics-invariants.md`.

## Prerequisites

- CUDA 13 driver, `pip install -e ".[gpu,multigpu,dev,kernelcheck]"`.
- The GPUs must be free. On the development machine (8× RTX PRO 6000 Blackwell) the kairyu
  vLLM containers occupy all GPUs; stop them only with the owner's agreement and restart them
  afterwards:
  `docker stop gemma-4-31b-8gpu-gemma-{0,1,2,3}-1 gemma-4-31b-8gpu-embed-{0,1}-1`, later
  `docker start` with the same names.
- Check: `nvidia-smi --query-compute-apps=pid,process_name --format=csv` shows no other jobs.

## Procedure

1. Commit the change (the gates record the commit; a dirty tree is flagged).
2. `python scripts/gpu_gates.py --gpus 8` (about 30-40 minutes):
   - `gpu-tests`: `pytest tests/gpu`, including the 1-vs-2-GPU equivalence tests;
   - `two-stream:{ei,ep,eip,ei_implicit,buneman}` on 8 GPUs: growth rate vs theory, Gauss law,
     div B, implicit energy;
   - the analytic-solution examples (`<example>:<config>`, criteria in `docs/design/output-analysis.md`
     ANA-D2); `particle_orbits` runs on one GPU at Δt and Δt/2 (convergence).
   While iterating use `--quick`, `--configs ei --skip-tests` or `--examples landau_damping:landau`
   (not recorded).
   A hung run is killed with its whole process group (no ranks are left spinning on the GPUs).
3. On ALL PASSED the script writes `bench/results/<date>-<commit>-gpu-gates.json` and updates
   `bench/results/index.json`. Commit both (CI validates the index), and update
   `PROGRESS.md` → Current Status if the gate status changed.
4. On failure: the JSON lists `failed_checks` and a log tail. Never relax a criterion to pass;
   a criterion changes only by a design amendment recorded in `docs/design/` and `PROGRESS.md`.

## Diagnostics

Wall time and ns/particle/step are recorded per gate but never decide a verdict. For scaling
studies run `examples/two_stream/run.py` with `--set cells_per_wavelength=64 --set t_end=4`
on 1/2/4/8 GPUs.

## Production runs (reported, not gated)

Long or large example configurations measure results that are not analytic, such as power-law indices, γ_max(t) and reconnection rates. Run them on a clean commit, then record the analysis:

```bash
python -m pypic.launch -n 8 examples/relativistic_shock/run.py --config examples/relativistic_shock/configs/pair_2d_long.yaml
python scripts/record_example.py output/relativistic_shock/pair_2d_long
```

`record_example.py` writes `bench/results/<date>-<commit>-example-<example>-<config>.json`. The commit, dirty flag and GPU count come from the run's `run.json`. It also adds a `type: example` entry to `index.json`, and CI validates both record types. Commit the record together with any summary in the example's README.

