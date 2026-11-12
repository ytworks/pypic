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

1. Commit the change (the gates record the commit). The record is flagged dirty when the code differs from the commit, including untracked files under `pypic/`, `examples/`, `scripts/` or `tests/` (`pypic.provenance`).
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

## Mixed-precision gather (PIC-D11, Proposed)

`python scripts/gpu_gates.py --gpus 8 --precision mixed` runs the explicit gates (two-stream ei/ep/eip/buneman
and every analytic example) with `precision: mixed` and records them as
`bench/results/<date>-<commit>-gpu-gates-mixed.json` (`precision: mixed` in the file and its index entry). The
GPU tests are not repeated (they cover both precisions and are part of the fp64 record); the implicit
configuration is excluded (mixed is explicit-only). Run it after the fp64 gates on the same commit.

## Compact particle storage (PIC-D19, Proposed)

`python scripts/gpu_gates.py --gpus 8 --precision compact` runs the same explicit gates with
`precision: compact` and records them as `bench/results/<date>-<commit>-gpu-gates-compact.json`
(`precision: compact` in the file and its index entry), under the same criteria. Run it after the fp64 and
the mixed gates on the same commit.

## Diagnostics

Wall time and ns/particle/step are recorded per gate but never decide a verdict. For scaling
studies run `scripts/bench_step.py --example two_stream --config examples/two_stream/configs/ei.yaml
--set cells_per_wavelength=64 --set t_end=4` on 1/2/4/8 GPUs (through `pypic.launch` for several GPUs): it
times the steps without output, and `--phases` gives the time per phase of the step.
With `--set t_end=4` (226 steps) the outputs requested as a count per run come far more often than in
a production run (phase-space histograms every 2 steps instead of 49, `diag_interval` 1 instead of 12): time
the cost of the output on the production schedule (no `t_end` override, `--steps N --with-output DIR`).

## Production runs (reported, not gated)

Long or large example configurations measure results that are not analytic, such as power-law indices, γ_max(t) and reconnection rates. Run them on a clean commit, then record the analysis:

```bash
python -m pypic.launch -n 8 examples/relativistic_shock/run.py --config examples/relativistic_shock/configs/pair_2d_long.yaml
python scripts/record_example.py output/relativistic_shock/pair_2d_long
```

`record_example.py` writes `bench/results/<date>-<commit>-example-<example>-<config>.json`. The commit, dirty flag and GPU count come from the run's `run.json`. It also adds a `type: example` entry to `index.json`, and CI validates both record types. Commit the record together with any summary in the example's README.

- An output directory that already holds a run is refused. To repeat a run in the same place, add `--overwrite`, which removes that run's output files first.
- `record_example.py` runs only the `analyze.py` of an example in this repository, chosen by the run's `meta.example`.
- Each run has its own record, named with the run id from `run.json`. Separate runs of a configuration are kept; recording a run again replaces its record.
- Do not switch the checkout or edit `examples/` or `pypic/` while `record_example.py` runs: it refuses to record when the analysis code changed during the analysis.
- Unmeasured metrics are written as `null`.

