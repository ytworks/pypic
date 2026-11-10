# Numerical and parallel invariants

The engine's value is its guarantees. Every change to `pypic/` must keep these contracts;
the GPU gates (`scripts/gpu_gates.py`, procedure in `docs/gpu-runbook.md`) check them.

## Conservation contracts (gate criteria)

| Quantity | Explicit | Implicit |
|---|---|---|
| Gauss law `max|div E - rho|` (checked region, below) | ≤ 1e-10 | ≤ 1e-10 |
| `max|div B|` | ≤ 1e-12 | ≤ 1e-12 |
| Total energy budget (particles + fields + boundary outflow + absorbed) | monitored only | relative ≤ 1e-8 |
| Two-stream growth rate of the k_max mode vs linear theory | ≤ 10 % | ≤ 10 % |
| 1-GPU vs N-GPU history (same seed) | relative ≤ 1e-10 | relative ≤ 1e-10 |

- Checked region of the Gauss law: the physical nodes without the conducting faces. At an open
  face it also leaves out the face node and the nodes that the charge left there reaches:
  `passes + compensate` nodes with a current filter, one more with the Marder zone (PIC-D15,
  PIC-D16; `Simulation.gauss_retreat`).
- Charge conservation comes from the deposition, never from a Poisson correction: keep the
  Esirkepov / orbit-split deposition consistent with the CIC charge shape. The Marder zone at
  open faces (PIC-D16) changes E only on edges between nodes outside the checked region; no
  divergence cleaning may reach a checked node.
- Implicit energy exactness needs the gather and the deposit to use the same path weights,
  `vbar = (u0+u1)/(g0+g1)`, and the per-particle energy closure. Never commit an implicit step
  whose Picard iteration did not converge (`ImplicitStepError` leaves the state at step n).
- Boris rotation must keep `|u|` exact; changes to the pusher keep `test_boris_*` green.

## Change gate

- Any change under `pypic/` (kernels, grid, particles, simulation, parallel) is verified on a
  GPU machine before merge: `pytest tests/gpu` and `python scripts/gpu_gates.py --gpus N`
  on a clean commit. Commit the result file and index entry (`bench/results/`) and note a
  status change in `PROGRESS.md`. CI cannot run these; say so in the PR if not yet done.
- A new `#define` used by a kernel goes into `pypic/kernelsrc.py`, and `scripts/check_kernels.py`
  gets a case that exercises it.

## Multi-GPU rules (each one was a real bug)

- Every collective (`allreduce*`, `exchange`, `total_particles`, gauss/energy diagnostics,
  migration, output diagnostics and their gathers) must be reached by every rank. Never call one
  inside a rank-0-only branch such as `if self.verbose:` or `if rank == 0:` — compute it first on
  all ranks; only the file writing itself is rank 0 only. Per-rank error flags are reduced before
  any rank raises.
- Buffers handed to NCCL inside a group (including temporary contiguous copies) must stay
  referenced until `groupEnd()` returns.
- Halo exchange and migration reach only the adjacent rank: the ghost width must not exceed the
  cells per rank (enforced at initialisation).
- Particle loading depends only on (seed, global cell, slot) — never on the rank layout — so
  results are decomposition-independent.

## Engine vs examples

`pypic/` stays experiment-agnostic. An experiment (its species mix, beam set-up, planner
targets, analysis) lives in `examples/<name>/`; move a helper into `pypic/setup/` only when a
second example reuses it.
