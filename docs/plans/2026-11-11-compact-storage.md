# Compact particle storage: closing R3 and R5

Status: In progress
Refs: roadmap R3, R5 (R10 stays open); PIC-D19 (new); PIC-D11 (default to be decided by the owner)

## Context

Owner request (2026-11-11): close the roadmap rows R3, R5, R6 and R8, each independent piece of
work in its own pull request. This plan is the R3 + R5 piece (branch `perf/compact-storage`, one
pull request; the stages below are its commits). The other two pieces have their own plans
(`2026-11-11-more-examples.md`, `2026-11-11-wakefield.md`) and do not touch this file.

State of the code (main at e883591; gates on 17f1b5f, FP64 and mixed, ALL PASSED):

- A particle is eight rows of float64 (64 B per slot); its position is the physical value in the
  whole box. `precision: mixed` (PIC-D11, Proposed) changes only the gather to float32.
- The push kernel is FP64-bound: 71 % of the 8-GPU step, 94 % on 1 GPU. Production runs are
  bound by memory (bow shock `full`: 79.8 of 103 GB per GPU).
- The periodic wrap is done in FP64 after the deposit (`pypic/kernels/particles.cu`, the block
  "periodic wrap on non-decomposed axes"), so the wrapped position is rounded: today's Gauss
  error is round-off (7.6e-12 against the limit 1e-10), not zero.
- The weight must stay float64 (the neutrality check of the loader compares to 1e-12; `PW == 0`
  marks a dead slot). Ids stay below 2^53.

Owner decisions (2026-11-11):

1. R3 closes with R5 and the decision on the default precision (PIC-D11). R10 (fixed-point
   deposit) stays a row of its own. The storage uses uint32 fixed-point positions, which does not
   stand in the way of R10.
2. float32 momenta in the new mode are allowed under two conditions fixed here: all explicit
   gates pass in that mode, and in the gyration test the relative drift of |u| is at most
   6e-8·√N after N steps with no secular trend. If either fails the momenta stay float64 (56 B).

## What closes the rows

| Row | Closed when |
|---|---|
| R5 | `precision: compact` (PIC-D19) is Accepted: the explicit gates pass on one commit in FP64, mixed and compact; the peak device memory of a production-size run is within 5 % of the planner's figure; the planner uses the bytes per slot of the mode. Or: rejected by the criteria of A0, with the numbers recorded |
| R3 | every remaining item of the arithmetic is adopted with gate evidence or measured against a criterion fixed before the measurement and dropped; the owner decides the default precision from the three gate records (PIC-D11 is no longer Proposed); the final breakdown of `push_deposit` is recorded; M6 is Complete and what remains (R10) is one sentence with a measured number in the roadmap |

## Design: `precision: compact`

Explicit solver only (refused with the implicit solver, as `mixed` is); the gather is the one of
`mixed`. One slot holds

| Field | Type | Bytes |
|---|---|---|
| weight | float64 | 8 |
| id | int64 | 8 |
| local linear cell | int32 | 4 |
| position within the cell, three axes | uint32, unit 2^-32 dx | 12 |
| momentum, three components | float32 | 12 |
| total | | 44 |

Fallback: float64 momenta, 56 B.

Rules that keep the deposit exact:

- Every end point that is stored is first snapped to the lattice of 2^-32 dx. The wall loop,
  `deposit_cut` and the Esirkepov deposit are computed in FP64 from these exact doubles, so the
  deposited charge is the charge of the stored positions.
- The periodic wrap and the reflection are integer operations.
- A migration record carries the global cell per axis (a rank at a PML end has another local
  size than its neighbour).

Out of scope: `compact` with the implicit solver; `compact` with the moving window (the pull
request that merges later adds the refusal and its test).

## Stages

All criteria below are fixed here, before the measurement or the run they judge. A criterion that
fails is recorded as failed and reported to the owner; a limit or a definition changes only by an
amendment written into this file with its reason (stop rule).

Full GPU gates are taken twice in this pull request: at the commit of A1 (null test of the
refactoring) and at the final head. Everything else during development is a targeted run.

### A0 — measurement (no product code changes)

Benchmark: 50 M particles, 256×64×64, on 1 GPU and on 8 GPUs, `scripts/bench_step.py --phases`,
alternating pairs. The breakdown of `push_deposit` is obtained by differences (mixed, fp64, zero
charge). Prototypes that are thrown away (periodic boundaries only, not part of the pull
request): (a) compact positions + float64 momenta, (b) float32 momentum storage + FP64 Boris,
(c) FP32 Boris arithmetic.

| Criterion | Limit |
|---|---|
| A speed item is adopted | ≥ 2 % of the step against `mixed` (the criterion of round 3) |
| (c) is preferred to (b) | ≥ 2 %; below that the Boris arithmetic stays FP64 |
| Compact storage goes on to A1 | not more than 2 % slower than `mixed` |
| Gauss residual of the prototype | ≤ 1e-10 (reported against 7.6e-12) |

Stop rule: if compact is slower than the limit, stop and report to the owner (its benefit is
memory; the limit changes only by an amendment). The pull request is then the record of the
measurement and the close-out of R3 / R5.

The timings need the GPUs to themselves (no other job on the devices that are measured).

### A1 — accessor layer

- `pypic/kernels/particles.cu`: the raw macros `PX .. PID` become accessors. The cold kernels
  (`deposit_rho`, `deposit_moment`, `gather_at_particles`, `boris_kick`, `kinetic`, body, flags)
  no longer depend on the mode.
- `pypic/particles.py`: the accessors of `Species` return FP64 physical values. Every reader of
  `.P[` moves to them (`simulation.py`: profile weights and perturbations;
  `output/diagnostics.py`; tests, examples, scripts).
- Null test of the refactoring: at this commit the FP64 and mixed gates reproduce the values of
  the records of 17f1b5f bit for bit.

### A2 — the compact mode (PIC-D19)

- `pypic/kernelsrc.py`: `COMPACT_STORE`, `GNX / GNY / GNZ`, and `OXI ..` also in compact.
  `scripts/check_kernels.py`: compact on one periodic rank; compact split, with open, reflecting
  and paired faces and a body.
- `#ifdef` bodies: `push_deposit`, `displace_deposit` (the displacement in an array of its own),
  `load_particles` (position within the cell = the upper 32 bits of the same random number),
  `inject_face`, `inject_partners`, `migrate.cu` (records of 7 words; classification and wrap in
  integers), `Species.sort` (argsort of the cell row), `compact_and_migrate`.
- `scripts/gpu_gates.py --precision compact`, record `-gpu-gates-compact.json`;
  `scripts/verify_bench_results.py`, `tests/cpu/test_evidence.py`, `docs/gpu-runbook.md` follow.
- One constant of bytes per slot per mode in `pypic/config.py`, used by
  `pypic/setup/planner.py`, `examples/bow_shock/plan.py`, `examples/blast_wave/plan.py`. The
  temporary memory of the sort is measured again.
- Tests (GPU): the integers of the loader agree between 1 rank and 2×2 ranks and lie within
  2^-32 dx of FP64; compact versions of `mgpu_migrate.py` and `mgpu_open_flow.py` (1 vs N
  ≤ 1e-10, Gauss); the gyration test with the limit of owner decision 2; compact + implicit is
  refused; particle output and track output are FP64 physical values. CPU: the planner's bytes
  per mode.
- `.claude/rules/numerics-invariants.md`: the Boris line gets "in `compact` up to the FP32
  rounding of the stored momentum" (only if 44 B is adopted). The precision text of the README
  follows.
- Reported: the drift of the energy budget of `bow_shock:gate` and `blast_wave:gate` against
  FP64.

| Criterion | Limit |
|---|---|
| All explicit gates with `--precision compact` | the limits of the FP64 gates, unchanged |
| Gyration test, relative drift of \|u\| after N steps (float32 momenta) | ≤ 6e-8·√N, no secular trend; otherwise float64 momenta (56 B) |
| 1 vs N GPUs, compact | relative ≤ 1e-10 |
| Loader, 1 rank vs 2×2 ranks | integers equal; within 2^-32 dx of FP64 |
| Peak device memory of a production-size run | within 5 % of the planner's figure |

### A3 — close-out

Decision of the owner (not taken here): the default among `fp64`, `mixed`, `compact`. The three
gate records and the production-size comparison are shown in the pull request; the last commit
applies the decision. Then: the `Status:` lines of PIC-D11 and PIC-D19, roadmap R3 and R5 to
Done, M6 Complete with the final breakdown, and the list of what was measured and dropped, with
numbers.

## Invariants

- The gate criteria are the same in every precision. `compact` changes one line of the contract
  (|u| up to the FP32 rounding of the stored momentum), by owner decision 2. FP64 and mixed keep
  their bits (null test of A1).
- Charge conservation comes from the deposit: stored end points are lattice points and the
  deposit is computed from them. No Poisson correction.
- Loading depends only on (seed, global cell, slot); the position within the cell in compact is
  a function of the same random number, so 1 vs N GPUs stays ≤ 1e-10.
- New `#define`s (`COMPACT_STORE`, `GNX ..`) go into `kernelsrc.py` with a case in
  `check_kernels.py`.

## Verification

- `ruff check .`, `python scripts/check_kernels.py`, `pytest tests/cpu` on every commit.
- Full gates (FP64, then mixed) at the commit of A1; full gates in the three precisions at the
  final head; the records and `index.json` are committed and `PROGRESS.md` follows.

## Results

### A0 — baseline (2026-11-12, e883591, no other job on the GPUs)

`scripts/bench_step.py --example two_stream --config examples/two_stream/configs/ei.yaml --set
cells_per_wavelength=64 --set t_end=4 --set overrides.diag_interval=10`: 256×64×64, 50,331,648
particles, 100 steps after 20 of warm-up, three runs each. `kernel_push` is the phase of
`--phases` (mean over the ranks; those runs synchronise the device at each phase and are 0.2 ms
(1 GPU) and 0.3–0.4 ms (8 GPUs) slower per step).

| GPUs | precision | ms / step, three runs | ns / particle / step | `kernel_push` ms |
|---|---|---|---|---|
| 1 | fp64 | 22.02, 22.01, 22.01 | 0.437 | 20.87 |
| 1 | mixed | 18.18, 18.19, 18.18 | 0.361 | 16.94 |
| 8 | fp64 | 3.61, 3.66, 3.72 | 0.573–0.592 (GPU time) | 2.76 |
| 8 | mixed | 3.42, 3.29, 3.42 | 0.523–0.544 (GPU time) | 2.31 |

The figures of the roadmap are reproduced (22.0 ms and 0.44 / 0.36 ns on 1 GPU; 3.4 ms on
8 GPUs is the `mixed` figure). The other phases of the 8-GPU `mixed` step: `migrate` 0.49–0.59,
`reduce_ghosts` 0.22–0.23, `field_push` 0.19–0.23, `record` 0.14, `sort` 0.10, not attributed
0.19 ms. With a record at every step (`diag_interval=1`, measured first by mistake) the 8-GPU
step is 4.5 ms (fp64) and 4.1–4.4 ms (mixed): the figures above need `diag_interval=10`.

Amendment (2026-11-12, before any prototype is measured) — how the 2 % criteria are judged.
Fact: the three 8-GPU runs of the whole step differ by 4 % (3.29–3.42 ms), more than the limit
of 2 %; the spread comes from `migrate`, `reduce_ghosts` and `field_push`, which a prototype
does not touch. `kernel_push` repeats within 0.1 % (2.305, 2.307, 2.308 ms) and the 1-GPU step
within 0.1 %. So a prototype is judged by the change of `kernel_push` (plus `sort`, which (a)
changes) in ms, relative to the whole `mixed` step of the same number of GPUs (18.18 ms on
1 GPU, 3.38 ms on 8 GPUs, the mean of the three runs). The limits (2 %) are unchanged; both
GPU counts are reported and the 8-GPU figure decides, since that is the production layout.

### A0 — prototypes: method (2026-11-12, before the timings)

The three prototypes are one kernel with compile-time variants, compiled behind the particle
kernels of the engine with the same constants (so it uses the same float32 gather and the same
atomics). Storage: local linear cell (int32), position within the cell (uint32, 2^-32 dx) and
momentum (float64 in (a), float32 in (b) and (c)); the weight row is the one of the engine. The
gather weights come from the stored integers (no `floor`); the end point is snapped to the
lattice before the Esirkepov deposit, which takes the cell and the two lattice points per axis;
the wrap on the axes that are not decomposed is an integer operation and the other axes flag
their migrants. Periodic boundaries only; no walls, no injection, no body.

- Timing ("shadow" run): in the real `mixed` benchmark run (same configuration, 1 and 8 GPUs),
  at every push of the engine the particle rows are converted to the compact arrays of each
  variant (not timed); then the push of the engine and the prototypes run one after the other on
  the same particle state and the same fields, each between two device synchronisations, in an
  order that rotates with the step. Only the push of the engine acts on the run; the prototypes
  write their current into scratch arrays. So every kernel sees the states of the real run on
  every rank, under the same load. The push of the engine timed this way is the `kernel_push` of
  the baseline and is compared with it. A run of the engine kernel with zero charge on a copy of
  the rows, and the prototypes compiled without the deposit, give the breakdown.
- Sort: at a sort of the engine, the sort of the compact layout (argsort of the int32 cell row,
  one `take` per row) is timed on the same state.
- End points: at one step the positions and momenta of each prototype are compared with those of
  the engine (largest difference over all particles and ranks).
- Gauss residual: a one-rank periodic run without sort in which the prototype is the push; the
  compact arrays live through the run, the current goes into the arrays of the engine, and the
  charge density of the check is the CIC density of the stored lattice points.

First check of the driver (64×16×16, 786,432 particles, 200 steps, the GPU shared with another
job, so no timing): largest Gauss residual 1.4e-13 with the engine and 4.7e-15, 4.3e-15 and
4.6e-15 with (a), (b) and (c); the electric energy of (a) follows the engine to 1e-10 relative.
With (c) the kinetic energy of the record differs by 3.6e-4 relative at step 0: the prototype
takes γ of the record from the float32 arithmetic, which cannot resolve γ − 1 of the cold ions.
A product version would compute the record from `u² / (γ + 1)` in FP64; the timings do not
depend on it (one step in ten records).
