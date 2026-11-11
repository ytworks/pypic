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

(filled in stage by stage)
