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

### A0 — prototypes: results (2026-11-12, engine code of e883591, no other job on the GPUs)

Shadow runs of the `mixed` benchmark (256×64×64, 50,331,648 particles, 100 steps after 20 of
warm-up), three runs on 8 GPUs and three on 1 GPU; ms per step, mean over the ranks, smallest
to largest of the three runs. "of the step" is relative to the whole `mixed` step of the
baseline (3.38 ms on 8 GPUs, 18.18 ms on 1 GPU), as the amendment above fixes.

| kernel | 8 GPUs, ms | change, of the step | 1 GPU, ms | change, of the step |
|---|---|---|---|---|
| push of the engine (`mixed`) | 2.220–2.223 | — | 16.945–16.948 | — |
| (a) compact position, f64 momentum | 1.893 | −9.7 % | 14.587–14.592 | −13.0 % |
| (b) f32 momentum stored, FP64 Boris | 1.949–1.954 | −7.9 to −8.0 % | 15.155–15.167 | −9.8 to −9.9 % |
| (c) f32 momentum, FP32 Boris | 1.298–1.310 | −27.0 to −27.3 % | 10.128–10.131 | −37.5 % |
| sort of the engine (per step) | 0.096–0.097 | — | 0.439–0.440 | — |
| sort, layout of (a) | 0.032–0.034 | −1.9 % | 0.281 | −0.9 % |
| sort, layout of (b) and (c) | 0.018 | −2.3 % | 0.246–0.247 | −1.1 % |

- The push of the engine timed in the shadow run is 16.95 ms on 1 GPU (baseline `kernel_push`
  16.94) and 2.22 ms on 8 GPUs (baseline 2.31; the 0.09 ms were not traced). The changes above
  are differences inside one shadow run, on the same particle states and fields.
- (b) against (a): 0.057–0.061 ms slower on 8 GPUs (+1.7 to +1.8 % of the step) and 0.57 ms
  slower on 1 GPU (+3.1 %). Storing float32 momenta and doing the Boris arithmetic in FP64 costs
  time; it would be a memory option only.
- (c) against (b): 0.64–0.65 ms faster on 8 GPUs (−19 % of the step), 5.03 ms on 1 GPU (−28 %).
- Breakdown (medians): without the deposit the engine kernel takes 1.144 ms on 8 GPUs (8.92 ms
  on 1 GPU), (a) 1.096 (8.41) and (c) 0.472 (3.40); the rest, the deposit, is 1.078 (8.02),
  0.797 (6.18) and 0.833 (6.73) ms. So (a) gains in the deposit (no `floor`, lattice end points)
  and (c) in the gather and the Boris push. (b) without the deposit was not timed.
- End points against the engine after one push (largest over all particles and ranks):
  (a) 3.5e-10 cells (the lattice is 2.3e-10) and 8.4e-12 in u; (b) 7.3e-9 cells and 2.1e-8;
  (c) 3.4e-8 cells and 4.3e-8.
- Gauss residual with the prototype as the push (one rank, 256×64×64, 50,331,648 particles,
  400 steps, no sort): engine 5.9e-13; (a) 4.1e-15, (b) 4.2e-15, (c) 4.4e-15. div B 8e-16 in all
  four. Electric energy at the end against the engine: (a) 7e-12 relative, (b) 4e-8, (c) 2.5e-7.

Against the criteria fixed in A0:

| criterion | limit | measured (8 GPUs decides) | outcome |
|---|---|---|---|
| a speed item is adopted | ≥ 2 % of the step | compact position (a): 9.7 %, plus 1.9 % from the sort | adopted |
| (c) preferred to (b) | ≥ 2 % | 19 % | FP32 Boris arithmetic, not FP64 on stored float32 |
| compact goes on to A1 | not more than 2 % slower than `mixed` | 9.7 % to 27 % faster | goes on |
| Gauss residual of the prototype | ≤ 1e-10 | 4.4e-15 at most | passes |

What the prototypes do not show: they are periodic only (no wall loop, injection, body,
migration record) and one species pair; the product kernels of A2 are measured again with the
gates. The float32 momentum of (c) stays under the two conditions of owner decision 2 (all
explicit gates in the mode; gyration drift of |u| ≤ 6e-8·√N without a secular trend), which A2
tests; if one fails the mode keeps float64 momenta — layout (a), 56 B — and (b) is not used.
The energy record of a float32 mode takes `u² / (γ + 1)` in FP64 (see the method section).

## A1 — the accessor layer as written

- Kernels (`pypic/kernels/particles.cu`): a kernel takes the storage as `PSTORE` / `PSTORE_RO`
  (`PSTORE_AS(Pe)` for the partner species of the injection), reads a particle with
  `PX..PID`, writes it with `SET_POS`, `SET_MOM`, `SET_W`, `SET_ID`, and takes the index
  coordinate from `XI`, `YI`, `ZI`. No kernel of the file indexes the array itself, with one
  exception: the remaining displacement of `displace_deposit` (three more rows during the set-up).
- Python (`pypic/particles.py`): `Species.store` (the arguments of `PSTORE`), the readers `pos`,
  `mom`, `weight`, `ids`, `rows` (FP64 physical values), the writers `scale_weight`, `add_mom`,
  `append`, and the set-up displacement (`set_displacement`, `displacing`, `displacement_max`,
  `end_displacement`). `simulation.py`, `output/diagnostics.py` and the GPU tests use only these.
- Left on the array, on purpose: `implicit.cu` and its two launches (`compact` refuses the
  implicit scheme, like `mixed`), and `migrate.cu` with `particles.migrate`, which copy whole
  rows and get their own record in A2.
- Null test of the refactor, two parts. (1) Bits: the particle rows after one step from the same
  loaded state, on the tree of 17f1b5f and on the A1 commit, FP64 and mixed; one step, because
  the order of the atomic deposits differs between two runs of the same code (ANA-D3), so later
  steps agree only to the last bits. (2) The FP64 and mixed gates on the A1 commit, every value
  beside its value in the records of 17f1b5f.
- Result of (1) (2026-11-12, one GPU; base = tree of e883591, which has the `pypic/` of 17f1b5f,
  run twice; A1 = e321bd2, which has the `pypic/` of c44d9e6; also a7497f5). Six set-ups: periodic,
  current filter, set-up perturbation, paired injection at an open face, absorbing body with test
  particles, shell of the Gauss check. Rows of all particles sorted by id and the six fields, at
  steps 0, 1 and 20.
  - Five set-ups (all but the perturbation): the loaded state is the same bits in every run, and
    the particle rows after one step are the same bits on base, A1 and a7497f5, FP64 and mixed
    (31,554 to 196,608 particles each). The fields after one step differ between the two runs of
    the base by 4.9e-19 to 1.7e-17, and by the same size between base and A1 (4.9e-19 to 2.1e-17).
    After 20 steps: mixed, the particle rows are still the same bits in all runs; FP64, they
    differ by 8.9e-16 to 3.6e-15 between the two runs of the base and by 8.9e-16 to 3.6e-15
    between base and A1.
  - Perturbation: two runs of the base already differ at step 0 (E_x by 5.6e-17 in FP64 and
    6.9e-17 in mixed, u_x of the electrons by 2.8e-17: the field of the set-up displacement is
    made by atomic deposits and the first half kick takes it), so this set-up cannot show bits.
    After one step the rows differ by 4.4e-16 between the two runs of the base and by 2.8e-17
    (FP64) and 4.4e-16 (mixed) between base and A1.
  (1) holds where the base reproduces itself. Result of (2): pending.
- Beside the null test, without a GPU: the PTX that NVRTC makes of the six cases of
  `scripts/check_kernels.py` (target `compute_120`), on the tree of e883591 and on the A1 tree.
  Field modules: the same text in all six cases. Particle modules: the text differs; with the
  registers and labels renamed in order of appearance 80 to 106 lines per case differ, and the
  sequence of instructions without the integer address arithmetic (`add.s64`, `sub.s64`,
  `shl.b64`, `mul.lo.s64`, `mul.wide.s32`, `mul.wide.u32`, `cvt.s64.s32`, `cvt.u64.u32`,
  `cvta.to.global.u64`) is the same in all six (16,112 to 19,168 instructions). The accessors
  form the addresses of the partner stores of `inject_face` in another order; no other kind of
  instruction was added, dropped or moved. This is not the criterion: the bits of (1) are.

## A2 — the compact mode as designed

Written before the code; what the code ends up doing differently is noted here with its reason.

- Arrays of one species (44 B per slot): `W` (2, cap) float64 = weight, id; `C` (cap) int32 =
  local linear cell `IDX(i, j, k)`; `F` (3, cap) uint32 = position within the cell in 2^-32 of
  the cell; `U` (3, cap) float32 = momenta. The id stays the integer-valued double of PIC-D8
  (the same 8 bytes as the int64 of the design table; no kernel or reader changes its type).
  `PSTORE` is these four pointers; `Species.store` returns them.
- Index coordinate: `XI(p) = (double)i + (double)F * 2^-32` is exact while `|i| < 2^21` (a
  double holds 53 bits, the fraction takes 32), so the kernels that work from `XI / YI / ZI`
  (`deposit_rho`, `deposit_moment`, `gather_at_particles`, `boris_kick`, `cell_index`) need no
  compact body, and the charge density of the Gauss check is the CIC density of the stored
  lattice points.
  Correction (2026-11-12): this item first read "exact (at most 31 + 32 bits)", which a double
  cannot hold. The limit is 21 bits of cell index, and it also covers the global index
  coordinate of the moving kernels and the points one step outside the box. `Simulation`
  therefore refuses `precision: compact` when the global cells of an axis plus two ghost widths
  exceed 2^20, or when the cells of one rank (with ghosts) do not fit the int32 of `C`.
- Physical position: `PX(p) = (XI(p) + OXI) * DX` — the global index first (exact), then one
  rounding, so the value does not depend on the rank layout. Used by the readers, the body test
  and the diagnostics; never by the deposit.
- Changed while writing the code: the moving kernels (`push_deposit`, `displace_deposit`, the
  injection) work in the *global* index coordinate (global cell plus fraction, exact in double),
  not the local one first planned. Reason: the walls (0 and `GN`), the points where a path
  meets a wall, the snap and the wrap `± GN` are then the same arithmetic on every rank layout,
  which the 1-GPU vs N-GPU contract needs; in local coordinates a wall hit point would be rounded
  at a different magnitude on each layout. The local cell enters only in the stored `C` and in
  the indices of the deposit and the gather.
- `push_deposit` (in the global index coordinate, see the item above): start
  `c0 + q0 * 2^-32`; the end point is snapped before anything uses it,
  `s = round((w0 + dX) * 2^32)`; walls are the integers `0` and `GN`; the reflection
  `2 * wall - x` and the periodic wrap `x ± GN` are exact in double for lattice points; a
  returned particle of a paired face gets a snapped end point and
  `wall - 2^-32` in place of `nextafter`; `deposit_cut` clips in index coordinates. A path
  without wall and cut takes the integer shape functions of the A0 prototype; the general path
  hands the two exact end points to the Esirkepov core.
- Boris rotation in FP32 (A0: 27 % of the step); the energy record takes `u² / (γ + 1)` in FP64.
  Kept only under owner decision 2 (all explicit gates in the mode, gyration drift).
- Loader: the cell is the cell of the loop and `F` the upper 32 bits of the same 53-bit random
  number, written directly — not through the physical position, whose rounding could land in
  the neighbouring cell. Face injection and its partners snap their end point the same way.
  The set-up displacement keeps its remaining part in an array of its own (three float64), which
  migrates with the particle during the set-up only.
- Migration record: 7 words of 8 bytes (weight, id, three words of global cell and fraction,
  two words of three float32). Classification by the global cell against the integer range of
  the rank, the periodic wrap `± GN` in integers; the receiver turns the global cell into its
  local linear cell.
  Changed while writing the code (2026-11-12): the record is 8 doubles, as in `fp64` — the
  global index coordinate (exact, see above), the three momenta as the doubles of the stored
  floats, weight, id — and 3 more while the set-up displacement is applied. Reason: the buffer
  kernels `mig_classify` and `mig_pack` and the exchange code then serve both layouts unchanged
  (the limits and the period are passed in cells), only `mig_take`, `mig_fill` and `mig_append`
  get a compact body, and the `fp64` module keeps its PTX. Cost: 64 bytes per migrating particle
  in the message in place of 56; the store itself is not affected.
- New defines: `COMPACT_STORE`, `GNX / GNY / GNZ` (global cells per axis); `OXI ..` also in
  compact.

## A2 — what has run

- 2026-11-12, 08daae9, 8 GPUs, no other job on them: `pytest tests/gpu -k compact` — 8 passed
  (55 s). These are the tests of the mode written in b064425: the loader against the `fp64`
  loader, the sort and the memory pool after it, the gyration bound of owner decision 2, 1 GPU
  against the decompositions (closed box and open flow), the refusal with the implicit solver,
  the output in FP64 physical values. The margins to the bounds were not printed in this run;
  a run that prints them follows the A1 gates (one GPU job at a time).
- Bytes by precision for the planners: `SLOT_BYTES` and `SORT_BYTES` in `pypic/config.py`, used
  by `pypic/setup/planner.py`, `examples/bow_shock/plan.py` and `examples/blast_wave/plan.py`.
  Not measured yet: the temporaries of the compact sort (32 bytes per particle is an estimate
  from the widths of the arrays) and the grid term of the two example planners in `mixed` and
  `compact` (31 arrays, measured in `fp64`; the float32 arrays of the gather are not counted
  separately). Both are part of the exit measurement of R5 (peak device memory within 5 % of
  the planner's figure).
