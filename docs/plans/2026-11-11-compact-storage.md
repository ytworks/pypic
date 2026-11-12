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

Amendment (2026-11-12): the fallback applies. Float32 momenta did not meet the gyration
condition, so the slot holds three float64 momenta and takes 56 B ("A2 — result" below).

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

Outcome (2026-11-12): the gyration criterion is not met with float32 momenta (a secular trend;
"A2 — result" below), so the mode has float64 momenta and 56 B. The gyration test of the
float32 mode goes with it, the Boris gyration test of `fp64` runs in `compact` too, and
`.claude/rules/numerics-invariants.md` is not changed.

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
  Amendment (2026-11-12): with float64 momenta ("A2 — result") `compact` changes no line of
  the contract.
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

Outcome (2026-11-12): the gyration condition fails with float32 momenta, in layout (b) by a
trend measured on the product kernel and in layout (c) by a change of |u| per step whose cause
is known; the mode keeps float64 momenta, layout (a). See "A2 — result".

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
  (1) holds where the base reproduces itself.
- Result of (2) (2026-11-12, 8 GPUs, e321bd2): FP64 ALL PASSED (93 GPU tests, five two-stream
  gates, 17 example gates) and mixed ALL PASSED (four and 17); records
  `bench/results/2026-11-12-e321bd238-gpu-gates.json` and `…-gpu-gates-mixed.json`. The table
  gives per gate the number of metrics compared (timings left out) and, for A1 against the record
  of 17f1b5f, how many differ and the largest relative difference with its metric. Beside it the
  same for the record of 9bfa1a5 against that of 17f1b5f: two runs of one `pypic/` tree (the
  analysis of the blast wave changed between them, so that row of the base pair holds this
  change too).

  | Gate | Metrics | FP64: A1 | FP64: base pair | mixed: A1 | mixed: base pair |
  |---|---|---|---|---|---|
  | `gpu-tests` | 3 | 0 | 0 | — | — |
  | `two-stream:ei` | 6 | 3, 2.6e-03 `gauss_max` | 3, 1.8e-02 `divB_max` | 3, 4.1e-03 `gauss_max` | 3, 1.5e-01 `divB_max` |
  | `two-stream:ep` | 6 | 3, 1.1e-01 `divB_max` | 3, 4.1e-02 `divB_max` | 3, 1.3e-01 `divB_max` | 3, 3.4e-02 `gauss_max` |
  | `two-stream:eip` | 6 | 3, 8.3e-02 `divB_max` | 3, 1.6e-01 `divB_max` | 3, 4.6e-02 `divB_max` | 3, 7.2e-02 `divB_max` |
  | `two-stream:ei_implicit` | 6 | 3, 4.0e-02 `divB_max` | 3, 1.1e-01 `divB_max` | — | — |
  | `two-stream:buneman` | 6 | 3, 6.2e-02 `divB_max` | 3, 1.7e-02 `divB_max` | 5, 4.1e-02 `divB_max` | 5, 1.0e-01 `divB_max` |
  | `plasma_waves:unmagnetized` | 18 | 8, 2.7e-02 `divB_max` | 8, 1.5e-01 `divB_max` | 8, 3.4e-02 `divB_max` | 8, 2.8e-02 `divB_max` |
  | `plasma_waves:perpendicular` | 38 | 16, 3.0e-03 `divB_max` | 16, 3.0e-02 `divB_max` | 16, 8.9e-02 `divB_max` | 16, 1.4e-01 `divB_max` |
  | `landau_damping:landau` | 9 | 1, 2.0e-03 `gauss_max` | 1, 1.9e-03 `gauss_max` | 1, 4.3e-04 `gauss_max` | 1, 4.1e-03 `gauss_max` |
  | `landau_damping:dawson` | 4 | 1, 5.9e-02 `gauss_max` | 1, 4.4e-03 `gauss_max` | 2, 2.3e-02 `gauss_max` | 1, 1.3e-01 `gauss_max` |
  | `landau_damping:wavebreaking` | 5 | 1, 8.9e-02 `gauss_max` | 1, 3.3e-02 `gauss_max` | 1, 3.0e-02 `gauss_max` | 1, 1.7e-01 `gauss_max` |
  | `landau_damping:relativistic` | 6 | 1, 9.4e-03 `gauss_max` | 1, 1.8e-02 `gauss_max` | 2, 2.5e-02 `gauss_max` | 1, 3.5e-02 `gauss_max` |
  | `filamentation:pair` | 56 | 2, 3.6e-02 `gauss_max` | 2, 1.6e-01 `gauss_max` | 2, 1.0e-02 `gauss_max` | 2, 1.0e-01 `gauss_max` |
  | `particle_orbits:gyration` | 14 | 0 | 0 | 0 | 0 |
  | `particle_orbits:exb` | 14 | 0 | 0 | 0 | 0 |
  | `particle_orbits:hyperbolic` | 14 | 0 | 0 | 0 | 0 |
  | `particle_orbits:e_gt_b` | 14 | 0 | 0 | 0 | 0 |
  | `particle_orbits:plane_wave` | 16 | 0 | 0 | 0 | 0 |
  | `plasma_expansion:isothermal_hydrogen` | 10 | 6, 8.4e-03 `electron_temperature_final_over_initial` | 6, 3.6e-03 `velocity_profile_median_rel_error` | 5, 1.8e-02 `velocity_profile_median_rel_error` | 6, 3.8e-03 `E_max_median_rel_error` |
  | `relativistic_shock:pair_2d` | 72 | 11, 2.3e-01 `compression_rel_error` | 12, 1.1e-01 `compression_rel_error` | 12, 2.6e-01 `compression_rel_error` | 13, 2.4e-01 `compression_rel_error` |
  | `reconnection:pair_sigma10` | 13 | 7, 6.1e-01 `E_parallel_share_top1pct` | 7, 1.2e-01 `E_parallel_share_top1pct` | 7, 1.0e+00 `E_parallel_share_top1pct` | 7, 7.1e-01 `E_parallel_share_top1pct` |
  | `bow_shock:gate` | 30 | 23, 2.3e-02 `gauss_max` | 23, 4.6e-02 `gauss_max` | 23, 1.1e-01 `gauss_max` | 23, 1.0e-01 `gauss_max` |
  | `blast_wave:gate` | 93 | 31, 1.7e-01 `gauss_max` | 54, 2.0e-02 `slice_y` | 56, 8.0e-02 `divB_max` | 59, 1.1e-01 `divB_max` |

  - The same in all digits, FP64 and mixed: the five `particle_orbits` gates (72 values); in
    FP64 the counts of `gpu-tests` (93 passed). Every other gate differs in some values, and so
    does the base pair: the gates do not show bits (the order of the atomic deposits, as in (1)).
  - Most of the largest differences are in `gauss_max` and `divB_max`, residuals of rounding
    (Gauss maximum of the record 7.61e-12 in FP64 and 7.65e-12 in mixed; 17f1b5f: 7.61e-12 and
    7.68e-12). The others, as 9bfa1a5 / 17f1b5f / A1: compression of the shock 3.068 / 3.076 /
    3.060 (FP64) and 3.071 / 3.055 / 3.072 (mixed); `E_parallel_share_top1pct` of the
    reconnection 0.0044 / 0.0038 / 0.0015 (FP64) and 0.0008 / 0.0026 / −0.0001 (mixed); final
    electron temperature of the expansion 1.0393 / 1.0378 / 1.0291 (FP64) and 1.0372 / 1.0383 /
    1.0427 (mixed). The growth rates of the two-stream gates are equal in all digits in the
    three records, except `buneman` in mixed (0.050969 / 0.051226 / 0.051226).
  - The A1 figure is the smaller one in 8 of the 17 gates that differ in FP64 and in 10 of 16 in
    mixed. One base pair per precision is not a spread; the table is a comparison of sizes.
  - In mixed two metrics differ on A1 where the base pair agrees in all digits:
    `freq_rel_error_median` of `landau_damping:dawson` (0.0030547052120 against 0.0030547064732)
    and of `landau_damping:relativistic` (0.00074694638 against 0.00074694771). Repeated for
    `dawson`, mixed, 8 GPUs, twice on each tree: e883591 gives 0.0030547065849 and
    0.0030547064732, e321bd2 gives 0.0030547064621 and 0.0030547064732. The metric changes in
    these digits between runs of the base, and one run of the base has the bits of the A1 record.
    `relativistic`, mixed, three times on each tree: e883591 gives 0.00074694771 once and
    0.00074694638 twice, e321bd2 gives 0.00074694771 three times; the base has both values.
  (2) shows no difference between A1 and the base that two runs of the base do not show.
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
  Changed (2026-11-12): `U` is float64 and the slot takes 56 B ("A2 — result").
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
  Changed (2026-11-12): not kept. The rotation and the energy records are those of `fp64`
  ("A2 — result").
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
  Changed (2026-11-12): with float64 momenta the record carries the momenta as they are
  stored; the message takes 64 bytes per migrating particle and the store 56 per slot.
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
  Measured (2026-11-12, 403eeb6, one GPU, no other job on it; `two_stream/ei` with
  `cells_per_wavelength=64`, 33,554,432 electrons and 16,777,216 ions): the memory pool emptied,
  then one sort; the growth of the pool and of the used device memory up to the end of the
  sort, per particle of the sorted species, twice per species with equal results.

  | precision | pool, bytes per particle | device, bytes per particle |
  |---|---|---|
  | fp64 | 40.22 | 40.25 |
  | mixed | 40.22 | 40.25 |
  | compact | 24.19 | 24.25 |

  `SORT_BYTES` of `compact` is 24 (the estimate was 32). The 40 of the bow-shock gate is
  reproduced. The grid term of the example planners is still open.

## A2 — result (2026-11-12): float32 momenta fail the gyration condition; momenta stay float64

New fact: with the product kernel of f322fe2 (float32 stored momenta, FP64 Boris arithmetic:
layout (b) of A0) the mean of |u| over particles of one |u| that gyrate in a uniform field
changes with the number of steps — by up to 26 standard errors after 2.3e6 steps, and at a
constant rate per step in the group that shows it most clearly. That is a secular trend, so the
second condition of owner decision 2 is not met. By the rule fixed with the A0 results ("if one
fails the mode keeps float64 momenta — layout (a), 56 B — and (b) is not used") `compact` stores
float64 momenta, 56 bytes per slot. The bound 6e-8·√N and the two conditions are not changed.
The owner can amend the decision; what is known for that is at the end of this section.

### Method

`precision: compact`, one GPU, 8³ cells with `dx` 0.5, uniform `B0` along z, one species with
q/m = −1 and weight 1e-20, no sort, no diagnostics. 1,024 particles per value of |u|; the phase
of each is uniform in [0, 2π) (`numpy.random.default_rng(seed)`), the momentum lies in the plane
across B and is set as the float32 of u·cos, u·sin. After n steps the relative change of |u| of
every particle is taken against its own start. In the tables "mean/se" is the mean of a group
over its standard error, "mean/n" the mean per step, and "max" the largest |change| among all
particles of the run over 6e-8·√n. The script and its outputs are not in the repository.

### Float32 momenta, product kernel of f322fe2

B = 0.5, dt = 0.05 (562 steps per turn at |u| = 2), |u| = 1, 2, 5, 11.3, seed 1:

| steps | mean/se | mean/n | max |
|---|---|---|---|
| 2,248 | +0.11 +2.71 −7.51 −4.04 | +1.2e-12 +3.1e-11 −1.6e-10 −6.4e-11 | 1.570 |
| 35,968 | −1.45 +2.41 −4.86 −3.31 | −3.9e-12 +6.3e-12 −1.8e-11 −1.3e-11 | 1.446 |
| 143,872 | −5.23 +3.63 −7.06 −0.12 | −7.1e-12 +4.7e-12 −1.4e-11 −2.3e-13 | 1.326 |
| 575,488 | −9.54 −3.05 −19.50 +3.60 | −6.8e-12 −2.2e-12 −1.4e-11 +3.5e-12 | 1.043 |
| 2,301,952 | −19.58 +21.46 −15.88 +16.41 | −6.5e-12 +5.1e-12 −6.3e-12 +4.3e-12 | 0.914 |

The same with seed 2:

| steps | mean/se | mean/n | max |
|---|---|---|---|
| 2,248 | +1.04 +4.04 −7.34 −2.35 | +1.1e-11 +4.6e-11 −1.6e-10 −3.6e-11 | 1.593 |
| 35,968 | −1.12 +2.65 −1.40 −1.90 | −3.1e-12 +7.4e-12 −5.4e-12 −7.0e-12 | 1.280 |
| 143,872 | −4.36 +5.58 −7.85 +0.09 | −5.9e-12 +7.2e-12 −1.5e-11 +1.8e-13 | 1.297 |
| 575,488 | −8.44 −2.00 −21.15 +2.98 | −6.1e-12 −1.4e-12 −1.5e-11 +2.8e-12 | 1.069 |
| 2,301,952 | −19.95 +25.09 −16.03 +15.99 | −6.9e-12 +5.6e-12 −6.8e-12 +4.2e-12 | 0.932 |

Values of |u| without short binary fractions, |u| = 0.37, 0.83, 1.7, 2.9, 4.6, 5.3, 6.7, 11.3
(B = 0.5, dt = 0.05, seed 1):

| steps | mean/se | max |
|---|---|---|
| 2,248 | −0.20 −0.98 −3.39 +0.86 +0.20 −0.35 +1.30 −3.35 | 1.675 |
| 35,968 | +0.68 −0.38 −0.65 +2.09 +3.06 −0.50 +4.18 −2.96 | 1.398 |
| 143,872 | −1.16 −4.93 −0.46 −1.24 +3.17 −2.97 +3.72 +3.27 | 1.605 |
| 575,488 | −9.81 −13.51 +4.88 +8.40 +7.31 −8.67 +9.70 +6.74 | 1.367 |
| 2,301,952 | −20.55 −16.26 +9.23 +16.71 +0.74 −26.00 +23.01 +17.21 | 0.859 |

B = 2, dt = 0.2 (35 steps per turn at |u| = 2), |u| = 0.5, 1, 2, 3, 5, 8, seed 1:

| steps | mean/se | max |
|---|---|---|
| 560 | +1.30 +0.16 −0.59 +0.18 +0.16 −1.52 | 1.323 |
| 2,240 | +1.51 +1.14 +2.80 +1.17 −0.22 −1.79 | 1.460 |
| 35,840 | +2.67 +2.23 +3.93 +2.09 −0.38 +3.53 | 1.279 |
| 143,360 | −8.74 +5.15 +2.33 +4.16 +4.92 +4.71 | 1.147 |
| 573,440 | −22.57 +13.60 −1.17 +1.01 +45.69 +7.09 | 0.744 |
| 2,293,760 | −24.52 +13.65 −1.07 +1.08 +49.97 +6.94 | 0.530 |

What the tables show:

- A trend. In the |u| = 1 group the mean per step is the same from 143,872 to 2,301,952 steps
  (−7.1e-12, −6.8e-12, −6.5e-12 with seed 1; −5.9e-12, −6.1e-12, −6.9e-12 with seed 2): the
  mean grows in proportion to n and ends 19.6 and 20.0 standard errors from zero. With the
  eight other values of |u|, seven groups end 9.2 to 26.0 standard errors from zero. With
  B = 2, dt = 0.2 four of six groups are 7 to 46 standard errors from zero at 573,440 steps.
- It needs a long run and a group to be seen: at 2,248 steps (2,240 with B = 2) 16 of the 22
  groups are within three standard errors (|u| = 5 is not: −7.5 and −7.3), and the test of this
  mode written in b064425 (one particle, 16 turns, |change| ≤ 6e-8·√n at the end of each turn)
  passed.
- Not explained, and not examined further: the |u| = 2 group changes sign twice (positive up
  to 143,872 steps, negative at 575,488, +21 and +25 standard errors at 2,301,952; both seeds);
  the |u| = 4.6 group is +7.3 standard errors from zero at 575,488 steps and +0.7 at 2,301,952;
  with B = 2, dt = 0.2 the mean of the |u| = 5 group stops growing after 573,440 steps (mean
  per step +2.1e-11, then +5.4e-12) while in |u| = 0.5 the mean and its standard error both
  grow about fourfold.
- The largest single change is 0.53 to 1.68 of 6e-8·√n at the marks: above 1 at 16 of the 17
  marks before 1e6 steps, below 1 at the four after. Read for every particle the bound is
  therefore exceeded by some particle, apart from any trend: in a CPU emulation of the same
  arithmetic with 16,384 particles of random |u| (three settings of B and dt, up to 35,968
  steps) the spread of the change is 2.15e-8·√n to 2.20e-8·√n (the bound is about 2.8 of
  these), 0.61 to 0.85 % of the particles lie above the bound at each mark, and the mean stays
  within 2.25 standard errors. Whether the bound of the decision is meant for every particle
  or for the spread was not fixed; the outcome here rests on the trend.

### Controls

- Float64 storage, same arithmetic (`precision: mixed`, same script, seed 1, |u| = 1, 2, 5,
  11.3): mean/se −0.26 +1.80 −0.49 −0.06 at 2,248 steps, +0.51 −0.11 −0.39 −0.95 at 35,968,
  −0.03 +1.00 −0.21 −0.99 at 143,872 and +0.79 +1.62 −0.32 +0.28 at 575,488; the mean per step
  is 4e-20 or less in size, max 0.000. The trend comes with the float32 storage.
- `compact` with float64 momenta (the tree of this commit, one GPU): the digits of the `mixed`
  control at 2,248, 35,968 and 143,872 steps (the run ends there), max 0.000.
- A NumPy emulation of the float32 storage with the FP64 rotation gives the printed digits of
  the product kernel at both seeds up to 143,872 steps. Later the |u| = 2 group differs at
  575,488 steps (−3.48 against −3.05 with seed 1, −2.36 against −2.00 with seed 2) and the
  |u| = 2 and 11.3 groups differ by 0.16 at most at 2,301,952 (+21.35 against +21.46 and
  +16.33 against +16.41 with seed 1); |u| = 1 and 5 agree throughout. The reason for the late
  differences was not examined.
- Where the change is made (emulation, |u| = 5, 4,096 particles, 2,248 steps; mean change
  −4.221e-7, −18.05 standard errors): the rotation in FP64 contributes −7.8e-17 and the
  rounding of its result to float32 the rest. Sorted by the angle of the momentum in 5° bins,
  the eight bins at 0°, 85°, 90°, 175°, 180°, 265°, 270° and 355° hold −334.6e-9 of the
  −422.1e-9 per particle; the other 64 bins lie between −8.3e-9 and +5.2e-9. The mechanism is
  not established.

### FP32 Boris arithmetic, layout (c): not used

Independent of the storage, the rotation in float32 changes |u| by a fixed amount per step.
With the rotation factors t and s = 2t / (1 + t²) rounded to float32, `u' = u + (u + u×t)×s`
multiplies u² by (1 − s·t)² + s², which is 1 only for the exact s. A float32 emulation of the
rotation alone (B = 2, dt = 0.2, 512 particles per group, |u| = 0.5, 1, 2, 3) against this
factor summed along the run, as the mean relative change of |u| per step:

| steps | measured | standard error | predicted |
|---|---|---|---|
| 35 | +4.301e-9 −1.091e-9 −8.643e-10 +4.285e-10 | 1.6e-10 1.4e-10 1.3e-10 1.7e-10 | +4.196e-9 −1.116e-9 −1.075e-9 +3.994e-10 |
| 140 | +3.818e-9 −1.130e-9 −9.447e-10 +4.067e-10 | 1.0e-10 6.3e-11 6.3e-11 7.6e-11 | +3.798e-9 −1.098e-9 −1.065e-9 +3.965e-10 |
| 560 | +2.925e-9 −8.795e-10 −9.364e-10 +3.947e-10 | 7.0e-11 2.6e-11 2.9e-11 3.8e-11 | +2.976e-9 −8.948e-10 −1.000e-9 +3.888e-10 |
| 2,240 | +1.458e-9 −4.290e-10 −7.822e-10 +3.388e-10 | 3.3e-11 8.2e-12 1.2e-11 1.9e-11 | +1.453e-9 −4.432e-10 −7.921e-10 +3.494e-10 |

Measured and predicted agree within 2.2 standard errors. In the emulation of layout (c) the
group means are +44.36 −52.27 −65.79 +17.51 −8.35 −2.01 standard errors from zero after 2,240
steps (|u| = 0.5, 1, 2, 3, 5, 8), where layout (b) has 2.16 at most (+26.76 for |u| = 0.5
after 35 steps already). So (c), the 27 % of A0, fails the same condition within a few turns,
and its cause is known.

### What changed in the code

- `CompactSpecies.U` is float64; the accessors `S_UX ..` and `S_SET_MOM` read and write
  doubles; the migration record carries the momenta as they are. `ROUND_U` and `GAMMA_M1`
  (rounding to the stored float, and γ − 1 from it) are removed; the energy records are those
  of `fp64` again. `SLOT_BYTES["compact"]` is 56.
- The PTX of the `fp64`, `mixed` and implicit builds is the one of A1; only the two compact
  particle modules change.
- Tests: `test_boris_gyration_conserves_gamma_and_period` runs in `fp64` and `compact`; its
  `compact` case fails on the float32 code and passes now (checked once). The gyration-bound
  test of the float32 mode is deleted with the feature. The loader and output tests compare
  the momenta with `fp64` exactly. On 4 GPUs `pytest tests/gpu -k "compact or boris"`: 9 passed.
- `.claude/rules/numerics-invariants.md` is not changed: with float64 momenta `compact` keeps
  the Boris line of the contract as it stands.

### For the owner's decision

The mode costs 12 bytes per slot more than planned (56 against 44; `fp64` has 64). Measured in
A0 with prototypes: 56 B with FP64 momenta is 9.7 % faster than the `mixed` step; float32
storage with FP64 arithmetic was 7.9 to 8.0 % (no gain in speed over float64 storage, only the
12 bytes). Size of the trend where it was measured: 1e-11 per step or less at 2,301,952 steps
with B = 0.5 (|u| = 1: 1.5e-5 in all), against a bound of 9.1e-5 there. Allowing float32
momenta needs an amendment of owner decision 2 by the owner; nothing here assumes it.

## A2 — speed of the product kernels

### First measurement (2026-11-12, 5edc7f9, no other job on the GPUs): slower than `mixed`

The benchmark of the A0 baseline (`two_stream/ei`, `cells_per_wavelength=64`, 256×64×64,
50,331,648 particles, 100 steps after 20 of warm-up, `diag_interval=10`), three runs per
precision; `kernel_push` and `sort` are the phases of the `--phases` runs (mean over the ranks),
smallest to largest of the three.

| GPUs | precision | ms / step, three runs | `kernel_push` ms | `sort` ms |
|---|---|---|---|---|
| 1 | fp64 | 23.16, 22.01, 22.02 | 20.866–20.869 | 0.442 |
| 1 | mixed | 18.75, 18.20, 18.20 | 16.944–16.947 | 0.442 |
| 1 | compact | 20.19, 19.95, 19.95 | 18.859–18.862 | 0.283–0.284 |
| 8 | fp64 | 3.71, 3.69, 3.65 | 2.762–2.778 | 0.102–0.107 |
| 8 | mixed | 3.32, 3.29, 3.31 | 2.306–2.310 | 0.103–0.111 |
| 8 | compact | 3.41, 3.43, 3.44 | 2.477–2.483 | 0.071–0.080 |

Judged as A0 fixed it (change of `kernel_push` plus `sort`, medians, relative to the `mixed`
step of the baseline, 18.18 ms on 1 GPU and 3.38 ms on 8 GPUs): +1.915 − 0.158 = +1.757 ms on
1 GPU, +9.7 % of the step; +0.172 − 0.028 = +0.144 ms on 8 GPUs, +4.3 % of the step. The
product kernel of 5edc7f9 is slower than `mixed` by more than the 2 % of the A0 criterion
("compact goes on: not more than 2 % slower than `mixed`"). Prototype (a) of A0, the same
layout, was 9.7 % (8 GPUs) and 13.0 % (1 GPU) faster in the push. The difference between the
prototype and the product kernel is therefore 0.5 ms (8 GPUs) and 4.3 ms (1 GPU) of
`kernel_push`; where it comes from is not established yet.

What differs in the code (read, not measured): the prototype decodes the linear cell once,
takes the gather weights from the stored integers, snaps the end point with one integer
rounding, deposits with integer shape functions, wraps in integers and stores cell and
fraction directly. The product `push_deposit` of 5edc7f9 goes through the global index
coordinate in double on every particle: three decodes of the cell through the accessors, Yee
weights with `floor`, the snap with `rint`, the Esirkepov core with `floor` on both ends, the
wrap in double, and `floor` again in the store (the wall loop and the cut are not in a
periodic build). The design above ("a path without wall and cut takes the integer shape
functions of the A0 prototype") was not carried out in the code of 5edc7f9.

Next: give the path without wall and cut the integer route of the prototype and measure again.
The criterion and its limit are not changed.

### Integer route, step 1: the Esirkepov core takes the shapes

`esirkepov` is split into `esirkepov_shapes` (the deposit from the three `AxisShape` and the
"this axis moved" flags) and a wrapper that builds the shapes from the two end points, so that
the integer route can hand over shapes made from the stored integers without a second copy of
the core. This changes the particle PTX of every build, and the statement under "What changed
in the code" above ("the PTX of the `fp64`, `mixed` and implicit builds is the one of A1")
holds up to 5edc7f9 only. Compared for the eight cases of `scripts/check_kernels.py` (823aa67
against the split): the field modules are identical; each particle module has the same number
of lines; in the two periodic cases the lines are the same set in another order; in the six
others two register names are exchanged in 6 or 16 lines and the rest is the same set in
another order. No instruction is added, removed or changed in kind.

### Integer route, step 2: the step that meets nothing (b830153)

In the compact build `push_deposit` starts from the stored integers and leaves them only where
a step meets something.

- The cell is decoded once (`load_cell`). The gather weights come from the cell and the
  fraction q (`yee_weights_cell`): q · 2^-32 from the node and (q + 2^31 mod 2^32) · 2^-32
  from the half-integer node, both exact in a double and rounded once to float.
- The end point is put on the lattice by one integer rounding per axis,
  `s1 = rn((w0 + u dt / (γ dx)) · 2^32)` with w0 the stored fraction: the change of cell is
  `s1 >> 32` and the new fraction the lower 32 bits. It depends on the stored fraction and on
  u only, not on the cell or the rank layout. Both paths below take this end point (the
  general path had `WSNAP` of the end point in the global index coordinate before).
- A step is "inner" when its end cell lies inside the box on every axis that is not periodic
  and its end point lies outside the body. Such a step meets no wall, no cut and no removal.
  Its shapes are built from the integers (`axis_shape_cells`) and deposited by
  `esirkepov_shapes`; the cell is wrapped on a periodic axis that is not decomposed; the
  particle is listed for migration when its end cell leaves the rank on a decomposed axis;
  cell and fraction are stored directly (`store_cell`).
- Every other step takes the general code as it was (wall loop, `deposit_cut`, removal at
  open faces, the body), in the global index coordinate.
- No new `#define`: `pypic/kernelsrc.py` and the cases of `scripts/check_kernels.py` are
  unchanged.

Checked (2026-11-12):

- Against a build with the route switched off (`inner = false`: every step takes the general
  path), one GPU, `ppc = 8`, electrons in two beams and ions, three cases: periodic, 16×8×8;
  reflecting on all faces, 16×8×8; open in x, reflecting in y, periodic in z with an absorbing
  sphere, 20³. The stored cells, fractions, momenta and weights and the particle counts have
  the same bits after 1 step and after 31 steps in all three. J differs by 3.9e-16 at most
  (largest |J| 0.45), E by 1.7e-16 at most, the history by 3.6e-14 at most (the budget, 20.3).
  Whether two runs of one build differ as much was not measured.
- `python scripts/check_kernels.py` passes, `pytest tests/cpu`: 72 passed; on 4 GPUs
  `pytest tests/gpu -k "compact or boris or reduced_precision or esirkepov"`: 12 passed.
- PTX of the eight cases, b830153 against cee6abe: the field modules are identical. The six
  particle modules without `COMPACT_STORE` have the same number of lines and the same set of
  lines, 12 to 42 of them at another position (the one change outside `COMPACT_STORE` is that
  `axis_shape` calls the new `axis_shape_cells`). The two compact particle modules grow from
  18,952 to 18,988 and from 23,343 to 23,756 lines.

### Second measurement (2026-11-12, b830153, no other job on the GPUs): faster than `mixed`

The benchmark, the script and the order of the runs of the first measurement.

| GPUs | precision | ms / step, three runs | `kernel_push` ms | `sort` ms |
|---|---|---|---|---|
| 1 | fp64 | 23.16, 22.01, 22.01 | 20.868–20.869 | 0.442 |
| 1 | mixed | 19.08, 18.20, 18.20 | 16.945–16.946 | 0.441–0.442 |
| 1 | compact | 15.40, 15.32, 15.33 | 14.238–14.246 | 0.283–0.284 |
| 8 | fp64 | 3.65, 3.61, 3.66 | 2.763–2.784 | 0.106–0.109 |
| 8 | mixed | 3.33, 3.25, 3.32 | 2.288–2.292 | 0.102–0.108 |
| 8 | compact | 2.86, 2.97, 2.99 | 1.950–1.958 | 0.077–0.078 |

Judged as before (change of `kernel_push` plus `sort`, medians, relative to the `mixed` step of
the baseline): −2.707 − 0.159 = −2.866 ms on 1 GPU, −15.8 % of the step; −0.338 − 0.028 =
−0.366 ms on 8 GPUs, −10.8 % of the step. The A0 criterion ("compact goes on: not more than
2 % slower than `mixed`") is met by the product kernel of b830153. The push alone is −14.9 %
(1 GPU) and −10.0 % (8 GPUs) of the step; prototype (a) of A0 had −13.0 % and −9.7 %.

Whole step without `--phases`, medians: 15.33 ms against 18.20 (`mixed`) and 22.01 (`fp64`) on
1 GPU, which is 0.305, 0.362 and 0.437 ns per particle and step; 2.97 ms against 3.32 and 3.65
on 8 GPUs. That is 15.8 % and 10.4 % below `mixed`, 30.4 % and 18.6 % below `fp64`.

Limits of this result:

- One benchmark, periodic and without a body: every step is inner. The speed of a build with
  walls, open faces or a body, where the test for the inner step is paid and some steps take
  the general path, is not measured.
- The `kernel_push` of `mixed` on 8 GPUs is 0.017 ms lower than in the first measurement (the
  ranges 2.288–2.292 and 2.306–2.310 do not overlap; on 1 GPU and in `fp64` they do). The
  `mixed` particle PTX of the two commits has the same lines in another order. Not examined.
- The correction of one figure of the first table: the third `mixed` run on 8 GPUs was
  3.31 ms (3.3146), not 3.32.

## A2 — the open-boundary gate configurations in the three precisions

2026-11-12, f34ac5a, 8 GPUs (rank grid (1, 2, 4)), no other job on them; one run per row, in
the order of the rows. Both configurations have six open faces with injection, an absorbing
sphere and the current filter (4 passes, compensated): the test for the inner step is paid for
every particle, and the particles that meet a face or the body take the general path. Neither
has a reflecting wall or a periodic axis.

- `bow_shock/gate`: 240×320×320 cells, 391,940,798 particles at the start, 5760 steps,
  capacity factor 2.
- `blast_wave/gate`: 384³ cells, 943,088,238 particles at the start, 1760 steps, capacity
  factor 1.5.

Device memory: `memory.used` of `nvidia-smi`, sampled once per second on every GPU during the
run. "Peak" is the largest sample of the GPU that had the largest one, "lowest GPU" the largest
sample of the GPU that had the smallest one. The planner's figure is `memory_per_gpu_GB` of the
example planner for 8 GPUs (1 GB = 1e9 bytes). The wall time is `Simulation.wall_time`: the
step loop with its diagnostics and output.

| run | precision | planner GB | peak GB | lowest GPU GB | peak − planner | wall s | ms / step | wall against fp64 |
|---|---|---|---|---|---|---|---|---|
| bow shock | fp64 | 11.00 | 10.87 | 10.36 | −1.2 % | 264.0 | 45.83 | — |
| bow shock | mixed | 11.00 | 10.78 | 10.12 | −2.0 % | 235.8 | 40.94 | −10.7 % |
| bow shock | compact | 9.42 | 9.36 | 8.79 | −0.7 % | 219.0 | 38.02 | −17.0 % |
| blast wave | fp64 | 21.45 | 22.21 | 17.67 | +3.6 % | 164.2 | 93.30 | — |
| blast wave | mixed | 21.45 | 22.15 | 19.74 | +3.2 % | 149.4 | 84.89 | −9.0 % |
| blast wave | compact | 18.52 | 18.78 | 17.40 | +1.4 % | 137.9 | 78.35 | −16.0 % |

Budget change: `(budget_end − budget_start) / budget_start` of the scalar history (the `dE/E0`
of the run log). Every run passed the criteria of its gate (`analyze`).

| run | precision | budget change | minus fp64 | Gauss max | div B max | particles at the end |
|---|---|---|---|---|---|---|
| bow shock | fp64 | +2.59608e-2 | — | 3.8e-14 | 7.0e-14 | 433,705,825 |
| bow shock | mixed | +2.59793e-2 | +1.9e-5 | 4.5e-14 | 7.4e-14 | 433,718,684 |
| bow shock | compact | +2.59088e-2 | −5.2e-5 | 3.6e-14 | 6.8e-14 | 433,766,996 |
| blast wave | fp64 | +1.94570e-2 | — | 7.1e-15 | 9.1e-15 | 936,673,613 |
| blast wave | mixed | +1.94545e-2 | −2.5e-6 | 7.3e-15 | 9.8e-15 | 936,670,065 |
| blast wave | compact | +1.94541e-2 | −2.9e-6 | 7.6e-15 | 8.4e-15 | 936,671,523 |

What the runs show:

- Memory: the peak is within 5 % of the planner's figure in all six runs. In the blast wave it
  is above the figure the planner calls "up to" in every precision, by 3.6 %, 3.2 % and 1.4 %.
  `compact` lowers the peak by 1.51 GB per GPU in the bow shock (13.9 %; the planner's
  difference is 1.57 GB) and by 3.43 GB in the blast wave (15.5 %; planner 2.93 GB).
- The GPUs of one run differ: in the blast wave the largest samples of the eight GPUs lie
  between 17.67 and 22.21 GB in `fp64` and between 17.40 and 18.78 GB in `compact`. Which
  allocation makes the difference between the GPUs was not examined.
- Time: the wall time of `compact` is 7.1 % (bow shock) and 7.7 % (blast wave) below `mixed`.
  One run each; the variation between runs of the same build was not measured, and the wall
  time holds the diagnostics and the output, which do not depend on the precision of the push.
- Energy budget: its change over the run differs from `fp64` by at most 5.2e-5 of the budget
  (`compact`, bow shock); the change itself is 2.6e-2 and 1.9e-2 (explicit scheme with open
  faces: monitored, no limit). The `energy_budget_change` of the blast-wave analysis, which is
  relative to the initial kinetic energy of the debris, is 0.109002, 0.108987 and 0.108985.
- The three precisions do not follow the same trajectories (the particle numbers at the end
  differ by up to 1.4e-4), so the differences above are not a rounding bound of the mode.

Not covered by these runs: a build with reflecting walls or periodic axes at this size (the
tests and the benchmark above cover them at their sizes), and a production-size run.
