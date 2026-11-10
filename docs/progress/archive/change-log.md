# Archived Change Log

Older `PROGRESS.md` Change Log entries, newest first. Frozen: insert only at the marker
(procedure: `docs/progress/archiving.md`).

<!-- ARCHIVE-INSERT-POINT -->

### 2026-11-10 — [progress] Performance round 3 gated: ALL PASSED on 158bc20, fp64 and mixed
- What: gates on the stack #27–#37: growth error +3.0 / −1.7 / −2.3 / +3.5 / −4.7 %, Gauss max 7.6e-12, implicit energy 1.9e-10, 15 example gates, 71 GPU tests; mixed: 19 explicit gates, Gauss max 7.7e-12. Recorded without a criterion: `div B` of `plasma_waves:perpendicular` 1.00e-12 (0.95e-12 mixed), in its earlier range. 3e013df → 158bc20, ms/step (two alternating runs): two-stream 1/2/4/8 GPUs 37.9/19.7/10.2/5.7 → 22.0/11.5/6.0/3.4 (−42/−41/−41/−40 %), mixed 33.5/17.5/9.3/5.3 → 18.2/9.9/5.3/3.1, shock 8 GPUs 10.9 → 6.3 (−43 %), implicit 1 GPU 513 → 488 (−4.9 %) and 8 GPUs 101–103 → 76–78 (−22 to −27 %, 60-step windows), production run with the full output set 40.6 → 23.0 (1 GPU) and 6.4–6.5 → 3.8 (8 GPUs).
- Why: evidence for the round-3 changes on one clean commit (CI cannot run the gates). M6 stays Partial: the push kernel is 71 % of the 8-GPU step (94 % on 1 GPU) and FP64-bound, and the open items change the arithmetic (R10, R5, PIC-D11).
- Refs: #4, #27–#37; `bench/results/2026-11-10-158bc20b8-gpu-gates.json`, `…-gpu-gates-mixed.json`; `docs/roadmap.md` R3/R4/R9, `README.md`
### 2026-11-10 — [progress] Performance round 3: items measured and not adopted
- What: below the criterion (2 % of the step time of the target case, alternating pairs) and discarded: phase-space histogram in one kernel (ceiling 1.0 % of the production run on 1 GPU, 1.3 % on 8), track selection in a kernel (0.6 %), field dumps in one gather (0.2 %), writer thread on rank 0 (−4.0 / −1.5 / −1.3 %, one pair of three; OUT-D2 unchanged), hot/cold split of `implicit_push` (`T(maxit=100) − T(4)` = 0.3 ms per launch, 0.5 % of the kernel), quadratic Picard start (start error 0.45–0.88 against 0.24–0.58 linear: worse), Anderson / plasma preconditioning (Picard contraction ≈ 1.8e-4 per sweep: no sweep to save), fusing `cg_scalars` (ceiling 4.3 µs per iteration), CUDA Graph for the CG (2 ranks: local kernels of an iteration 44 µs, 43 µs replayed, of 159 µs; a captured NCCL exchange does not replay), paired filter passes on a decomposed axis (exchanges 5 → 3, bitwise equal; shock −0.1 / −0.3 / −1.7 %). Closed without a measurement: tiled deposition (deposit atomics are 0.2 % of the PC samples, round-2 `ncu`), tuning for `diag_interval=1` (the production value is 12).
- Why: every task of the round is closed as adopted or as measured and rejected, so that none is proposed again without a new fact.
- Refs: #4 (measurements: issuecomment-6763298161); `docs/roadmap.md` R3/R4/R9
### 2026-11-10 — [progress] Correction: the round-2 output and `diag_interval=1` numbers hold for the bench schedule only
- What: corrects "Diagnostics step: field energy on the device, one reduction per record (issue #23)" and "Performance round 2 gated: ALL PASSED on 5d7015a, fp64 and mixed" (both archived). `diag_interval=1` is not the planner default: the planner gives 12 for the production two-stream run and 1 only under the bench override `t_end=4` (226 steps), where the outputs requested as a count per run also come every 2–37 steps instead of 49–824. Production run (4936 steps, full output set) against the same run without output: 3e013df +2.8 % on 1 GPU (1.1 ms/step) and +12–14 % on 8 GPUs (0.74 ms/step); 158bc20 +4.0 % (0.89 ms/step) and +13–14 % (0.45 ms/step).
- Why: the override shortens the run, not the number of outputs per run; the cost of the output is timed on the production schedule from now on.
- Refs: #4; `docs/gpu-runbook.md`, `docs/design/output-analysis.md` (Measured), `docs/roadmap.md` R9
### 2026-11-10 — [progress] Current filter: all passes along a self-periodic axis in one launch (issue #37)
- What: along an axis that wraps onto the rank the N passes (+ compensator) are one composite stencil (`binomial_wrap3`, indices modulo the period, exact dyadic coefficients, evaluated symmetrically), with no ghost fill; decomposed and wall axes are unchanged. Owned cells differ from the passes one by one at rounding level; the ghosts of a wrapped axis now hold the periodic images. Shock (`pair_2d_long`, 8 GPUs), 3 alternating pairs vs 7fe5815: −9.0 / −10.3 / −8.8 % (6.87 → 6.26 ms/step); synchronised `filter_J` 1.61 → 1.01 ms/step.
- Why: 10 of the 15 sweeps of the shock's filter were along y and z, each a ghost fill and a pass; the filter was 21 % of the step.
- Refs: #37, #4; `pypic/kernels/fields.cu`, `pypic/grid.py`, `docs/design/engine.md` (PIC-D10); `tests/gpu/test_engine.py` (`test_current_filter_kernels_match_the_array_expressions`)
### 2026-11-10 — [progress] Particle sort without a second copy of the particle arrays (issue #36)
- What: `Species.sort` gathers the rows one at a time into the key array (free once the order is known) instead of building `P[:, order]`. Same particles in the same order. One species of 25 million particles: extra pool memory of a sort 2.43 → 0.94 GiB (−61 %), 21.7–23.4 ms per sort before and after; synchronised `sort` phase 0.439 → 0.439 ms/step (two-stream, 1 GPU) and 0.077 → 0.069 (shock, 8 GPUs).
- Why: the copy was the peak of the device memory and bounds the particles per GPU; judged by the memory, with the time as a check.
- Refs: #36, #4; `pypic/particles.py`; `tests/gpu/test_engine.py` (`test_sort_orders_by_cell_and_keeps_every_particle`)
### 2026-11-10 — [progress] Migration lists ordered by a counting kernel instead of a radix sort (issue #35)
- What: the two lists of a migration (leaving indices, destination keys) are ordered by `mig_rank32` / `mig_rank64` — one launch, each key to the slot given by the number of keys before it — up to `RANK_SORT_MAX` = 16384 keys, by the radix sort above. Same values in the same order: the migrated particles are unchanged. 3 alternating pairs vs 6061977 on 8 GPUs: two-stream fp64 −1.8 / −3.8 / −4.0 %, mixed −2.9 / −5.1 / −4.4 %, shock −2.3 / −2.3 / −1.0 %; synchronised `migrate` 0.54 → 0.45 (two-stream) and 0.54 → 0.41 ms/step (shock). Not adopted: the count inside the data message (ceiling 99 µs; the overflow round must be a collective, net ≈ 1.4 % of the step).
- Why: the two sorts were 181 µs of the 606 µs of `migrate` per step; a radix sort of 2 500–14 000 keys costs 24–86 µs, the counting kernel 14–60 µs.
- Refs: #35, #4; `pypic/kernels/migrate.cu`, `pypic/particles.py` (`ascending`); `tests/gpu/test_engine.py` (`test_migration_lists_are_ordered_like_a_sort`)
### 2026-11-10 — [progress] Explicit field update: ghost layer of B computed locally, 3 → 1 exchanges per step (issue #34)
- What: `update_B` runs on the update range and one ghost layer around it (from the two exchanged ghost layers of E); `push_B` no longer exchanges. Its curl products are rounded one by one (`__dmul_rn`, explicit `fma`): the kernel variants with and without the absorbing layer were contracted differently, which left the ghosts 1 ulp off between an edge rank and a middle rank on an open decomposed axis. B changes at rounding level. 3 alternating pairs vs ac1b7ff on 8 GPUs: two-stream fp64 −1.0 / −3.0 / −4.6 % (3.66 → 3.55 ms/step), mixed −2.4 / −2.2 / −2.5 %, shock −0.8 / −1.0 / −1.4 %; synchronised `field_push` 0.36 → 0.17 ms/step.
- Why: the two exchanges of B were 198 µs per step on 8 GPUs, 55 % of the field update.
- Refs: #34, #4; `pypic/kernels/fields.cu`, `pypic/grid.py`, `docs/design/engine.md` (PIC-D5); `tests/gpu/test_engine.py` (`test_locally_computed_ghost_layer_of_B_equals_the_neighbours_values`)
### 2026-11-10 — [progress] Explicit push: Esirkepov and gather with 46 % fewer executed instructions (issue #33)
- What: `esirkepov` takes the per-axis CIC weights from integer cell conditions (no `cic()` compares), builds each current from the prefix sum over two faces (the second only when the move crosses a cell face) and factorises the transverse weight (`W = S0' A + dS' B`); the gather reads the 8 corners at constant offsets from one index. Scheme and FP64 unchanged, rounding order only. 3 alternating pairs vs a2eca7a on 8 GPUs: two-stream fp64 5.73 → 3.63 ms/step (−35 to −37 %), mixed 5.21 → 3.33, shock 10.79 → 6.95; 1 GPU fp64 37.9 → 22.0, mixed 33.5 → 18.2. Gauss error unchanged (7.3e-14 → 7.4e-14 periodic, 4.7e-15 → 5.3e-15 walls).
- Why: `ncu` — a warp executes the union of all 27 node paths per component, each with FP64 compares and a 7-operation weight, and a 64-bit 3-D index per access (5726M → 3116M executed thread-instructions per launch of 4.19M particles).
- Refs: #33, #4; `pypic/kernels/particles.cu`, `pypic/kernels/common.cuh`
### 2026-11-10 — [progress] Implicit CG: matrix product in one launch on periodic grids (issue #32)
- What: new kernel `cg_matvec` — where no rank has a wall or an absorbing layer, `w = a u + (dt²/4) curl curl u` and the three inner products are one launch instead of `curl_E`, `curl_B`, `cg_matvec_end` and a temporary field; with walls the kernel sequence is unchanged. Agrees with the sequence to rounding. Implicit two-stream, 4 alternating pairs vs a2eca7a on 8 GPUs: −1.8 / −4.2 / −5.3 / −3.3 % (`cg` 18.60 → 17.05 ms/step); 1 GPU −0.1 %.
- Why: on the per-rank grid of 8 GPUs the three launches cost 42 µs per product (74 products per step), 12 µs of it the last kernel alone.
- Refs: #32, #4; `pypic/kernels/fields.cu`, `pypic/simulation.py`; `tests/gpu/test_engine.py` (`test_cg_matrix_product_in_one_launch_matches_the_kernel_sequence`)
### 2026-11-10 — [progress] Implicit half-step fields and Picard change in dedicated kernels (issue #31)
- What: `half_E` builds `Eh = (E0 + E1)/2` on every point, ghosts included (the ghosts of `E0` and `E1` are already those of the neighbours), so `_half_fields` drops its exchange of `Eh`; `half_B` gives `Bh = B0 - dt/2 curl Eh` in one launch; `picard_change` takes the block maxima of `|Enew - E1|` and `|Enew|` (NaN kept) instead of four whole-grid temporaries. Bitwise the previous values. Implicit two-stream, 4 alternating pairs vs 3c2532b on 8 GPUs: −4.8 / −3.0 / −5.8 / −3.1 % (`_half_fields` 2.82 → 0.99 ms/step); 1 GPU −0.4 %.
- Why: five half-step builds and Picard checks per step were launches, temporaries and one NCCL exchange each, on a per-rank grid 8× smaller than on 1 GPU.
- Refs: #31, #4; `pypic/kernels/fields.cu`, `pypic/simulation.py`; `tests/gpu/test_engine.py` (`test_half_step_fields_and_picard_change_match_the_array_expressions`)
### 2026-11-10 — [progress] Implicit particle map: gamma sweeps stop by a contraction bound (issue #30)
- What: in `particle_map` the sweeps over `g1` also stop when `L = 2|b||av|/(g0+g1)^2 <= 1/2` (`b = 2 eps B`, `av = u0 + eps E`) and `L dg <= 1e-16 |u1|`: the remaining sweeps then change `u1` and `g1` by at most `2 L dg`. The sweep is `ub = (m^2 av + m av x b + (av.b) b)/(m^2 + b^2)`, one division instead of two. Sweeps per map evaluation: one 29 → 75 %, two 67 → 25 %, three 4 → 0 %. Implicit two-stream, alternating pairs vs 50f41cf: 8 GPUs −3.2 / −4.8 / −6.2 % (87.7 → 84.8, 87.9 → 83.7, 86.9 → 81.5 ms/step), 1 GPU −4.5 / −4.4 % (513 → 490). Scheme, tolerances, energy closure and Picard counts unchanged.
- Why: `implicit_push` is 90 % of the 1-GPU implicit step; two of three map evaluations spent a second sweep (two FP64 divisions) only to confirm that `g1` no longer changes.
- Refs: #30, #4; `pypic/kernels/implicit.cu`; `tests/gpu/test_engine.py` (`test_implicit_momentum_is_the_fixed_point_of_the_midpoint_rotation`)
### 2026-11-10 — [progress] Ghost fill in three launches (issue #29)
- What: `slab_pack3` packs / unpacks both faces of an axis in one launch (one buffer, two NCCL messages, same order); new `slab_wrap3` fills the ghosts of all axes that wrap onto the rank in one launch from the periodic images, after the exchanges — the values of the axis-by-axis sequence, bitwise. `fill_ghosts` 8 → 3 launches per rank with one decomposed axis; `reduce_ghosts` shares the two-face pack.
- Why: one `fill_ghosts` per CG matrix product (74 per step). Implicit 8 GPUs −3.0 / −3.9 / −4.0 / −3.2 % (4 alternating pairs; CG 21.1 → 18.6 ms/step); shock −1.2 %, explicit within noise, 1 GPU −0.3 %.
- Refs: #29, #4; `pypic/grid.py`, `pypic/kernels/fields.cu`; `tests/gpu/test_engine.py` (`test_fill_ghosts_gives_the_periodic_images`)
### 2026-11-10 — [progress] Implicit CG: one reduction per iteration (issue #28)
- What: `_cg` uses the Chronopoulos–Gear recurrence: the last kernel of the matrix product (`cg_matvec_end`) also sums `r·u`, `w·u`, `r·r` per thread block, one `allreduce` carries the three, `cg_scalars` computes α and β on the device and `cg_step` does the four vector updates and the preconditioning (before: 2.25 reductions and 3 vector kernels per iteration). Implicit two-stream, 8 GPUs, two alternating pairs: 97.5 → 89.3 and 98.0 → 89.7 ms/step (−8.4 / −8.5 %); 1 GPU 513.3 → 515.2 / 514.9 (+0.3 %: one more matrix product per solve). Iteration counts and the true residual equal the old CG (periodic / reflecting / open, `dt` 0.5–2).
- Why: 163 NCCL `AllReduce` per step were 25 ms of the 8-GPU implicit step (`nsys`); the CG took longer on 8 GPUs than on 1.
- Refs: #28, #4; `pypic/simulation.py`, `pypic/kernels/fields.cu`, `pypic/kernelsrc.py` (`CG_BLOCK`); `tests/gpu/test_engine.py`
### 2026-11-10 — [progress] Time-averaged lineout and modes share one profile per step (issue #27)
- What: `FieldSource.profile(comp, axis)` sums a component over the two other axes once per step and shares it; `LineoutDiag` and `ModesDiag` accumulate that real profile in place, and the mode transform is applied once per output (linear, so the same result). Two-stream production run with output (4936 steps, 8 GPUs, two alternating pairs): 6.62 → 6.17 and 6.42 → 6.19 ms/step (−6.9 / −3.6 %); the per-step accumulation 0.364 → 0.118 ms/step. The 13 output files agree with 3e013df within 1e-12.
- Why: the averaged diagnostics accumulate every step; each did its own reduction of the whole grid with temporaries and a complex transform per step — the largest part of the output cost of a production run (8–10 % of the run on 8 GPUs).
- Refs: #27, #4; `pypic/output/diagnostics.py`
### 2026-11-10 — [progress] Correction: `div B` of `plasma_waves:perpendicular` is within its earlier range
- What: the entry below ("Performance round 2 gated") says the value was 0.85–0.95e-12 "in the three records before". That is wrong: the 18 earlier fp64 gate records span 0.85–1.17e-12, eight of them ≥ 1.0e-12 (77f26df, two records before, has 1.17e-12). The 1.02e-12 on 5d7015a is inside that range; `push_B` and the `div B` diagnostic are unchanged, and no regression is observed. The case has no `div B` criterion.
- Why: found in the review of PR #26; the wrong comparison read as a drift towards the 1e-12 contract.
- Refs: PR #26 review; `bench/results/*-gpu-gates.json` (`plasma_waves:perpendicular`, `divB_max`)
### 2026-11-10 — [progress] Performance round 2 gated: ALL PASSED on 5d7015a, fp64 and mixed
- What: gates on the stack #19–#25: growth error +3.0 / −1.7 / −2.3 / +3.5 / −4.7 %, Gauss ≤ 7.6e-12, implicit energy 1.9e-10, 15 example gates, 58 GPU tests; mixed: 19 explicit gates. Baseline f31e7f9 → 5d7015a, ms/step (200-step windows, two alternating runs): two-stream 1/2/4/8 GPUs 43.1/25.1/14.6/9.5 → 37.9/19.7/10.2/5.6 (mixed 39.0/9.5 → 33.5/5.3 on 1/8), `diag_interval=1` 8 GPUs 10.9 → 6.6, shock 14.9 → 10.9, implicit 597/137 → 513/94 on 1/8. Measured without gain and not adopted: compiling out the loss reductions, `sort_interval` tuning, sort rework, paired filter passes (−0.7 %). Recorded without a criterion: `div B` of `plasma_waves:perpendicular` (B0 = 0.5) is 1.02e-12, 0.85–0.95e-12 in the three records before; the gated two-stream cases are ≤ 1.7e-14.
- Why: evidence for the round-2 changes on one clean commit (CI cannot run the gates), and the same-window before/after numbers for the roadmap.
- Refs: #4, #19–#25; `bench/results/2026-11-10-5d7015a34-gpu-gates.json`, `…-gpu-gates-mixed.json`; `docs/roadmap.md` R3/R4
### 2026-11-10 — [progress] Launcher binds each rank to the NUMA node of its GPU (issue #25)
- What: `pypic.launch` sets each rank's CPU affinity to the CPUs of the node its GPU is attached to (PCI address → `numa_node` → `cpulist`), within the CPUs the launcher may use; `--no-bind` turns it off; unknown node → no binding. 8 GPUs, default vs `--no-bind`: implicit two-stream 96.98 → 94.65 ms/step (−2.4 %, 4 pairs; −2.2 / −2.8 % in two earlier batches with external pinning), `diag_interval=1` 6.48 → 6.44 (−2.7 % with external pinning), two-stream −1.5 %, shock −0.5 %.
- Why: 4 NUMA nodes with two GPUs each; the host thread that drives a GPU ran wherever the scheduler put it. The gain is on the cases with the most host work per step; bound runs also spread less (implicit 94.5–94.8 vs 95.8–99.2 ms).
- Refs: #25, #4; `pypic/launch.py`; `tests/gpu/test_binding.py`
### 2026-11-10 — [progress] Explicit step queues the field update before the migration (issue #24)
- What: `_step_explicit` calls `push_B`/`push_E`/`push_B` before `_migrate()`. 8 GPUs: two-stream 6.19 → 5.90 ms/step (−4.7 %, 4 pairs), mixed 5.54 → 5.34, `diag_interval=1` 6.69 → 6.44, shock 11.02 → 10.83 (−1.8 %). A short explicit history agrees with the previous commit to rounding.
- Why: the migration starts with a host read that waits for the push kernels; the field update (J and fields only, independent of the particles the migration moves) was prepared after it, serial with the push instead of behind it.
- Refs: #24, #4; `pypic/simulation.py`
### 2026-11-10 — [progress] Diagnostics step: field energy on the device, one reduction per record (issue #23)
- What: `YeeGrid.field_energy_local` sums `E^2` and `B^2` with one reduction kernel each and leaves the result on the device; `_record` reduces it together with the per-species diagnostics in its one `allreduce` and one host read (before: six host reads, two `allreduce`). Two-stream with `diag_interval=1` (the planner default), 8 GPUs: 7.57 → 7.09 ms/step (−6.4 %, 6 pairs); 1 GPU 41.46 → 41.09. Recorded energies agree with the previous commit to rounding (rel ≤ 2e-16).
- Why: with the history recorded every step the field energy and `_record` were 14 % of the 8-GPU step, almost all of it host synchronisation.
- Refs: #23, #4; `pypic/grid.py`, `pypic/simulation.py`, `scripts/bench_step.py`
### 2026-11-10 — [progress] implicit_push with reciprocals instead of FP64 divisions (issue #22)
- What: the particle map of the implicit kernel computes `eps / gbar`, `1/(1+t^2)` and `DT/(g0+g1)` once each, the Picard residual norms multiply by `1/h`, and `eps = q dt / 2m`, `q/m` come from the host. Implicit two-stream 1 GPU 586 → 513.5 ms/step (−12.4 %), 8 GPUs 109.5–110.6 → 99.7–101.0. Picard iterations (145 in 30 steps) and the energy budget error (3.03e-11) are unchanged.
- Why: same reason as #21 — the kernel is FP64-pipe bound and the map runs once per sub-iteration of every particle. Rounding-level change, scheme unchanged.
- Refs: #22, #4; `pypic/kernels/implicit.cu`, `pypic/simulation.py`
### 2026-11-10 — [progress] push_deposit with reciprocals instead of FP64 divisions (issue #21)
- What: `1/gm`, `2/(1+t^2)` and `DT/g1` are one reciprocal each (9 divisions → 3), `h = q dt / 2m` comes from the host, the periodic wrap is a comparison (same value for moves below one period), `g0` only on diagnostics steps. Two-stream 1 GPU 43.0–43.9 → 38.0–38.1 ms/step, 8 GPUs 6.6–6.8 → 6.1–6.3 (mixed 6.1–6.2 → 5.5), shock 8 GPUs 12.0–12.1 → 10.9–11.0. Not adopted (no gain): compiling out the open-boundary loss reductions.
- Why: the kernel is FP64-pipe bound and 75–96 % of the explicit step; divisions are its most expensive operations. Rounding-level change, scheme unchanged.
- Refs: #21, #4; `pypic/kernels/particles.cu`; `test_boris_gyration_conserves_gamma_and_period`
### 2026-11-10 — [progress] Implicit field solve: one ghost exchange per CG matvec (issue #20)
- What: `curl_E` is also taken on the ghost layer that `curl_B` reads, so the matvec exchanges ghosts once (bitwise the same on the owned unknowns); the CG vector updates are three fused kernels (`vec_add_scaled`, `cg_update`, `cg_direction`); the unused `itmax` atomic of `implicit_push` is gone (no measurable effect). Implicit two-stream, 8 GPUs: `_cg` 58–60 → 31 ms, step 136–145 → 108–110 ms/step; 1 GPU 597 → 586. CG iterations unchanged.
- Why: on 8 GPUs the CG was 42 % of the implicit step, host-side sequencing of two exchanges and ~10 small launches per iteration, not grid work.
- Refs: #20, #4; `Simulation._matvec`, `Simulation._cg`; `pypic/kernels/fields.cu`
### 2026-11-10 — [progress] Particle migration by dedicated kernels (issue #19)
- What: `particles.migrate` takes, fills, classifies, packs and appends with five kernels (`pypic/kernels/migrate.cu`) on a particle-major buffer; one host sync fewer and no allreduce; particle order unchanged (old vs new arrays bitwise equal on 1/4/8 ranks). 8 GPUs: migrate 3.8 → 0.6 ms; two-stream 9.7–10.0 → 6.7–6.8 ms/step (mixed 9.2–9.4 → 6.2–6.5), shock 14.7–15.0 → 11.8.
- Why: migration was 23–40 % of the 8-GPU step: host time of ~100 small CuPy launches for a few thousand particles (0.3 ms on 1 GPU).
- Refs: #19, #4; `scripts/bench_step.py`; `tests/gpu/test_engine.py::test_migration_across_faces_and_corners_keeps_every_particle`
### 2026-11-10 — [design] PIC-D11 mixed-precision field gather, Proposed (issue #10)
- What: `precision: mixed` gathers E, B for the explicit push from float32 copies with rank-independent weights; positions, pusher and deposit stay FP64; default `fp64`. Push 18.4 → 14.2 ms; two-stream 1 GPU 42.7 → 38.6 ms/step; shock 8 GPUs 14.4 → 13.4. All 19 explicit gates pass in mixed mode.
- Why: `push_deposit` is FP64-pipe bound; the gather is the part that can leave FP64 without touching Gauss (deposit) or 1 vs N GPU (global weights).
- Refs: #10, #4; PIC-D11 in `docs/design/engine.md`; `test_mixed_gather_is_decomposition_independent`
### 2026-11-10 — [amendment] Esirkepov deposit on the 3-node CIC support (issue #6)
- What: the explicit deposit loops over the 3 nodes per axis that a sub-cell CIC move touches (2 prefix sums per component), and the gather shares its weights. Electron `push_deposit` 37.9 → 30.8 ms; two-stream 1 GPU 51.6 → 42.7 ms/step, 8 GPUs 10.7 → 9.2–10.2; shock 8 GPUs 16.0 → 14.6–14.9.
- Why: `push_deposit` is FP64-pipe bound (85 %) on SM120; same scheme, fewer FP64 operations.
- Refs: #6, #4; PIC-D2 review amendment in `docs/design/engine.md`; `test_esirkepov_current_moment_equals_the_displacement`
### 2026-11-10 — [progress] Ghost exchange with fused slab kernels (issue #8)
- What: ghost copy / add / zero / pack / unpack of up to 3 arrays per launch, one message per face; a self-wrapping axis is copied in place. 8 GPUs: two-stream 11.7 → 10.5–10.9 ms/step (diag every step 13.7 → 11.6–12.4), shock 17.8 → 16.0 ms/step.
- Why: the 8-GPU step was bound by host-side sequencing of ~30 small operations per ghost call, not by bandwidth; no compute/communication overlap (the field kernels take < 0.1 ms per step, nothing to hide behind).
- Refs: #8, #4
### 2026-11-10 — [progress] Batched particle migration (issue #7)
- What: `particles.migrate` moves every species at once: one message group per decomposed axis, sort-based partition, 3 host syncs per step (was ~20 for two species), one allreduce for the diagnostics record. 8-GPU two-stream: NCCL calls 16.2 → 14.2 per step; step time unchanged within noise (11.3–12.2 → 11.7 ms).
- Why: the issue's premise; the profile after it shows the 8-GPU step is bound by host-side sequencing (GPU idle ~5 of 14 ms under nsys), not by the syncs' count.
- Refs: #7, #4; `tests/gpu/test_engine.py::test_migration_across_faces_and_corners_keeps_every_particle`
### 2026-11-10 — [progress] Current filter as fused ping-pong pass kernels (issue #9)
- What: one kernel per pass filters up to 3 arrays, passes alternate with work buffers; bitwise the former result. Shock `pair_2d_long` on 8 GPUs: filter cost 2.1 → ~1.6 ms of a 17.8 → ~17.6 ms step; ~250 → ~35 launches per step.
- Why: the filter was ~250 small launches and two sweeps per pass; it is memory-bound (~1 TB/s per pass), and the halo exchange per pass (15 per step) is now most of what remains (follow-up: two passes per sweep).
- Refs: #9, #4; PIC-D10; `tests/gpu/test_engine.py::test_current_filter_kernels_match_the_array_expressions`
### 2026-11-10 — [progress] CG inner product in two stages (issue #5)
- What: `Simulation._dot` reduces over the inner axes, then sums; implicit two-stream (50.3M particles) 1017 → ~570 ms/step on 1 GPU, 151 → 144 ms/step on 8 GPUs.
- Why: a CuPy ReductionKernel to one scalar runs as one thread block (3.1 ms per call, 49 % of the implicit GPU time).
- Refs: #5, #4; `tests/gpu/test_engine.py::test_cg_inner_product_matches_the_owned_unknowns`
### 2026-11-10 — [progress] Particle row offsets are 64-bit (issue #11)
- What: the particle kernels take `cap` as 64-bit, so `row * cap` no longer overflows above ~2.7e8 particles of one species per GPU (it was an illegal address). `Species` rejects a capacity above 2**31 - 1 (the kernels' int particle index).
- Why: a 96 GB GPU holds ~1e9 particles; the int32 offsets capped a species at ~1/4 of that (R11 needs large runs). The regression test failed on the old code.
- Refs: #11, #4; `tests/gpu/test_engine.py::test_particle_rows_beyond_int32_offsets`

### 2026-11-10 — [progress] GitHub wiki articles of the examples (EN + JA)
- What: one page per example and language: conditions table, GIF animations and comparison charts with theory or the published results, all generated from the runs (`scripts/build_wiki.py run|collect|render`, templates and numbers in `docs/wiki/`). New charts in the examples (growth ratio, branch errors, convergence, jump conditions, p vs σ, reconnection field lines with X/O points); `pypic.analysis.animate` lays out with the first frame's titles and spans the spectrum axis over all frames. README numbers corrected.
- Why: owner request; the numbers come from runs (fresh gate configs at a8caa8f, recorded production runs), so the pages cannot drift from the data.
- Refs: `docs/wiki/`, `tests/cpu/test_wiki.py`, `bench/results/2026-11-10-d03ab601a-gpu-gates.json`
### 2026-11-10 — [amendment] PR #2 fifth review: rows of a failed HDF5 flush
- What: `SeriesWriter._write` keeps the rows buffered until the HDF5 flush succeeds, so a retry repeats the write and the flush.
- Why: the review found that after a failed flush the next flush() returned without flushing, and SWMR readers saw no row. The new CPU test failed on that code.
- Refs: review amendments in `docs/design/output-analysis.md`; `tests/cpu/test_output.py`
### 2026-11-10 — [amendment] PR #2 fourth review: scalar rows of a failed write
- What: scalar history rows count as written once the scalars writer has taken them (rank 0), not before the write. A failed write is retried, and rows already taken stay buffered.
- Why: the review found a regression from the third-review change: after a failed write, a retry wrote nothing. The new test failed on that code.
- Refs: review amendments in `docs/design/output-analysis.md`; `tests/gpu/test_output.py`
### 2026-11-10 — [progress] Second recorded realizations of the production runs
- What: `pair_sigma10_large` (run e879cfd6) and `pair_2d_long` (run 40c0ca85) were run again and recorded, one record per run.
  - Reconnection, two recorded runs: p = 2.10 / 2.10 and rate 0.064 / 0.074 v_A.
  - Shock, two recorded runs: β_sh −5.9 / −7.6 %, n₂/n₁ +6.9 / +8.1 %, tail p 3.91 / 3.85, γ_max ∝ t^0.36 / t^0.41.
- Why: the 2026-11-09 entry "Production runs of the acceleration examples" quoted a first realization of each that was never recorded (second review). Those values are withdrawn: p = 2.15 and rate 0.080; −4.2 %, +7.8 %, p 3.93 and t^0.42.
- Refs: `bench/results/*-run-e879cfd6.json`, `*-run-40c0ca85.json`; `examples/{reconnection,relativistic_shock}/README.md`; ANA-D3
### 2026-11-10 — [amendment] PR #2 third review: rank-0 file errors, untracked code
- What:
  - Every rank-0 file operation of the output shares its outcome with all ranks (`OutputError` on every rank).
  - Provenance has one definition (`pypic.provenance`), used by run.json, gate records and example records. Untracked code files make a commit dirty; changes outside the code paths do not.
- Why: the review reproduced a rank-0-only `PermissionError` that left the other ranks in a collective, and an untracked analysis helper that was missing from `analysis_commit`. Each fix has a test that failed on the old code.
- Refs: OUT-D2, ANA-D3 and review amendments in `docs/design/output-analysis.md`

### 2026-11-09 — [amendment] PR #2 second review: recorded results
- What:
  - `run.json` carries a run id and a start time, both fixed at start-up. Example records are one per run, with the run id in the file name and index entry, so realizations no longer overwrite each other.
  - The analysis commit is taken before the analysis, and nothing is recorded if the code changes during it.
  - Evidence is strict JSON (NaN → null; the five existing files were normalized), and the verifier enforces it.
- Why: the review reproduced 3 issues; each fix has a test that failed on the old scripts.
- Refs: ANA-D3 and review amendments in `docs/design/output-analysis.md`; `scripts/bench_records.py`
### 2026-11-09 — [amendment] PR #2 review fixes
- What:
  - The density-perturbation inverse is a bisection with a residual check reduced over ranks. Newton stalled at α → 1.
  - A test species' virtual move deposits `qdep`.
  - One run per output directory (`output.overwrite` / `--overwrite`).
  - Also fixed: empty tracks, typed `output.*` overrides, `fit_two_wave` (finite weights, coefficients at t = 0), the work split at B = 0, the drift-frame orbit theory up to t_end, and `record_example.py` (repository examples only).
- Why: the review reproduced 10 bugs (1 P1, 9 P2). Each fix has a regression test that failed on the old code.
- Refs: review amendments in `docs/design/engine.md` and `docs/design/output-analysis.md`; PR #2
### 2026-11-09 — [progress] Production runs of the acceleration examples
- What:
  - pair_3d: 3-D jump conditions, β_sh −3.7 %, n₂/n₁ +3.5 %.
  - pair_2d_long: early Fermi tail p ≈ 3.9 at t = 2000, γ_max ∝ t^0.36–0.42.
  - isothermal_hydrogen_long: E_max −12 % up to τ = 30.
  - pair_sigma10_large: p = 2.10–2.15 vs Werner+2016's 2.12; rate 0.06–0.08 v_A (two realizations: the atomic deposit order is not deterministic).
  - All recorded as typed example evidence (`scripts/record_example.py`, validated in CI).
- Why: the reported (non-analytic) results planned for PR #2.
- Fixes: the shock speed is fitted over 0.3–0.7 t_valid (the precursor reaches the plasma edge later); run.json takes the git state once at start-up.
- Refs: `bench/results/*-example-*.json`; `examples/{relativistic_shock,reconnection,plasma_expansion}/README.md`

### 2026-11-09 — [design] Analytic-solution examples and their set-up features (PR #2, phases 2-4)
- What: examples plasma_waves, landau_damping (+ Dawson, wave breaking, relativistic), filamentation, particle_orbits, plasma_expansion, relativistic_shock, reconnection, each gated against theory (ANA-D2); engine PIC-D9 (Gauss-exact perturbations, explicit half kick, test species and particle lists, E0 / exact plane waves / field profiles, density profiles) and PIC-D10 (binomial current filter, Gauss- and implicit-energy-exact); `pypic.theory` waves/streaming.
- Why: owner request (one example per analytic problem); each feature keeps the conservation contracts (GPU tests: Gauss <= 1e-12, implicit energy <= 1e-11, 1 vs 2 GPUs).
- Findings: Bernstein needs theta_e = 1e-3 (relativistic shift at 1e-2); PIC expansion ~15 % below Mora's isothermal E_max (kinetic electrons; gate set above the measured offset); Buneman needs remove_mean_current.
- Refs: PIC-D9, PIC-D10, ANA-D2; `examples/*/README.md`; plan `docs/plans/2026-11-09-output-analysis-examples.md`

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
