# The ids of the entering particles: one block per step of the size of the layers

Status: Implemented and gated on 58166d1 (PR #53; not merged)
Refs: PIC-D13 (amended: "Ids"), PIC-D25 (amended: `meta.json`), PIC-D8

## Context

The research topic `qpar-electron-conditions` (plan approved by the owner on 2026-11-14) runs long quasi-1-D
piston shocks with one inflow face. Step 4 of that plan names the route for an engine shortfall: a branch from
main, one pull request, the GPU gate, a Codex review with no findings, the owner's instruction to merge.

The shortfall (main 70202df): the particles that enter through open faces (PIC-D13) take one block of ids per
step, and `ids_per_step` was `2 · streams · ppc · Π (N_a + 2)`: one id per slot of the whole extended grid, per
stream, and once more for the partners of a paired face. Only the cells of the layers outside the injecting faces
are ever used. For a box that is long along the normal of its inflow face the block is larger than the layers by
the length of the box:

| run | cells | ppc | streams | steps | ids needed / 2^53 | stops at step |
|---|---|---|---|---|---|---|
| the topic's V1 | 270000 × 2 × 2 | 250 | 2 | 5,000,000 | 2.40 | 2,084,984 |

(`loaded = 5.4e8`, `ids_per_step = 4,320,032,000`.) The engine stops such a run on every rank when the next block
would reach 2^53 (the step check of PIC-D13), so the run cannot be done.

## Approach

**Block = the planes of the injecting faces.** The plane of a face on axis `a` is the extended grid (cells −1..N)
of the two other axes. It contains the layer of the face with the edge and corner cells that the face fills
(`injection.face_slabs` gives an edge or corner cell to the first of its faces in the order x, y, z).
- `injection.id_planes(cfg)`: per injecting face the first cell of its plane in the sub-block, and the cells of
  all planes. `ids_per_step = 2 · streams · ppc · (cells of all planes)`; 0 without injection.
- id = `n_ids + step · ids_per_step + (stream sub-block) + (first cell of the face's plane + cell in the plane) · ppc
  + slot`; the partners of a paired face use the second half of the block as before.
- Kernel `inject_face`: `lin` is the index of the layer cell in the plane of the face (it was the index in the
  extended grid). No new argument: `_inject` adds the first id of the plane to `id_step`. No new `#define`.
- The key of the random numbers (`gcell`, from the global cell) is not touched: positions, momenta, weights, the
  bank and the fields of a run keep their bits; only the ids of the particles that entered change.
- For V1: one plane of 4 × 4 cells, `ids_per_step = 16,000`, ids needed 8.9e-6 of 2^53.

**Checkpoints.** A checkpoint written with the larger block cannot be continued with the smaller one: the ids of
step n would start at `n_ids + n · (new block)`, below ids given before the checkpoint. `checkpoint.save` stores
`ids_per_step` in `meta.json`; `check_restart` refuses a run with injection when the stored value is missing
(a checkpoint written before this change) or differs. A run without injection is read as before. The format
name stays `pypic-checkpoint-1`: nothing else in a checkpoint changed, and a bump would refuse every old
checkpoint, also those without injection.

Rejected: a block of exactly the layer cells (no unused ids). It needs the per-rank slab extents in the kernel's
id, saves less than a factor (N + 2) / N per axis against the planes, and the planes already remove the factor
that made runs impossible (the length along the face normal).

## Invariants

- Ids (PIC-D8) stay unique and independent of the decomposition: the id depends on (step, face, global cell in
  the plane, stream, slot); each layer cell belongs to one face on every layout.
- The step check of PIC-D13 (every rank raises together, from the configuration and the step) is unchanged and
  uses the new block.
- Gauss, div B, energy: not touched (no deposit, push or field code changes).
- Engine stays experiment-agnostic: nothing of the topic enters `pypic/`.

## Verification

Fixed before the first run of each test.
- CPU, `tests/cpu/test_injection.py`: the V1-shaped configuration stays below 2^53 over 5e6 steps (fails on
  70202df: `810000000 + 5000000 · 6480048000`); the block holds at least one id per stream, slot and layer cell
  of `face_slabs`; 0 without injection.
- CPU, `tests/cpu/test_checkpoint.py`: a restart with injection accepts the stored `ids_per_step`, refuses a
  missing and a different one; a run without injection reads a checkpoint without the entry (fails without the
  check: "DID NOT RAISE").
- GPU, `tests/gpu/test_engine.py`: faces on three axes of a 6 × 8 × 10 box, one paired, a listed test species,
  20 steps: no id twice, none behind the blocks of the steps taken, every species entered. Checked once that it
  fails without the face offsets (5309 unique of 5425).
- Unchanged and still required: the exact id sums 1 GPU against N GPUs of `tests/gpu/mgpu_open_flow.py` and the
  restart tests (`mgpu_restart.py`).
- `ruff check .`, `scripts/check_kernels.py`, `pytest tests/cpu`, `pytest tests/gpu`, then one full
  `scripts/gpu_gates.py --gpus 8` on the final head; the step time of the open-flow gates is compared with the
  record on fef8bde (diagnostic; the kernel does the same work per entering particle).

## Results

- Full gate on 58166d1 (2026-11-15, 8 GPUs, clean tree): ALL PASSED in `fp64` (32 gates), `mixed` (30) and
  `compact` (28); `pytest tests/gpu` 165 passed. Records: `bench/results/2026-11-15-58166d1ad-gpu-gates.json`,
  `-mixed.json`, `-compact.json`.
- Against the record on fef8bde (main before this change): the growth errors, the implicit energy error (1.9e-10)
  and the Gauss maximum of `fp64` and `mixed` (7.6e-12) are the same to the digits recorded in `PROGRESS.md`.
  Different: the `buneman` growth error in `mixed` (-4.7 %, was -4.2 %) and in `compact` (-5.1 %, was -4.7 %; the
  record on 9666c49 had -5.1 %), with `fp64` at -4.7 % in all three records; the Gauss maximum of `compact`
  (2.5e-13, was 2.1e-13). Both are inside their criteria (10 %, 1e-10). Which of these differences come from this
  change and which from run to run was not measured (no second gate on the same commit).
- Step time (diagnostic): `alfven_inflow:gate` 11.5 / 11.6 / 12.1 ms per step in `fp64` / `mixed` / `compact`
  (fef8bde: 11.9 / 12.0 / 12.0); the wall times of the three `bow_shock` gates are within 0.3 % of fef8bde's.
- The run that showed the shortfall (270000 x 2 x 2 cells, 250 per cell, 5e6 steps) is accepted by the planner on
  this branch: its ids take 8.9e-6 of 2^53.
