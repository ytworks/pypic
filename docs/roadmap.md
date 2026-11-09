# Roadmap

Product: a GPU PIC engine that runs plasma experiments end to end — choose a physical set-up,
get numerically trustworthy (conservation-guaranteed) 3-D results on one to many GPUs, and an
analysis against theory. Decisions: `docs/design/engine.md`.

## Done (2026-11-09)

Explicit and exactly energy-conserving implicit schemes, three boundary types, multi-GPU NCCL
decomposition, e-i / e-p / e-i-p neutral mixtures, two-stream example with planner and analysis,
CPU CI, GPU gates with recorded evidence.

## Next (in priority order)

| ID | Item | Notes |
|---|---|---|
| R1 | Open-boundary particle injection (thermal flux) | in the original plan; open faces currently only remove particles |
| R2 | Field / particle output and restart (HDF5 or npz checkpoints) | needed for long runs; was dropped from v0.1 |
| R3 | Performance: shared-memory tiled deposition, fewer host syncs, comm/compute overlap | explicit 8-GPU strong scaling is 4.0× at 50M particles |
| R4 | Implicit solver cost: CUDA graphs or fused kernels for CG, better plasma preconditioning | ~20× explicit per step; Picard 5–6 × CG ~15 |
| R5 | float32 / mixed-precision particle storage | must keep the gate criteria |
| R6 | More examples: Buneman, filamentation/Weibel (transverse), beam–plasma with B0 | each in `examples/<name>/` |
| R7 | Multi-node runs (NCCL over network) | launcher is single-node today |
