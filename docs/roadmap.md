# Roadmap

Product: a GPU PIC engine that runs plasma experiments end to end — choose a physical set-up,
get numerically trustworthy (conservation-guaranteed) 3-D results on one to many GPUs, and an
analysis against theory. Decisions: `docs/design/engine.md`.

## Done (2026-11-09)

Explicit and exactly energy-conserving implicit schemes, three boundary types, multi-GPU NCCL
decomposition, e-i / e-p / e-i-p neutral mixtures, two-stream example with planner and analysis,
CPU CI, GPU gates with recorded evidence. Research output (HDF5/openPMD + reduced diagnostics,
cadence from physics), particle ids, CPU analysis/visualisation package with movies, Buneman
config (PR #2, `docs/design/output-analysis.md`).

Analytic-solution examples (PR #2): plasma waves, Landau damping / Dawson / wave breaking,
filamentation, test-particle orbits, plasma expansion (Mora), relativistic shock, reconnection,
with the set-up features they need (initial perturbations, test species, initial fields, density
profiles, current filter; PIC-D9, PIC-D10).

## Next (in priority order)

| ID | Item | Notes |
|---|---|---|
| R1 | Open-boundary particle injection (thermal flux) | in the original plan; open faces currently only remove particles |
| R2 | Restart checkpoints (float64 full state) | research output done in PR #2; restart still needed for long runs |
| R3 | Performance: shared-memory tiled deposition, fewer host syncs, comm/compute overlap | explicit 8-GPU strong scaling is 4.0× at 50M particles |
| R4 | Implicit solver cost: CUDA graphs or fused kernels for CG, better plasma preconditioning | ~20× explicit per step; Picard 5–6 × CG ~15 |
| R5 | float32 / mixed-precision particle storage | must keep the gate criteria |
| R6 | More examples | PR #2 set done; later: magnetised shocks, beam–plasma with B0, long production runs (shock tail, large reconnection boxes) |
| R8 | Laser / beam drivers and a moving window (wakefield acceleration) | needs field injection at a boundary |
| R9 | Asynchronous output writing | rank 0 writes synchronously today (~17 % of wall time for the full two-stream output set) |
| R7 | Multi-node runs (NCCL over network) | launcher is single-node today |
