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

Open-boundary particle injection (R1, 2026-11-11): `flux` and `paired` injecting faces (PIC-D13),
gated on 84f5808. The null test of a magnetised open flow (M2 of the bow-shock plan) waits for the
filter with open faces and the Marder zone (its Phase 3).

## Next (in priority order)

| ID | Item | Notes |
|---|---|---|
| R2 | Restart checkpoints (float64 full state) | research output done in PR #2; restart still needed for long runs |
| R3 | Performance of the explicit step | rounds 1–3 (#5–#37): 8-GPU strong scaling 6.4× at 50M particles (3.4 ms/step; 22.0 ms on 1 GPU, 0.44 ns/particle/step), shock 6.3 ms/step on 8 GPUs. The push kernel is 71 % of the 8-GPU step (94 % on 1 GPU) and FP64-bound; what remains changes the arithmetic: a fixed-point deposit (R10), float32 storage (R5), `precision: mixed` by default (PIC-D11). Measured and dropped: tiled deposition (deposit atomics are 0.2 % of the PC samples), CUDA graphs (< 1 % on the local kernels, a captured NCCL exchange does not replay), the migration count inside the data message (1.4 %), paired filter passes on a decomposed axis (< 2 %) |
| R4 | Implicit solver: robustness for ωp·dt > 2 | cost after rounds 2–3 (#20, #22, #28–#32): 488 ms/step on 1 GPU at 50M particles (~22× explicit), 76–78 ms on 8 GPUs; `implicit_push` is 92 % of the 1-GPU step and 73 % of the 8-GPU step, the CG 19 % on 8 GPUs. At the tested ωp·dt the Picard iteration contracts by ≈ 1.8e-4 per sweep, so Anderson acceleration or a plasma preconditioner saves no sweep there; they are open only for larger time steps. Measured and dropped: CUDA graphs, a hot/cold split of `implicit_push` (0.5 % of the kernel), a quadratic Picard start (worse start error) |
| R5 | float32 / mixed-precision particle storage | must keep the gate criteria |
| R6 | More examples | PR #2 set done; later: magnetised shocks, beam–plasma with B0, long production runs (shock tail, large reconnection boxes) |
| R8 | Laser / beam drivers and a moving window (wakefield acceleration) | needs field injection at a boundary |
| R10 | Deterministic (order-independent) current deposition option | atomics make long turbulent runs realization-dependent (ANA-D3); needed for bitwise reproducibility |
| R11 | Converged Fermi tail of the pair shock | t >~ 8000 / omega_pe, ~1e9 particles (quasi-2-D); needs R5/R10 performance work or many nodes |
| R9 | Asynchronous output writing | rank 0 writes synchronously. Production two-stream run with the full output set: +4.0 % on 1 GPU (0.89 ms/step), +13–14 % on 8 GPUs (0.45 ms/step). A writer thread on rank 0 was measured at −4.0 / −1.5 / −1.3 % of that run on 8 GPUs (one pair of three above the 2 % criterion) and not adopted; reopen when the write time of a production case is measured above 2 % |
| R7 | Multi-node runs (NCCL over network) | launcher is single-node today |
