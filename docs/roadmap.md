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
gated on 84f5808. Current filter with open faces and the Marder zone at open faces (PIC-D15,
PIC-D16, 2026-11-11): with them a magnetised wind through one paired and four flux faces stays
uniform for three crossings at the bow-shock gate size (null test M2 of the plan, passed). Absorbing
sphere (PIC-D14, 2026-11-11): an insulating body that removes the particles ending a step inside it
and keeps their charge; with the gate dipole its charge levels within three crossings.

Bow shock of a magnetised planet (R6, 2026-11-11): `examples/bow_shock` puts these features
together (wind as external fields, five injecting faces, dipole, absorbing sphere). Its gate
configuration (240×320×320, 3.9e8 particles, 8 GPUs, 264 s) gives a magnetopause at 0.92 to 0.93 of
the pressure-balance distance and a shock at 1.38 (Spreiter et al. 1966: 1.37); gate M3 of the plan.

Blast wave around a magnetised body (R6, 2026-11-11): `examples/blast_wave` follows the debris
clouds of Le et al. 2021 and Keenan et al. 2022 in full PIC: a shell with a radial drift (radial
beams, PIC-D17) expands into a magnetised plasma at rest in a box with six injecting faces, around
an absorbing sphere with a dipole. Gate configuration (384³, 9.4e8 particles, 8 GPUs, 164 s) and
production configuration (576³, 3.1e9 particles, 788 s): front at 2.1 fast magnetosonic speeds,
debris extent 21 to 22 % above the scaling of Le et al. 2021 (a 2-D fit used with D = 3).

Laser and beam drivers and a moving window (R8, 2026-11-12; plan `docs/plans/2026-11-11-wakefield.md`):
a laser antenna (the current of a magnetisation on one layer of cells, PIC-D20), a bunch driver
loaded on a heavy partner of the opposite charge (PIC-D21), lattice places for the load (PIC-D23)
and a box that follows the driver at c one cell at a time (PIC-D22, output in laboratory places
OUT-D5); all explicit scheme only, charge conservation by the deposition alone. `examples/wakefield`
has two gates in the moving window: the wake of an a₀ = 0.1 pulse (amplitude −1.3 % of the linear
theory, wavenumber −0.02 %) and of a rigid bunch (+0.5 %, −0.25 %), equal to the box at rest on the
same interval. One recorded run with an a₀ = 4 pulse (`laser_bubble`, 1536×320², 2.1e9 to 2.4e9
particles, 172 c/ω_pe of plasma, 52 min on 8 GPUs): cavity radius 4.69 c/ω_pe (scaling of Lu et al.
2007: 4), slope of E_x on the axis 0.446 (0.5), trapped electrons up to u_x = 457. Not included:
ionisation, a boosted frame, a window with injection, a body or the implicit scheme, a bare bunch
with a Poisson start.

## Next (in priority order)

| ID | Item | Notes |
|---|---|---|
| R2 | Restart checkpoints (float64 full state) | research output done in PR #2; restart still needed for long runs |
| R3 | Performance of the explicit step | rounds 1–3 (#5–#37): 8-GPU strong scaling 6.4× at 50M particles (3.4 ms/step; 22.0 ms on 1 GPU, 0.44 ns/particle/step), shock 6.3 ms/step on 8 GPUs. The push kernel is 71 % of the 8-GPU step (94 % on 1 GPU) and FP64-bound; what remains changes the arithmetic: a fixed-point deposit (R10), float32 storage (R5), `precision: mixed` by default (PIC-D11). Measured and dropped: tiled deposition (deposit atomics are 0.2 % of the PC samples), CUDA graphs (< 1 % on the local kernels, a captured NCCL exchange does not replay), the migration count inside the data message (1.4 %), paired filter passes on a decomposed axis (< 2 %) |
| R4 | Implicit solver: robustness for ωp·dt > 2 | cost after rounds 2–3 (#20, #22, #28–#32): 488 ms/step on 1 GPU at 50M particles (~22× explicit), 76–78 ms on 8 GPUs; `implicit_push` is 92 % of the 1-GPU step and 73 % of the 8-GPU step, the CG 19 % on 8 GPUs. At the tested ωp·dt the Picard iteration contracts by ≈ 1.8e-4 per sweep, so Anderson acceleration or a plasma preconditioner saves no sweep there; they are open only for larger time steps. Measured and dropped: CUDA graphs, a hot/cold split of `implicit_push` (0.5 % of the kernel), a quadratic Picard start (worse start error) |
| R5 | float32 / mixed-precision particle storage | must keep the gate criteria |
| R6 | More examples | PR #2 set, the bow shock and the blast wave done; later: the bow shock with a southward or oblique IMF and at a larger Dp, the blast wave at a larger R_m / d_i and with a Gaussian cloud, magnetised shocks, beam–plasma with B0, long production runs (shock tail, large reconnection boxes) |
| R10 | Deterministic (order-independent) current deposition option | atomics make long turbulent runs realization-dependent (ANA-D3); needed for bitwise reproducibility |
| R11 | Converged Fermi tail of the pair shock | t >~ 8000 / omega_pe, ~1e9 particles (quasi-2-D); needs R5/R10 performance work or many nodes |
| R9 | Asynchronous output writing | rank 0 writes synchronously. Production two-stream run with the full output set: +4.0 % on 1 GPU (0.89 ms/step), +13–14 % on 8 GPUs (0.45 ms/step). A writer thread on rank 0 was measured at −4.0 / −1.5 / −1.3 % of that run on 8 GPUs (one pair of three above the 2 % criterion) and not adopted; reopen when the write time of a production case is measured above 2 % |
| R7 | Multi-node runs (NCCL over network) | launcher is single-node today |
