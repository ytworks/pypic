# Engine design

Status: **Accepted** (2026-11-09, PR #1; PIC-D8..D10 PR #2); PIC-D11 **Proposed** (2026-11-10, #10).
Applies to: `pypic/` (engine), `pypic/setup/` (set-up layer), `examples/two_stream/`.

## Requirements (owner, 2026-11-09)

GPU-accelerated, multi-GPU PIC; relativistic Buneman–Boris pusher with numerically guaranteed
energy behavior; field scheme preserving the conservation laws; cubic macro-particles;
periodic / open / reflecting boundaries; everything 3-D; two-stream example; electron, ion
(mass ratio configurable) and positron species in e-i, e-p, e-i-p mixtures, charge neutral by
default; time step and plasma parameters chosen automatically for the experiment.

## Decisions

### PIC-D1 — Python + CuPy RawKernel, NVRTC-compiled CUDA C specialised by `#define`
Hot loops are CUDA C compiled at run time with the local grid as compile-time constants
(`pypic/kernelsrc.py`); orchestration, set-up and analysis are Python.
Why: no nvcc dependency, constant-folded indexing, one source for every grid layout; the same
source builder lets CPU CI compile-check every kernel with NVRTC (`scripts/check_kernels.py`).

### PIC-D2 — Explicit scheme: Yee + relativistic Buneman–Boris + Esirkepov, CIC shape
Leapfrog with B split in two half steps (fields stored at integer times); Esirkepov deposition
for CIC (cubic) particles; staggered CIC gather.
Why: Gauss law and div B to round-off; Boris keeps |u| exact in the magnetic rotation. Total
energy is not exactly conserved (finite-grid effects) — that is what PIC-D3 is for.

### PIC-D3 — Implicit scheme with exact energy and charge conservation (owner choice: both, switchable)
Crank–Nicolson Maxwell (preconditioned CG on `1 + dt^2/4 curl D curl`) + relativistic implicit
midpoint (`vbar = (u0+u1)/(g0+g1)`), orbit split at cell faces with the same path weights for
gather and deposit (Chen–Chacón type). Outer Picard on E; per-particle Picard → Newton (trust
region, line search) → 1-D bracketed solve; per-particle energy closure; B path-averaged across
walls. A step that does not converge is not committed.
Why: the owner asked for numerically guaranteed energy conservation; this is the scheme that
delivers it for the coupled particle–field system, at ~20× the explicit cost per step.

### PIC-D4 — Boundaries per face
Periodic; open = CPML (explicit) / matched conductive layer (implicit, dissipation counted
exactly), particles removed with their energy and charge recorded; reflecting = perfect
conductor + specular reflection with the deposit split at the wall.
Why: wall splitting keeps charge conservation; the implicit layer keeps the energy budget exact.

### PIC-D5 — 3-D Cartesian decomposition, one process per GPU, NCCL
Halo fill/reduce axis by axis (edges and corners by composition), particle migration axis by
axis, launcher `python -m pypic.launch -n N` (spawns ranks, distributes the NCCL id).
The explicit update of B computes one ghost layer from the exchanged ghosts of E instead of
exchanging B (issue #34); its products are rounded one by one, so ranks with and without an
absorbing layer give the same bits.
Why: no MPI dependency on the target machine; results equal to round-off for any rank layout.

### PIC-D6 — Neutral, decomposition-independent loading
All beams of all species are placed on the same positions (counter-based RNG keyed by seed,
global cell and slot); weights follow the densities, so rho = 0 and E = 0 exactly at t = 0.
Drifting Maxwell–Jüttner momenta: exact rejection sampler for theta < 1, Sobol for theta ≥ 1.
Why: neutrality without a Poisson solve; identical initial state on any number of GPUs.

### PIC-D7 — Experiment planner from linear theory
The two-stream planner solves the warm relativistic multi-beam dispersion relation for
k_max / gamma_max / omega_r and derives box, grid, dt (CFL, omega_pe dt, omega_ce dt, growth and
oscillation resolution) and run length, reporting each margin; every value is overridable.
Why: the owner asked that time step and plasma parameters be set optimally for the experiment.

### PIC-D8 — Persistent particle ids (2026-11-09, PR #2)
Particle storage gains an 8th row: an integer-valued float64 id, assigned by the loader as
`((beam - 1) · N_cells + global cell) · ppc + slot` (linear global cell index), so ids are unique and
independent of the rank layout; they move with the particle through sort, migration and removal.
Initialisation refuses runs with `n_beams · N_cells · ppc ≥ 2^53` (exact float64 integers).
Why: tracked orbits and particle samples for acceleration studies (OUT-D1/D2) must follow the
same physical macro-particles on any number of GPUs.

### PIC-D9 — Initial-state extensions (2026-11-09, PR #2)
- **Perturbations** (`Beam.perturbation`):
  - Velocity modes are added directly.
  - Density modes (exact inverse of the cumulative density, by bisection on a bracket that holds for every particle when |α| < 1; the residual is checked on all ranks before anything moves) and displacement modes move the loaded particles virtually. The move goes in rounds of less than one cell; each round's current is deposited with Esirkepov and the particles migrate between rounds (the number of rounds is agreed by all ranks). Then E −= Δt J, so div E = ρ holds to round-off from a neutral start, on any number of GPUs.
- **Explicit half step**: with non-zero initial E or B, a Boris step of −Δt/2 turns the loaded momenta into u^{−1/2}.
- **Test species** (`SpeciesConfig.test`, `particles`):
  - They use separate force and deposit charges (`qdep = 0`), also in the virtual move of their perturbations.
  - They are excluded from the budget and from Gauss's law.
  - They are loaded after the plasma species, so they never change the plasma's random numbers.
  - Explicit particle lists get ids after all beam ids.
- **Initial fields**:
  - Uniform `E0`. It is rejected if it is tangential to a conducting face.
  - Plane waves that are the exact discrete eigenmode of each solver (Yee with the stored-B half-step average, or Crank–Nicolson).
  - 1-D profiles (`pypic.profiles`) of a component along another axis, so the divergence is zero bitwise.
- **Density profiles** (`Beam.profile`): the profile sets the particle weights; shapes below 1e-10 are cut and those particles removed. Equal profiles keep the plasma exactly neutral, and neutrality is checked slot by slot.

Why: the analytic-solution examples (Landau and Dawson, orbits, Mora expansion, Harris sheets) need these initial states without giving up any conservation contract.

### PIC-D10 — Binomial current filter (2026-11-09, PR #2)
- `current_filter: {passes: N, compensate: bool}`: N passes of (1, 2, 1)/4 along every axis, plus an optional (−N/4, 1 + N/2, −N/4) compensator. It is applied after the ghost reduction and before the wall boundary conditions.
- **Wall rules**: node-type data are mirrored about the face node; wall-normal staggered data have a zero ghost; the hi face node is filtered too. With these rules the filter is symmetric and commutes with the discrete divergence.
- The Gauss check filters ρ in the same way.
- Implicit particles gather F·E_h, so the energy identity Σ E_h·(F J) = Σ (F E_h)·J stays exact.
- Along an axis that wraps onto the rank itself (periodic, not decomposed) the passes are applied as one composite stencil with indices modulo the period (issue #37): the same operator, still symmetric, equal to the passes one by one to rounding.
- Periodic and reflecting faces only; open faces are rejected.

Why: the numerical Cherenkov instability of relativistic flows (shocks, reconnection outflows). Verified: transfer function to 1e-12; Gauss ≤ 1e-12 and implicit energy ≤ 1e-11 with walls; 1 vs 2 GPUs with a decomposed wall axis.

### PIC-D11 — Mixed-precision field gather (2026-11-10, #10) — Status: Proposed
- `precision: "fp64" | "mixed"` (default `fp64`). With `mixed`, the explicit push reads E and B from float32 copies made once per step. The six trilinear interpolations run in FP32.
- **Unchanged (FP64):** positions, momenta, the Boris rotation, the Esirkepov current deposit, ρ and the field solver.
- **Why the deposit stays FP64:** Gauss's law to 1e-10 needs the deposited current to satisfy the discrete continuity equation against the FP64 CIC charge. FP32 shape factors break that telescoping sum at ~1e-7 per step.
- **Rank independence:** the gather weights come from global index coordinates `x / dx`, computed in FP64 and rounded once to FP32. The local index is the global cell minus the rank's integer offset (`OXI..`). So the weights do not depend on the rank layout, and 1 vs N GPUs stays equal to round-off.
- **Scope:** explicit solver only. The implicit energy identity needs the gather and the deposit on the same FP64 path weights, so `mixed` with `implicit` is rejected.

Why: `push_deposit` is FP64-pipe bound on SM120, where FP64 runs at 1/64 of FP32. The gather is the part of the push that can leave FP64 without touching a conservation contract.

Measured (2026-11-10, 8× RTX PRO 6000; diagnostic only):
- **Speed:** the push without deposit falls from 18.4 to 14.2 ms per step. The 1-GPU two-stream step falls from 42.7 to 38.6 ms, the 8-GPU two-stream step from 9.4 to 9.1 ms, and the 8-GPU shock step from 14.4 to 13.4 ms.
- **Gates:** all 19 explicit gate configurations pass with `precision: mixed`: two-stream ei/ep/eip/buneman and the 15 analytic-example gates, including the particle-orbit Δt convergence. Growth errors and Gauss (≤ 7.7e-12) match FP64.
- **Tests:** 1 GPU and a 2×2 rank grid agree to 1e-10. A small plasma differs from FP64 by ~2e-11 in kinetic energy.

Status stays Proposed: the default remains `fp64` until the owner accepts `mixed` as a default.


## Review amendments

- 2026-11-09 (PR #1 review): CG guards against 0/0 after exact convergence; exact MJ sampling;
  no commit of non-converged implicit steps; ghost width validated against cells per rank;
  nested planner overrides converted through `SimConfig.from_dict`; launcher polls all ranks.
- 2026-11-09 (PR #1 second review): implicit orbit folding chooses the wall from the unfolded
  midpoint (charge conservation for |d| > L/2); typed `SpeciesConfig` overrides accepted.
- 2026-11-09 (PR #2): implicit-step error flags are reduced over ranks before raising (a
  rank-local raise left the other ranks in a collective); a particle exactly on the upper
  reflecting wall of a decomposed axis stays with the edge rank (migration used to drop it).
  The in-memory mode / phase-space diagnostics and `Simulation.save()` were replaced by the
  output system (`docs/design/output-analysis.md`).
- 2026-11-09 (PR #2 review):
  - The density-perturbation inverse used Newton's method, which did not converge where the density nearly vanishes. With α = 0.99, 170 of 8192 particles were left, and the virtual move then asked for ~1e16 rounds. It is now a bisection with a residual check reduced over ranks.
  - The virtual move deposited a test species' real charge. It now uses `qdep`.
- 2026-11-10 (#6, PIC-D2 implementation): the explicit Esirkepov deposit iterates over the 3-node
  CIC support of each axis (from min(floor X0, floor X1)) and the 2 non-zero prefix sums of each
  current component, instead of 4 nodes and 3 prefix terms with zero tests; the push shares the
  node / half-cell gather weights of each axis. Same scheme and per-term arithmetic (round-off
  order only), fewer FP64 operations: the kernel is FP64-pipe bound on SM120 (FP64 = 1/64 FP32).
