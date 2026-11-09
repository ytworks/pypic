# Engine design

Status: **Accepted** (2026-11-09, PR #1; PIC-D8 PR #2).
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
