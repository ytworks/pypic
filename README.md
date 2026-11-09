# pypic

**English** | [日本語](README.ja.md)

pypic is a GPU-accelerated 3-D, special-relativistic, electromagnetic Particle-in-Cell (PIC) engine.
The hot loops are CUDA C kernels, JIT-compiled through CuPy/NVRTC (no `nvcc` needed); set-up and analysis are in Python.
On several GPUs it uses a 3-D domain decomposition with NCCL, one process per GPU.

## Features

| | Explicit `solver: explicit` (default) | Implicit `solver: implicit` (exact energy conservation) |
|---|---|---|
| Particle push | Relativistic Buneman–Boris (`u = γv`); the magnetic rotation conserves \|u\| exactly | Relativistic implicit midpoint; `v̄ = (u⁰+u¹)/(γ⁰+γ¹)` makes `m(γ¹−γ⁰) = q Δt v̄·E` hold exactly |
| Fields | Yee-lattice FDTD (leapfrog); ∇·B = 0 to machine precision | Crank–Nicolson Maxwell, solved with preconditioned CG |
| Current deposition | Esirkepov; continuity equation, and therefore Gauss's law, holds to machine precision | Chen–Chacón-type deposition, with the orbit split at cell faces. Field gather and current deposit use adjoint weights, so charge and total energy are conserved together |
| Macro-particle shape | Cube (CIC = 1st-order B-spline) | Same |
| Nonlinear solver | — | Outer Picard iteration on the field. The per-particle orbit falls back from Picard → Newton (trust region + line search) → 1-D bracketed root solve |

- **Boundary conditions**: `periodic`, `open` or `reflecting`, chosen per axis and per face.
  - `periodic`: periodic boundary.
  - `open`: fields are absorbed by a CPML (explicit) or a matched conductive layer (implicit). Particles that leave are removed and their energy and charge are recorded.
  - `reflecting`: perfect-conductor wall with specular particle reflection. The orbit is split at the wall for the current deposit, so charge conservation holds.
- **Species**: electrons, ions (mass ratio `m_i/m_e` and charge number `Z` are free parameters), positrons.
- **Mixtures**: `e-i`, `e-p` or `e-i-p`, always generated charge-neutral (`n_e = Z n_i + n_p`). All beams share the same positions, so the initial ρ is exactly 0 and the initial E is 0.
- **Automatic parameter planner**: finds the fastest-growing mode from the linear dispersion relation. From it, it chooses box size, grid, time step, run length, particles per cell, memory and rank layout, and reports the margin against every constraint. Every value can be overridden.
- **Initial distribution**: drifting Maxwell–Jüttner (Zenitani 2015). The random numbers are counter-based, so the initial particles are identical whatever the number of GPUs.
- **Diagnostics**:
  - kinetic energy per species, electric and magnetic field energy, momentum
  - total energy budget, including boundary outflow and dissipation in absorbing layers
  - max\|∇·E−ρ\| and max\|∇·B\|
  - complex Fourier modes of E_x
  - (x, u_x) phase space

Units are normalised: `c = ε0 = 1`, length `c/ω_pe`, time `1/ω_pe`, charge `e`, mass `m_e`, density `n0` (the total electron density).
Temperatures are given in `m_e c²` (`T_e`), or in eV via `T_e_eV` etc.

## Installation

```bash
python3 -m venv .venv && .venv/bin/pip install -e ".[multigpu,test]"
```
You need a CUDA 13 driver and CuPy (`cupy-cuda13x`). Multi-GPU runs also need `nvidia-nccl-cu13`.
Tested on Blackwell (sm_120).

## Quick start: two-stream instability

```bash
# show the parameter plan only
python examples/two_stream/run.py --config examples/two_stream/configs/ei.yaml --plan-only
# run on 1 GPU
python examples/two_stream/run.py --config examples/two_stream/configs/ei.yaml
# run on 8 GPUs
python -m pypic.launch -n 8 examples/two_stream/run.py --config examples/two_stream/configs/ep.yaml
# analyse: growth rate vs linear theory, conservation, phase space -> output/<name>/summary.png
python examples/two_stream/analyze.py output/ep
```

Provided configurations:

| File | Mixture | Set-up |
|---|---|---|
| `configs/ei.yaml` | electron–ion | two counter-streaming electron beams over ions at rest (`mass_ratio` adjustable) |
| `configs/ep.yaml` | electron–positron | neutral counter-streaming pair beams (`mode: symmetric`) |
| `configs/eip.yaml` | electron–ion–positron | positrons carry 30% of the positive charge (`positron_fraction`) |
| `configs/ei_implicit.yaml` | electron–ion | same physics as e-i with the exactly energy-conserving implicit solver |

Override YAML values with `--set key=value`, e.g. `--set u_drift=0.3 --set solver=implicit`.
Any `SimConfig` field can be overridden with `--set overrides.ppc=64`.

### Main parameters (`TwoStreamParams`)

- **Plasma**
  - `mixture`: `e-i` / `e-p` / `e-i-p`
  - `mass_ratio`, `Z`, `positron_fraction`
  - `mode`
    - `electron`: counter-streaming electron beams, positive species at rest
    - `symmetric`: every species split into two neutral counter-streaming beams
    - `buneman`: electrons drift relative to the positive species
  - `u_drift` (γv/c), `beam_fraction`
  - `T_e`, `T_i`, `T_p` (or `*_eV`)
  - `B0`
- **Numerics**
  - `solver`, `n_wavelengths`, `cells_per_wavelength`, `transverse_ratio`
  - `ppc` (positions per cell; each beam gets one particle per position)
  - `cfl`, `wpe_dt_max`, `growth_times` (run length in units of `1/γ_max`)
  - `boundaries`, `max_cells`

### What the planner chooses

1. It solves the multi-component warm-fluid relativistic dispersion relation `1 = Σ ω_ps²/γ_s³ / ((ω − k v_s)² − 3k²v_ts²)` for `k_max`, `γ_max` and `ω_r`.
2. It sets `L_x` to `n_wavelengths` fastest-growing wavelengths and rounds cell counts to multiples of 8.
3. It sets `Δt = min(` CFL (explicit: `cΔt/Δx = 0.95/√3`, implicit: 1)`,` `ω_pe Δt ≤ 0.2` (implicit 0.5)`,` `ω_ce Δt ≤ 0.2``,` `γ_max Δt ≤ 0.05``,` `|ω_r| Δt ≤ 0.2 )`.
4. For the explicit solver with `Δx > π λ_De`, it warns about finite-grid heating and recommends the implicit solver.
5. It estimates the particle count, memory per GPU and rank layout.

## Using it as a library

```python
from pypic import SimConfig, SpeciesConfig, Beam, Simulation
cfg = SimConfig(n_cells=(64, 32, 32), dx=(0.1,) * 3, dt=0.05, n_steps=1000, solver="explicit",
                boundaries=(("open", "open"), ("periodic", "periodic"), ("reflecting", "reflecting")),
                species=[SpeciesConfig("electron", -1, 1, [Beam(1.0, (0.1, 0, 0), 1e-3)]),
                         SpeciesConfig("ion", 1, 1836, [Beam(1.0, (0, 0, 0), 1e-6)])], ppc=16)
sim = Simulation(cfg)      # for several GPUs, start with: python -m pypic.launch -n N
sim.run(progress=100)
sim.save()                 # output/history.npz
```

## Tests

```bash
pytest tests        # the multi-GPU consistency tests also run when 2 or more GPUs are available
```
The tests cover:
- relativistic gyration (\|u\| conservation, rotation angle) and E×B drift
- Gauss's law and ∇·B (explicit, all 3 boundary types)
- implicit total energy (relative error < 1e-11, all 3 boundary types)
- Yee numerical dispersion of a vacuum wave; CPML absorption (residual energy < 1e-4)
- particle handling at reflecting and open boundaries
- charge neutrality of the mixtures
- dispersion solver vs analytic results
- 3-D two-stream growth rate vs theory (within 15%)
- agreement between 1 and 2 GPUs (relative difference < 1e-10)

## Verification (8× RTX PRO 6000 Blackwell)

These are the 3-D two-stream examples run on 8 GPUs with each configuration's default parameters.
γ is fitted to the linear growth phase of the mode at k_max. Because e-i-p has ω_r ≠ 0, its fit includes both counter-propagating growing waves.

| Run | Grid / particles | γ_fit | γ_theory | Error | max \|ΔE_tot\|/E₀ | max\|∇·E−ρ\| | max\|∇·B\| |
|---|---|---|---|---|---|---|---|
| e-i (explicit) | 128×32×32 / 6.3M | 0.351 | 0.342 | +2.6% | 2.2e-3 | 2.1e-12 | 4.6e-15 |
| e-p (explicit) | 128×32×32 / 8.4M | 0.453 | 0.468 | −3.2% | 2.4e-3 | 9.2e-13 | 9.4e-15 |
| e-i-p (explicit) | 128×32×32 / 8.4M | 0.320 | 0.332 | −3.7% | 3.5e-4 | 3.4e-12 | 8.7e-15 |
| e-i (implicit) | 64×16×16 / 0.8M | 0.349 | 0.342 | +2.1% | **1.9e-10** | 3.4e-13 | 1.6e-15 |

- **Multi-GPU agreement**: 1 GPU vs 2 GPUs (2×1×1, 1×2×1, 1×1×2) and 8 GPUs (4×1×2) agree to a relative difference of about 1e-15, for both solvers and all three boundary types.
- **Implicit total energy**: conserved to about 1e-13 with periodic, reflecting and open boundaries (the open-boundary budget includes dissipation in the absorbing layer).
- **Explicit strong scaling** (256×64×64 cells, 50M particles):

  | GPUs | 1 | 2 | 4 | 8 |
  |---|---|---|---|---|
  | Wall time | 11.6 s | 7.0 s | 4.1 s | 2.9 s |
  | Speed-up | 1.0× | 1.7× | 2.8× | 4.0× |

  A single GPU runs at 1.0 ns per particle per step.
- **Implicit cost**: a step costs about 20× an explicit one (5–6 Picard iterations × about 15 CG iterations). In exchange it allows Δx ≫ λ_De and cΔt/Δx ≥ 1.

## Layout

```
pypic/
  kernels/common.cuh     interpolation, RNG, warp reduction
  kernels/fields.cu      Yee updates (CPML), curl, div
  kernels/particles.cu   Boris + Esirkepov (wall splitting), rho, loading (Maxwell-Juttner)
  kernels/implicit.cu    implicit push (orbit-split deposition, Newton / bracketing, energy closure)
  grid.py                Yee grid, ghost exchange, boundary conditions, CPML, field diagnostics
  particles.py           SoA particles, sorting, inter-rank migration
  simulation.py          explicit / implicit main loop, CG, diagnostics, output
  parallel/              NCCL communication, 3-D domain decomposition
  launch.py              multi-GPU launcher
  setup/                 mixtures (species), dispersion relation, automatic parameter planner
examples/two_stream/     run.py, analyze.py, configs/*.yaml
tests/
```

## Notes and limitations

- **Explicit energy**: Boris conserves \|u\| exactly in the magnetic rotation, but the total energy of the system is not exactly conserved (typically about 1e-3). Use the implicit solver for machine-precision total energy. Its error is set by the Picard tolerance (`implicit.picard_tol`), about 1e-10.
- **Implicit non-converged particles**: a particle grazing a node plane or wall may not converge in its individual iteration. In that case \|u¹\| is corrected to match the work done along the deposited path, so total energy and charge conservation still hold. The count is available as `sim.nonconverged`.
- **Explicit open boundaries**: the Poynting-flux term of the energy budget is approximate.
- **Implicit open boundaries**: the matched conductive layer accounts for dissipated energy exactly, but absorbs less well than CPML.
- **`mode: buneman`**: carries a net current, so the k=0 mode oscillates at ω_pe.
- **Precision**: double precision only for now. Planned next: a float32 mode, shared-memory current deposition, and overlapping communication with computation.
