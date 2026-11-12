# pypic

**English** | [日本語](README.ja.md) | [Wiki: examples with animations and comparison charts](https://github.com/ytworks/pypic/wiki)

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
- **Diagnostics**: kinetic energy per species, field energies, momentum, the total energy budget (including boundary outflow and absorbing-layer dissipation), max\|∇·E−ρ\| and max\|∇·B\|.
- **Research output** (HDF5, see [Output and analysis](#output-and-analysis)): openPMD field and particle dumps, plus light reduced diagnostics — lineouts, slices, Fourier modes, k spectra, phase-space and energy histograms, tracked particle orbits. The output cadence is derived from the physics of the run, and the files are identical whatever the number of GPUs.

Units are normalised: `c = ε0 = 1`, length `c/ω_pe`, time `1/ω_pe`, charge `e`, mass `m_e`, density `n0` (the total electron density).
Temperatures are given in `m_e c²` (`T_e`), or in eV via `T_e_eV` etc.

## Installation

```bash
python3 -m venv .venv
.venv/bin/pip install -e ".[gpu,multigpu,dev,viz]"   # GPU machine (CuPy, NCCL, pytest, ruff, MP4 movies)
.venv/bin/pip install -e ".[dev,kernelcheck,viz]"    # CPU only: set-up, analysis, tests/cpu, kernel check
```
The engine needs a CUDA 13 driver and CuPy (`cupy-cuda13x`, extra `gpu`); multi-GPU runs need
`nvidia-nccl-cu13` (extra `multigpu`). The set-up layer (`pypic.config`, `pypic.setup`) is pure
Python and runs anywhere. Tested on Blackwell (sm_120).

## Quick start: two-stream instability

```bash
# show the parameter plan only
python examples/two_stream/run.py --config examples/two_stream/configs/ei.yaml --plan-only
# run on 1 GPU
python examples/two_stream/run.py --config examples/two_stream/configs/ei.yaml
# run on 8 GPUs
python -m pypic.launch -n 8 examples/two_stream/run.py --config examples/two_stream/configs/ep.yaml
# analyse: growth rate vs linear theory, conservation, phase space -> output/two_stream/ep/summary.png
python examples/two_stream/analyze.py output/two_stream/ep --movie   # + two_stream.mp4
```

`pypic.launch` starts one process per GPU and binds each to the CPUs of the NUMA node its GPU is
attached to; `--no-bind` leaves the CPU affinity alone.

Provided configurations:

| File | Mixture | Set-up |
|---|---|---|
| `configs/ei.yaml` | electron–ion | two counter-streaming electron beams over ions at rest (`mass_ratio` adjustable) |
| `configs/ep.yaml` | electron–positron | neutral counter-streaming pair beams (`mode: symmetric`) |
| `configs/eip.yaml` | electron–ion–positron | positrons carry 30% of the positive charge (`positron_fraction`) |
| `configs/ei_implicit.yaml` | electron–ion | same physics as e-i with the exactly energy-conserving implicit solver |
| `configs/buneman.yaml` | electron–ion | Buneman instability: all electrons drift through ions at rest, γ ≈ (√3/2)(m_e/2m_i)^{1/3} ω_pe |

Override YAML values with `--set key=value`, e.g. `--set u_drift=0.3 --set solver=implicit`.
Any `SimConfig` field can be overridden with `--set overrides.ppc=64`, and any output setting with
`--set output.budget_GB=5` or `--set output.codec=deflate`.

A run writes to `output/<example>/<config>/` by default, or to the directory given with `--out DIR`.
A directory that already holds a run is refused, so the files of two runs never mix. `--overwrite`
replaces that run: its output files are removed first, and nothing else in the directory is touched.

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

## Examples with analytic or known solutions

Each example has the same structure:
- `examples/<name>/` holds a planner (physics → grid, time step, run length, output cadence), configs, `theory.py` (the solution) and `analyze.py` (verdict JSON, figures, `--movie`).
- The gate configs run in `scripts/gpu_gates.py`, with the criteria given in each README.

| example | physics | compared with | measured (gate config) |
|---|---|---|---|
| [`two_stream`](examples/two_stream) | electrostatic two-stream (e-i, e-p, e-i-p, implicit), Buneman | warm-fluid multi-beam dispersion | growth −2.3 … +3.5 %, Buneman −4.7 % |
| [`plasma_waves`](examples/plasma_waves) | ω–k of thermal noise: Langmuir, EM, X/O, Bernstein, whistler, ion acoustic | kinetic / cold magnetised dispersion | Langmuir 1.4 %, EM 0.3 %, X/O 0.07 %, Bernstein 0.5–0.7 % (medians) |
| [`landau_damping`](examples/landau_damping) | Landau damping, Dawson oscillation, wave breaking, relativistic oscillation | Z-function root, exact sheet solutions | γ +4.6 %, ω −0.1 %; Dawson 0.3 %; crossing time −1.2 %; relativistic period 0.08 % |
| [`filamentation`](examples/filamentation) | current filamentation (Weibel) | cold relativistic dispersion | growth 2.2 % (median over modes) |
| [`particle_orbits`](examples/particle_orbits) | gyration, E×B, hyperbolic, E > B, plane wave (figure 8) | exact relativistic orbits | errors ≤ 0.3 %, second-order convergence (ratio 3.76–4.00) |
| [`plasma_expansion`](examples/plasma_expansion) | TNSA-type ion acceleration into vacuum | Mora 2003 self-similar solution | velocity profile −8 %, E_max −15 % (kinetic electrons, see README) |
| [`relativistic_shock`](examples/relativistic_shock) | γ₀ = 15 pair shock, Fermi acceleration | relativistic Rankine–Hugoniot (Synge gas) | β_sh −5.4 %, n₂/n₁ +2.2 % |
| [`reconnection`](examples/reconnection) | σ = 10 pair reconnection, plasmoids, power law | Harris equilibrium; published rate and index (empirical) | equilibrium 1.1 % of B₀; rate ≈ 0.1–0.2; p = 2.4 (small box) |
| [`bow_shock`](examples/bow_shock) | 3-D bow shock of a magnetised planet in the solar wind (open box, dipole, absorbing sphere) | pressure-balance magnetopause; Spreiter et al. (1966) stand-off; perpendicular MHD shock (reported) | magnetopause 0.92 to 0.93 R_mp; shock 1.38 R_mp (1.37); upstream within 0.1 % |
| [`blast_wave`](examples/blast_wave) | 3-D blast wave: a debris shell expanding from a magnetised body into a magnetised plasma at rest (open box, radial beams, dipole, absorbing sphere) | fast magnetosonic speed; debris extent of Le et al. (2021), Eq. 4 with D = 3; perpendicular MHD shock (reported) | front at 2.1 fast speeds, B_z up to 2.7 B₀; debris extent 1.24 R_m (1.02); plasma ahead within 0.05 % |

New set-up features these examples use (all engine-level; see `docs/design/engine.md` PIC-D9, PIC-D10):

- **Initial perturbations** (`Beam.perturbation`): density, displacement or velocity modes. They are applied by a charge-conserving virtual move, so Gauss's law holds exactly.
- **Density profiles** (`Beam.profile`).
- **Test species with explicit particle lists** (`SpeciesConfig.test`, `particles`).
- **Initial fields**: uniform `E0`, and `initial_fields` for exact plane waves and 1-D profiles.
- **Binomial current filter** (`current_filter`), against numerical Cherenkov noise in relativistic flows.
- **`remove_mean_current`** for net-current set-ups such as Buneman.

Open-flow features that `bow_shock` and `blast_wave` use (explicit scheme only; PIC-D12 to PIC-D16):

- **External fields** (`external_fields`): a uniform E and B and regularised dipoles that the particles feel. The grid holds the difference to them, so open faces absorb only the difference; outputs hold the total fields.
- **Injection through open faces** (`injection`): `flux` and `paired` faces keep a drifting plasma entering; the charge is conserved by the deposit of the entering path.
- **Absorbing sphere** (`body`): removes the particles that enter it and keeps their charge.
- **Current filter with open faces, and the Marder zone** (`marder_zone`) that removes the charge left at the open faces, outside the region where Gauss's law is checked.
- **`capacity_factor`** sizes the particle arrays for a run that gains particles; the history key `fill` reports how full they are.

`blast_wave` also uses **radial beams** (`Beam.origin`, `Beam.drift_radial`, profile axis `"r"`; PIC-D17): a density shape of the distance from a point and a drift along the radius, loaded only into the cells the shape reaches. Both schemes load them; only the explicit one was run with them.

## Using it as a library

```python
from pypic import SimConfig, SpeciesConfig, Beam, Simulation
cfg = SimConfig(n_cells=(64, 32, 32), dx=(0.1,) * 3, dt=0.05, n_steps=1000, solver="explicit",
                boundaries=(("open", "open"), ("periodic", "periodic"), ("reflecting", "reflecting")),
                species=[SpeciesConfig("electron", -1, 1, [Beam(1.0, (0.1, 0, 0), 1e-3)]),
                         SpeciesConfig("ion", 1, 1836, [Beam(1.0, (0, 0, 0), 1e-6)])], ppc=16)
cfg.output = OutputConfig(dir="output/my_run", timescales={"plasma_period": 6.283}, diagnostics=[
    {"name": "lineout", "kind": "lineout", "fields": ["Ex"], "axis": "x"},          # auto cadence
    {"name": "phase", "kind": "histogram", "quantities": ["x", "ux"], "bins": [256, 128], "every": {"count": 100}},
    {"name": "fields", "kind": "fields", "fields": ["E", "B", "n:electron"], "coarsen": [2, 2, 2], "every": {"t": 10.0}},
])
sim = Simulation(cfg)      # for several GPUs, start with: python -m pypic.launch -n N
sim.run(progress=100)
sim.close()                # flush and close the output files
```
(`from pypic.config import OutputConfig`.) Without `output.dir` nothing is written; the energy
history stays in `sim.history`.

## Output and analysis

Design: [`docs/design/output-analysis.md`](docs/design/output-analysis.md).

```
output/<run>/run.json                      configuration, schedule, units, time centring, theory (meta)
            scalars.h5                     energy / momentum / conservation history
            reduced/<diag>.h5              lineout, slice, modes, kspectrum, histogram (time series)
            fields/<diag>/<diag>_%08d.h5   openPMD meshes (E, B, J, rho, species moments)
            particles/<diag>/<diag>_%08d.h5   openPMD particles (sorted by id)
            tracks/<diag>_<species>.h5     tracked orbits (x, u, and E, B at the particle)
```

| kind | content |
|---|---|
| `fields` | 3-D fields (`E`, `B`, `J`, `rho`, `n:<sp>`, `flux:<sp>`, `ekin:<sp>`), block-averaged with `coarsen`, cropped with `region` |
| `slice` | a 2-D plane of fields over time |
| `lineout` | 1-D profiles averaged over the other two axes (x–t and ω–k diagrams) |
| `modes` | complex Fourier amplitudes of chosen `(mx, my, mz)` (growth / damping rates in 3-D) |
| `kspectrum` | shell-averaged \|F(k)\|² |
| `histogram` | 1–3-D weighted histograms of x y z ux uy uz u gamma ekin (energy spectra, phase space) |
| `particles` | particle samples (`stride` on the id, `select: {gamma_min: ...}`, `region`) |
| `tracks` | a fixed set of particles (`stride`), optionally with E and B at the particle |

- **Cadence from physics** — `every` takes one of these forms:
  - `25` (steps), `{t: 0.5}`, `{count: 100}`
  - `{per: growth_time, samples: 10}`, `{nyquist: shortest_period, oversample: 2}`
  - `auto`: applies every requirement of the diagnostic kind whose time scale the planner set, and takes the strictest. For spectral diagnostics this means the Nyquist rate of the fastest oscillation (time-averaged over the interval, as an anti-aliasing filter) and at least 10 samples per growth time.

  The planner report prints the resulting interval, number of outputs and size of every diagnostic, and warns about aliasing. `budget_GB` caps the total size by widening only the `auto` cadences.
- **Lightweight files**:
  - The default is float32, with byte-shuffle and zstd lossless compression (`hdf5plugin`).
  - `codec: deflate` is readable by any HDF5 tool.
  - `keep_bits: n` is lossy: it rounds the mantissa, with a relative error of at most 2^-(n+1).
  - Most of the saving comes from writing reduced data instead of full dumps.
- **Analysis** (`pypic.analysis`, CPU only):

  ```python
  from pypic.analysis import open_run, Movie
  run = open_run("output/two_stream/ei"); print(run.summary())       # cadence, sizes, compression
  lo = run.series("lineout"); lo["Ex"], lo.t, lo.x("Ex")             # (n_t, n_x) + coordinates
  f = run.fields("fields"); fd = f.get(f.iterations[-1], "Ez")       # data + coords + time
  Movie(run).phase("phase", "electron").xt("lineout", "Ex").history().save("movie.mp4")
  ```

  Modules:
  - `spectral`: ω–k diagrams, ridge extraction, growth and damping fits.
  - `particles`: spectra, power-law fits, Maxwell–Jüttner temperature, cut-off energy, work by each E component, parallel and perpendicular work along tracks.
  - `fields`: flux function, fronts.
  - `plot`, `animate`, `xdmf`: figures, movies, and ParaView/VisIt 3-D viewing.

  The same is available on the command line:

  ```bash
  python -m pypic.analysis info  output/two_stream/ei
  python -m pypic.analysis plot  output/two_stream/ei
  python -m pypic.analysis movie output/two_stream/ei --panel phase:phase:electron --panel xt:lineout:Ex -o ei.mp4
  python -m pypic.analysis xdmf  output/two_stream/ei --diag fields
  ```

## Tests and CI

```bash
pytest                       # tests/cpu always; tests/gpu when CuPy and a CUDA device are present
python scripts/check_kernels.py   # NVRTC compile check of every kernel (no GPU needed)
python scripts/gpu_gates.py --gpus 8   # GPU verification gates, recorded in bench/results/
```
- `tests/cpu` (no GPU) covers:
  - charge neutrality of the mixtures
  - the dispersion solver vs analytic results
  - nested planner overrides
  - the launcher's handling of a failed rank
- `tests/gpu` covers:
  - relativistic gyration (\|u\| conservation, rotation angle) and E×B drift
  - Gauss's law and ∇·B (explicit, all 3 boundary types)
  - implicit total energy (relative error < 1e-11, all 3 boundary types)
  - Yee numerical dispersion; CPML absorption
  - particle handling at reflecting and open boundaries
  - Maxwell–Jüttner moments vs analytic values
  - CG, Picard non-convergence and ghost-width regressions
  - two-stream growth rate vs theory
  - 1-vs-2-GPU agreement
  - the launcher's CPU binding of a rank
- **CI** (GitHub Actions, `.github/workflows/ci.yml`, CPU only) runs:
  - ruff
  - the `PROGRESS.md` budget check
  - validation of the recorded GPU evidence (`bench/results/index.json`)
  - the NVRTC kernel compile check (sm_80, sm_120)
  - `tests/cpu`

  Work that executes kernels is verified on a GPU machine by `scripts/gpu_gates.py`
  ([GPU runbook](docs/gpu-runbook.md)).

## Development harness

The repository is set up for continuous development across sessions, by people and coding agents:

| File | Role |
|---|---|
| `CLAUDE.md` / `AGENTS.md` | Entry point for Claude Code / Codex: session start, test policy, where things live |
| `PROGRESS.md` | Cross-session memory: Product, Current Status, Change Log. A SessionStart hook and CI enforce its size budget; old entries are archived to `docs/progress/archive/` |
| `.claude/rules/` | `progress-log.md` (how to keep PROGRESS.md); `numerics-invariants.md` (conservation contracts, gate criteria, multi-GPU rules) |
| `docs/design/engine.md` | Design decisions with IDs (PIC-D1..D8) and review amendments |
| `docs/design/output-analysis.md` | Output and analysis decisions (OUT-D1..D4, ANA-D1) |
| `docs/roadmap.md`, `docs/plans/` | Next work items (R1..R7) and per-change implementation plans |
| `docs/gpu-runbook.md`, `bench/results/` | GPU gate procedure and the recorded evidence |

## Verification (8× RTX PRO 6000 Blackwell)

These are the 3-D two-stream examples run on 8 GPUs with each configuration's default parameters, as recorded by the GPU gates in `bench/results/2026-11-10-dded54e12-gpu-gates.json`.
γ is fitted to the linear growth phase of the mode at k_max. Because e-i-p has ω_r ≠ 0, its fit includes both counter-propagating growing waves.

| Run | Grid / particles | γ_fit | γ_theory | Error | max \|ΔE_tot\|/E₀ | max\|∇·E−ρ\| | max\|∇·B\| |
|---|---|---|---|---|---|---|---|
| e-i (explicit) | 128×32×32 / 6.3M | 0.352 | 0.342 | +3.0% | 2.2e-3 | 2.1e-12 | 4.5e-15 |
| e-p (explicit) | 128×32×32 / 8.4M | 0.460 | 0.468 | −1.7% | 2.4e-3 | 8.5e-13 | 9.9e-15 |
| e-i-p (explicit) | 128×32×32 / 8.4M | 0.325 | 0.332 | −2.3% | 2.9e-4 | 3.5e-12 | 8.4e-15 |
| e-i (implicit) | 64×16×16 / 0.8M | 0.354 | 0.342 | +3.5% | **1.9e-10** | 5.8e-13 | 1.5e-15 |
| Buneman (explicit) | 128×8×8 / 0.26M | 0.0512 | 0.0537 | −4.7% | 4.4e-3 | 7.5e-12 | 1.5e-14 |

- **Multi-GPU agreement**: 1 GPU vs 2 GPUs (2×1×1, 1×2×1, 1×1×2) and 8 GPUs (4×1×2) agree to a relative difference of about 1e-15, for both solvers and all three boundary types.
- **Implicit total energy**: conserved to about 1e-13 with periodic, reflecting and open boundaries (the open-boundary budget includes dissipation in the absorbing layer).
- **Explicit strong scaling** (256×64×64 cells, 50M particles):

  | GPUs | 1 | 2 | 4 | 8 |
  |---|---|---|---|---|
  | Step time | 22.0 ms | 11.5 ms | 6.0 ms | 3.4 ms |
  | Speed-up | 1.0× | 1.9× | 3.6× | 6.4× |

  A single GPU runs at 0.44 ns per particle per step (0.36 with `precision: mixed`, 0.30 with
  `precision: compact`). Measured with `scripts/bench_step.py` (no output, history every 10
  steps).
- **Implicit cost**: a step costs about 22× an explicit one on one GPU (5–6 Picard iterations × about 15 CG iterations). In exchange it allows Δx ≫ λ_De and cΔt/Δx ≥ 1.

## Layout

```
pypic/
  kernels/common.cuh     interpolation, RNG, warp reduction
  kernels/fields.cu      Yee updates (CPML), curl, div
  kernels/particles.cu   Boris + Esirkepov (wall splitting), rho, loading (Maxwell-Juttner)
  kernels/implicit.cu    implicit push (orbit-split deposition, Newton / bracketing, energy closure)
  kernelsrc.py           kernel source assembly + compile-time constants (pure Python)
  grid.py                Yee grid, ghost exchange, boundary conditions, CPML, field diagnostics
  particles.py           SoA particles, sorting, inter-rank migration
  simulation.py          explicit / implicit main loop, CG, energy / conservation diagnostics
  output/                cadence planning (schedule), HDF5/openPMD writers, GPU diagnostics, manager
  analysis/              reader, spectral / particle / field analysis, plots, movies, XDMF, CLI (CPU only)
  parallel/              NCCL communication, 3-D domain decomposition, gathers to rank 0
  launch.py              multi-GPU launcher
  setup/                 mixtures (species), dispersion relation, planners, example runner
examples/two_stream/     run.py, analyze.py, configs/*.yaml
tests/cpu/, tests/gpu/   CPU-only tests (CI) / GPU tests
scripts/                 check_kernels, gpu_gates, verify_bench_results, check_progress_size
docs/                    design/, roadmap.md, plans/, gpu-runbook.md, progress/archive/
bench/results/           recorded GPU gate evidence (index.json)
```

## Notes and limitations

- **Explicit energy**: Boris conserves \|u\| exactly in the magnetic rotation, but the total energy of the system is not exactly conserved (typically about 1e-3). Use the implicit solver for machine-precision total energy. Its error is set by the Picard tolerance (`implicit.picard_tol`), about 1e-10.
- **Implicit non-converged particles**: a particle grazing a node plane or wall may not converge in its individual iteration. In that case \|u¹\| is corrected to match the work done along the deposited path, so total energy and charge conservation still hold. The count is available as `sim.nonconverged`.
- **Explicit open boundaries**: the Poynting-flux term of the energy budget is approximate.
- **Implicit open boundaries**: the matched conductive layer accounts for dissipated energy exactly, but absorbs less well than CPML.
- **`mode: buneman`**: carries a net current, so the k=0 mode oscillates at ω_pe.
- **Particle ids**: every particle carries an id (8th storage row, +14 % particle memory).
- **Precision** (`precision:` in the configuration; explicit solver): `fp64`, the default, is double precision throughout. `mixed` gathers the fields at the particles from float32 copies (PIC-D11). `compact` does the same and stores a position as a cell and a 32-bit fraction of the cell, 56 bytes per particle slot in place of 64 (PIC-D19); momenta, weights, the Boris push, the current deposit and the field solver stay in double precision in every mode. Float32 momenta were measured and not adopted: the mean of \|u\| of gyrating particles drifts. The implicit solver takes `fp64` only.
