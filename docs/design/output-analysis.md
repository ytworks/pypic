# Output and analysis design

Status: **Accepted** (2026-11-09, PR #2). Plan: `docs/plans/2026-11-09-output-analysis-examples.md`.
Decisions: OUT-D1..D4, ANA-D1..D3.
Applies to: `pypic/output/` (schedule, writers, GPU diagnostics, manager), `pypic/analysis/`,
`pypic/setup/runner.py`, every `examples/<name>/`.

## Requirements (owner, 2026-11-09)

Save results so that they can be analysed later as research data; write at the step interval the
research needs, derived from the configuration; keep the files as small as possible (lossless by
default, lossy optional); general analysis/visualisation tools for plasma physics and particle
acceleration, with animations; one example per classic problem with an analytic / known solution.

## Decisions

### OUT-D1 — HDF5: openPMD 1.1 (ED-PIC) dumps + appendable reduced series
- Field and particle dumps: openPMD 1.1 with the ED-PIC extension, file-based
  (`fields/<diag>/<diag>_%08d.h5`, `particles/<diag>/<diag>_%08d.h5`). Readable by
  openPMD-viewer, yt, ParaView/VisIt; validated with `openPMD_check_h5` in CPU tests.
  `unitSI` is SI when `output.n0_m3` is set, else 1 (normalised units, recorded in `run.json`).
  Time levels in `timeOffset`: explicit momenta −dt/2, J −dt/2 (both schemes).
  Particle `weighting` is the physical number per macro-particle (openPMD requires unitSI 1):
  `w · n0 (c/ω_pe)^3`, i.e. `w` itself for normalised runs.
- Reduced diagnostics (lineout, modes, slice, kspectrum, histogram, scalars, tracks): one HDF5 file
  each, datasets `t`, `step` and one per quantity of shape (n_out, …), SWMR (readable while running).
  Rows are buffered and written a whole chunk at a time: rewriting a partly filled compressed chunk
  leaks its old file space (measured: 11 MB instead of 1.3 MB for frequent track rows).
- One run per directory: a directory that already holds run output (`run.json`, `scalars.h5`,
  `reduced/`, `fields/`, `particles/`, `tracks/`) is refused on every rank. `output.overwrite`
  (`--overwrite`) removes exactly those entries first. Otherwise dumps a re-run does not rewrite
  would be read as part of it (review amendment, PR #2).
- `run.json`: configuration, resolved schedule, layout of every diagnostic, time centering, units,
  ED-PIC attributes, git commit, run id and start time (both fixed at start-up), rank grid, and
  experiment metadata (`meta`: planner theory, run statistics). The reader (`pypic.analysis.open_run`) relies only on these files.
Why: openPMD is the community standard for PIC data, so other tools read it; reduced
diagnostics have no standard, so their schema is documented here and in `run.json`.

### OUT-D2 — Collective diagnostics, rank 0 writes, decomposition-independent files
Every diagnostic is computed on all ranks at the same step (start of `step()`, state at time
level n; once more after the last step) and gathered to rank 0 over NCCL (`gather_blocks`,
`gather_columns`: counts first, then data, buffers kept until `groupEnd`). Fields are reduced on
each rank to block sums of the coarse cells it touches (any alignment of ranks, region and
coarsening); particles are gathered and sorted by id. Files are therefore identical (float32
round-off) on any number of GPUs (`tests/gpu/test_output.py`).
Why: no MPI on the target machine (PIC-D5); gathering reduced data is cheap; identical files
make results comparable across runs and GPUs.

### OUT-D3 — Lightweight by construction
float32 by default; byte-shuffle + zstd (hdf5plugin) lossless, `codec: deflate` for tools without
the plugin; optional lossy `keep_bits: n` rounds the mantissa on the GPU (round to nearest even,
relative error ≤ 2^-(n+1), decoded by any HDF5 reader). Most of the saving comes from writing
reduced data: lineouts, slices, Fourier modes, histograms, block-averaged (`coarsen`) or cropped
(`region`) fields, id-strided particle samples and fixed track sets.

### OUT-D4 — Cadence from physics
`every` is a step count, `{t}`, `{count}`, `{per: <timescale>, samples}`, `{nyquist: <timescale>,
oversample}` or `auto`. Planners publish named time scales (`plasma_period`, `shortest_period`,
`growth_time`, …). `auto` applies every requirement of the kind whose time scale is known and takes
the strictest (spectral kinds: Nyquist of the shortest period ×2 *and* 10 samples per growth time,
time-averaged over the interval as an anti-aliasing filter; tracks: 16 per shortest period;
dumps: a fixed number). Aliasing / under-sampling warnings and a storage table are part of the
planner report; `budget_GB` widens only `auto` cadences and refuses fixed ones that do not fit.

### PIC-D8 — Persistent particle ids (see `docs/design/engine.md`)

### ANA-D1 — Analysis package on CPU only
`pypic.analysis` (reader, spectral, particles, fields, plot, animate, xdmf, CLI) needs no CuPy,
so runs are analysed anywhere and tested in CPU CI. General dispersion relations used by several
examples live in `pypic/theory/`; problem-specific exact solutions stay in `examples/<name>/`.
Figures follow one style: categorical series in a fixed validated order, single-hue sequential
maps for magnitudes, a blue–gray–red diverging map for signed fields, legends for every
multi-series plot, theory dashed black. Movies fix colour limits over all frames.

### ANA-D2 — Example gates against theory
Every example's `analyze.py` returns `passed` against its theory; `scripts/gpu_gates.py` runs the
gate configs and also requires Gauss ≤ 1e-10. Criteria (fixed before the gate evidence was recorded):

| gate | criterion |
|---|---|
| plasma_waves | each gated branch found at ≥ half of the k points, ≥ 90 % of them within max(3 %, 2Δω/ω) |
| landau_damping:landau | γ within 10 %, ω within 2 % of the kinetic root |
| landau_damping:dawson / relativistic | sheet frequency within 1 % (amplitude 2 %) / exact relativistic period within 0.5 % |
| landau_damping:wavebreaking | first crossing time within 5 % of arcsin(ω_pe / k v₀) |
| filamentation:pair | median growth-rate error ≤ 10 % over m = 2…10 (common linear window) |
| particle_orbits:* | position error ≤ 2 % of the orbit size and error ratio Δt / (Δt/2) in [3, 5] |
| plasma_expansion:isothermal_hydrogen | self-similar ion velocity within 10 %, E_max within 20 % (median) |
| relativistic_shock:pair_2d | shock speed and compression within 10 % of the Synge-gas jump conditions |
| reconnection:pair_sigma10 | B_x(y) within 5 % of B0 of the Harris sheet before the onset |

Why: agreement with analytic and linear theory is the test policy's observable result.
- Where the model and collisionless PIC differ physically, the tolerance sits above the measured, explained offset. The case is Mora's isothermal electrons; see examples/plasma_expansion/README.md.
- Non-analytic results are reported, not gated: power-law indices and reconnection rates.

### ANA-D3 — Production results are recorded, not gated
- Long or large example runs measure results that are not analytic: Fermi tails, γ_max(t), reconnection power laws and rates. `scripts/record_example.py` records them as `type: example` evidence in `bench/results/`. The record holds the run's code commit, the commit of the analysis code and the GPU count.
- One record per run: the run id from `run.json` is part of the record's file name and index entry. Separate runs of a configuration are kept side by side, and recording a run again replaces only its own record.
- The analysis commit is taken before the analysis code is loaded. If HEAD or the code under `examples/` and `pypic/` changes while the analysis runs, nothing is recorded.
- Evidence files (gates and examples) are strict JSON: a metric that could not be measured is `null`, never NaN. `verify_bench_results.py` rejects anything else in CI.
- Each record is one realization. The atomic current deposit has no fixed order, so turbulent runs from the same seed diverge at round-off level. Two runs of `pair_sigma10_large` gave p = 2.15 and 2.10, and the READMEs quote both.

Why: these numbers are research results with statistical spread, not contracts.

## Measured (2026-11-09, `examples/two_stream` ei, 8 GPUs)

Full research output (modes, lineout, phase space ×100, spectra, 6 field dumps, slices, 2000
tracks): 32 MB on disk for 62 MB raw; wall time 18.3 s vs 15.6 s without output.

## Review amendments

- 2026-11-09 (PR #2 review):
  - **Run directory**: one run per directory (OUT-D1). An occupied directory is refused unless `output.overwrite` is set.
  - **Empty tracks**: a series with an empty dimension, such as a tracks diagnostic that selects no particle of a species, is written with an unlimited chunked dimension.
  - **Runner overrides**: `output.*` values from the YAML file or `--set` are validated through `OutputConfig`, so diagnostics given as dicts become typed objects.
  - **`fit_two_wave`**:
    - The relative weights 1/|z| are computed in float64 and floored at 1e-6 of the largest |z|. Exact zeros, as in complex64 data, kept infinite weights before.
    - a and b are returned for t = 0 rather than for the window start, so the documented model reproduces the signal.
    - An envelope-based weight was tried and rejected: it moved the Buneman growth fit from −6.7 % to −13.4 %.
  - **`parallel_perpendicular_work`**: computed in float64. Where B = 0 at the particle there is no parallel direction, so E counts as perpendicular there and W_par + W_perp equals the total work.
  - **`record_example.py`**: accepts only the name of an example in this repository and loads that example's `analyze.py` by its path.
- 2026-11-09 (PR #2 second review), all on recorded results (ANA-D3):
  - **Separate runs**: two runs of one configuration recorded on one day overwrote each other. Records now carry the run id.
  - **Analysis commit**: it was read after the analysis, so a checkout during the analysis was recorded as the code that produced the metrics. It is now taken before the analysis and checked after.
  - **Strict JSON**: unmeasured metrics were written as NaN, which is not JSON. They are now null. The five existing evidence files that held NaN were normalized (NaN → null, nothing else changed).
