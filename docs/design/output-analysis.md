# Output and analysis design

Status: **Accepted** (2026-11-09, PR #2). Plan: `docs/plans/2026-11-09-output-analysis-examples.md`.
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
- `run.json`: configuration, resolved schedule, layout of every diagnostic, time centering, units,
  ED-PIC attributes, git commit, rank grid, and experiment metadata (`meta`: planner theory,
  run statistics). The reader (`pypic.analysis.open_run`) relies only on these files.
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

## Measured (2026-11-09, `examples/two_stream` ei, 8 GPUs)

Full research output (modes, lineout, phase space ×100, spectra, 6 field dumps, slices, 2000
tracks): 32 MB on disk for 62 MB raw; wall time 18.3 s vs 15.6 s without output.
