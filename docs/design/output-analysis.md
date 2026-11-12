# Output and analysis design

Status: **Accepted** (2026-11-09, PR #2). Plan: `docs/plans/2026-11-09-output-analysis-examples.md`.
OUT-D5 **Accepted** (2026-11-12, wakefield plan).
Decisions: OUT-D1..D5, ANA-D1..D3.
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
- Field convention (`run.json: field_convention`): `E` and `B` of every diagnostic (dumps, slices, lineouts,
  modes, spectra, fields of the tracks) are the fields the particles feel. With external fields
  (PIC-D12) that is Yee field + external field, while the energies `electric` / `magnetic` of the
  scalar history are those of the Yee fields alone.
- Scalar history, particles crossing open faces (PIC-D13): `lost_kinetic`, `lost_charge`,
  `injected_kinetic`, `injected_charge` and `n_injected` are cumulative sums over the run;
  `fill` is the largest particles / capacity ratio over ranks and species at the record. It
  describes the ranks' arrays, so it is the one history quantity that depends on the decomposition.
- Scalar history, particles absorbed by the body (PIC-D14): `absorbed_kinetic`, `absorbed_charge`
  and `n_absorbed` are cumulative sums over the run. `n_absorbed` counts every species; test
  species carry no charge and no energy into the other two. Particles left out at load are not counted.
Why: openPMD is the community standard for PIC data, so other tools read it; reduced
diagnostics have no standard, so their schema is documented here and in `run.json`.

### OUT-D2 — Collective diagnostics, rank 0 writes, decomposition-independent files
Every diagnostic is computed on all ranks at the same step (start of `step()`, state at time
level n; once more after the last step) and gathered to rank 0 over NCCL (`gather_blocks`,
`gather_columns`: counts first, then data, buffers kept until `groupEnd`). Fields are reduced on
each rank to block sums of the coarse cells it touches (any alignment of ranks, region and
coarsening); particles are gathered and sorted by id. Files are therefore identical (float32
round-off) on any number of GPUs (`tests/gpu/test_output.py`).
Every file operation of rank 0 (preparing the directory, run.json, step writes, flush, close)
shares its outcome with all ranks by one scalar reduction. On a failure every rank raises
`OutputError` together, so no rank is left waiting in a collective (review amendment, PR #2).
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

### OUT-D5 — Moving window: diagnostics are placed in the window, the files give laboratory places (2026-11-12) — Status: Accepted
Plan: `docs/plans/2026-11-11-wakefield.md` (Stage 4). Accepted with the gate records of d523e3b (FP64 and `precision: mixed`; `bench/results/2026-11-12-d523e3bbb-gpu-gates.json`, `-mixed.json`). The window itself: PIC-D22 (`docs/design/engine.md`).
- **Placed in the window:** with `window:` the `region` of a diagnostic, the place of a slice plane and the x range and x bins of a histogram are coordinates of the window: the diagnostic follows the driver. A particle region and a histogram filter compare `x − window_shift · dx`.
- **Written:** `window_shift`, the cells that the window has moved at the output (the record of a step precedes the shift at its end, PIC-D22). Field dumps: `gridGlobalOffset[0] = (region_lo + window_shift) dx` per iteration (openPMD). Reduced series (lineout, modes, slice, kspectrum, histogram): a dataset `window_shift` (int64) with one value per row, only in a run with a window. `scalars.h5`: `window_shift` per row. `run.json`: `window: {start, margin_cells, front_cell, window_shift_end}`, or null. Particle and track positions are laboratory positions as they are stored. The J of an output is the current of the step before it (at t − dt / 2); the shift at the end of that step moves it one plane down with the fields (PIC-D22), so it lies at the laboratory places that the offset gives (review of PR #42: it lay one cell too far up, and was zero after a shift that loads a warm plane).
- **Reader:** `Run.window`; `Series.window_shift` (zeros without a window); `Series.x(comp, i)`, `plane_axes(comp, i)` and `extent(comp, i)` give the laboratory x of the output i, and positions in the window without i. The edges of an x histogram stay those of the window (`edges + window_shift[i] dx` in the laboratory). `FieldSeries.get` and the XDMF index take the offset of the iteration.

Why: what a window run looks at moves with the driver, so a region fixed in the window is the useful one and the arrays keep one shape; one integer per row gives the laboratory place back exactly, where an axis of doubles per row would cost more and round.

Not established: a time-averaged output (`average`, the `auto` cadence of spectral kinds) averages at fixed places of the window over an interval in which the window moves, and its row carries the shift of the output step; the phase of a Fourier mode along x refers to the window. Neither is used by a recorded run.

Measured (2026-11-12, RTX PRO 6000; diagnostic only): the laser run of W2 with a window on 2×2×1 GPUs against the box at rest on 1 GPU, regions at the same laboratory places, outputs at 85 and 200 shifts: dump coordinates within 1.8e-15, dump E_y within 6.0e-15 of its maximum, slice 2.2e-15, lineout of the density 9.5e-15, x histogram 2.2e-16 (relative), the same 245,716 particles of a region with positions within 5.3e-15. With J moved by the shift (0e34252; the same two runs, the two dumps): J_y within 4.4e-13 and 9.4e-13 of its maximum in the window, E_y within 3.0e-15 and 5.0e-15 of the largest E_y of the box at rest.

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
| filamentation:pair_b0 | flow-aligned B0 at half the stabilising field: median growth-rate error ≤ 10 % over m = 4…10 against the relation with Ω (`pypic.theory.streaming.filamentation_growth`) |
| filamentation:pair_b0_stable | B0 at 1.5 times the stabilising field: Σ\|B_m\|² (m = 2…10) at the end ≤ e² times its value at 10 % of the run |
| particle_orbits:* | position error ≤ 2 % of the orbit size and error ratio Δt / (Δt/2) in [3, 5] |
| plasma_expansion:isothermal_hydrogen | self-similar ion velocity within 10 %, E_max within 20 % (median) |
| relativistic_shock:pair_2d | shock speed and compression within 10 % of the Synge-gas jump conditions |
| relativistic_shock:pair_2d_sigma | σ = 0.3, field out of the plane: shock speed and compression within 10 % of the magnetised jump conditions; mean downstream B_z / B1 within 10 % of the compression |
| reconnection:pair_sigma10 | B_x(y) within 5 % of B0 of the Harris sheet before the onset |
| bow_shock:gate | upstream density and u_x within 5 % of the wind over the last crossing; sheath compression ≥ 2; magnetopause within 25 % of the pressure-balance distance; shock upstream of it (definitions in docs/plans/2026-11-10-bow-shock.md) |
| bow_shock:south | the criteria of bow_shock:gate; B_z reverses once between the shock and the body, at x_rev with \|x_rev / R_mp − 1\| ≤ 0.25 |
| bow_shock:oblique | the criteria of bow_shock:gate with the upstream state taken at z − z_c ≥ 1.5 R_mp (definitions of the three in docs/plans/2026-11-11-more-examples.md) |
| blast_wave:gate | front faster than the fast magnetosonic speed and largest ring mean of B_z above 1.5 B0 while the front is between 0.3 and 0.7 R_B; front radii on the axes and the diagonals within 5 %; largest extent of the debris ions within 25 % of Eq. 4 of Le et al. 2021 (D = 3) and reached before the end; ambient densities and \|B\| ahead of the front within 5 % of step 0; div B ≤ 1e-12 (definitions in docs/plans/2026-11-11-blast-wave.md) |
| wakefield:laser_linear | fit of the on-axis E_x behind the pulse, c A_th(x) cos(k x + φ) + b with the linear theory A_th at the place: \|c − 1\| ≤ 5 %, \|k v_g − 1\| ≤ 2 %; E_r / E_x at r = w / 2 within 10 % (definitions in docs/plans/2026-11-11-wakefield.md) |
| wakefield:beam_linear | the same fit over five wavelengths behind the bunch: \|c − 1\| ≤ 5 %, \|k β − 1\| ≤ 2 %; amplitude at r = σ_r over the one on the axis within 10 % of R(σ_r) / R(0) |

Why: agreement with analytic and linear theory is the test policy's observable result.
- Where the model and collisionless PIC differ physically, the tolerance sits above the measured, explained offset. The case is Mora's isothermal electrons; see examples/plasma_expansion/README.md.
- Non-analytic results are reported, not gated: power-law indices and reconnection rates.

### ANA-D3 — Production results are recorded, not gated
- Long or large example runs measure results that are not analytic: Fermi tails, γ_max(t), reconnection power laws and rates. `scripts/record_example.py` records them as `type: example` evidence in `bench/results/`. The record holds the run's code commit, the commit of the analysis code and the GPU count.
- One record per run: the run id from `run.json` is part of the record's file name and index entry. Separate runs of a configuration are kept side by side, and recording a run again replaces only its own record.
- The analysis commit is taken before the analysis code is loaded. If the code changes while the analysis runs, nothing is recorded.
- One definition of provenance (`pypic.provenance`) serves run.json, the gate records and the example records. A commit is "-dirty" when a tracked file differs from it, or an untracked, not ignored file exists, under the code paths (`pypic/`, `examples/`, `scripts/`, `tests/`, `pyproject.toml`). Documentation and recorded evidence are not code.
- Evidence files (gates and examples) are strict JSON: a metric that could not be measured is `null`, never NaN. `verify_bench_results.py` rejects anything else in CI.
- Each record is one realization. The atomic current deposit has no fixed order, so turbulent runs from the same seed diverge at round-off level. The two recorded runs of `pair_sigma10_large` gave p = 2.10 both times, with reconnection rates 0.064 and 0.074. The READMEs quote recorded runs only.

Why: these numbers are research results with statistical spread, not contracts.

## Measured (2026-11-09, `examples/two_stream` ei, 8 GPUs)

Full research output (modes, lineout, phase space ×100, spectra, 6 field dumps, slices, 2000
tracks): 32 MB on disk for 62 MB raw; wall time 18.3 s vs 15.6 s without output.

The wall times are those of the code of that date. On 158bc20 (2026-11-10) the production run
(`cells_per_wavelength=64`, 4936 steps, 50M particles) with the full output set costs 3.80–3.83
ms/step against 3.36–3.37 without output on 8 GPUs (+13–14 %, 0.45 ms/step), and 23.03 against
22.15 on 1 GPU (+4.0 %). A bench run shortened with `t_end` keeps the number of outputs per run and
overstates this cost (`docs/gpu-runbook.md`). A writer thread on rank 0 was measured at −4.0 / −1.5 /
−1.3 % of the 8-GPU run and not adopted: OUT-D2 (synchronous write, errors raised on every rank at
the output) is unchanged.

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
  - **Unrecorded realization**: the READMEs compared two runs each of `pair_sigma10_large` and `pair_2d_long`, but only one run of each had been recorded. Each configuration was run again and recorded (2026-11-10). The READMEs now quote the two recorded runs, and the unrecorded values are dropped.
- 2026-11-10 (PR #2 third review):
  - **Rank-0 file errors**: a failure to create or clean the output directory raised on rank 0 only, and the other ranks went on into the next collective. Every rank-0 file operation of the output now shares its outcome with all ranks (OUT-D2).
  - **Untracked code**: an untracked helper imported by an analysis was outside the provenance check, so the record looked clean and editing the helper during the analysis passed. Untracked code files are now part of the provenance (`pypic.provenance`). This applies to run.json and the gate records too. For tracked changes the check is now limited to the code paths, so an uncommitted evidence index no longer marks the next run dirty.
- 2026-11-10 (PR #2 fourth review): the third-review change marked the scalar history rows as written before the write, so a failed write lost them on retry (a regression). A row now counts as written only once the scalars writer has taken it into its buffer, and only rank 0 tracks this. Whether to write is still decided on every rank from the history length.
  - If creating the file or its datasets fails, no row is taken and the next write retries them.
  - A failed chunk write keeps the rows buffered, and rewriting is repeatable.
- 2026-11-10 (PR #2 fifth review): `SeriesWriter` emptied its buffer before the HDF5 flush. After a failed flush, the next flush() returned without flushing, and SWMR readers did not see the rows. The rows now leave the buffer only after the flush succeeds, and a retry repeats the write and the flush.
