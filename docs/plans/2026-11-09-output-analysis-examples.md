# Research output, analysis/visualisation tools and analytic-solution examples
Status: In progress (PR #2) — Phase 1 done; Phases 2–4 planned
Refs: roadmap R2 (output part), R6; OUT-D1..D4, ANA-D1, PIC-D8..D10

## Context

The engine kept only an in-memory energy history, Ex modes and (x, u_x) histograms and wrote
them to `history.npz` from rank 0. The owner asked (2026-11-09) for: file output that can be
analysed later as research data; output cadence derived from the configuration at the interval the
research needs; storage-light files (lossless by default, lossy optional); general analysis and
visualisation tools for plasma physics and particle acceleration, including animations; and one
example per classic problem with an analytic or well-known solution. The owner chose all eight
proposed problems: plasma-wave dispersion, Landau damping + Dawson oscillation / wave breaking,
filamentation (Weibel), Buneman, test-particle orbits, plasma expansion into vacuum (Mora),
relativistic pair shock (Fermi acceleration), relativistic reconnection.

## Approach

Phase 1 — output foundation (done):
- `pypic/output/`: cadence resolver (`schedule.py`, pure Python, used by planners and engine),
  HDF5/openPMD writers (`h5.py`), GPU diagnostics (`diagnostics.py`: fields, slice, lineout, modes,
  kspectrum, histogram, particles, tracks), manager (collective compute at the start of each step,
  rank 0 writes). `SimConfig.output`; `Simulation.save()`, `mode_diag`, `phase_space_*` removed.
- Particle ids (PIC-D8); `gather_blocks` / `gather_columns`; kernels `deposit_moment`,
  `gather_at_particles`.
- `pypic/analysis/`: reader, spectral (omega-k, ridge, growth and damping fits), particles
  (spectra, power laws, Maxwell–Jüttner fits, work along tracks), fields, plot, animate, xdmf, CLI.
- `pypic/setup/runner.py` (shared example CLI); two-stream ported (planner emits diagnostics and
  time scales), `buneman.yaml` added.
- Bug fixes: implicit error flags reduced over ranks; particle on the upper reflecting wall kept.

Phase 2 — waves and kinetics: `pypic/theory/` (Z function, kinetic Langmuir, Bohm–Gross, EM with
Yee correction, ion acoustic, Stix cold magnetised branches, Bernstein, streaming instabilities
incl. filamentation); Gauss-exact initial perturbation (virtual Esirkepov displacement in rounds
with migration; explicit backward half kick); examples `plasma_waves`, `landau_damping`,
`filamentation`.

Phase 3 — acceleration I: test species (`q`/`qdep` split, excluded from budgets, no RNG
streams), explicit particle lists, uniform `E0`, Yee/CN-exact plane waves, 1-D field and density
profiles; examples `particle_orbits`, `plasma_expansion`.

Phase 4 — acceleration II: binomial current filter (+ compensator) consistent with Gauss's law at
walls and with implicit energy (gather from F·E); examples `relativistic_shock`, `reconnection`.

## Invariants

All collectives (diagnostics and gathers) are reached by every rank; only writing is rank-0.
Output must not change the physics: Gauss / div B / implicit energy / 1-vs-N-GPU gates unchanged.
New engine features (Phases 2–4) each get a Gauss ≤ 1e-12 and implicit-energy ≤ 1e-11 GPU test.

## Verification

- CPU: schedule forms and budget, writer/reader round trip + `openPMD_check_h5`, mantissa rounding
  bound, omega-k / damping / power-law fits on synthetic data, GIF movie.
- GPU: files equal the in-memory state; every output file equal on 1 and 2 GPUs; the two
  regression tests (each fails on the pre-fix code, checked once).
- Gates: `python scripts/gpu_gates.py --gpus 8` per phase, plus one analytic gate per example.
