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

More examples (R6, 2026-11-12; plan `docs/plans/2026-11-11-more-examples.md`): the seven items
that were left.
- Bow shock with a southward IMF (`bow_shock/south`, gate): B_z reverses once, 0.93 to 1.00 R_mp
  from the centre in three runs with the same seed; with an oblique IMF (`oblique`, gate; 45° on
  the subsolar line): magnetopause 0.96 to 0.97 R_mp, shock 1.25 to 1.28 R_mp; at twice the
  magnetosphere (`dp10`, R_mp = 10 d_i, 480×640², 3.1e9 particles, one recorded run of 1.0 h):
  magnetopause 0.94 R_mp, shock 1.44 R_mp, compression 3.8.
- Blast wave at R_m = 12 d_i (`rm12`, 576³, 3.2e9 particles, recorded): debris extent 15 % above
  Eq. 4 of Le et al. 2021 (21 % at 8 d_i); as a Gaussian cloud without the body (`cloud`,
  reported): the gate criteria on the front hold (1.8 fast speeds, B_z 2.1 B0).
- Beam–plasma with B0 along the drift (gates): `two_stream/ei_b0` (growth +2.2 %),
  `filamentation/pair_b0` (median error 5.7 %) and `pair_b0_stable` (above the stabilising field:
  the mode energy stays within the e² limit); the growth rate is in `pypic/theory/streaming.py`.
- Magnetised piston shock (`relativistic_shock/pair_2d_sigma`, σ = 0.3, gate): speed, density and
  field jump within 7.2 % of the magnetised jump conditions; it uses a uniform `E0` tangential to
  reflecting faces (PIC-D18).
- Long runs (one recorded run each): `reconnection/pair_sigma10_xl` (800 × 400 c/ω_pe to
  t = 4000): index 2.04, converged in time, against 2.10 in the half-size box; the energy that
  10⁻⁴ of the particles exceed is 17 to 32 % higher. `relativistic_shock/pair_2d_t4000` and
  `pair_2d_tail` (t = 8000, 3.3e8 particles): the slope above the thermal peak stays at 3.9 to 4.0
  and the largest energy rises by 13 % from t = 4000 to 8000 — no extended tail (see R11).

Performance of the explicit step and compact particle storage (R3, R5, 2026-11-12; plan
`docs/plans/2026-11-11-compact-storage.md`; milestone M6 Complete).
- Rounds 1–3 (#5–#37): 8-GPU strong scaling 6.4× at 50M particles, shock 6.3 ms/step on 8 GPUs.
  Measured and dropped there: tiled deposition (deposit atomics are 0.2 % of the PC samples),
  CUDA graphs (< 1 % on the local kernels, a captured NCCL exchange does not replay), the
  migration count inside the data message (1.4 %), paired filter passes on a decomposed axis
  (< 2 %).
- `precision: compact` (PIC-D19; explicit solver): a position is stored as its cell and a 32-bit
  fraction of the cell, 56 bytes per particle slot in place of 64; momenta, weights, the Boris
  push, the current deposit and the field solver stay FP64. The explicit gates pass on one
  commit in `fp64`, `mixed` (PIC-D11) and `compact` under the same criteria. The default stays
  `fp64` (owner decision); both are options.
- Two-stream benchmark (50,331,648 particles, medians of three): 3.55 / 3.17 / 2.84 ms per step
  on 8 GPUs in `fp64` / `mixed` / `compact`, 22.0 / 18.2 / 15.3 ms on 1 GPU (0.44 / 0.36 / 0.30 ns
  per particle per step). With the phases timed, `kernel_push` is 70 / 63 / 60 % of the 8-GPU
  step. In the bow-shock and blast-wave gate configurations `compact` is 10.5 % and 12.0 % below
  `mixed` and 20 % below `fp64`.
- Device memory (`blast_wave/full`, 3.1e9 particles, 8 GPUs): peak 53.72 GB per GPU in `compact`
  against 61.78 in `fp64`, below the planner's figure (56.57 and 66.16), since the particle sort
  has a memory pool of its own (PIC-D13, amended); before that the peaks were 10.5 % and 13.6 %
  above the figure. The planners take the bytes per slot of the mode.
- Measured and dropped: float32 stored momenta (44 bytes per slot; the mean |u| of gyrating
  particles shows a secular trend, up to 26 standard errors after 2.3e6 steps), FP32 Boris
  arithmetic (the rotation changes |u| by a fixed amount per step), emptying the default memory
  pool after every sort (+7.8 % and +8.3 % of the step).
- What remains is the deposit (R10, below).

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
ionisation, a boosted frame, a window with injection, a body, the implicit scheme or the compact
storage (PIC-D19), a bare bunch with a Poisson start.

## Next (in priority order)

| ID | Item | Notes |
|---|---|---|
| R2 | Restart checkpoints (float64 full state) | research output done in PR #2; restart still needed for long runs |
| R4 | Implicit solver: robustness for ωp·dt > 2 | cost after rounds 2–3 (#20, #22, #28–#32): 488 ms/step on 1 GPU at 50M particles (~22× explicit), 76–78 ms on 8 GPUs; `implicit_push` is 92 % of the 1-GPU step and 73 % of the 8-GPU step, the CG 19 % on 8 GPUs. At the tested ωp·dt the Picard iteration contracts by ≈ 1.8e-4 per sweep, so Anderson acceleration or a plasma preconditioner saves no sweep there; they are open only for larger time steps. Measured and dropped: CUDA graphs, a hot/cold split of `implicit_push` (0.5 % of the kernel), a quadratic Picard start (worse start error) |
| R10 | Deterministic (order-independent) current deposition option | atomics make long turbulent runs realization-dependent (ANA-D3); needed for bitwise reproducibility. Its share of the step (R3, final breakdown): the part of the push that is skipped with zero charge — the shape functions and the Esirkepov deposit — is 0.876 ms of the 2.837 ms `compact` step on 8 GPUs (30.9 %) and 6.323 of 15.317 ms on 1 GPU (41.3 %) |
| R11 | Converged Fermi tail of the pair shock | the run to t = 8000 / omega_pe (R6, `pair_2d_tail`: 3.3e8 particles, 62 c/omega_pe across, 8 particles per cell, 4-pass filter) has no extended tail: index 3.9–4.0 from t = 2000 on, largest energy +13 % from t = 4000 to 8000 (t^1/2: +41 %). Time alone is not what is missing; the transverse size, the particles per cell and the filter were not varied and are not separated. A run at ~1e9 particles needs performance work (R10; `precision: compact` of R5 has not been run on this configuration) or many nodes |
| R9 | Asynchronous output writing | rank 0 writes synchronously. Production two-stream run with the full output set: +4.0 % on 1 GPU (0.89 ms/step), +13–14 % on 8 GPUs (0.45 ms/step). A writer thread on rank 0 was measured at −4.0 / −1.5 / −1.3 % of that run on 8 GPUs (one pair of three above the 2 % criterion) and not adopted; reopen when the write time of a production case is measured above 2 % |
| R7 | Multi-node runs (NCCL over network) | launcher is single-node today |
