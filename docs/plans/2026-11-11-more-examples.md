# More examples: the seven open items of R6

Status: In progress
Refs: roadmap R6 (R11 stays open); PIC-D18 (new); ANA-D2, ANA-D3; PIC-D12..D17 (used as they are)

## Context

Owner request (2026-11-11): close the roadmap rows R3, R5, R6 and R8, each independent piece of
work in its own pull request. This plan is the R6 piece (branch `examples/r6`, one pull request;
the stages below are its commits). The other two pieces have their own plans
(`2026-11-11-compact-storage.md`, `2026-11-11-wakefield.md`) and do not touch this file.

R6 names seven items that are not done: the bow shock with a southward or oblique IMF and at a
larger Dp; the blast wave at a larger R_m / d_i and with a Gaussian cloud; magnetised shocks;
beam–plasma with B0; long production runs (shock tail, large reconnection boxes). One of them
needs a change of the engine: a magnetised piston shock needs an initial electric field tangential
to its reflecting wall, which `Simulation._set_initial_fields` refuses today.

Owner decision (2026-11-11): the shock-tail run goes up to 3.3e8 particles, once. The converged
tail at ~1e9 particles stays R11, after R2.

Defaults chosen without a question, stated to the owner with the plan and not refused: the larger
Dp of the bow shock is 10 (8 if the memory exceeds the planner's bound); the Gaussian cloud has no
body and is reported, not gated.

## Approach

All criteria below are fixed here, before the first run of each configuration. A criterion that
fails is recorded as failed and reported to the owner; a limit or a definition changes only by an
amendment written into this file with its reason (stop rule).

Every formula that is new here gets a CPU test of its limits in `tests/cpu/test_theory.py` and is
compared with the cited paper before the first run; the outcome of that comparison is recorded in
the section of the item.

| Roadmap item | Configuration | Gate / record |
|---|---|---|
| Bow shock, southward IMF | `bow_shock/south` (gate size) | gate |
| Bow shock, oblique IMF | `bow_shock/oblique` (`imf: [1, 0, 1]`) | gate if O1 passes, reported otherwise |
| Bow shock, larger Dp | `bow_shock/dp10` (480×640², 3.15e9 particles, 78.6 GB per GPU, about 1.1 h); falls back to Dp = 8 if the device memory at step 100 exceeds the planner's bound | record |
| Blast wave, larger R_m / d_i | `blast_wave/rm12` (576³, 68.6 GB per GPU, about 14 min) | record |
| Blast wave, Gaussian cloud | `blast_wave/cloud` (r_s = 0, width 0.33 R_m, no body; a width of 0.1 R_m is out of reach: peak density 752 n0) | reported |
| Beam–plasma with B0 | `two_stream/ei_b0`, `filamentation/pair_b0`, `filamentation/pair_b0_stable` | gate |
| Magnetised shock | `relativistic_shock/pair_2d_sigma` (σ = 0.3, field out of the plane) | gate |
| Long runs, reconnection | `reconnection/pair_sigma10_xl` (3200×1600×4, 32,000 steps) | record |
| Long runs, shock tail | `relativistic_shock/pair_2d_t4000`, `pair_2d_tail` (32,512×160×4, 3.3e8 particles, 41,300 steps) | record |

### Stage 1 — bow shock: `south`, `oblique`, `dp10`

`BowShockParams.planet` (bool, default true) leaves out the sphere and its dipole: the null test.

O1 (null test of the oblique wind, before `oblique`): the wind of the gate with `imf = (1, 0, 1)`
and no planet, four crossings, the limits of M2 of `docs/plans/2026-11-10-bow-shock.md` (Gauss
≤ 1e-10 at every check; density of each species and region within 5 % over the last crossing;
mean velocity per component within 0.05 |v|; T_e of the interior at the end ≤ 1.5 T0; electric and
magnetic energy of the difference fields, last crossing over the one before, ≤ 1.5). The slope of
the magnetic energy is reported (in M2 the oblique case had not levelled after three crossings).
Stop rule: if O1 fails, `oblique` is reported and is not a gate.

Definitions of O1 that M2 left to its analysis script (which is not in the repository), fixed
here before the run and implemented as `null_test` in `examples/bow_shock/analyze.py`:

- Samples: the 3-D output (every 480 steps, blocks of 8 cells). Regions: the interior (two blocks
  or more from every face) and the slab of two blocks at each of the six faces. References: the
  sample of step 0. "Last crossing": the samples of the last 1920 steps.
- Density: Σn of the region over its value at step 0. Velocity: Σflux / Σn per component minus the
  wind velocity, over |v|. Both limits apply to the samples of the last crossing (those with a step
  ≥ N − 1920, as in the upstream check of the gate); the largest velocity deviation over all
  samples is reported.
- T: (2/3) (Σekin / Σn − m (γ_u − 1)) with γ_u of the mean velocity |Σflux / Σn|, over its value
  at step 0 (for a drifting Maxwellian this is proportional to T at a fixed drift).
- Energies: the mean of the scalar-history records `electric` and `magnetic` with a step in
  (N − 1920, N], over the mean of those in (N − 3840, N − 1920], N the last step. Reported slope:
  the linear fit of `magnetic` over the last crossing, per 480 steps, over its mean.

Definition of the reversal of `south`: the profile of B_z on the subsolar line is taken like the
density profile of the gate (slice z = z_c, band of ±0.1 R_mp, running mean over 0.1 R_mp, mean of
the outputs of the last half crossing). Its changes of sign are counted between the shock and the
surface of the body; the criterion needs exactly one, and x_rev is its linear interpolation.

Reported ±z asymmetry: on the line x = x_c of the plane y = y_c (same band, running mean and
outputs), the distance from the centre of the outermost point with n_i ≥ 1.5 n_0 towards +z and
towards −z; asymmetry = (plus − minus) / (plus + minus).

The compression of an oblique fast shock goes into `examples/bow_shock/theory.py`; at θ = 90° it
must reproduce `pypic.setup.mhd.compression`. `analyze.py` gains: the point x_rev on the subsolar
line where B_z changes sign; the signed B_z / |B_imf| in the movie; the choice of the upstream
region by z.

| Criterion | Limit | |
|---|---|---|
| Gauss; sheath compression; shock upstream of the magnetopause | 1e-10; ≥ 2; holds | gate, both |
| Magnetopause by density | within 25 % of R_mp | gate, both |
| `south`: B_z reverses once between the shock and the body | \|(x_c − x_rev) / R_mp − 1\| ≤ 0.25 | gate |
| Upstream n and u_x (`oblique`: the region z − z_c ≥ 1.5 R_mp) | 5 % | gate |
| Side of the foreshock, compression against the oblique jump, stand-off of Spreiter, ±z asymmetry, erosion x_rev − x_mp | — | reported |

`dp10` is run once and recorded (ANA-D3) with the criteria of `bow_shock:gate` evaluated and
reported.

### Stage 2 — blast wave: `cloud`, `rm12`

Both are run once; the criteria of `blast_wave:gate` are evaluated and reported. The result is
Δ98 against Eq. 4 of Le et al. 2021: whether the +21 % of the gate depends on R_m / d_i and on the
shape of the debris.

### Stage 3 — beam–plasma with B0

`ei_b0`: B0 along the drift, Ω_ce = ω_pe. The electrostatic mode along B0 is unchanged, so the
growth rate of the k_max mode stays within 10 % of the unmagnetised theory.

Filamentation with a flow-aligned B0 (symmetric counter-streaming beams, Ω = q B0 / (γ0 m)):
Γ² = ½ [−(s + Ω²) + √((s − Ω²)² + 4 k² A)], s = k² + b, with A and b the coefficients of the
unmagnetised relation of `pypic/theory/streaming.py`. For Ω² ≥ A every k is stable. Implemented as
`filamentation_growth(..., omega_c=0.0)`. `filamentation/analyze.py` takes the linear interval
from the energy of the modes (the total magnetic energy contains B0).

| Criterion | Limit | |
|---|---|---|
| `ei_b0`: growth rate of the k_max mode | 10 % | gate |
| `pair_b0` (Ω = 0.5 Ω_c): median error of the growth rates of m = 4…10 | ≤ 10 % | gate |
| `pair_b0_stable` (Ω = 1.5 Ω_c): Σ\|B_m\|² (m = 2…10), end over 10 % of the run | ≤ e² | gate |
| Oblique two-stream modes, saturation | — | reported |

Ω_c is the field of marginal stability, Ω_c² = A.

Definitions fixed before the runs (2026-11-11): `pair_b0_stable` has no growth time of its own and
runs as long as `pair` (25 / Γ_max of the case without the field, 1128 steps); its criterion takes
Σ|B_m|² at the last output over that at the first output from 10 % of the run. The fits of
`pair_b0` use the modes from `gate_m_min: 4` (m = 1 is stable there and m = 2 marginal, Γ = 0.055).

Check of the relation (2026-11-11): no paper text is at hand, so it was not compared with a
printed formula. It was derived from the cold relativistic fluid equations and is compared in
`tests/cpu/test_theory.py` with the largest eigenvalue of the linearised fluid + Maxwell system
(pair and electron–ion beams, Ω / Ω_c = 0, 0.5, 0.9, 1.5, three k each): agreement to 1e-7. The
field of marginal stability, ω_ce = β √γ ω_p, is the known one (Godfrey et al. 1975; Bret et al.).

### Stage 4 — PIC-D18 and the magnetised shock

PIC-D18: an initial field `E0` tangential to a *reflecting* face is allowed (explicit solver; an
open face keeps refusing it). `_set_initial_fields` already applies the walls after loading the
field, which is the piston set-up: the wall zeroes the tangential field on its nodes and the
plasma that carries E = −v × B runs into it. `external_fields.E` tangential to a reflecting face
is refused (today nothing refuses it, and the wall acts only on the difference field). Files:
`pypic/simulation.py`, `pypic/initial_fields.py`, `pypic/external.py`. GPU test: a cold plasma in
E × B drift between reflecting walls feels no force in the interior, Gauss ≤ 1e-10 (the code
before the change raises).

Example: `ShockParams.sigma`, `field`. The magnetised jump conditions go into
`examples/relativistic_shock/theory.py` (σ → 0 reproduces `jump_conditions`; γ0 → ∞ gives the form
of Kennel & Coroniti 1984).

| Criterion (`pair_2d_sigma`) | Limit | |
|---|---|---|
| β_sh and n2 / n1 against the magnetised jump | 10 % | gate |
| Downstream ⟨B_z⟩ / B1 against r | 10 % | gate |
| Gauss | 1e-10 | gate |
| Precursor, in-plane field, σ = 0.1 and 1 | — | reported |

### Stage 5 — long runs

Each run starts with a probe of 200 steps (ms per step, fill, device memory). If the estimate
exceeds three hours the run is not started and is reported. Recorded with
`scripts/record_example.py`, FP64. Reported: power-law index, reconnection rate, γ_max(t), the
cut-off against the box size.

## Invariants

- PIC-D18 changes which set-ups are accepted, not the step: the Gauss law, div B and the 1-GPU vs
  N-GPU history keep their limits, and the existing gates must reproduce their recorded values.
- No Poisson correction: the piston set-up starts from a neutral plasma with a divergence-free E0.
- `pypic/` stays experiment-agnostic: the only theory added to `pypic/theory/` is the shared
  filamentation relation; the oblique and magnetised jump conditions live in their examples.

## Verification

- CPU: `ruff check .`, `python scripts/check_kernels.py`, `pytest tests/cpu`.
- GPU, on the clean head commit: `pytest tests/gpu`, then `python scripts/gpu_gates.py --gpus 8`
  and `--precision mixed`, with the new example gates in the list; record and index entry
  committed, `PROGRESS.md` updated.
- The wall test of PIC-D18 is run once on the code before the change and must fail there.
- Long runs: `scripts/record_example.py`; results in the README of each example and in the wiki.
- End state: R6 moved to Done in `docs/roadmap.md`; PIC-D18 Accepted.

## Results

(Filled in stage by stage.)

### Stage 1 — O1, null test of the oblique wind (2026-11-11): passed

`bow_shock/oblique` with `planet: false`, 8 GPUs, FP64, 7680 steps (four crossings), 3.93e8
particles, 41 ms per step. Evaluated by `null_test` of `examples/bow_shock/analyze.py`
(17 samples, 5 of them in the last crossing).

| Criterion | Limit | Measured | |
|---|---|---|---|
| Gauss, every check | 1e-10 | 2.7e-15 | pass |
| Density, last crossing, worst region | 5 % | e 0.42 %, i 0.41 % (slab y+, last sample) | pass |
| Mean velocity per component, worst region | 0.05 \|v\| | e 8.2e-4, i 7.6e-4 (slab z+); the same over all samples | pass |
| T_e of the interior at the end | 1.5 T0 | 1.003 (ions 1.001) | pass |
| Electric energy, last crossing over the one before | 1.5 | 0.99987 | pass |
| Magnetic energy, last crossing over the one before | 1.5 | 0.99962 | pass |

Reported: slope of the magnetic energy over the last crossing −3.8e-4 of its mean per 480 steps
(levelled; it had not in M2 after three crossings); largest T_e of any region at the end 1.021
(slab y+); the energy monitor of the explicit scheme drifts by +8.3e-3 over the run (monitored
only); particles 393,216,311 → 393,126,101. `oblique` is therefore a gate.

### Stage 1 — `bow_shock/south` (2026-11-11): gate passed

8 GPUs, FP64, 5760 steps (three crossings), 3.92e8 particles at the start, 46 ms per step
(wall 264 s). Evaluated by `analyze` of `examples/bow_shock/analyze.py`.

| Criterion | Limit | Measured | |
|---|---|---|---|
| Gauss, every check | 1e-10 | 6.1e-14 | pass |
| Upstream density and u_x, last crossing | 5 % | n: e 0.13 %, i 0.12 %; u_x: e 0.03 %, i 0.04 % | pass |
| Sheath compression on the subsolar line | ≥ 2 | 2.48 (at 28.25 d_e; MHD 2.97) | pass |
| Magnetopause by density | within 25 % of R_mp | 0.910 R_mp | pass |
| Shock upstream of the magnetopause | holds | shock at 1.178 R_mp | pass |
| B_z reverses once between shock and body | \|(x_c − x_rev) / R_mp − 1\| ≤ 0.25 | one reversal, 0.005 (x_rev at 1.005 R_mp from the centre) | pass |

Reported: the reversal lies 0.094 R_mp (1.9 d_e) upstream of the magnetopause by density
(erosion x_rev − x_mp = −0.094 R_mp: the density falls back to the wind value inside the field
reversal, not outside it). Against the northward gate of the same wind (README: magnetopause
0.92 to 0.93 R_mp, shock 1.38 R_mp, largest density 3.8 n0, sheath over magnetopause distance
0.49 to 0.50): the magnetopause is at the same place, the shock is 0.20 R_mp closer (Spreiter
1.371 R_mp), the largest density is 2.48 n0 and the sheath is 0.295 of the magnetopause distance
(Spreiter 0.371). Magnetopause by magnetic pressure 0.883 R_mp. Flank (line x = x_c, plane
y = y_c): +z 2.236 R_mp, −z 2.300 R_mp, asymmetry −0.014. Body charge +432; absorbed particles
per step 1138, 403, 311 in the three crossings; particles 3.919e8 → 4.300e8; largest fill ratio
0.57; energy monitor +2.7e-2 (monitored only).
