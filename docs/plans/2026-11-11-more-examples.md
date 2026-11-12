# More examples: the seven open items of R6

Status: Done (PR #41) — the seven items have a gated or a recorded configuration; gate records on c54e5d1 in FP64 and with `precision: mixed`, ALL PASSED (Results, Gate records on c54e5d1); PIC-D18 Accepted
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

Check of the oblique jump (2026-11-11): no text of Priest (1982) is at hand, so the relation was not
compared with the printed formula. It was derived from the conservation of mass, of both momentum
components, of the tangential electric field and of the energy across a plane shock, and
`tests/cpu/test_theory.py` checks it with the equation the derivation of the downstream state
leaves over (the energy flux, to 1e-9, four winds and five angles) and with its limits: θ = 90°
gives `pypic.setup.mhd.compression` exactly, θ = 0 the gas-dynamic shock of the sonic Mach number
or the switch-on shock X = M_A².

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

Amendment (2026-11-11, before the implementation): the refusal of `external_fields.E` tangential
to a reflecting face is dropped. New fact: PIC-D12 (Accepted) allows an external field tangential
to any face, and `tests/gpu/mgpu_migrate.py` (variant `dipole`) runs with an external E tangential
to its reflecting z faces; the refusal would have removed accepted behaviour. PIC-D18 states the
difference instead: `E0` is a Yee field and the wall zeroes the total tangential field; an
external E is not shorted by a wall. `pypic/initial_fields.py` and `pypic/external.py` are not
changed (plane waves and 1-D profiles keep refusing a tangential E). The test region is the
interior farther than c t + 2 from both walls: each wall sends a front into the plasma at c.

Example: `ShockParams.sigma`, `field`. The magnetised jump conditions go into
`examples/relativistic_shock/theory.py` (σ → 0 reproduces `jump_conditions`; γ0 → ∞ gives the form
of Kennel & Coroniti 1984).

| Criterion (`pair_2d_sigma`) | Limit | |
|---|---|---|
| β_sh and n2 / n1 against the magnetised jump | 10 % | gate |
| Downstream ⟨B_z⟩ / B1 against r | 10 % | gate |
| Gauss | 1e-10 | gate |
| Precursor, in-plane field, σ = 0.1 and 1 | — | reported |

Check of the magnetised jump (2026-11-11): no text of Kennel & Coroniti (1984) is at hand, so the
relation was not compared with the printed formula. The two conservation equations of `theory.py`
were derived in the downstream frame (particle flux, momentum flux with the Maxwell stress, energy
flux with the Poynting flux, flux freezing). `tests/cpu/test_theory.py` checks σ → 0 against
`jump_conditions` and γ0 → ∞ in 3-D against the Kennel–Coroniti expression for the downstream
four-velocity in the shock frame as remembered (agreement to 4e-7 at γ0 = 1e3); a misremembered
coefficient there would fail that test, an error common to both would not be seen.

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
0.49 to 0.50): the magnetopause is 0.01 to 0.02 R_mp closer, the shock is 0.20 R_mp closer (Spreiter
1.371 R_mp), the largest density is 2.48 n0 and the sheath is 0.295 of the magnetopause distance
(Spreiter 0.371). Magnetopause by magnetic pressure 0.883 R_mp. Flank (line x = x_c, plane
y = y_c): +z 2.236 R_mp, −z 2.300 R_mp, asymmetry −0.014. Body charge +432; absorbed particles
per step 1138, 403, 311 in the three crossings; particles 3.919e8 → 4.300e8; largest fill ratio
0.57; energy monitor +2.7e-2 (monitored only).

### Stage 1 — `bow_shock/oblique` (2026-11-11): gate passed

8 GPUs, FP64, 5760 steps (three crossings), 3.92e8 particles at the start, 47 ms per step
(wall 269 s). IMF along (1, 0, 1) / √2 (θ_Bn = 45° on the subsolar line); the upstream region of
the gate is z − z_c ≥ 1.5 R_mp.

| Criterion | Limit | Measured | |
|---|---|---|---|
| Gauss, every check | 1e-10 | 3.2e-14 | pass |
| Upstream density and u_x (z − z_c ≥ 1.5 R_mp), last crossing | 5 % | n: e 0.11 %, i 0.12 %; u_x: e 0.09 %, i 0.02 % | pass |
| Sheath compression on the subsolar line | ≥ 2 | 3.26 (at 28.25 d_e) | pass |
| Magnetopause by density | within 25 % of R_mp | 0.960 R_mp | pass |
| Shock upstream of the magnetopause | holds | shock at 1.251 R_mp | pass |

Reported: the compression is 2.4 % above the oblique MHD jump (3.19); the shock stands at
1.251 R_mp against 1.345 R_mp of Spreiter (sheath over magnetopause distance 0.303 against 0.345).
The mirrored region (z − z_c ≤ −1.5 R_mp, the side the field lines of the shock reach upstream)
deviates from the wind by 0.09 % in density and 0.08 % in u_x: no foreshock is seen in the moments
at this size and length (22.5 / Ω_ci). The flank is not symmetric: the shock on the line x = x_c of
the plane y = y_c is at +z 2.382 R_mp and −z 1.574 R_mp (asymmetry +0.204; `south` −0.014).
Magnetopause by magnetic pressure 0.958 R_mp. Body charge +522; absorbed particles per step 1529,
491, 357 in the three crossings; particles 3.919e8 → 4.200e8; largest fill ratio 0.59; energy
monitor +1.4e-2 (monitored only).

### Stage 3 — beam–plasma with B0 (2026-11-11): three gates passed

Targeted runs on 8 GPUs, FP64 (the working tree of 3ab48e0; `pypic/` differs from it only by the
PIC-D18 check of `E0`, which these configurations do not reach). Each run takes 5 to 8 s.

| Gate | Criterion | Limit | Measured | |
|---|---|---|---|---|
| `two_stream/ei_b0` | growth rate of the k_max mode | 10 % | +2.2 % (fit 0.3491, theory 0.3415) | pass |
| | Gauss, div B | 1e-10, 1e-12 | 1.9e-12, 4.9e-13 | pass |
| `filamentation/pair_b0` (Ω = 0.5 Ω_c) | median growth-rate error, m = 4…10 | ≤ 10 % | 5.7 % (largest 18 %, m = 4) | pass |
| | Gauss | 1e-10 | 6.3e-13 | pass |
| `filamentation/pair_b0_stable` (Ω = 1.5 Ω_c) | Σ\|B_m\|² (m = 2…10), end over 10 % of the run | ≤ e² = 7.39 | 6.02 | pass |
| | Gauss | 1e-10 | 3.5e-14 | pass |
| `filamentation/pair` (no field, existing gate) | median growth-rate error, m = 2…10 | ≤ 10 % | 2.2 % (unchanged analysis path) | pass |

Reported:
- **Oblique two-stream modes** (`ei_b0` against `ei` of the same working tree, from the 3-D field
  outputs): the growth rate of the k_max mode is the same (`ei` +3.0 %), the saturation is not.
  Without the field 82 % of Σ|E_k|² is in modes with k_⊥ ≠ 0 at t = 29.3 / ω_pe and the electric
  energy peaks at 5.3 % of the initial electron kinetic energy (t = 21.1); with Ω_ce = ω_pe along the
  drift 2.0 % is in k_⊥ ≠ 0 at the same time (17 % at the end, t = 87.8), E_x holds 99 % of |E|², the
  peak is 18.8 % (t = 37.3) and the amplitude of the (k_max, 0, 0) mode peaks 8.8 times higher
  (0.0326 against 0.0037). The electrons lose 4.7 % of their kinetic energy by the end (0.65 % without
  the field).
- **Filamentation, saturation:** the magnetic energy of the perturbation peaks at 19.5 % of the
  initial kinetic energy with and without the field (`pair_b0` at t = 33.9, `pair` at t = 31.7): at
  Ω = 0.5 Ω_c the field lowers the rates (m = 6: fit 0.410 against 0.497; theory 0.417 against 0.510)
  and not the level. In `pair_b0` the modes m = 1 (stable in the linear theory) and m = 2 (Γ = 0.055)
  rise at 0.134 and 0.248 in the fit window: they are driven by the faster modes, which is why the
  gate starts at m = 4.
- **Stable case:** the mode energy grows 6.0-fold between 10 % of the run and its end with no
  unstable k (the unstable case would give e^(2 Γ t) with Γ t ≈ 22). The margin to the limit e² is
  small (19 %).

### Stage 2 — `blast_wave/cloud` (2026-11-11): reported, all criteria of the gate hold

The gate configuration without the body and its dipole, the debris a Gaussian ball (shell radius 0)
of width 0.33 R_m = 10.56 d_e, peak density 20.9 n0. 384³ cells, 8 GPUs, FP64, 1760 steps,
1.185e9 particles, 109 ms per step (wall 192 s, 0.78 ns per particle and step). Launched from the
working tree of 23b348e with the PIC-D18 check of `E0` uncommitted (this configuration does not
reach it), so the run is reported and not recorded.

| Criterion of `blast_wave:gate` | Limit | Measured | |
|---|---|---|---|
| Gauss, div B | 1e-10, 1e-12 | 7.4e-15, 8.6e-15 | holds |
| Front speed between 0.3 and 0.7 R_B | above the fast speed (0.0774 c) | 0.1371 c = 1.77 fast speeds = 0.73 V_d | holds |
| Largest ring mean of B_z / B0 there | above 1.5 | 2.08 | holds |
| Front radii, axes against diagonals | 5 % | 0.45 % | holds |
| Δ98 of the debris ions against Eq. 4 of Le et al. 2021 (D = 3), reached before the end | ± 25 % | 39.0 d_e = 1.22 R_m at step 572 of 1760; Eq. 4: 32.6 d_e; +19.5 % | holds |
| Ambient n_e, n_i, \|B\| ahead of the front | 5 % | 0.045 % | holds |

Reported: the shell of the gate gives Δ98 = 39.5 d_e (+21 %) with the same R_m, M_A and seed, so
the excess over Eq. 4 changes by 1.5 points between a shell at 0.33 R_m around a body with a dipole
and a ball at the centre without them. The front is slower than that of the shell (1.77 against
2.10 fast speeds) and the ring mean is lower (2.08 against 2.69; perpendicular MHD shock at the
speed of the front: 1.95): the ball starts 0.33 R_m further in and the 98 % extent of the debris is
17.5 d_e at step 0. Front at the end at 77.0 d_e = 1.16 R_B (shell: 1.16 R_B); cavity 34.25 d_e =
0.51 R_B (shell 0.50 R_B); distance between front and debris 10.7 d_e at step 660; debris rms extent
36.8 d_e along B0 and 14.9 d_e across; the debris ions keep 21.5 % of their energy and the ambient
ions gain 38.3 % of it; energy monitor +1.9e-2 (budget change with the boundary terms 0.107;
monitored only); particles 1.185e9 → 1.115e9; largest fill ratio 0.670. 200-step probe before the
run: 123 ms per step, peak device memory 29.3 GiB = 31.5 GB on the largest GPU against the
planner's 32.4 GB.

### Stage 2 — `blast_wave/rm12` (2026-11-11): recorded, all criteria of the gate hold

The gate configuration with R_m = 12 d_i = 48 d_e (R_B = 99.8 d_e, debris gyroradius 0.25 R_m
instead of 0.375 R_m; shell and body scale with R_m). 576³ cells, 8 GPUs, FP64, 2640 steps
(10.3 / Ω_ci), 3.183e9 particles, 286 ms per step (wall 757 s, 0.73 ns per particle and step).
Launched from the clean tree of c58f533. Record:
`bench/results/2026-11-11-c58f533c5-example-blast_wave-rm12-run-031fb299.json`.

| Criterion of `blast_wave:gate` | Limit | Measured | |
|---|---|---|---|
| Gauss, div B | 1e-10, 1e-12 | 7.3e-15, 1.3e-14 | holds |
| Front speed between 0.3 and 0.7 R_B | above the fast speed (0.0774 c) | 0.1616 c = 2.09 fast speeds = 0.86 V_d | holds |
| Largest ring mean of B_z / B0 there | above 1.5 | 2.87 | holds |
| Front radii, axes against diagonals | 5 % | 0.16 % | holds |
| Δ98 of the debris ions against Eq. 4 of Le et al. 2021 (D = 3), reached before the end | ± 25 % | 56.4 d_e = 1.17 R_m at step 2376 of 2640; Eq. 4: 48.9 d_e; +15.2 % | holds |
| Ambient n_e, n_i, \|B\| ahead of the front | 5 % | 0.082 % | holds |

Reported:
- **The excess over Eq. 4 is 15 % at R_m = 12 d_i, against 21 % at 8 d_i** (same M_A, m_i/m_e,
  shell shape in units of R_m and seed). One run each with one seed; repeated runs of the gate
  with that seed agree in Δ98 to 10 digits or more, another seed was not run.
- **The time history of Δ98 is not that of the gate.** Δ98 (41 samples, every 66 steps) has a
  first maximum of 44.4 d_e = 0.93 R_m at step 528 (9 % *below* Eq. 4), falls to 44.1 d_e at
  step 594 and then rises without a further dip to 56.4 d_e at step 2376; the last four samples
  fall by 0.13 d_e (0.2 %), and the nine samples from step 2112 on lie within 0.25 d_e. The
  criterion "reached before the end" holds by these four samples: the value is a plateau at the
  end of the run and no fall like that of the gate (39.5 → 34.3 d_e) was seen. In the gate and in
  `cloud` the first maximum is the largest one (`cloud`: 39.0 d_e at step 572, second maximum
  38.2 d_e at step 1056).
- Front: 2.09 fast speeds (gate 2.10); ring mean 2.87 (gate 2.69; perpendicular MHD shock at the
  speed of the front 2.25); at the end at 117.5 d_e = 1.18 R_B (gate 1.16 R_B), 17.9 d_e ahead of
  Δ98 when it leaves the window (step 1056). Cavity 50.25 d_e = 0.50 R_B = 1.05 R_m (gate
  0.50 R_B; the static dipole alone 12.1 d_e).
- Debris rms extent at the end 52.7 d_e along B0 and 23.4 d_e across; the debris ions keep 13.5 %
  of their energy (gate 20 %) and the ambient ions gain 46.5 % of it (gate 41 %); energy monitor
  +2.7e-2 (budget change with the boundary terms 0.152 of the debris energy; monitored only).
- Body: 306,982 particles absorbed, charge +33.1 at the end (largest 65.6); potential at two
  dipole radii 0.08 T_e / e across B0 and 0.41 T_e / e along it (gate 4.9 and 1.3).
- Particles 3.183e9 → 3.161e9; largest fill ratio 0.671.
- **Device memory: the run needs more than the planner's figure.** 200-step probe: 293 ms per
  step, 64,297 MiB = 67.4 GB on the fullest GPU. During the run the fullest GPU reached
  74,489 MiB = 78.1 GB (the others 46 to 71 GiB), against 68.6 GB of the planner for full particle
  arrays: 14 % above it, with the arrays 67 % full. A 200-step probe therefore does not bound the
  memory of a run; `dp10` below is judged with this in mind.

### Stage 1 — `bow_shock/dp10` (2026-11-11): recorded, all criteria of the gate hold; no fallback to Dp = 8

The gate wind with R_mp = 10 d_i = 40 d_e = 80 cells (gate: 5 d_i). 480×640×640 cells, 8 GPUs,
FP64, 11,520 steps (3 crossings = 45 / Ω_ci), 3.136e9 particles at the start, 313 ms per step at
the end (293 at the start; wall 3607 s, 0.72 ns per particle and step). Launched from the clean
tree of 4d55109. Record:
`bench/results/2026-11-11-4d55109b5-example-bow_shock-dp10-run-4b073ac2.json`.

| Criterion of `bow_shock:gate` | Limit | Measured | |
|---|---|---|---|
| Gauss (div B, limit 1e-12) | 1e-10 | 4.7e-14 (1.4e-13) | holds |
| Upstream n and u_x of both species, last crossing | 5 % | 0.054 % | holds |
| Ion density compression in the sheath on the subsolar line | ≥ 2 | 3.77 | holds |
| Magnetopause by density | within 25 % of R_mp | 0.940 R_mp | holds |
| Shock upstream of the magnetopause | holds | shock at 1.439 R_mp | holds |

Reported:
- **Fallback rule of the table above**: the 200-step probe (295 ms per step, fill 0.50) reached
  69,457 MiB = 72.8 GB on the fullest GPU, below the planner's 78.6 GB, so the run stayed at
  Dp = 10. **During the run the fullest GPU reached 75,017 MiB = 78.66 GB: 0.08 % above the
  planner's figure for full particle arrays, with the arrays 57 % full** (the value of
  `nvidia-smi`, which includes the blocks that the memory pool holds; the live part was not
  measured). As for `rm12`, the planner's figure is not a bound on the device memory of a run;
  24 GB of the 102.6 GB (97,887 MiB) of a device were left.
- Shock at 1.439 R_mp against 1.371 of Spreiter et al. (1966) with the MHD compression: +4.9 %
  (gate: 1.38 against 1.37). Magnetopause by magnetic pressure 0.973 R_mp. Sheath thickness over
  magnetopause distance 0.530 (Spreiter 0.371; gate 0.49 to 0.50): 19.9 d_e = 5.0 d_i, against
  9 d_e ≈ 2 d_i in the gate.
- Largest density on the subsolar line 3.77 n0 (MHD 2.97; gate 3.8): the maximum does not come
  down to the MHD value with a sheath of 5 d_i. One run; whether a plateau follows the maximum
  was not analysed.
- Flank (line x = x_c, plane y = y_c): shock at 2.354 R_mp on +z and 2.341 R_mp on −z, 0.26 %
  apart (northward field, symmetric as expected).
- **Not stationary at the end**: the particle number grows by 10.5 % over the run (3.136e9 →
  3.465e9), by 0.43 % over the last 2400 steps and by 0.13 % over the last 960. The body absorbs
  5752, 2294 and 926 particles per step in the three crossings (3.44e7 in all); its charge is
  +1668 at the end (largest +2739 at t = 1160, smallest −478).
- Largest fill ratio 0.567, no reallocation. Energy (monitored only): total +2.3 %, budget with
  the boundary terms +3.3 %.

### Stage 4 — PIC-D18 and `relativistic_shock/pair_2d_sigma` (2026-11-11): gate criteria passed

The wall test of PIC-D18 (`test_exb_drifting_plasma_between_reflecting_walls_keeps_its_drift`)
passes on the branch and was run once on the code before the change, where the construction
raises `ValueError: E0 component y is tangential to the conducting faces of axis x`.

`pair_2d_sigma` (γ0 = 15, σ = 0.3, B1 = 3.0 along z, E_y = −β0 B1; 1792×80×4 cells, 2066 steps,
4.6e6 particles), targeted run on 8 GPUs, FP64, on ec0ce42 (7 s). It is the first run of the
configuration; the limits were fixed before it.

| Criterion | Limit | Measured | |
|---|---|---|---|
| β_sh against the magnetised jump (2-D gas: 0.6448) | 10 % | 0.5984, −7.2 % | pass |
| n2 / n1 against the jump (2.547) | 10 % | 2.666, +4.7 % | pass |
| Downstream ⟨B_z⟩ / B1 against r (2.547) | 10 % | 2.651, +4.1 % | pass |
| Gauss | 1e-10 | 1.2e-12 | pass |

The existing gate `pair_2d` (no field, same analysis code) in the same session: β_sh −5.3 %,
n2 / n1 +1.9 %, Gauss 2.4e-12 (README: −5.1 %, +1.7 %; the runs are independent realizations).
The gate record is taken with the full gates at the final head of the pull request.

Reported (one run each, 8 GPUs, FP64, `--set` on `pair_2d_sigma`; the gate criteria are evaluated
and are not gates here):

| Run | Jump used | β_sh | n2 / n1 | ⟨B⟩ / B1 downstream | Precursor δB / B1 (rms, max) | Density peak at the front |
|---|---|---|---|---|---|---|
| σ = 0.3, B along z (gate) | 2-D gas: 0.645, 2.547 | 0.598 (−7.2 %) | 2.666 (+4.7 %) | 2.651 (+4.1 %) | 0.71, 1.9 | 11.4 n1 |
| σ = 0.3, B along y (`field=in`) | 3-D gas: 0.549, 2.817 | 0.571 (+4.0 %) | 2.800 (−0.6 %) | 2.746 (−2.5 %) | 0.59, 2.2 | 8.1 n1 |
| σ = 0.1, B along z | 2-D gas: 0.561, 2.780 | 0.522 (−6.8 %) | 2.995 (+7.7 %) | 2.972 (+6.9 %) | 1.08, 3.5 | 5.9 n1 |
| σ = 1, B along z | 2-D gas: 0.783, 2.274 | 0.809 (+3.3 %) | 2.262 (−0.6 %) | 2.281 (+0.3 %) | 0.34, 0.73 | 12.5 n1 |

- **Flux freezing holds**: the mean downstream field over B1 is within 2 % of the measured
  density compression in the four runs (−0.6 %, −2.0 %, −0.8 %, +0.8 %).
- **In-plane field**: the particles gyrate out of the simulation plane, the gas has three degrees
  of freedom and the run follows the 3-D jump (−0.6 % in n2 / n1), not the 2-D one (2.547: +10 %).
- **Precursor**: an electromagnetic wave train runs ahead of the shock in all four runs
  (line-out of the field component along B1, measured at t ≈ 0.7 of the validity time between
  5 skin depths ahead of the shock and the farthest point that deviates by 5 % from B1). Its
  leading edge is at 278 to 285 c/ω_pe at t = 294 to 302, i.e. it has travelled at c since the
  start of the run: it contains the front that the wall launches when it shorts E_y at t = 0
  (PIC-D18) and the emission of the forming shock; the two are not separated. The amplitude falls
  with σ between 0.1 and 1. The density in the wave train deviates by 1 % rms. No comparison
  with a published precursor amplitude was made.
- The deviations of β_sh and n2 / n1 have opposite signs in every run, as mass conservation asks
  (β_sh = β0 / (r − 1) in the wall frame); the pair measured at σ = 0.3 satisfies it within 0.1 %,
  the others within 2 to 4 % (the speed is a fit over a time window, the density a median at its
  end).
- Downstream spectrum at t = 397: no power-law tail in any of the runs (the fit finds none);
  temperature of the Maxwell–Jüttner fit 3.9 (σ = 0.3), 4.9 (in-plane), 3.7 (σ = 0.1), 3.7
  (σ = 1) against 6.8, 4.3, 7.1, 6.4 of the jump; the gate without a field gives 4.2 against 7.4
  with the same fit, so the fit of the quasi-2-D runs is low by the same factor with and without
  the field and is not a statement about the heating.

### Stage 5 — long runs (2026-11-11): three runs recorded; the shock tail does not develop by t = 8000

All three on 8 GPUs, FP64, launched from the clean tree of 2945145, one run each. Each started
with the 200-step probe; no estimate came near the three-hour limit.

| Run | Cells | Particles | Probe: ms / step, fill, estimate | Run: steps, wall, ms / step at the end | Energy budget | Gauss |
|---|---|---|---|---|---|---|
| `reconnection/pair_sigma10_xl` | 3200×1600×4 | 3.68e8 | 70.4, 0.29, 38 min | 32,000, 1314 s, 41.1 | −0.34 % | 2.8e-12 |
| `relativistic_shock/pair_2d_t4000` | 16256×160×4 | 1.66e8 | 37.5, 0.80, 13 min | 20,656, 628 s, 30.2 | +0.29 % (step 20,000) | 4.8e-15 |
| `relativistic_shock/pair_2d_tail` | 32512×160×4 | 3.33e8 | 71.1, 0.78, 49 min | 41,312, 2466 s, 59.4 | +0.56 % (step 40,000) | 8.9e-15 |

Records: `bench/results/2026-11-11-29451458d-example-reconnection-pair_sigma10_xl-run-0e8c5abd.json`,
`…-relativistic_shock-pair_2d_t4000-run-014c11db.json`, `…-relativistic_shock-pair_2d_tail-run-b012ba98.json`.
The gate criteria of the two examples hold in the three runs (`passed: true`).

Device memory of the fullest GPU (`nvidia-smi`, blocks held by the memory pool included): 23,129,
20,107 and 39,631 MiB in the probes (sampled every 2 s), 16,235, 9,839 and 18,665 MiB in the runs
(sampled every 10 s). The probes are higher because the largest use is a transient during the
load, which the coarser sampling of the runs missed; the probe values are the ones to plan with.

Reported — reconnection, the cut-off against the box size (`large` is 400 × 200 c/ω_pe to
t = 2000 with two recorded runs, `xl` is 800 × 400 to t = 4000):

| | `large`, t = 2000 | `xl`, t = 1000 | `xl`, t = 2000 | `xl`, t = 3000 | `xl`, t = 4000 |
|---|---|---|---|---|---|
| index p over 5 ≤ γ−1 ≤ 50 | 2.101 ± 0.027, 2.104 ± 0.026 | 2.050 ± 0.029 | 2.048 ± 0.021 | 2.036 ± 0.018 | 2.044 ± 0.017 |
| γ−1 that 10⁻³ of the particles exceed | 90.4, 86.9 | 77.1 | 90.4 | 97.9 | 97.9 |
| the same for 10⁻⁴ | 140.1, 157.8 | 129.4 | 157.8 | 184.9 | 184.9 |
| the same for 10⁻⁵ | 225.5, 285.9 | 200.2 | 264.2 | 274.8 | 285.9 |

- The index is converged in time in the `xl` box from t = 1000 on and is 2.04 at the end, −3.7 %
  from 1.9 + 0.7/√σ = 2.12 (Werner et al. 2016); the half-size box gives 2.10 in both runs. The
  difference of 0.06 is 1.8 times the combined fit error of one run of each box.
- The energy that 10⁻⁴ of the particles exceed is 17 to 32 % higher with twice the box and twice
  the time, not twice as high; at equal time the boxes agree within the spread of the two `large`
  runs; between t = 3000 and 4000 the 10⁻³ and 10⁻⁴ values do not move. The spectrum has 160
  logarithmic bins (steps of 4 %), which is why equal values recur. No cut-off form was fitted and
  there is one `xl` run, so a scaling law with the box is not claimed. These energies were read
  with `pypic.analysis.particles.energy_cutoff` from the `spectrum` output and are not in the
  records.
- Rate: 0.066 v_A on average over the active phase, peak 0.168 (`large`: 0.064 and 0.074, peak
  0.19). Magnetic energy converted 77.7 %. E∥ share of the 1 % most energetic tracked electrons
  −0.004.

Reported — shock tail (`long` is t = 2000 with two recorded runs; the box doubles with the time):

| | `long` (two runs) | `t4000` | `tail` |
|---|---|---|---|
| β_sh against the 2-D jump | −5.9 %, −7.6 % | −7.5 % | −7.0 % |
| n2 / n1 | +6.9 %, +8.1 % | +5.1 % | +3.6 % |
| spectrum taken at t | 1995 | 3919 | 7838 |
| tail index, automatic range | 3.95 on γ 27–93 (first run) | 4.04 ± 0.04 on γ 32–126 | none found (NaN) |
| tail index, window 3 × peak to 0.7 × cut-off | 3.91, 3.85 on γ 40–101 | 4.00 ± 0.06 on γ 40–121 | 3.86 ± 0.10 on γ 43–136 |
| non-thermal fraction | 15–16 % | 16.3 % | 18.3 % |
| γ that 10⁻⁴ of the downstream particles exceed, at the end | 144 | 172 | 195 |
| exponent of its growth over the second half | 0.36, 0.41 | 0.21 | 0.10 |
| magnetic energy at the end over the initial kinetic energy | 0.0079, 0.0067 | 0.0040 | 0.0021 |

- **The tail of the literature is not reached.** t = 8000 ω_pe⁻¹ is 2070 relativistic plasma
  times, near the 2250 at which Spitkovsky (2008b) measured p ≈ 2.4. Here the slope above the
  thermal peak stays at 3.9 to 4.0 from t = 2000 to 8000, and the cut-off rises by 13 % from
  t = 4000 to 8000 (172 for t = 3360–4320, 183 for t = 4480–7520, 195 from t = 7680 in `tail`;
  `t4000` also gives 172 at its end), where t^1/2 gives 41 %. The histogram has 120 logarithmic
  bins up to γ = 1500 (steps of 6.3 %), so it does not cap the values.
- The magnetic energy in the box does not grow with the shocked volume: in units of the kinetic
  energy of the `long` box it is 0.0079 and 0.0067 at t = 2000, 0.0081 at t = 4000 and 0.0085 at
  t = 8000.
- **What this changes**: the README of the example said that a converged index needs
  t ≳ 8000 ω_pe⁻¹ and about 10⁹ particles. The time was reached with 3.3e8 particles
  (L = 1.57 c t) and neither the index nor the cut-off approaches the published values, so time
  alone is not what is missing in this set-up. Not varied, and therefore not separated: the
  transverse size (62 c/ω_pe = 16 relativistic skin depths in all runs), 8 particles per cell,
  the 4-pass current filter, the 4 cells along z, and the region of the spectrum (0.3 to 0.8 of
  the shock position). Roadmap R11 (the converged tail) should start from a scan of the
  transverse size at t = 4000, which costs 10 min per run at the present width, before a larger
  particle number is spent on a longer run.
- These are reported values; no criterion was fixed for them, so none failed.

Gate runs of this pull request: the full gates (`pytest tests/gpu`, `scripts/gpu_gates.py --gpus 8`
in FP64 and mixed) are taken once, on its final head (owner instruction of 2026-11-11); the stages
above were verified with targeted runs.

### Gate records on c54e5d1 (2026-11-12): ALL PASSED in FP64 and mixed

`bench/results/2026-11-12-c54e5d104-gpu-gates.json` (29 gates: 94 GPU tests, 6 two-stream, 22
examples) and `…-gpu-gates-mixed.json` (27 gates: 5 two-stream, 22 examples), 8 GPUs, clean
commit. The largest Gauss residual is 7.6e-12 in FP64 and 7.7e-12 in mixed (`two-stream:buneman`).
The gates this plan adds:

| Gate | Criterion | Limit | FP64 | mixed |
|---|---|---|---|---|
| `two-stream:ei_b0` | growth rate of k_max | 10 % | +2.2 % | +2.2 % |
| `filamentation:pair_b0` | median growth error, m = 4…10 | ≤ 10 % | 5.7 % | 5.7 % |
| `filamentation:pair_b0_stable` | Σ\|B_m\|² (m = 2…10), end over 10 % of the run | ≤ e² = 7.39 | 6.02 | 6.02 |
| `relativistic_shock:pair_2d_sigma` | β_sh | 10 % | −7.2 % | −7.2 % |
| | n₂/n₁ | 10 % | +4.2 % | +4.4 % |
| | ⟨B_z⟩ / B₁ against r | 10 % | +4.0 % | +4.0 % |
| `bow_shock:south` | magnetopause by density | within 25 % of R_mp | 0.883 | 0.846 |
| | reversals of B_z; \|(x_c − x_rev)/R_mp − 1\| | 1; ≤ 0.25 | 1; 0.009 | 1; 0.066 |
| | sheath compression | ≥ 2 | 3.34 | 3.35 |
| | upstream n and u_x | 5 % | 0.13 % | 0.14 % |
| `bow_shock:oblique` | magnetopause by density | within 25 % of R_mp | 0.965 | 0.970 |
| | sheath compression | ≥ 2 | 3.11 | 2.99 |
| | upstream n and u_x (z − z_c ≥ 1.5 R_mp) | 5 % | 0.12 % | 0.12 % |

Gauss's law of the five example gates: ≤ 1.3e-12 (limit 1e-10; the bow-shock runs ≤ 7.3e-14); the
shock stands upstream of the magnetopause in the four bow-shock runs.

**The bow-shock values of the gate runs are not those of the targeted runs of Stage 1.** Three
runs of each configuration with the same seed (R_mp units; reported values below the line of
dashes):

| | `south`, Stage 1 | `south`, FP64 record | `south`, mixed record | `oblique`, Stage 1 | `oblique`, FP64 record | `oblique`, mixed record |
|---|---|---|---|---|---|---|
| magnetopause | 0.910 | 0.883 | 0.846 | 0.960 | 0.965 | 0.970 |
| largest density (n₀) | 2.48 | 3.34 | 3.35 | 3.26 | 3.11 | 2.99 |
| reversal of B_z | 1.005 | 0.991 | 0.934 | — | — | — |
| - - - | | | | | | |
| shock | 1.178 | 1.169 | 1.077 | 1.251 | 1.280 | 1.267 |
| sheath / magnetopause distance | 0.295 | 0.323 | 0.273 | 0.303 | 0.326 | 0.306 |
| erosion x_rev − x_mp | −0.094 | −0.108 | −0.088 | — | — | — |
| flank, +z and −z | 2.236, 2.300 | 2.305, 2.278 | 2.276, 2.244 | 2.38, 1.57 | 2.385, 1.546 | 2.383, 1.570 |

- The analysis did not change between the runs: the present `analyze.py` gives the Stage 1
  values again on the Stage 1 outputs (`south`: 2.4799, 0.9103, 1.1785, 1.0048; `oblique`: 3.263,
  0.9601, 1.2513). The differences are those between runs with the same seed (the order of the
  deposit, roadmap R10; the third run also has the float32 gather), amplified by the flow over
  5760 steps, as for the northward gate (`examples/bow_shock/README.md`).
- They are larger for `south` than for the northward gate (magnetopause 0.85 to 0.91 against
  0.92 to 0.93; largest density 2.48 to 3.35): the layer at the field reversal is not at rest.
  Every run passes every limit; the closest is the magnetopause of the mixed `south` run, 15.4 %
  inside R_mp against the limit of 25 %.
- Consequence for the text: the README of the example and the roadmap now give the ranges of the
  three runs for `south` and `oblique` instead of the values of the first run. No limit and no
  definition changed.

PIC-D18 is Accepted with these records (`relativistic_shock:pair_2d_sigma` on 8 GPUs).

### Review fixes (2026-11-12, after the gates on c54e5d1)

Four findings of the review of PR #41, all in `examples/` (no file under `pypic/` changes):

- `filamentation/plan.py`: `stable` is now Ω² ≥ A (every k), not "no analysed mode grows". A field
  below the marginal one whose unstable modes all lie above `n_modes` is refused (`b0 = 0.98` with
  16 modes: growth starts at m = 17).
- `bow_shock/analyze.py`: the null test fails when a species was not measured (a nan in the ions
  was skipped by `max`).
- `relativistic_shock/analyze.py`: the field ahead of the shock is nan when the shock leaves no
  room ahead (`t_end = 50`), instead of an exception.
- `filamentation/analyze.py`: the plot of the stable case draws only the modes the run analysed
  (`n_modes = 4`).

Checks: a CPU test per finding, each failing on the code before the fix. The plans of the four
shipped filamentation configs are equal before and after (1 and 8 GPUs), and the recorded runs
`relativistic_shock/pair_2d_sigma`, `bow_shock/oblique_null`, `filamentation/pair_b0` and
`pair_b0_stable` give the same analysis results with the old and the new code. The two cases of
the review ran end to end on 1 GPU. Seen there and not examined: `pair_b0_stable` with
`n_modes = 4` (a longer run, three modes in the sum) gives a mode-energy ratio of 9.9, above the
limit e².

### Review fix (2026-11-12, second round): the direction of the field

One finding of the review of 00c07ee, in `examples/filamentation/plan.py`: the stability test
of the first round compared the signed Ω with the marginal one, so a field against the drift
(`b0 < 0`) at or above the marginal strength was taken as unstable and the planner raised
(`b0 = -1`: ZeroDivisionError; `b0 = -1.5`: math domain error; both planned before 00c07ee).
`stable` is Ω² ≥ A now, as in `filamentation_growth`. CPU test: `b0` and `-b0` give the same
verdict, growth rates and run length and opposite `B0`, for 0.5, 1 and 1.5 times the marginal
field (fails on the code before the fix, checked once). No file under `pypic/` changes.

### Amendment (2026-11-12): one pull request

Owner instruction (2026-11-12, after this review): the three pull requests become one, the
fixes go into it and it is gated once. The branch `wakefield` (PR #42) holds this branch by the
merges 031851c and 2263870; PR #41 is closed without a merge of its own, its commits and review
comments stay. The gate records on c54e5d1 are those of this piece alone; the gates of the one
pull request are in `docs/plans/2026-11-11-wakefield.md` (last sections).
