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
