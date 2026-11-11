# Wakefield: laser and beam drivers, the moving window, `examples/wakefield`

Status: In progress
Refs: roadmap R8; PIC-D20, PIC-D21, PIC-D22, OUT-D5 (new); PIC-D17 (extended); ANA-D2, ANA-D3

## Context

Owner request (2026-11-11): close the roadmap rows R3, R5, R6 and R8, each independent piece of
work in its own pull request. This plan is the R8 piece (branch `wakefield`, one pull request;
the stages below are its commits). The other two pieces have their own plans
(`2026-11-11-compact-storage.md`, `2026-11-11-more-examples.md`) and do not touch this file.

The engine has no field source, no moving window and no Poisson solver, and it refuses a bare
bunch that is not neutral (PIC-D17). R8 asks for laser and beam drivers and a moving window; it
is closed when PIC-D20 (antenna), PIC-D21 (bunch) and PIC-D22 (window) are Accepted, the gates
`wakefield:laser_linear` and `wakefield:beam_linear` are in a gate record, and one non-linear
production run is recorded.

Defaults chosen without a question, stated to the owner with the plan and not refused: the beam
driver is a bunch with a co-located, very heavy partner of the opposite charge (an initial field
from a Poisson solve is out of scope); `precision: compact` together with the moving window is out
of scope, and whichever of the two pull requests lands second adds the refusal and its test.

## Approach

All criteria below are fixed here, before the first run of each test and configuration. A
criterion that fails is recorded as failed and reported to the owner; a limit or a definition
changes only by an amendment written into this file with its reason (stop rule).

Every formula that is new here gets a CPU test of its limits and is compared with the cited paper
before the first run; the outcome of that comparison is recorded in the section of the item.

### Stage 1 — laser antenna (PIC-D20)

The current of a magnetisation M given on one plane, `J = curl M`: its discrete divergence is zero,
it uses global indices only, and it is built like `external.dipole_B`. It is added in
`_step_explicit` between `apply_bc_J` and `push_B`. No kernel change and no `#define`. New module
`pypic/antenna.py`. Schema `lasers: [{plane_x, polarization, a0, omega, waist, focus_x, center,
envelope, phase}]`. The energy record gains the work term `antenna`. Refused: the implicit
solver; a plane inside the PML, the retreat region of an open face or a body.

Test A1 (vacuum):

| Criterion | Limit |
|---|---|
| Amplitude and waist at the focus and at one Rayleigh length against the paraxial beam | 2 % |
| Forward flux against backward flux | equal within 1 % |
| Work of the antenna against the Poynting flux | 1 % |
| Gauss, div B | 1e-10, 1e-12 |
| 1 GPU against 2×2 GPUs with the plane on a rank boundary | history ≤ 1e-10 |

Definition fixed before the first run of A1: the paraxial beam of the comparison is the one of
the grid. A wave with a small transverse wavenumber q has `k_x = k - q^2 / (2 k_d)` on the Yee
grid, `k_d = sin(k dx) / dx`, so the antenna and the test use `z_R = k_d w0^2 / 2` (k_d is 0.95 k
at the 12 cells per wavelength of the test, and tends to k on a fine grid). The flux and the
amplitudes are averages over one period of the steady beam of a flat envelope.

CPU: the divergence of curl M is round-off; a sub-block equals the whole box bit for bit; the
amplitude factor of the 1-D Yee grid; k_d against the dispersion relation of the grid.

### Stage 2 — bunch driver (PIC-D21, extends PIC-D17)

`Beam.profile` accepts a list (a product of 1-D shapes); `_beam_block` uses the box of the support
on each axis. The driver is the bunch species plus a co-located partner of the opposite charge and
a very large mass. Every position is neutral, so the initial E = 0 is consistent and Gauss stays
at round-off. The field builds up over a run-up of γ · r_max in vacuum.

Test B1: E_r(r) in the mid-plane of the bunch within 3 % of the field of the Gaussian bunch, and
B_θ = β E_r. An initial field from a Poisson solve is out of scope (reconsidered only if B1 fails
at the run-up within reach).

Definitions fixed before the first run of B1 (2026-11-11):

- The field of the comparison is the steady field of the rigid Gaussian bunch in vacuum, in its
  mid-plane: `bunch_field` in `examples/wakefield/theory.py`. The formula was derived here (the
  potential of a Gaussian at rest in the frame of the bunch); no paper was consulted. It is checked
  on the CPU against its three limits: a point charge, the Gaussian sphere at rest, the line charge.
- Run: vacuum, all faces open, 448×192×192 cells, σ_x = σ_r = σ = 8 cells, dt = dx / 2, γ = 10,
  n_b = 0.05, 8 particles per cell, the bunch and its partner with the mass 1e6 (rigid), the
  partner at rest, the axis of the bunch on a line of nodes, its centre 8 σ from the x− face at
  t = 0. The comparison is made after 640 steps (the bunch has travelled 40 σ).
- E_r at the radii (m + 1/2) dx, m = 4 … 15 (0.56 σ to 1.94 σ): the mean of the four samples of
  E_y and E_z on the two lines through the axis along y and z, interpolated along x to the centre
  of the bunch with the parabola through the three nearest nodes. B_θ in the same way from B_z and
  B_y. Criteria: `|E_r / E_bunch − 1| ≤ 3 %` at every radius, `|B_θ / (β E_r) − 1| ≤ 1 %` at every
  radius, Gauss ≤ 1e-10, div B ≤ 1e-12.
- Correction of the run-up stated above. The mode solution of the start on the partner
  (`bunch_startup`, same file) shows that the part of the start that travels with the bunch loses
  its phase at the rate `k_perp² / (2 k_x) + k_x / (2 γ²)`: the field is formed after a time of the
  order 2 γ² σ_x (200 σ at γ = 10), not γ · r_max. Predicted at the 40 σ of the test:
  E_r is above the steady field by 0.8 % (0.56 σ) to 1.6 % (1.94 σ), and `B_θ / (β E_r) − 1` is
  0.07 % to 0.15 %; at 20 σ the excess is 0.9 % to 1.3 %. The criterion stays the one above (the
  steady field, 3 %); the ratio to the mode solution is reported. Other estimates: the charge shape
  and the grid below 0.3 %, the dispersion of the grid about 0.1 %, the open walls at 12 σ below
  0.5 % at 2 σ.
- Stage 3 takes this into account: the criteria of `wakefield:beam_linear` are evaluated with the
  field the bunch has at the time of the measurement, and the run-up is chosen from the mode
  solution.

### Stage 3 — `examples/wakefield`, static box

N1 (cold plasma, no driver): |E_x| on the axis at most 1 % of the predicted wake, density within
1 %. N2: the same with anisotropic cells. Gates (the formulas and the group-velocity correction of
the Yee grid are in `examples/wakefield/theory.py`):

| Gate | Parameters | Criterion |
|---|---|---|
| `wakefield:laser_linear` | a0 = 0.1, k_p L = √2, k_p w = 4, ω0 / ω_p = 10; 3200×128², 8.4e8 particles, about 4 min | amplitude of E_x on the axis 5 %, wavelength 2 %, E_r / E_x at r = w / 2 10 % |
| `wakefield:beam_linear` | n_b / n_0 = 0.05, k_p σ_x = k_p σ_r = 1, rigid bunch with γ ≈ 10; 900×192², about 2 min | amplitude 5 %, wavelength 2 %, R(σ_r) / R(0) 10 % |

### Stage 4 — moving window (PIC-D22, OUT-D5)

`window: {start, margin}`. Positions stay in the laboratory frame. Under `#ifdef MOVING_WINDOW`
the kernels subtract a run-time integer shift from the index coordinate (exact) and use run-time
x bounds, computed once in Python as doubles (the same values that the migration classifier
gets). The shift comes after the second `push_B`: E, B and the twelve psi arrays move by one plane
(scratch array and swap), a new x-plane exchange of B and psiB, the kernel `window_drop`,
`_migrate`, the load of the new laboratory plane (ids in a block per shift, with the 2^53 check),
the energy term `window_out`. x− is open, x+ is reflecting behind a vacuum margin. Refused: the
implicit solver, a periodic x, injection, a body, the dipole, `B0` / `E0`, `compact`. Output:
`gridGlobalOffset` per iteration, slices relative to the window, `window_shift` in `scalars.h5`
and `run.json`. `check_kernels.py`: a window split along x, the window with the mixed gather.

| Test | Criterion |
|---|---|
| W0 | without a window: the bits of `mgpu_regress` and the assembled kernel source do not change |
| W1 | thermal plasma, 400 shifts: density 2.5 %, T_e ≤ 1.05, mean velocity per x slab < 0.02 v_te; Gauss, div B |
| W2 | the laser run with the window against the static box: after 200 shifts E agrees within 1e-6 max\|E\| in the region farther than 32 cells from x− (bit equality is not expected) |
| W3 | 1 against (2,2,1) against (1,2,2) GPUs: history ≤ 1e-10, the loaded planes equal by (id, position) |
| W4 | a pulse leaving through x−: reflected energy ≤ 1e-4; reported to the owner if it fails |

### Stage 5 — window runs

The gate configurations become window runs (laser 1152×128², about 1.5 min). Production (a0 = 4,
1536×320², 2.5e9 particles, 19,000 steps, about 80 min, 25 GB per GPU): bubble radius, peak E_x,
trapped charge and spectrum are reported (ANA-D3). Wiki page, README, roadmap R8.

## Invariants

- Charge conservation stays with the deposition: the antenna current is a discrete curl, the
  bunch is neutral at load, the window subtracts an integer from the index coordinate. No Poisson
  correction is added.
- Loading stays a function of (seed, global laboratory cell, slot) only; the planes loaded by the
  window keep 1 against N GPUs ≤ 1e-10.
- Every new collective (the shift exchange of the window, the number of dropped particles, the
  work of the antenna) is reached by every rank.
- The new `#define` (`MOVING_WINDOW`) goes into `pypic/kernelsrc.py`, with a case in
  `scripts/check_kernels.py`.
- Without `lasers` and `window` the bits of the existing runs do not change (W0).
- `pypic/` stays experiment-agnostic: the wake theory lives in `examples/wakefield/`.

## Verification

- `ruff check .`, `python scripts/check_kernels.py`, `pytest tests/cpu`.
- Targeted GPU runs of the new tests and configurations while developing; one full gate run
  (`pytest tests/gpu`, `python scripts/gpu_gates.py --gpus 8`, then `--precision mixed`) on the
  clean final head of the pull request, record and `index.json` committed, `PROGRESS.md` updated.
- Production run: `scripts/record_example.py`; results in the README of the example and the wiki.

## Out of scope

An initial field from a Poisson solve; the window with injection, a body, the implicit solver or
a speed other than c; the window with `precision: compact`; boosted frame; ionisation.

## Results

(Filled in stage by stage.)
