# Wakefield: laser and beam drivers, the moving window, `examples/wakefield`

Status: Done (PR #42) — PIC-D20, PIC-D21, PIC-D22 and OUT-D5 Accepted with the gate records on d523e3b in FP64 and with `precision: mixed`; the one pull request of R3, R5, R6 and R8 (this branch with PR #43 and PR #41 merged into it, and the review fixes) is gated on 300c1a4 in `fp64`, `mixed` and `compact`, ALL PASSED (Results, last section); PIC-D23 Proposed
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

Amendment (2026-11-12) — the particle places of the gate runs. New fact: with the configuration
as fixed above both linear gates fail on the amplitude (−6.90 % and −8.28 % against 5 %), and the
test fixed in the Results (the same `beam_linear` with the particles of a cell on 2×2×2 lattice
places instead of random places, same number of particles, same cost) gives +0.51 %: the random
places of the load are the cause of the loss inside the gate window. What changes: the four runs
of this stage (`laser_linear`, `beam_linear`, N1, N2) load their particles on a lattice of 2×2×2
places per cell. What does not change: every limit, every measurement and every other parameter
of the table. The alternative that needs no new code, 32 particles per cell with random places,
gives −2.78 % at four times the memory and time (2.1e9 particles for `beam_linear`), and the laser
gate would need 3.0e9; it is not taken.

This needs an engine option, which is a design decision of its own: **PIC-D23** (Proposed; the
owner can veto it, the gates then stay failed or go to 32 per cell):

- `SimConfig.lattice: [a, b, c]` (default: none = the random places of today, bit for bit).
  `a b c` must equal `ppc`. The slot `s` of a cell sits at `((s mod a) + 0.5) / a`,
  `(((s / a) mod b) + 0.5) / b`, `((s / (a b)) + 0.5) / c` of the cell. Momenta are drawn as
  today (seed, stream, global cell, slot).
- The places depend only on (global cell, slot), so the load stays independent of the rank layout
  and species loaded on the same places stay neutral place by place (PIC-D21).
- In the kernel: `#ifdef LOAD_LATTICE` in `load_particles` with the three numbers as defines
  (`pypic/kernelsrc.py`, a case in `scripts/check_kernels.py`); without the define the source of
  the kernel is the one of today. The planes that the moving window loads use the same kernel.
- Not changed: the particles injected at an open face keep their random places (a lattice has no
  meaning for a flux through a face).
- GPU test: the charge density of one species of a uniform beam loaded on a lattice with an even
  number of places per axis is uniform to round-off (with random places it scatters by
  1 / sqrt(particles per cell)); the places agree between 1 rank and 2×2 ranks.
- Stage order: the commit of the `#ifdef MOVING_WINDOW` macros (Stage 4, W0) comes first, so the
  comparison of the compiled kernels against 477223a is taken before `load_particles` changes.

Formulas fixed before the first run of the gates (2026-11-11; `examples/wakefield/theory.py`, the
limits are CPU tests in `tests/cpu/test_wakefield.py`):

- Pulse (`laser_wake`): `phi'' + phi = a_env² / 4` point by point, so behind a pulse with the field
  envelope `a0 exp(-(t / tau)²) exp(-r² / w²)` the amplitude of E_x is
  `a0² (tau / 4) sqrt(pi / 2) exp(-tau² / 8) exp(-2 r² / w²)` and E_r / E_x is `4 r v_g / w²`.
  Compared with Schroeder, Benedetti, Esarey and Leemans, "Nonlinear pulse propagation and phase
  velocity of laser-driven plasma waves" (2011, OSTI 1013067): their Eq. (6),
  `x'' = (gamma_perp² / x² - 1) / 2` with `gamma_perp² = 1 + a² / 2`, is this equation for
  `x = 1 + phi` to first order, and their resonant Gaussian pulse, `k_p L = 1` with
  `a = a0 exp(-zeta² / (4 L²))`, is `tau = 2`, where `laser_wake` is largest. The `k_p L = sqrt(2)`
  of the table is the same pulse written as `a ~ exp(-zeta² / (2 L²))`: the antenna gets
  `duration = 2`. The radial field (the same equation at every radius, `E_r = -d phi / dr`) is
  derived here and was not compared with a paper.
- Bunch (`beam_wake`): `E_x = -q n_b sqrt(2 pi) sigma_x exp(-k_w² sigma_x² / 2) R(r) cos(k_w xi)`
  with `R(r) = k_w² int r' n_perp(r') I_0(k_w r_<) K_0(k_w r_>) dr'`. Compared with Muggli,
  "Beam-Driven Systems, Plasma Wakefield Acceleration" (CERN Accelerator School 2019,
  arXiv:2007.05226), Eqs. (1) and (3): the same field and the same R(r) for beta = 1. Here
  `k_w = 1 / beta` in the cosine and in R (derived here from the pole of 1 / eps): 0.5 % in k at
  gamma = 10.
- Predicted for the parameters of the gates: pulse, E_x on the axis at the focus 3.801e-3, E_r / E_x
  at r = w / 2 `0.5 v_g`; bunch, E_x on the axis 0.03505, `R(sigma_r) / R(0)` = 0.7516,
  `k_w` = 1.00504. Wavelengths `2 pi v_g` and `2 pi beta`.
- New fact for the pulse gate, diffraction: with `k_p w = 4` and `omega0 / omega_p = 10` the Rayleigh
  length is 80 c/omega_p and the box is 63 c/omega_p long, so `a0²` on the axis changes by up to
  13 % along the box. The wake follows the local `a0(x)²` and `w(x)` (`beam_at`), and the gate
  compares with the theory at the place of the measurement. The limits of the table stay.
- Group velocity on the grid (`yee_group_velocity`): for `dx = lambda0 / 32`, `dt = dx / 2` it is
  0.99142 (0.99499 without the grid), 0.86 % below c against the 2 % of the wavelength criterion.
- Correction of the pulse formula (2026-11-11, found while writing the example, before the first
  run): the first item gives the amplitude of phi and of E_x for a pattern that moves with c. With
  the pattern `phi(t - x / v_g)` the fields are `E_x = phi / v_g` and `E_r = (4 r / w²) phi`
  (`laser_wake` had `phi` and `(4 r v_g / w²) phi`: both low by the factor v_g, the ratio right);
  CPU test: E_x from its own equation, `E_x'' + E_x = (1 / v) d/dt (a² / 4)`. And a0 is the
  amplitude inside the plasma: through the ramp the pulse keeps its energy flux, `a² = a0² / eta`,
  `eta = sqrt(1 - omega_p² / omega0²)` (`amplitude_in_plasma`; CPU test: the wave equation through
  a ramp of 20 wavelengths). Predicted E_x on the axis at the focus: 3.853e-3 (3.801e-3 above:
  +0.87 % from 1 / v_g, +0.50 % from eta); E_r / E_x at r = w / 2 stays `0.5 v_g` = 0.4957. Neither
  term was taken from a paper.

Set-up and measurement of the gates, fixed before the first run (2026-11-11). Both runs are
explicit, with electrons and ions (mass 1836) at 8 particles per cell each, open x faces, and a
density that rises from zero to n0 along x over a ramp (`table` profile) and falls again before
the far face, so no plasma touches an open face.

- `laser_linear`: 3200×128² cells, `dx = lambda0 / 32`, `dy = dz = 24 / 128`, `dt = dx / 2`;
  transverse faces periodic. Antenna in the vacuum before the ramp, `t_peak = 4 tau`, focus at the
  centre of the box; the run ends when the centre of the pulse is at `x_p` = 54 (the travel time
  is the integral of 1 / v_g over the local density, with `yee_group_velocity`). Measured on the
  slice through the axis that is normal to the polarisation (E_z is the radial field there, the
  laser field has no component in it), at the last output, in the window from one c/omega_p
  behind the start of the plateau to `x_p - 3 tau`.
  - Amplitude and wavelength: least-squares fit of the on-axis `E_x(x) = c A_th(x) cos(k x + phi)
    + b`, with `A_th(x)` the amplitude of the theory at x (`beam_at` for the local a0 and w,
    `amplitude_in_plasma`, `laser_wake` with the speed v_g); free c, k, phi, b; k starts at the
    peak of the spectrum. Criteria: `|c - 1|` ≤ 5 %, `|k v_g - 1|` ≤ 2 %.
  - E_r / E_x at r = w / 2: E_x and E_z interpolated to the two lines at the distance w / 2 from
    the axis (E_z taken with its sign away from the axis, the two sides averaged); with k and phi
    of the on-axis fit, `E_x = c_x A_x,th(x) cos(k x + phi) + b` and
    `E_r = -c_r A_r,th(x) sin(k x + phi) + b'`. Criterion: `|c_r / c_x - 1|` ≤ 10 %.
- `beam_linear`: 900×192² cells, `dx = 1 / 8` in all directions (sigma = 8 cells), `dt = dx / 2`;
  transverse faces reflecting: the partner leaves a net charge that oscillates at its place, and
  between periodic transverse faces that charge gives a uniform E_x of about 2 % of the wake,
  which the conducting walls remove. The bunch (q = -1, mass 1e6, gamma = 10) and its partner
  (q = +1, mass 1e6, at rest) start inside the plateau at x = 12, centred on the axis; the run
  ends when the bunch is at `x_b` = 100. Window: from `max(42, x_b - 5 sigma_x - 5 * 2 pi beta)`
  to `x_b - 5 sigma_x` (five wavelengths, at least 30 c/omega_p from the start of the bunch).
  - Amplitude and wavelength: fit of the on-axis `E_x = c A_th cos(k (x - x_b) + phi) + b`,
    `A_th = beam_wake(r = 0)`. Criteria: `|c - 1|` ≤ 5 %, `|k beta - 1|` ≤ 2 %.
  - `R(sigma_r) / R(0)`: the amplitude of E_x on the lines at the distance sigma_r from the axis
    in the slice (8 cells, on nodes; both sides averaged) with k and phi of the on-axis fit, over
    the on-axis amplitude, against 0.7516. Criterion: 10 %.
  - Reported, not a criterion: the on-axis amplitude at the output before the last (is the wake
    stationary?), the phase phi.
- N1 is the box of `beam_linear` and N2 the box of `laser_linear` (cells with dx ≠ dy), each run
  once with `driver: none` and recorded here: max |E_x| on the axis ≤ 1 % of the predicted
  amplitude of that gate, mean electron density over the plateau within 1 % of n0. A cold plasma
  loaded with its ions at the same places has no field at all, so the expected result is zero up
  to round-off; the runs check that the ramps, the faces and the anisotropic cells add nothing.
- Both gates also need Gauss ≤ 1e-10 in the checked region (the gate driver's criterion).
- The time estimates of the table are replaced by the plans of the configurations: about 6,200
  steps for `laser_linear` and 1,420 for `beam_linear`.
- Notes written with the example (2026-11-11, before the first run; no limit and no measured
  quantity changes):
  - The null runs are the gate configurations with `--set drive=false`, not a third value of
    `driver`: the box, the run time and the predicted amplitude that the 1 % refers to are those
    of the gate, and `driver` names it. Max |E_x| is taken over all outputs.
  - The mean electron density is taken on the plane through the axis, from one c/omega_p inside
    each end of the plateau, without the two nodes next to each transverse face (a node on a
    wall holds half a cell).
  - `A_th(x)` uses the Rayleigh length of the grid in vacuum. Inside the plasma the wavenumber is
    0.5 % smaller, which moves the focus by 0.13 c/omega_p and changes `a0²` on the axis by
    less than 0.2 % in the window; not included.
  - The measurement is checked on planted wakes before the first run
    (`tests/cpu/test_example_analysis.py`): a wake with 1.03 times the amplitude, 1.004 times the
    wavenumber and 0.95 times the radial field of the theory returns these three numbers.

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

Amendment (2026-11-11, before any Stage 4 code) — W0. The kernel files get `#ifdef
MOVING_WINDOW` branches and macros for the x bounds, so the assembled text of every configuration
changes and the second half of W0 cannot hold as written. It is replaced by what it stood for:
without a window the macros expand to the tokens that are there now, and the code that NVRTC
compiles from the cases of `scripts/check_kernels.py` (the PTX) equals the one of 477223a. This
is compared once, with a script outside the repository; the outcome goes into Results. The first
half (the bits of `mgpu_regress`) stays as it is.

Design fixed before coding (the reasons are in PIC-D22):

- Shift: after a step that ends at the time t the window has moved by
  `s = max(0, floor((t - start) / dx))` cells, at most one more than before (dt < dx). Every rank
  computes s from the configuration and the step number.
- Values of the run: the particle module holds `window_state` (s, the two faces of the window
  along x, the two bounds of the rank along x), written by the kernel `window_set`; no kernel gets
  a new argument. The Python doubles are the ones that the migration classifier gets. Two
  simulations of one process must not share the module: the defines carry a number of the
  instance.
- Exactness: the index coordinate along x is `(x - OX) / DX - s`. For 0 <= s <= X < 2^53 the
  difference is exact, so a shift lowers the index coordinate of every particle by exactly 1 and
  the charge on the grid moves with the fields.
- Fields: E, B and the psi arrays are copied one plane down into a scratch array, which then
  takes the place of the array. Before the copy psiE is exchanged along x (its ghost planes are
  not kept up to date otherwise); after it E, B and psiB are exchanged along x. This is wider
  than "B and psiB" above: the plane that enters a rank from its neighbour is needed for all of
  them.
- The x+ wall moves one cell: its surface charge moves with it. The current that carries it, on
  the edges between the old and the new wall node, changes E_x there by the value of E_x on the
  edge before, so `E_x(W - 1) = E_x(W - 2)`; every other field of the new cell is zero. Gauss's
  law at the node that stops being a wall then holds by the continuity of this current. No
  Poisson solve, no correction of a field that is already there.
- x−: after the shift `window_drop` is the only kernel that fills the flag list of this step: dead
  slots, particles below the new face (their kinetic energy, charge and number are summed, their
  weight is set to 0, nothing is deposited) and particles outside the bounds of the rank. The
  charge of a dropped particle stays on the nodes -1 .. reach of the new window, which the Gauss
  check leaves out at an open face already (PIC-D15).
- New plane: the plasma is loaded up to the window cell `n - 1 - m`, `m = round(margin / dx)`;
  after the shift to s the laboratory plane `s + n - 1 - m` is loaded, by the loader of the
  initial load (same seed, global cell and slot), ids from a block per shift. A beam without a
  drift along x is loaded as it was at the shift before and moved for the time since then, with a
  reflection at the upper face of the plane: with this a thermal plasma is a Maxwellian behind
  the front at every load (a plane that is loaded where it stands counts the particles that
  crossed into it from the plane before twice, a mean velocity of v_te^2 / c per slab). A beam
  with a drift along x is loaded where it stands. No half kick: the fields at the front are those
  of the vacuum margin.
- Energy: `window_out` is the sum of the field energy that the copy removes, the kinetic energy
  of the dropped particles and minus the kinetic energy of the loaded planes. It joins the one
  reduction of `_record`; the recorded value is the one at the time of the record (before the
  shift of the step), like the work of the antenna.
- Antenna: its two rows follow the shift. It is switched off when its cell reaches the nodes
  next to x− that the Gauss check leaves out; the configuration is refused unless the envelope
  is Gaussian and has fallen to exp(-25) by then.
- Refused with a window, besides the list above: `initial_fields`, `external_fields`, a
  `perturbation` of a beam, x boundaries other than (`open`, `reflecting`), a margin of fewer
  than 4 cells plus the reach of the current filter, or one that leaves no plasma.

Amendment (2026-11-12, before the kernels of Stage 4) — the load of a plane, the energy term,
the x+ wall. Three facts that the design above did not have:

1. The beams of a plane are loaded on the same positions (this is what makes the load neutral),
   and the move "for the time since the shift before" gives each species its own displacement.
   After it the charge on the grid is not zero, and E holds none of it: Gauss's law would fail at
   checked nodes. So the move is a move of the engine like any other: it deposits its current
   (Esirkepov, two halves so that no segment is longer than a cell; the reflection at the upper
   face of the plane splits a segment there), the current is reduced, filtered and gets the wall
   condition like the current of a step, and `E -= dt J`. The move is at most `dx + dt` long, less
   than two cells: with the ghost width 2 the current of a particle of the first cell of a rank
   stays inside the arrays and the particle migrates to the adjacent rank.
2. "A beam with a drift along x is loaded where it stands" is not the state of such a beam. The
   particles of the plane before have moved into the new plane (or away from it) by the drift
   since the shift before: the new plane would have to be loaded without the part they fill, or
   with a part of the plane ahead. Both take or add particles of one species only at places
   that the other species of the plane keep, which ends the neutrality position by position that
   Gauss's law rests on here; and the profile of such a beam is not at its place of t = 0 when
   the plane enters. No run of R8 needs it (the bunch and its partner lie in the cells of the
   first load). The line is replaced by a refusal, below.
3. The field energy that a shift takes out of the physical cells is not only the plane that the
   copy moves into the absorbing layer: the move of fact 1 and the rule of the wall change E too.

What changes in the design:

- Beams: every beam has a range of laboratory cells along x that it loads, from its profile on
  the cells `0 .. n - m + s_end` (`s_end`: the shift after `n_steps` steps; no profile along x:
  all of them). The first load takes the cells below `n - m`, the shift to s the plane
  `s + n - 1 - m`. A beam whose range ends at or below the cell `n - m` is loaded once and may be
  of any kind (a drifting bunch, a radial beam). A beam that enters later planes must have no
  drift along x, no radial drift and no radial profile: refused otherwise. A run that goes on
  beyond `n_steps` stops with an error at the shift `s_end + 1` (a function of the configuration
  and of the step, so every rank raises).
- A beam with a temperature that enters later planes: the plane is loaded as at the shift before
  and moved as in fact 1. This is the exact state of a plasma that is uniform across the plane
  and symmetric in v_x, except for the particles faster than `dx / T` along x (T: the time since
  the shift before, at most `dx + dt`; with dt = dx / 2 that is 2 c / 3), which the reflection at
  the upper face of the plane would have to hold back in the plane before as well. With a density
  that varies across the plane it is as good as the first load is. For the first shift the time
  is counted from `start`: a thermal plasma that waits for the window (start > 0) has expanded
  into the margin by then, as a plasma with a free surface does. W1 uses `start = 0`. A cold beam
  is loaded where it stands, without the move.
- Neutrality of a plane: the same check as for the first load (net charge position by position
  over the beams with the same cells across the plane), reduced over the ranks at every shift
  that loads a plane; every rank raises if it fails.
- Ids of a plane: `n_ids + ((p - (n - m)) n_streams + stream - 1) ny nz ppc` plus the number of
  the particle in the plane, p the laboratory plane. Checked once at the start for the last
  plane: below 2^53.
- Energy: `window_out` is (field energy + kinetic energy of the particles) before the shift minus
  after it: the field energy of the physical cells before the copy minus the one after all that
  the shift does to the fields (copy, wall rule, the E of the move), plus the kinetic energy of
  the dropped particles, minus the kinetic energy of the loaded ones.
- x+ wall: the rule `E_x(W - 1) = E_x(W - 2)` is the rule for a wall node without the charge of a
  particle. A particle within 1 + (reach of the current filter) cells of the x+ wall at a shift
  is an error of the configuration (the driver or the plasma reached the wall before the window
  left): `window_drop` counts them, the count is reduced with every record and at the end of a
  run, and every rank raises. The margin of at least 4 cells plus the reach keeps the loaded
  plasma away; nothing that starts behind the front catches up with a wall that moves at c.

Amendment (2026-11-12, with the kernels of Stage 4, before the shift is coded) — the move of a
plane goes in rounds. New fact: the arrays of a rank hold the nodes −2 .. n + 1 of each axis
(ghost width 2). A particle of the cell c reaches the nodes c and c + 1, so the deposit of a
push, which starts in an owned cell and is shorter than a cell, ends at the node n + 1 at most.
The move of fact 1 above is up to `dx + dt` long. Along x upwards the upper face of the plane
bounds it, and downwards two ghost nodes are enough; but across the plane a particle of the last
cell of a rank can go on for more than a cell (1.5 cells with cubic cells and dt = dx / 2) and
its second half would be deposited on the node n + 2, outside the array. The sentence "with the
ghost width 2 the current ... stays inside the arrays" of fact 1 holds for the downward
directions only. What changes:

- The move is done in `R = floor(T / min(dx, dy, dz)) + 1` rounds of the time `T / R` each (two
  rounds in the runs of R8), every round shorter than a cell on every axis, with a migration
  after each round: the way the initial perturbation is applied (PIC-D9,
  `Simulation._apply_perturbations`). The current of all rounds is collected and then reduced,
  filtered and applied once, as fixed in fact 1. The kernel is `window_move`.
- Order within a shift, so that the particles of the new plane are the last ones of a rank's
  arrays while they are moved: fields (copy, wall rule, exchanges) → `window_set` → `window_drop`
  on the particles of before the shift and their migration → load of the plane → the rounds of
  the move, each followed by a migration (`migrate` appends what stays and what arrives at the
  end of the arrays, so the plane stays one block) → the current of the move into E.
- The record of a step is taken before the shift at its end: the history of the step n holds the
  shifts up to the one that made the state at t^n, like the fields and the particles it records.
- `window_drop` and `window_move` sum into a buffer of their own per species (kinetic energy,
  charge and number of the dropped particles; the count of charged particles within the reach of
  the x+ wall; kinetic energy of the loaded particles), reduced with the record. A particle of a
  plane that leaves through an open face across the window during its move is counted as loaded
  and as lost, like a particle lost in a push.

What does not change: the state that a plane is loaded in, the reflection at its upper face, the
criteria of W1 to W4.

### Stage 5 — window runs

The gate configurations become window runs (laser 1152×128², about 1.5 min). Production (a0 = 4,
1536×320², 2.5e9 particles, 19,000 steps, about 80 min, 25 GB per GPU): bubble radius, peak E_x,
trapped charge and spectrum are reported (ANA-D3). Wiki page, README, roadmap R8.

Design (2026-11-12, fixed before the first window run of the example).

- **Planner:** `window_lead` (a length; none: the box at rest of Stage 3) and `window_margin`.
  With it the x faces are (open, reflecting) and the window starts at the time at which the driver
  is `window_lead` behind the x+ wall, so the driver keeps that place in the window (it falls back
  by (1 − v) per unit of path). The ramp, the focus, the start of the bunch, `x_end` and the run
  time are laboratory values as before.
- **Gate configurations** (same names, same driver, plasma, cells, lattice, limits):
  `laser_linear` 1152×128² (22.6 c/ω_pe; the focus stays at x = 10π, the centre of the box of
  Stage 3), `window_lead` 8; `beam_linear` 384×192² (48 c/ω_pe), `window_lead` 8.
- **Measurement:** as in Stage 3, on the laboratory x of the output (OUT-D5), with one rule added
  for every run: the interval of the fit and of the density starts no closer than 1 c/ω_pe to the
  lower face of the box (for a box at rest this changes nothing: the intervals start at x = 9 and
  beyond 42). The density of the plateau is the mean up to 1 c/ω_pe behind the front of the
  plasma. Consequence: the laser fit has 7 to 8 c/ω_pe (1.2 wavelengths) instead of 39; the bunch
  fit keeps its five wavelengths.
- **Check of the short interval, done before the window run** on the Stage 3 run in the box at
  rest (68a3a29, last output): the fit on the last 7.3 c/ω_pe before the pulse gives the amplitude
  −1.33 % (whole plateau: −1.91 %), the wavelength −0.02 % (+0.07 %), E_r / E_x +0.06 % (+0.04 %);
  on 8.3, 12 and 16 c/ω_pe the same within 0.02 % (amplitude), 0.01 % (wavelength) and 0.2 %
  (ratio). The short interval measures the same wake.
- **Criteria of the four window runs:** the gate limits (amplitude 5 %, wavelength 2 %, radial
  10 %; without a driver |E_x| ≤ 1 % and density within 1 %), Gauss ≤ 1e-10. Reported: the
  difference to the fit of the Stage 3 run on the same laboratory interval.
- **Production** (`laser_bubble`; reported, ANA-D3): a0 = 4, ω0 = 10 ω_pe, duration 2, waist 4
  focused at the top of the ramp (x = 8), 1536×320² cells (30.2 × 40² c/ω_pe, dx = λ0 / 32,
  dy = dz = 0.125), periodic across, 8 per cell on the 2 × 2 × 2 lattice, ions of the mass 1836,
  `window_lead` 8, to `x_end` = 180 (about 19,000 steps). Outputs: the plane through the axis
  (E_x, E_y, E_z, electron density), the electron histograms (x, u_x) and u_x. Reported at the
  last output: the radius of the cavity (largest half-width across of the region with electron
  density below 0.5 n0 behind the pulse), max and min of E_x on the axis behind the pulse and its
  slope at the zero inside the cavity (Lu et al. 2007: k_p R ≈ 2 √a0 = 4, slope 1 / 2), the charge
  and the largest u_x of the electrons with u_x > 20 (twice the γ of the wake), the spectrum. A
  200-step probe first (ms per step, device memory); an estimate above 3 h is reported instead of
  run.

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

### Stage 1 — antenna, test A1 (2026-11-11): the criteria of the table hold

`tests/gpu/test_engine.py::test_laser_antenna_emits_the_paraxial_beam` (336×96², dx = λ/12,
dy = dz = λ/6, dt = λ/16, six open faces, w0 = 2λ, the plane 12.6 λ before the focus, steady
beam; one GPU, FP64). E along y and E along z give the same numbers to 15 digits.

| Criterion | Limit | Measured |
|---|---|---|
| Width at the focus; 1.05 z_R behind it | 2 % | +0.13 %; +0.74 % |
| Amplitude at the focus; 1.05 z_R behind it | 2 % | −0.21 %; −1.44 % |
| Forward flux against backward flux | 1 % | −0.076 % |
| Work of the antenna against the flux through the two planes | 1 % | −0.005 % |
| Gauss, div B | 1e-10, 1e-12 | 1.5e-14, 1.7e-14 |
| 1 GPU against 2×2 GPUs (2,2,1 for E along y, 2,1,2 for E along z), plane on a rank face | history ≤ 1e-10 | passed |

A check that is not in the table, written into the test before its first run, failed in that
run: `|budget| < 1 % of the work of the antenna` after 35 periods gave 1.33 % (645.1 of 48,325).
It is not a leak. The budget of the explicit scheme (monitored only) keeps an offset that does
not grow: from 35 to 80 periods the work rises from 48,281 to 114,650 and the budget goes from
645 to 611. Of the 645, 447 (0.92 % of the work) are the definition of the recorded magnetic
energy: the record takes B at t^n, and `B_n²` is above the `B^(n−1/2) · B^(n+1/2)` of the
leapfrog energy by `((dt / 2) curl E)²`, which is `tan²(ω dt / 2)` = 4.0 % of the magnetic energy
of a wave with ω dt = 2π / 16 (computed in the run from `curl E`, not estimated). The other 199
(0.41 %) appear during the first 10 periods, while the beam is switched on and the flux through
the open faces is not yet steady, and stay between 141 and 199 afterwards. The check was replaced
by the statement that does hold exactly in the steady state: the budget does not change over one
period once the box is full, within 1 % of the work of that period (measured: 5.8e-5). It
protects the sign and the place of the term `antenna` in the budget.

### Stage 2 — bunch, test B1 (2026-11-11): failed as first defined; amended measurement passes

First run, as fixed in the test before the run (448×192², σ = 8 cells, γ = 10, n_b = 0.05,
8 particles per cell, 640 steps = 40 σ; E_r and B_θ as the mean of the four lines ±y, ±z through
the axis, at the radii (m + 1/2) dx, m = 4…15; limits 3 % for E_r against the steady field and
1 % for B_θ against β E_r): **failed**. `E_r / E_steady − 1` = +5.62, +6.07, +4.48, +0.70,
+0.12, −0.20, +1.14, +0.89, +1.41, +1.60, +1.67, +1.70 %; `B_θ / (β E_r) − 1` up to 1.97 %.
Gauss 1.7e-15, div B 6.6e-16.

Cause: the shot noise of the bunch particles in the measurement, which the prediction of the
start-up (+0.83 % at the smallest radius to +1.63 % at the largest, B / βE + 0.07 to 0.15 %) did
not contain. Shown by the same run with other seeds and particle numbers (one GPU each):

| Particles per cell | Seed | Four lines: max E dev. | Circle mean: max E dev. | Circle mean: max B dev. |
|---|---|---|---|---|
| 8 | 12345 | 6.07 % | 3.67 % | 0.79 % |
| 8 | 2024 | 3.63 % | 3.95 % | 0.92 % |
| 8 | 777 | 2.10 % | 3.31 % | 0.68 % |
| 27 | 12345 | — | 1.87 % | 0.36 % |
| 27 | 777 | — | 2.26 % | 0.63 % |
| 64 | 12345 | 2.07 % | 1.77 % | 0.79 % |
| 64 | 2024 | — | 1.93 % | 0.37 % |
| 64 | 777 | 2.28 % | 1.93 % | 0.43 % |

("Circle mean": the mean over 64 points of the circle of each radius, bilinear in the plane.
With 8 particles per cell single lines differ from each other by up to 14 %. With 64 the four
lines still give a B deviation of 1.78 % for the seed 12345.) With 64 particles per cell the
largest deviation of E is at the largest radius for every seed, 1.77 to 1.93 % against the
1.63 % of the start-up.

**Amendment (2026-11-11) to the measurement of B1.** E_r and B_θ are the means over the circle
of each radius (64 points, bilinear in the mid-plane), and the bunch is loaded with 64 particles
per cell. Reason: the criterion is about the field of the bunch; with 8 particles per cell on
four lines the measurement is limited by the shot noise of the bunch (table above), not by that
field. Not changed: the radii, the run length, the seed, and the limits 3 % and 1 %.

Result with the amended measurement (seed 12345, one GPU, 55 s): `E_r / E_steady − 1` = +0.91,
+1.05, +1.13, +1.09, +0.86, +0.96, +1.13, +1.34, +1.61, +1.67, +1.70, +1.77 % (limit 3 %; the
start-up alone gives +0.83 to +1.63 %); `B_θ / (β E_r) − 1` = +0.79 % at the smallest radius,
within ±0.23 % at the others (limit 1 %); Gauss 1.6e-15, div B 5.6e-16. Passed.
`test_bunch_is_decomposition_independent` (1 GPU against 2×2): passed.

Not shown: the 0.79 % of B at the smallest radius is not the start-up (0.07 %) and differs from
seed to seed (0.37 to 0.79 % with 64 particles per cell); it was not followed further. An initial
field from a Poisson solve stays out of scope (B1 does not fail at the run-up within reach).

### Stage 3 — static-box gates (2026-11-11/12): N1 and N2 pass; both linear gates fail on the amplitude

First runs of the configurations as fixed above (8 GPUs, FP64, 8 particles per cell; commit
848517e for `beam_linear`, 3823577 for the others — the two differ only in the progress line of
`Simulation.run`).

| Run | Criterion | Limit | Measured | |
|---|---|---|---|---|
| N1 (`beam_linear`, `drive=false`) | max \|E_x\| on the axis / predicted; mean plateau density | 1 %; 1 % | 0.0; +0.013 % | passed |
| N2 (`laser_linear`, `drive=false`) | the same | 1 %; 1 % | 0.0; +0.006 % | passed |
| `beam_linear` | amplitude `c − 1` | 5 % | **−6.90 %** | **failed** |
| | wavelength `k beta − 1` | 2 % | +0.59 % | |
| | `R(sigma_r) / R(0)` against 0.7516 | 10 % | −0.44 % | |
| `laser_linear` | amplitude `c − 1` | 5 % | **−8.28 %** | **failed** |
| | wavelength `k v_g − 1` | 2 % | +1.15 % | |
| | `E_r / E_x` at r = w / 2, `c_r / c_x − 1` | 10 % | −0.22 % | |

Gauss: 7.8e-15 (`beam_linear`), 8.5e-13 (`laser_linear`), 2.6e-15 (both null runs). Sizes:
`beam_linear` 900×192², 5.3e8 particles, 1,416 steps, 52 s, 11.4 GB per GPU at the peak;
`laser_linear` 3200×128², 7.6e8 particles, 6,160 steps, 297 s, 15.0 GB per GPU. The driver is
where the plan puts it (bunch centre 100.02 against 100.06, pulse centre 53.98 against 54.04).
Both gates stay failed; no limit and no definition is changed here.

What fails: the wake is right behind the driver and loses amplitude with the distance behind
it. Fits as in the gate, but over one wavelength each (`amp / theory`; `beam_linear` at
t = 88.5, wavelengths counted from `x_b − 5 sigma_x`; the gate window is the first five):

| Wavelength behind the driver | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `beam_linear` (8 per cell, dx = 0.125, width 24) | 0.986 | 0.967 | 0.934 | 0.904 | 0.855 | 0.830 | 0.807 | 0.793 | 0.758 | 0.760 | 0.713 | 0.712 |
| 32 per cell (`ppc=32`) | 1.007 | 0.999 | 0.981 | 0.957 | 0.914 | 0.887 | 0.855 | 0.863 | 0.817 | 0.825 | 0.784 | 0.774 |
| width 48 (`width=48 n_perp=384`) | 0.987 | 0.965 | 0.941 | 0.888 | 0.866 | 0.830 | 0.817 | 0.768 | 0.759 | 0.751 | 0.715 | 0.717 |
| dx = 0.25, 256 per cell (`dx=0.25 nx=450 n_perp=96 ppc=256`, t = 89.0) | 1.000 | 0.997 | 0.991 | 0.976 | 0.967 | 0.940 | 0.933 | 0.896 | 0.813 | 0.789 | 0.728 | 0.692 |

`laser_linear` at t = 60.5 (one wavelength each, the delay behind the pulse centre in 1 / omega_p):
0.962 at 9.2, 0.957 at 15.5, 0.939 at 21.8, 0.919 at 28.0, 0.898 at 34.3, 0.872 at 40.6.

Diagnostic runs (not gates; each one changes the `--set` values named in the table, 8 GPUs):

| Run | Particles | Gate fit: amplitude, `k beta`, `R` ratio, rms residual | Gauss | Wall, peak per GPU |
|---|---|---|---|---|
| `beam_linear` | 5.3e8 | 0.9310, 1.0059, 0.7483, 0.224 | 7.8e-15 | 52 s, 11.4 GB |
| 32 per cell | 2.12e9 | 0.9722, 1.0016, 0.7491, 0.191 | 2.9e-14 | 197 s, 42.4 GB |
| width 48 | 2.06e9 | 0.9310, 1.0058, 0.7493, 0.223 | 9.1e-15 | 191 s, 38.4 GB |
| dx = 0.25, 256 per cell | 2.12e9 | 0.9824, 0.9912, 0.7564, 0.069 | 2.2e-14 | 139 s, 53.0 GB |

What the runs show:

- The loss follows the age of the oscillation, not the place or the time: in `beam_linear` the
  first wavelength behind the bunch has 0.986, 0.981, 0.986 and the fifth 0.865, 0.843, 0.855 at
  t = 66.4, 77.4, 88.5, although these are different places of the box.
- It is not the transverse walls: twice the width gives the same table within the scatter and
  the same gate value (0.93104 against 0.93101).
- It is not the dispersion of the grid (an error of order dx²): with the same number of particles
  per volume the cell of twice the size loses less (0.967 against 0.914 at the fifth wavelength).
- It becomes smaller with the particles per cell: `−ln(amp / theory)` at the fifth wavelength is
  0.157 with 8, 0.090 with 32 (both dx = 0.125) and 0.034 with 256 (dx = 0.25). It is not
  proportional to 1 / (particles per cell) (8 → 32: a factor 1.7, not 4).
- The last four columns of the dx = 0.25 run (x < 45) lie where the bunch field was still
  building up (start at x = 12, see B1); their residuals are 0.34 to 0.82 of the amplitude, and
  they are not used.
- The null runs have no field at all: the two species are loaded on the same places, so the
  shot noise of the load has no charge until a driver displaces the electrons.

Not shown — the cause. Working explanation: the load is random (`urand` positions in each cell),
so the electron density that the wake displaces varies from cell to cell by about
1 / sqrt(particles per cell); the oscillation then does not have one frequency, and the wave
loses its phase coherence with its age. This fits the dependence on the particles per cell and
the null runs, and it was not tested directly. Other causes that depend on the particles per
cell are not excluded by the runs above.

Next test, fixed before it is run: `beam_linear` unchanged (8 per cell, the same grid, seed and
measurement) with the particles of each cell on the same lattice of 2×2×2 places instead of
random places, as a patch that is not committed. If the load is the cause, the amplitude of the
table above stays at the value of the first wavelength for all wavelengths (within the 2 % that
the 256-per-cell run reaches), and the gate value is within the 5 %. The outcome is recorded
here before anything is changed in the engine or in a gate definition.

Outcome of that test (2026-11-12; 8 GPUs, FP64, code of 3823577 with the uncommitted patch: the
slot `s` of a cell is put at `((s & 1) + 0.5) / 2`, `(((s >> 1) & 1) + 0.5) / 2`,
`(((s >> 2) & 1) + 0.5) / 2` of the cell instead of three random numbers; momenta, weights, ids,
grid, seed and measurement unchanged):

| Run | Criterion | Limit | Random places (above) | Lattice 2×2×2 | |
|---|---|---|---|---|---|
| `beam_linear` | amplitude `c − 1` | 5 % | −6.90 % | **+0.51 %** | within the limit |
| | wavelength `k beta − 1` | 2 % | +0.59 % | −0.25 % | |
| | `R(sigma_r) / R(0)` against 0.7516 | 10 % | −0.44 % | −0.77 % | |
| | rms residual of the fit / amplitude | — | 0.224 | 0.042 | |
| the same, `drive=false` | max \|E_x\| on the axis / predicted; mean plateau density | 1 %; 1 % | 0.0; +0.013 % | 0.0; 0.000 % | |

Gauss 5.2e-15 (2.5e-16 without the driver); 52 s and 10.9 GB per GPU at the peak (random
places: 52 s, 11.4 GB). The gate value at the output before the last (t = 77.4) is +0.42 %.

One wavelength each, as in the table above (`amp / theory`, lattice load):

| Wavelength behind the driver | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| t = 66.4 | 1.013 | 1.013 | 1.016 | 1.005 | 1.003 | 1.002 | 0.952 | 0.968 | 0.937 | | | |
| t = 77.4 | 1.013 | 1.012 | 1.012 | 1.002 | 0.988 | 0.992 | 0.976 | 0.953 | 0.936 | 0.920 | 0.892 | |
| t = 88.5 | 1.013 | 1.013 | 1.011 | 1.002 | 0.992 | 0.975 | 0.977 | 0.971 | 0.927 | 0.931 | 0.878 | 0.855 |

Against what was fixed before the run:

- The gate value is within the 5 %: yes (+0.51 %).
- The amplitude stays at the value of the first wavelength for all wavelengths within 2 %: **no**.
  It does over the first four wavelengths (within 1.1 %; the fifth is 2.1 % below the first at
  t = 88.5), then it falls: 3.6 to 4.2 % below at the wavelengths 6 to 8, 8 to 16 % at 9 to 12.
- So the random places of the load cause most of the loss in the gate window (0.855 → 0.992 at
  the fifth wavelength, 0.712 → 0.855 at the twelfth) and not all of the loss far behind the
  bunch. What remains is smaller than the loss of the 256-per-cell run with random places
  (0.967 at the fifth wavelength, 0.896 at the eighth). Its cause is not identified. The last
  four columns at t = 88.5 lie at x < 45, where the dx = 0.25 run had residuals of 0.34 to 0.82;
  here the residuals are 0.04 there as everywhere, so the build-up of the bunch field does not
  explain them.
- The measurement is quieter: the rms residual of the gate fit is 0.042 instead of 0.224 of the
  amplitude.

Both linear gates stay recorded as failed in the configuration fixed above (random places, 8 per
cell). See the amendment of 2026-11-12 in Stage 3 for what follows from this outcome.

Found on the way:

- The E_x map of the analysis figure was drawn transposed (822bc9a); no number depended on it.
- A run that starts without energy (a null run) raised `ZeroDivisionError` in the progress line
  (3823577, with a regression test).
- `pypic.launch` does not notice a rank that raises and then hangs in the interpreter shutdown:
  after the error above the launcher waited without end (the chain used `timeout`). Not fixed
  here.
- The budget of `laser_linear` has a constant offset of 0.039 on a work of the antenna of 62.75
  (0.06 %): the `tan²(omega dt / 2)` = 0.24 % of the magnetic energy 15.5 that the record takes
  at t^n (see A1).

The lattice option in the engine (2026-11-12, PIC-D23, Proposed): `SimConfig.lattice`, the defines
`LATTICE_X / _Y / _Z`, the `#ifdef` in `load_particles`, one case in `scripts/check_kernels.py`.
Null check of the change: the six cases without a lattice and without a window give the same
kernel text after the C preprocessor (`cpp -P`) as before it, both units. Their PTX is the same in
23 of the 24 compilations; in one (absorbing sphere, mixed gather, `compute_80`) the number of one
virtual register differs (`%fd873` / `%fd874`, three bytes). The unchanged commit 477223a shows
the same difference between two working directories of the compiling process, so the PTX hash of
that case is not a reproducible null test; the preprocessed text is (this also qualifies the
"byte for byte" of W0 below: it held for the runs that were compared).

The four runs with the lattice places (2026-11-12, commit 68a3a29, 8 GPUs, FP64, 2×2×2 places per
cell through `lattice: [2, 2, 2]` in the two gate configurations; everything else as fixed in the
table of the stage): **all four pass.**

| Run | Criterion | Limit | Measured | |
|---|---|---|---|---|
| N1 (`beam_linear`, `drive=false`) | max \|E_x\| on the axis / predicted; mean plateau density | 1 %; 1 % | 0.0; 0.0 | passed |
| N2 (`laser_linear`, `drive=false`) | the same | 1 %; 1 % | 0.0; 0.0 | passed |
| `beam_linear` | amplitude `c − 1` | 5 % | +0.51 % | passed |
| | wavelength `k beta − 1` | 2 % | −0.25 % | |
| | `R(sigma_r) / R(0)` against 0.7516 | 10 % | −0.77 % | |
| `laser_linear` | amplitude `c − 1` | 5 % | −1.91 % | passed |
| | wavelength `k v_g − 1` | 2 % | +0.07 % | |
| | `E_r / E_x` at r = w / 2, `c_r / c_x − 1` | 10 % | +0.04 % | |

Reported: Gauss 5.2e-15 (`beam_linear`), 1.1e-13 (`laser_linear`), 2.5e-16 and 2.3e-15 (null
runs); rms residual of the fit 0.042 and 0.023 of the amplitude; amplitude at the output before
the last +0.42 % (`beam_linear`, t = 77.4) and −2.13 % (`laser_linear`, t = 52.9); the centre of
the pulse is measured at x = 53.98 against 54.04 of the group velocity of the grid, the bunch at
100.056 (measured and predicted). Wall time 66 s and 307 s (null runs 52 s and 264 s); peak
device memory 10.9 and 12.4 GB per GPU. The `beam_linear` figures repeat those of the lattice
test above (the uncommitted patch on 3823577) in all digits shown. One run each; the spread
between runs of one seed was not measured here.

Not shown by these runs: the amplitude far behind the bunch (the fall of 8 to 16 % in the
wavelengths 9 to 12 found in the lattice test stays unexplained; it lies outside the gate
window), and the same gates in a gate record — they enter `scripts/gpu_gates.py` with the stage
that closes the track, and the record is taken once at the final head.

### Stage 4 — moving window

W0, the compiled kernels (2026-11-12, the commit that adds the `#ifdef MOVING_WINDOW` macros
`XIDX`, `FLO`, `FHI`, `RX0`, `RX1` to `common.cuh` and uses them in `particles.cu`): the six cases
of `scripts/check_kernels.py` without a window, both units, both architectures (24 compilations)
give the same PTX as the commit 477223a, byte for byte (sha256 of the PTX of NVRTC 13.5; script
outside the repository, as fixed in the amendment). The two new cases with a window compile. The
other half of W0 (`mgpu_regress` bits) is taken with the GPU tests of the stage.

W0, the bits (2026-11-12): `tests/gpu/test_engine.py` with the stored bits of `mgpu_regress`
passes on the commit that puts the shift into the step (26 targeted GPU tests, no window in any
of the old ones). W0 holds in both halves.

W1 to W4 (2026-11-12, `tests/gpu/mgpu_window.py`, `tests/gpu/test_engine.py`; RTX PRO 6000):
**all pass.**

| Test | Criterion | Limit | Measured |
|---|---|---|---|
| W1 | density per slab of eight cells, 400 shifts | 2.5 % | 0.45 % (0.9955 to 1.0042; front slab 0.7497 against 0.75) |
| | mean u_x per slab | 0.02 v_te | 0.012 v_te |
| | ⟨u²⟩ / 3θ | 1.05 | 1.019 |
| | Gauss; div B | 1e-10; 1e-12 | 3.2e-14; 3.4e-16 |
| W2 | E against the box at rest, farther than 32 cells from x−, 200 shifts | 1e-6 max\|E\| | 6.9e-15 (1 GPU), 6.9e-15 (2,2,1), 7.1e-15 (1,2,2) |
| W3 | histories, 1 GPU against (2,1,1), (2,2,1), (1,2,2) | 1e-10 | kinetic 1.1e-14, electric 9.4e-16, magnetic 6.3e-15, `window_out` 5.4e-13, budget 1.2e-13 |
| | loaded planes | equal | 995,303 particles in every layout; test particles equal by a hash of (id, position, momentum) |
| W4 | energy behind the forward pulse that the box at rest does not have, over the pulse energy | 1e-4 | 5.3e-10 (1 GPU and (2,1,1)) |

W4 as measured. The plan fixed the limit, not the measurement. Set-up: a plane wave in vacuum
(a0 = 0.1, ω = 10, duration 0.3), periodic across, 128×8×8 cells, window from t = 1.5, 200 shifts;
the pulse that the antenna emits towards −x meets the lower face around t = 2.7. Measured is the
field energy in the 28 planes of the window behind the forward pulse, at every output from
t = 3.8 on, in units of the energy of the forward pulse (equal to the one of the pulse towards −x:
1.00009 passed through these planes). First measurement: 2.4e-11 at t = 3.8, but 1.5e-6 at t = 9
and 9.0e-5 at the end (t = 11.5), a rise that is not a reflection: the same run in a box at rest
has the same energy at the same laboratory planes (9.045e-5 at the end). It is the short waves of
the forward pulse (three periods long, so its spectrum reaches ω = 25, where the group velocity
of the grid is 0.82 c), which fall behind. The test therefore takes the difference to the box at
rest: at most 5.3e-10. A reflection would add to the window run only.

Two defects found by the tests and fixed in the same commits: a run split along x hung at the
first shift (`YeeGrid.shift_x` called the exchange of psi only on ranks with an absorbing layer: a
collective inside a rank-dependent branch); a run with a window and no species raised in the
check at its end.

OUT-D5 (2026-11-12): the window run of W2 on 2×2×1 GPUs and the box at rest on 1 GPU write the
same diagnostics at the same laboratory places (the regions are given in the window and in the
box). Read through `pypic.analysis` at 85 and 200 shifts: dump coordinates within 1.8e-15 and E_y
within 6.0e-15 of its maximum; slice 2.2e-15; density lineout 9.5e-15; x histogram 2.2e-16;
particles of a region: the same 245,716, positions within 5.3e-15. `scalars.h5` holds the shift
of every record (the last one, step 504: 196; the output at step 512: 200). The test fails
without the change of the diagnostics (checked once).

Not shown by Stage 4: a production-size window run, the cost of a shift, the window with the
mixed gather on a GPU (compiled only) — all three come with Stage 5 and the gates.

### Stage 5 — the gate configurations as window runs (2026-11-12, 40edad6, 8 GPUs, FP64): all four pass

| Run | Box | Amplitude (5 %) | Wavelength (2 %) | Radial (10 %) | Gauss (1e-10) |
|---|---|---|---|---|---|
| `laser_linear` | 1152×128², ranks (8,1,1), 2,024 shifts | −1.33 % | −0.02 % | +0.06 % | 6.6e-13 |
| `beam_linear` | 384×192², ranks (4,1,2), 482 shifts | +0.51 % | −0.25 % | −0.77 % | 9.1e-14 |

| Run without the driver | max \|E_x\| on the axis / predicted (1 %) | Density of the plateau (1 %) | Gauss |
|---|---|---|---|
| `laser_linear` | 0 | 1 exactly | 2.3e-15 |
| `beam_linear` | 0 | 1 exactly | 2.2e-16 |

Difference to the Stage 3 runs in the box at rest (68a3a29), last output, the same laboratory
interval (laser: x = 40.74 to 48.04, 372 samples; bunch: x = 63.80 to 95.06, 250 samples; the
outputs are float32): the E_x on the axis of the laser run differs by at most 1.9e-9 of its
maximum, the one of the bunch run in no stored bit; the three fitted numbers of each run agree to
the eight digits printed. The values of the laser run are the ones that the check of the short
interval gave before the run (−1.33 %, −0.02 %, +0.06 %); the bunch keeps the numbers of Stage 3.
The window changes nothing that the gates measure.

Cost: 235 M particles at the load in both; laser 6,160 steps in 122 s (17 to 20 ms per step),
bunch 1,416 steps in 27 s. The box at rest of Stage 3 holds 2.8 (laser) and 2.3 (bunch) times
the cells.

### Stage 5 — production set-up `laser_bubble`: formulas and analysis, before the first run (2026-11-12)

Formulas (`theory.blowout`), compared with the text of Lu et al., Phys. Rev. ST Accel. Beams 10,
061301 (2007) (arXiv physics/0612227): matched spot k_p R ≃ k_p w0 = 2 √a0 (Eq. 1); the field on
the axis e E_z / (m c ω_p) ≃ ξ / 2, peak √a0; etching speed c ω_p² / ω0²; phase velocity of the
wake c [1 − 3 ω_p² / (2 ω0²)]; dephasing length (2 / 3) (ω0² / ω_p²) R (Eq. 4); energy gain
(2 / 3) m c² (ω0 / ω_p)² a0 (Eq. 6). They agree with the numbers of the design above
(k_p R ≈ 4, slope 1 / 2); for a0 = 4, ω0 = 10 ω_pe: R = 4, peak E_x 2, dephasing length 267,
energy gain 267 m c². The run ends after 172 c/ω_pe of plasma, before the dephasing length.

The planner gets `regime: blowout` (laser only): no linear theory and no gate; the plane through
the axis as before, the electron histograms over (x, u_x) (two cells, Δu_x = 2) and over u_x
(Δu_x = 1) on −20 ≤ u_x ≤ 400 (a particle outside is not counted; the report says whether the
last bin is occupied).

Definitions of the report, fixed here:

- Pulse: centre = centroid of E_y² on the axis; a0 = max |E_y| / ω0 there.
- Cavity: the largest connected region of the plane with electron density below 0.5 n0 between
  1 c/ω_pe above the lower face (or the top of the ramp) and the centre of the pulse, after the
  density is averaged along x over one spacing across (6 samples, so that a sample is as long as
  wide). Radius = half of its largest extent across.
- Field: max and min of E_x on the axis in the same interval; the zero at which E_x rises inside
  the cavity that is nearest to its widest place, and the slope of a line fitted within
  1 c/ω_pe of it.
- Electrons: those of the bins at and above u_x = 2 ω0 / ω_pe = 20: number in n0 (c/ω_pe)³, mean
  and largest u_x (upper edge of the last occupied bin), the most populated bin; converted for a
  laser of 0.8 µm (n0 = 1.74e19 cm⁻³; 3.60e7 electrons, 5.76 pC per unit).
- The same numbers at every output (history).

CPU test on planted data (`test_wakefield_blowout_report_recovers_a_planted_cavity_field_and_bunch`):
a sphere without electrons behind a pulse, a smaller hole behind it, vacuum ahead, E_x with
several rising zeros, a histogram with a bin at the threshold.

### Stage 5 — probe of `laser_bubble` (2026-11-12, 61f454b, 8 GPUs)

`laser_bubble` with `x_end=12`, `n_out=2`: 1,842 steps to t = 18.08 (before the window starts),
2,099,609,600 particles at the load, scalars every 3 steps.

| | Probe | Design above |
|---|---|---|
| ms / step | 134–141 | 250 (80 min for 19,000 steps) |
| Peak device memory per GPU | 38.3 GB | 25 GB |
| Estimate for the 19,120 steps (up to 2.43e9 particles) | about 52 min | about 80 min |

The memory figure of the design was too low by a half; it stays far below the device (97 GB).
The estimate is below the limit of 3 h, so the production run is started.

What the output of a real run showed (analysis from end to end, no change to the definitions):

- The sum of the electron histogram is 37,052.647 n0 (c/ω_pe)³ at the three outputs (equal to
  3e-8) and is the sum of the profile over the loaded lattice planes, 4 ≤ x ≤ L_x − 1 (37,052.6).
- t = 18.1: a0 = 3.97 (4), centre of the pulse 12.019 (planned 12.010), cavity radius 3.81 while
  it is still on the ramp, slope of E_x 0.458, no electron above u_x = 20. Gauss 7.5e-12.
- Energy budget (monitored only): work of the antenna 1.004e5; the budget is off by at most
  7.6e2 while the pulse enters and by 63 (6e-4 of the work) at the end.
- Figure: the axis of E_x moved to the left (it lay on the colour bar of the map), with the line
  of slope 1 / 2 through the measured zero.

### Stage 5 — first production run of `laser_bubble` (2026-11-12, 61f454b, 8 GPUs, FP64)

Record: `bench/results/2026-11-12-61f454b3f-example-wakefield-laser_bubble-run-1c8c65d7.json`.
19,120 steps to t = 187.71 in 3,133 s (164 ms per step, 0.54 ns per particle and step), 8,117
shifts of the window; 2.100e9 particles at the load, 2.433e9 at the end.

| Last output (t = 187.7) | Run | Lu et al. 2007 |
|---|---|---|
| Radius of the cavity | 4.69 | 4 |
| Slope of E_x at its zero in the cavity | 0.446 | 0.5 |
| a0 of the pulse | 3.87 (largest of the outputs 5.56 at t = 84) | 4 at the focus |
| Centre of the pulse | 178.30 (planned 180.18; 0.980 c from t = 18.8 to the end) | |
| Electrons with u_x > 20 | 54.4 n0 (c/ω_pe)³; for 0.8 µm 1.96e9 electrons, 313 pC | |
| u_x of these electrons | mean 146, most populated bin 95.5, largest ≥ 400 (limit of the histogram) | energy gain 267 at the dephasing length |

The extremes of E_x on the axis behind the pulse are 4.42 and −8.66 at x = 174.25 and 174.31:
a spike 0.06 c/ω_pe wide, not the field of the cavity (which the slope measures; the peak of the
cavity without load is 2).

History (every output from the one at which the cavity exists):

| t | Centre of the pulse | a0 | Radius | Slope | Electrons with u_x > 20 | Largest u_x | Most populated bin |
|---|---|---|---|---|---|---|---|
| 18.8 | 12.7 | 3.98 | 3.94 | 0.453 | 0 | – | – |
| 28.2 | 22.0 | 4.03 | 3.94 | 0.427 | 0.000539 | 24 | 20.5 |
| 37.5 | 31.4 | 4.32 | 3.81 | 0.438 | 7.32 | 70 | 23.5 |
| 46.9 | 40.6 | 4.78 | 3.94 | 0.431 | 13 | 100 | 20.5 |
| 56.3 | 49.7 | 5.10 | 4.06 | 0.418 | 17.9 | 129 | 36.5 |
| 65.7 | 59.0 | 5.42 | 4.06 | 0.415 | 21.8 | 165 | 32.5 |
| 75.1 | 68.3 | 5.53 | 4.06 | 0.416 | 24.9 | 220 | 50.5 |
| 84.5 | 77.5 | 5.56 | 4.19 | 0.432 | 27.3 | 235 | 68.5 |
| 93.9 | 86.7 | 5.49 | 4.31 | 0.441 | 29.7 | 250 | 84.5 |
| 103.2 | 95.9 | 5.30 | 4.31 | 0.442 | 32.2 | 275 | 99.5 |
| 112.6 | 105.1 | 5.09 | 4.44 | 0.439 | 35 | 299 | 107.5 |
| 122.0 | 114.3 | 4.92 | 4.56 | 0.436 | 37.5 | 322 | 119.5 |
| 131.4 | 123.7 | 4.78 | 4.56 | 0.440 | 40.2 | 352 | 74.5 |
| 140.8 | 132.9 | 4.60 | 4.56 | 0.436 | 42.5 | 373 | 147.5 |
| 150.2 | 142.0 | 4.25 | 4.56 | 0.442 | 44.5 | 387 | 63.5 |
| 159.6 | 151.2 | 4.26 | 4.56 | 0.446 | 46.8 | 395 | 84.5 |
| 168.9 | 160.5 | 4.07 | 4.69 | 0.442 | 49.3 | 400 | 88.5 |
| 178.3 | 169.6 | 3.78 | 4.69 | 0.461 | 51.8 | 400 | 92.5 |
| 187.7 | 178.3 | 3.87 | 4.69 | 0.446 | 54.4 | 400 | 95.5 |

The pulse steepens to a0 = 5.6 at t = 84 and falls back to 3.9; the cavity widens from 3.9 to 4.7
with it; the slope stays between 0.41 and 0.46. Electrons are trapped from t = 28 on, all through
the run (the number still grows by 2.6 per output at the end). The last bin of the u_x histogram is
occupied from t = 168.9: the largest u_x of the last three outputs is not measured.

Diagnostics of the run:

- Gauss residual 2.7e-15 at step 0 and 1.69e-10 at step 14,516, the only two checks of the run:
  above the 1e-10 of the contract.
- div B 2.17e-12 at most (contract 1e-12).
- Energy (monitored only): work of the antenna 1.004e5; the budget is off by 162 at the end
  (1.6e-3 of the work) and by 743 at most; `window_out` 4,965.
- Device memory: 97,139 MiB at the peak on one GPU (the probe of 1,842 steps: 38.3 GB).

Four things of this list are not properties of the wake. Each has its own section below: the
device memory, the two Gauss checks, the Gauss residual and div B, the limit of the histogram.

### Stage 5 — device memory of a run whose ranks refill (e5c2117)

`Species.sort` gave its temporaries back to the device only after a sort of more particles than
any sort before. A rank of a window run loses particles at the lower face and refills: its n rises
from sort to sort but stays below its largest value, so the freed temporaries of each sort were too
small for the next one and stayed in CuPy's pool beside the new ones. The pool is now freed after a
sort of more particles than the last sort. `test_sort_of_a_growing_species_does_not_grow_the_memory_pool`
now starts from an empty pool, adds a fall followed by a rise below the largest n, and fails on the
old code.

A/B on a reduced `laser_bubble` (`n_perp=160 width=20 x_end=100`: 1536×160² cells, 5.25e8
particles, 10,900 steps, 8 GPUs), device memory of every GPU every 2 s:

| | ed83bb4 (before) | e5c2117 |
|---|---|---|
| Largest device memory of a GPU | 63,549 MiB | 12,859 MiB |
| Wall time, ms per step at the end | 448.5 s, 41.1 | 449.4 s, 41.2 |

### Stage 5 — the Gauss check of the planner (2e7911e)

The engine records the Gauss check only at a step of the history that is also a multiple of
`gauss_interval`. The planner set the two intervals independently (`n_steps // 500` and
`n_steps // 25`); in the production run (38 and 764 steps) only steps 0 and 14,516 are multiples
of both. The Gauss interval of the wakefield planner is now 20 history intervals.
`test_wakefield_plan_records_the_gauss_check_through_the_whole_run` fails on the old planner for
`laser_linear` (13 checks) and `laser_bubble` (2 checks).

Not changed here (outside this track, reported to the owner): the planners of `filamentation`,
`landau_damping`, `plasma_waves`, `plasma_expansion`, `relativistic_shock` and `reconnection` set
the two intervals independently in the same way. How many checks their runs keep was not counted.

### Stage 5 — Gauss residual at the planes between ranks: the index coordinate of the window (07c668c)

**Where the residual is.** A diagnostic script (not in the repository) ran the reduced
`laser_bubble` above on the code of the first production run and gave, every 440 steps, the
residual on every rank with the place of its largest value. Three places appear:

| t | Shift | Largest residual | Place |
|---|---|---|---|
| 8.6 to 21.6 | 0 | 6.9e-12 to 7.7e-12 | plane of the antenna (x = 1.98), no charge there |
| 25.9 | 0 | 2.3e-11 | on the axis, window cell 577: one node above the plane between ranks 2 and 3 (576) |
| 30.2 to 43.2 | 97 to 757 | 3.1e-11 to 3.2e-11 | on the axis at x = 15.10 (laboratory cell 769; 768 is the plane between ranks 3 and 4 before the first shift) |
| 47.5 to 56.2 | 977 to 1,417 | 5.4e-12 to 2.1e-11 | on the axis, x = 32.8, 38.4, 42.5 |
| 60.5 to 73.4 | 1,637 to 2,297 | 3.3e-11 to 3.5e-11 | on the axis, x = 45.3, then 47.5 at three checks |
| 77.7 to 86.4 | 2,517 to 2,957 | 2.4e-11 to 1.6e-11 | on the axis, x = 49.5, 54.2, 58.6 |
| 90.7, 95.0 | 3,177, 3,397 | 1.08e-10 | on the axis, laboratory cell 3,833, which was window cell 768 (the plane between ranks 3 and 4) at shift 3,065, between the checks at t = 86.4 and 90.7 |
| 99.4 to 107.0 | 3,617 to 4,007 | 1.30e-10 to 1.74e-10 | on the axis, window cell 769, in electrons of −550 to −595 n0 |

A residual stays at its laboratory place once it is there (nothing removes it) until that place
leaves the window, so a maximum is seen at the same place at several checks.

**It is reproducible at the planes between ranks.** Three runs of the same code and seed (A, B,
C; `x_end=40`), residual on the line of the axis at step 2,640 (before the first shift): within two
nodes of a plane between ranks the rms is 2.6e-12 in each and the three runs agree (correlation
+0.991 for each pair; node 576 +2.27e-11, +2.53e-11, +2.29e-11, node 577 −2.22e-11, −2.38e-11,
−2.33e-11); elsewhere the rms is 2.6e-13 to 2.8e-13 and the runs do not agree (correlation −0.01
to +0.06), as round-off of the deposit sums in a different order does (R10). A pair of opposite
signs on the two nodes of a plane between ranks is charge that moved by a part of a cell without
a current.

**Cause.** With a window the index coordinate was `(x − OX) / dx − s` with the origin OX of the
rank, and the compiler contracts the product and the subtraction of s into one fused operation.
`x − OX` is rounded with the origin of the rank, so two ranks give one position index coordinates
that differ by more than their whole cells; a particle that changes rank starts its next deposit at
another place than the one at which its last deposit ended. The two forms were compiled with the
options of the engine and run for 4,096 positions within one cell of the plane between two ranks,
with the cells of the production run (dx = 0.019635, origins of the ranks at cells 382 and 574):
without a shift the two ranks differ for 1,495 positions (by up to 2.8e-14 of a cell), at shift
917 for 1,582, at shift 8,117 (the end of the production run) for all 4,096 with a mean of
1.1e-12 of a cell: the difference grows with x, and late in the run it has one sign. With cells of
0.1 (the tests of Stage 4) no position differs: the inverse is the round number 10.

**Fix.** With a window the index coordinate is now `(x · (1 / dx) + XLOW) − (OXI + XLOW + s)`: the
product is rounded by itself (`__dmul_rn`), XLOW (the cells that the arrays of the first rank hold
below the box, the same number on every rank) keeps the sum from being negative, and the integer
taken away is not negative on any rank, so the subtraction is exact. Two ranks and two shifts now
differ by whole cells only, bit for bit. New define `XLOW` (`kernelsrc.particle_defines`,
`Decomposition.x_low`), exercised by the window cases of `check_kernels.py`. The kernels of a box
at rest are unchanged (W0 holds). A first form without XLOW was not exact on the first rank while
`OXI + s` is negative (found by an exact emulation before any run; not committed).

**Test.** The tracers case of `mgpu_window.py` now uses cells of 0.0999 and deposits 32 single
particles on every rank that holds the nodes around the middle plane;
`test_moving_window_is_decomposition_independent` asks that the ranks give every node the same
charge bit for bit. On the old form the largest difference is 1.3e-11 on (2,2,1); on the new form
it is 0 on (2,2,1) and (1,2,2), with the history hash of one GPU.

**Result on the reduced run.** Run D (`x_end=40`) and run E (`x_end=100`) used two earlier forms
of a rank-independent index coordinate (D: the fused form with the global cell; E: the product
rounded by itself, without XLOW), not the committed one; the second production run below is the
run of the committed form. At step 2,640 the nodes at the planes between ranks are like the others
(rms within two nodes 2.4e-13 and 2.6e-13, elsewhere 2.6e-13 and 2.7e-13; node 576 −6e-14 and
+1.3e-12; correlation with A, B, C at most 0.04). Largest residual of run E against the old code:

| t | Old code | Run E |
|---|---|---|
| 30.2 to 43.2 | 3.1e-11 to 3.2e-11 | 1.4e-11 to 1.6e-11 |
| 60.5 to 82.1 | 3.3e-11 to 3.5e-11, then 2.4e-11 and 2.0e-11 | 2.7e-11 to 2.8e-11 |
| 90.7 | 1.08e-10 | 1.7e-11 |
| 99.4 | 1.30e-10 | 2.4e-11 |
| 107.0 (end) | 1.74e-10 | 3.4e-11 |

**What remains, and what is not established.**

- From t = 30 on the largest residual of runs D and E is on the axis (y = z = 10), where the
  largest density of the run is 230 to 1,570 n0, at places that differ between the two runs
  (t = 34 to 43: D 1.9e-11 at x = 17.75, E 1.6e-11 at x = 17.51). In run E it rises from 1.7e-11
  to 3.4e-11 between t = 91 and 107. Its cause is not established. Whether it stays below 1e-10
  to t = 188 is read from the second production run.
- The residual in the plane of the antenna has the same size with the fix (7.4e-12 to 8.0e-12;
  old code 6.9e-12 to 7.7e-12).
- div B is unchanged by the fix: 1.1e-11 to 1.7e-11 in the plane of the antenna while the pulse
  enters (|B| up to 362 there), 2.1e-12 to 2.9e-12 in the pulse afterwards (|B| 42 to 56). With
  dx = 0.0196 one rounding unit of B = 50 over dx is 3.6e-13 (2.9e-12 for B = 362), so these
  values are six to eight units; the places differ from run to run. The 1e-12 of the contract is an absolute number and
  is not met by a field of this size on these cells. Reported to the owner; the contract is not
  changed here.
- The index coordinate of a box at rest, `(x − OX) / dx` (every axis, also y and z of a window
  run), has the same dependence on the origin of the rank. It is not changed in this track: it
  would change the bits of every run on `main`. Reported to the owner.

### Amendment (2026-11-12): range of the u_x histograms of `laser_bubble`

The set-up above fixed −20 ≤ u_x ≤ 400. The first run reached the limit at t = 168.9, so the
largest u_x of its last three outputs is not measured. `laser_bubble.yaml` now sets
`u_range: [-20, 520]` (da4e27e); the width of the bins is unchanged (1, and 2 over (x, u_x)), the
default of the planner stays 400. 520 is above the continuation of the fastest growth of the
largest u_x in the first run (2.6 per unit of time from t = 103 to 141: 495 at the end). No other
definition of the report changes.

### Stage 5 — second production run of `laser_bubble` (2026-11-12, da4e27e, 8 GPUs, FP64)

Record: `bench/results/2026-11-12-da4e27e2c-example-wakefield-laser_bubble-run-b738687c.json`.
The run of the committed code: the index coordinate of the window path (07c668c), the Gauss check
every 760 steps (2e7911e), the sort that frees its pool (e5c2117), histograms to u_x = 520.
19,120 steps to t = 187.71 in 3,098 s (162 ms per step, 0.53 ns per particle and step), 8,117
shifts; 2.100e9 particles at the load, 2.433e9 at the end.

**The questions left open above.**

- Gauss residual: 26 checks, 2.4e-15 at step 0, 3.72e-11 at most (step 18,240), 3.71e-11 at the
  last check. Below the 1e-10 of the contract through the whole run. It rises through the run
  (largest of each quarter of the checks: 2.0e-11, 2.9e-11, 3.6e-11, 3.7e-11); the cause of the
  residual on the axis is still not established.
- div B: 9.3e-12 to 1.6e-11 at the checks of steps 760 to 3,040 (the pulse enters through the
  plane of the antenna), 1.7e-12 to 2.7e-12 from step 3,800 on. Above the 1e-12 of the contract
  at every check after step 0, as in the first run; the section on the index coordinate gives the
  size of one rounding unit. Reported to the owner; the contract is not changed here.
- Device memory: 47,067 MiB at the peak on one GPU (every GPU every 20 s; first run 97,139 MiB).
- Histogram: the largest u_x is 457 at the end, below the limit of 520 in every output.

**The wake.** The last output gives the numbers of the first run: radius of the cavity 4.69
(Lu et al. 4), slope 0.446 (0.5), a0 3.87, centre of the pulse 178.30, extremes of E_x on the axis
4.42 and −8.66. The histories of a0, radius, slope and the minimum of E_x agree with the first run
in every digit of the analysis; the centre of the pulse differs by 1.8e-6 at most.

| t | Electrons with u_x > 20 | Largest u_x | first run |
|---|---|---|---|
| 150.2 | 44.5 | 387 | 387 |
| 159.6 | 46.8 | 405 | 395 |
| 168.9 | 49.3 | 424 | 400 (limit) |
| 178.3 | 51.8 | 441 | 400 (limit) |
| 187.7 | 54.4 | 457 | 400 (limit) |

Electrons with u_x > 20 at the end: 54.38 n0 (c/ω_pe)³, for 0.8 µm 1.955e9 electrons and 313 pC;
mean u_x 146.3, most populated bin 95.5. The largest u_x grows by 1.9 per unit of time from
t = 150.2 to the end. The energy gain of Lu et al. at the dephasing length is 267; the pulse has
crossed 172 of the 267 c/ω_pe of that length, and the largest u_x is above 267 from t = 103.2 on
(275).

Correction to the section of the first run: its histogram left out the electrons above 400, it
did not collect them in the last bin. Its number with u_x > 20 is equal to the one of this run to
3e-8 up to t = 150.2 and smaller by 6.1e-6 at t = 159.6 and by 0.0046 at the end; its largest u_x
at t = 159.6 (395) was already not the largest of the run (405), one output before the last bin
was occupied.

Energy (monitored only): work of the antenna 1.004e5; the budget is off by 162 at the end (1.6e-3
of the work) and by 743 at most; `window_out` 4,965.

### Gate records on d523e3b (2026-11-12, 8 GPUs): FP64 and `precision: mixed` ALL PASSED

One full gate run on the final code of the pull request, as fixed under Verification:
`python scripts/gpu_gates.py --gpus 8`, then `--precision mixed`, on the clean commit d523e3b.

| Record | Gates | Passed | GPU tests | Gauss max of all gates |
|---|---|---|---|---|
| `2026-11-12-d523e3bbb-gpu-gates.json` (FP64) | 25 | 25 | 110 passed, 0 failed, 0 skipped | 7.6e-12 |
| `2026-11-12-d523e3bbb-gpu-gates-mixed.json` | 23 | 23 | (not part of this mode) | 7.6e-12 |

The 25 are the GPU tests, five two-stream gates and 19 example gates; the mixed record has the
four explicit two-stream gates and the same 19. Two-stream growth against linear theory: +3.0,
−1.7, −2.3, +3.5 (implicit), −4.7 % (limit 10 %); implicit energy 1.9e-10 (limit 1e-8); these
figures equal those of the record of 17f1b5f on main in the digits given.

| Gate | Record | Amplitude (5 %) | Wavelength (2 %) | Radial (10 %) | Gauss (1e-10) |
|---|---|---|---|---|---|
| `wakefield:laser_linear` | FP64 | −1.33 % | −0.02 % | +0.06 % | 1.1e-13 |
| `wakefield:laser_linear` | mixed | −1.33 % | −0.02 % | +0.06 % | 1.1e-13 |
| `wakefield:beam_linear` | FP64 | +0.51 % | −0.25 % | −0.77 % | 3.9e-15 |
| `wakefield:beam_linear` | mixed | +0.51 % | −0.25 % | −0.77 % | 3.9e-15 |

The three measured ratios repeat the window runs of 40edad6 (Stage 5 above) in the digits given.
The Gauss residual does not: 6.6e-13 and 9.1e-14 there. Between the two commits the index
coordinate of the window (07c668c) and the steps of the Gauss check (2e7911e) changed; which of
them accounts for the difference was not examined.

With these records PIC-D20, PIC-D21, PIC-D22 and OUT-D5 are Accepted, and the close condition of
R8 (Context) is met. PIC-D23 (lattice places) stays Proposed: it entered by an amendment and
waits for the owner's decision; both wakefield gates use it. The window with the mixed gather
has no test of its own; the two wakefield gates of the mixed record are its GPU evidence.

### Review fixes (2026-11-12, after the gates on d523e3b)

Four findings of the review of PR #42; the first changes the engine, the second the analysis
package, the last two the example:

- `pypic/grid.py`, `pypic/simulation.py` (PIC-D22, OUT-D5): a shift of the window left J, the
  current of the step that the outputs read, one cell too far up, and zero after a shift that
  loads a plane with a temperature (the move of that plane wrote its current into J). `shift_x`
  now moves J one plane down with the fields, after an exchange of its x ghost planes, and the
  current of the move of a warm plane has arrays of its own, allocated at the first such plane.
- `pypic/analysis/animate.py`: `PhasePanel` shows the laboratory place of x on whichever axis of
  the histogram x is (only a histogram with x on its first axis moved with the window).
- `examples/wakefield/analyze.py`: a run with one output (`n_out = 1`) raised an IndexError: the
  output before the last is the one at t = 0, with no interval behind the driver to fit. Its
  amplitude is reported only where that interval holds a wavelength of the wake (nan otherwise);
  the verdict is that of the last frame.
- `examples/wakefield/analyze.py`: a blowout box with `drive: false` raised a KeyError (the report
  of a cavity read E_y, which such a run does not write). It is the null test now, against the
  largest E_x of the scaling, sqrt(a0); the movie of a run without the driver has no driver panel.

Checks: a test per finding, each failing on the code before the fix. GPU: the current before and
after a shift for a warm and for a cold plane (1 GPU), and J_y of the window dump on 2×2×1 GPUs
against the box at rest. CPU: the movie of a histogram over (u_x, x); `analyze()` on a planted
wake with two outputs; `analyze()` on `laser_bubble.yaml` with `drive: false`. Measured on
0e34252 (diagnostic): J_y of the two dumps within 4.4e-13 and 9.4e-13 of its maximum in the
window, E_y within 3.0e-15 and 5.0e-15 of the largest E_y of the box at rest.

Not measured: the cost of a shift with the one more exchange and the three more arrays that go
down. A run that loads warm planes holds three more grid arrays from the first such plane on; the
wakefield configurations are cold and do not allocate them.

### Amendment (2026-11-12): stacked on PR #43 and PR #41; `compact` refused with a window

Owner instruction (2026-11-12): the three pull requests must merge without conflicts. They
touched the same files (`particles.cu`, `kernelsrc.py`, `gpu_gates.py`, `PROGRESS.md`,
`bench/results/index.json`, the design files), so they are stacked by merge commits, nothing is
rebased: PR #41 into PR #43 (031851c), PR #43 into this branch (2263870). The merge order is
fixed: #41, #43, #42; the base of each stays `main`.

What the merge changed in this piece:

- The window kernels `window_drop` and `window_move` read and write the particles through the
  accessor layer of PIC-D19. In the float64 layout the index coordinate along x of that layer is
  `XIDX` and the range of the rank along x is `RX0` to `RX1`, so every kernel has the window.
- `precision: compact` with a moving window is refused (Context: the later of the two pull
  requests adds the refusal; that is this one): `window.check_config`, with a case of the CPU
  test of the refusals that fails without the check, and an `#error` in `particles.cu`.
- Lattice places (PIC-D23) with the compact storage: the loader sets the place through `WPOINT`;
  a case of `scripts/check_kernels.py` and a compact case of the lattice GPU test (2×2×1 GPUs).
- Verification, amended: the one full gate run of this pull request is taken on the head that
  holds both merges, in three precisions (FP64 with `pytest tests/gpu`, `mixed`, `compact`),
  because PR #43 adds the third. `scripts/gpu_gates.py --precision compact` leaves out the two
  wakefield gates, which run in a window; the window has no compact evidence and needs none.

The records on d523e3b are those of this piece before the review fixes and the merge. The
`pypic/` tree of the merged head is the tree of no record: its GPU tests and gates are open.

### Amendment (2026-11-12): one pull request for R3, R5, R6 and R8, gated once

Owner instruction (2026-11-12, after the review of 00c07ee and a078127): the three pull
requests become one, the review fix goes into it, and the GPU runs are cut down: what is on the
branch now and the fix are verified together. That replaces "three pull requests, each with
its own full gate run" of the owner's plan and the merge order of the amendment above:

- PR #42 (this branch, base `main`) is the one pull request. It holds PR #41 and PR #43 by the
  merges 2263870 and 969e73d. PR #41 and PR #43 are closed without a merge of their own.
- The gate run that had started on 969e73d was stopped in its first gate; it left no record.
- The review fix (the direction of the field in the filamentation planner) is described in
  `docs/plans/2026-11-11-more-examples.md` (last sections). It is in `examples/`.
- One full gate run on the head that holds the fix: `fp64` (with `pytest tests/gpu`), `mixed`
  and `compact`, 8 GPUs. No GPU run before it for the fix, which is a planner on the CPU.

### Gate records on 300c1a4 (2026-11-12, 8 GPUs): `fp64`, `mixed` and `compact` ALL PASSED

The one full gate run of the one pull request (amendment above), on the head that holds the three
pieces, both merges and the review fixes of the three: `scripts/gpu_gates.py --gpus 8` in `fp64`
(with `pytest tests/gpu`), then `--precision mixed`, then `--precision compact`, on the clean
commit 300c1a4, no other job on the GPUs. The commits after it change no file under `pypic/`,
`scripts/`, `tests/` or `examples/`.

| record (`bench/results/2026-11-12-300c1a498-…`) | precision | gates | result |
|---|---|---|---|
| `gpu-gates.json` | fp64 | 124 GPU tests; six two-stream gates (with `ei_implicit`); 24 example gates | ALL PASSED |
| `gpu-gates-mixed.json` | mixed | five two-stream gates; 24 example gates | ALL PASSED |
| `gpu-gates-compact.json` | compact | five two-stream gates; 22 example gates (without the two window gates) | ALL PASSED |

| Gate | Record | Amplitude (5 %) | Wavelength (2 %) | Radial (10 %) | Gauss (1e-10) |
|---|---|---|---|---|---|
| `wakefield:laser_linear` | fp64 | −1.33 % | −0.02 % | +0.06 % | 1.2e-13 |
| `wakefield:laser_linear` | mixed | −1.33 % | −0.02 % | +0.06 % | 1.2e-13 |
| `wakefield:beam_linear` | fp64 | +0.51 % | −0.25 % | −0.77 % | 3.5e-15 |
| `wakefield:beam_linear` | mixed | +0.51 % | −0.25 % | −0.77 % | 3.7e-15 |

The three measured ratios of both gates equal those of the records of d523e3b in the digits given,
in `fp64` and in `mixed`; the Gauss residual was 1.1e-13 / 3.9e-15 there (`fp64`).

| value | fp64 | mixed | compact | limit |
|---|---|---|---|---|
| growth error of `ei` / `ep` / `eip` | +2.95 / −1.68 / −2.26 % | +2.95 / −1.68 / −2.26 % | +2.95 / −1.68 / −2.26 % | 10 % |
| growth error of `buneman` / `ei_b0` | −4.66 / +2.23 % | −5.14 / +2.23 % | −4.66 / +2.23 % | 10 % |
| growth error and energy error of `ei_implicit` | +3.52 %, 1.9e-10 | not run | not run | 10 %, 1e-8 |
| largest Gauss residual of all gates | 7.6e-12 (`two-stream:buneman`) | 7.6e-12 (`two-stream:buneman`) | 2.6e-13 (`reconnection:pair_sigma10`) | 1e-10 |
| Gauss residual of the five explicit two-stream gates | 8.5e-13 to 7.6e-12 | 8.5e-13 to 7.6e-12 | 1.0e-14 to 4.4e-14 | 1e-10 |
| largest div B of the two-stream gates | 5.5e-13 (`ei_b0`) | 5.1e-13 (`ei_b0`) | 5.0e-13 (`ei_b0`) | 1e-12 |
| div B of `plasma_waves:perpendicular` | 1.07e-12 | 9.48e-13 | 9.64e-13 | reported (no criterion in this gate) |
| `bow_shock:gate`: magnetopause / shock, in R_mp | 0.925 / 1.382 | 0.916 / 1.380 | 0.923 / 1.381 | analysis of the example |
| `bow_shock:south`: magnetopause / shock, in R_mp | 0.867 / 1.116 | 0.872 / 1.147 | 0.878 / 1.161 | analysis of the example |
| `bow_shock:oblique`: magnetopause / shock, in R_mp | 0.960 / 1.249 | 0.949 / 1.221 | 0.962 / 1.233 | analysis of the example |
| largest Gauss residual of the three `bow_shock` gates | 7.0e-14 | 6.6e-14 | 6.3e-14 | 1e-10 |
| `relativistic_shock:pair_2d_sigma`: error of β_sh / n2/n1 / B2/B1 | −7.2 / +4.3 / +4.0 % | −7.2 / +4.3 / +4.1 % | −7.2 / +4.1 / +4.1 % | 10 % |
| `filamentation:pair_b0`: median growth error, m = 4…10 | +5.7 % | +5.7 % | +5.7 % | 10 % |
| `filamentation:pair_b0_stable`: mode-energy ratio | 6.02 | 6.02 | 6.02 | e² = 7.39 |
| `blast_wave:gate`: front speed in fast speeds | 2.102916 | 2.102916 | 2.102916 | analysis of the example |
| `blast_wave:gate`: `energy_budget_change` | 0.109002 | 0.108992 | 0.108973 | reported |
| `blast_wave:gate`: Gauss residual | 7.4e-15 | 6.2e-15 | 7.7e-15 | 1e-10 |

- Against the records of a213c5d (PR #43 with PR #41), these rows differ in the digits given: growth
  error of `buneman` / `ei_b0`: −4.66 / +2.23 % , −4.66 / +2.23 % , −5.14 / +2.23 % there; largest
  Gauss residual of all gates: 7.6e-12 (`two-stream:buneman`) , 7.7e-12 (`two-stream:buneman`) ,
  2.0e-13 (`filamentation:pair`) there; Gauss residual of the five explicit two-stream gates:
  8.4e-13 to 7.6e-12 , 8.4e-13 to 7.7e-12 , 1.2e-14 to 4.8e-14 there; largest div B of the
  two-stream gates: 5.0e-13 (`ei_b0`) , 4.8e-13 (`ei_b0`) , 6.2e-13 (`ei_b0`) there; div B of
  `plasma_waves:perpendicular`: 1.05e-12 , 1.01e-12 , 1.06e-12 there; `bow_shock:gate`: magnetopause
  / shock, in R_mp: 0.919 / 1.377 , 0.915 / 1.378 , 0.916 / 1.380 there; `bow_shock:south`:
  magnetopause / shock, in R_mp: 0.881 / 1.141 , 0.875 / 1.132 , 0.876 / 1.170 there;
  `bow_shock:oblique`: magnetopause / shock, in R_mp: 0.953 / 1.229 , 0.945 / 1.234 , 0.964 / 1.243
  there; largest Gauss residual of the three `bow_shock` gates: 6.9e-14 , 6.6e-14 , 6.9e-14 there;
  `relativistic_shock:pair_2d_sigma`: error of β_sh / n2/n1 / B2/B1: −7.2 / +4.7 / +4.1 % , −7.2 /
  +4.5 / +4.1 % , −7.2 / +4.3 / +4.1 % there; `blast_wave:gate`: `energy_budget_change`: 0.109002 ,
  0.108977 , 0.108972 there; `blast_wave:gate`: Gauss residual: 6.4e-15 , 6.4e-15 , 7.2e-15 there.
  The other rows are equal.
- The limit of div B is a criterion of the two-stream gates; the largest value there is that of
  `ei_b0`, the gate with B0. `plasma_waves:perpendicular` reports a div B of 1.07e-12 / 9.48e-13 /
  9.64e-13 (`fp64` / `mixed` / `compact`); the gate has no criterion on it. It is above 1e-12 in
  `fp64` here, and above 1e-12 in 22 of the 53 earlier full gate records that hold the gate
  (8.52e-13 to 1.17e-12). Not examined.
- `buneman` fits 0.051226 in `fp64`, 0.050969 in `mixed` and 0.051226 in `compact` here; a213c5d had
  0.051226, 0.051226 and 0.050969, and 2c8aa5f 0.051226, 0.050969 and 0.050969. Which of the two
  values a `mixed` or `compact` run gives was not examined.

What the records do not cover: the window with the compact storage (refused); the cost of a shift
with the exchange of J; a window run that loads warm planes at gate size (its GPU test is the
evidence). The 124 GPU tests run once, with the `fp64` gates.
