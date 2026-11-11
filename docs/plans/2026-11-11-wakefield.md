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
