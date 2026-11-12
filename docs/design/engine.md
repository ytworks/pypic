# Engine design

Status: **Accepted** (2026-11-09, PR #1; PIC-D8..D10 PR #2; PIC-D12..D16 2026-11-11, bow-shock plan;
PIC-D17 2026-11-11, blast-wave plan); PIC-D11 **Proposed** (2026-11-10, #10); PIC-D19 **Proposed**
(2026-11-12, compact-storage plan).
Applies to: `pypic/` (engine), `pypic/setup/` (set-up layer), `examples/two_stream/`.

## Requirements (owner, 2026-11-09)

GPU-accelerated, multi-GPU PIC; relativistic Buneman–Boris pusher with numerically guaranteed
energy behavior; field scheme preserving the conservation laws; cubic macro-particles;
periodic / open / reflecting boundaries; everything 3-D; two-stream example; electron, ion
(mass ratio configurable) and positron species in e-i, e-p, e-i-p mixtures, charge neutral by
default; time step and plasma parameters chosen automatically for the experiment.

## Decisions

### PIC-D1 — Python + CuPy RawKernel, NVRTC-compiled CUDA C specialised by `#define`
Hot loops are CUDA C compiled at run time with the local grid as compile-time constants
(`pypic/kernelsrc.py`); orchestration, set-up and analysis are Python.
Why: no nvcc dependency, constant-folded indexing, one source for every grid layout; the same
source builder lets CPU CI compile-check every kernel with NVRTC (`scripts/check_kernels.py`).

### PIC-D2 — Explicit scheme: Yee + relativistic Buneman–Boris + Esirkepov, CIC shape
Leapfrog with B split in two half steps (fields stored at integer times); Esirkepov deposition
for CIC (cubic) particles; staggered CIC gather.
Why: Gauss law and div B to round-off; Boris keeps |u| exact in the magnetic rotation. Total
energy is not exactly conserved (finite-grid effects) — that is what PIC-D3 is for.

### PIC-D3 — Implicit scheme with exact energy and charge conservation (owner choice: both, switchable)
Crank–Nicolson Maxwell (preconditioned CG on `1 + dt^2/4 curl D curl`) + relativistic implicit
midpoint (`vbar = (u0+u1)/(g0+g1)`), orbit split at cell faces with the same path weights for
gather and deposit (Chen–Chacón type). Outer Picard on E; per-particle Picard → Newton (trust
region, line search) → 1-D bracketed solve; per-particle energy closure; B path-averaged across
walls. A step that does not converge is not committed.
Why: the owner asked for numerically guaranteed energy conservation; this is the scheme that
delivers it for the coupled particle–field system, at ~20× the explicit cost per step.

### PIC-D4 — Boundaries per face
Periodic; open = CPML (explicit) / matched conductive layer (implicit, dissipation counted
exactly), particles removed with their energy and charge recorded; reflecting = perfect
conductor + specular reflection with the deposit split at the wall.
Why: wall splitting keeps charge conservation; the implicit layer keeps the energy budget exact.

### PIC-D5 — 3-D Cartesian decomposition, one process per GPU, NCCL
Halo fill/reduce axis by axis (edges and corners by composition), particle migration axis by
axis, launcher `python -m pypic.launch -n N` (spawns ranks, distributes the NCCL id).
The explicit update of B computes one ghost layer from the exchanged ghosts of E instead of
exchanging B (issue #34); its products are rounded one by one, so ranks with and without an
absorbing layer give the same bits.
Why: no MPI dependency on the target machine; results equal to round-off for any rank layout.

### PIC-D6 — Neutral, decomposition-independent loading
All beams of all species are placed on the same positions (counter-based RNG keyed by seed,
global cell and slot); weights follow the densities, so rho = 0 and E = 0 exactly at t = 0.
Drifting Maxwell–Jüttner momenta: exact rejection sampler for theta < 1, Sobol for theta ≥ 1.
Why: neutrality without a Poisson solve; identical initial state on any number of GPUs.

### PIC-D7 — Experiment planner from linear theory
The two-stream planner solves the warm relativistic multi-beam dispersion relation for
k_max / gamma_max / omega_r and derives box, grid, dt (CFL, omega_pe dt, omega_ce dt, growth and
oscillation resolution) and run length, reporting each margin; every value is overridable.
Why: the owner asked that time step and plasma parameters be set optimally for the experiment.

### PIC-D8 — Persistent particle ids (2026-11-09, PR #2)
Particle storage gains an 8th row: an integer-valued float64 id, assigned by the loader as
`((beam - 1) · N_cells + global cell) · ppc + slot` (linear global cell index), so ids are unique and
independent of the rank layout; they move with the particle through sort, migration and removal.
Initialisation refuses runs with `n_beams · N_cells · ppc ≥ 2^53` (exact float64 integers).
Why: tracked orbits and particle samples for acceleration studies (OUT-D1/D2) must follow the
same physical macro-particles on any number of GPUs.

### PIC-D9 — Initial-state extensions (2026-11-09, PR #2)
- **Perturbations** (`Beam.perturbation`):
  - Velocity modes are added directly.
  - Density modes (exact inverse of the cumulative density, by bisection on a bracket that holds for every particle when |α| < 1; the residual is checked on all ranks before anything moves) and displacement modes move the loaded particles virtually. The move goes in rounds of less than one cell; each round's current is deposited with Esirkepov and the particles migrate between rounds (the number of rounds is agreed by all ranks). Then E −= Δt J, so div E = ρ holds to round-off from a neutral start, on any number of GPUs.
- **Explicit half step**: with non-zero initial E or B, a Boris step of −Δt/2 turns the loaded momenta into u^{−1/2}.
- **Test species** (`SpeciesConfig.test`, `particles`):
  - They use separate force and deposit charges (`qdep = 0`), also in the virtual move of their perturbations.
  - They are excluded from the budget and from Gauss's law.
  - They are loaded after the plasma species, so they never change the plasma's random numbers.
  - Explicit particle lists get ids after all beam ids.
- **Initial fields**:
  - Uniform `E0`. It is rejected if it is tangential to a conducting face.
  - Plane waves that are the exact discrete eigenmode of each solver (Yee with the stored-B half-step average, or Crank–Nicolson).
  - 1-D profiles (`pypic.profiles`) of a component along another axis, so the divergence is zero bitwise.
- **Density profiles** (`Beam.profile`): the profile sets the particle weights; shapes below 1e-10 are cut and those particles removed. Equal profiles keep the plasma exactly neutral, and neutrality is checked slot by slot.

Why: the analytic-solution examples (Landau and Dawson, orbits, Mora expansion, Harris sheets) need these initial states without giving up any conservation contract.

### PIC-D10 — Binomial current filter (2026-11-09, PR #2)
- `current_filter: {passes: N, compensate: bool}`: N passes of (1, 2, 1)/4 along every axis, plus an optional (−N/4, 1 + N/2, −N/4) compensator. It is applied after the ghost reduction and before the wall boundary conditions.
- **Wall rules**: node-type data are mirrored about the face node; wall-normal staggered data have a zero ghost; the hi face node is filtered too. With these rules the filter is symmetric and commutes with the discrete divergence.
- The Gauss check filters ρ in the same way.
- Implicit particles gather F·E_h, so the energy identity Σ E_h·(F J) = Σ (F E_h)·J stays exact.
- Along an axis that wraps onto the rank itself (periodic, not decomposed) the passes are applied as one composite stencil with indices modulo the period (issue #37): the same operator, still symmetric, equal to the passes one by one to rounding.
- Periodic and reflecting faces. Open faces were rejected until PIC-D15 (2026-11-11), which gives the rules there.

Why: the numerical Cherenkov instability of relativistic flows (shocks, reconnection outflows). Verified: transfer function to 1e-12; Gauss ≤ 1e-12 and implicit energy ≤ 1e-11 with walls; 1 vs 2 GPUs with a decomposed wall axis.

### PIC-D11 — Mixed-precision field gather (2026-11-10, #10) — Status: Proposed
- `precision: "fp64" | "mixed"` (default `fp64`). With `mixed`, the explicit push reads E and B from float32 copies made once per step. The six trilinear interpolations run in FP32.
- **Unchanged (FP64):** positions, momenta, the Boris rotation, the Esirkepov current deposit, ρ and the field solver.
- **Why the deposit stays FP64:** Gauss's law to 1e-10 needs the deposited current to satisfy the discrete continuity equation against the FP64 CIC charge. FP32 shape factors break that telescoping sum at ~1e-7 per step.
- **Rank independence:** the gather weights come from global index coordinates `x / dx`, computed in FP64 and rounded once to FP32. The local index is the global cell minus the rank's integer offset (`OXI..`). So the weights do not depend on the rank layout, and 1 vs N GPUs stays equal to round-off.
- **Scope:** explicit solver only. The implicit energy identity needs the gather and the deposit on the same FP64 path weights, so `mixed` with `implicit` is rejected.

Why: `push_deposit` is FP64-pipe bound on SM120, where FP64 runs at 1/64 of FP32. The gather is the part of the push that can leave FP64 without touching a conservation contract.

Measured (2026-11-10, 8× RTX PRO 6000; diagnostic only):
- **Speed:** the push without deposit falls from 18.4 to 14.2 ms per step. The 1-GPU two-stream step falls from 42.7 to 38.6 ms, the 8-GPU two-stream step from 9.4 to 9.1 ms, and the 8-GPU shock step from 14.4 to 13.4 ms.
- **Gates:** all 19 explicit gate configurations pass with `precision: mixed`: two-stream ei/ep/eip/buneman and the 15 analytic-example gates, including the particle-orbit Δt convergence. Growth errors and Gauss (≤ 7.7e-12) match FP64.
- **Tests:** 1 GPU and a 2×2 rank grid agree to 1e-10. A small plasma differs from FP64 by ~2e-11 in kinetic energy.

Status stays Proposed: the default remains `fp64` until the owner accepts `mixed` as a default.

### PIC-D12 — External static fields: the Yee arrays hold the difference field (2026-11-11) — Status: Accepted
Plan: `docs/plans/2026-11-10-bow-shock.md` (D-A, approved by the owner). Explicit solver only.
- `external_fields: {E: [..], B: [..], dipoles: [{moment, center, core}]}` (`ExternalFields`, `Dipole`). The particles feel the Yee fields plus the external fields; the Yee arrays hold the difference to them.
- **Unchanged:** the Maxwell update, the CPML, the wall conditions, the deposit. The external fields are static, so the difference fields obey the same equations. An external field may therefore be normal or tangential to any face, which `B0` / `E0` cannot (a conducting wall zeroes the normal B and the tangential E of the Yee fields). At an open face the CPML absorbs only the difference, which relaxes the fields there to the external values.
- **External E is uniform only** (div = curl = 0), e.g. the motional field −v×B of a flowing plasma. Gauss's law and div B hold for the Yee fields as before, with the same contracts.
- **Dipole:** B is the discrete curl, on the Yee faces, of `A = m × r / (r² + a²)^(3/2)` on the Yee edges (`pypic/external.py`). Its discrete divergence vanishes to round-off. For r ≫ a it is `(3 (m·n) n − m) / r³`; at the centre `2m / a³`. The current that sustains it (the discrete curl of this B, concentrated in the core) is a fixed external current: it is not deposited, so Ampère's law for the difference field is unchanged.
- **Rank independence:** the dipole is computed from global node indices, so a rank holds the same bits as the same faces of the undecomposed grid.
- **No periodic images:** along a periodic axis the dipole is not continuous across the box edge. It is meant for open axes; this is not rejected.
- **Gather:** the kernel `ext_sum` writes Yee + external array + uniform value once per step. FP64: only for the components that have an external part. `mixed`: it replaces the float32 copy of PIC-D11. Runs without external fields keep their bits.
- **Total field elsewhere:** the initial half kick, the field outputs (`Ex..Bz`) and the fields of the tracks use the total field (`Simulation.total_field`); `run.json` records it in `field_convention`.
- **Energies:** the history keys `electric`, `magnetic` and the Poynting outflow are those of the Yee (difference) fields. The work of the external E and of the dipole's current on the plasma is not recorded, so the monitored energy budget is not closed with external fields.
- **Rejected at construction:** the implicit solver; `B0` or `E0` together with `external_fields` (they are initial values of the Yee fields; the uniform external fields go into `external_fields`). A dipole without a body that covers its core (PIC-D14): inside the core the field is as strong as 2m/a³ and is not that of a dipole.

Why: a planet's dipole and the solar wind's field and motional E cross every face of an open box. As Yee fields they would be changed by the wall behind the absorbing layer and absorbed by the layer itself. OSIRIS, WarpX and Smilei also give the dipole as a field felt by the particles only.

Measured (2026-11-11, 1× RTX PRO 6000; diagnostic only):
- **Equivalence:** 50 test particles over 200 steps in external uniform E and B hold the same bits as in the same fields given by `B0` / `E0`, in FP64 and `mixed`.
- **Dipole:** discrete div B ≤ 8.9e-16 for max |B| = 1.7; within 0.7 % of the regularised formula beyond 6 core radii at `dx = 0.3..0.5` core radii; 1 GPU and a 2×2 rank grid agree to 1e-10 (`mixed`, with uniform fields).
- **Cost:** 120×160×160 cells, 4.9e7 particles (the per-rank size of the bow-shock gate), uniform E and B and one dipole against no external field, 3 alternating pairs of 20 steps: FP64 22.97–23.00 → 23.12–23.18 ms/step (+0.5 to +0.9 %), `mixed` 19.09–20.17 → 19.16–20.24 (+0.3 to +0.4 %).

### PIC-D13 — Plasma entering through open faces: `flux` and `paired` injection (2026-11-11) — Status: Accepted
Plan: `docs/plans/2026-11-10-bow-shock.md` (D-B and D-F, approved by the owner; the two modes come from its Phase 0 measurements). Closes roadmap R1. Explicit solver only.
- `injection: {faces: {"x-": "paired", "y-": "flux", ...}, species: null, neutralizer: "e"}` (`Injection`). `species: null` means every plasma species with beams. A species enters with the distribution of its own beams.
- **Layer:** every step, after the push of all species and before the ghost reduction of J, the one cell layer outside each injecting face is sampled like the loader: `ppc` slots per cell and beam, a uniform position, a drifting Maxwell–Jüttner momentum. Each sample moves one step in a straight line and enters if it ends inside the box. One cell is enough because a particle moves less than a cell per step. Kernel `inject_face`, one thread per (cell, slot, beam); the rules and the layers are in `pypic/injection.py`.
- **Keys:** the random numbers depend on (seed, step, global cell, slot) only, with a salt apart from the loader's, so what enters does not depend on the decomposition. Species and beams share the position key: electrons and ions are sampled at the same points.
- **Edges and corners:** a cell outside two or three faces is sampled once, and only if all of those faces inject. It belongs to the first of its faces in the order x, y, z and follows the mode of that face.
- **`flux`:** every injected species is launched on its own. The whole path from the layer into the box is deposited (Esirkepov), so Gauss's law holds at the face as it does for a particle that leaves.
- **`paired`** (the inflow face of a flow): one species is the neutraliser (the electrons), every other injected species is a carrier.
  - A neutraliser that reaches the face from inside is not removed. It is returned from the wall point with a momentum drawn from the inward flux distribution `v_n f` of its beams and finishes its step.
  - The neutraliser is not launched from the layer. Each entering carrier brings a neutraliser particle born at the same point, with the weight that cancels its charge. The pair appears at the end point of the carrier and deposits nothing. The partner's momentum is drawn from `v_n f` (kernel `inject_partners`, one launch per step).
  - **Bank:** a carrier that leaves through the face is removed, and its charge is added to the debt of the face column of its position at the start of the step. The next carriers entering in that column pay it: each enters without a partner and deposits its path. The payers are the first entering samples of the column's cell in the order (slot, beam), decided from the debt before the launch, so neither the thread order nor the decomposition changes the result.
  - The charge that has crossed a paired face is the unpaid rest: `injected_charge − lost_charge = −sign(q_carrier) · Σ bank`.
- **Compile-time mask:** `PAIRED_FACES` (bit `2·axis + side`, `pypic/kernelsrc.py`) tells `push_deposit` which faces return the neutraliser. With 0 the compiler removes that path, so a run without a paired face keeps its step time.
- **Ids** (PIC-D8): the injected particles take one block of ids per step behind the loaded ones, from (step, cell, stream, slot). A run stops on every rank before an id reaches 2^53.
- **Capacity:** `capacity_factor` (default 1.3) sizes the particle arrays against the loaded number. Before a step the arrays are grown if what the step can inject would not fit. Growing holds the old and the new array at once, so a large run sets the factor and watches `fill`. The planner's memory estimate counts 8 rows per particle times the factor.
  - **Sort** (amended 2026-11-11, Phase 5 of the plan): after a sort larger than any before, `Species.sort` returns its temporaries (40 bytes per particle of the sorted species) to the device; with a constant particle number the pool keeps them and the sort costs the same as before (an unconditional release cost 0.9 ms per sort of 3.1e6 particles). CuPy's pool keeps a freed block only for requests that fit, so with a growing particle number every sort allocated new blocks and kept the old ones: the bow-shock gate held up to 82 GB on one GPU for 10 GB of arrays. The bow-shock planner's estimate counts the sort, the ghost and absorbing layers, and the arrays an output step holds.
  - **Sort, amended 2026-11-12** (PR #43; plan `docs/plans/2026-11-11-compact-storage.md`, "R5 exit"): the temporaries of the sort are allocated from a memory pool of their own (`pypic.particles.sort_memory`). A sort of more particles than any before returns the blocks of that pool and the free blocks of the default pool to the device first; a sort of fewer particles is served from the blocks of the largest one, so the pool holds the temporaries of the largest sort and nothing else. The rule above freed only when a species passed its own largest number; under it the bytes held by the pool grew at the argsort of the sort while the bytes in use stayed flat (blast wave at gate size, every rank; at production size, not instrumented, the peak was 75.16 GB per GPU against the planner's 66.16), and the regression test reproduces the growth with a number that falls and rises again below its largest. Emptying the default pool after every sort was measured and not adopted (+7.8 % of the step in the two-stream benchmark); with the pool of the sort the benchmark shows no difference from the rule above. With it (ac90f90) the peak device memory of `blast_wave/full` is 61.78 GB per GPU in `fp64` (75.16 before; planner 66.16), and the peaks of the blast wave and the bow shock at gate size are below the planner's figure too. The pool of the sort holds 40.2 bytes per particle of the largest species on the rank (24.2 in `compact`; its number at the end of the run), and at the peak the default pool has at most 0.71 GB of free blocks. Tables: the plan, "Device memory with the pool of the sort".
- **History:** new keys `lost_charge`, `injected_kinetic`, `injected_charge`, `n_injected` (cumulative) and `fill` (the largest particles / capacity over ranks and species; by its definition it depends on the decomposition). `budget = total + lost_kinetic − injected_kinetic + boundary flux + dissipated` (plus `absorbed_kinetic` with a body, PIC-D14). A returned neutraliser adds its old energy minus its new energy to `lost_kinetic`.
- **Rejected at construction:** the implicit solver; an injecting face that is not open or that meets a reflecting face; test species; beams with a profile or a perturbation; a paired face without a neutraliser, without a carrier of the opposite charge, or out of which a species drifts; a neutraliser without a paired face; paired faces on more than one axis.
  - **Paired faces lie on one axis** (amended 2026-11-11, review of PR #39): a returned neutraliser takes a whole momentum of the inflow distribution, so near an edge where two paired faces meet it can cross the other face in the rest of its step, and again after that return. The wall loop of `push_deposit` takes three hits; a neutraliser still outside after them was removed like a particle at an open face, against the rule that a paired face returns it. Measured (1 GPU, 16×16×96 cells, dx 0.5, dt 0.2, 32 ppc, θ = 0.01, no drift, periodic z, 2000 steps): with x± and y± paired 713 neutralisers were lost (Gauss 1.3e-14, exact as before); with x± paired, none. The configuration is refused, not handled in the kernel: no example, test or gate has two paired axes, and the edge cells of two paired faces follow one face's mode (Approximations), which no null test has checked.
- **Approximations:** at a paired face the neutraliser enters only as returned particles and partners, not with its own thermal flux. An edge or corner cell follows one face's mode, and a carrier entering from such a cell of a paired face always brings a partner (it pays no debt).

Why: the Phase 0 probe (plan, run table) showed that neither mode serves every face. Flux injection of a wind with u / v_te = 2.5 settled in a source sheath (density 1.51 at the face), where the paired face stayed stationary for 6e4 steps. Where the normal velocity is thermal the order was reversed: flux was clean and the paired rule gained 9 % density in 2e4 steps. Paying for a lost ion with the next entering ions left a debt of 611 macro-charges over 1024 columns, against 2233 when an electron was let out instead. TRISTAN injects pairs at one point; iPIC3D and MHD-AEPIC refill the boundary cells every step.

Null test M2 (a magnetised flow through one paired and four flux faces stays uniform): taken with the filter and the Marder zone, which open flow at the gate parameters needs (plan, Phase 0), and passed; the numbers and what it leaves open are under PIC-D16, Measured.

Measured (2026-11-11, RTX PRO 6000; diagnostic only):
- **Tests:** the number entering through a flux face is within 5σ of `n E[max(v_n, 0)] · area · time` (quadrature) for each species; a cold diagonal beam through two faces enters with the swept volume, fills the edge cells once and creates no field (electric energy < 1e-20); at a paired face no neutraliser is lost, the bank identity holds to 1e-9, no column is paid twice, and the partners have the mean energy of the inflow distribution within 5σ. Gauss's law holds to 1e-10 in all three.
- **Decomposition:** a magnetised oblique flow through six open faces (one paired, four flux), 1 GPU against a 2×2 rank grid: the same ids, the same injected counts and the same bank; histories within 1e-10. The same holds with `precision: mixed` (run by hand, not a test: histories within 2.7e-13).
- **Cost:** 120×160×160 cells, 4.9e7 particles (the per-rank size of the bow-shock gate), FP64, windows of 20 steps, against the commit before. Without injection: periodic 23.06–24.16 → 23.05–24.13 ms/step, six open faces 25.61–26.21 → 25.62–26.29. With one paired and four flux faces: 29.29–29.53 ms/step, of which the injection launches take 2.0 ms (synchronised timers, 7 %). The open run without injection loses 8 % of its particles in the measurement, so it is not a baseline at equal particle number.
- **Measured and replaced:** the paired mask as a run-time argument cost +11 % (periodic) and +5 % (open) without any injection; one thread per layer cell took 4.7 ms/step for the same launches.

### PIC-D14 — Absorbing sphere: an insulator that keeps the charge of the particles it takes (2026-11-11) — Status: Accepted
Plan: `docs/plans/2026-11-10-bow-shock.md` (D-D; `body.relax` was built by the rule of its Phase 4 measurement). Explicit solver only.
- `body: {shape: sphere, center, radius, relax: 0.0}` (`Body`); the rules, the kernel constants and the masks are in `pypic/body.py`.
- **Absorption:** a particle that ends a step inside the sphere (`|x − center| < radius`) is removed after its whole path was deposited, by the removal path of the open faces (`push_deposit`, under `#ifdef BODY_SPHERE`). Test species are removed too.
- **Charge:** the charge of an absorbed particle stays at the nodes its shape covered at its end point. The body is an insulator: it charges up, and nothing neutralises its charge. The current filter spreads that charge `reach = passes + compensate` nodes further (PIC-D15).
- **Gauss check:** leaves out the nodes closer to the centre than `radius + (1 + reach) · |dx|`, `|dx| = sqrt(dx² + dy² + dz²)` (`body.gauss_radius`, `body.gauss_mask`). Over those nodes `Σ (div E − ρ) dV` is the absorbed charge of the history, to round-off. At every other node the contract is the one before.
- **Load:** a particle loaded inside the sphere gets weight zero (`mark_body`) and is removed by one collective migration after the set-up, without a second copy of the arrays. Species loaded at the same points stay neutral. These particles are not counted as absorbed.
- **`body.relax`** (rate in ω_pe; 0 = off): every step, after the update of E, the Yee E is multiplied by `exp(−relax · dt)` on every edge whose two nodes are inside the sphere (`body.relax_edges`).
  - Such an edge changes `div E` only at nodes inside the sphere, all outside the Gauss check, and takes from one of them what it gives to the other: the charge of the body moves towards its surface, its total is unchanged, and no checked node is touched.
  - The Yee field is the difference field (PIC-D12), so the inside becomes an equipotential of the difference field; the external E is not shorted.
  - The factor depends on the global position of the edge only, so ghosts keep the values of their owners: no exchange and no collective. The field energy it removes is not booked in `dissipated`.
  - Not for use with a dipole: with the gate dipole the insulating body passes the charge criteria of the measurement below and the body with `relax: 0.05` fails them (cases B and Br). Without a dipole it does not level the charge either (Ar). It is kept as an option that is off by default.
- **Box:** the ball of the left-out nodes must lie inside the box on every axis, and `reach + 1` cells clear of an open face (the Marder zone and the nodes the filter fills from the face, PIC-D15, PIC-D16).
- **Dipole:** a dipole of `external_fields` needs a body that covers its core, `|dipole.center − body.center| + core ≤ radius` (PIC-D12).
- **History:** new keys `absorbed_kinetic`, `absorbed_charge`, `n_absorbed` (cumulative; `sdiag` has 12 slots, the same reduction as before). `budget = total + lost_kinetic + absorbed_kinetic − injected_kinetic + boundary flux + dissipated`.
- **Compile-time constants:** `BODY_SPHERE`, `BODY_CX`, `BODY_CY`, `BODY_CZ`, `BODY_R2` (`pypic/kernelsrc.py`). Without a body none is defined and the kernels are compiled from the source as before: those runs keep their bits.
- **Rejected at construction:** the implicit solver; a negative `relax`; a shape other than `sphere`; a ball of left-out nodes that leaves the box or comes too close to an open face; a dipole whose core the body does not cover.

Why: an absorbing surface is what the kinetic models of small magnetospheres use for the planet (Lavorenti et al. 2022). A charge-conserving deposit leaves the charge of a removed particle on the grid (the same fact as at an open face, PIC-D15); keeping it and leaving the nodes out of the check costs nothing per step and needs no Poisson correction. With the dipole the charge of the body levels within three crossings (below), so the bow-shock example runs the insulating body.

Not established:
- Without a dipole the charge of the body does not level within three crossings, with or without `body.relax` (cases A and Ar below). The mechanism is not identified. No zone or neutralisation at the sphere was built: an absorber without a dipole is outside what was shown to reach a stationary state.
- A use of `body.relax`: in neither measured configuration did it bring the charge of the body to a stationary state. One rate (0.05 ω_pe), fixed before the runs and not tuned; other rates were not run.
- One size, one seed, three crossings.
- The step cost of the body was not measured on its own (the sphere runs below take 42.4 to 44.9 ms/step on 8 GPUs with their probes, against 40.7 for the null test M2 of the same size).

Measured (2026-11-11, RTX PRO 6000; diagnostic only):
- **Tests:** a drifting plasma on a sphere, without and with `CurrentFilter(1, True)`: no particle is left inside after the load and the load is neutral (1e-13); the number absorbed equals the number missing; the charge at the left-out nodes equals `absorbed_charge` (1e-12) and Gauss's law holds elsewhere (1e-12). One electron and one test particle flying in are each counted once, with the electron's energy and charge. With `relax = 30` the sum of E² on the relaxed edges is 5.6e-4 of the value without, with the same two identities. On the CPU: the left-out ball covers every node an absorbed particle can charge, and a field on the relaxed edges has no divergence outside the sphere and none in total.
- **Decomposition:** the open flow of PIC-D13 with a sphere centred on the edge shared by four ranks and `relax = 1`, 1 GPU against a 2×2 rank grid: the same ids, counts and bank; histories within 1e-10.
- **Sphere in the gate wind** (plan, Phase 4; criteria fixed before the runs; scratch scripts, not in the repository): the wind of the null test M2 (`aligned`) on a sphere of radius 13.33 at 50 from the inflow face, 240×320×320 cells, 3.9e8 particles, 2×2×2 ranks, `CurrentFilter(4, True)`, `marder_zone: 0.125`, 5761 steps (three crossings of 1920), FP64. A: no dipole. B: the dipole of the bow-shock gate (moment 5680 along −z, core 4, magnetopause planned at 20). Ar, Br: the same with `relax: 0.05`. Criteria: Gauss ≤ 1e-10; the charge Q of the body changes over the last crossing by at most 0.10 max|Q|; the rms charge per node of the body at the end is at most 1.10 of the value one crossing before; upstream density and u_x within 5 % of the wind over the last crossing.

  | Case | Gauss max | Change of Q, last crossing / max\|Q\| | Node rms, end / one crossing before | Upstream, worst | Verdict |
  |---|---|---|---|---|---|
  | A | 3.6e-14 | 0.299 | 1.374 | 0.08 % | failed |
  | B | 4.0e-14 | 0.043 | 1.056 | 0.11 % | passed |
  | Ar | 3.2e-14 | 0.386 | 1.534 | 0.08 % | failed |
  | Br | 4.1e-14 | 0.112 | 1.229 | 0.11 % | failed |
  - B: Q rises to +750 and stays between 751 and 759 from step 4320; 1479 / 621 / 427 particles are absorbed per step in the three crossings (4.9e6 in all). The surface on the stagnation line is at +145 to +191 T_e/e, the ion ram energy being 50 T_e. Electric energy of the difference field 0.97, magnetic 1.05 of the crossing before.
  - A: Q grows almost linearly to +851, by 64 per 480 steps at the end. The sphere takes 4500 particles per step (2.5e7 in all), with 0.4 % more ion than electron charge. The potential at the centre of the body (zero at the inflow face) is between −21 and −4 T_e/e at the six probed steps. The magnetic energy of the difference field still grows (1.24 of the crossing before).
  - Ar: Q reaches +903; the largest node charge is 6.7 against 2.0 in A (the relaxation gathers the charge at the surface).
  - Br: Q peaks at +521 (step 2401) and falls to +263, by 58 over the last crossing; the rms charge per node grows throughout (0.64 at the end against 0.18 in B; largest node 9.3 against 2.5). The surface on the stagnation line is at +21 to +75 T_e/e. 1409 / 632 / 387 particles are absorbed per step; the profile on the stagnation line is that of B (density 3.6 n_0 at 24 upstream of the centre).
  - 42.4 ms/step (A), 42.5 (Ar) and 44.9 (B and Br, 4.34e8 particles at the end, `fill` 0.56), probes every 480 steps included.

### PIC-D15 — Current filter with open faces: the deposit is cut at the face, which is the filter's wall (2026-11-11) — Status: Accepted
Plan: `docs/plans/2026-11-10-bow-shock.md` (D-E as amended by its Phase 0 measurements). Amends PIC-D10. Explicit solver only.
- **Cut of the deposit** (every open face, with or without a filter): the deposited path is clipped to the box. A particle that leaves deposits up to the face and leaves its charge on the nodes of the face, around the point where it crossed; a particle entering through a `flux` face (PIC-D13) deposits from the face on. Nothing is deposited in the absorbing layer. Before, the whole path was deposited and the charge stayed on the layer nodes around the end point, one cell behind the face. Kernel helper `deposit_cut` (`push_deposit`, `inject_face`).
- **Wall of the filter:** on an open axis the wall rules of PIC-D10 apply at the face of the box (nodes 0 and N), not at the conducting wall behind the layer, and the mirror rows in the layer are zeroed after the passes of the axis. The current in the layer is zero, and the filter stays symmetric and commutes with the divergence inside the box.
- **Gauss check:** the charge left on the face nodes belongs to no particle, so `div E − ρ` is not zero there, and the filter spreads it `reach = passes + compensate` nodes into the box. The check leaves out the face node (as before) and the next `reach` nodes at every open face, `reach + 1` with the Marder zone (PIC-D16). Without a filter and a zone the checked region is the one before.
- **Rejected at construction:** a filter with an open face in the implicit solver (its deposit is not cut).
- An open flow with a filter needs the Marder zone to stay stationary (plan, Phase 0: runs `prf4w`, `prf4mw`); the filter alone is accepted, the docstring of `CurrentFilter` says what it needs.
- **Bits:** the cut changes the results of runs with an open face that a charged particle crosses. Periodic and reflecting runs keep their bits.

Why: Phase 0 of the plan. With the filter's wall at the conducting wall the charge left in the layer is spread into the box with every pass; with a filter, neither the zone alone nor the face wall alone was stable, and the uncut deposit was worse (`prf4m`, `prf4w`, `prf4mn`). The cut puts the charge on one plane of nodes, where the zone takes it.

### PIC-D16 — Marder zone at open faces: the charge left at a face diffuses to the wall (2026-11-11) — Status: Accepted
Plan: `docs/plans/2026-11-10-bow-shock.md` (D-C; the Marder zone is the alternative the plan named, chosen by its Phase 0 measurements). Explicit solver only; needs an open face.
- `marder_zone: ν` (`SimConfig.marder_zone`), the Marder number, 0 < ν ≤ 1/6; 0 switches it off.
- **Zone:** at every open face, the nodes of the absorbing layer, the face node and the nodes up to `reach + 1` inside the box (`reach` of PIC-D15; 0 without a filter). Membership depends on the global node index only.
- **Step:** after the particles have moved and migrated, `F = div E − ρ` at the zone nodes, with ρ the charge density of the new positions, filtered like the current; F = 0 at the other nodes and on the conducting walls. Then `E_a += ν h² (F(+a) − F) / d_a` on every edge whose two nodes are in the zone, with `h² = 3 / Σ d_a⁻²`, followed by the wall condition and the ghost exchange of E. Kernels `marder_residual`, `marder_apply`; `Grid.marder_zone`.
- **What it does:** F diffuses inside the zone, `F ← F + ν h² ∇²F`. The conducting wall behind the layer holds F = 0 and takes it; no edge leads out of the zone into the box, so nothing of F moves inwards. ν ≤ 1/6 keeps the explicit diffusion stable.
- **Contract:** E changes only on edges between two zone nodes. Those nodes are outside the Gauss check (PIC-D15), so in the checked region E is advanced by the Maxwell update alone and Gauss's law holds by the deposit, as before. The zone is not a Poisson correction of the plasma: it removes the field of charge that no particle carries any more. `div B` is untouched.
- **ρ of the zone:** `deposit_rho` takes a shell width and skips the particles farther than `2 (reach + 1)` cells from every open face. That is exactly the set that reaches a zone node through the filter.
- **Collectives:** one ghost exchange of F and one of E per step, on every rank.
- **Rejected at construction:** ν outside [0, 1/6]; the implicit solver; no open face.

Why: the CPML is driven by the derivatives of the fields, so it does not relax the static field of the charge that removal and injection leave on the face nodes. In the Phase 0 probe, relaxing E in the layer changed nothing at any rate, and the zone did; ν = 0.125 and 0.03 gave the same state. TRISTAN uses slightly conducting boundaries for the same purpose.

Not established: the depth of the zone was not varied. There is no zone at the absorbing sphere: with the dipole its charge levels without one, without a dipole it does not level (PIC-D14, Not established).

Measured (2026-11-11, RTX PRO 6000; PIC-D15 and PIC-D16; diagnostic only):
- **Tests:** two beams leaving through open x and y faces with `CurrentFilter(1, True)`, with and without the zone: Gauss's law to 1e-12 in the checked region while charge is lost, and the ρ of the shell equals the full ρ at every zone node (1e-12). Without particles, the field of a charge on one face node: one zone step gives the value of the discrete diffusion at that node (1e-12), `div E` inside the box stays exactly zero, and after 1500 steps `max|div E|` is below 1e-6 of its first value.
- **Tight:** with one node less of retreat the Gauss check fails, and with one cell less of shell the ρ of the zone is wrong (both tried once by hand; neither is a test).
- **Decomposition:** the open flow of PIC-D13 with `CurrentFilter(1, True)` and the zone in a box of 8 cells along x (the zone ends at the rank face and the filter reaches across), 1 GPU against a 2×2 rank grid: the same ids, injected counts and bank; histories within 1e-10.
- **Cost:** 120×160×160 cells, 4.9e7 particles, six open faces with one paired and four flux faces, FP64, 1 GPU, windows of 20 steps: no filter 31.5–32.6 ms/step; zone alone 32.3–33.4 (+3 %); `CurrentFilter(4, True)` 34.4–35.8; filter and zone 40.0–41.2 (+13 % over the filter). 42 % of the particles are in the shell of this box (23 % of the bow-shock gate's box). With ρ deposited from every particle the zone cost +12 % alone and +17 % over the filter (ρ deposit 3.5 ms, its filter 2.2 ms, the zone kernels and exchanges 0.45 ms). The run without a filter took 29.3–29.5 ms/step in the measurement of PIC-D13; whether the cut of the deposit accounts for the difference was not separated.
- **Null test M2** (plan, Phase 3; criteria fixed before the run; scratch scripts, not in the repository): a uniform magnetised wind at the bow-shock gate parameters, 240×320×320 cells, 3.9e8 particles, 2×2×2 ranks, x− paired, y and z faces flux, x+ removing, external `E = −v×B`, `CurrentFilter(4, True)`, `marder_zone: 0.125`, 5760 steps (three crossings), FP64. Two cases: `aligned` (v along x, B along z) and `oblique` (v = (0.225, 0.0875, 0.065), B = (0.075, −0.1, 0.2165)). Both pass all five criteria.
  - Gauss max 3.0e-14 / 3.2e-14 (criterion 1e-10); `div B` 1.5e-15 / 3.2e-15.
  - Density of each species in the interior and in the slab of two 8-cell blocks at each face, over the last crossing, against the sample of step 0: at most 0.70 % / 1.25 % off (criterion 5 %). The largest values are at faces the wind leaves through or runs along (y+ of `aligned`; z+ of `oblique`).
  - Mean velocity, same samples: within 0.09 % / 0.84 % of |v| per component (criterion 5 %).
  - Electron temperature in the interior at the end: 1.003 / 1.024 T0 (criterion 1.5), level over the last crossing. In the face slabs at most 1.044 (y+, `aligned`, still rising by 0.001 per 480 steps) and 1.12 (z+, `oblique`; 1.14 at step 2400).
  - Energy of the difference fields, mean of the last crossing over the one before: electric 1.003 / 1.000, magnetic 1.002 / 1.043 (criterion 1.5).
  - 40.7 ms/step on 8 GPUs (0.83 GPU-ns per particle and step), 234 s for the run; `fill` 0.77 with `capacity_factor` 1.3. The `oblique` run shared the GPUs with the test suite, so its time is not recorded here.
  - Not established by M2: the magnetic energy of `oblique` had not levelled after three crossings (+0.7 % per 480 steps at the end; `aligned` is level from step 1440). Nothing was run beyond 5760 steps, at another seed or at another size, and the edges where the paired face meets a flux face were not looked at apart from the face slabs.

### PIC-D17 — Radial beams: a shape of the distance from a point, a drift along the radius, loaded into the cells the shape reaches (2026-11-11) — Status: Accepted
Plan: `docs/plans/2026-11-11-blast-wave.md` (Phase 1). Accepted with the gate record of `examples/blast_wave` (Phase 3; ae4bb69, FP64). Both solvers load this way; only the explicit one was run.
- `Beam.origin` (a point, c/ω_pe), `Beam.drift_radial` (a proper velocity along the unit vector from the origin; no drift for a particle at the origin itself), and the profile axis `"r"`: the shape is evaluated at the distance from the origin. Each of the two can be used without the other; a beam is radial when it has an origin.
- **Loader:** `load_particles` takes the radial drift and the origin as run-time arguments (no new `#define`), builds the drift vector from the global position of the particle and gives it to `boost_drift`. The random numbers stay keyed by (seed, global cell, slot) and a beam without a radial drift takes the path it took before.
- **Block:** a beam with a radial profile is loaded only into the global cells that the ball of its support reaches, `[floor((o − r_s) / dx), ceil((o + r_s) / dx))` per axis, cut to the box and then to the rank (`Simulation._beam_block`, `_local_block`; `Species.load(block=)`). `r_s = profiles.support(shape, r_max)` with `r_max` the distance from the origin to the farthest corner of the box: the shape is sampled at `r_max / 2^16` and `r_s` is one sample beyond the last one that reaches the cut 1e-10 (`profiles.CUT`, the weight below which a particle is not kept). A part of a shape narrower than a sample can be missed. Every other beam loads the box, as before. The id of a particle is the one a load of the whole box gives it.
- **Capacity:** a species has room for `capacity_factor · ppc ·` the cells its beams load on the rank; for a block smaller than the box a rank counts at least its even share of the block, `ceil(block cells / ranks)`. The share is there because `fill` is the largest particles / capacity over ranks and species (PIC-D13): a rank that loads nothing of the block would hold the particles that reach it in an array of the minimum size, grown to 1.25 of what it holds, and its `fill` near 1 would hide that of the species that fill the box.
- **Neutrality:** beams that load the same block share positions; the plasma beams of one block must cancel position by position (`max|Σ q w| ≤ 1e-12 max w`, one reduction over the ranks for all blocks). A block smaller than the box that does not cancel is an error on every rank: its charge has a field that the initial E = 0 does not hold. For the box the warning of before stays.
- **Removal of the particles with weight 0,** by species (the choice depends on the configuration only, every rank makes the same collective calls): a species with a profile is compacted (`compact_and_migrate`); of a species without one only the particles inside a body go, by the flag list and one migration (PIC-D14). Before, one profile in the configuration compacted every species and the flag path was taken only without any profile. Reason: the flag list holds `capacity / 8` entries and a species that overflows it is taken out as a whole (a copy), and a shell loaded into its block loses about half of its particles; compaction copies only what it keeps and keeps the order. The species that fill the box (the large ones) stay on the path without a copy when another species has a profile. A configuration with a body and a profile changes the path of its species without a profile (no example has both before `blast_wave`); their order changes, their particles do not.
- **Output:** the default momentum range of a species counts `|drift_radial|` like `|drift|`.
- **Rejected at construction:** a radial profile or `drift_radial` without an origin; an origin that neither uses; `drift_radial` together with `drift`; a perturbation on a radial beam; a radial beam in an injected species (PIC-D13: what enters through a face is the uniform beam).

Why: a blast wave is driven by a shell that expands from a point (Keenan et al., arXiv:2112.12114; Dieckmann et al. 2018), and a beam had one drift for all particles and a shape along one axis. A beam loaded the whole box with capacity for it: at 576³ a debris beam would take about 12 GB per GPU (plan, Context).

Not established: no run of the implicit solver with a radial beam; at the size of the gate the `fill` of the debris species and the time of the load were not measured (the gate gives the largest `fill` over ranks and species); the mean of `u_r` was tested for θ ≤ 0.01 only (the expectation used is first order in θ).

Measured (2026-11-11, RTX PRO 6000; diagnostic only):
- **Tests:** a Gaussian shell (centre 3, width 0.6, 2 n_0, θ_e = 0.01, radial drift 0.1) of electrons and ions in 32³ cells, fixed seed. The weight of each species is within 5σ of `n_0 ∫ f dV` (a bound of the sampling error); the mean of `u_r` is within 5σ of `u_d (1 + 5θ/2)`, the mean proper velocity of a drifting Maxwellian along its drift to first order in θ; the ids are those of a load of the whole box where the shape is not cut; Gauss's law to 1e-10 at the load and over 40 steps (2.3e-15). A shell of electrons with uniform ions is refused. In a run of the same set-up by hand: weight 245.81 against 245.36 (0.19σ), electron `u_r` 0.100184 against 0.100025 (1.7σ), 54,364 of the 110,592 loaded particles kept per species, `fill` 0.37.
- **Tight:** with the block one cell short on the high side the id comparison fails (528 particles of 54,364 missing; tried once by hand, not a test). The weight alone does not see it.
- **Decomposition:** an ambient plasma entering through six `flux` faces, an absorbing sphere and a drifting shell of two more species (two beams each) in 16×24×12 cells; with the 2×2 rank grid the block lies across the rank face along x, the sphere on it, and the two ranks of the high half of y load no shell particle. 1 GPU against 2×2: the same ids, the same counts (9748 shell particles loaded, 1264 of them in the high half after 40 steps, 44 absorbed, 5002 injected), histories within 1e-10, Gauss 1.9e-15.
- **Gate** (`blast_wave:gate`, ae4bb69, 8 GPUs, FP64; `bench/results/2026-11-11-ae4bb69c6-gpu-gates.json`): a debris shell of two species with four beams each (about 37e6 of the 943e6 particles) in 384³ cells with six `flux` faces and the absorbing sphere, 1760 steps: Gauss 7.9e-15, div B 9.0e-15, largest `fill` 0.67 with `capacity_factor` 1.5; all criteria of the example passed, and the other 22 gates of the record with it.
- **Gate with `precision: mixed`** (PIC-D11; 9bfa1a5, 8 GPUs; `bench/results/2026-11-11-9bfa1a537-gpu-gates-mixed.json`, after the FP64 record of the same commit): the same configuration passes all criteria: Gauss 6.9e-15, div B 8.6e-15, largest `fill` 0.67; the speed of the front and the debris extent have 7 digits or more in common with FP64.

### PIC-D19 — Compact particle storage: the cell and the position within it in 2^-32 of a cell; momenta stay float64 (2026-11-12) — Status: Proposed
Plan: `docs/plans/2026-11-11-compact-storage.md` (A0 to A2). The mode is written and compiles (`scripts/check_kernels.py`, two cases). Gates: the 21 explicit gates ALL PASSED with `precision: compact` on 2c8aa5f, 8 GPUs, beside the records in `fp64` and `mixed` of the same commit (2026-11-12, `bench/results/2026-11-12-2c8aa5f8e-gpu-gates-compact.json`). The decision stays Proposed until the owner decides on the memory criterion of R5 (under Measured) and on the default precision. What is below is the design; under "Measured" are the prototypes of A0, which are not the product kernels, the gyration runs that decided the momenta, and the speed, the device memory and the gates of the product kernels.
- `precision: "compact"` (default stays `fp64`): explicit solver only, refused with the implicit one like `mixed`. The push gathers from the float32 field copies of `mixed` (PIC-D11). The Yee arrays, the current deposit, ρ and the field solver stay FP64.
- **Storage, 56 bytes per slot in place of 64:** weight and id (two float64; the id stays the integer-valued double of PIC-D8), the cell of the particle in the arrays of its rank (int32, `IDX`), the position within that cell per axis in units of 2^-32 of a cell (uint32), the momenta (float64). `CompactSpecies` holds four arrays. A kernel reaches a particle only through the accessors at the head of `particles.cu`; `COMPACT_STORE` selects their compact bodies, and the builds of `fp64`, `mixed` and the implicit solver are unchanged by it.
- **A position is a lattice point.** The index coordinate, cell + fraction · 2^-32, is exact in a double while the cell index stays below 2^21. Every end point is put on the lattice before anything uses it; the Esirkepov deposit, the wall loop and `deposit_cut` are then computed in FP64 from the two exact end points, so the deposited current is the one of the stored positions (the charge conservation of PIC-D2 comes from the deposit, as before; no correction is added).
- **Moving kernels work in the global index coordinate** (cells from the lower corner of the box): the walls, the period and the range of the rank are integers there, so a reflection, a periodic wrap and a change of rank keep a lattice point exact, and the arithmetic of a path is the same on every rank layout. The local cell enters only in the stored cell and in the indices of the gather and the deposit.
- **Integer route of the push:** the end point of a step is one integer rounding per axis of (stored fraction + u dt / (γ dx)) · 2^32, so it depends on the fraction and on u only, not on the cell or the rank layout. A step whose end cell lies inside the box on every axis that is not periodic, and whose end point lies outside the body, meets no wall, no cut and no removal: its gather weights and Esirkepov shapes are built from the stored integers, the periodic wrap and the change of rank are integer tests, and cell and fraction are stored directly. Every other step takes the general path in the global index coordinate from the same end point. With the route on and off the particles had the same bits after 31 steps in three small cases (periodic; reflecting; open, reflecting and periodic with a body); the currents differed by 3.9e-16 at most.
- **Readers** (`Species.pos`, `mom`, `rows`; outputs, diagnostics, the body test) get FP64 physical values: the global index coordinate times the cell size, one rounding, independent of the rank layout.
- **Loader and injection:** the cell is the cell of the loop and the fraction the upper 32 bits of the random number that `fp64` takes for the position, written directly (not through the physical position, whose rounding could land in the neighbouring cell). The random numbers stay keyed by (seed, global cell, slot) (PIC-D6).
- **Migration:** the record is the 8 doubles of `fp64` — global index coordinate, the momenta, weight, id — and 3 more while the set-up displacement is applied, which is kept in an array of its own. The limits and the period are passed in cells; only `mig_take`, `mig_fill` and `mig_append` have a compact body.
- **Sort:** `argsort` of the cell array; the rows are gathered one at a time through one buffer.
- **Momenta:** float64, with the Boris rotation and the energy records of `fp64`; the rotation keeps `|u|` as it does there (the Boris gyration test runs in `compact` too). Float32 momenta (44 bytes per slot) were allowed by the owner (2026-11-11) under two conditions fixed before the code: all explicit gates pass in the mode, and in the gyration test the relative drift of `|u|` stays within 6e-8·√N without a secular trend; if one fails, the momenta stay float64 (56 bytes). The second is not met (2026-11-12, see Measured), so the momenta stay float64. FP32 Boris arithmetic is not used either: with the rotation factors rounded to float32 the rotation changes `|u|` by a fixed amount per step.
- **Refused at construction:** the implicit solver; a rank whose arrays (with ghosts) have 2^31 cells or more; a box with more than 2^20 cells along an axis including two ghost widths.

Why: production runs are bound by memory (bow shock `full`: 79.8 GB per GPU) and the push kernel by FP64 (71 % of the 8-GPU step). A position as cell and fraction takes 16 bytes in place of 24 and makes the walls and the wrap integer operations; the roadmap item R10 (a fixed-point deposit) can build on the same lattice.

Measured (2026-11-12, A0, 8× RTX PRO 6000; prototypes of the push, periodic only, one species pair; diagnostic only):
- **Speed, of the whole `mixed` step (3.38 ms on 8 GPUs, 50,331,648 particles in 256×64×64):** compact positions with float64 momenta −9.7 %; float32 momenta stored with FP64 Boris −7.9 to −8.0 %; float32 momenta with FP32 Boris −27.0 to −27.3 %. The sort of the compact layouts takes 1.9 to 2.3 % of the step less.
- **Gauss residual** over 400 steps on one rank: 4.4e-15 at most with the prototypes, 5.9e-13 with the engine.
- **Gyration with float32 stored momenta and FP64 Boris arithmetic** (2026-11-12; product kernel of f322fe2, one GPU, uniform B, 1,024 particles with random phases per value of `|u|`): the group mean of the relative change of `|u|` leaves zero as the steps go on. At 2,301,952 steps (B = 0.5, dt = 0.05) it is −19.6, +21.5, −15.9, +16.4 standard errors from zero for `|u|` = 1, 2, 5, 11.3 (a second seed: −20.0, +25.1, −16.0, +16.0), and for seven of eight other values of `|u|` 9.2 to 26.0; for `|u|` = 1 the mean per step is −7.1e-12, −6.8e-12, −6.5e-12 at 143,872, 575,488 and 2,301,952 steps. With float64 storage (`mixed`, and `compact` after the change) the same runs stay within 1.8 standard errors, 5.2e-20 per step at most. The mechanism is not established. Tables: the plan, "A2 — result".
- **Speed of the product kernels** (2026-11-12, b830153; the benchmark of A0, three runs per precision, medians): `kernel_push` plus `sort` against `mixed` −0.366 ms on 8 GPUs (−10.8 % of the 3.38 ms step) and −2.866 ms on 1 GPU (−15.8 % of 18.18 ms); the whole step 2.97 ms against 3.32 (`mixed`) and 3.65 (`fp64`) on 8 GPUs, 15.33 against 18.20 and 22.01 on 1 GPU. Before the integer route (5edc7f9) the same figure was +4.3 % (8 GPUs) and +9.7 % (1 GPU), above the 2 % of the A0 criterion. The benchmark is periodic without a body: every step takes the integer route. Tables: the plan, "A2 — speed of the product kernels".
- **Device memory at production size** (2026-11-12; `blast_wave/full`, 576³ cells, 3.13e9 particles, 8 GPUs, `pypic/` of b830153, one run per precision; `nvidia-smi` once per second, the GPU with the most): 62.48 GB in `compact` against 75.16 in `fp64` (−12.67 GB, 16.9 %), wall time 648 s against 788 s; both runs pass the criteria of the example. Against the planner (56.57 and 66.16 GB) the peaks are +10.5 % and +13.6 %: **the criterion of R5 (within 5 %) is not met, in either precision.** In instrumented runs of `blast_wave/gate` (bad2803) every rise of the bytes that CuPy's pool holds after the set-up, on all eight ranks and in both precisions, came from `cp.argsort` of the particle sort, while the bytes in use stayed at their early maximum (`fp64`: 17.11 GB in use, 22.06 held, planner 21.45; `compact`: 14.91, 18.37, 18.52). The production-size runs were not instrumented. Tables: the plan, "R5 exit".
- **Device memory with the pool of the sort** (2026-11-12; PIC-D13 as amended, `pypic/` of ac90f90, 8 GPUs, one run per precision; `nvidia-smi` once per second and `cudaMemGetInfo` at every allocation of the two pools, peaks agree to 0.001 GB): `blast_wave/full` 53.72 GB per GPU in `compact` against 61.78 in `fp64` (−8.06 GB, 13.0 %); `blast_wave/gate` 17.56 against 20.04; `bow_shock/gate` 9.01 against 10.28. No peak exceeds the planner's figure (production size 56.57 and 66.16 GB): they are 5.03 % and 6.63 % below it, 4.3 % to 6.6 % at gate size. **By its letter the criterion of R5 (within 5 %) is still not met; the difference is now on the low side.** The planner counts the sort of full arrays (`SORT_BYTES` per slot) and the runs fill 0.57 to 0.67 of their slots: with the sort term at the largest fill of the run the planner's formula is within 1.6 % of the peak in all six runs. The planner, `SORT_BYTES` and the limit are unchanged; an amendment of the criterion is proposed to the owner in the plan ("R5 memory criterion") and not applied. Measured particle arrays: 56.5 bytes per slot (64.5 in `fp64`) against the 56 and 64 of `SLOT_BYTES`; the half byte is `Species.flag_list`. Tables: the plan, "Device memory with the pool of the sort".
- **Speed with open faces, injection, a body, the filter and the Marder zone** (2026-11-12, b2d9e76, 8 GPUs; `bow_shock/gate` with 3.9e8 and `blast_wave/gate` with 9.4e8 particles, 100 steps after 20, three runs per precision, medians): the whole step is 34.10 ms against 38.10 (`mixed`) and 42.63 (`fp64`) in the bow shock, 72.59 against 82.44 and 90.99 in the blast wave: 10.5 % and 12.0 % below `mixed`, 20.0 % and 20.2 % below `fp64`. `kernel_push` is 16 % below that of `mixed` in both (18.32 against 21.79 ms, 42.89 against 51.03). `deposit_rho` of the Marder zone finds the particles near the open faces from the stored cell; with the test on the physical position (91de7e4) it took 2.54 and 5.40 ms per step against 1.59 and 3.12 in `mixed`, now 1.16 and 2.17. Tables: the plan, "A3 — speed in the open-boundary gate configurations".
- **Gates and final speed** (2026-11-12, 2c8aa5f, 8 GPUs): the growth errors of the two-stream gates are those of `mixed` (`ei` / `ep` / `eip` +2.95 / −1.68 / −2.26 % as in `fp64` too; `buneman` −5.14 %, `fp64` −4.66 %). The largest Gauss residual of the 21 gates is 1.9e-13 (`fp64` and `mixed`: 7.6e-12); it is 11 to 940 times below that of `fp64` in ten gates and 0.93 to 2.3 times in the other six with a residual. Two-stream benchmark at the same `pypic/` tree (medians of three): 2.84 ms per step on 8 GPUs against 3.17 (`mixed`) and 3.55 (`fp64`), 15.3 against 18.2 and 22.0 on 1 GPU; of the 1.94 ms of the push on 8 GPUs, 0.88 ms are the shape functions and the deposit (skipped with zero charge). Tables: the plan, "A3 — final breakdown" and "A3 — gates at the final head".

Not established: the speed in a build with reflecting walls at such a size; which arrays make the grid term of the blast-wave planner larger than the measured one (by 0.5 GB at gate size and 1.0 GB at production size, in both precisions; in the bow shock by 0.02 to 0.06 GB).

## Review amendments

- 2026-11-09 (PR #1 review): CG guards against 0/0 after exact convergence; exact MJ sampling;
  no commit of non-converged implicit steps; ghost width validated against cells per rank;
  nested planner overrides converted through `SimConfig.from_dict`; launcher polls all ranks.
- 2026-11-09 (PR #1 second review): implicit orbit folding chooses the wall from the unfolded
  midpoint (charge conservation for |d| > L/2); typed `SpeciesConfig` overrides accepted.
- 2026-11-09 (PR #2): implicit-step error flags are reduced over ranks before raising (a
  rank-local raise left the other ranks in a collective); a particle exactly on the upper
  reflecting wall of a decomposed axis stays with the edge rank (migration used to drop it).
  The in-memory mode / phase-space diagnostics and `Simulation.save()` were replaced by the
  output system (`docs/design/output-analysis.md`).
- 2026-11-09 (PR #2 review):
  - The density-perturbation inverse used Newton's method, which did not converge where the density nearly vanishes. With α = 0.99, 170 of 8192 particles were left, and the virtual move then asked for ~1e16 rounds. It is now a bisection with a residual check reduced over ranks.
  - The virtual move deposited a test species' real charge. It now uses `qdep`.
- 2026-11-10 (#6, PIC-D2 implementation): the explicit Esirkepov deposit iterates over the 3-node
  CIC support of each axis (from min(floor X0, floor X1)) and the 2 non-zero prefix sums of each
  current component, instead of 4 nodes and 3 prefix terms with zero tests; the push shares the
  node / half-cell gather weights of each axis. Same scheme and per-term arithmetic (round-off
  order only), fewer FP64 operations: the kernel is FP64-pipe bound on SM120 (FP64 = 1/64 FP32).
