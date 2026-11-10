# Engine design

Status: **Accepted** (2026-11-09, PR #1; PIC-D8..D10 PR #2; PIC-D12, D13, D15, D16 2026-11-11, bow-shock plan);
PIC-D11 **Proposed** (2026-11-10, #10).
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
- **Rejected at construction:** the implicit solver; `B0` or `E0` together with `external_fields` (they are initial values of the Yee fields; the uniform external fields go into `external_fields`). A dipole needs a body that covers its core (particles inside the core see a field of 2m/a³); that rule comes with the body (PIC-D14).

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
- **History:** new keys `lost_charge`, `injected_kinetic`, `injected_charge`, `n_injected` (cumulative) and `fill` (the largest particles / capacity over ranks and species; by its definition it depends on the decomposition). `budget = total + lost_kinetic − injected_kinetic + boundary flux + dissipated`. A returned neutraliser adds its old energy minus its new energy to `lost_kinetic`.
- **Rejected at construction:** the implicit solver; an injecting face that is not open or that meets a reflecting face; test species; beams with a profile or a perturbation; a paired face without a neutraliser, without a carrier of the opposite charge, or out of which a species drifts; a neutraliser without a paired face.
- **Approximations:** at a paired face the neutraliser enters only as returned particles and partners, not with its own thermal flux. An edge or corner cell follows one face's mode, and a carrier entering from such a cell of a paired face always brings a partner (it pays no debt).

Why: the Phase 0 probe (plan, run table) showed that neither mode serves every face. Flux injection of a wind with u / v_te = 2.5 settled in a source sheath (density 1.51 at the face), where the paired face stayed stationary for 6e4 steps. Where the normal velocity is thermal the order was reversed: flux was clean and the paired rule gained 9 % density in 2e4 steps. Paying for a lost ion with the next entering ions left a debt of 611 macro-charges over 1024 columns, against 2233 when an electron was let out instead. TRISTAN injects pairs at one point; iPIC3D and MHD-AEPIC refill the boundary cells every step.

Null test M2 (a magnetised flow through one paired and four flux faces stays uniform): taken with the filter and the Marder zone, which open flow at the gate parameters needs (plan, Phase 0), and passed; the numbers and what it leaves open are under PIC-D16, Measured.

Measured (2026-11-11, RTX PRO 6000; diagnostic only):
- **Tests:** the number entering through a flux face is within 5σ of `n E[max(v_n, 0)] · area · time` (quadrature) for each species; a cold diagonal beam through two faces enters with the swept volume, fills the edge cells once and creates no field (electric energy < 1e-20); at a paired face no neutraliser is lost, the bank identity holds to 1e-9, no column is paid twice, and the partners have the mean energy of the inflow distribution within 5σ. Gauss's law holds to 1e-10 in all three.
- **Decomposition:** a magnetised oblique flow through six open faces (one paired, four flux), 1 GPU against a 2×2 rank grid: the same ids, the same injected counts and the same bank; histories within 1e-10. The same holds with `precision: mixed` (run by hand, not a test: histories within 2.7e-13).
- **Cost:** 120×160×160 cells, 4.9e7 particles (the per-rank size of the bow-shock gate), FP64, windows of 20 steps, against the commit before. Without injection: periodic 23.06–24.16 → 23.05–24.13 ms/step, six open faces 25.61–26.21 → 25.62–26.29. With one paired and four flux faces: 29.29–29.53 ms/step, of which the injection launches take 2.0 ms (synchronised timers, 7 %). The open run without injection loses 8 % of its particles in the measurement, so it is not a baseline at equal particle number.
- **Measured and replaced:** the paired mask as a run-time argument cost +11 % (periodic) and +5 % (open) without any injection; one thread per layer cell took 4.7 ms/step for the same launches.

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

Not established: the depth of the zone was not varied; whether the filter needs a zone at the absorbing sphere is open until PIC-D14.

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
