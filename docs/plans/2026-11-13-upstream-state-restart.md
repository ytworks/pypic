# Frozen-in initial states, a turbulent inflow face, and restart checkpoints

Status: In progress
Refs: roadmap R2 (restart), R12 (new: frozen-in fluctuation states and a turbulent inflow); PIC-D24, PIC-D25,
PIC-D26, OUT-D6 (new, Proposed); PIC-D9, PIC-D13, PIC-D18 (extended)

## Context

The research topic `shock-alfvenic-turbulence` (branch `research/shock-alfvenic-turbulence`, README "Approach")
needs a perpendicular shock (reflecting-wall piston) driven into large-amplitude Alfvénic fluctuations that are
frozen into the flowing plasma, and runs of 20–160 h. The owner's instructions (2026-11-13):
- engine changes needed by research go through a pull request from main;
- the PR is merged only after a review by Codex (model gpt-6.1-sol, reasoning effort ultra) reports no findings;
- the merge waits for the owner's instruction; until then the work does not stop.
The plan was approved with the research Approach on 2026-11-13.

What the engine cannot do today (source read on 15e5c2c):
1. **B with a normal component that varies along the wall.** A 1-D `profile` B_x(y) is divergence-free, but the
   reflecting wall zeroes normal B on its plane (`grid.apply_bc_B`), so cell 0 then has div B ≈ B_x/dx, which
   breaks the 1e-12 contract. A field that vanishes near the wall needs B_x(x, y) = ∂_y A_z, with B_y = −∂_x A_z.
2. **E tangential to a reflecting wall that is not uniform.** `initial_fields._check_tangential_E` refuses it;
   PIC-D18 allows a uniform `E0` only. The motional field of a fluctuating flow, E = −v × B, is such a field.
3. **Perturbations with a non-periodic factor.** `Beam.perturbation` takes sinusoidal modes along a periodic axis
   only, so δv cannot be confined to the upstream region.
4. **Restart** (R2): nothing exists. Long runs on GPUs shared with other jobs need it.

All three changes the research named go into this one pull request (owner, 2026-11-13: "出し惜しみするな"):
- the frozen-in initial state (E1, stage 1);
- the turbulent inflow face (E2, stage 3);
- restart (E3, stage 2), onto the same or a different number of GPUs.

The research's main runs use a unidirectional circularly polarized R-mode slab, which a loaded upstream can carry
(stage 1). The inflow face (stage 3) holds the age of the fluctuations at the shock fixed. It is needed for:
- the linear-polarization contrast;
- spectra with L-mode or counter-propagating content;
- runs longer than the decay times.

## Approach

All criteria below are fixed here, before the first run of each test. A criterion that fails is recorded as failed
and reported; a limit changes only by an amendment written into this file with its reason.

### Stage 1 — PIC-D24: frozen-in initial states

Three generic primitives of the initial state. They are not specific to an experiment.

1. **Vector potential** (`initial_fields` kind `potential`):
   - Form: `{"kind": "potential", "component": "Az", "shape": {"x": {...}, "y": {...}}, "scale": s}`.
   - A_c is s times the product of 1-D shapes (`pypic.profiles`) along the listed axes, evaluated on the c-edges (the
     positions of E_c). It is added to B as the discrete curl.
   - Example: for A_z, B_x += (A_z[j+1] − A_z[j]) / dy and B_y −= (A_z[i+1] − A_z[i]) / dx.
   - So div B stays zero bitwise. Several terms add up.
   - A is evaluated from global coordinates on every local index, ghosts included, so no exchange is needed and the
     result is independent of the decomposition.
   - Refused:
     - an A component that is not zero (|A| ≤ 1e-14 · max|A|) on the plane of a non-periodic face it is tangential
       to (the wall would zero the normal B and leave div B);
   - Amended 2026-11-13, before the first run: the first draft also refused a shape along an axis of one cell and
     the implicit solver with a non-periodic face. Both are dropped:
     - a shape along a one-cell axis is evaluated at the single place of that cell, the same on every rank;
     - the other initial-field kinds are allowed with the implicit solver at walls too.
2. **Shaped perturbations** (`Beam.perturbation` gains the key `shape`):
   - Form: `{"shape": {"y": {...}, "x": {...}}, "velocity": u1, "velocity_component": "x"}` or
     `{"shape": {...}, "displacement": ξ, "displacement_component": "y"}`.
   - The perturbation is the amplitude times the product of the 1-D shapes at the loaded position.
   - The axes may be periodic or not. The old form (axis, mode, phase) stays as it is.
   - Displacements use the existing virtual move of PIC-D9 (rounds below one cell, Esirkepov deposit, E −= Δt J),
     so div E = ρ holds to round-off.
   - Refused:
     - a displacement that is not zero for particles within one cell plus the largest displacement of a
       non-periodic face (the wall or an open face would take them);
     - `density` together with `shape` (the inverse of the cumulative density exists only for the periodic cosine).
3. **E from a charge displacement:** no new key.
   - A target E field is loaded by displacing one species by ξ = −E / (q n) with primitive 2: moving the charge
     deposits the current whose −Δt J is that E. This holds for any E, with div E = ρ exact.
   - The uniform part stays in `E0` (PIC-D18).
   - This is documented in PIC-D24 with the formula. The engine does not compute E = −v × B itself.

Configuration checks in `pypic/config.py` (CPU-testable). No kernel change and no new `#define`: the shapes are
evaluated in Python/CuPy, and the virtual move is the existing kernel `k_displace`.

### Stage 2 — PIC-D25 and OUT-D6: restart checkpoints (R2)

Order of work: stage 1, stage 2, stage 3. The research pilot needs stage 1 and stage 2 first; stage 3 continues
in the same pull request while the pilot is prepared.

- **Configuration:** `checkpoint: {every: N, keep: 2}` writes a checkpoint every N steps into
  `<output dir>/checkpoint/step-<n>/`. `run.py --restart <step dir>` or `restart:` in the configuration continues
  from it.
- **Files:** one directory per checkpoint.
  - Every rank writes its own arrays (`rank-<r>.npz`): fields with their ghosts, CPML ψ, the implicit arrays,
    species stores (all rows, every precision) with n and cap, the injection bank, `sdiag`, window and antenna state.
  - Rank 0 writes `meta.json`: format version, commit, the configuration, the layout, step, time, and the scalar
    state (history, kinetic_species, flux_out, dissipated, antenna work, output manager state).
  - The directory is completed by a rename after every rank has written (`step-<n>.tmp` → `step-<n>`), with a
    barrier, so a crash never leaves a partial checkpoint that looks complete.
- **Restore rules:**
  - The same grid and physical configuration are required; the number of ranks and the dims may differ (below).
  - The configuration must be equal except for an allow-list: the run length, output cadences, the checkpoint block.
  - The output directory is opened in append mode (OUT-D6 amends OUT-D1's refusal of an existing run directory):
    - histories continue from the stored step;
    - averaged diagnostics restore their accumulators;
    - tracks keep their ids.
- **Collectives:** every collective in write and restore is reached by all ranks. Only file writing is
  rank-local or rank 0.
- **Restart onto a different layout** (another number of GPUs or other dims):
  - Particles are written in global form: positions in global physical coordinates, compact rows expanded by
    their stored cell. On restore every rank reads all particle files in chunks and keeps the particles inside its
    domain, so each particle is kept exactly once (its cell decides).
  - Fields are written per rank with their global index ranges. On restore every rank assembles its owned block
    from the blocks that overlap it, then fills the ghosts by the usual exchange.
  - The injection bank is written per face column with its global column index, and restored by column.
  - Per-rank cumulative sums (`sdiag`) are reduced before writing and restored on rank 0, zero elsewhere. The
    reductions then give the same totals.
  - The same-layout restore is the fast path; the other one is checked to give the same result.

### Stage 3 — PIC-D26: a turbulent inflow face (E2)

A face that injects a flowing plasma (PIC-D13, `paired` or `flux`) can carry a prescribed fluctuation
δB(r⊥, t), δv(r⊥, t) of the inflow. The fluctuation is a sum of modes given in the frame of the inflowing plasma
(wave vector, polarization, amplitude, phase, and ω from the cold Hall-MHD dispersion or given). It is a function
of the coordinates along the face and of time only: the plasma carries the pattern in.
- **Particles:**
  - The beam table of the injected carriers gains a velocity δv(r⊥, t) evaluated at the birth point and step.
  - Paired partners and returned neutralisers take the same δv at their own birth point. Today they draw from the
    uniform table `rb`, which would leave a current n e (δv_i − δv_e) at the face.
  - The rejection bound `vcap` includes max |δv|.
  - Randomness stays keyed by (seed, step, global cell, slot), so injection remains decomposition-independent.
- **Fields:**
  - The fluctuation field enters through a total-field / scattered-field surface a few cells inside the face.
    Outside it the CPML absorbs only the scattered field. On it the incident fields (δB and E = −v × B of the
    incoming plasma) are added to the Yee update as equivalent surface currents.
  - The incident field is evaluated on the staggered positions from global coordinates, without exchanges.
  - The surface sits on unchecked nodes (the retreat of PIC-D15 and PIC-D16 is extended by its width), so no
    change of E reaches a Gauss-checked node.
  - The current that the injected δv carries must equal curl δB − ∂_t δE of the incident field at the surface.
    This is the consistency condition of the mode table: δv and δB are not independent.
  - The table is checked against the cold dispersion relation at set-up and refused beyond a tolerance.
- **Uniform background:** the uniform motional field E0 = −v0 × B0 with a reflecting face on the same axis as an
  open inflow face. Today it is refused, because the axis has an open face. It is extended: E0 tangential to the
  reflecting face is allowed when the open face of that axis injects the flow that carries it. In the TF/SF
  region the incident field holds it.
- **Energy record:** the incident Poynting flux through the TF/SF surface and the injected kinetic energy of δv
  enter the budget as source terms. The budget is monitored for the explicit solver.
- **Refused:** the implicit solver; fluctuations on a face that is not an inflow face; a window.

## Invariants

- **Gauss ≤ 1e-10 and div B ≤ 1e-12:** stage 1 keeps them by construction (discrete curl; charge moved by the
  deposit). Tests check both on 1 and 2 GPUs with reflecting walls.
- **1-GPU vs N-GPU ≤ 1e-10:** the shapes use global coordinates only; tested.
- **Restart:** a continued run must agree with an uninterrupted one. The deposit uses atomics, so two uninterrupted
  runs already differ in their last bits (R10). The criterion is relative ≤ 1e-10 on histories and fields after a
  short quiet run (no turbulent amplification), the same bound as 1 vs N GPUs.
- **Stages 1 and 2:** no change to kernels, the pusher, the deposit, the filter or the boundaries.
- **Stage 3:** it changes the injection kernel (`inject_face`, `inject_partners`, the neutraliser return in
  `push_deposit`) and adds the TF/SF field source.
  - A new `#define` goes into `pypic/kernelsrc.py` with a case in `scripts/check_kernels.py`.
  - Without a fluctuation table the kernels compute what they compute today, bit for bit. The existing injection
    tests and gates check that.

## Verification

- **CPU** (`tests/cpu`):
  - schema and refusals of `potential` and shaped perturbations;
  - the discrete curl of a product A gives div B = 0 on a NumPy grid;
  - `meta.json` round trip and the allow-list of the configuration comparison.
  - Each refusal test fails on the code before the change (checked once).
- **GPU** (`tests/gpu`):
  - (a) Potential A_z = a(y) h(x) with reflecting x walls and h = 0 near them: div B ≤ 1e-12, normal B zero on the
    wall, B equal to the analytic curl within 1e-12 relative.
  - (b) Shaped velocity and electron displacement that load E_z(x, y) = −(v × B)_z of an Alfvénic slab over a uniform
    drift in a reflecting-wall box (`tests/gpu/mgpu_frozen.py`): E_z equal to the target within 1e-3 of the
    perturbation; Gauss ≤ 1e-10 and div B ≤ 1e-12 over 240 steps; outside the light cones of the walls each species
    keeps its mean drift within 1 % (E × B balance of the uniform part). Amended 2026-11-13, before the first run of
    the test: the first wording asked for 2 Ω_ci⁻¹, which the light cones of the walls do not leave room for in a
    test box. Measured on the first run: E_z within 5.6e-5, drift within 0.2 %. The mode amplitude of δE_z changes
    by 8 % over 0.75 Ω_ci⁻¹ (a linear polarization is a sum of R and L waves of different speed at k d_i = 0.39).
    That is wave physics, not a test of the engine. The test state carries the Ampère current J_z = −∂_y δB_x on
    the ions; without it δE_z grew by 30 % (the probe of 2026-11-13).
  - (c) 1 vs 2 GPUs for (a)+(b): relative ≤ 1e-10.
  - (d) Restart: 2N steps against N + checkpoint + restart + N, for the two-stream benchmark configuration on 1 and
    2 GPUs, an open-flow configuration with paired injection (bank state), `precision: compact`, a window
    configuration, and the implicit solver. Histories and fields relative ≤ 1e-10; particle numbers and ids equal.
  - (e) A partial checkpoint (crash before the rename) is ignored; a changed configuration is refused.
  - (f) Restore onto another layout: a checkpoint from 2 GPUs restored on 1 and on 4 GPUs, against the same-layout
    restart. Relative ≤ 1e-10, particle numbers and ids equal.
  - (g) Turbulent inflow, open-flow box without a shock:
    - the fluctuation that enters and crosses the box equals the incident mode within 2 % in amplitude and 1 % in
      phase speed after one crossing;
    - Gauss ≤ 1e-10 on the checked region;
    - no current sheet at the face (|J| on the face layer within 5 % of the incident current);
    - the scattered field outside the TF/SF surface below 1 % of the incident amplitude;
    - 1 vs 2 GPUs ≤ 1e-10.
  - (h) Inflow into a reflecting-wall piston: the shock forms, and the upstream fluctuation at a fixed distance
    ahead of it keeps the amplitude of (g) within 5 % over the run.
- **Gates:** one full `scripts/gpu_gates.py --gpus 8` run on the final head of the pull request (FP64 and
  `mixed`; `compact` as in the last records), after the Codex review reports no findings. The gate set is
  extended by the turbulent-inflow check (g) as a gate of an example configuration if it fits the gate time,
otherwise as a GPU test; the research pilot is recorded in the topic, not here.
- **Review:** `~/.local/bin/codex review --base main -c model="gpt-6.1-sol" -c model_reasoning_effort="ultra"`,
  repeated after each round of fixes until it reports no findings. Each round and its outcome are recorded in the
  PR.
