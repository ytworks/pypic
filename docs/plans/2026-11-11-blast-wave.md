# Blast wave around a magnetised body (`examples/blast_wave`) and radial beams (PIC-D17)
Status: In progress — Phase 2 (the example). Phase 1 done: radial beams (PIC-D17, Proposed) gated on 574c904. Phase 0 done: N1 passed, N2 failed the charge criterion, which the owner then made a reported value (Phase 0, Amendment)
Refs: roadmap R6; PIC-D12..D16 (used as they are), PIC-D17 (new, Phase 1)

## Context

Owner request (2026-11-11): an example of a three-dimensional spherical shock with a magnetised
body at the centre; the way the shock is made follows prior work as closely as possible; the engine
gets the boundary conditions and obstacles the example needs.

Survey: no PIC or hybrid work was found that drives a flow onto a magnetised sphere from all
directions. The spherical shocks of the literature are driven outwards from the centre:

- Keenan et al. (arXiv:2112.12114; 3-D hybrid) and Le et al. 2021 (JGR Space Physics 126,
  e2021JA029125; 3-D hybrid): a debris cloud with a radial drift V_d in a uniform field. Parameters:
  the Alfvén Mach number M_A = V_d / v_A, the equal-mass radius R_m (the ambient mass inside R_m
  equals the debris mass), the debris gyroradius ρ_d = V_d / Ω_ci, and the magnetic stopping radius
  R_B = R_m M_A^(2/3) (debris energy = field energy inside R_B). The shock leaves the debris when
  ρ_d < R_m.
- Dieckmann et al. 2018 (Phys. Plasmas 25; 2-D full PIC): a dense shell expanding into a magnetised
  ambient plasma; the kind of shock depends on the angle to the field.
- Nersisyan & Osipyan (arXiv:0905.0051): an explosion in a dipole field.

Owner decisions (2026-11-11): the shock is driven by a cloud expanding from the centre; the field is
a uniform field plus a dipole at the centre.

Differences to the prior work, to be stated in the README: full PIC instead of hybrid; a reduced
mass ratio; a shell around the body instead of a Gaussian cloud; a smaller R_m / ρ_i; the dipole
and the body at the centre are additions of this example.

## Approach

Units: c = 1, lengths in d_e, times in 1/ω_pe, the ambient density n0 = 1. B0 along z.

Four species: ambient electrons and ions at rest, injected through the six open faces (`flux`,
PIC-D13); debris electrons and ions, loaded as a shell around the body with a radial drift V_d and
not injected. The uniform field B0 and the dipole (moment along B0, inside the absorbing sphere,
PIC-D14) are external fields (PIC-D12). Current filter and Marder zone as in the bow shock
(PIC-D15, PIC-D16). The debris species split their density over K = 4 beams with the same
positions and independent momenta, which lowers the weight ratio to the ambient particles from
16.5 to 4.1 without an engine change.

Common: dx = 0.5, dt = 0.25, 8 particles per cell and beam, ω_pe / Ω_ce = 4 (B0 = 0.25),
v_te = 0.1 (β = 0.64), T_i = T_e, M_A = 3, R_m = 8 d_i (ρ_d / R_m = 0.375, the value of Keenan et
al.), shell radius 0.33 R_m and width 0.1 R_m (Gaussian in r; peak density 16.5 n0), r_dip =
0.2 R_m (where the equatorial dipole field equals B0), body radius r_dip / 1.5, dipole core
0.2 r_dip, filter 4 + compensation, Marder number 0.125, 12 CPML cells, run 1.25 R_B / V_d.

| | gate | production |
|---|---|---|
| m_i / m_e, d_i | 16, 4 | 25, 5 |
| V_d = 3 v_A | 0.1875 | 0.15 |
| R_m, R_B | 32, 66.6 | 40, 83.2 |
| grid, half-width | 384³, 1.44 R_B | 576³, 1.73 R_B |
| r_dip, body radius | 6.4, 4.27 | 8, 5.33 |
| ambient particles | 9.1e8 | 3.1e9 |
| steps | 1760 | 2800 |
| estimate from the bow-shock records (not measured) | 3 min, 20 GB per GPU | 15 min, 58 GB per GPU |

What the engine lacks is the load of the debris: a density that depends on the radius and a drift
along the radius (PIC-D17, Phase 1). A beam is loaded into every cell of the box with capacity for
all of them; the debris beams would take about 12 GB per GPU each at 576³, so a radial beam is
loaded only into the block of cells that holds its shell.

Phases: 0 this plan, the planner skeleton and the null tests (no change of `pypic/`); 1 PIC-D17;
2 the example; 3 the gate; 4 the production run and the documents.

## Phase 0 — null tests (criteria fixed here before the first run)

Both runs: the gate configuration without the debris (`--set debris=false`), 8 GPUs, FP64.
N1 without body and dipole (`--set body=false`), N2 with them.

Samples: the 3-D output on blocks of 6 cells at 11 times (every 176 steps). Regions: the interior
(blocks four or more block layers from every face) and, for each of the six faces, the slab of the
two block layers at the face and the slab of the next two. N2 leaves out the blocks whose centre is
closer than 2 r_dip to the centre of the body or to the axis through it along B0 (the flux tube the
body empties). Moments of a region are sums over its blocks; references are the sample of step 0.
T = (2/3) (Σekin / Σn − m |Σflux / Σn|² / 2).

| Criterion | Limit |
|---|---|
| Gauss, every check | ≤ 1e-10 |
| div B, every check | ≤ 1e-12 |
| Density, each species, region and sample | within 2.5 % of step 0 |
| T_e, each region and sample | ≤ 1.05 of step 0 |
| Mean velocity \|Σflux / Σn\|, each species, region and sample | < 0.02 v_te |
| N2: potential of the body against the plasma at 2 r_dip, mean of the slice samples of the last third, in the plane z = centre and in the plane y = centre | ≤ 5 T_e / e each |
| N2: change of the body charge over the last third of the run | ≤ 0.10 max\|Q\| |

Reported without a limit: the lowest T_e ratio, T_i, the fill of the particle arrays, the cost per
step, the body charge.

Potential of the body: φ = ∫ ⟨E_r⟩ dr from the body radius to 2 r_dip, ⟨E_r⟩ the mean over a ring
of one cell in the slice (E_x cos + E_y sin of the samples of each component, at their own
positions). Change against the approved plan, which had |Q| / (4π R_b): the charge of the body is
shielded by a sheath of a few cells, so the vacuum formula overstates the potential by the factor
1 + R_b / (sheath thickness), 5 to 44 here; the bow-shock gate (Q = +744 on a radius of 13.3) would
read 445 T_e while ions of 50 T_e kept reaching the body. The limit of 5 T_e is unchanged.

Rules fixed with the criteria: if N2 fails, stop and report to the owner. If N1 fails, the fallback
is open faces without injection; the region of measurement of the gate then shrinks to
r < 1.17 R_B. A limit or definition changed after a run is recorded here as an amendment with its
reason.

### Results (2026-11-11; working tree on 7e89c55, no change of `pypic/`; 8 GPUs, FP64, one run each)

Both runs: 384³ cells, 9.06e8 particles, 1760 steps, rank grid (1, 2, 4), fill 0.667.

| Criterion | Limit | N1 | N2 |
|---|---|---|---|
| Gauss | ≤ 1e-10 | 8.3e-16 | 8.7e-16 |
| div B | ≤ 1e-12 | 4.5e-16 | 6.0e-16 |
| Density, worst region and sample (e / i) | 2.5 % | 0.44 / 0.42 % | 0.43 / 0.41 % |
| T_e, worst region and sample | ≤ 1.05 | 1.030 | 1.030 |
| Mean velocity (e / i), in v_te | < 0.02 | 0.0013 / 0.0006 | 0.0013 / 0.0006 |
| Potential of the body, planes z / y, in T_e / e | ≤ 5 | — | −0.33 / −0.87 |
| Change of the body charge over the last third | ≤ 0.10 max\|Q\| = 3.18 | — | **10.16 (failed)** |
| Verdict | | **passed** | **failed** |

Reported without a limit: T_e of the interior at the end 1.015 (both runs), lowest T_e ratio 0.997;
T_i at most 1.002; the worst density and T_e are in the slab at the y− face at the last sample;
total energy +2.2 % at the end of N1 (grid heating of the electrons, Debye length 0.2 dx);
85 ms/step (N1) and 89 ms/step (N2), 0.75–0.79 ns per particle and step × 8 GPUs; N2 measures
172,592 interior blocks (N1 175,616).

Body charge of N2 (particle weight 0.0156; 188,280 particles absorbed): −1.8 at step 0, −31.8 at
step 66 (the largest |Q|: the electrons of the flux tube arrive first), back to −16 by step 220,
then a steady rise: zero near step 880, +5.6 at step 1166, +15.75 at the end; +0.021 per step in
the middle third, +0.018 in the last. There the ions are absorbed at 52.3 particles per step, the
electrons at 51.2 (2.2 % fewer; 2.6 % in the middle third). The shot noise of the net count over
the last third (59,188 absorbed) is 1σ = 3.8 in charge: the limit (3.18) lies below it, but the
rise is not noise (the same sign and nearly the same slope from step 220 to the end).

Diagnostic run, not a re-test (N2 with `--set run=3.75`, 5360 steps): the charge does not level:
+68.5 at the end, +45.2 at two thirds; the slope falls only from +0.021 per step (steps 1340 to
2010) to +0.012 (the last 940 steps); both absorption rates fall (ions 54 → 43, electrons 54 → 42
particles per step from the first third to the last) with the ions 2 % ahead in the last two
thirds. The potential stays inside its limit (−0.67 / −1.14 T_e), as do Gauss (1.3e-15), div B
(1.4e-15), density (1.1 %) and mean velocity (0.0025 v_te). T_e exceeds its limit at this length
(1.060 in the slab at the y+ face, 1.040 in the interior: about +1.5 % per 1760 steps), which the
criteria, fixed for 1760 steps, do not cover.

What the runs show: outside the excluded region the plasma does not see the body (the worst
density and T_e of N2 differ from those of N1 by 0.01 % of n0 and 0.02 % of T0), and the potential
of the body against the plasma at 2 r_dip stays within 1.2 T_e / e and does not follow Q (plane y:
−0.82, −0.78, −0.87 over the thirds of N2 while Q goes from −32 to +16, and −1.14 in the last third
of the long run at Q = +45 to +68): the charge is shielded close to the surface (in the last third
of N2 the ring means in the first d_e from the surface are n_e / n_i = 0.52 / 0.46 n0 in the plane
z and 0.59 / 0.53 in the plane y; both recover to n0 by about 14 d_e). Not shown: why the ions are
ahead (the likely reading, not tested: the electrons reach the body only along the flux tube and
are mirrored by the dipole, the ions also across the field, thermal gyroradius 1.6 d_e against a
body radius of 4.27); why the potential stays negative under a positive charge; a level of the
charge; the charge and the potential under the debris, where the cavity removes the shielding
plasma (an unshielded Q = +16 on the body is Q / (4π R_b) = 0.29 = 29 T_e / e; the kinetic energy
of a debris ion is 0.28).

How the criterion came about: the approved plan asked for a potential within 5 T_e whose change
over the last third is below 10 %, with the potential taken as |Q| / (4π R_b); when the potential
was redefined before the runs (above), the second half stayed a criterion on Q.

### Amendment (2026-11-11, owner decision after the results)

The criterion "change of the body charge over the last third ≤ 0.10 max|Q|" is no longer a limit of
N2; the charge, its largest value and its change over the last third are reported. N2 stays recorded
as failed against the criteria fixed before the run; the six criteria that remain were met in it.
The work goes on to Phase 1.

Reason: what the criterion was to protect is measured directly and holds at one and at three times
the run length (the plasma outside 2 r_dip as in N1, the potential of the body within 1.2 T_e / e
and not following Q); the limit was a remainder of the vacuum formula for the potential and lies
below the shot noise of the absorbed count. Not decided by this: a level of the charge was not
shown, and the charge and the potential of the body under the debris are not known; the gate
reports the potential of the body without a limit, as in the approved plan.

## Phase 1 — radial beams (PIC-D17)

- `Beam.origin`, `Beam.drift_radial` (a constant proper velocity along the unit vector from the
  origin), profile `axis: "r"` (the distance from the origin). Rejected at construction: a radial
  profile or drift without an origin, `drift_radial` together with `drift`, a perturbation on a
  radial beam, a radial beam in an injected species.
- `profiles.support(spec, rmax)`: the largest radius at which the shape is 1e-10 or more.
- `load_particles` builds the unit vector from the global position and passes it to `boost_drift`;
  the arguments are run-time values (no new `#define`); the random numbers stay keyed by
  (seed, global cell, slot).
- A radial beam is loaded into the block of cells that holds its support, cut to the rank; the
  capacity of a species counts the cells its beams load. The species that share a block must load
  it neutral: checked per block, an error otherwise. The zero-weight particles are removed by one
  path.
- `default_umax` of the output includes `drift_radial`.
- Tests. GPU: mass and radial momentum of a shell against its shape (fixed seed, sigma bounds);
  Gauss ≤ 1e-10 at step 0; a block that straddles rank faces; 1 against 2×2×1 GPUs, history within
  1e-10 (`tests/gpu/mgpu_radial.py`). CPU: the rejections, `support`.
- Acceptance: `pytest tests/gpu`, `scripts/check_kernels.py`, the gates of the existing examples
  unchanged.

### Changes while building (2026-11-11)

- Removal of the zero-weight particles: not one path. A species with a profile is compacted, a
  species without one loses only the particles inside the body, by the flag list. Fact: the flag
  list holds capacity / 8 entries and a species that overflows it is copied as a whole; a shell
  loaded into its block loses about half of its particles (54,364 of 110,592 in the test set-up
  were kept), so the flag path would copy every shell species, and compaction of every species (the
  path of a configuration with a profile before) would copy the ambient species, 9.1e8 particles at
  the gate.
- Capacity: for a block smaller than the box a rank counts at least its even share of the block.
  Fact: `fill` is the largest particles / capacity over ranks and species, and an array grown on
  demand is grown to 1.25 of what it must hold; without the share the debris arrays of the ranks that
  load no block cell would set `fill`, and the gate could not show the fill of the ambient species,
  which decides the production size (Risks: capacity).
- `Beam.origin` without a radial profile or drift is refused too (nothing would use it).

### Result (2026-11-11)

Commit 574c904, 8 GPUs, FP64 (`bench/results/2026-11-11-574c90452-gpu-gates.json`): 93 GPU tests
pass (the new ones: mass and radial momentum of a shell within 5σ, ids, Gauss 2e-15 at step 0, 1
against 2×2 GPUs with a block across a rank face); `scripts/check_kernels.py` passes; all 21 gates
pass. The growth rates of the five two-stream gates equal those of 9b901c1 bit for bit;
`bow_shock:gate` gives magnetopause 0.922 and shock 1.381 R_mp (9b901c1: 0.917, 1.381; the recorded
runs of one seed spread over 0.917–0.922). No gated example uses a radial beam: the gates show the
existing paths unchanged, the radial path is covered by the tests above and by the gate of Phase 3.
`--precision mixed` was not run on this commit (Phase 3 runs both).

## Phase 2 — the example

`examples/blast_wave/{plan.py, theory.py, analyze.py, README.md, configs/gate.yaml, configs/full.yaml}`.

- Planner rejections: Ω_ce dt > 2 at the surface of the body; ω_pe dt > 1.2 at the peak of the
  shell; particle ids; capacity; a shell that overlaps the body; a body radius below 6 cells; a
  half-width below 1.3 R_B.
- `theory.py`: peak density of the shell, R_B, the debris gyroradius. The cavity radius of Le et al.
  2021 is used only after its formula is copied from the paper with its equation number.
- `alfven_speed`, `plasma_beta`, `fast_mach`, `compression` move to `pypic/setup/mhd.py` (second
  user); the gate numbers of the bow shock must not change.
- CPU tests: the engine accepts the planner's configuration; each rejection; the cavity radius of
  synthetic data is recovered; nan fails the gate.

## Phase 3 — gate (criteria fixed here before the first run of the gate configuration)

Measured in the plane z = centre (perpendicular to B0), B_z of the total field.
Front r_f(t, direction): the outermost radius on a ray at which B_z / B0, averaged over one d_i
along the ray, exceeds 1.25. Rays: ±x, ±y (axes) and the four diagonals. Window: the slice samples
with the mean r_f between 0.3 and 0.7 R_B. Cavity r_c(t): the outermost radius below the front at
which the ring mean of B_z / B0 is below 0.5.

| Criterion | Limit | |
|---|---|---|
| Gauss, div B, every check | ≤ 1e-10, ≤ 1e-12 | gate |
| Ambient plasma ahead of the front (3-D blocks with r > largest r_f + 0.2 R_B, outside the flux tube of the body), every sample: density of each ambient species and \|B\| | within 5 % of step 0; no such block at the last sample = fail | gate |
| Largest cavity radius over the run | 0.5 to 1.1 R_B (replaced by the formula of Le et al. ± 25 % once it is copied) | gate |
| Speed of the front (fit of the mean r_f over the window) | > the fast magnetosonic speed √(v_A² + c_s²) | gate |
| Largest ring mean of B_z / B0 in the window | > 1.5 | gate |
| Mean r_f of the axes against the diagonals, every sample of the window | within 5 % | gate |
| Compression against the perpendicular MHD shock, distance of front and debris, extent of the debris along and across B0, potential and charge of the body, energy budget, fill | — | reported |

`("blast_wave", "gate")` joins `scripts/gpu_gates.py`; the gates run on a clean commit with 8 GPUs,
FP64 and `--precision mixed`, and are recorded in `bench/results/`.

## Phase 4 — production run and documents

`full.yaml` once, recorded with `scripts/record_example.py` (reported, not gated). README with the
citations and the differences to them; wiki page; `docs/roadmap.md` (R6); `PROGRESS.md`.

## Risks and the measurement that decides each

| Risk | Decided by |
|---|---|
| `flux` faces with a magnetised plasma at rest (tested so far with a wind only) | N1 |
| Charging of the body in a plasma at rest and under the debris | N2; potential at the gate |
| The cubic box and the grid break the symmetry | the axes-against-diagonals criterion |
| The shock does not leave the debris (R_m / ρ_i = 2.67; Keenan et al. 8) | distance of front and debris: reported at the gate, judged at the production run |
| Grid heating in the dense shell (dx / λ_D = 20 at the peak) | a 128³ trial, K = 1 against K = 4 |
| Capacity | the fill measured at the gate, before the production size is fixed |

## Invariants

- Gauss and div B: unchanged contracts and checked region (open faces with filter and zone, the
  ball of the body). The debris is loaded neutral cell by cell (shared positions), so div E = ρ = 0
  holds at step 0 without a field solve.
- Loading depends only on (seed, global cell, slot): the block of a radial beam is given in global
  cells and cut to the rank.
- Collectives: the neutrality check and the removal of zero-weight particles are reached by every
  rank; error flags are reduced before any rank raises.
- No divergence cleaning is added. `pypic/` stays experiment-agnostic: radial beams are a general
  load, the shell and its numbers live in the example.

## Verification

- CPU: `ruff check .`, `python scripts/check_kernels.py`, `pytest tests/cpu`.
- GPU: `pytest tests/gpu`, `python scripts/gpu_gates.py --gpus 8` (FP64 and mixed) on a clean commit.
- Example: `python examples/blast_wave/run.py --config … --plan-only`; the gate run and
  `analyze.py`.

## Out of scope

Collisions and ionisation of the debris; a conducting or emitting body; an oblique dipole or a
moving body; the implicit scheme (external fields and injection are explicit only).
