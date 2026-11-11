# Blast wave around a magnetised body (`examples/blast_wave`) and radial beams (PIC-D17)
Status: Done (2026-11-11). Phase 4 done: the production run is recorded and the wiki page is built locally; the wiki and the branch are not pushed (Phase 4, Result). Phase 3 done: the gate records passed on 9bfa1a5 in FP64 and with `precision: mixed` (Phase 3, Result of the records on 9bfa1a5). The FP64 gate record passed on ae4bb69 with `blast_wave:gate` in it (Phase 3, Result of the gate record) and PIC-D17 is Accepted. The first run of the gate configuration passed all seven criteria (Phase 3, Result of the first run). Phase 2 done: the example with planner, analysis and CPU tests. Phase 1 done: radial beams (PIC-D17, Proposed) gated on 574c904. Phase 0 done: N1 passed, N2 failed the charge criterion, which the owner then made a reported value (Phase 0, Amendment)
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

Correction (2026-11-11, after reading the papers again for the documents of Phase 3; the list above
is left as it was written):

- Le et al. 2021 present 2-D hybrid simulations, not 3-D ones. Keenan et al. is Phys. Plasmas 29,
  012107 (2022) and is 3-D.
- Le et al. 2021 do not say that the shock leaves the debris when ρ_d < R_m. They give ρ_d < R_m as
  the regime in which the debris couples to the ambient plasma through the fields and as the range in
  which their Eq. 4 applies, and ρ_d > R_m as the regime in which it decouples.
- Their runs differ from this example in three more points than the list of differences below
  names: a debris-to-ambient ion mass ratio of 3 (here 1), a Gaussian cloud of width 0.1 R_m (Keenan
  et al.: 4 d_i), and R_m of 10 to 150 d_i (Keenan et al.: 40 d_i; here 8 d_i).
- No decision and no gate criterion rests on the corrected sentences: the configuration keeps
  ρ_d / R_m = 0.375 of Keenan et al., and whether the front leaves the debris is measured (reported).

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

### Changes while building (2026-11-11)

- `fast_speed` moves to `pypic/setup/mhd.py` with the four functions named above (the front
  criterion needs it). The planner report of the bow shock is the same before and after the move,
  for `gate.yaml` and `full.yaml`.
- `MAX_CAPACITY` (2³¹ − 1 particle slots) moves from `pypic/particles.py` to `pypic/config.py`: the
  planner compares against the engine's constant instead of a copy, and `pypic/particles.py` needs
  CuPy.
- The synthetic test recovers the front and the debris extent (the gated quantities); the cavity
  radius is a reported value since the amendment of Phase 3 and has no test of its own.
- Planner figures for 8 GPUs: gate 21.5 GB per GPU, production 66.2 GB per GPU and 2800 steps (the
  table of the approach estimated 20 and 58 GB from the bow-shock records).

### Trial of the beams per debris species (2026-11-11; working tree on b4e2c7c, 8 GPUs, FP64, one run each)

`gate.yaml` with `rm=4 body=false half_width=1.3` (R_m = 16 d_e = 4 d_i, R_B = 33.3 d_e, 176³
cells, 880 steps; the risk table said 128³, the planner's smallest box for these scales is 176³),
`debris=false`, `debris_beams=1` and `debris_beams=4`. With the uniform external field only, the
energy budget has no source term, so its change is numerical.

| | no debris | K = 1 | K = 4 |
|---|---|---|---|
| change of the budget over the run | 289.7 | 372.1 | 315.1 |
| excess over the run without debris / kinetic energy of the debris at step 0 (5808) | — | 1.42 % | 0.44 % |
| electric field energy at the end | 34.3 | 49.3 | 44.0 |
| Gauss's law, largest | — | 2.7e-15 | 5.8e-15 |

K = 4 stays (the value of the approach): it lowers the numerical energy gain that comes with the
debris by a factor 3.2 and costs 4 % more particles.

The gate analysis ran on both runs as a first exercise on real data. It is not a gate: the
configuration has ρ_d = 0.75 R_m (the gate 0.375 R_m). Values of K = 4, with those of K = 1 in
brackets where they differ: front speed 0.166 c (fast speed 0.077 c), largest ring mean of B_z
2.19 B0, axes against diagonals 0.8 % (0.5 %), ambient plasma ahead of the front within 0.4 %. The
largest ∆98 of the debris ions is 34.1 d_e (33.9) = 2.13 R_m at step 682, against 16.3 d_e from
Eq. 4: this criterion fails in the trial.
∆98 follows the free expansion at V_d up to t ≈ 100 and is larger than the mean front radius
between t ≈ 50 and 185, so the debris ions are not held behind the front there. The free gyration
of a debris ion in B0 carries it 2 ρ_d = 1.5 R_m outwards in the trial and 0.75 R_m at the gate.
The criterion of the gate stays as fixed; the trial does not show what the gate configuration gives.

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

### Amendment before the first run with debris (2026-11-11)

The formula of Le et al. 2021 (arXiv:2109.00583) was copied from the paper, as the criterion on
the cavity asked:

- Eq. 3, the extent of the debris at a time: `∫_0^∆98 n_d dz = 0.98 ∫_0^∞ n_d dz`, n_d the density
  of the debris ions on the cut through the centre perpendicular to the background field; "a metric
  that is based on the debris density rather than the magnetic field profile". `∆max` is the
  largest ∆98 over time.
- Eq. 4: `∆max / R_m ∼ (1 + c0 M_A) / (1 + c1 / M_A²)^(1/D)`, D the number of dimensions, with
  `c0 = 0.03`, `c1 = 2` fitted to 2-D hybrid runs (D = 2, m_d / m_b = 3, Z_d / Z_b = 1, Gaussian
  cloud of the size 0.1 R_m); it "only applies for systems with ρ_d < R_m". For 3-D the paper gives
  no fit; it expects coefficients "of similar magnitude in 2D and 3D" (Section 7).

Three facts change the criterion "largest cavity radius 0.5 to 1.1 R_B":

1. ∆max is defined by the density of the debris ions, not by the field. A limit taken from Eq. 4
   has to be compared with that quantity.
2. The coefficients are a 2-D fit for another mass ratio and cloud; in 3-D only the form is given.
   The band of ± 25 % stays as fixed; the comparison is a check against an extrapolated fit, which
   the README has to say.
3. Eq. 4 with D = 3 gives 1.02 R_m = 32.6 d_e = 0.49 R_B at the gate (M_A = 3, R_m = 32 d_e): the
   placeholder 0.5 to 1.1 R_B began where the formula puts the value.

Replaced criterion (gate): the largest ∆98 of the debris ions over the slice samples lies within
25 % of Eq. 4 with D = 3, c0 = 0.03, c1 = 2 (24.5 to 40.8 d_e at the gate), and it is not the value
of the last sample (then the run did not show a maximum). ∆98 is taken in the plane z = centre from
the ring mean of the debris-ion density, as the line integral of Eq. 3 from the centre to the edge
of the plane (the ring mean replaces the single cut of the 2-D runs; the shell at step 0 gives
15 d_e). The magnetic cavity radius r_c (definition above, unchanged) becomes a reported value;
with the dipole along B0 the static field alone is below 0.5 B0 on the equator inside
2^(1/3) r_dip = 8.1 d_e, which is reported with it.

Definitions fixed while building the analysis, before the first run with debris:

- Rays start at the surface of the body (at the centre without a body) and end at the largest
  radius the plane holds on the axes; samples every dx, bilinear in the slice; the running mean is
  over round(d_i / dx) samples. A ray without a front counts r_f = 0, so a missing front fails the
  symmetry criterion once the window has begun.
- The speed of the front needs two samples in the window; an empty window fails the speed, the
  compression and the symmetry criterion (nan). The ring means of B_z start at the surface of the body.
- "Outside the flux tube of the body": distance from the axis through the centre along B0 of
  2 r_dip or more (the exclusion of the null test N2). |B| is the magnitude of the block mean of the
  total field. The densities and |B| are compared as sums over the blocks ahead of the front,
  against the sums over the same blocks at step 0. "Largest r_f" at a 3-D output is the largest over
  the rays and over the slice samples up to that output: where the front has passed, the plasma does
  not count as undisturbed again when the front is lost later.
- Planner rejections (Phase 2): "capacity" is the particle slots of one species on one rank against
  the kernels' limit 2³¹ − 1, with the engine's formula (PIC-D17); "a shell that overlaps the body"
  is more than 1 % of the shell's particles inside the body radius (0.03 % at the gate); ω_pe at the
  peak of the shell is √(1 + n_pk) (ambient and debris electrons).

### Result of the first run (2026-11-11; 1223a21, clean tree, 8 GPUs, FP64, one run)

`gate.yaml` as committed: 943,087,492 particles at the start, 1760 steps, 165 s (93.5 ms per step,
0.80 ns per particle and step, times 8 GPUs). All seven criteria pass.

| Criterion | Measured | Limit |
|---|---|---|
| Gauss, div B | 6.9e-15, 9.0e-15 | 1e-10, 1e-12 |
| Ambient plasma ahead of the front | within 0.050 % (n_e; n_i 0.044 %, \|B\| 0.015 %); 258,678 blocks at step 0, 145,526 at the end | 5 % |
| Speed of the front (14 samples in the window) | 0.1627 c = 2.10 fast speeds = 0.87 V_d | > 0.0774 c |
| Largest ring mean of B_z / B0 in the window | 2.69 | > 1.5 |
| Axes against diagonals | 0.21 % | 5 % |
| Largest ∆98 of the debris ions | 39.5 d_e = 1.24 R_m at step 616 (of 1760) | 24.5 to 40.8 d_e (Eq. 4: 32.6) |

- ∆98 is 21 % above Eq. 4 with D = 3, 4 points inside the band. This is one run, and the bow-shock
  gate showed differences in the third digit between runs with the same seed; how far ∆98 moves
  between runs has not been measured. After its maximum ∆98 falls to 34.3 d_e at the end.
- The front leaves the debris: it is 10.5 d_e ahead of ∆98 at the end of the window (step 660) and
  43 d_e at the end of the run (front 77.4 d_e = 1.16 R_B, ∆98 34.3 d_e). The largest ring mean of B_z
  falls from 2.6 B0 (step 352) to 1.5 B0 (end).
- Compression: the perpendicular MHD shock gives 2.51 at M_A = 3 and 2.26 at the measured speed of
  the front (M_A = 2.60); the measured 2.69 is the largest ring mean of B_z, not a mean behind a
  plane shock.
- Magnetic cavity: largest radius 33.25 d_e = 0.50 R_B = 1.04 R_m, constant over the last 350 steps
  (the static field alone: 8.1 d_e).
- Debris at the end (3-D density of the debris ions): rms distance 36.8 d_e along B0, 14.8 d_e per
  direction across.
- Energy: the budget changes by +1.95 % of the energy at step 0, which is 10.9 % of the kinetic
  energy of the debris (N1, the ambient plasma alone: +2.2 %). The debris ions keep 20 % of their
  kinetic energy; the ambient ions gain 41 % of the kinetic energy of the debris.
- Body: potential against the plasma at two dipole radii 4.9 T_e / e in the plane across B0 and
  1.3 T_e / e in the plane along it; charge +22.5 at the end (largest 32.9); 106,480 particles
  absorbed. Reported without a limit (Phase 0, Amendment).
- Fill of the particle arrays 67 % at the most (capacity factor 1.5); the particle number falls
  from 943.1e6 to 936.7e6.

### Result of the gate record (2026-11-11; ae4bb69, clean tree, 8 GPUs, FP64)

`python scripts/gpu_gates.py --gpus 8`: ALL PASSED
(`bench/results/2026-11-11-ae4bb69c6-gpu-gates.json`): 93 GPU tests, the five two-stream gates with
the growth rates of 574c904 (+3.0 / −1.7 / −2.3 / +3.5 / −4.7 %, Gauss max 7.6e-12, implicit energy
1.9e-10) and 17 example gates, with `bow_shock:gate` (magnetopause 0.916, shock 1.381, upstream
within 0.11 %) and `blast_wave:gate`. PIC-D17 is Accepted with this record. `precision: mixed` was
not run on this commit.

`blast_wave:gate` in the record is a second run of the configuration of the first run (same seed;
the code of the run is the same, ae4bb69 only registers the gate). Against the first run:

| Quantity | First run (1223a21) | Gate record (ae4bb69) |
|---|---|---|
| Speed of the front, largest ring mean of B_z / B0, axes against diagonals | 0.16274968 c, 2.6865894, 0.21272969 % | the same in all 16 digits |
| Largest ∆98 and its step | 39.5455585 d_e, step 616 | the same in all 16 digits |
| Ambient plasma ahead of the front | 0.0500307 % | 0.0500287 % |
| Gauss, div B | 6.9e-15, 9.0e-15 | 7.9e-15, 9.0e-15 |
| Front at the end | 77.4205320 d_e | 77.4205321 d_e |
| Debris rms distance along B0 | 36.809888 d_e | 36.809882 d_e |
| Body potential (plane across B0) | 4.93716 T_e / e | 4.93719 T_e / e |
| Body charge at the end, particles absorbed | +22.4937968, 106,480 | +22.4937968, 106,480 |
| Particles at the end | 936,674,000 | 936,673,989 |
| Time per step | 93.5 ms | 93.2 ms |

- The sentence of the first-run result "how far ∆98 moves between runs has not been measured" is
  answered for FP64 and this seed: not at all in these two runs. The 21 % above Eq. 4 are therefore
  not a difference between runs. Not shown: another seed, `precision: mixed`, more than two runs.
- The sums over the blocks ahead of the front (relative to step 0) agree in all digits up to step
  704 and differ in the 13th digit at step 880 and in the 8th at step 1760.
- The count "943,087,492 particles at the start" of the first-run result is not the history value
  at step 0, which is 943,088,238 in both runs (175,432 injected and 175 absorbed by then); where
  the first number came from was not traced.

### Result of the records on 9bfa1a5 (2026-11-11; clean tree, 8 GPUs, FP64 and then mixed)

9bfa1a5 adds the movie of the example and two reported ratios; `pypic/analysis/animate.py` changed,
no code of the step. The runbook runs the mixed gates after the FP64 gates on the same commit, so
both were run: ALL PASSED (`bench/results/2026-11-11-9bfa1a537-gpu-gates.json`: 93 GPU tests, five
two-stream gates, 17 example gates; `…-gpu-gates-mixed.json`: four two-stream gates, 17 example
gates). `blast_wave:gate`, same seed as before:

| Quantity | FP64 (ae4bb69) | FP64 (9bfa1a5) | mixed (9bfa1a5) |
|---|---|---|---|
| Speed of the front | 0.16274968378855556 c | the same | 0.16274967880523641 c |
| Largest ring mean of B_z / B0 | 2.6865894164931907 | the same | 2.6865894211811967 |
| Axes against diagonals | 0.21272968616 % | the same | 0.21272970359 % |
| Largest ∆98, its step | 39.5455585300089 d_e, 616 | 39.54555853000368 d_e, 616 | 39.54555829703839 d_e, 616 |
| Ambient plasma ahead of the front | 0.0500287 % | 0.0500293 % | 0.0500666 % |
| Gauss, div B | 7.9e-15, 9.0e-15 | 6.2e-15, 8.8e-15 | 6.9e-15, 8.6e-15 |
| Energy budget change / debris kinetic energy | 10.899889 % | 10.899890 % | 10.897725 % |
| Body charge at the end, particles absorbed | +22.4937968, 106,480 | +22.4937968, 106,480 | +21.7170596, 106,358 |
| Largest `fill` | 0.6701040 | 0.6701040 | 0.6701177 |
| Time per step | 93.2 ms | 93.2 ms | 84.8 ms |

- Three FP64 runs with the same seed: the front criteria are equal in all digits, ∆98 in 12 digits.
  This replaces "the same in all 16 digits" of the two-run comparison above for ∆98.
- `precision: mixed` passes the seven criteria; the front speed, the ring mean and ∆98 have 7 digits
  or more in common with FP64. The 21 % above Eq. 4 are in both precisions. Not shown: another seed.
- `bow_shock:gate` in these records: magnetopause 0.923 (FP64) and 0.915 (mixed), shock 1.381 and
  1.382; the ae4bb69 record had 0.916 and 1.381 in FP64.

## Phase 4 — production run and documents

`full.yaml` once, recorded with `scripts/record_example.py` (reported, not gated). README with the
citations and the differences to them; wiki page; `docs/roadmap.md` (R6); `PROGRESS.md`.

### Result (2026-11-11; 3752d32, clean tree, 8 GPUs, FP64, one run)

Record: `bench/results/2026-11-11-3752d3211-example-blast_wave-full-run-ea70ccc5.json`. m_i/m_e = 25,
576³, 3,130,149,797 particles at step 0, 2800 steps, 788 s (281.5 ms per step), fill 0.668.

| quantity | production | gate configuration (FP64) |
|---|---|---|
| front speed in the window | 0.13001 c = 2.100 fast speeds | 0.16275 c = 2.103 |
| largest ring mean of B_z / B0 | 2.685 | 2.687 |
| axes against diagonals | 0.38 % | 0.21 % |
| largest ∆98 | 49.74 d_e = 1.243 R_m at step 980; 1.220 × Eq. 4 | 39.55 d_e = 1.236 R_m at step 616; 1.212 × Eq. 4 |
| ambient ahead of the front | 0.045 % | 0.050 % |
| Gauss, div B | 8.0e-15, 1.3e-14 | 6e-15 to 8e-15, 9e-15 |
| cavity, largest radius | 41.75 d_e = 0.502 R_B | 33.25 d_e = 0.500 R_B |
| front at the end | 98.2 d_e = 1.18 R_B | 77.4 d_e = 1.16 R_B |
| front ahead of ∆98 at the end of the window | 13.6 d_e (step 1050) | 10.5 d_e (step 660) |
| debris rms along / across B0 at the end | 49.2 / 18.3 d_e | 36.8 / 14.8 d_e |
| debris ions keep / ambient ions gain (of the debris energy) | 22.7 % / 39.6 % | 20.4 % / 40.9 % |
| budget change, of the debris energy | 25.9 % (2.82 % of the energy at step 0) | 10.9 % (1.95 %) |
| body: absorbed, charge at the end (largest) | 248,006, +25.8 (53.9) | 106,480, +22.5 (32.9) |
| body potential, across / along B0 | 6.1 / 1.1 T_e / e | 4.9 / 1.3 T_e / e |

- All seven criteria of the gate hold; they are reported here, not gated.
- Risk "the shock does not leave the debris" (table below), judged here as the plan says: the front
  leaves the debris. At the end it is at 98.2 d_e while ∆98 has fallen to 42.9 d_e; the ring mean of
  B_z has fallen to 1.5 B0 by then.
- The excess of ∆98 over Eq. 4 is 22 % (gate 21 %). The configurations differ in mass ratio, V_d / c
  and box width and share R_m / d_i, the shell shape and the seed. Not shown: the excess at another
  R_m / d_i, shell shape or seed.
- The budget change is a larger fraction of the debris energy than at the gate. The energy at step 0
  is 9.2 times the debris energy (gate 5.6) and the run has 1.6 times the steps. No run without
  debris was made at this size, so the numerical heating of the ambient plasma is not separated.
- Wiki: `docs/wiki/{en,ja}/Blast-Wave.md` and the Home rows; data and images from a fourth FP64 run
  of the gate configuration (3752d32, run 4f8a4615): front speed, ring mean and symmetry equal to
  the gate records in all digits, ∆98 equal to the ae4bb69 record in all digits, other values within
  4e-5. Rendered into the local wiki checkout; not published.
- The count "943,087,492 particles at the start" of the first-run result (Phase 3) is the number on
  the launcher's start line. It is not the first history value (gate: 943,088,238; production:
  3,130,150,010 on the start line, 3,130,149,797 in the history); the differences were not traced.

## Risks and the measurement that decides each

| Risk | Decided by |
|---|---|
| `flux` faces with a magnetised plasma at rest (tested so far with a wind only) | N1 |
| Charging of the body in a plasma at rest and under the debris | N2; potential at the gate |
| The cubic box and the grid break the symmetry | the axes-against-diagonals criterion |
| The shock does not leave the debris (R_m / ρ_i = 2.67; Keenan et al. 8) (corrected 2026-11-11: R_m / ρ_d = 2.67 here and in Keenan et al.; what differs is R_m, 8 d_i here and 40 d_i there) | distance of front and debris: reported at the gate, judged at the production run |
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
