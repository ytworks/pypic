# Exploratory runs: three conditions as close to a supernova remnant as the engine allows

Topic `spherical-shock-injection`. Designed 2026-11-12. **Status: designed and checked with the planner; not run.**
Japanese version: `exploratory-conditions.ja.md` (this English file is the original).

**Corrected 2026-11-13** after the review of these conditions (issue #48, kept in
`reviews/2026-11-13-exploratory-conditions.md`): Sects. 2.2 and 2.3 now state the speed and the electron
temperature behind each remnant value, and only β of run L is listed as matched. The physical conditions of
L, M and H are unchanged.

**Superseded in part, 2026-11-13:** Sect. 6 (checks 2 and 6) and Sect. 8 are replaced by
`exploratory-protocol.md`, which is fixed before the first run. The run time of Sect. 7 was an estimate; the
16-step check measured 0.9–1.1 s per step, which is about 45–50 minutes per run. The sections stay below as
designed.

Citekeys refer to `research/literature/`; "survey" is `survey-snr-electron-injection.md`. "Our arithmetic" marks
values derived here. "Estimate" marks an extrapolation that no run has tested.

## 1. Purpose and status

- **Purpose.** These runs are trials to build hypotheses about how electrons are accelerated at a spherical
  shock. They are not evidence for a claim. After their analysis the research question of the topic is to be
  narrowed (the owner's decision).
- **Brief from the owner (2026-11-12):**
  - conditions as close to a supernova remnant (SNR) as possible;
  - a large box, few particles, a short run time;
  - two or three runs;
  - nothing is run at this stage.
- **Decisions of the owner (2026-11-12):**
  - only what the existing example `examples/blast_wave` can express, so no new code;
  - the runs differ in the field strength only;
  - 768³ cells with 2 particles per cell and species.
- **What exists now:** three configuration files in `configs/` that the planner of the example accepts
  (Sect. 10), and this record. No run has been made, on a GPU or otherwise.
- **The most important limit:** the front of such a blast is slower than the debris and slows down. Run H is
  nominally above the threshold for surfing, but its front is expected to be at or below it (Sect. 4.3).

## 2. How close to a supernova remnant

### 2.1 What limits the choice
With the units of the engine (c = 1, lengths in d_e = c/ω_pe) and the parameters of the example:

- M_A = (V_d/c) · (ω_pe/Ω_ce) · (m_i/m_e)^(1/2), and β = 4 (v_te/c)² (ω_pe/Ω_ce)² for T_i = T_e
  (`examples/blast_wave/plan.py`).
- **Debye length.** The cell is 0.5 d_e, and the engine's practice is a Debye length of 0.2 cells, that is
  v_te = 0.1 c. Colder electrons would heat numerically (Sect. 4.3). With v_te fixed, a higher M_A at the same
  flow speed means a weaker field and so a higher β.
- **Box.** The planner needs a half-width of at least 1.3 R_B, with R_B = R_m M_A^(2/3) the magnetic stopping
  radius. The debris couples to the ambient plasma through the field only if its gyroradius M_A d_i does not
  exceed the equal-mass radius R_m (Le et al. 2021, as quoted in `examples/blast_wave/README.md`). At 768³
  cells and m_i/m_e = 16 this allows M_A up to about 8.7.
- **Run length.** A run lasts 1.25 R_B/V_d, which is 1.25 (R_m/d_i) M_A^(−1/3) ion gyro-times: the higher M_A,
  the fewer gyro-times.

### 2.2 The numbers against a remnant

The Mach numbers of the runs are nominal: they use the drift speed V_d of the debris, not the speed of the front.

| Quantity | Supernova remnant | L | M | H |
|---|---|---|---|---|
| Alfvén Mach number M_A | about 250–300 with the unamplified field, about 15 with the amplified field (5000 km/s) | 3.5 | 6.0 | 8.7 |
| M_A / (m_i/m_e)^(1/2); whistler-mediated injection needs ≲ 1 | 6–7 unamplified, 0.35 amplified | 0.88 | 1.5 | 2.2 |
| M_A / [1.2 (m_i/m_e)^(2/3)]; surfing needs ≥ 1 | 1.4–1.7 unamplified, 0.08 amplified | 0.46 | 0.79 | 1.14 |
| (M_A/√β_e) / [(m_i/m_e)^(1/2)/2] = √2 × (flow speed / electron thermal speed); whistlers from reflected electrons need ≳ cos θ_Bn | about 18 at 5000 km/s and 10⁴ K; 3.2 at 3000 km/s and T_e = 10 eV | 3.5 | 3.5 | 3.5 |
| Sonic Mach number | about 440 (Tycho, 4990 km/s) | 5.5 | 5.5 | 5.5 |
| Fast Mach number | — | 2.95 | 4.05 | 4.64 |
| β | 0.4–0.8 at 10⁴ K; about 5–10 at 10 eV | 0.49 | 1.44 | 3.03 |
| Shock speed / c | 0.007–0.02 | 0.25 | 0.25 | 0.25 |
| Part of the sphere that is superluminal | under 2 % | 25 % (θ_Bn > 76°) | 25 % | 25 % |
| Flow speed / electron thermal speed | about 13 at 5000 km/s and 10⁴ K; 2.3 at 3000 km/s and T_e = 10 eV | 2.5 | 2.5 | 2.5 |
| m_i/m_e | 1836 | 16 | 16 | 16 |
| ω_pe/Ω_ce | about 350 | 3.5 | 6.0 | 8.7 |
| Run length in ion gyro-times | — | 13.2 | 7.7 | 5.3 |
| Debris gyroradius / R_m | — | 0.22 | 0.54 | 1.00 |

Sources for the remnant column:
- M_A: `morlino2012strong` Sect. 2.3; `riquelme2011electron` Sect. 5; survey Sect. 2. The values are for
  about 5000 km/s.
- Thresholds: `riquelme2011electron` Sect. 3.1 (whistlers); `matsumoto2012electron` Eq. 8 (surfing);
  `amano2010critical` Eqs. 2, 3 (whistlers from reflected electrons). Eq. 3 gives
  M_A/√β_e ≃ 68 (V_s / 3000 km/s) (T_e / 10 eV)^(−1/2), and the threshold at the real mass ratio is 21.4. So
  the row is 3.2 for the paper's 3000 km/s and 10 eV, and about 18 for 5000 km/s and 10⁴ K (0.86 eV), the
  values behind the rows for β and for the flow speed (our arithmetic).
- This row and the row of the flow speed over the electron thermal speed are one quantity: with
  β_e = 2 (m_e/m_i) (v_te/v_A)² and v_te = (kT_e/m_e)^(1/2), M_A/√β_e = (m_i/m_e)^(1/2) V/(√2 v_te), and so
  the row is √2 times that ratio at any mass ratio (our arithmetic).
- Sonic Mach number: `morlino2012strong` Sect. 2.3 (M_0 = 443 at 4990 km/s).
- β: our arithmetic for ions and electrons at 10⁴ K, 0.8 for Tycho (survey Sect. 2) and 0.4 for the density
  and field of SN 1006. At 10 eV it is 11.6 times larger.
- ω_pe/Ω_ce: our arithmetic in the survey, Sect. 2.
- Shock speed: `vink2012supernova` Sect. 11.3. Superluminal part: survey Sect. 2.
- Flow speed over electron thermal speed: our arithmetic. v_te is 390 km/s at 10⁴ K and 1330 km/s at 10 eV.

The values for L, M and H are the planner's (Sect. 10) and our arithmetic from them.

### 2.3 What is matched and what is not
- **Matched:** β in run L, for a remnant with an upstream temperature of 10⁴ K. Of the remnant's numbers this
  is the only one that is matched.
- **Alike in kind, not a matched number:**
  - the initial geometry: a spherical drive in a uniform field, so that the angle between the radial
    direction and the field runs from 0° to 90°. This describes the set-up, not a shock: no shock is expected
    near the poles (Sect. 4.2), and θ_Bn has to be taken from the measured front where a shock exists;
  - the side on which the nominal M_A lies against the whistler limit and the surfing threshold. L lies below
    the whistler limit, as a remnant with an amplified upstream field does. H lies nominally above the surfing
    threshold, as a remnant with the unamplified field does. M lies between the two. The fronts are slower
    than the nominal values (Sect. 4.3).
- **Not matched, in every run:**
  - the Mach numbers themselves, the sonic Mach number above all;
  - the mass ratio and the shock speed, with a superluminal band that is more than ten times too wide;
  - the ratio of flow speed to electron thermal speed, and with it M_A/√β_e against the threshold of
    `amano2010critical`, which is the same quantity (Sect. 2.2). The runs have 2.5 (3.5 in units of the
    threshold). A remnant at 5000 km/s and 10⁴ K, the temperature at which its β is 0.4–0.8, has about 13
    (18). It has 2.3 (3.2) only at 3000 km/s and T_e = 10 eV, where its β is about ten times larger. The value
    is the same in all three runs because the flow speed and the temperature are fixed;
  - the size of the system in units of the shock thickness.
- **Not matched in M and H:** β, which rises to 1.4 and 3.0 because the temperature is fixed. The scan of M_A
  is therefore also a scan of β.

## 3. The three conditions

### 3.1 Common to all three
- **Driver:** the debris shell of `examples/blast_wave`. It is a Gaussian shell at 0.33 R_m, of width 0.1 R_m and
  peak density 16.5 n₀, with a radial drift V_d = 0.25 c, loaded as 4 beams per species.
- **Ambient plasma:** at rest, T_i = T_e, v_te = 0.1 c, m_i/m_e = 16, injected through the six open faces.
- **Field:** uniform B₀ along z, as an external field. No absorbing sphere and no dipole (`body: false`): a
  remnant has no magnetised body.
- **Grid:** 768³ cells of 0.5 d_e, so 384 d_e = 96 d_i across. dt = 0.25/ω_pe; 2960 steps.
- **Particles:** 2 per cell and species: 1.81 × 10⁹ ambient particles, plus 1.2 × 10⁷ to 7.5 × 10⁷ of the debris.
- **Numerics:** explicit solver, `precision: compact`, 4 filter passes with compensation, Marder zone 0.125,
  seed 12345. On 8 GPUs the rank grid is (1, 2, 4).
- **R_B** is about 148 d_e = 37 d_i in all three (set through `rm`), so the box and the run length in 1/ω_pe are the
  same.

### 3.2 What differs

| Run | Config | `mach_alfven` = `wpe_wce` | `rm` (d_i) | R_m (d_e) | v_A / c | β | Run length (1/Ω_ci) |
|---|---|---|---|---|---|---|---|
| L | `configs/explore-L.yaml` | 3.5 | 16.0 | 64 | 0.0714 | 0.49 | 13.2 |
| M | `configs/explore-M.yaml` | 6.0 | 11.2 | 44.8 | 0.0417 | 1.44 | 7.7 |
| H | `configs/explore-H.yaml` | 8.7 | 8.7 | 34.8 | 0.0287 | 3.03 | 5.3 |

### 3.3 Why a scan of the field strength
- The drift speed and the temperature are the same in all three, so the sonic Mach number is the same and only
  the field changes. This is the difference between a remnant shock with and without an amplified upstream
  field (survey Sect. 4.1: whistlers need M_A ≲ 20, surfing needs M_A above the trapping threshold).
- The survey lists as a disagreement which Mach number SNR shocks need for injection (Sect. 8.1, item 2). The
  three runs are placed on each side of the two thresholds and between them, by their nominal Mach numbers.
- Colder electrons at higher M_A would have kept β fixed. That was not chosen, because the Debye length would
  fall to 0.1 cells and the upstream temperature could not be held (Sect. 4.3).

## 4. What the runs can and cannot show

### 4.1 Can show (first look, for hypotheses)
- The shape of the front and of the field around the sphere, and the local θ_Bn from the measured front normal.
- Where on the front energetic electrons appear, as a function of θ_Bn and time, and how this changes from L
  to H.
- Electron spectra by position in the meridional and equatorial planes.
- Along tracked electrons: the energy history and the work of the field parallel and perpendicular to B.
- Whether electrons run ahead of the front, and at which θ_Bn.

### 4.2 Cannot show
- **Surfing on Buneman waves.** The waves have about 6 cells per wavelength (2π · 2V_d/ω_pe ≈ 3.1 d_e; our
  arithmetic), are partly damped by the filter, and sit in the noise of 2 particles per cell. Whether such waves
  appear can be seen; how well they trap electrons cannot. An absence of surfing in these runs says nothing
  (survey Sect. 7).
- **A quasi-parallel shock.** Ion-driven waves at quasi-parallel shocks need tens of ion gyro-times, and a
  foreshock about 125 or more (`morris2023preacceleration` Sect. 4.5); the runs last 13, 7.7 and 5.3. Along B₀
  the debris streams freely (`examples/blast_wave/README.md`) and reaches the z faces near the end. Near the
  poles the runs show a beam running into a plasma, not a shock.
- **Diffusive shock acceleration.** The box is far too small for it (the thresholds are in the survey, Sect. 2).
- **The dependence on mass ratio, shock speed or β at fixed M_A.** Only the field strength varies, and β
  varies with it.
- **Electron reflection where θ_Bn > 76°.** That band is superluminal at 0.25 c. It narrows as the front slows.

### 4.3 Risks
- **The front is slower than V_d and slows down.** In the gated run of the example the front moved at 0.87 V_d
  between 0.3 and 0.7 R_B (`docs/plans/2026-11-11-blast-wave.md`, results of the gate), and more slowly beyond.
  - At 0.87 V_d the front of H has M_A = 7.6, which is the surfing threshold and not above it (our arithmetic).
  - The front slows on the scale of R_m, and R_B is 4.2 R_m in H against 2.1 in that gated run. So the front of H
    is expected to be well below the threshold for most of the run; an extrapolation gives a front M_A of about
    3 for L, 4–5 for M and at most about 5 for H (estimate).
  - Each run is to be labelled by the measured Mach number of its front, not by the nominal one.
- **Coupling in H.** The debris gyroradius equals R_m. Le et al. 2021 give gyroradius < R_m as the condition for
  coupling. In a trial of the example at 0.75 R_m the debris ions were not held behind the front
  (`docs/plans/2026-11-11-blast-wave.md`, trial with one and four beams). In H, energetic electrons inside
  about 2 R_m = 70 d_e may come from debris running through the ambient plasma, not from a shock.
- **Numerical heating of the upstream electrons.** The null run of the example, with 8 particles per cell, ended
  after 1760 steps with T_e up by 1.5 % in the interior and 3.0 % in the worst slab at a face, and the total
  energy up by 2.2 % (same plan, null test N1). If the heating scales as the inverse of the particle number,
  these runs end with T_e up by about 10 % in the interior and 13–20 % near the faces (estimate). It is
  measured in each run (Sect. 6).
  - The energy gained numerically by the ambient plasma is then comparable to the energy of the debris, or
    larger, in M and H (estimate). The budget figure of the example, which is in units of the debris energy,
    will look large for that reason.
- **Noise.** With 2 particles per cell the fields are noisy. Tails of spectra need counts (Sect. 6).
- **The last tenth of the run.** A front at 1.17 R_B, where the gated front ended, is 19 d_e from the faces
  (our arithmetic). The end of each run is affected by the boundary.
- **Untested paths of the output.** The selection of particles by energy (`select`) has no test and no example
  in the repository, and a dump that selects nothing, as at step 0, writes empty datasets. A 16-step check on
  the GPU machine comes before the first run (Sect. 8).

## 5. Outputs

The planner's own output has no particle data, so each config gives the list of diagnostics itself. The first
three repeat the planner's diagnostics under the same names, because `examples/blast_wave/analyze.py` reads them.

| Name | Kind | Content | Every (steps) | Outputs | Size, uncompressed |
|---|---|---|---|---|---|
| `coarse` | fields | moments of the four species, E, B, ρ in blocks of 13 / 9 / 7 cells (L / M / H) | 296 | 11 | 0.26 / 0.76 / 1.6 GB |
| `slice_z` | slice | equatorial plane at full resolution: densities, ion and electron flux, electron energy density, E, B, J | 74 | 41 | 1.9 GB |
| `slice_y` | slice | meridional plane (contains B₀) at full resolution, same fields | 74 | 41 | 1.9 GB |
| `spec_meridional` | histogram | all four species in a slab of 16 d_e around the meridional plane: (x, z, kinetic energy), 192 × 192 × 96 bins, energy from 10⁻⁴ to 10² m_e c², 16 bins per decade | 148 | 21 | 1.2 GB |
| `spec_equatorial` | histogram | the same around the equatorial plane: (x, y, kinetic energy) | 148 | 21 | 1.2 GB |
| `spec_all` | histogram | energy spectra of the whole box, all four species | 74 | 41 | under 1 MB |
| `spec_up` | histogram | energy spectra of the ambient species in a corner column that the blast never reaches | 74 | 41 | under 1 MB |
| `hot_e` | particles | every electron (ambient and debris) with γ ≥ 1.5: position, momentum, weight, id | 148 | 21 | see below |
| `tracks_e` | tracks | 1.9 × 10⁵ ambient electrons (ids divisible by 4733) followed from the start, with E and B at the particle | 4 | 741 | 6.8 GB |

- **`hot_e`:** γ = 1.5 is a kinetic energy of 0.5 m_e c², the ram energy of an ion (½ m_i V_d²) and 50 times the
  upstream electron temperature.
  - One mirror reflection of an electron at rest gives 0.13, 0.29, 0.67 and 2.3 m_e c² at θ_Bn = 0°, 45°, 60° and
    70° (our arithmetic, relativistic). So γ ≥ 1.5 marks one reflection at 56° or more, or several interactions
    at smaller angles.
  - The planner's table gives 16 GB per output. That is its flat planning figure of one tenth of all
    particles. The expected number is between 10² and a few 10⁶ per dump, depending on the downstream electron
    temperature (estimate).
  - All selected particles are gathered on one GPU. If a dump exceeds about 10⁷ particles, the threshold is
    raised to γ ≥ 2.
- **`tracks_e`:** the set is fixed at the first output.
  - The stride 4733 is a prime that does not line the tracked electrons up on the grid; a stride that shares a
    factor with 768 would put them in planes.
  - About 7 × 10⁴ of the tracked electrons lie where the front passes (our arithmetic). If one in a thousand is
    accelerated, about 70 tracks show it. `hot_e` is the complete census.
  - One sample every 4 steps gives about 5, 9 and 14 samples per electron gyration (L, M, H) in a field
    compressed fourfold (our arithmetic).
  - Its cost in run time is not known and is to be measured in the first run.
- **`spec_up`:** the column is at a cylindrical radius of 204 d_e or more from the axis along B₀ and holds
  5.8 × 10⁶ electrons. It gives the upstream electron temperature against time and shows any numerical tail.
- Scalars (energies, Gauss's law, div B) are written every 37 steps.
- The analysis of the electron data is not written yet. It will be `analyze.py` in the topic folder.

## 6. Checks fixed before the runs

These decide whether a run can be used for building hypotheses. They are not physics claims. Names in
`code font` are results of `examples/blast_wave/analyze.py`.

1. **Conservation contracts of the engine:** `gauss_max` ≤ 10⁻¹⁰ and `divB_max` ≤ 10⁻¹².
2. **A shock formed at the equator:** between 0.3 and 0.7 R_B there are at least two samples of the front
   (`window_samples`), `front_speed` exceeds `fast_speed`, and `compression` exceeds 1.5. These are the criteria
   of the example's gate. A run that fails is reported as "no shock"; its electron data are not read as shock
   acceleration.
3. **The front is ahead of the debris:** `front_debris_distance` > 0. If not, the run counts as debris running
   through the ambient plasma.
4. **The upstream is undisturbed:** `ambient_deviation` ≤ 0.05.
5. **Symmetry in the equatorial plane:** `asymmetry` ≤ 0.05.
6. **Numerical heating:** the mean electron energy in `spec_up` stays within 1.25 of its value at step 0. A run
   that exceeds this is not used for statements that depend on β or on the electron thermal speed.
7. **Particle arrays:** `fill_max` < 0.95.
8. **Counts:** the number of electrons with γ ≥ 1.5 is reported per 10° band of θ_Bn. A band with fewer than 100 is
   not used for a spectrum above that energy.
9. **Reported, not a criterion:**
   - the speed of the front in units of V_d and of v_A (the regime of the run);
   - the largest extent of the debris. The example's own verdict `passed` includes a band for that extent, which
     is not expected to hold for M and H;
   - `energy_budget_change`.

## 7. Cost

- **Device memory:** the planner gives up to 57.9 / 55.9 / 55.3 GB per GPU (L / M / H) on 8 GPUs with
  `precision: compact`, and 63.7 / 61.5 / 60.8 GB with `fp64`.
  - The planner's figure does not count the diagnostics of the config. The histograms need about 7 GB of
    temporary arrays while they are computed (estimate), which is within the planner's own allowance for output.
  - The GPUs have 103 GB each (`docs/plans/2026-11-10-bow-shock.md`), and the largest recorded use is 79.8 GB.
    Both precisions fit; `compact` is chosen for speed.
- **Run time:** about 15–20 minutes per run on 8 GPUs. This is an estimate, not a measurement.
  - It scales the phases of the example's gate run (384³, 9.4 × 10⁸ particles, 73–91 ms per step) to 8 times
    the cells and 1.9 times the particles, which gives about 285 ms per step in `compact` and 320 ms in `fp64`.
  - It then adds 10–25 % for the output.
  - The grid costs more than the particles here, because there are only 2 particles per cell.
- **Disk:** about 13–15 GB per run uncompressed, without `hot_e`.
- **The GPUs are shared.** They must be free for a run; see `docs/gpu-runbook.md` for the containers that have to
  be stopped with the owner's agreement.

## 8. How to run later

Not to be done at this stage. When the owner starts the runs:

1. Work on a clean commit of the topic branch. The run-recording tools of the research harness
   (`research/tools/run.py`, `record_run.py`) are planned and not built. Until they exist, a run is not
   recorded as evidence; `run.json` in the output directory still holds the commit and the run id.
2. Check that the GPUs are free (`docs/gpu-runbook.md`).
3. Make a 16-step check first. It exercises the start at 768³ and every diagnostic at step 0, the empty
   particle dump included. If that dump fails, add `--set output.codec=none`.
4. Run L first. It is the nearest to the gated regime of the example, and it gives the first measurements of
   run time, memory, numerical heating and the size of `hot_e`. Then H, then M.

```bash
# 16-step check
python -m pypic.launch -n 8 examples/blast_wave/run.py \
  --config research/topics/spherical-shock-injection/configs/explore-L.yaml \
  --set overrides.n_steps=16 --out output/research/spherical-shock-injection/check-L
# the run and the example's analysis
python -m pypic.launch -n 8 examples/blast_wave/run.py \
  --config research/topics/spherical-shock-injection/configs/explore-L.yaml \
  --out output/research/spherical-shock-injection/explore-L
python examples/blast_wave/analyze.py output/research/spherical-shock-injection/explore-L
```

The plan alone, on any machine without a GPU (it writes nothing; it needs numpy, pyyaml and scipy):

```bash
PYPIC_WORLD_SIZE=8 PYTHONPATH=. python examples/blast_wave/run.py \
  --config research/topics/spherical-shock-injection/configs/explore-L.yaml --plan-only
```

The planner checks the names of the diagnostics and of their keys. The histogram settings and the `select` key
are checked only when a run starts; they were checked here by hand against `pypic/output/diagnostics.py`.

## 9. Left out, and what it would need

- **A front that keeps its speed (piston mode).** A much larger debris mass would hold the Mach number of the
  front. The planner ties the box to R_B, which then becomes too large; a planner of the topic's own would be
  needed. This is the direct answer to the first risk of Sect. 4.3.
- **M_A of 20–30.** Possible in a thin, quasi-two-dimensional box with the field in the plane: a cylindrical
  front that still has every θ_Bn. It also needs a planner of its own, and it leaves the sphere.
- **Cold electrons** (β and the electron thermal Mach number of a remnant): cells of 0.2 d_e or less, which a
  three-dimensional box of this size cannot afford.
- **Other mass ratios and shock speeds.** These are needed before any statement about scaling (survey Sect. 7).
- **Planar runs at the same local parameters**, to separate what the curved front adds.
- **Test electrons placed where they are needed**, in place of the stride subset of `tracks_e`; spectra binned in
  radius and angle; a selection by region and energy together.

## 10. The planner's output for the three configs

Produced on 2026-11-12 with the last command of Sect. 8 (CPU only; the plan is for 8 GPUs).

**L**

```text
=== pypic blast-wave planner (debris shell around a magnetised body) ===
ambient: m_i/m_e=16, B0=0.285714 along z, beta=0.49, v_te=0.1 c (Debye length 0.20 dx), v_A=0.0714 c, fast speed 0.0848 c
debris (on): M_A=3.5 (fast Mach 2.95), V_d=0.2500 c, R_m=64 d_e = 16 d_i, R_B=147.5 d_e, gyroradius 14 (0.219 R_m); shell at 21.12, width 6.40, peak density 16.5 in 4 beam(s) per species (omega_pe dt <= 1.05, 0.00 % inside the body); extent after Le et al. 2021, Eq. 4 with D = 3: 67.2 d_e
body: OFF (no sphere, no dipole)
grid (768, 768, 768), dx=0.5, L=384 d_e (half-width 1.30 R_B), dt=0.25, 2960 steps = 1.25 R_B/V_d = 13.2 / Omega_ci
1,811,939,328 ambient particles (2 per cell and species), about 74,672,289 of the debris (block of 8,998,912 cells), capacity factor 1.5; memory up to 57.9 GB/GPU on 8 GPU(s), rank grid (1, 2, 4); ids 6.0e-03 of 2^53
faces: six open faces with flux injection of the ambient species; filter 4+comp, Marder zone 0.125, npml 12
output: 3-D moments in blocks of 13 cells every 296 steps, slices z and y through the centre every 74 steps, history every 37 steps
output plan (uncompressed sizes; zstd typically saves 10-50 % more):
  diagnostic       kind         every    dt_out per period outputs      each     total  rule
  coarse           fields         296        74        0.1      11   23.3 MB    257 MB  296 steps
  slice_z          slice           74      18.5        0.3      41   47.2 MB   1.93 GB  74 steps
  slice_y          slice           74      18.5        0.3      41   47.2 MB   1.93 GB  74 steps
  spec_meridional  histogram      148        37        0.2      21   56.6 MB   1.19 GB  148 steps
  spec_equatorial  histogram      148        37        0.2      21   56.6 MB   1.19 GB  148 steps
  spec_all         histogram       74      18.5        0.3      41   3.58 kB    147 kB  74 steps
  spec_up          histogram       74      18.5        0.3      41   1.79 kB   73.5 kB  74 steps
  hot_e            particles      148        37        0.2      21   16.3 GB    342 GB  148 steps
  tracks_e         tracks           4         1        6.3     741   9.19 MB   6.81 GB  4 steps
  scalars every 37 steps; total ~356 GB
```

**M**

```text
=== pypic blast-wave planner (debris shell around a magnetised body) ===
ambient: m_i/m_e=16, B0=0.166667 along z, beta=1.44, v_te=0.1 c (Debye length 0.20 dx), v_A=0.0417 c, fast speed 0.0618 c
debris (on): M_A=6 (fast Mach 4.05), V_d=0.2500 c, R_m=44.8 d_e = 11.2 d_i, R_B=147.9 d_e, gyroradius 24 (0.536 R_m); shell at 14.78, width 4.48, peak density 16.5 in 4 beam(s) per species (omega_pe dt <= 1.05, 0.00 % inside the body); extent after Le et al. 2021, Eq. 4 with D = 3: 51.9 d_e
body: OFF (no sphere, no dipole)
grid (768, 768, 768), dx=0.5, L=384 d_e (half-width 1.30 R_B), dt=0.25, 2960 steps = 1.25 R_B/V_d = 7.7 / Omega_ci
1,811,939,328 ambient particles (2 per cell and species), about 25,607,222 of the debris (block of 3,112,136 cells), capacity factor 1.5; memory up to 55.9 GB/GPU on 8 GPU(s), rank grid (1, 2, 4); ids 6.0e-03 of 2^53
faces: six open faces with flux injection of the ambient species; filter 4+comp, Marder zone 0.125, npml 12
output: 3-D moments in blocks of 9 cells every 296 steps, slices z and y through the centre every 74 steps, history every 37 steps
output plan (uncompressed sizes; zstd typically saves 10-50 % more):
  diagnostic       kind         every    dt_out per period outputs      each     total  rule
  coarse           fields         296        74        0.1      11   68.7 MB    756 MB  296 steps
  slice_z          slice           74      18.5        0.3      41   47.2 MB   1.93 GB  74 steps
  slice_y          slice           74      18.5        0.3      41   47.2 MB   1.93 GB  74 steps
  spec_meridional  histogram      148        37        0.2      21   56.6 MB   1.19 GB  148 steps
  spec_equatorial  histogram      148        37        0.2      21   56.6 MB   1.19 GB  148 steps
  spec_all         histogram       74      18.5        0.3      41   3.58 kB    147 kB  74 steps
  spec_up          histogram       74      18.5        0.3      41   1.79 kB   73.5 kB  74 steps
  hot_e            particles      148        37        0.2      21   16.3 GB    342 GB  148 steps
  tracks_e         tracks           4         1        6.3     741   9.19 MB   6.81 GB  4 steps
  scalars every 37 steps; total ~356 GB
```

**H**

```text
=== pypic blast-wave planner (debris shell around a magnetised body) ===
ambient: m_i/m_e=16, B0=0.114943 along z, beta=3.03, v_te=0.1 c (Debye length 0.20 dx), v_A=0.0287 c, fast speed 0.0539 c
debris (on): M_A=8.7 (fast Mach 4.64), V_d=0.2500 c, R_m=34.8 d_e = 8.7 d_i, R_B=147.2 d_e, gyroradius 34.8 (1.000 R_m); shell at 11.48, width 3.48, peak density 16.5 in 4 beam(s) per species (omega_pe dt <= 1.05, 0.00 % inside the body); extent after Le et al. 2021, Eq. 4 with D = 3: 43.5 d_e
body: OFF (no sphere, no dipole)
grid (768, 768, 768), dx=0.5, L=384 d_e (half-width 1.30 R_B), dt=0.25, 2960 steps = 1.26 R_B/V_d = 5.3 / Omega_ci
1,811,939,328 ambient particles (2 per cell and species), about 12,002,231 of the debris (block of 1,481,544 cells), capacity factor 1.5; memory up to 55.3 GB/GPU on 8 GPU(s), rank grid (1, 2, 4); ids 6.0e-03 of 2^53
faces: six open faces with flux injection of the ambient species; filter 4+comp, Marder zone 0.125, npml 12
output: 3-D moments in blocks of 7 cells every 296 steps, slices z and y through the centre every 74 steps, history every 37 steps
output plan (uncompressed sizes; zstd typically saves 10-50 % more):
  diagnostic       kind         every    dt_out per period outputs      each     total  rule
  coarse           fields         296        74        0.1      11    144 MB   1.58 GB  296 steps
  slice_z          slice           74      18.5        0.3      41   47.2 MB   1.93 GB  74 steps
  slice_y          slice           74      18.5        0.3      41   47.2 MB   1.93 GB  74 steps
  spec_meridional  histogram      148        37        0.2      21   56.6 MB   1.19 GB  148 steps
  spec_equatorial  histogram      148        37        0.2      21   56.6 MB   1.19 GB  148 steps
  spec_all         histogram       74      18.5        0.3      41   3.58 kB    147 kB  74 steps
  spec_up          histogram       74      18.5        0.3      41   1.79 kB   73.5 kB  74 steps
  hot_e            particles      148        37        0.2      21   16.3 GB    342 GB  148 steps
  tracks_e         tracks           4         1        6.3     741   9.19 MB   6.81 GB  4 steps
  scalars every 37 steps; total ~357 GB
```
