# Pilot: flow-through control of V1's plasma to 40 / Omega_ci (no shock), second of two, with the velocity histogram

Kind: plan step 4 (engine check). Not a run of the plan's table; it tests no hypothesis. The companion of
`2026-11-15-pilot-through-pml192.md`: the same control with another seed and one more diagnostic, the histogram
of the electrons in (u_x, u_y, u_z) in the slab 20-30 d_i, which criterion 2 of the notebook (entry of
2026-11-15, "Pass criteria of the flow-through control") needs. The four criteria are applied to this run on its own.

| | |
|---|---|
| Date | 2026-11-15 (05:00 JST), run id `7c26e111` |
| Commit | d7b01dd50a8ae5c97921503de58083986b5a89cc (clean tree; engine = main 70202dfab) |
| Command | `python run.py --config configs/pilot-through-velocity.yaml --out <dir> --overwrite --progress 100000` (single process) |
| Config | `configs/pilot-through-velocity.yaml`: as `pilot-through.yaml` (M_A 20 (r4), mass ratio 100, theta_Bn 30 deg, v0 0.1 c, beta_e = beta_i = 0.5, `through: true`, lx 40 d_i, t_end 40 / Omega_ci, 125 particles per cell and species, 32 filter passes), plus `velocity_slab: [20, 30]`, 40 outputs |
| Grid | 4000 x 2 x 2 cells, dx 0.1 d_e, dt 0.045, 1,333,334 steps, `compact`, npml 192 |
| Seed | 23456 |
| GPUs | 1 (RTX PRO 6000 Blackwell) |
| Output | `/mnt/nvme/pypic/qpar/pilot/through-velocity`, 38.5 MB (not in git) |
| Analysis | `pilot_check.py <dir> --json` at commit 1f0d14a (clean tree) -> `2026-11-15-pilot-through-velocity.json`; `control_look.py <dir> --pair <through-pml192> --json` at commit 8287f8d (clean tree) -> `2026-11-15-pilot-through-velocity-look.json` |

## Numbers (from the two JSON files)

- Wall time 17,478.0 s, 13.11 ms/step, 3.28 ns per particle and step. Particles 4,000,013 -> 3,997,716;
  59,983,832 entered through the face; largest fill of the storage 0.769.
- Gauss max 2.7e-16 and div B max 6.9e-17 are the values of step 0 only (same defect as in the first control;
  fixed in commit 7d5598b after this run had started). The later values of this run were not measured.
- **Criterion 1, slope:** zone 30-37 d_i, t = 4-40, 180 outputs: slope of <B_y>/B_0y = +1.08e-5 +- 1.03e-5 per
  1 / Omega_ci; slope + 2 standard errors = 3.15e-5 against the limit 3.3e-4: **pass**.
  <B_y>/B_0y in the window: 0.99923 at its start, 1.00094 at its end, smallest 0.99773, largest 1.00442.
- **Criterion 1, halves:** t = 4-22: +3.09e-5 +- 2.69e-5; t = 22-40: -5.85e-5 +- 3.05e-5. Difference 8.95e-5
  against the tolerance 8.14e-5 (two standard errors of the difference): **fail**.
- **Criterion 2, incoming electrons** (velocity histogram, slab 20-30 d_i, frame of the flow, 36 outputs from
  t = 4, against the first output): T_par 0.9953-1.0056 (mean 1.0006), T_perp 0.9962-1.0034 (mean 0.9998); every
  output within 2 %: **pass**. Shares above 2, 3, 4 v_th,e, window mean against the first output, in units of
  the counting error: 0.26300 against 0.26267 (0.52 sigma), 0.029597 against 0.029560 (0.15 sigma), 0.0012455
  against 0.0011880 (1.16 sigma); all within 3 sigma: **pass**.
- **Criterion 3, drift** (slab 20-30 d_i, 180 outputs from t = 4), mean +- standard error in units of n_0 v_0:
  x -3.8e-6 +- 4.9e-5, y +1.5e-5 +- 3.5e-5, z +2.0e-5 +- 3.0e-5, against the limit 2e-3: **pass**.
- **Criterion 4 (reported, no threshold):** delta B_perp / B_0 in the slab = 0.00581 (B_y 0.00407, B_z 0.00413);
  first control 0.00584; N1 level 0.00469 (its three differences are in the first control's record).
- Whole box, not a criterion: moment temperatures against the loaded ones, 10 outputs from t = 4 to 40:
  electrons 0.9993-1.0010, ions 0.9984-1.0004. Slope of <B_y>/B_0y over 3-37 d_i from t = 1: +1.4e-6 +- 4.6e-6;
  second half +5.4e-6 +- 1.4e-5 (first control: -6.2e-5 +- 1.2e-5).
- Diagnostic look (`control_look.py`; not a criterion): rms of the residuals about the straight line 1.43e-3.
  Autocorrelation at lags of 1, 2, 5, 10, 20, 40 outputs: 0.626, 0.185, 0.073, -0.131, -0.062, 0.115.
  Means of <B_y>/B_0y - 1 over blocks of 6 / Omega_ci from t = 4: +8.6e-5, +5.6e-4, +5.2e-4, +1.23e-3, +3.3e-4,
  +6.3e-4. Largest spectral lines: shares 9.7 %, 7.6 %, 7.0 % at periods 11.9, 8.9, 4.0 / Omega_ci.
- The two controls against each other (`--pair`):
  - Residuals, output by output: correlation +0.098. Without a common part the standard error is 0.108
    (86 effective outputs of 180). The rms of the half sum is 1.06e-3 and of the half difference 0.96e-3; two
    independent series of this rms (1.43e-3) give 1.01e-3 for both.
  - Block means (six blocks of 6 / Omega_ci): correlation +0.757. In both runs the largest block is t = 22-28
    and the next one (t = 28-34) is lower by 8.3e-4 (first) and 9.1e-4 (second). The last block falls further in
    the first run (-2.0e-4) and rises in the second (+6.3e-4).

## Readings and limits

- The halves part of criterion 1 fails on this run too (difference 1.10 times the tolerance; first control 2.09
  times). Criteria 2 and 3 pass. By the fixed rule V1 does not start on these two results. No criterion was
  changed after the runs.
- In both runs the first half rises and the second falls. The readings of the first control's record stand, and
  the two runs add this:
  - Output by output (0.2 / Omega_ci) the two series share nothing measurable (0.098 +- 0.108): the fast part of
    the fluctuation depends on the seed. Its neighbouring outputs correlate (0.56 and 0.63), so the fit errors of
    `pilot_check.py`, which treat outputs as independent, are too small by about 1.9 and 2.1 (AR(1) form assumed).
    With these factors the tolerances would be 1.49e-4 and 1.70e-4 against the differences 1.65e-4 and 8.95e-5:
    the first run stays outside, the second would be inside.
  - On the scale of 6 / Omega_ci the two runs resemble each other (0.757 from six block means; for six
    independent pairs a value this large or larger in magnitude has a probability of about 0.08). A common slow
    part of the size 5e-4 to 1e-3, caused by the boundaries or the start, is neither shown nor excluded by two
    runs of 40 / Omega_ci.
- What neither run answers: the field's change between t = 40 and 150 (V1's length), and whether the slow part
  grows with the box length or the particle number.
- Size, as a fact and not as a criterion: every output of both runs lies within 0.9955-1.0051, and every block
  mean within -6.3e-4 to +1.23e-3, of the loaded B_0y. The tolerance of V1's Q monitor is 5 %.
- Limits: two seeds; 125 particles per cell (V1: 250); a box of 40 d_i (V1: 2700 d_i); 40 of V1's 150 / Omega_ci;
  no shock.

## Correction (added later, 2026-11-15)

The numbers above stand. Two statements are corrected by `notebook.md` (2026-11-15, "controls to 150", facts
2 and 3):

- "too small by about 1.9 and 2.1 (AR(1) form assumed)": the AR(1) form is contradicted by the lag-2
  correlations (AR(1) gives 0.31 and 0.39; measured 0.137 and 0.185). The outputs share plasma while the flow
  crosses the zone (0.467 / Omega_ci): correlation 0.571, 0.143, 0 at lags 1, 2, >= 3, and the fit errors are
  1.56 times those of independent outputs in both runs. With 1.56 the differences of the halves are 2.7 and
  1.4 standard errors. The result by the fixed rule is not changed by this.
- "Largest spectral lines": the periods 11.9, 8.9 and 4.0 / Omega_ci are the harmonics 3, 4 and 9 of the
  36 / Omega_ci window; they are not lines.
