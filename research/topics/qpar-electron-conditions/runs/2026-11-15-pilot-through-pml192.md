# Pilot: flow-through control of V1's plasma to 40 / Omega_ci (no shock), first of two

Kind: plan step 4 (engine check). Not a run of the plan's table; it tests no hypothesis. Its purpose: V1's plasma
flows through a box with V1's two boundaries (open x-, paired face at x+ with the absorbing layer) and no shock,
and the pass criteria fixed in the notebook before the run (entry of 2026-11-15, "Pass criteria of the
flow-through control") are applied to it. The companion run (`through-velocity`, another seed, with the
velocity histogram that criterion 2 needs) has its own record.

| | |
|---|---|
| Date | 2026-11-15 (04:08 JST), run id `b2d6ea61` |
| Commit | 5bf6998bfd8ca2d44dfc69f6fea888182defea01 (clean tree; engine = main 70202dfab) |
| Command | `python run.py --config configs/pilot-through.yaml --out <dir> --overwrite --progress 100000` (single process) |
| Config | `configs/pilot-through.yaml`: M_A 20 (r4), mass ratio 100, theta_Bn 30 deg, v0 0.1 c, beta_e = beta_i = 0.5, `through: true`, lx 40 d_i, t_end 40 / Omega_ci, 125 particles per cell and species, 32 filter passes |
| Grid | 4000 x 2 x 2 cells, dx 0.1 d_e, dt 0.045, 1,333,334 steps, `compact`, npml 192 |
| Seed | 12345 |
| GPUs | 1 (RTX PRO 6000 Blackwell) |
| Output | `/mnt/nvme/pypic/qpar/pilot/through-pml192`, 35.9 MB (not in git) |
| Analysis | `pilot_check.py <dir> --json` at commit 7d5598b (clean tree) -> `2026-11-15-pilot-through-pml192.json`; `control_look.py <dir> --json` at commit d5ad433 (clean tree) -> `2026-11-15-pilot-through-pml192-look.json` |

## Numbers (from the two JSON files)

- Wall time 17,369.6 s, 13.03 ms/step, 3.26 ns per particle and step. Particles 3,999,994 -> 3,997,985;
  59,975,682 entered through the face; largest fill of the storage 0.769.
- Gauss max 3.3e-16 and div B max 6.9e-17 are the values of step 0 only: with 1,333,334 steps the run's two
  diagnostic intervals had no common multiple inside the run (defect fixed in commit 7d5598b, after this run
  had started). The later values of this run were not measured.
- **Criterion 1, slope:** zone 30-37 d_i, t = 4-40, 180 outputs: slope of <B_y>/B_0y = +4.43e-6 +- 1.03e-5 per
  1 / Omega_ci; slope + 2 standard errors = 2.50e-5 against the limit 3.3e-4: **pass**.
  <B_y>/B_0y in the window: 0.99958 at its start, 1.00028 at its end, smallest 0.99554, largest 1.00510.
- **Criterion 1, halves:** t = 4-22: +9.29e-5 +- 2.19e-5; t = 22-40: -7.25e-5 +- 3.28e-5. Difference 1.65e-4
  against the tolerance 7.89e-5 (two standard errors of the difference): **fail**.
  Each half on its own: |slope| + 2 standard errors = 1.37e-4 and 1.38e-4 (below 3.3e-4; not a criterion).
- **Criterion 2:** needs the velocity histogram of the companion run; not applicable to this run. Moment
  temperatures of the whole box against the loaded ones, 10 outputs from t = 4 to 40: electrons 0.9984-0.9992,
  ions 0.9985-1.0003. No weight above 5 or 8 thermal momenta in the (x, |u|) histograms at the end.
- **Criterion 3, drift** (slab 20-30 d_i, 180 outputs from t = 4), mean +- standard error in units of n_0 v_0:
  x -7.9e-7 +- 3.7e-5, y +5.3e-6 +- 3.1e-5, z +8.7e-6 +- 3.1e-5, against the limit 2e-3: **pass**.
- **Criterion 4 (reported, no threshold):** delta B_perp / B_0 in the slab = 0.00584 (B_y 0.00408, B_z 0.00416).
  N1 level: 0.00469 (B_y 0.00329, B_z 0.00333) from `n1-full` with `pilot_check.py --slab 5 15` (195 outputs from
  t = 4 to 150). The two differ in three ways: N1 has 250 particles per cell and species (here 125), its slab is
  5-15 d_i of a periodic box of 20 d_i, and it runs to t = 150. `n1-full` is a diagnostic pilot of the notebook,
  not a clean-tree record.
- Whole box (3-37 d_i), not a criterion: slope of <B_y>/B_0y from t = 1: +2.1e-7 +- 4.6e-6; second half
  -6.2e-5 +- 1.2e-5. <E_z>/E_0z from t = 1: -4.2e-6 +- 4.3e-5 (rms about the fit 6.7e-3).
- Diagnostic look at the zone's series (`control_look.py`; written after the numbers above had been read; not a
  criterion): rms of the residuals about the straight line 1.43e-3. Autocorrelation of the residuals at lags of
  1, 2, 5, 10, 20, 40 outputs (0.2 to 8 / Omega_ci): 0.560, 0.137, 0.240, 0.134, -0.163, 0.088.
  Means of <B_y>/B_0y - 1 over blocks of 6 / Omega_ci from t = 4: -6.3e-4, +5.4e-4, +6.3e-4, +8.5e-4, +2.7e-5,
  -2.0e-4. Largest spectral lines of the residuals: shares 5.4 %, 4.7 %, 4.2 % at periods 35.8, 7.2 and
  11.9 / Omega_ci.

## Readings and limits

- By the rule fixed before the run ("V1 starts only if all of the following pass"), the halves part of
  criterion 1 fails on this run, so V1 does not start on this result. No criterion was changed after the run.
- Why the halves disagree is open. Two readings, which this one run does not separate:
  (a) The outputs are not independent. The fit errors treat them as independent (a detail of `pilot_check.py`
  fixed before the run), and the residuals of neighbouring outputs correlate with 0.56. An AR(1) process with
  that value has fit errors 1.88 times larger; with that factor the tolerance would be 1.49e-4 against the
  difference 1.65e-4, which is still outside. The AR(1) form is an assumption: the correlation at lag 5 (0.24)
  is larger than 0.56^5 = 0.06.
  (b) The rate changes with time: the block means rise from t = 4 to 28 and fall afterwards. In that case
  36 / Omega_ci do not bound the field's change up to t = 150.
- The companion run has another seed. A rise and fall at the same times there supports (b); a different pattern
  supports a fluctuation that depends on the seed. Two runs do not settle it.
- The slope over the whole window passes with a margin of 13 (limit over slope + 2 standard errors). With the
  errors enlarged by 1.88 it would be 4.43e-6 + 2 x 1.94e-5 = 4.3e-5, still below the limit. This is not a
  criterion of the notebook.
- Limits: one seed; 125 particles per cell (V1: 250); a box of 40 d_i (V1: 2700 d_i); 40 of V1's 150 / Omega_ci;
  no shock, so nothing reaches the face from downstream.

## Correction (added later, 2026-11-15)

The numbers above stand. Two statements are corrected by `notebook.md` (2026-11-15, "controls to 150", facts
2 and 3):

- Reading (a), "AR(1), fit errors 1.88 times larger": the AR(1) form is contradicted by the lag-2 correlation
  (AR(1) gives 0.31; measured 0.137). The outputs share plasma while the flow crosses the zone (0.467 /
  Omega_ci): correlation 0.571, 0.143, 0 at lags 1, 2, >= 3, and the fit errors are 1.56 times those of
  independent outputs. With 1.56 the difference of the halves is 2.7 standard errors: the fail stands.
- "Largest spectral lines": the periods 35.8, 7.2 and 11.9 / Omega_ci are the harmonics 1, 5 and 3 of the
  36 / Omega_ci window; they are not lines.
