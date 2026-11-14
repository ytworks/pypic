# Pilot: the oblique piston of V1's plasma in a short box, paired face, one seed species

Kind: plan step 4 (engine check). Not a run of the plan's table; it tests no hypothesis.

| | |
|---|---|
| Date | 2026-11-15 (02:07 JST), run id `10c436a9` |
| Commit | ba98d45bd8aa22589a23003ade7d2807ffd6501b (clean tree; engine = main 70202dfab) |
| Command | `python run.py --config configs/pilot-piston.yaml --out <dir> --overwrite` (single process) |
| Config | `configs/pilot-piston.yaml`, no `--set`: M_A 20 (r4), mass ratio 100, theta_Bn 30 deg, v0 0.1 c, beta_e = beta_i = 0.5, one seed species (theta 0.02), `face: paired` (default), lx 40 d_i, t_end 2 / Omega_ci |
| Grid | 4000 x 2 x 2 cells, dx 0.1 d_e, dt 0.045, 66,667 steps, 32 filter passes, `compact`, npml 24, Marder zone 0.125 |
| Seed | 12345 |
| GPUs | 1 (RTX PRO 6000 Blackwell) |
| Output | `/mnt/nvme/pypic/qpar/pilot/piston-paired`, 37.5 MB (not in git) |
| Analysis | `pilot_check.py <dir> --json` -> `2026-11-15-pilot-piston-paired.json` |

## Numbers (from the JSON)

- Wall time 868.6 s, 13.0 ms/step, 1.29 ns per particle and step.
- Gauss max 2.3e-14, div B max 5.7e-13.
- Particles 6.00 M -> 10.12 M; 4.79 M entered through the face; largest fill of the storage 0.365.
- Shock (largest x with n_i > 2, mean over 0.5 d_i), per 0.2 / Omega_ci:
  0, 2.66, 5.51, 7.88, 10.43, 12.49, 14.48, 15.85, 17.12, 18.80, 18.67 d_i.
  Speed over the second half 0.0457 c; compression between 0.2 and 0.8 of the shock position 2.52.
  MHD jump of the planner: 0.0337 c and 3.96.
- Upstream at the end (21.7-37.0 d_i): n_e 1.046, v_x / v0 0.945, B_y / B0y 1.043, E_z / E0z 1.034.
  The moment temperatures there are 6.5 (electrons) and 42.8 (ions) times the loaded ones.
- Histogram totals (weight in the box): electrons and ions 16.0 -> 28.0, seeds 16.0 -> 24.96.
  Share of the seeds above u0 + 5 sqrt(theta): 1.35e-5 -> 6.7e-3; above 8: 0 -> 6.7e-6.

## Limits

- 2 / Omega_ci: the shock is still forming (faster and less compressed than the MHD jump). Nothing about
  acceleration can be read from this run.
- The upstream region holds reflected particles, so its moment temperatures are not temperatures of the
  incoming plasma.
- B_y and E_z in the upstream region are 4 % and 3 % above the loaded values at the end. A rise of B_y at the
  paired face is also present without a wall (notebook, 2026-11-15); its size at later times was not known when
  this run was made.
- One seed, one resolution, one particle number.
