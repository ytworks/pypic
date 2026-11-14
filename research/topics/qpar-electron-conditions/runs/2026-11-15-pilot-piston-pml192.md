# Pilot: the oblique piston of V1's plasma in a short box, absorbing layer of 192 cells

Kind: plan step 4 (engine check). Not a run of the plan's table; it tests no hypothesis. It repeats
`2026-11-15-pilot-piston-paired.md` with the planner's new depth of the absorbing layer of the x+ face
(192 cells instead of 24); nothing else differs.

| | |
|---|---|
| Date | 2026-11-15 (03:55 JST), run id `a53ee494` |
| Commit | 92686fb7f7740a9c2a2ee1fe4c72ad5f21762fe0 (clean tree; engine = main 70202dfab) |
| Command | `python run.py --config configs/pilot-piston.yaml --out <dir> --overwrite --progress 10000` (single process) |
| Config | `configs/pilot-piston.yaml`, no `--set`: M_A 20 (r4), mass ratio 100, theta_Bn 30 deg, v0 0.1 c, beta_e = beta_i = 0.5, one seed species (theta 0.02), `face: paired` (default), lx 40 d_i, t_end 2 / Omega_ci |
| Grid | 4000 x 2 x 2 cells, dx 0.1 d_e, dt 0.045, 66,667 steps, 32 filter passes, `compact`, npml 192, Marder zone 0.125 |
| Seed | 12345 |
| GPUs | 1 (RTX PRO 6000 Blackwell) |
| Output | `/mnt/nvme/pypic/qpar/pilot/piston-pml192`, 37.5 MB (not in git) |
| Analysis | `pilot_check.py <dir> --json` at commit 5bf6998 -> `2026-11-15-pilot-piston-pml192.json` |

## Numbers (from the JSON; in brackets the run with 24 cells)

- Wall time 864.4 s, 13.0 ms/step, 1.28 ns per particle and step (868.6 s).
- Gauss max 2.1e-14, div B max 5.6e-13 (2.3e-14, 5.7e-13).
- Particles 6.00 M -> 10.10 M; 4.79 M entered through the face (the same count, 4,787,763, as with 24 cells);
  largest fill of the storage 0.365.
- Shock (largest x with n_i > 2, mean over 0.5 d_i), per 0.2 / Omega_ci:
  0, 2.66, 5.50, 8.01, 10.23, 12.47, 14.52, 16.10, 17.89, 18.27, 19.85 d_i.
  Speed over the second half 0.0479 c; compression between 0.2 and 0.8 of the shock position 2.56
  (0.0457 c, 2.52). MHD jump of the planner: 0.0337 c and 3.96.
- Upstream at the end (22.9-37.0 d_i): n_e 1.027, v_x / v0 0.971, B_y / B0y 1.020, E_z / E0z 0.996
  (21.7-37.0 d_i: 1.046, 0.945, 1.043, 1.034). The moment temperatures there are 6.2 (electrons) and 27.7
  (ions) times the loaded ones.
- Histogram totals (weight in the box): electrons and ions 16.0 -> 28.0, seeds 16.0 -> 24.79.
  Share of the seeds above u0 + 5 sqrt(theta): 1.35e-5 -> 6.4e-3; above 8: 0 -> 3.9e-6.

## Limits

- 2 / Omega_ci: the shock is still forming (faster and less compressed than the MHD jump). Nothing about
  acceleration can be read from this run.
- The upstream region holds reflected particles, so its moment temperatures are not temperatures of the
  incoming plasma, and its B_y (2 % above the loaded value) is not the state that the face sets alone. That
  state is measured without a wall (notebook, 2026-11-15; `configs/pilot-through.yaml`).
- The two runs differ in the upstream zone that `pilot_check.py` takes (it follows the shock position), so the
  bracketed upstream numbers are not a like-for-like comparison of the two depths.
- One seed, one resolution, one particle number.
