# Pilot: the oblique piston of V1's plasma to 10 / Omega_ci in a box of 80 d_i

Kind: plan step 4 (engine check). Not a run of the plan's table; it tests no hypothesis. Its purpose: `spectra.py`
(V1's analysis) on a formed shock, and the response of the Q monitor and the edge criterion when the precursor
reaches the face. The box is 34 times shorter than V1's, so every window, region and zone below is scaled and none
of the verdicts is V1's.

| | |
|---|---|
| Date | 2026-11-15 (05:06 JST), run id `ab3ae170` |
| Commit | d7b01dd50a8ae5c97921503de58083986b5a89cc (clean tree; engine = main 70202dfab) |
| Command | `python run.py --config configs/pilot-piston.yaml --set lx=80 --set t_end=10 --out <dir> --overwrite` (single process) |
| Config | `configs/pilot-piston.yaml`: M_A 20 (r4), mass ratio 100, theta_Bn 30 deg, v0 0.1 c, beta_e = beta_i = 0.5, one seed species, `face: paired`; lx 80 d_i, t_end 10 / Omega_ci |
| Grid | 8000 x 2 x 2 cells, dx 0.1 d_e, dt 0.045, 333,334 steps, `compact`, npml 192 |
| Seed | 12345 |
| GPUs | 1 (RTX PRO 6000 Blackwell) |
| Output | `/mnt/nvme/pypic/qpar/pilot/piston-t10`, 74.7 MB (not in git) |
| Analysis | at commit 07b3bfd (clean tree): `pilot_check.py <dir> --json` -> `2026-11-15-pilot-piston-t10.json`; `spectra.py <dir> --floor <n1-full> --smooth 1 --early 5 10 --blocks 5 10 1 --region 1.0 0.25 --clear 0.1 --q-zone 1 10 --q-block 1 --edge 5 7.5 10 --edge-zone 10 --json` -> `2026-11-15-pilot-piston-t10-spectra.json` |

## Numbers (from the two JSON files)

- Wall time 5685.9 s, 17.1 ms/step, 0.55 ns per particle and step. Gauss max 3.9e-16, div B max 6.9e-17.
- Correction (added later on 2026-11-15): the two maxima above are the values of step 0. With 333,334 steps the
  run's two diagnostic intervals (333 and 3333) have no common multiple inside the run, and the engine evaluates
  Gauss and div B only at such steps. The later values of this run were not measured (notebook, entry of
  2026-11-15 on the Gauss / div B samples).
- Particles 12.00 M -> 30.74 M; 23.77 M entered through the face; largest fill of the storage 0.590.
- Shock (largest x with n_i > 2, mean over 0.5 d_i), per 1 / Omega_ci:
  0, 12.77, 19.95, 19.52, 23.28, 26.41, 22.44, 29.01, 31.6, 30.43, 33.85 d_i.
- t = 5-10 (`spectra.py`, mean over 1 d_i): shock speed 0.0156 +- 0.0008 c against v0 / 3 = 0.0333 c
  (ratio 0.47; the MHD jump gives 0.0337 c); compression 5.45 (4.73-5.94), MHD jump 3.96. Both scaled checks of the early window fail.
- Upstream at the end (36.9-77.0 d_i): n_e 1.26, v_x / v0 0.68, B_y / B0y 0.97, E_z / E0z 0.93; moment
  temperatures 8.6 (electrons) and 220 (ions) times the loaded ones.
- Blocks of 1 / Omega_ci over t = 5-10, region x_sh - 20 to x_sh - 5 d_i:
  electrons p_th 0.0296-0.0307 m_i v_pt, eta_nt 1.8e-4 to 3.3e-4, phi_e = 0 and no slope is returned;
  ions p_th 0.395-0.430 m_i v_pt, eta_nt 1.4e-4 to 1.8e-3.
  sigma_pt of the last block: p_th 1.5 % (electrons), 3.1 % (ions); eta_nt 0.25 (electrons), 4.8 (ions) of its value.
- Q monitor (zone 68.1-77.1 d_i), t = 0 to 9: 1.001, 1.003, 1.002, 1.009, 1.017, 0.995, 0.926, 0.725, 0.835;
  largest deviation 0.275: fail.
- Edge (zone 68.1-78.1 d_i): delta B_perp / B_0 0.100, 0.301, 0.191 at t = 5, 7.5, 10, which is 23, 66 and 42 times
  the level of the run without a shock (N1, `n1-full`): fail.
- Histogram totals (weight in the box): electrons and ions 32.0 -> 90.67, seeds 32.0 -> 64.55.

## Readings and limits

- The shock is slower and more compressed than the MHD jump, and the upstream state differs from the loaded one
  by 4-32 % (n_e, v_x, B_y, E_z). In this box the shock is 34-56 d_i from the face zone, and reflected particles and their waves fill
  the whole upstream within a few 1 / Omega_ci. Two readings are open and this run does not separate them: the
  shock runs into a flow that its own precursor has already slowed and heated; or the face does not hold the inflow once the precursor reaches it. The Q monitor and the edge criterion
  exist to flag such a state: the edge criterion fails at all three times, and Q stays within 2 % of 1 to t = 6
  and is 0.93, 0.72, 0.84 afterwards.
- V1's box is 2700 d_i; whether its face zone stays quiet to t = 150 is not answered here. It is read from V1's
  own Q monitor and edge criterion.
- The non-thermal fractions are those of 5 / Omega_ci of a disturbed shock, with a region of 15 d_i and blocks of
  1 / Omega_ci; sigma_pt of eta_nt of the ions is 4.8 times its last value. Nothing about acceleration can be read from this run.
- `spectra.py` ran through all its parts on a formed shock (shock position, blocks, sigma_pt, sensitivity
  regions, Q monitor, edge). No slope was returned (phi_e = 0), so the slope fit is not exercised here.
- One seed, one resolution, one particle number.
