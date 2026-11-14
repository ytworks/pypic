# Pilot: the particle dumps of N12 and N6 (a) on a real run (V1's plasma, box of 60 d_i, to 6 / Omega_ci)

Kind: plan step 4 (engine check). Not a run of the plan's table; it tests no hypothesis. Its purpose: the mechanics
of V1's dumps (`origin_e`, `wall_e`, `edge_e`, `edge_i`, `sample_e*`) and of `n12_n6.py` on dumps the engine wrote.
The box is 45 times shorter than V1's and the run 25 times shorter, so every length and time below is scaled and no
verdict below is V1's. Nothing here is a statement about acceleration.

| | |
|---|---|
| Date | 2026-11-15 (07:02 JST), run id `d02a11ee` |
| Commit | d74ba7c03cc320ca372c15220bd1fedc48e16f65 (clean tree; engine = main 70202dfab) |
| Command | `python run.py --config configs/pilot-dumps.yaml --out <dir> --overwrite --progress 25000` (single process) |
| Config | `configs/pilot-dumps.yaml`: M_A 20 (r4), mass ratio 100, theta_Bn 30 deg, v0 0.1 c, beta_e = beta_i = 0.5, `face: paired`; lx 60 d_i, t_end 6 / Omega_ci; dumps: origin (x 0-30 d_i, t >= 4, at the spectra outputs), wall (x <= 5 d_i, every 0.6), edge (zone of 10 d_i, every 0.75), sample (every 101st id, x 5-55 d_i, t = 1, 2, 4) |
| Grid | 6000 x 2 x 2 cells, dx 0.1 d_e, dt 0.045, 200,000 steps, ppc 125, `compact`, npml 192 |
| Seed | 12345 |
| GPUs | 1 (RTX PRO 6000 Blackwell) |
| Output | `/mnt/nvme/pypic/qpar/pilot/dumps`, 89.4 MB (not in git): `particles/` 33.2 MB, `reduced/` 56.0 MB |
| Analysis | at commit 0c62c11 (clean tree): `pilot_check.py <dir> --json` -> `2026-11-15-pilot-dumps.json`; `n12_n6.py <dir> --smooth 1 --region 0.75 0.25 --clear 0.1 --blocks 4 6 0.4 --xc 2.5 2 3 --slab 1 10 --edge 3 4.5 6 --margin 5 --json` -> `2026-11-15-pilot-dumps-n12n6.json` |

A first attempt from commit 43dfa1c was stopped by the session's time limit near step 125,000-140,000. It is not
used (`/mnt/nvme/pypic/qpar/pilot/dumps-stopped`, not a record).

## Numbers

Run (`2026-11-15-pilot-dumps.json`, the log, the scalars):

- Wall time 2648.8 s, 13.2 ms/step, 0.88 ns per particle and step. dE/E0 at the end -1.58e-4.
- Particles 6.00 M -> 14.98 M; 8.99 M entered through the face; largest fill of the storage 0.520.
- Gauss and div B, 100 samples (every 2000 steps): Gauss max 6.2e-14. div B max 2.0e-12 (step 190,000); it
  grows with the step count: 5.2e-14, 1.9e-13, 5.5e-13, 7.2e-13, 9.6e-13, 8.9e-13, 1.3e-12, 1.4e-12, 1.8e-12 at
  steps 20,000 to 180,000 in steps of 20,000; first above 1e-12 at step 124,000.
- Shock (largest x with n_i > 2), every 0.6 / Omega_ci from t = 0:
  0, 7.89, 14.92, 18.69, 19.53, 20.07, 21.54, 22.7, 24.1, 22.14, 24.12 d_i; compression 5.45.

Dumps written (electrons or ions per dump: first / middle / last; total; bytes per particle):

| Dump | Dumps | Particles per dump | Total | Bytes | Bytes per particle |
|---|---|---|---|---|---|
| `origin_e` (steps 136,000-200,000, every 4000) | 17 | 42,263 / 48,278 / 60,534 | 852,501 | 24,675,775 | 28.9 |
| `wall_e` (steps 0-200,000, every 20,000) | 11 | 0 / 167 / 32,619 | 127,779 | 4,045,800 | 31.7 |
| `edge_e` (every 25,000) | 9 | 0 / 20 / 357 | 565 | 292,696 | (about 30 kB per file) |
| `edge_i` (every 25,000) | 9 | 0 / 9 / 17,466 | 20,253 | 817,295 | 40.4 |
| `sample_e0`, `_e1`, `_e2` (steps 33,333, 66,667, 133,333) | 1 each | 29,613, 37,449, 48,721 | 115,783 | 3,401,197 | 29.4 |

A dump file is written also when no particle is selected.

`n12_n6.py` (`2026-11-15-pilot-dumps-n12n6.json`). Scaled arguments: region x_sh - 15 to x_sh - 5 d_i (D_fix);
X_c = 50 (40, 60) d_i in a box of 60 d_i; wall zone x <= 5 d_i; slab 48.08-57.08 d_i, dump zone 48.08-60 d_i;
dump interval over slab crossing time 12.5 (V1: 0.577).

- Dump against histogram (README check, criterion 1 %): electrons above p_d,e in D_fix, whole bins at and above
  p_d,e, 17 outputs: largest relative difference 1.8e-8. Pass. Beside: the dump holds 0.962-0.987 of the count
  that the histogram gives when the bin cut by p_d,e is split by its logarithmic fraction.
- Samples (every 101st electron id, no threshold; counts in macro-particles):

  | t | x_sh | In D_fix | Below p_d,e and left out by (L) | Share | p_th,e | p_th,e ratio under (L) |
  |---|---|---|---|---|---|---|
  | 1.0 | 12.50 | 2520 | 2520 | 1.000 | - | - |
  | 2.0 | 19.30 | 11,835 | 11,577 | 0.978 | 0.03025 | 0.866 (changes) |
  | 4.0 | 22.15 | 19,930 | 9572 | 0.481 | 0.02572 | 0.984 |

  Ids in the same samples (not judged): loaded electrons upstream of x_sh + 5 d_i against x_0 - v0 t: 13,661
  electrons, residual mean 0.92, rms 5.2 d_i (t = 1); 3725, 14.6, 17.1 d_i (t = 2); 2151, 71.6, 18.9 d_i (t = 4).
  Injected electrons: 4943, 12,410, 27,472, smallest x 31.98, 11.93, 5.01 d_i against the ballistic front of the
  inflow at 45.0, 30.0, 0.0 d_i.
- N12 blocks: five blocks of 0.4 / Omega_ci over t = 4-6 with 3, 3, 4, 3, 4 outputs; no dump missing; no output
  with the region near the wall or outside the dump. Ids in the wall dumps up to each block: 46,945, 46,945,
  55,781, 62,536, 74,437. Bins that turn negative after the subtraction: at most 4 per variant, each above
  -1e-5 macro-particles.
- D_fix, electrons above p_0,e per block (macro-particles), variants (2) none / (1) L / (3) L at 40 / (3) L at 60 /
  (4) E / (5) L and E:
  752.5 / 5.7 / 24.5 / 0.3 / 338.5 / 5.7; 1108.8 / 13.1 / 48.2 / 0.9 / 473.0 / 12.9;
  2657.0 / 14.1 / 65.3 / 2.9 / 553.0 / 11.1; 1695.9 / 12.2 / 43.3 / 2.1 / 268.6 / 8.6;
  3831.3 / 70.4 / 202.3 / 29.5 / 1060.9 / 62.8.
- Last block, D_fix: 11.78 M macro-particles in 4 outputs, p_th,e 0.02942, p_0,e 0.1471 (above p_d,e = 0.1).
  eta_nt,e: (2) 3.25e-4 with sigma_pt 3.9e-5; (1) 5.97e-6; (3) 1.72e-5 and 2.50e-6; (4) 9.00e-5; (5) 5.33e-6.
  Rule: "wall effect detected", range 0.992 of the value (above the 30 % limit). phi_e is 0 in all variants and
  no slope q is returned (q not judged). eta_e above the fixed momentum: (2) 5.2e-8, the others 0, sigma_pt 1.4e-7.
  Region from D_fix to the shock: eta_nt,e (2) 3.89e-4, (1) 9.9e-6, (4) 1.35e-4, (5) 9.3e-6; same verdict.
- N6 (a), sums of 9 dumps per species, intervals regular:

  | t | Species | p_0 | E_nt | N+ | N | U_far | Verdict |
  |---|---|---|---|---|---|---|---|
  | 3.0 | e | 0.172 | 2.08 | 0 | 0 | 0 | pass |
  | 3.0 | i | 3.25 | 0 | - | - | - | none (nothing above p_0 in the box) |
  | 4.44 | e | 0.135 | 3727.5 | 7.45 | 7.45 | 2.09 | pass (R 0.0020) |
  | 4.44 | i | 2.12 | 6570.5 | 465.3 | 465.3 | 1673.0 | not decided (R 0.071, lower -0.18) |
  | 6.0 | e | 0.144 | 3627.7 | 185.3 | 183.9 | 30.6 | not decided (R upper 0.0511, R 0.0507, lower 0.042) |
  | 6.0 | i | 1.97 | 145,084 | 11,757.7 | 11,757.7 | 6515.1 | not decided (R 0.081, lower 0.036) |

  Energies in macro-particle units of m_e c^2. Output: `pass` false, first time that is not a pass 3.0. The engine's
  cumulative kinetic energy of removed particles at the three times: 1962, 8887, 21,436.

## V1's dump sizes from these numbers (arithmetic, not a measurement)

29-32 bytes per dumped particle; about 30 kB for an empty file. V1 has 250 particles per species and cell and
writes 17 `origin_e`, 251 `wall_e`, 201 `edge_e`, 201 `edge_i` and 3 sample dumps.

- If every particle of the zone were above the threshold (bound): `origin_e` about 9 GB per dump (0-760 d_i
  downstream at compression 4), `wall_e` about 0.6 GB per dump, each edge dump 0.6-0.8 GB (zone of 202 d_i at the
  loaded density): about 150 + 160 + 240 to 330 GB in total. Free space on `/mnt/nvme`: 2.5 TB.
- With the shares of this pilot's last dumps (1.0 % of the electrons of the origin zone, 3.9 % of the wall zone,
  from the histogram at step 200,000): `origin_e` about 90 MB per dump (1.6 GB in total), `wall_e` about 25 MB
  per dump (6 GB in total). V1's shares at t = 100-150 are not known.
- Samples: at most about 90 MB each (x = 150-900 d_i all at compression 4).

## Readings and limits

- What the run checks: the engine writes the five kinds of dumps at the planned steps, with ids; `n12_n6.py` reads
  them and reproduces the engine's histogram counts above p_d,e to 1.8e-8; the id split
  into loaded and injected electrons and all exclusions run; the N6 (a) branches "pass", "not decided" and "none"
  occur on real dumps ("fail" did not occur; it is covered only by the synthetic check).
- The N12 verdict "wall effect detected" is the verdict of this box, where D_fix lies 7-19 d_i from the wall and
  X_c = 50 d_i covers most of the box: 98 % of the non-thermal electrons of D_fix were loaded at x_0 <= X_c and
  72 % have their id in a wall dump. It shows that the rule returns "detected" in a case where the wall's start
  supplies the non-thermal electrons. It says nothing about V1.
- Sample check, "expected zero": downstream plasma at x holds material loaded near x_0 = r x (r = compression), so
  in this box D_fix holds material loaded at x_0 <= 50 d_i (share 1.0, 0.98, 0.48). In V1 the material of D_fix at
  t = 100-150 comes from x_0 of about 800 d_i or more against X_c = 400 d_i (arithmetic with r = 4; beyond
  2700 d_i it has entered through the face). The pilot does not test that expectation; it shows that the check
  reports the share and the change of p_th,e.
- x_0 from the id: at t = 1 the residual (mean 0.92, rms 5.2 d_i) is of the size of the thermal spread along B
  (sqrt(theta_e) cos 30 deg t = 4.3 d_i per 1 / Omega_ci) and the selection is cut at x_sh + 5 d_i; it is not a
  sharp test (the exact arithmetic of the id is tested in `2026-11-15-n12-n6-synthetic`). At t = 2 and 4 the loaded
  electrons still upstream are those that did not move with the flow, so the residual there does not test the id.
  Injected electrons ahead of the ballistic front: by 13 d_i at t = 1, which is three thermal spreads.
- N6 (a): with 12.5 slab crossings per dump interval the sums are not estimates of the escaped energy; only the
  bookkeeping is exercised. The ions reach the face zone of this box by t = 4.4-6. A species with nothing above
  p_0 in the box gets no verdict, and the code then counts that time as "not a pass" (t = 3.0 here).
- div B: the engine's gate criterion (1e-12) is passed in the gate cases; this run goes above it from step 124,000.
  Runs of the same plasma: without a shock (`n1-full`, periodic, 5,000,000 steps) 7.3e-13 at the end, growing as
  the square root of the step count between 1.25 M and 4.95 M steps; with a shock 5.6e-13 at step 65,934
  (`piston-pml192`). Two readings, not separated: round-off of the field update that accumulates faster where the
  fields are larger (behind the shock |B| is several B_0); or a contribution of the wall, the face or the layer.
  Where the maximum sits was not measured. V1 has 25 times as many steps; its value is not known.
- Gauss and div B are sampled only at steps that are multiples of both `diag_interval` and `gauss_interval`. With
  `plan.py` at this commit that is every 2000 steps here, but only step 0 in runs whose step count is not round
  (`piston-t10`: its record's "Gauss max, div B max" are the values of step 0).
