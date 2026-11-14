# Notebook — qpar-electron-conditions

Dated lab notes, newest first. Write down what was tried, what was seen, and what was decided, failed attempts included.
Never rewrite a past entry. To correct one, add a new entry that refers to it.

## 2026-11-15 — Step 4, second part: pilots on the merged engine; paired face, absorbing layer of 192 cells

Step 4 of the plan, after the engine change. Unless noted, the runs of this entry are scratch runs (diagnostic,
not recorded in `runs/`): 1 GPU, V1's plasma (`configs/pilot-piston.yaml`: M_A 20, mass ratio 100, θ_Bn 30°,
v0 0.1 c, β_e = β_i = 0.5), 4000×2×2 cells (40 d_i), dx 0.1 d_e, dt 0.045, 32 filter passes, `compact`,
`npml = 24` (unless the run's line gives another depth), seed 12345, engine tree 68b4926 (`pypic/` identical
to main 70202df). Outputs:
`/mnt/nvme/pypic/qpar/pilot/<name>`. Times are in 1/Ω_ci (15,000 ω_pe⁻¹ = 33,333 steps), fields in units of
the loaded B0y and E0z = v0 B0y.

### Engine
- PR #52 merged on the owner's instruction: main 70202df (`external_fields` with `B0` / `E0`, the refusal of a
  `B0` normal to a wall, injected test species). Taken into trunk (0df1f2e) and into this topic (d5b2162) with
  `research-sync-main`. The topic code now runs on the topic branch itself.

### N1 in a short box: the grid does not heat this plasma in 5 / Ω_ci when the current is filtered
`configs/pilot-n1.yaml`: the same plasma, periodic, 2000×2×2 cells, 2000 particles per column, 166,667 steps.
T is the moment temperature of `pilot_check.py`, averaged over the box, relative to the first output.

| run | change | ms/step | Gauss | div B | T_e end | T_i end |
|---|---|---|---|---|---|---|
| n1-f32 | as configured (32 passes) | 1.55 | 2.7e-15 | 1.0e-13 | 0.99999 | 0.99986 |
| n1-f4 | 4 passes | 1.49 | 6.0e-15 | 1.1e-13 | 0.99982 | 0.99986 |
| n1-f0 | no filter | 1.47 | 9.8e-14 | 1.1e-13 | 1.02563 | 1.00148 |
| n1-f32c | 32 passes + compensator | 1.56 | 3.8e-15 | 1.1e-13 | 0.99990 | 0.99997 |
| n1-f32-fp64 | precision `fp64` | 1.90 | 7.1e-15 | 1.2e-13 | 0.99991 | 0.99990 |
| n1-f32-ppc125 | 1000 per column | 0.91 | 2.7e-15 | 1.2e-13 | 0.99995 | 0.99991 |

- Without the filter T_e is 1.2 % up at half time and 2.6 % up at the end (still rising). With 4 or 32 passes
  every temperature stays within 2e-4 of the loaded one.
- Share of the electron histogram above u0 + 5 √θ at the end: ≤ 5e-7 in every run (one histogram count); above
  8 √θ: 0.
- Limits: 5 / Ω_ci of V1's 150; one seed; a periodic box (no face, no wall, no shock). The full length is below.

### N1 at full length (150 / Ω_ci)
`n1-full`: `configs/pilot-n1.yaml` with `--set t_end=150` (2000×2×2 cells, periodic, 4.0 M particles, 1 GPU):
5,000,000 steps in 7699.4 s (1.54 ms/step, 0.385 ns per particle and step), output 18 MB. A diagnostic pilot
like the others of this entry, not a clean-tree record. `pilot_check` summary: scratch `n1-full.json`.

- Gauss max 1.8e-14, div B max 7.5e-13, dE/E0 −2.6e-7.
- T_e / first output at t = 0, 15, …, 150: 1, 0.99995, 0.99991, 0.99986, 0.99990, 0.99994, 0.99992, 0.99998,
  0.99994, 0.99993, 0.99988.
- T_i / first output: 1, 1.00004, 1.00009, 0.99998, 0.99999, 0.99986, 0.99983, 0.99985, 0.99987, 0.99983, 0.99985.
- Share of the histograms above u0 + 5 √θ and above 8 √θ: 0 for electrons and ions, at the start and at the end.
- The electron kinetic energy stays between 0.99916 and 1.00012 of its first value, without a trend (the mean
  electron v_x / v0 moves within −0.99982 … −1.00038, the mean E_z within ±1 % of E0z, the mean B_y is constant).
  The ion kinetic energy stays within 1e-5.
- div B: 1.0e-13 after 166,667 steps (short box above), 7.5e-13 after 5,000,000. The contract is 1e-12. V1 has
  the same number of steps. Two points do not give the growth law; a run with more steps may come close to the
  contract.
- Limits: one seed; a periodic box (no face, no wall, no shock).

### Measured step times
- 1 GPU, periodic (N1): 1.55 ms/step for 4.0 M particles (0.39 ns per particle and step), 0.91 ms for 2.0 M.
- 1 GPU, wall or open faces (the pilots below): 12.5–14.3 ms/step for 4, 6, 8 and 16 M particles at t = 0. The
  step time of these runs does not follow the particle number; what it follows was not looked at. (One run
  without the Marder zone ran at 9.7 ms/step.)
- 8 GPUs, V1's box (270,000×2×2 cells, wall, face, `-n 8`, uniform plasma at t = 0, no seeds):
  540 M particles (2000 per column) 37.06 ms/step; 2.16 G particles (8000 per column) 136.79 ms/step.
  Two points give: step ≈ 3.8 ms + 0.493 ns × (particles on one rank).
- Own estimates for V1 (5.0 × 10⁶ steps) from that line. Assumptions: the wall rank fills at compression 4 behind
  a shock at the MHD speed; uniform decomposition; no seed species; the line holds for a compressed rank
  (not measured).

  | V1 variant | days |
  |---|---|
  | as configured (2700 d_i, 2000 per column) | 6.6 |
  | 1000 per column | 3.4 |
  | box 1350 d_i | 3.7 |
  | as configured, if the load were balanced | 2.9 |

  The plan has 0.9 days and the entry of 2026-11-14 (second review) 1.8 days. One seed species adds 50 % to the
  particles. Not decided: the box length, the particles per column, the seeds.

### Observation: with the flux face the tangential field at the face decays from 0.4–0.6 / Ω_ci on
The planner of eef89a4 used a `flux` face at x+ (the recommendation of PIC-D26 for a wave table). In the piston
pilot the upstream field was gone at the end: B_y / B0y = −0.005 and E_z / E0z = 0.06 over 22.8–37.0 d_i at
t = 2, with Gauss 2.1e-14 and div B 5.6e-13 (the conservation contracts hold; they do not see this).

B_y / B0y − 1, mean over the 0.5 d_i next to the face:

| t | piston | thru | thru-ns | wall-ns | wall-80 | wall-pml48 | wall-f0 | wall-m0 | flux-ppc500 | flux-th8 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.2 | −0.000 | −0.003 | −0.003 | −0.003 | +0.009 | −0.006 | +0.026 | −0.011 | −0.009 | +0.003 |
| 0.4 | +0.023 | +0.030 | +0.030 | +0.030 | +0.026 | +0.024 | −0.523 | +0.024 | +0.017 | +0.021 |
| 0.6 | −0.502 | −0.507 | −0.507 | −0.512 | −0.337 | −0.389 | −0.458 | −0.225 | −0.586 | −0.835 |
| 0.8 | −0.653 | −0.744 | −0.790 | −0.713 | −0.777 | −0.801 | −0.438 | −0.259 | −0.820 | −0.865 |
| 1.0 | −0.699 | −0.622 | −0.743 | −0.781 | −0.792 | −0.799 | −0.707 | −0.570 | −0.762 | −0.807 |
| 1.2 | −0.791 | −0.865 | −0.867 | −0.655 | −0.732 | −0.768 | −0.797 | −0.624 | −0.550 | −0.725 |

Runs: `piston` as configured (wall, one seed species); `thru` x− open instead of the wall (no shock);
`thru-ns` the same without seeds; `wall-ns` wall, no seeds; `wall-80` box 80 d_i; `wall-pml48` absorber 48
cells; `wall-f0` no current filter; `wall-m0` no Marder zone; `flux-ppc500` four times the particles
(no seeds); `flux-th8` θ_Bn 8° (no seeds).

In the same zone and the same runs:
- n_i − n_e rises from 0.001 (t = 0.2) to 0.02–0.09 n0 (t = 1.2): an electron deficit at the face.
- The electron energy per particle rises from 0.0068 (loaded) to 0.025–0.054.
- In `piston`, within 0.5 d_e of the face, E_x grows from t ≈ 0.40 on; the front of the changed field moves
  inward at about the flow speed.

What the table excludes as the cause: the wall and the shock (`thru`), the seed species (`thru-ns`,
`wall-ns`), the box length (`wall-80`), the absorber depth 24 → 48 (`wall-pml48`), the current filter
(`wall-f0`; without it the decay starts earlier), the Marder zone (`wall-m0`; without it the decay is slower
but reaches −0.62), the field angle (`flux-th8`), particle noise at 1000 → 4000 per column (`flux-ppc500`;
one run).

Readings that remain, none tested further:
1. The source sheath that PIC-D13 records for a flux face when the flow is supersonic for the electrons
   (there: u / v_te = 2.5, density 1.51 at the face). Here u / v_te = 3.0. The charge layer and E_x at the face
   fit this reading; PIC-D13 does not mention a loss of the tangential field.
2. An interaction of that layer with the incident-field surface of PIC-D26 (the tangential field enters through
   the face plane; a current layer on it changes what enters). Own thought.
3. Something specific to this configuration (external B_x together with `B0` / `E0`). The PR #52 test of the
   oblique piston runs 71 steps, so it does not reach this time.

### The paired face holds the field; B_y at the face rises slowly
The same pilots with `faces: {x+: paired}`, `neutralizer: e` (the inflow face of a wind in PIC-D13).
B_y / B0y − 1 in the 0.5 d_i next to the face:

| t | paired | paired-ns | paired-th8 | thru-paired | paired-80 | piston-paired |
|---|---|---|---|---|---|---|
| 0.2 | +0.009 | +0.009 | +0.020 | +0.009 | +0.003 | +0.009 |
| 0.4 | +0.004 | +0.004 | +0.006 | +0.005 | +0.007 | +0.004 |
| 0.6 | +0.013 | +0.011 | +0.009 | +0.011 | +0.014 | +0.019 |
| 0.8 | +0.002 | +0.004 | −0.028 | +0.005 | +0.015 | +0.003 |
| 1.0 | +0.019 | +0.024 | +0.021 | +0.020 | +0.014 | +0.025 |
| 1.2 | +0.031 | +0.027 | +0.039 | +0.026 | +0.021 | +0.028 |

Runs: `paired` wall and one seed species; `paired-ns` wall, no seeds; `paired-th8` θ_Bn 8°, no seeds;
`thru-paired` x− open (no wall, no shock), no seeds; `paired-80` wall, box 80 d_i, no seeds; `piston-paired`
the committed planner (below). The electron energy per particle at the face at t = 1.2 is 0.0071 (loaded 0.0068).

At t = 2, the 0.5 d_i next to the face:

| quantity | paired-ns (wall, 40 d_i) | thru-paired (no wall) | paired-80 (wall, 80 d_i) |
|---|---|---|---|
| B_y / B0y − 1 | +0.049 | +0.043 | +0.035 |
| B_z / B0y | −0.097 | +0.000 | −0.005 |
| electron energy per particle | 0.0090 | 0.0068 | 0.0069 |
| n_i − n_e | +0.001 | −0.001 | −0.001 |

- B_y at the face against time (paired-ns / thru-paired / paired-80): 0.6: +0.011 / +0.011 / +0.014;
  1.0: +0.024 / +0.020 / +0.014; 1.4: +0.030 / +0.024 / +0.022; 1.8: +0.051 / +0.038 / +0.031;
  2.0: +0.049 / +0.043 / +0.035. At 2–3 d_i from the face at t = 2: +0.059 / +0.036 / +0.038; at 10–11 d_i:
  +0.031 / +0.022 / +0.017; at 20–21 d_i: (shock) / +0.004 / +0.014.
- Electron energy per particle at t = 2, 10–11 d_i from the face: 0.0122 / 0.0067 / 0.0071; at 20–21 d_i:
  0.0366 / 0.0067 / 0.0074.
- `thru-paired` at t = 2, over the 64 d_e next to the face: B_y / B0y 1.035–1.045, E_z / E0z 0.985–1.031,
  the electron and ion v_x within 0.7 % and 0.05 % of v0, density 0.976–1.016 n0 in the first cells and
  1.000–1.001 beyond 8 d_e, electron energy 6.70–6.81e-3.

What these numbers show:
- B_z and the electron heating at the face appear only with the wall, and later in the longer box. They are
  absent without a shock. So they come from the shock side (particles or waves of the shock reaching the
  face), which is the matter of the box length (N6), not of the face rule.
- The rise of B_y (about +0.02 B0y per 1/Ω_ci up to t = 2, E_z with it, spreading inward) is there without a
  wall. It belongs to the face side (the paired rule, the incident-field surface or the absorbing layer).
- These runs end at t = 2 and do not say whether the rise saturates. V1 runs to 150 / Ω_ci. The next section
  follows the rise to t = 8.

### The rise of B_y follows the depth of the absorbing layer
Flow-through runs (scratch `thru_paired.py`: x− open, paired face at x+, no wall, no shock, no seeds; 40 d_i
unless noted; the same seed in all). Slope of B_y / B0y − 1 per 1/Ω_ci (straight-line fit), mean over the
10 d_i next to the face:

| run | npml | change | fit t = 1–4 | fit t = 1–8 | fit t = 4–8 | last value |
|---|---|---|---|---|---|---|
| tp-pml12 | 12 | | +0.3080 | | | +1.128 at t = 4 |
| tp-long | 24 | | +0.0193 | +0.0189 | +0.0189 | +0.145 at t = 8 |
| tp-m0 | 24 | no Marder zone | +0.0193 | | | +0.069 at t = 4 |
| tp-ppc500 | 24 | 4 times the particles | +0.0191 | | | +0.070 at t = 4 |
| tp-80 | 24 | box 80 d_i | +0.0180 | | | +0.069 at t = 4 |
| tp48-long | 48 | | +0.0014 | +0.0010 | +0.0010 | +0.0075 at t = 8 |
| tp96-long | 96 | | +0.0003 | −0.0001 | −0.0001 | −0.0013 at t = 8 |
| tp192-long | 192 | | +0.0002 | −0.0002 | −0.0001 | −0.0018 at t = 8 |

- The whole box (0–37 d_i from the face), fit over t = 2–8: +0.0189 (24), +0.0010 (48), −0.0001 (96),
  −0.0002 (192); rms about the fit 0.0003. The zone 30–37 d_i, t = 4–8: +0.0186, +0.0007, −0.0004, −0.0005.
  The rise is the same over the whole box.
- E_z / E0z − 1 has the same slopes within its errors (rms about the fit 0.005–0.014): E_z = v0 B_y holds.
- B_y / B0y − 1 in the 0.5 d_i next to the face at t = 0.62, 1.23, …, 8.00 (13 outputs):
  - 24: +0.011, +0.021, +0.039, +0.042, +0.059, +0.079, +0.083, +0.100, +0.100, +0.121, +0.124, +0.140, +0.156
  - 48: +0.000, −0.001, +0.006, −0.002, +0.004, +0.014, +0.006, +0.013, +0.003, +0.011, +0.004, +0.008, +0.014
  - 96: −0.001, −0.003, +0.004, −0.005, +0.001, +0.010, +0.001, +0.007, −0.003, +0.004, −0.004, −0.000, +0.005
  - 192: −0.001, −0.003, +0.004, −0.005, −0.000, +0.010, +0.001, +0.007, −0.004, +0.004, −0.005, −0.001, +0.003

  The fluctuations are common to the runs (same seed). Differences from these lists: (48 − 96) grows by about
  +0.0011 per 1/Ω_ci; (192 − 96) reaches −0.0013 at t = 8.
- In the four runs to t = 8: B_z within ±0.007 B0y, electron energy per particle 0.0067–0.0068 (loaded
  0.0068), n_i − n_e within ±0.0016, the same at every depth.
- Step time: 13.1, 12.8 and 13.1 ms/step with 48, 96 and 192 cells (266,667 steps in 3500, 3423 and 3503 s),
  as with 24. The depth has no measurable cost.

What the numbers show, and what they do not:
1. With 24 cells the rise is linear to t = 8 and does not saturate. If it stayed linear it would be about
   2.8 B0y over V1's 150 / Ω_ci. 24 cells are not usable for V1.
2. The rate does not depend on the Marder zone, on the particle number (×4) or on the box length (40 → 80 d_i).
3. The rate depends on the depth: a factor 16 from 12 to 24 cells and a factor 14–19 from 24 to 48, about
   npml⁻⁴.
4. With 96 and 192 cells the fitted slopes lie between −0.0002 and +0.0003 per 1/Ω_ci for the 10 d_i next to
   the face and for the whole box, with either sign depending on the window; in the smaller zone 30–37 d_i over
   t = 4–8 they are −0.0004 and −0.0005. The errors of these fits were not computed. In 8 / Ω_ci the rate is not
   distinguished from zero at the level of 0.0003–0.0005. A linear rate of 0.0003 would be 4.5 % over
   150 / Ω_ci, one of 0.0005 would be 7.5 %. The npml⁻⁴ line from 24 cells would give 7e-5 (1 % over
   150 / Ω_ci) with 96 and 5e-6 (0.07 %) with 192; that is an extrapolation, not a measurement.
5. With 48 cells the remaining +0.0010 per 1/Ω_ci would be 15 % over 150 / Ω_ci, if linear.
6. The mechanism is not identified. An own, untested thought: the layer of the engine has no frequency shift
   (`grid._init_pml`: coefficients `exp(-σ dt)`, `exp(-σ dt) − 1`), and PIC-D16 says that it does not relax a
   static field; a uniform (δB_y, δE_z) scattered into the region beyond the face would then stay and add up.
   This reading does not by itself give the dependence on the depth.
7. Limits: V1's plasma only; no wall; 8 / Ω_ci of V1's 150. How the rate depends on the plasma, the mass ratio
   or dt is not known, so another run group may need its own control.
8. `docs/plans/2026-11-14-oblique-piston-test-injection.md` (main) says of this drift "gone with `npml ≥ 24`".
   That held for the 71 steps of the PR #52 test; over 4–8 / Ω_ci it does not hold with 24 cells. The file is
   main's and is not changed here.

### Planner change: the face rule is a parameter, default `paired` (ba98d45)
- `TopicParams.face`: `paired` (default) or `flux`; with `paired` the electrons are the neutraliser; seeds enter
  as at a flux face with either rule (PIC-D13 as amended by PR #52).
- Limits of `paired`, stated with it:
  - PIC-D13: "at a paired face the neutraliser enters only as returned particles and partners, not with its own
    thermal flux." The electrons that enter are not an independent draw of the upstream distribution.
  - PIC-D26 recommends `flux` for a wave table, because a paired face ties the electrons' birth rate to the
    ions' flux. The topic's table has zero amplitude; a run with upstream waves (not in the plan's table) would
    meet this.
  - PIC-D13 records a 9 % density gain with the paired rule where the normal velocity is thermal. Here the flow
    is 3 v_te and the density at the face stays within about 2 %.
  - The rise of B_y is not a limit of the rule as far as measured: its rate follows the depth of the absorbing
    layer (section above). Whether it would also appear with a flux face was not tested (that face loses the
    field first).

### Recorded pilot: the committed planner (`runs/2026-11-15-pilot-piston-paired.md`)
- Commit ba98d45, `configs/pilot-piston.yaml` as it is (wall, paired face, one seed species), t = 2, 1 GPU,
  868.6 s. Gauss 2.3e-14, div B 5.7e-13. Particles 6.00 M → 10.12 M.
- Shock at 18.7 d_i at the end; speed over the second half 0.0457 c and compression 2.52 against 0.0337 c and
  3.96 of the MHD jump. In every piston pilot of this entry (1.2–2 / Ω_ci) the shock is faster and less
  compressed than the MHD jump; at this age it is still forming, and nothing more is read from it.
- Upstream (21.7–37.0 d_i) at the end: n_e 1.046, v_x / v0 0.945, B_y / B0y 1.043, E_z / E0z 1.034. Reflected
  particles are in that region (moment "temperatures" 6.5 and 42.8 times the loaded ones).
- Seeds: the weight in the box goes 16.0 → 24.96 (plasma species 16.0 → 28.0), so seeds enter with the flow and
  the upstream keeps them. Share of the seeds above u0 + 5 √θ: 1.35e-5 → 6.7e-3.

### Planner changes after these pilots
- 92686fb: the absorbing layer of the x+ face is 192 cells deep by default (`TopicParams.npml`; before: 24).
  The depth is a config value of the planner, so this needs no engine change. Chosen between the two depths
  that the data do not separate (96, 192) because the depth has no measurable cost. The docstring of `plan.py`
  gives the rate with 96 and 192 cells as "within ±0.0003"; that is the range of the fits next to the face and
  over the whole box (point 4 above has the far zone too). It is replaced when the control below has its
  number. `runs/index.json` added
  (`check_branch.py` asks for it; the earlier record is listed).
- 5bf6998: `through: true` opens x− instead of the wall (no shock), so that the control of the upstream state
  is a planner configuration; `pilot_check.py` gives the field deviations against time and their slopes;
  `configs/pilot-through.yaml` (V1's plasma, 40 d_i, 40 / Ω_ci).

### Recorded pilot with 192 cells (`runs/2026-11-15-pilot-piston-pml192.md`)
- Commit 92686fb, `configs/pilot-piston.yaml` as it is, t = 2, 1 GPU, 864.4 s (13.0 ms/step; 868.6 s with 24
  cells). Gauss 2.1e-14, div B 5.6e-13. The same number of particles entered (4,787,763).
- Shock at 19.9 d_i at the end; speed over the second half 0.0479 c, compression 2.56 (24 cells: 18.7 d_i,
  0.0457 c, 2.52).
- Upstream (22.9–37.0 d_i) at the end: n_e 1.027, v_x / v0 0.971, B_y / B0y 1.020, E_z / E0z 0.996 (24 cells,
  21.7–37.0 d_i: 1.046, 0.945, 1.043, 1.034). The zones differ and both hold reflected particles, so this is
  not the measurement of the face's own state; the flow-through runs above are.

### Control started: the upstream state to 40 / Ω_ci with 192 cells
- `configs/pilot-through.yaml` at commit 5bf6998 (clean tree), 1 GPU, 1,333,334 steps, output
  `/mnt/nvme/pypic/qpar/pilot/through-pml192`. Started 2026-11-15 04:09 JST; about 4.9 h at the measured
  13.2 ms/step.
- Fixed before the run: the quantity is the slope of B_y / B0y − 1 (mean over the box without 3 d_i at each
  end) per 1/Ω_ci with its standard error, fitted over t ≥ 1 and over the second half, and the value at the
  end (`pilot_check.py`, block `through`). The run measures the rate; it has no pass threshold. Its result
  times 150 is the stated limit of V1's upstream field.
- The result goes into `runs/` and a new entry.

### What step 4 has checked, against the list of the plan
| item of the plan | state |
|---|---|
| piston shock with a field component normal to the wall | runs on main 70202df; Gauss and div B hold; a shock forms |
| injection in a long quasi-1-D box, `compact` | flux face: the upstream field decays; paired face with 192 cells: fitted rate of B_y between −0.0005 and +0.0003 per 1/Ω_ci to t = 8 (no wall); control to t = 40 running |
| test species | enter through the face and stay in the upstream (PR #52) |
| numerical heating (N1) | with the filter none within 2e-4 in 5 / Ω_ci and in 150 / Ω_ci; no tail |
| measured ms/step | above; own estimate for V1: 6.6 days as configured, 3.4 days with 1000 per column (plan: 0.9 days) |

### Open after this entry
- The control to 40 / Ω_ci (above).
- V1's cost against the plan: the box length, the particles per column (2000 in `configs/v1.yaml`, about 1000
  in the plan) and the seeds are not decided.
- div B after 5.0e6 steps is 7.5e-13 against the contract 1e-12 (N1 at full length).
- The flux face: the cause of the decay is not identified (three readings above). The topic does not use it.
- The mechanism of the rise of B_y, and its rate for the other run groups (limit 7 above).
- Unchanged from the entry of 2026-11-14: the second definition of the injection rate and the length of the
  downstream region before V1; the single T_e before P0.

## 2026-11-14 — Step 4, first part: the engine falls short in three places (PR #52); the topic planner

Later on the same day as the entry below. Step 4 of the plan checks the engine with short pilots. All runs of this
entry are scratch runs (diagnostic, not recorded in `runs/`): 1 GPU, 1500×2×2 cells, mass ratio 25,
ω_pe/Ω_ce 4, flow 0.18 c, 20000 steps, engine main fa09726 unless noted.

### What the pilots found on main fa09726
1. **A `B0` normal to a reflecting wall breaks div B and is not refused.** max |div B| = 1.1 at the wall against
   the contract 1e-12. Normal to an open face it is fine (1.6e-14).
2. **The piston with an oblique field cannot be configured.** The normal component must be an external field (the
   wall zeroes a normal Yee B); the tangential component and the motional E must be Yee fields (a wall does not
   short an external E); the face must carry the uniform fields in (PIC-D26 table of zero amplitude). The engine
   refuses `external_fields` together with `B0` / `E0`. With that one refusal bypassed the run holds Gauss
   ≤ 1.2e-14 and div B ≤ 2.2e-13 at 30° and at 8°, and a shock forms.
3. **Test species cannot enter through an injecting face.** A test species loaded in the box is swept into the
   wall and is not replaced: after one crossing time the upstream holds no seed electrons.

### The engine change (plan step 4, the route named there)
- Branch `engine/qpar-piston` from main, draft PR #52, head fef8bde: `external_fields` together with `B0` / `E0`
  (PIC-D12, PIC-D18 amended); a `B0` component normal to a reflecting face is refused; a test species listed in
  `injection.species` enters as at a flux face, with no deposit (PIC-D13, PIC-D9 amended). No kernel change.
  Detail and test results: `docs/plans/2026-11-14-oblique-piston-test-injection.md` on that branch.
- State when this entry was written: CPU suite and the targeted GPU tests pass; Codex review round 1 reports no
  findings; the full GPU gate on fef8bde is running. Not merged. Recorded runs wait for the merge and
  `research-sync-main`.

### Facts for the set-up of the runs
- **Absorber depth at the PIC-D26 face.** With the default `npml = 12` the upstream state drifts next to the face
  within 20000 steps: momentum 2.2e-4 u, B_y 7.4e-4, E_z 1.1e-3 (relative to the initial values). The same amounts
  appear with an open x− instead of the wall; with the whole field in `B0` (90°) the momentum is off by 5.1e-4 u.
  With `npml = 24`: 1.4e-5 u, 4.9e-5, 5.7e-5. The planner sets `npml = 24`. Not checked: longer runs, deeper layers.
- **A seed species costs as many particles as the electrons.** `SimConfig.ppc` is one number for every beam of
  every species. One seed species in a run adds 50 % to the particles. The plan's cost estimates (§6.2) do not
  include this. Not decided: fewer seeds through a second, shorter run, or the cost as it is.
- **No load balancing.** The decomposition is uniform, so the rank at the wall holds up to the compression times
  the upstream particles and sets the step time. The plan's cost formula does not include this either. Measured
  timings follow with the pilots.

### An observation at the wall (30°, scratch pilot)
- Measured: shock speed 0.078–0.079 c in the box and density compression 3.3. The fast shock of ideal MHD
  (γ = 5/3) for the same inflow: 0.0725 c and 3.48.
- Readings, none tested:
  1. A conducting wall forces E_z = 0 at the wall, while the oblique jump leaves a tangential flow and a
     tangential field downstream. A layer between the wall and the fast shock would then take part of the jump
     (own reasoning).
  2. The shock is young (20000 steps, 300 d_e of box) and has not reached a steady state.
  3. The downstream of a kinetic shock at mass ratio 25 with 2 cells across need not have γ = 5/3 or the MHD
     pressure balance (energetic particles, anisotropy).
- `gupta2024electron` uses a reflecting piston with the field at 30° as well; what its downstream looks like next
  to the piston was not looked up.
- Consequence for the runs: the measured shock speed and compression of each run are recorded and used; the MHD
  values of the planner's report are only the expectation.

### Conventions used for V1 (run B1 of `gupta2024electron`)
- From the source: M_A and M_s are defined on v_sh = 4 v_pt / 3; M_s = v_sh / √(kT/m_i); β = 4 (M_A/M_s)² is the
  total of both species with T_e = T_i; dx = d_e/10, c dt/dx = 0.45, 200 particles per cell per species, 5 cells
  across, 32 passes of the current filter.
- Not found in the source: ω_pe/Ω_ce, and whether v_A holds the ion mass only.
- Own working numbers (convention `r4` of `plan.py`; v_A with the ion mass only): b = Ω_ce/ω_pe = 0.0667,
  v_A = 0.00667 c, Ω_ci = 6.67e-4 ω_pe, T_e/m_e c² = 1.1e-3, Debye length 0.33 Δx, electron sonic Mach number 4.0.
  The MHD jump at these numbers gives r = 3.96 and M_A = 20.06 on its own shock speed: the two conventions differ
  by 0.3 % here.
- The Debye length is a third of a cell. How much the grid heats this plasma over 5 × 10⁶ steps is what N1
  measures; the pilot N1 runs come next.

### Topic code
- `plan.py`, `run.py`, `configs/v1.yaml`, `configs/pilot-piston.yaml`, `configs/pilot-n1.yaml` (commit eef89a4).
  They run only on an engine with PR #52.
- V1 at 2000 particles per column in a fixed box of 2700 d_i: 2.7 × 10⁵ cells along x, 5.4 × 10⁸ particles at
  t = 0, 5 × 10⁶ steps. The box length is not fixed; it follows the measured timings.

## 2026-11-14 — Second review applied: inputs corrected, criteria and thresholds fixed

Later on the same day as the entry below. The second review by `research-reviewer` (README, `inputs.md`, plan) was
reported to the owner with its proposals and the open points.

### Decision
- The owner, 2026-11-14: proceed with the recommendations. Read as: the recommendations of the second review are
  approved.
- Where the review gives no recommendation, the approved plan stays as it is. No new fact requires a change there.
- Where the review requires an item to be settled and names no value, no value is invented. The item is marked
  "not yet fixed" with the facts and a point by which it must be fixed.

### The nine corrections to `inputs.md`, each checked in its source before it was written
1. **PSP electron spectrum.** δ ~2.7 from 50 keV to 2.5 MeV, the ~6.7 MeV and the isotropy are in
   `jebaraj2024direct`, not in `jebaraj2024acceleration`, which says "to 6 MeV". It is a 10-minute average upstream
   (07:03-07:13 UT). Above 2.5 MeV δ ~4.18 with "increased uncertainties due to the instrument's dynamic threshold
   mode". The isotropy is read from two telescope directions with similar intensities.
2. **Saturn speed.** ~450 km/s in `masters2017situ` is a typical value quoted for the quasi-perpendicular crossing.
   No speed is given for 2007-02-03. The 500 km/s is the U_1 of the test-particle model of `guo2015acceleration`.
3. **Saturn crossing time.** ~01:05 UT in `masters2017situ`, ~01:10 UT in `masters2013electron`. The connection at
   the spacecraft is ≥ 60-65 min, ω_ci t ≥ 34-37. The reviewer's "01:08" was not found in our check; only ~01:10
   is written.
4. **PSP compression.** The ~4 of `jebaraj2024acceleration` follows from 184 ± 35 over 46.5 ± 15 cm⁻³ (own
   propagation: 4.0 ± 1.5). In `kouloumvakos2026wave` 4 is the limit "expected at M_A ≃ 7.5", not a measurement.
5. **M_A 4.68 ± 3.45.** It rests on the SPAN-i density, which `chiappetta2026evolution` calls "a lower limit
   estimation, at least for the upstream side". It is biased low, not an equal third estimate.
6. **Sizes.** The 15-25 R_E and 50-60 R_S of `oka2025maximum` are adopted characteristic sizes. Their kind is
   "assumed"; ω_ci t from size over speed is a crossing-time proxy.
7. **Earth M_A.** `lalti2022database` writes "the highest counts between M_A = 5 and 15", for all crossings. The
   same work has one crossing with θ_Bn 87.3° from the model and 23° from local data, which bears on Saturn's
   model normal.
8. **Earth "no".** "only noise > 12 keV" describes three example distributions of `wilson2016relativistic`. Of the
   10 enhancements 5 are at SLAMS, 2 at hot flow anomalies, 3 at foreshock bubbles.
9. **Kartavykh's second table value** is the midrange, not the median: M_A 5.8 = (3.2 + 8.4)/2 against a median of
   6.35; θ_Bn 65° = (88° + 42°)/2.

The notes `kartavykh2025statistical`, `lalti2022database` and `wilson2016relativistic` carry corrections 7-9 too.

### Added to `inputs.md` beside the corrections
- **Convention:** M_se = u_1/√(kT_e/m_e), as in `gupta2024electron` and in the plan's formula. Own values: PSP 1.28
  (T_e 20 eV) and 0.90 (40 eV) at 2400 km/s; Earth 0.27-0.41; Saturn 0.34 at T_e 10 eV and 1.07 at T_e ≈ 1.0 eV,
  the value that M_s 25 implies with γ 5/3 and T_e = T_i.
- **PSP densities:** 46.5 ± 15 and 184 ± 35 cm⁻³ from the plasma frequency; SPAN-i ~22 cm⁻³, "a lower limit".
- **DSA and the measured index (own):** for relativistic electrons q = δ + 2 = 3r/(r − 1); δ 2.7 gives r ≈ 2.8.
  The authors: "DSA predicts δ=2".
- **Scale of the PSP precursor (own):** ~57 d_i per second at the spacecraft; a box of 10³ d_i is the last ~17 s.
- **Saturn:** e V B L at the event's 0.1 nT is 87-210 keV (own), below the observed ~1 MeV; the shock motion of
  100 km/s is an assumed value; the signal is ≥ 10 times the seeds below 100 keV; 28 of the 31 crossings with a
  signal have it below 100 keV only; the authors hold that injection "may have been taking place at all these
  crossings".
- **1 au:** the four shocks left out of the averages of `kartavykh2025statistical`, among them the fast
  quasi-parallel one of 2022-09-06 (M_A 7.4, θ_Bn 30°).
- **Backgrounds:** split into what is collected as statements and what is not collected.

### Fixed in the README (Approach), from the review's table, with its reasons
- Reproduction: V1 slope within ±0.3 of 4; D_1,e in [0.3, 1.2] D_B0 at t = 135 and p = m_i v_pt; D_2/D_1 in
  [0.05, 0.2]. V2: no slope criterion; p_max and η_nt,e within a factor 2.
- σ_pt from the seed pair and ≥ 5 time blocks of ≥ 20 ω_ci⁻¹.
- N1 at each v_sh/c, signal ≥ 5 times the floor and ≥ 100 counts per bin. N2 ≤ 30 % and ≤ 0.15 in slope. N3 and
  N6 ≤ 2σ_pt. N4: residual ≤ 2σ_pt, "decided" only beyond twice the extrapolated σ. N7 within 10 % and a log-slope
  of 1 ± 0.1. N8: slope < 0.3 gives "stationary within the run" and no transfer beyond the run's age.
- Five tests added: N9 frozen against tracked (factor 1.5), N10 test-electron domain doubled, N11 without the wave
  electric field, N12 electrons that touched the piston region left out, N13 V1 at Gupta's particle number.
- Confrontation: H1 "yes" at ≥ 3 times the background at a named energy, "no" at ≤ 1/3, otherwise "not decided";
  the observed side per case; H2/H3 orders only above 3σ_pt; H6 in three classes (2/3 and 1/3); "differ strongly"
  for 1-D against 2-D as a factor 2 or 3σ_pt.
- Own checks of the review's numbers: the N4 factor is 2.24 at 0.008 c and 3.6-4.1 at 0.0020-0.0013 c (review: 2.2
  and 3.9); Δx/λ_D is 0.75, 1.49 and 3.03 at 0.133, 0.067 and 0.033 c for M_se 1; the downstream grows by M_A/4 d_i
  per ω_ci⁻¹; Eq. 1 and Fig. 2 of `gupta2025speed` give the per-cycle gain that H6 uses.

### Changes against the approved plan, each with the fact that requires it
| Fact | Change |
|---|---|
| At the events' M_se (0.27-1.28) a large share of the thermal electrons already moves faster than the shock | The injection rate is no longer "the share passing p_inj". First definition: η_nt,e of `gupta2024electron`, p ≥ 5 p_th,e downstream |
| `gupta2024electron` uses 200 particles per cell per species over 5 cells, 2000 per column; the plan has ~1000 | V1 runs at 2000 per column. Its estimate goes from 0.9 to 1.8 days (cost linear in the particle number) |
| Δx/λ_D changes from 0.75 to 3.0 along the v_sh/c scan | N1 is run at each speed, not once |
| Carrying D̃ from t = 200 to 1.1 × 10⁵ within a factor 2 needs a slope ≤ 0.11; a run resolves ~0.3-0.5 | N8 labels a value "stationary within the run" only. No value is transferred to the age of the PSP shock on N8 alone |
| The published M_A of the PSP crossing (4.68 from a lower-limit density, 7-7.5, 9.1) lie at or below the M_A ≳ 10 of `gupta2024electron` | The PSP result is stated "for M_A ≤ 7.5". A run at M_A 9-10 is decided after P0 |
| Run B2_II (M_A 5, M_se 0.9, β 0.06, mass ratio 1836) is close to the PSP conditions | Novelty part 1 names the step from B2_II: M_A 5 → 7.5, θ_Bn 30° → 8°, β_i 0.03 → 0.2 |
| The observed "no" is an upper limit in every case (Saturn: channel backgrounds; Earth: three examples, bursts at SLAMS; 1 au: proton contamination) | Novelty part 2 and the confrontation are worded as intensities and upper limits at named energies |
| The ~10 % and 5-6 % are values of the 58-minute precursor in 0.05-1 Hz; a box covers the last ~17 s | δB/B of a run is compared in the same window and band only |
| Without n and T_e no absolute intensity follows at Saturn | S1-S3 decide orderings, not an intensity, until those inputs have a source |
| W stands for an average over all geometries | W states a result for an M_A 3 shock at 30° only |

### Kept as planned (no recommendation in the review, no new fact)
- The frame (downstream), with N1 as its check.
- ~1000 particles per column in the runs other than V1, the same at each speed of the scan.
- H4 with passive seeds; the waves that seeds would drive are named as a limit.
- The run list; the border regimes of `simulation-claims.md` §6.3 outside it are named as limits.
- P0 at M_A 7.5 and θ 8°.

### Not yet fixed (the review requires them and names no value)
- The second definition of the injection rate, and the length of the fixed downstream region. Facts: Gupta's
  1000 d_i "at late times" is not reached (750 d_i in V1 at t = 150, 375 d_i in P0 at t = 200). Before V1.
- The single T_e of P0. Facts: 20 ± 5 eV adopted in `jebaraj2024acceleration`, "≲ 40 eV, unresolved" in
  `jebaraj2026wave`; M_se 1.28-0.90. Before P0.
- Which B enters H1b at Saturn. Facts: 87-210 keV at 0.1 nT (own), 0.2-2 MeV for 0.2-1.0 nT (`oka2025maximum`).
  Before the confrontation.
- Thresholds of N10 and N11: none given; the measured change is reported.
- Whether wave-growth-limited confinement is a separate candidate for electrons (README, Research question).

### Slips of our own, caught while applying
- The DSA relation was first written as q = 2δ + 2, the non-relativistic form. Corrected to q = δ + 2.
- "δB/B ≥ 0.5 close to the ramp" was written before the check. The source says it of the shock transition
  "preceding the ramp-overshoot region"; the text now quotes it.
- "the 477 crossings without a signal" was our arithmetic, not a number of the source. Removed.

### Not checked, carried from the review (can wait until stage 3 or the write-up)
- Sources: `dresing2016efficiency`, `talebpoursheshvan2023role`, `yang2019electron`, `tsurutani1985acceleration`
  against the inputs; the M_A 9.4 / 33° example of `lalti2022database`; the full texts of `masters2011electron`,
  `went2011new`, `achilleos2006orientation`.
- Inputs: Saturn's density and T_e; a connection time at Earth; the quasi-parallel subset near 1 au; the Solar
  Orbiter backgrounds.

## 2026-11-14 — Step 3, second part: what a simulation can state, and the boundary areas

Later on the same day as the entry below. The first part of step 3 collected the observed side only. The research
is done by simulation, so the question "what can a simulation state here, and what not" has to be surveyed with the
same care, together with the boundary areas. This entry records that survey. Result file: `simulation-claims.md`.

### What was done
- **Extraction from what was already read.** The notes and stored texts of the fully kinetic works on quasi-parallel
  shocks (`gupta2024electron`, `gupta2025speed`, `park2015simultaneous`, `kato2015particle`, `shalaby2022mechanism`,
  `zekovic2025slams`, `fang2019early`, all on `research/spherical-shock-injection`; `fang2022kinetic`,
  `yu2024nonthermal`, `bessho2023electron` here) were gone through for four things: what each states, on what basis,
  which limits its authors name, and how it is tied to an observation.
- **A second search, 26 keyword queries on OpenAlex** (10-12 results each, with abstracts):
  - methods and numerics: "particle-in-cell collisionless shock simulation mass ratio dependence electron
    acceleration"; "numerical heating particle-in-cell shock simulation number of particles per cell electron
    acceleration"; "shock simulation injection method reflecting wall piston comparison particle-in-cell";
    "three-dimensional particle-in-cell quasi-parallel shock electron"; "finite grid instability numerical heating
    drifting cold plasma particle-in-cell"; "simulation box transverse size cosmic ray current filamentation shock
    convergence"; "semi-implicit particle-in-cell collisionless shock realistic mass ratio";
  - neighbouring methods: "hybrid simulation test particle electrons quasi-parallel shock acceleration"; "MHD-PIC
    shock electron acceleration injection"; "magnetohydrodynamic-particle-in-cell method cosmic rays shock Bai";
    "Vlasiator electron foreshock hybrid-Vlasov electrons"; "eVlasiator electrons hybrid-Vlasov global simulation";
    "Vlasov-Fokker-Planck shock precursor cosmic ray electron injection simulation"; "continuum Vlasov-Maxwell
    collisionless shock electrons simulation field-particle correlation"; "dHybridR hybrid relativistic ions shock
    simulation"; "focused transport electron shock acceleration interplanetary model relativistic electrons"; "free
    escape boundary diffusive shock acceleration maximum energy time-dependent finite size electrons heliospheric
    shock"; "minimal model ion injection quasi-parallel shocks reformation";
  - runs against events: "global hybrid simulation foreshock transient electron acceleration test particle";
    "comparison particle-in-cell simulation MMS observations quasi-parallel shock electrons"; "interplanetary shock
    particle-in-cell simulation observed parameters electron acceleration Solar Orbiter Wind"; "Saturn bow shock
    particle-in-cell simulation Cassini high Mach number";
  - laboratory and border regimes: "laboratory magnetized collisionless shock electron acceleration experiment";
    "shock obliquity dependence electron acceleration particle-in-cell quasi-parallel oblique transition"; "weak
    quasi-parallel shock electron preacceleration intracluster particle-in-cell critical Mach number"; "electron
    injection shock speed dependence scaling non-relativistic shock simulation extrapolation realistic shock
    velocity".
- **Screening and reading.** Abstracts of about 90 works were kept and screened. 29 notes were written:
  - read in the body (arXiv PDF text; the read range is in each note): `zhou2025planar`, `marcowith2020multi`,
    `simon2026faux`, `bai2015magnetohydrodynamic`, `vanmarle2022diffusive`, `haggerty2019dhybridr`,
    `caprioli2018diffusive`, `caprioli2015simulations`, `orusa2025criteria`;
  - read as abstracts: `swisdak2023comparison`, `kropotina2023weibel`, `zhang2024kinetic`, `tsubouchi2004full`,
    `johlander2022quasi`, `liu2017fermi`, `nyberg2026electron`, `john2024effects`, `battarbee2021vlasov`,
    `lezhnin2021kinetic`, `kang2019electron`, `reville2013universal`, `bessho2025electron`,
    `liseykina2019electron`, `kylivnyk2026hybrid`, `gingell2017mms`, `schaeffer2019direct`,
    `vanthieghem2024electron`, `lindberg2025fermi`, `fiuza2020electron`.
- **Sources that did not work:** no open full text for `swisdak2023comparison`, `zhang2024kinetic`,
  `tsubouchi2004full`; our PDF text extractor gave nothing for arXiv:2310.13779 (`kropotina2023weibel`); the journal
  page of `nyberg2026electron` could not be read. All were read as abstracts and are marked so.
- **Passages re-read in the stored texts:** `kato2015particle`, the drag on energetic electrons in 1-D (their Eq. 4
  and the particle numbers at which acceleration fails); `gupta2024electron`, the definition of M_A, the relation of
  piston and shock speed, the frame, the dependence of spectra and efficiency on the shock speed; `gupta2025speed`,
  p_max in time.

### What the survey shows (detail and sources in `simulation-claims.md`)
- No published simulation in what was read measures D̃_e or an electron injection rate at the conditions of our four
  cases, or confronts one of the events with numbers from self-consistent electron kinetics. This repeats the
  result of the first survey with a wider base; it is still not exhaustive.
- The limits that authors themselves name recur in eight kinds: dimension, mass ratio, v_sh/c, run time, box and
  boundaries (with the launching method and the frame), particle noise, extrapolation by a formula, and ties to
  observations that are indirect (a global ratio at a remnant, a slope, a mechanism seen in situ).
- Published ways of tying a run to an observation: a global ratio at remnants; a synthetic probe through the box;
  hybrid runs at event parameters; test electrons in prescribed or hybrid fields; a model in prescribed fields
  fitted to one event; a mechanism compared with in-situ data; transport from data; laboratory shocks; Monte Carlo
  with prescribed inputs.
- The methods next to PIC either have no kinetic electrons or take the electron injection (or the fields at
  electron scales) as an input.
- Facts that bear on statements of the plan:
  - `gupta2024electron` already reports the electron injection efficiency for v_sh ≈ 0.067-0.267 c and mass ratios
    up to 1836: it "does not show much dependence on shock speed and m_R", and drops by a factor ~5 at low Mach
    number. The energetic electrons become fewer at high speed; the authors state that 1-D exaggerates this.
  - `gupta2024electron` runs in the upstream rest frame, chosen against numerical heating; its efficiency is a
    factor ≲ 2 above the downstream-frame run of `park2015simultaneous`. Our set-up is in the downstream frame.
  - `kato2015particle`: in 1-D, energetic electrons lose energy at a rate set by the number of electrons per
    electron skin depth; with too few, no acceleration appears.
  - `caprioli2018diffusive`: seeds drive waves through their current and can switch on the injection of thermal
    ions. Our seeds are a passive species.
  - `vanmarle2022diffusive`: "injection is, fundamentally, a 3D process" (ions, oblique shocks).

### Our calculations (own; not results)
- **Numerical drag against a Bohm-rate gain.** Loss from Kato's Eq. 4, gain from the formula of `gupta2025speed`
  with compression 4 and D_2 = D_1/10. Table by run in `simulation-claims.md` §5: loss/gain is 0.02-0.11 for W, V2,
  P0, P2; 0.3-0.5 for P5a, V1, E; above 1 for S1, P5b, P6, S2, at 5000 electrons per d_e and D̃ = 1.
- **A correction of our own first number.** The first evaluation for Gupta's run B1 used the piston speed (0.1 c) in
  place of the shock speed. `gupta2024electron` defines M_A with v_sh and sets v_sh = v_pt R/(R−1), so v_sh is
  0.133 c. With it the formula gives 2.1e-4 m_e c² ω_pe, against a mean rate of 2.0e-4 from the published p_max at
  300 ω_ci⁻¹. The earlier value (8.9e-5, never written to a file) is void.
- **A tension that is not resolved.** The same formula gives loss/gain ≈ 0.24 for Gupta's own B1, while that work
  reports agreement with its p_max formula and no sensitivity to the particle number. Either the 1-D loss formula
  overstates the drag in a filtered quasi-1-D run, or the effect is inside their agreement. Not checked; a short
  test run with the particle number varied would decide it.

### Changes of the README, each with its fact
| Fact | Change |
|---|---|
| The second search and its 29 notes | Related work: "Second search" added, with a pointer to `simulation-claims.md` |
| `gupta2024electron` reports the injection efficiency against shock speed and mass ratio (quotation above) | "Not new": that item added; what remains open is named |
| `liseykina2019electron`: test electrons in 1-D hybrid quasi-parallel shocks | "Not new": added beside `hanusch2020electron` |
| `kato2015particle` Eq. 4 and our estimate | "Strongest objections": numerical drag |
| `gupta2024electron` on the frame, against `park2015simultaneous` | "Strongest objections": the frame |
| `caprioli2018diffusive` | "Strongest objections": passive seeds cannot drive waves |
| `vanmarle2022diffusive` | "Strongest objections": injection as a 3-D process |
| `inputs.md`: v_sh/c is 0.0013-0.0020 at Earth and 0.0015 at Saturn | "Out of reach": the range corrected (it said 0.0013-0.0015 for both) |
| The points of `simulation-claims.md` §8 | Approach: a line that names them as open, for the owner; set-up, definitions and criteria unchanged |
| The second review is done; this survey is done | Status: "Done" and "Next" updated |

`literature.md`: five Map bullets for the 29 works and a pointer to `simulation-claims.md`.

### Second review by `research-reviewer`
Done on README, `inputs.md` and the plan. It lists nine corrections to `inputs.md` (attributions, a crossing time,
the kind of some values), finds that the published M_A values of the PSP event straddle M_A ≈ 10, that three of the
four "no" classes are limited by detection according to their own authors, and that part of the dependence on shock
speed is published; it proposes thresholds for the criteria. None of it is applied in this commit: the corrections
are to be checked in their sources first, and the thresholds are proposals until the owner decides.

### Not done
- The search is keyword-based and not exhaustive; no forward-citation pass was made for the 29 new works.
- 20 of the 29 works were read as abstracts only.
- The engine's implicit scheme as a route to larger cells was not examined (`zhou2025planar` is the one work read
  on semi-implicit PIC for shocks).
- The drag estimate is not tested by a run.
- The reviewer's corrections to `inputs.md` and the instrument backgrounds as intensities.

### Next
The owner's decisions on the open points (`simulation-claims.md` §8) and on the criteria; the corrections to
`inputs.md`; then the engine pilot (step 4).

## 2026-11-14 — Step 3, first part: the inputs of the four cases from primary sources

Later on the same day as the entry below.

### What was done
- The values the runs and the transport model need were taken from the stored full texts and from five further
  works, and written to `inputs.md` with source and kind (measured, derived, model, assumed, own).
- Passages re-read for this: `jebaraj2024acceleration` (A), `jebaraj2026wave`, `kouloumvakos2026wave`,
  `jebaraj2026diffusive`, the table of `chiappetta2026evolution` (PSP); `masters2013electron`, `masters2017situ`,
  `masters2016suprathermal` (Saturn); `wilson2016relativistic`, `shi2025compound`, `liu2019relativistic`,
  `liu2026efficient` (A), the table of `oka2025maximum` (Earth, Saturn); table 1 of `kartavykh2025statistical`
  (now added to its note), `dresing2016efficiency`, `talebpoursheshvan2023role` (1 au).
- New notes: `lalti2022database` (body, arXiv PDF text); `dresing2022role`, `masters2011electron`,
  `achilleos2006orientation`, `went2011new` (abstracts; full texts not obtained).
- This closes the item "Open from the plan's list of unread works" of the entry below: Dresing et al. 2022 is
  identified (ApJL 925, L21) and read as an abstract.
- Looked at, no note: Sulaiman et al. 2015 (PRL 115, 125001; quasi-perpendicular high-M_A shocks at Saturn; abstract).

### What the sources say about the kind of each input
- **PSP:** B, density and T_i are measured. T_e is not ("≲ 40 eV, unresolved"; 20 ± 5 eV is adopted from earlier
  orbits). M_A is 7-7.5, 9.1 ± 1.35, 4.68 ± 3.45 or ~4. Reasons given by the authors: a coarser estimate (9.1); an
  underestimated upstream density from SPAN-i (4.68; the compression rises from 1.39 to ~2.54 with the Langmuir-wave
  density). The hybrid model with M_A ~4 and ~1350 km/s was read as an abstract only.
- **Saturn:** only B (~0.1 nT), the times and the electron spectra are measured. M_A ~100 comes from a propagated
  solar wind model, θ_Bn ~20° from a model normal, M_s ~25 from typical temperatures, the flow speed (~450-500 km/s)
  and the shock motion (100 km/s) are typical or assumed values. No density, no β_e.
- **Earth:** the M_A range (most crossings at 5-15) is from 2797 MMS crossings. The limit without transients is the
  noise above 12 keV. No source read gives a connection time.
- **1 au:** averages of 51 Solar Orbiter shocks over all geometries (M_A 3.3 ± 0.3, gas compression 2.0 ± 0.1,
  β 0.84 ± 0.12). No values for quasi-parallel shocks alone.

### Our calculations (unverified; not claims)
- ω_ci: 7.66 rad/s at 80 nT; 9.58 × 10⁻³ rad/s at 0.1 nT; 0.479 and 0.766 rad/s at 5 and 8 nT.
- PSP: 4 h at 80 nT is ω_ci t ≈ 1.1 × 10⁵; β_e 0.04-0.08 for T_e 20-40 eV.
- Saturn: the quasi-parallel connection at the spacecraft lasts ≥ 60 min before the crossing (00:05 to 01:05 UT):
  ω_ci t ≥ 34. Size over speed (50-60 R_S, 300-600 km/s) gives 46-112. M_A 100 at 0.1 nT and 450 km/s implies a
  density of 0.24 cm⁻³. M_A 100 with M_s 25 implies a total β of ~19.
- Earth: size over speed (15-25 R_E, 400-600 km/s) at 5-8 nT gives ω_ci t 76-305. β_e ≈ 1 for 12 eV, 5 cm⁻³, 5 nT.
- 1 au: 1-3 days at 5 nT is ω_ci t 4 × 10⁴ - 1.2 × 10⁵.

### Changes of the README, each with its fact
No decision of the plan (question, hypotheses, runs, proposed criteria) was changed.
- No source read gives Mach numbers for the 15 quasi-parallel shocks of `yang2019electron`; the measured values are
  M_ms 1.25-2.85 in other samples and M_A 3.3 ± 0.3 for 51 shocks of all geometries. → The case table no longer
  states "M_A ~1.5-3" for that subset.
- `lalti2022database`: most crossings of Earth's bow shock have M_A 5-15; `oka2025maximum` lists 400-600 km/s.
  → The Earth row carries the ranges; ~8 stays as the representative value of run E.
- Saturn's M_A and speed are model and typical values (above). → Said so in the case table and in the list of
  objections.
- The connection estimates of the entry below (Earth 120-240, Saturn 30-70) are replaced in the README by the values
  above (76-305; ≥ 34 and 46-112). The planned run times (E: 240, S1-S3: 50) lie inside them.

### Not done
- Instrument backgrounds as intensities for each case (needed for H8 and for the confrontation by intensity).
- Still not found: the Gupta et al. 2025 conference contribution; the journal reference of `jebaraj2026diffusive`; the
  exact title of `galishnikova2026self`.
- A caution for later reading: each entry of table 1 of `kartavykh2025statistical` has a second value that is not the
  median of the four listed M_A in our check; only averages and standard errors are used.

### Next
- Re-review of the README, `inputs.md` and the plan by `research-reviewer`; the owner's decision on the criteria and
  thresholds.

## 2026-11-14 — Topic started: decisions, survey, differences from the approved plan

### Decisions (owner)
- Theme: electron acceleration at strong quasi-parallel shocks (Saturn/Cassini, Parker Solar Probe, SN 1006).
- Question: what decides whether a quasi-parallel shock accelerates electrons to relativistic energy.
- The plan (question in operational form, hypotheses H1a-H8, runs, proposed criteria) was approved after one review
  by `research-reviewer`. The criteria and null tests stay proposals until the review of step 3 and the owner's
  approval.
- Earlier candidates of this session (four themes, and selection criteria set without the owner) were rejected by the
  owner as too small and driven by the simulation instead of by observations. They are not reused.

### Survey (plan step 1)
- **Forward citations** (OpenAlex and INSPIRE) of 15 seed works; 337 citing works after merging. Seeds:
  `masters2013electron`, `masters2016suprathermal`, `masters2017situ`, `jebaraj2024acceleration`,
  `jebaraj2024direct`, `jebaraj2026wave`, `kouloumvakos2026wave` (no citing work yet), `liu2019relativistic`,
  `wilson2016relativistic`, `raptis2025revealing`, `gupta2024electron`, `gupta2025speed`, `guo2015acceleration`,
  `shalaby2022mechanism`, `zekovic2025slams`.
- **Venue sweep** (OpenAlex, title and abstract, 2016-2026): ApJ, ApJL, ApJS, JGR Space Physics, GRL, MNRAS, A&A, PRL,
  Physics of Plasmas, Nature Communications. Eight queries: electron acceleration with "quasi-parallel" / with
  "parallel shock"; "electron injection" with shock; electron "diffusion coefficient" with shock and kinetic or
  particle-in-cell; "relativistic electrons" with bow shock or foreshock; energetic electrons with "interplanetary
  shock"; electron "diffusive shock acceleration" with particle-in-cell or hybrid; electron, foreshock and
  acceleration with Saturn, Jupiter, Mars, Mercury or Venus. 84 works not already among the citing works.
- **Reading:** all titles screened; 50 abstracts; the text of about 22 works in part or whole. Each note states its
  range in its `Read:` line. 37 notes were written (the files added to `research/literature/` together with this entry).
- **Result:** no work was found that measures the electron injection rate or D̃_e at the conditions of an observed
  event from self-consistent kinetics, or that confronts an event with a transport model built on such a measurement.
- **Limits of the search:** most of the 84 sweep works were screened by title only. ADS and Semantic Scholar did not
  answer (HTTP 405 and 429). Publisher pages of IOP, A&A (HTML), Elsevier and PubMed Central refused automated
  reading; where no arXiv or repository version existed, the abstract from OpenAlex or Crossref was read. Installing a
  PDF library was not permitted; PDF text was extracted with a small script of our own, which loses digits in some
  files (`grassi2023electron`: only statements checked in context are quoted).
- **Not found:** a 2025 conference contribution of Gupta et al. named in the plan; the journal reference of
  `jebaraj2026diffusive` (arXiv:2610.08034); the exact title of `galishnikova2026self` (arXiv:2611.15879), to be
  checked before citing. The journal abstract of `bebesi2019observations` (a conference abstract was read instead).
- **Looked at, no note** (off the question): arXiv:2303.10969 (Dresing et al. 2023), arXiv:2409.03174 (Guo et al.
  2024), arXiv:2007.09194; Do et al. 2025 (ApJ 979, 50; time-dependent acceleration and escape of protons at
  travelling shocks near the Sun; abstract); Burkholder et al. 2025 (10.3847/1538-4357/ada440; foreshock
  fluctuations against Mach number at planetary bow shocks; abstract); Yang et al. 2018 and Qin et al. 2023 (shock
  drift at interplanetary shocks; abstracts). `winner2020evolution` and `morris2022preacceleration` have notes on
  `research/spherical-shock-injection`.
- **Open from the plan's list of unread works:** the plan names "Dresing et al. 2022" (electron intensities against
  the Mach number of coronal shocks). No note was written and this notebook cannot say that it was read. To do at
  step 3.

### Facts from the survey that bear on the hypotheses
Stated without ranking. Sources and ranges are in the notes.
- H1a/H1b: for protons at the PSP event, a steady balance holds only near 1.3 MeV and confinement fails between 2.3
  and 6 MeV; the age limit (~8.9 MeV) and the size limit (~9.5 MeV) lie above (`jebaraj2026diffusive`, preprint).
- H1b: the observed maxima often reach the Hillas limit (`oka2025maximum`); a Bohm scaling with the size matches
  Earth, Saturn and Jupiter if dB/B_0 ≫ 1 and enough time are assumed (`raptis2026relativistic`).
- H2/H3: mass ratio 100 suppresses electron acceleration at M_A 5.3 and inflates it at M_A 21.3 in 1-D parallel shocks
  (`shalaby2024energy`). Runs at M_A ~17-20 with mass ratio 30 show electron tails (`fang2022kinetic`,
  `yu2024nonthermal`).
- H4: seeds before every interplanetary shock (`tsurutani1985acceleration`); > 500 keV at Earth only with 1-5 keV seeds
  (`raptis2025revealing`); seeds to ~100 keV at Saturn (`masters2013electron`).
- H5: at Earth, a model without transients stops at 10-20 keV and with them reaches ≥ 200 keV, at a
  quasi-perpendicular shock (`shi2025compound`); a planar interplanetary shock lacks mature SLAMS
  (`raptis2026compressive`); persistent SLAMS need M_A ≳ 20 (`hemler2025shock`, abstract).
- H6: the 4 electron responses among 51 shocks have θ_Bn 88°, 42°, 80°, 61° (`kartavykh2025statistical`); reconnection
  in the transition energizes electrons in a 2-D run at Earth's parameters (`bessho2023electron`, abstract).
- H7/H8: the PSP crossing is a flank crossing with a saturated low-energy instrument (`dresing2025reason`,
  `kouloumvakos2026wave`); its M_A is 9.1 ± 1.35, 7-7.5, 4.68 ± 3.45 or ~4 depending on the work
  (`jebaraj2024acceleration`, `jebaraj2026wave`, `chiappetta2026evolution`, `giacalone2025source`).
- On the long-term state: "largely insensitive to the details of CR injection" (`galishnikova2026self`, abstract).

### Differences between the README and the approved plan
Each line names the fact first. No decision of the plan (question, hypotheses, runs, proposed criteria) was changed.
- `wilson2016relativistic` is now read (the plan listed it as unread): outside the disturbances THEMIS shows "only
  noise > 12 keV". → The Earth "no" is written as an instrument limit.
- `shi2025compound` confronts an Earth event with a data-constrained transport model; `raptis2026relativistic` gives a
  size scaling over three planets; `hemler2025shock` states the strategy "extend PIC by test particles";
  `grassi2023electron` (the LLNL report of the plan) and `hanusch2020electron` measure or follow electrons in kinetic
  or hybrid fields. → All five are added to the "not new" list.
- `gupta2024electron` varied the piston speed (0.05, 0.1, 0.2 c) for spectra and efficiency; `gupta2025speed`
  measured D_e at one speed. → Gap 2 reads "D̃_e is measured at one shock speed", not "the dependence on v_sh/c is
  unmeasured".
- The statement "> 500 keV only with 1-5 keV seeds" and the Bohm length of 10-100 R_E are in `raptis2025revealing`
  (Nat. Commun. 16, 488), not in arXiv:2411.12815 (`raptis2025multi`). → Citations corrected; the size limit is cited
  to `oka2025maximum` and `raptis2026relativistic`.
- `kartavykh2025statistical`, body: the table has 58 forward shocks, the electron search covers 51, and 4 respond.
  → "4 of 51" replaces "about 8 % of 58".
- `giacalone2025source` (abstract) and `chiappetta2026evolution` give M_A ~4 and a speed of ~1350 km/s for the PSP
  crossing. → The case table carries the spread; fixing the input is part of step 3.
- `bessho2023electron` (abstract): reconnection in the shock transition. → Named under H6 as a further non-DSA
  process and under the limits of the quasi-1-D set-up.
- `jebaraj2026diffusive`: wave growth limits proton confinement at the PSP event. → Written as an open point under
  the question; no candidate was added.
- DOI 10.3847/2041-8213/aeaf56 belongs to `kouloumvakos2026wave` (arXiv:2610.10067).

### Estimates (ours, unverified; not claims)
The p_max formula of `gupta2025speed` with D_1 = D_Bohm and D_2 = D_1/10, compression 4 (the reviewer reproduced the
values). Conventions (1/3, 1/2, the non-relativistic part of the integral) leave a factor 2-4.

| Case | Input (u_1, ω_ci t) | Reached with Bohm | Observed | Reading |
|---|---|---|---|---|
| Earth | 400 km/s, 240 | ~40 keV | none above the noise without transients | inside the factor 2-4: not decided |
| PSP | 2400 km/s, ~10⁵ (4 h, 80 nT) | 6.7 MeV at ω_ci t ≈ 700 | ≥ 6.7 MeV | H1 suffices if D̃ ≲ ~150; the measured proton κ_eff is ~180 u r_L at 2 MeV (our calculation) |
| Saturn | 450 km/s, 30-70 | ~10 keV | ≥ ~1 MeV | two orders short with H1 alone |
| weak shock at 1 au | 150 km/s, ~10⁵ | a few MeV; ~20 keV at 100 × Bohm | mostly 2-10 keV | no conflict if D̃ is two orders larger; needs compression and intensity (H8) |

- The largest energies of the two "yes" cases both lie at p ≈ 1.0 m_i v_sh (Saturn 2.78 against 2.75 m_e c; PSP 14.1
  against 14.7 m_e c). Two readings: the instrument limits happen to lie there; or the stall scale p ≈ m_i v_pt of
  `gupta2024electron`.
- Connection times can change by a factor 2.5-3 with the geometry of the field-line sweep (the reviewer's judgement).

### Next
- Step 3: the input table from primary sources; re-review of the README and the plan; the owner's decision on the
  criteria and thresholds.
