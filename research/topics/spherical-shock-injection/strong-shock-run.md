# A strong shock on a sphere: the heavy-piston run

Topic `spherical-shock-injection`. Written 2026-11-13. Sections 1–6 are committed before the run starts;
Sect. 7 (the observation) is added after it.

## 1. Purpose and standing

- **Instruction.** The owner asked on 2026-11-13: why are the shocks still weak? Start from making a strong
  shock on a sphere. This instruction is the fact that changes the plan: the next run is designed for the
  strength of the shock, not for the length of the drive (S25) and not for the ion-scale run of
  `shock-forming-run.md` (Sect. 6), which stays withdrawn.
- **Standing.** An observation run, like S25. It supports no hypothesis and no claim about electron
  acceleration. The theme (the acceleration mechanism at a shock on a sphere) needs a strong shock on a
  sphere first; this run tries to make one and measures how strong and how long-lived it is.
- **Question of the run.** Can this box hold a spherical front that is a valid shock (protocol Sect. 4) with a
  fast Mach number of 3 or more, and for how long?

## 2. Why the shocks so far are weak

| | L | M | H | S25 |
|---|---|---|---|---|
| nominal M_A; nominal fast Mach | 3.5; 2.95 | 6.0; 4.05 | 8.7; 4.64 | 2.5; 2.27 |
| debris gyroradius ρ_d = M_A d_i [d_e] | 14 | 24 | 34.8 | 10 |
| ambient β | 0.49 | 1.44 | 3.03 | 0.25 |
| R_m / ρ_d | 4.57 | 1.87 | 1.0 | 8.04 |
| half-width of the box / ρ_d | 13.7 | 8.0 | 5.5 | 19.2 |
| fast Mach in the valid bins, 10 / 50 / 90 % | 1.55 / 1.82 / 2.24 | 1.57 / 1.75 / 2.13 | none valid | 1.54 / 1.74 / 2.03 |

- The box is 384 d_e = 96 d_i wide. The scale of the coupling is ρ_d = M_A d_i, so the box is only
  48 / M_A gyroradii from the centre to a face.
- The planner of the example ties the box to the stopping radius: half-width ≥ 1.3 R_B, R_B = R_m M_A^(2/3).
  With the box fixed, R_m / ρ_d = (R_B / d_i) M_A^(−5/3): a higher M_A gives a cloud that is small against its
  own gyroradius, and the debris passes through the ambient plasma instead of pushing it (M, H).
- A lower M_A couples well but is weak by construction (L, S25). So the measured fast Mach number stayed at
  1.5–2.5 whatever the nominal M_A was.

## 3. The lever: a heavy piston in a box smaller than 1.3 R_B

- **Idea.** Give up the stopping of the debris. A cloud with R_m of the order of the half-width of the box
  loses little speed inside the box, so the front stays near V_d, and M_A can be high while R_m / ρ_d stays
  high too. The run ends before the front reaches the faces.
- **Runner.** `run_piston.py` runs the planner of `examples/blast_wave` with its limit on the half-width
  lowered (1.3 → 0.3 R_B). Nothing else changes: same engine, same example, same open faces.
- **V_d is raised, not B lowered.** B₀, the grid, the temperatures and the seed are those of L, so β = 0.49 and
  `L-null` is the null run of this condition. M_A = 5 then means V_d = 0.357 c (L: 0.25 c).
- **What limits the design.**
  - The peak density of the shell is capped by ω_pe dt ≤ 1.2 (planner). With the shell of L (0.33 R_m, width
    0.1 R_m, peak 16.5 n₀) the debris edge (0.1 n₀) starts at 0.65 R_m. A heavier cloud is a larger one.
  - A valid front must stay inside r_outer − V_n / Ω_ci (protocol Sect. 4, criterion e): about
    half-width − 8 d_e − 0.85 ρ_d.
  - Validity has built-in delays in gyro-times: the lead over the debris edge (0.7 V_n / Ω_ci) and a
    persistence of 2 / Ω_ci. In L the front needed 2.2 ρ_d of travel from the initial debris edge.

| candidate | M_A | V_d | nominal fast Mach | R_m [d_e] | R_m / ρ_d | initial edge [d_e] | limit of (e) [d_e] | room [ρ_d] |
|---|---|---|---|---|---|---|---|---|
| L (run) | 3.5 | 0.250 | 2.95 | 64 | 4.57 | 42 | 172 | 9.3 |
| S25 (run) | 2.5 | 0.250 | 2.27 | 80 | 8.04 | 52 | 176 | 12.3 |
| M_A 4, rm 28 | 4.0 | 0.286 | 3.37 | 112 | 7.0 | 73 | 170 | 6.1 |
| M_A 5, rm 30 | 5.0 | 0.357 | 4.21 | 120 | 6.0 | 78 | 167 | 4.5 |
| M_A 5, rm 35 | 5.0 | 0.357 | 4.21 | 140 | 7.0 | 91 | 167 | 3.8 |
| M_A 5, rm 40 | 5.0 | 0.357 | 4.21 | 160 | 8.0 | 104 | 167 | 3.2 |

- The room is the distance from the initial debris edge to the limit of (e). The trade is visible: a heavier
  cloud keeps its speed for longer but starts further out and leaves less room.

## 4. What L and S25 say about the design

The two well-coupled runs, read in units of the debris (`analysis/front.npz`; medians over the polar angles
above 60°, where the front exists from the start):

| | L (R_m / ρ_d = 4.57) | S25 (R_m / ρ_d = 8.04) |
|---|---|---|
| first valid front: t [1/Ω_ci] | 3.30 | 3.24 |
| … radius: [ρ_d]; [R_m]; initial debris edge + [ρ_d] | 5.23; 1.14; + 2.3 | 7.10; 0.88; + 1.9 |
| … front ahead of the debris edge [ρ_d] | 0.66 | 0.72 |
| front speed over V_d at R / R_m = 0.83–0.90 | 0.93 | 0.94 |
| … at 1.0 | 0.85 | 0.89 |
| … at 1.14 | 0.80 | 0.86 |
| … at 1.2 | 0.77 | 0.85 |
| … at 1.37 | 0.72 | 0.79 |
| … at 1.5 | 0.68 | 0.75 |

- **The front is fast before it is valid.** It moves at 0.93–0.95 V_d as soon as it has a speed (L: 2.0 / Ω_ci,
  S25: 2.8 / Ω_ci). It becomes valid at 3.2–3.3 / Ω_ci in both runs, when it leads the debris edge by
  0.7 ρ_d: the lead criterion (b). The radius then is the initial debris edge plus 1.9–2.3 ρ_d.
- **The speed follows R / R_m.** The front keeps 0.85–0.9 V_d out to R_m and 0.7–0.8 V_d at 1.4 R_m, in both
  runs, although M_A and R_m / ρ_d differ. A weak shock was therefore not a failure of the coupling in L and
  S25. Their M_A was low, and the front was valid only late, beyond R_m.
- **What this gives for M_A = 5** (V_d = 0.357, ρ_d = 20 d_e, 1 / Ω_ci = 56, perpendicular fast speed 0.0848,
  nominal fast Mach 4.21; the limit of criterion (e) is at 167 d_e = 8.35 ρ_d):

| | rm 30 | rm 35 | rm 40 |
|---|---|---|---|
| R_m [d_e]; R_m / ρ_d | 120; 6.0 | 140; 7.0 | 160; 8.0 |
| first valid front [d_e]; [R_m] | 116–124; 0.97–1.03 | 129–137; 0.92–0.98 | 142–150; 0.89–0.94 |
| fast Mach there (equator) | 3.6–3.75 | 3.75–3.9 | 3.8–3.95 |
| limit of (e) [R_m]; fast Mach there | 1.39; 3.0–3.3 | 1.19; 3.25–3.6 | 1.04; 3.6–3.7 |
| valid for [1/Ω_ci] | 2.6–3.1 | 1.7–2.2 | 0.9–1.4 |
| front after 2 / Ω_ci of criteria a–c [d_e] | 153 | 167 | 182 |

- The last row is the persistence that criterion (d) asks for. With rm 40 it ends 2 d_e inside r_outer
  (184 d_e): no margin. With rm 30 the fast Mach number falls to about 3 at the end, which is the level the
  run asks about. **rm 35 is chosen**: the fast Mach number stays above 3.2 by either calibration, and the
  valid shock lasts about 2 / Ω_ci.
- **What is extrapolated.** The table assumes that the two regularities above hold at M_A = 5. M (M_A = 6,
  R_m / ρ_d = 1.9) did not couple; no run has both a high M_A and a high R_m / ρ_d. That is what this run tests.

## 5. The run

- **Config.** `configs/strong-P5.yaml`: `explore-L.yaml` with `mach_alfven` 3.5 → 5.0, `rm` 16 → 35 d_i,
  `debris_beams` 4 → 2, and `half_width`, `run` and the output cadences set for this box. The plasma, the grid,
  the cell, the particles per cell, the filter and the seed are those of L.
- **Planner** (`run_piston.py --plan-only` with the launcher on 8 GPUs). B₀ = 0.2857, v_A = 0.0714, fast speed
  0.0848, β = 0.49; V_d = 0.357, nominal fast Mach 4.21; R_m = 140 d_e, ρ_d = 20 d_e, R_B = 409 d_e (outside
  the box); shell at 46.2 d_e, width 14.0 d_e, peak density 16.5 n₀ in 2 beams per species (ω_pe dt ≤ 1.05);
  3.9 × 10⁸ debris particles (L: 0.75 × 10⁸); 768³ cells of 0.5 d_e, dt = 0.25; 68.1 GB per GPU.
- **Length.** 1440 steps, t = 360 = 6.4 / Ω_ci. By Sect. 4 the front is then at the edge of the usable box
  (184 d_e). A shorter run (t = 320) would leave only 0.4 / Ω_ci between the predicted first valid front plus
  the persistence of criterion (d) and the end of the run; this one leaves 1.1 / Ω_ci.
- **Outputs.** Slices and the two box spectra every 36 steps (9 time units = 0.16 / Ω_ci, 41 outputs; L:
  0.33 / Ω_ci); the spectra by position and `hot_e` every 72 steps; tracks every 2 steps as in L. The front
  speed is a slope over five outputs as before, here over 0.64 / Ω_ci (L: 1.3 / Ω_ci). The measured error of
  that slope in L and S25 (median 0.002 c) scales to about 0.004 c: 0.05 in the fast Mach number.
- **Name and place.** Run `P5`; output in `/mnt/nvme/pypic/output/research/spherical-shock-injection/P5`
  (linked from `output/research/spherical-shock-injection/P5`). The command of the protocol (Sect. 8) with
  `run_piston.py` in place of `examples/blast_wave/run.py`; the stop rules of the protocol. Expected: 12 minutes
  on 8 GPUs (S25: 0.44 s per step) and about 40 GB.
- **Null run.** `L-null`: the same plasma, grid, particles and seed without the debris; its outputs up to
  t = 360 give the floor.

## 6. What is observed (fixed before the run)

From `analyze.py front` (protocol Sect. 4, unchanged) and `analyze.py drive`, which now also counts the *high*
bins: fast Mach number 3 or more.

- **O1** first valid front: time in 1/Ω_ci; radius in d_e, R_m and ρ_d; its lead over the debris edge; the
  speed of the debris edge over the speed of the front.
- **O2** fast Mach number in the valid bins: 10 / 50 / 90 % over all of them and per polar angle.
- **O3** valid and high bins: how many, at how many polar angles, the longest duration at one polar angle; the
  same for the bins that are also driven.
- **O4** front speed over V_d against R / R_m (the table of Sect. 4 with a third column), and the exponent m of
  R ∝ t^m per polar angle.
- **O5** what ends the valid bins: criterion (e), the edge of the box, or criterion (c), the weakening.
- **O6** the range of θ_Bn of the valid bins, the compression, and the angle beyond which no de
  Hoffmann–Teller frame exists (cos θ_Bn < V_n / c).
- **O7** upstream electron and ion temperature at the start and the end, the conservation checks (Gauss law,
  div B, energy budget), `fill`. The checks of `examples/blast_wave/analyze.py` that use windows in R_B do not
  apply: R_B lies outside the box.
- **O8** census of the electrons with γ ≥ 1.5 and the tracks, as counts by band. No residual, no origin split
  (`census.SHELL_RADIUS` is the shell of L) and no mechanism is read from them in this document.

**Predictions** (Sect. 4, near the equator; they assume that the two regularities of L and S25 hold at
M_A = 5):

| | predicted | if the coupling fails as in M |
|---|---|---|
| E1 first valid front | t = 2.8–3.4 / Ω_ci (157–190); R = 125–140 d_e | none, or late and slow |
| E2 fast Mach number there | 3.6–3.9 | below 3 |
| E3 valid and high bins, longest at one polar angle | 1.7–2.3 / Ω_ci, ended by criterion (e) | none |
| E4 front speed over V_d at R = R_m; at 1.2 R_m | 0.85–0.89; 0.77–0.85 | below 0.7 at R_m |
| E5 exponent m | 0.8–1 (a driven front; L and S25: 0.53–0.54) | – |

- The range of E1 is wider than 3.24–3.30 / Ω_ci because L and S25 give the first output after the event
  (intervals of 0.33 and 0.46 / Ω_ci).
- **Verdict, fixed now.** "A strong shock on a sphere was made" if the valid and high bins last 1.5 / Ω_ci or
  longer at one polar angle and the median fast Mach number of all valid bins is 3 or more. "Not made" if no
  bin is valid and high, or they last less than 0.5 / Ω_ci. Undecided between.
- **If no bin is valid because the persistence (d) does not fit into the run**, the verdict is "not shown":
  the bins that meet (a)–(c) and (e) are then described with their fast Mach number, and the definition is not
  changed after the run.
- **Side observation.** The time of the first valid front bears on the readings of `shock-forming-run.md`
  (Sect. 5): reading T expects 2.7–3.3 / Ω_ci here too; reading R (a length of 72–73 d_e) expects no delay of
  its own, because the debris edge starts beyond that length.

**Limits, known before the run.**

- One condition, one seed, mass ratio 16, 2 particles per cell. "Strong" means a fast Mach number of 3–4:
  supercritical, far below the Mach numbers of a remnant.
- The shock is young and thin: valid for about 2 / Ω_ci, with a downstream layer of about one debris
  gyroradius between the front and the debris. It is driven all the time; the decay of a free blast wave is
  not in this run.
- V_d = 0.357 c. A de Hoffmann–Teller frame exists only below θ_Bn of about 70°; the planar baseline of the
  protocol (one reflection) applies only there.
- The debris has half the particles per unit of charge of L (2 beams), and the thermal electron gyroradius
  (0.35 d_e) is below the cell (0.5 d_e), as in L.
- The valid domain is a definition. Its lead criterion alone sets a lowest time of about 0.7 r / Ω_ci
  (r: compression), and its persistence needs 2 / Ω_ci of the run after that. So the run shows at most the
  first 2–3 / Ω_ci of a strong shock.
