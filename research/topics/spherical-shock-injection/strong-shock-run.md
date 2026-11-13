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
