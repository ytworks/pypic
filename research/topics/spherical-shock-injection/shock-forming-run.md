# A run in which the debris drives the front for longer (S25)

Topic `spherical-shock-injection`. Written 2026-11-13, before the run. Sections 1–6 were committed before the
run started; Sect. 7 (the observation) is added after it.

## 1. Purpose and standing

- **Instruction.** The owner asked on 2026-11-13: first run a condition in which a shock forms, and look at it.
- **Standing.** This is an observation run. It supports no hypothesis and no claim about electron acceleration.
  - `exploratory-results.md` (Sect. 8) recommended another order: read the two unread predecessors, narrow the
    question, then design. The owner chose to run first. The two readings and the narrowing stay conditions
    before any hypothesis.
  - No review was asked for before the run: the front and its valid domain are those of the protocol
    (Sect. 4), unchanged, and nothing is claimed. A review comes before anything is claimed from it.
- **Question of the run.** In L the front became a valid shock only when the debris had almost stopped. Does a
  longer drive give a shock that forms while the debris still pushes?

## 2. What L, M and H showed about the formation

`analyze.py drive` (new, `drive.py`) on the recorded outputs of L and M; H has no valid bin.

| | L | M | H |
|---|---|---|---|
| R_m / ρ_d | 4.57 | 1.87 | 1.0 |
| first valid front: t [1/Ω_ci] | 3.30 | 2.70 | – |
| … its radius [d_e; R_m; d_i; R_B] | 73.2; 1.14; 18.3; 0.50 | 71.6; 1.60; 17.9; 0.48 | – |
| … debris edge [R_m] | 1.00 | 1.34 | – |
| … debris-edge speed over front speed | 0.57 | −0.35 | – |
| valid bins; longest at one polar angle [1/Ω_ci] | 241; 8.59 | 180; 4.62 | 0 |
| valid and driven bins; longest [1/Ω_ci] | 35; 1.32 | 0 | 0 |
| fast Mach in the valid bins, 10 / 50 / 90 % | 1.55 / 1.82 / 2.24 | 1.57 / 1.75 / 2.13 | – |
| valid bins with fast Mach ≥ 2; longest [1/Ω_ci] | 72; 2.31 | 36; 1.16 | – |
| exponent m of R ∝ t^m, median | 0.54 | 0.51 | – |

- **Driven** means: the debris edge (0.1 n₀) moves outwards at half the radial speed of the front or faster.
  Both speeds are slopes over five outputs. The duration of k outputs in a row is k − 1 output intervals.
- **Which criterion delays validity.** Criterion (c), the speed, holds from the first output that has a speed
  (L: t = 111, M: t = 148). Criterion (b), the lead over the debris edge of 0.7 V_n/Ω_ci, is the last to hold.
  - A shock ahead of a steady piston gains the lead V_sh t / r (r: compression). The lead of (b) then needs
    t ≥ 0.7 r / Ω_ci: 1.4–2.1 / Ω_ci for r = 2–3, counted from the time the compression exists.
  - In L and M the lead grew mainly because the debris slowed down. So their shocks became valid as free,
    decaying fronts: in L the debris edge was at R_m and at 0.57 of the front speed, in M it moved inwards.
- **Readings that L and M cannot separate.** The first valid front lies at 3.3 and 2.7 / Ω_ci (a time), at
  73 and 72 d_e = 18 d_i = 0.5 R_B (a length), and at the stalling of the debris (L only).

## 3. The lever

- The debris gyroradius is ρ_d = V_d / Ω_ci = M_A d_i. The piston pushes at nearly V_d only inside about
  0.5–0.7 R_m (snowplough: V_d / (1 + (R/R_m)³)). The drive therefore lasts about R_m / ρ_d gyro-times.
- R_m / ρ_d = (R_B / d_i) M_A^(−5/3). The box of 768³ cells of 0.5 d_e holds R_B ≤ 148 d_e = 37 d_i (the
  planner's lowest half-width, 1.3 R_B). So in this box, at this cell size and mass ratio, **M_A is the only
  lever**, and a longer drive costs shock strength.

| M_A | R_m [d_e] | R_m / ρ_d | β | nominal fast Mach | 1/Ω_ci | run [1/Ω_ci] |
|---|---|---|---|---|---|---|
| 6.0 (M) | 44.8 | 1.9 | 1.44 | 4.05 | 96 | 7.7 |
| 3.5 (L) | 64.0 | 4.6 | 0.49 | 2.95 | 56 | 13.2 |
| 3.0 | 71.2 | 5.9 | 0.36 | 2.63 | 48 | 15.4 |
| **2.5 (S25)** | **80.4** | **8.0** | **0.25** | **2.27** | **40** | **18.5** |
| 2.2 | 87.6 | 9.9 | 0.19 | 2.04 | 35 | 21.0 |

- M_A = 2.5 is chosen. Below it the front is expected to fall under the fast Mach number 1.5 of criterion (c)
  (in L the front moved at 0.79 V_d when it became valid; 0.79 × 2.04 = 1.6 leaves no margin).
- A drive that is both long and strong (R_m / ρ_d ≈ 9 at M_A = 3.5) needs R_B ≈ 295 d_e: cells of 1 d_e in
  this grid. That is a second run (Sect. 6), with numerics that have not been measured.

## 4. The run

- **Config.** `configs/explore-S25.yaml`: `explore-L.yaml` with `mach_alfven` and `wpe_wce` 3.5 → 2.5 and
  `rm` 16.0 → 20.1. Everything else is the same: grid, box, steps, particles, seed, filter, output list.
- **Planner.** B₀ = 0.4, v_A = 0.1, fast speed 0.1099, V_d = 0.25; R_m = 80.4 d_e, R_B = 148.1 d_e, ρ_d = 10 d_e;
  shell at 26.5 d_e, width 8.0 d_e, peak density 16.5 (ω_pe dt ≤ 1.05); 2960 steps, t = 740 = 18.5 / Ω_ci;
  1.48 × 10⁸ debris particles per species (L: 0.75 × 10⁸).
- **Name and place.** Run `S25`; output in `/mnt/nvme/pypic/output/research/spherical-shock-injection/S25`
  (linked from `output/research/spherical-shock-injection/S25`); command and stop rules of the protocol
  (Sect. 8). Expected: 22 minutes on 8 GPUs and 62 GB, as measured on L.
- **No null run** for this condition unless the census shows something to explain. The floor of `L-null` (same
  grid, particles and temperature, weaker field) is quoted as an order of magnitude only.

## 5. What is observed (fixed before the run)

All from `analyze.py front` (protocol Sect. 4, unchanged) and `analyze.py drive`.

- **O1** first valid front: time in 1/Ω_ci; radius in d_e, R_m, ρ_d, d_i, R_B; debris edge in R_m; the speed
  of the debris edge over the speed of the front.
- **O2** valid bins and driven valid bins; their longest duration at one polar angle.
- **O3** fast Mach number in the valid and in the driven bins; how long it stays at 2 or above.
- **O4** exponent m of R ∝ t^m per polar angle (1: free expansion; 2/5: Sedov).
- **O5** upstream electron and ion temperature at the start and the end (numerical heating), conservation
  checks, and the nine checks of `examples/blast_wave/analyze.py`.
- **O6** census of the electrons with γ ≥ 1.5 (the registered threshold) by origin and band, as counts. No
  residual is claimed: this condition has no null run and no second seed.

**Predictions of the three readings of Sect. 2** (outputs are 18.5 apart = 0.46 / Ω_ci):

| reading | first valid front in S25 | driven valid duration |
|---|---|---|
| T: a time, 2.7–3.3 / Ω_ci | t = 108–132; R ≤ 0.8 R_m; debris edge at ≥ 0.7 of the front speed | about 5 / Ω_ci |
| R: a length, 72–73 d_e | t = 148–167 (3.7–4.2 / Ω_ci); R = 0.9 R_m | about 4 / Ω_ci |
| P: the stalling of the debris (edge at R_m, 0.57 of the front speed) | t ≥ 222 (5.5 / Ω_ci) | ≤ 2.5 / Ω_ci |

- The durations come from L's debris edge scaled by R_m (80.4 / 64): it falls below half the front speed at
  about t = 325 = 8.1 / Ω_ci. Under P the driven duration of L, 74 time units, becomes 93 = 2.3 / Ω_ci.
- **Verdict, fixed now.** "A driven shock formed" if the valid and driven bins last 4 / Ω_ci or longer at one
  polar angle; "no better than L" at 2.5 / Ω_ci or less; undecided between. A result that fits none of the
  three readings is reported as such.
- **Expected strength.** With the front at 0.8–0.9 V_d the fast Mach number is 1.8–2.0 next to the equator:
  a weak shock, below the first critical Mach number. Few or no strong bins are expected.

## 6. Limits and the next run

- One condition, one seed, mass ratio 16, 2 particles per cell. β = 0.25 differs from L (0.49), so S25 against
  L changes M_A, β and R_m / ρ_d together.
- The valid domain is a definition (protocol Sect. 4). Its lead criterion alone sets a lowest time of about
  0.7 r / Ω_ci, so reading T cannot be told from that built-in delay by this run alone.
- **Second run (not designed in detail, decided after S25).** The plasma of L (M_A = 3.5) with cells of 1 d_e:
  box 768 d_e, R_m = 128 d_e, R_m / ρ_d = 9.1, nominal fast Mach 2.95. It needs a wider shell (ω_pe dt ≤ 1.2
  at the peak), has a Debye length of 0.1 cells (L: 0.2) and does not resolve electron scales: an ion-scale
  observation whose upstream heating must be measured in the run itself.

## 7. The observation (added 2026-11-13, after the run)

Run `S25` from commit `66a3629`: 2960 steps in 1310 s on 8 GPUs, 60 GB. Record: `runs/` (name `S25`).
Analysis: `front`, `drive`, `census`, `tracks` and the nine checks of the example.

**O5, the run itself.** The nine checks pass (Gauss law 6.8 × 10⁻¹⁵, div B 2.5 × 10⁻¹⁴, `fill` 0.68, energy
budget +9.1 % of the initial energy, compression 3.0). Upstream T_e rises from 0.0100 to 0.0110 (10 %), T_i
stays at 0.0099.

**The front and its valid domain.** 738 bins (41 outputs × 18 polar angles); a front in 387, valid in 214,
all at polar angles of 37.5° to 87.5° (none nearer to the axis). Medians per polar angle in the valid bins:
V_n = 0.180–0.199 c, fast Mach number 1.71–1.81, compression of B 1.8–2.1 (bin) and 2.0–2.3 (cell), front ahead
of the debris edge by 12–33 d_e. No de Hoffmann–Teller frame beyond θ_Bn = 76–81°.

| | T: a time | R: a length | P: stalling | observed |
|---|---|---|---|---|
| O1 first valid front, t [1/Ω_ci] | 2.7–3.3 | 3.7–4.2 | ≥ 5.5 | **3.24** (output before: 2.78) |
| … radius | ≤ 0.8 R_m | 0.9 R_m (72–73 d_e) | – | **0.88 R_m** (71.0 d_e; output before: 66 d_e) |
| … debris edge: place; speed over front speed | –; ≥ 0.7 | – | R_m; 0.57 | **0.79 R_m; 0.66** |
| O2 valid and driven, longest [1/Ω_ci] | about 5 | about 4 | ≤ 2.5 | **2.78** |

- **O1.** At 6 polar angles at once. In other units: 7.1 ρ_d, 17.8 d_i, 0.48 R_B. Fast Mach number 2.1 there.
- **O2.** Valid: 214 bins at 11 polar angles, the longest 11.1 / Ω_ci (87.5°). Valid and driven: 66 bins at 11
  polar angles, the longest 2.78 / Ω_ci (42.5°), median 2.3. The speed of the debris edge over that of the
  front falls steadily: 0.66 at 3.2 / Ω_ci, 0.50 at 6.5, 0.34 at 8.3, below 0.15 after 11.5.
- **O3.** Fast Mach number, 10 / 50 / 90 %: 1.54 / 1.74 / 2.03 in the valid bins, 1.80 / 1.99 / 2.09 in the
  driven ones. At 2 or above: 29 bins at 8 polar angles, the longest 2.3 / Ω_ci; all of them are driven. None
  at 3 or above. Median over the bins with a front: 2.15 at 2.8 / Ω_ci, 1.99 at 4.6, 1.75 at 8.3, 1.39 at 17.6.
- **O4.** m = 0.53–0.55 at the 11 polar angles (median 0.536); L: 0.54.
- **What ends the valid bins.** Criterion (c) at every polar angle: the fast Mach number falls below 1.5,
  between 9.3 / Ω_ci (37.5°, R = 115 d_e) and 14.3 / Ω_ci (87.5°, R = 159 d_e). The box does not end them.

**Verdict by the rule of Sect. 5: undecided.** The valid and driven bins last 2.78 / Ω_ci, between 2.5 and 4.
S25 doubled the driven time of L (1.32 / Ω_ci) and did not reach "a driven shock formed".

**The three readings.**

- **P is rejected.** The front is valid at 3.24 / Ω_ci, while the debris edge is at 0.79 R_m and still moves at
  0.66 of the front speed. Validity does not wait for the stalling.
- **T and R are not separated.** The time fits T (3.24 against 3.30 in L); T's other predictions are not all
  met (radius 0.88 against ≤ 0.8 R_m; edge speed 0.66 against ≥ 0.7; driven 2.8 against about 5 / Ω_ci). The
  radius, 71 d_e, is again the length of R (72–73 d_e), but R's time (3.7–4.2 / Ω_ci) is not met. With outputs
  0.46 / Ω_ci apart the two cannot be told apart in this run.
- **Seen beside the readings** (used in `strong-shock-run.md` Sect. 4): in units of the debris, L and S25
  agree. The front becomes valid 1.9–2.3 ρ_d outside the initial debris edge, when it leads the edge by
  0.7 ρ_d, and its speed over V_d follows R / R_m. This is a description of two runs, not a tested rule.

**Strength, as expected:** weak. The fast Mach number is 2.1 at the first valid front and falls from there.

**O6, the census (counts; no residual is claimed).** 1,234,713 electron records with γ ≥ 1.25 in 41 outputs.

- **Ambient electrons with γ ≥ 1.5.** 57 at t = 37 and 533 at 55.5, before any front is found; 3245–3506 at
  t = 92.5–129.5 (2.3–3.2 / Ω_ci, the time up to the first valid front); 1231 at 222 (5.6 / Ω_ci); 109 at 370;
  11–33 after 444 (11 / Ω_ci). Nearly all are in the layer or downstream and are classed "swept by the front".
- **Upstream of the front** there are at most 13 of them in one output, against the registered signal
  threshold of 100: no band is measurable. Summed over the outputs (one electron counts once per output):
  9 "S0 at exit", 10 undetermined and 3 residual at a detected front, of which 7, 2 and 2 are in the valid
  domain; 116 undetermined where no front is found.
- **Debris electrons with γ ≥ 1.5** (weights; they are the driver, not the ambient plasma): 1337 at t = 74,
  2530 at 129.5, about 1000 at 370–407 and 3185 at the end; 70–570 of them are upstream of the front from
  t = 148 on.
- **Limits of the census here.** No null run and no second seed for this condition. `census.SHELL_RADIUS`
  (27.5 d_e) is the shell edge of L; in S25 the shell ends at 34.6 d_e, so ambient electrons born between the
  two are counted as "swept by the front". The numerical scattering rate is taken from this run
  (ν = 0.0029 at γ = 1.25–1.5).

