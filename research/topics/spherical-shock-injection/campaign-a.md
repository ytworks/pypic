# Campaign A: what sets where on the sphere the shock reflects and energises electrons

Designed 2026-11-13, before the runs that need the predictions (owner: "first A"; use the machine fully, analyse
and delete what can be re-run). Configs: `configs/campaign_a.py`. Queue: `run_queue_guarded.sh`; after each run
`post_run.sh` (analyses, then pruning of tracks and particle dumps; the record keeps commit, config and seed).

## 1. Background (P5, notebook 2026-11-13)

- P5 (M_A 5, v_te 0.10, m_i/m_e 16) made a strong shock on a sphere (fast Mach 3.0–3.6). Electrons that cross it
  once are heated, not accelerated beyond the heated distribution in the plasma frame.
- The shock reflects electrons back upstream at high latitude: returning with γ ≥ 1.25 for 1–4 % of the entries
  at polar angles 10–50° (θ_Bn 15–50°, V_HT 0.27–0.37 c), about none above 60°; energies at the level of one
  adiabatic reflection off a mirror moving at V_HT (γ_out over γ_HT²(1 + V_HT²): median 1.05).
- **Stage 0** (`P5/analysis/stage0.json`, `stage0_pref.json`; 261,352 entries of electrons born outside the
  initial debris shell): the probability that an entry returns upstream with a gain of more than 2 %,
  **P_ref**, falls with V_HT: 0.17 at 0.20–0.25 c, 0.022 at 0.4–0.5 c, 0.0008 at 0.8–1.0 c. Against x = V_HT / v
  (v the electron's own speed) it falls smoothly: 0.18 at x = 0.5–1, 0.026 at 2–2.5, 0.0039 at 4–5, 0.0003
  above 20.

## 2. Question

Is the latitude up to which the shock reflects and energises electrons set by V_HT against the electron speed
(so that hotter upstream electrons are reflected closer to the luminal angle and reach higher energies), and does
the Alfvén Mach number change it at a fixed ratio?

- **H_T:** P_ref depends on x = V_HT / v only.
- **H_0:** P_ref depends on V_HT only (temperature does not matter).
- **H_M:** at a fixed x, P_ref changes with M_A (A7 against A5).

## 3. Runs (all 8 GPUs; T_i = 0.01 m_e c² everywhere; seeds 12345 and 54321)

| set | M_A | v_te | cells [d_e] | box [d_e] | steps | β | nominal fast Mach |
|---|---|---|---|---|---|---|---|
| A5 (P5's box, field, debris) | 5 | 0.10 (P5, + seed 2), 0.14, 0.20, 0.28 | 0.5 | 384 | 1440 | 0.49–2.17 | 4.21–2.99 |
| A7 | 7.5 | 0.14, 0.20, 0.28 | 0.75 | 576 | 1600 | 1.63–4.87 | 4.88–3.33 |

Null runs (no debris, same plasma, field and grid): L-null for A5 at v_te 0.10; one per other temperature and
set. A7 uses R_m = 45 d_i = 6 ρ_d and a shell of width 0.2 R_m (ω_pe dt at its peak). The Debye length is
0.19–0.56 cells (P5: 0.20).

## 4. Measures (fixed now; `analyze.py campaign`)

- Entries: steps of an electron born outside the initial debris shell from upstream into the transition layer
  (sites of `sites.md`); outcome: back upstream or into the steep wall; γ_out read 8 samples after leaving.
  V_HT = V_n / cos θ_Bn with V_n the front speed at the polar angle and θ_Bn local (`sites.local_tbn`).
- **P_ref(V_HT)** in the bins 0.20 / 0.25 / 0.30 / 0.35 / 0.40 / 0.50 / 0.60 / 0.80 / 1.0 / 1.5 / 3.0 c, and
  **V_c**, the V_HT at which P_ref falls through 1 % (linear between bin centres).
- **P_ref(x)**, x = V_HT / v, in the bins of stage 0.
- Energy: γ_out of the returns (p50, p90, p99, maximum) and their share with γ_out ≥ 1.25 by polar angle.

## 5. Predictions (fixed now)

P_ref per V_HT bin (centres 0.225 … 2.25 c); H_0 is P5 as measured; H_T rescales the momenta of P5's entering
electrons by v_te / 0.1 and applies P5's P_ref(x):

| V_HT bin [c] | 0.20–0.25 | 0.25–0.30 | 0.30–0.35 | 0.35–0.40 | 0.40–0.50 | 0.50–0.60 | 0.60–0.80 | 0.80–1.0 | 1.0–1.5 | 1.5–3.0 |
|---|---|---|---|---|---|---|---|---|---|---|
| H_0 (P5) | 0.172 | 0.146 | 0.081 | 0.042 | 0.022 | 0.0078 | 0.0026 | 0.0008 | 0.0006 | 0.0005 |
| H_T, v_te 0.14 | 0.123 | 0.121 | 0.098 | 0.077 | 0.056 | 0.034 | 0.018 | 0.0083 | 0.0041 | 0.0016 |
| H_T, v_te 0.20 | 0.110 | 0.131 | 0.127 | 0.114 | 0.094 | 0.066 | 0.039 | 0.019 | 0.0088 | 0.0028 |
| H_T, v_te 0.28 | 0.087 | 0.118 | 0.140 | 0.142 | 0.128 | 0.101 | 0.070 | 0.040 | 0.018 | 0.0048 |

V_c: H_0 0.53 c (θ_Bn ≈ 59° at V_n = 0.27 c); H_T 0.86 c (0.14; θ_Bn ≈ 72°), beyond c for 0.20 and 0.28 (then
up to the luminal angle). Energy: under H_T the hotter runs return electrons where V_HT is larger, so p99 of
γ_out of the returns rises with v_te; under H_0 it does not.

## 6. Verdicts (fixed now)

- **Per A5 temperature (0.14, 0.20, 0.28), each seed:** misfit M = mean over bins with ≥ 2000 entries of
  [ln((P_meas + 0.002) / (P_pred + 0.002))]². "H_T favoured" if M_T < M_0 / 2 in both seeds; "H_0 favoured" if
  M_0 < M_T / 2 in both; else undecided. The campaign verdict on temperature: H_T if favoured at two of the three
  temperatures and H_0 at none; H_0 likewise; else undecided.
- **Seed check:** P5 against A5-T10-s2: P_ref per bin within a factor 1.5 in the bins with ≥ 2000 entries. If
  not, the misfit criterion above is not used and the result is reported as not reproducible at this size.
- **M_A (H_M):** at the same v_te, P_ref(x) of A7 against A5 in the x bins with ≥ 2000 entries in both: "M_A
  matters" if the ratio is outside 0.5–2 in at least half of those bins in both seeds; "no M_A effect" if inside
  in all; else undecided.
- **Nulls:** a null run with any electron of γ ≥ 1.5 upstream at more than 10 per output, or a P_ref-like count
  (none expected: there is no front), is reported, and the run of that temperature is read only where its
  counts exceed the null by 10 times.

## 7. Limits, known now

- The Buneman wave (shock surfing) has a wavelength of 3–4 cells in every set and is damped by the filter: that
  mechanism is not tested.
- Mass ratio 16; one geometry per set. Higher β lowers the fast Mach number (2.99 at A5-T28).
- The entry outcome uses site boundaries defined after seeing P5 (`sites.md`); the same definition is used in
  every run unchanged.
- Electrons placed inside the debris shell at t = 0 are excluded by their birth radius; an engine option to keep
  ambient plasma out of the shell is planned on main (B-1), not used here.
- Outputs are pruned after the analyses; a re-analysis needing the tracks needs a re-run.

## 8. Amendment after the design review (2026-11-13, before the runs that need it)

Review: `research-reviewer` on this document (advice, not approval). Taken: the return floor from boundary
flicker, the missing mirror ratio, the confounded M_A comparison, the energy measure, an absolute fit scale.
The configs do not change; only the measures and verdicts below replace those of Sects. 4–6 where they differ.

- **Two definitions of a reflection**, both scored; a verdict must hold under both:
  - *2 %*: the entry returns upstream with γ_out > 1.02 γ_in (Sect. 4);
  - *persistent*: the same, and after the return the electron stays upstream for one upstream gyro-period
    2π γ_in / Ω_ce0 or gets ahead of the front by d_i + 2 ρ_e (ρ_e = γ v / Ω_ce0).
- **Bins:** the verdicts use the V_HT bins below 1.0 c only. Returns beyond the luminal angle are reported apart.
- **H_A, an adiabatic mirror in the de Hoffmann–Teller frame without potential**, per entry from the run's own
  u_par, u_perp, V_HT and the compression B_max / B_1 of `front.npz` (`campaign.adiabatic_reflection`); no
  fitted constant. In P5 (persistent reflections against H_A): 0.066 / 0.25 at V_HT 0.20–0.25, 0.0052 / 0.066 at
  0.40–0.50, 0.0015 / 0.014 at 0.50–0.60, 0 / 0 above 0.8 c: the same cut-off, a level 4–12 times lower, and a
  ratio that falls with V_HT. By polar angle H_A extends to 60° where the measurement stops near 40°.
  **Verdict on H_A:** in each run and seed, the ratio measured / H_A per bin (bins with ≥ 2000 entries and
  H_A ≥ 0.005) against the ratio of P5 in that bin: "H_A carries the temperature trend" if within a factor 2 in at
  least two thirds of those bins in both seeds; "H_A misses a temperature-dependent factor" if outside in at
  least two thirds; else undecided.
- **Absolute scale:** "neither H_T nor H_0 fits" when both misfits exceed 4 times the misfit of P5 against
  A5-T10-s2 under the same definition (computed from those two runs before any other run is read).
- **Energy:** γ_out / γ_in of the persistent returns (p50, p90) and γ_out over the bound γ_HT²(1 + V_HT²) γ_in;
  against H_A's own γ_out for the same entries. Not γ_out alone.
- **A7 against A5:** reported as a comparison of two sets (β, ω_pe / Ω_ce, cell size and R_m differ too), not as
  an effect of M_A. Partial control at similar β: A7-T14 (β 1.63) against A5-T28 (β 2.17).
- **Counts:** entries per V_HT bin and per polar angle (`n_by_VHT_and_polar`) are reported with every result.

Predictions, persistent definition (P5 measured for H_0; P5's P(x) applied to rescaled momenta for H_T;
`P5/analysis/campaign_predictions.json` holds both definitions):

| V_HT bin [c] | 0.20–0.25 | 0.25–0.30 | 0.30–0.35 | 0.35–0.40 | 0.40–0.50 | 0.50–0.60 | 0.60–0.80 | 0.80–1.0 |
|---|---|---|---|---|---|---|---|---|
| H_0 | 0.066 | 0.054 | 0.027 | 0.014 | 0.0052 | 0.0015 | 0.0003 | 0 |
| H_T, v_te 0.14 | 0.038 | 0.037 | 0.028 | 0.021 | 0.015 | 0.010 | 0.0060 | 0.0035 |
| H_T, v_te 0.20 | 0.034 | 0.041 | 0.039 | 0.034 | 0.026 | 0.018 | 0.011 | 0.0066 |
| H_T, v_te 0.28 | 0.025 | 0.037 | 0.045 | 0.045 | 0.039 | 0.028 | 0.018 | 0.011 |

## 9. Seed check and absolute scale (2026-11-13, from P5 and A5-T10-s2, before any other run is read)

- P_ref of A5-T10-s2 over P5, bins 0.20–1.0 c: 2 % definition 0.93 / 0.92 / 1.04 / 0.99 / 1.22 / 1.25 / 1.04 /
  0.31; persistent 0.86 / 0.80 / 1.05 / 0.98 / 1.04 / 1.05 / 0.98 / (0 against 0.0001). Every bin has ≥ 2000
  entries, but the last one holds 16 against 6 reflections (2 %) and 0 against 2 (persistent).
- **Amendment of the seed check (Sect. 6):** the factor 1.5 applies to bins with at least 50 reflections in both
  runs; the rule as first written failed on a bin of a few events only. With it the check passes under both
  definitions. Recorded before any other run is read.
- **Misfit of the seeds** (Sect. 6 formula, all eight bins): 0.016 (2 %), 0.0087 (persistent). "Neither fits"
  therefore means both misfits above 0.065 (2 %) or 0.035 (persistent). Numbers in `P5/analysis/seed_misfit.json`.
