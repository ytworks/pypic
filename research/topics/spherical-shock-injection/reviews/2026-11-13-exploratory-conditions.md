# Review — exploratory conditions (2026-11-13)

Source: issue #48 of `ytworks/pypic`, "Research (spherical-shock-injection): findings of the review of the
exploratory conditions", opened by the owner on 2026-11-12 (UTC), https://github.com/ytworks/pypic/issues/48.
The issue says that the review itself is not stored in the repository; the text below is the body of the issue,
unchanged, kept as the record. It reviews `exploratory-conditions.md` and `configs/explore-{L,M,H}.yaml` as of
2c33f53. The answers are in `notebook.md` (2026-11-13) and in `exploratory-protocol.md`.

---

Findings of the research review of the three exploratory conditions designed in #47 (`research/topics/spherical-shock-injection/exploratory-conditions.md`, configs `configs/explore-{L,M,H}.yaml`). Nothing has been run. The review itself is not stored in the repository; this issue keeps what has to be acted on.

## A. An inconsistency in the record (to fix)

Sect. 2.2 and 2.3 list one quantity both as matched to a supernova remnant and as not matched.

- The row "(M_A/√β_e) / [(m_i/m_e)^(1/2)/2]" (threshold of `amano2010critical`) is algebraically √2 × (flow speed / electron thermal speed).
- The remnant value 3.2 in that row assumes T_e = 10 eV (3000 km/s). The row "flow speed / electron thermal speed ≈ 13" assumes 10⁴ K (about 0.9 eV, 5000 km/s). With 10⁴ K the first row is about 18 for a remnant, against 3.5 in the runs: not matched.
- To do:
  - [ ] state the speed and electron temperature assumed in each remnant row;
  - [ ] remove this row from "matched" in Sect. 2.3;
  - [ ] reword "θ_Bn runs from 0° to 90°": it holds for the initial geometry, not for a shock, since no shock forms near the poles (Sect. 4.2);
  - [ ] say that, of the remnant's numbers, only β of run L is matched;
  - [ ] apply the same changes to `exploratory-conditions.ja.md`.

## B. Design points raised by the review

1. **The expected signal is not written down.** The flow is only 2.5 electron thermal speeds, so thermal electrons can be mirror-reflected directly, and one reflection above θ_Bn = 56° already exceeds γ = 1.5. The dump `hot_e` will first of all map one-step reflection below the luminal angle, plus the downstream thermal tail. A prediction of reflection for the measured front speed, θ_Bn and compression should be written before the runs; hypotheses then come from the residual.
2. **Numerical artefacts are not excluded.**
   - There is no null run at 2 particles per cell. `spec_up` (5.8 × 10⁶ electrons) cannot see a tail of a few hundred electrons among 9 × 10⁸.
   - Numerical pitch-angle scattering could stand in for the scattering that stochastic shock drift acceleration needs; its rate is not estimated.
   - The heating bound of check 6 (1.25) is looser than its effect on the Maxwellian tail: a rise of T_e by 10–20 % multiplies the fraction beyond V/cos θ_Bn by 2–3 at 45° and by 3–9 at 60° (reviewer's arithmetic).
   - `hot_e` should be split by origin through the ids: debris, ambient from inside the shell, ambient swept by the front, injected at a face.
   - The most useful single addition: run L with `debris: false` and the same outputs (heating map, count of `hot_e`, numerical scattering rate from `tracks_e`).
3. **Check 2 does not show a shock.** It shows a fast compression at the equator. Its criteria were also met in the example's trial in which the debris led the front. Coupling is a risk for M as well as for H (debris gyroradius 0.54 R_m, between the trial's 0.75 and the gate's 0.375). H lasts 0.84 ion gyro-periods, M 1.2. Needed: per ray in `slice_y`, the front ahead of the debris ions at each time; a front definition that also holds where B is not compressed; a domain in (θ_Bn, time) outside which electron data are not read as shock.
4. **The scan is confounded.** The fronts of M and H are expected at nearly the same Mach number (4–5 and at most 5). β, coupling, shock age and the electron gyroradius all change with the nominal M_A. The more informative three runs would be L, the null of L, and M, with H kept as a contrast case.
5. **Smaller points.**
   - Test Σ qE·v Δt against Δγ on `tracks_e` before splitting the work; lower the stride after the first run.
   - Assign θ_Bn by the foot point of the field line on the front.
   - Take the luminal angle from the measured front speed.
   - Keep one energy threshold for all runs.

## C. Decisions for the owner

- [ ] Fix the record as in A?
- [ ] Change the set of runs? Options: replace H by the null of L (this changes the decision "three runs that differ in the field strength only"), or add the null of L as a fourth run.
- [ ] Write the reflection baseline (B.1) and the stricter shock definition (B.3) before the first run?

## Proposed order of the first analysis (from the review)

1. The numerical floor from the null run.
2. Front radius, speed, normal, Mach number, V/v_te and luminal angle per latitude and time; the valid domain.
3. The upstream electron spectrum ahead of the front.
4. The census of `hot_e` by origin and by foot-point θ_Bn, against the reflection baseline.
5. Tracks of the residual: magnetic moment, parallel and perpendicular work, number of encounters, change of the foot-point θ_Bn.

Refs: #47 (conditions), #46 (survey corrections), `research/topics/spherical-shock-injection/reviews/2026-11-12-purpose.md` (first review).
