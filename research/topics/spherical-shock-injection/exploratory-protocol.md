# Exploratory runs: protocol fixed before the runs

Topic `spherical-shock-injection`. Written 2026-11-13, before any run of L, M or H.
**Status: fixed before the first run; the numbers marked "proposal" may change only through the review R1
(Sect. 11), and only before the first run.**

This file answers the review of the exploratory conditions (issue #48, `reviews/2026-11-13-exploratory-conditions.md`).
It replaces Sect. 6 (checks 2 and 6) and Sect. 8 of `exploratory-conditions.md`. The physical conditions of
L, M and H are unchanged. Code: `shockkin.py` (kinematics, count model), `analyze.py` (commands).
Units: c = 1, m_e = 1, lengths in d_e, times in 1/ω_pe.

## 1. The question of this exercise and its three possible answers

- **Question (from the owner, 2026-11-13).** Do runs L, M, H, analysed, give a research hypothesis about the
  electron acceleration mechanism that is novel and certain enough to build a study on?
- **A hypothesis is built only from a residual:** what is left of the energetic electrons after (i) the
  numerical floor of a null run, (ii) the known one-step reflection (baseline S0, Sect. 3), and (iii) the
  regions that hold no shock (Sect. 4) are taken away.
- **The answer is one of three** (criteria in Sect. 8):
  - **obtained:** all of C1–C6 hold;
  - **candidate:** a residual with a positive signature exists, and the missing evidence is named;
  - **not obtained from these runs:** with the conditions that a run would need.
- **Honest expectation.** One mass ratio (16), 2 particles per cell and V/c ≈ 0.2: the best outcome is a
  hypothesis with its test designed, not a result.

## 2. Facts measured in the 16-step check of L (run 0, 2026-11-13)

The check ran `configs/check-L.yaml` on 8 GPUs at commit `3e97c6c`; outputs in
`output/research/spherical-shock-injection/check-L` (2.9 GB, not in git).

- All nine diagnostics wrote files. Gauss law 4 × 10⁻¹⁵, div B ≤ 1.3 × 10⁻¹⁶, fill 0.667.
- **Cost:** 0.90–1.06 s per step, so a run of 2960 steps takes about 45–50 minutes, not the 15–20 minutes of the
  conditions record (Sect. 7 there was an estimate). The whole run set is about 3–5 GPU-hours.
- **Selection of the dump works:** with `gamma_min: 1.05` the smallest γ in the dump is 1.05, and the dump
  holds 17,782,723 ambient electrons at step 0 (a fraction 0.019628).
- **The loaded electrons are Maxwell–Jüttner** with T = v_te² = 0.01 m_e c², f ∝ exp(−(γ−1)/T): the measured
  tail fractions agree with it to 0.03 % at γ ≥ 1.05 and to 3 % at γ ≥ 1.15, while a Gaussian in u gives
  0.0166 at γ ≥ 1.05. The baseline uses the Maxwell–Jüttner form.
- **Ids.** Ambient electrons have ids below 2N (N = 768³ cells); `id // 2` is the linear index of the birth
  cell (checked: every position lies in its cell). Electrons injected at the open faces have ids ≥ 20N.
  Debris electrons are the species `de`; their weights are small and vary, so they are counted by weight.
- **Tracks:** 191,416 ambient electrons, momentum half a step before position and fields.

## 3. The known signal: one-step reflection (baseline S0)

### 3.1 Physics
The flow is only about 2.2 electron thermal speeds, so thermal electrons can be mirror-reflected by the moving
compression of the field (the "fast Fermi" reflection of the bow-shock literature; shock drift acceleration in
the de Hoffmann–Teller frame). This is known physics. No novelty is claimed for it.

- The upstream is at rest with E = 0 and B = B₀ ẑ. A field line is a line of constant distance ρ from the axis.
- For an axisymmetric front R(ϑ, t) (ϑ from the axis), the foot point of a field line moves along it with
  **V_HT = R_t R / (R_ϑ sin ϑ + R cos ϑ)**, and cos θ_Bn = |R_ϑ sin ϑ + R cos ϑ| / (R² + R_ϑ²)^(1/2). For a sphere
  this is V_sh / cos ϑ and θ_Bn = ϑ.
- In the frame that moves with the foot point, an electron keeps its energy and its magnetic moment. With
  g = (1 − V_HT²)^(−1/2), u′ = g (u_∥ − V_HT γ), γ′ = g (γ − V_HT u_∥):
  - it meets the front if u_∥/γ < V_HT;
  - it is reflected if (γ′ + φ)² − 1 − b u_⊥² < 0, with b = B_max/B₀ and φ ≥ 0 the rise of the potential in
    that frame (in m_e c²/e);
  - after the reflection u_∥ = g (−u′ + V_HT γ′), γ = g (γ′ − V_HT u′). A cold electron gets
    γ = (1 + V_HT²)/(1 − V_HT²).
- Where V_HT ≥ 1 (the superluminal band, cos θ_Bn ≤ V_n) no reflection is possible.
- The energy after a reflection depends on V_HT and on the momentum before it, not on b and φ. b and φ decide
  only **how many** electrons are reflected. So the energies are a firm prediction and the counts are not.

### 3.2 Two uses
- **B1, per electron (the primary test).** For an energetic electron found ahead of the front:
  1. take its guiding centre (r + u × B / (q B²)) and its parallel velocity;
  2. follow it back along its field line as a free particle to the last time it was at the foot point;
  3. take V_HT and θ_Bn of the fitted front there, and undo the reflection (the map is its own inverse).

  The electron **is explained by S0** if a foot-point crossing exists, its kinetic energy before the reflection
  is at most E_th, and it was outside the loss cone for b_max and φ = 0 (the most permissive case).
  - E_th is the energy above which the Maxwell–Jüttner distribution at the measured upstream temperature holds
    fewer than 10 macro-electrons in the whole box: 0.20 m_e c² at T = 0.010, 0.245 at T = 0.012.
  - b_max is the largest measured compression over the detected front (bin averages, Sect. 4). A second pass
    uses the largest single-cell value, as a test of how much a rippled front could change the answer.
  - The distribution of the energies before the reflection is also compared with the flux-weighted thermal
    distribution, as a consistency check of the whole procedure.
- **B2, counts (a consistency check).** The number of reflected electrons per band of foot-point θ_Bn (10°) and
  per energy bin (γ 1.5–2, 2–3, 3–5, > 5) that are still in the box at each dump, from the flux of a
  Maxwell–Jüttner upstream on the fitted front.
  - Upper bound: b(ϑ, t) as measured, φ = 0, T_e as measured upstream.
  - Lower bound: the same with the initial T_e and b one standard error lower; with a measured φ if its error
    is under 50 %, otherwise the lower bound carries the note "φ not measured".
  - B2 is computed for the whole detected front and, separately, for the valid shock domain (Sect. 4). A moving
    compression reflects electrons whether or not it is a shock.
  - If the measured count of S0-explained electrons falls outside [lower/2, 2 × upper] in the bands that hold at
    least 100 electrons, the baseline or the front measurement is wrong, and no residual is discussed until
    this is understood.

### 3.3 Prediction before the runs
A nominal sphere R = 21.12 + 0.2175 t (the shell radius, and 0.87 V_d from the example's measured front), a
Maxwell–Jüttner upstream at T = 0.01, 16 macro-electrons per d_e³, box half-width 192 d_e, φ = 0
(`analyze.py predict`, numbers in `data/predict-nominal.json`). The luminal angle is 77.4°.
The nominal prediction is the same for L, M and H: they share V_d, v_te, the shell and the box.

Macro-electrons with γ ≥ 1.5, both hemispheres, at the last dump (t = 740):

| b | reflected in the run | still in the box | of these γ 1.5–2 | γ 2–3 | γ ≥ 3 | bands that hold 99 % |
|---|---|---|---|---|---|---|
| 2 | 1.0 × 10³ | 5.3 × 10² | 5.3 × 10² | 0 | 0 | 50–70° |
| 3 | 6.7 × 10⁴ | 3.3 × 10⁴ | 3.2 × 10⁴ | 69 | 0 | 50–70° |
| 4 | 4.0 × 10⁵ | 1.9 × 10⁵ | 1.9 × 10⁵ | 1.8 × 10³ | 0 | 50–70° |

What this fixes before the runs:
- **S0 gives γ ≥ 1.5 almost only for foot-point θ_Bn between 50° and 70°, and nothing at γ ≥ 3.** The cold
  limit alone would give 61°; the electrons below it are those with a large perpendicular momentum.
- **The count changes by a factor of about 400 between b = 2 and b = 4.** A count that exceeds S0 is therefore
  not evidence by itself; the per-electron test B1 is the primary one.
- The real front slows and is not a sphere; the prediction is recomputed with the measured front. The nominal
  table stays as the record of what was expected.
- **The null run:** a Maxwell–Jüttner plasma at T = 0.01 holds 2 × 10⁻¹² macro-electrons with γ ≥ 1.5 in the
  whole box, so every such electron in the null run is numerical or comes from the open faces.

### 3.4 Limits of the baseline
- The back-tracing assumes a uniform upstream field and free flight. It is valid only where the upstream is
  undisturbed (check 4 of the conditions record) and for flight times below 0.1/ν_num (Sect. 5, N3).
- The front is treated as axisymmetric and smooth on the scale of the bins. Reflection by ripples smaller than
  a bin counts as "beyond S0" here; the single-cell pass of B1 bounds this.
- The potential φ is hard to measure at 2 particles per cell. φ = 0 gives the largest reflected count, so the
  upper bound does not need it.

## 4. What counts as a shock: the front and its valid domain

All from `slice_y` (the meridional plane through the centre; it contains B₀), every 74 steps. Polar bins of
2 d_e in radius and 5° in polar angle about the centre; the two sides in x and the two hemispheres are four
samples of one front and are compared before they are folded. d_i = 4 d_e in all runs.

- **Front position R(ϑ, t).** The outermost radius at which the density of the ambient ions `n:i`, averaged over
  d_i in radius, is at least 1.5 n₀ in two neighbouring bins; the crossing is interpolated. This also works
  where B is not compressed (the quasi-parallel side). At the equator it is compared with the example's
  definition (B_z ≥ 1.25 B₀).
- **Front model.** A least-squares fit of the positions with {1, cos 2ϑ, cos 4ϑ} × {1, τ, τ², τ³}; V_HT, θ_Bn
  and the normal speed V_n follow analytically. The rms residual of the fit is reported; it has to be below
  d_i for B1 to be used.
- **Compression b(ϑ, t).** The largest bin average of |B|/B₀ between 3 d_i behind and 1 d_i ahead of the front.
- **Valid shock domain (proposal).** A bin (ϑ, t) belongs to it if all of these hold:
  - (a) a front position exists;
  - (b) the front is at least 2 d_i ahead of the debris ions (outermost radius with `n:di` ≥ 0.1 n₀);
  - (c) the normal speed from five neighbouring outputs is at least 1.5 times the fast speed for the local
    θ_Bn (L: 0.071–0.085 c; M: 0.046–0.062 c; H: 0.046–0.054 c);
  - (d) (a)–(c) hold in five outputs in a row;
  - (e) the front is between 0.3 R_B (44 d_e) and 172 d_e (5 d_i from the open face).
- **Outside the valid domain** the electron data are reported as "fast compression, no shock". They are not
  read as shock acceleration, in any run.
- **Luminal angle.** From the measured front: the foot points with V_HT ≥ 1. The band is marked in every
  table and figure.
- **θ_Bn of an electron** is the θ_Bn at the foot point of its field line: at the time of its reflection for
  B1, and at the time of the dump for the census.

## 5. Numerical floor: the null run of L

`configs/explore-L-null.yaml`: the same plasma, field, grid, seed, length and outputs as L, without the debris.

- **N1, counts.** Electrons with γ ≥ 1.5 per dump, by origin (ambient, injected at a face) and by position.
  The expectation is zero. A count in a run is a **signal** only if it is at least 100 and at least 10 times
  the count of the null run in the same class.
- **N2, heating.** T_e(x, z, t) from `spec_meridional` and T_e(t) from `spec_up` (fit to the bulk). The old
  bound of check 6 (1.25) is dropped: the baseline is computed with the initial and with the measured
  temperature, and a band in which the two differ by more than the residual is "explained by heating".
- **N3, numerical scattering.** The pitch-angle scattering rate ν_num from the tracks in the uniform field,
  ν = ⟨(Δμ)²/(1 − μ²)⟩/Δt at lags of 4, 16 and 64 time units, for all tracked electrons and for the fastest
  tenth. Uses:
  - back-tracing (B1) only for flight times below 0.1/ν_num;
  - a mechanism that needs scattering (stochastic shock drift acceleration) is not assigned unless
    ν_num τ < 0.1, with τ the crossing time of the front layer (its measured width over V_n).
- **N4, closure of the work.** On the tracks, Σ q E·v Δt against Δγ (the momentum is recorded half a step
  before the fields, so the velocity is centred first). The split into parallel and perpendicular work is used
  only if the median error is below 10 % for the electrons that gained Δγ ≥ 0.1. Otherwise orbits are
  classified by the magnetic moment and the number of encounters only.
  - In the null run no electron gains that much, so N4 is tested on run L itself, on the tracked electrons
    behind the front, before any split is read.
- A null run of M is made only if M has a valid domain.

## 6. The census of the energetic electrons

One energy threshold for all runs: **γ ≥ 1.5** (50 kT_e; the ram energy of an ion is 0.5 m_e c² at V_d).
Every electron of the dump `hot_e` (every 148 steps) gets:

- **origin**, from the id: debris (species `de`, counted by weight); injected at a face (id ≥ 20N); ambient
  born inside the shell (birth radius ≤ 27.5 d_e); ambient swept (the front has passed its birth cell);
  ambient not yet reached;
- **position class**, from its guiding centre and the fitted front: upstream (more than 1 d_i ahead of the
  front, or on a field line that the front has not reached), front layer (within 1 d_i ahead and 3 d_i behind),
  downstream; and whether the foot point is in the valid domain;
- **foot-point θ_Bn** (10° bands) and the luminal flag;
- for upstream electrons: the result of B1;
- **history across dumps**, by matching ids: the dump in which it first appears, and its class there.

Tables: counts per dump × origin × position class × θ_Bn band × energy bin, with the null-run count beside
each. A band with fewer than 100 electrons is "not measurable" and is not used for a hypothesis.

Downstream, the reference is the thermal tail: the count above γ = 1.5 expected from a Maxwell–Jüttner fit to
the bulk of the local spectrum (`spec_meridional`). Downstream electrons enter a hypothesis only with tracked
orbits or id histories that show how they gained their energy.

## 7. Candidate signatures, looked for in the residual only

The residual: ambient electrons with γ ≥ 1.5 in the valid domain that B1 does not explain (upstream), or that
exceed the thermal tail (downstream), above the floor N1.

- **S0** one-step reflection: known; the baseline.
- **K1, a path across latitudes.** An electron gains its energy at a high θ_Bn and leaves upstream along the
  field at a lower one: the foot-point θ_Bn at the gain and at the exit differ by more than 20°, and the energy
  exceeds the S0 maximum of the exit latitude. In the ideal axisymmetric picture the drifts are azimuthal, so
  the latitude changes only by motion along the draped field behind the front.
- **K2, several encounters.** Two or more crossings of the front, and an energy above one reflection.
- **K3, gain in the superluminal band,** where reflection is impossible, with the magnetic moment broken.
  **A missing K3 is not a result:** the Buneman wavelength is about 4 cells and the current filter damps it.
- **K4, anything else,** only with a positive signature and a control run that switches it off.

Measures on tracked orbits (the tracked set is fixed at the start; only ambient electrons are tracked):
magnetic moment, number of front crossings (with a hysteresis of 1 d_i), foot-point θ_Bn against time, and
the parallel and perpendicular work if N4 holds.

**Track stride.** If fewer than 300 tracked electrons fall in the residual of L, the stride is lowered for the
later runs (a prime that does not divide 768; track file at most 70 GB).

## 8. Criteria for "a hypothesis was obtained"

- **C1, the residual is real.** At least 100 electrons, at least 10 times the floor, at least 3σ (Poisson)
  above the upper bound of the reference; present in both hemispheres and in all four quadrants in azimuth;
  present in a second run with another seed (L-s2).
- **C2, the mechanism is shown on orbits:** tracked orbits or id histories show the signature (magnetic moment,
  encounters, change of the foot-point θ_Bn).
- **C3, it is not numerical and not the driver:** N1–N4 hold; the origin is ambient plasma; the electrons are in
  the valid domain; the late phase (front beyond 172 d_e) is excluded.
- **C4, it is testable:** a control parameter, a prediction that can fail, and the design of the run that tests it.
- **C5, it is new to the field:** stated against the direct predecessors on curved shocks, which have been read
  (Sect. 10). If a predecessor could not be obtained, the answer is at most "candidate".
- **C6, the research reviewer** does not return "needs rethinking" (review R2).

## 9. Runs

All on 8 GPUs, from a clean committed tree (`git status --porcelain` empty), with the GPUs free.

| Order | Name | Config | Purpose | Condition |
|---|---|---|---|---|
| 0 | check-L | `check-L.yaml` | all output paths, cost | done (Sect. 2) |
| 1 | L-null | `explore-L-null.yaml` | numerical floor N1–N3 | always |
| 2 | L | `explore-L.yaml` | main run | always |
| 3 | M | `explore-M.yaml` | second condition | after steps 1–4 of the analysis of L |
| 4 | M-null | (M without debris) | floor of M | only if M has a valid domain |
| 5 | H | `explore-H.yaml` | contrast: weak coupling | last; its electrons support no hypothesis |
| 6 | L-s2 | L with another seed, denser tracks | C1 and orbit statistics | only if L has a residual |

```bash
.venv/bin/python -m pypic.launch -n 8 examples/blast_wave/run.py \
  --config research/topics/spherical-shock-injection/configs/<config> \
  --out output/research/spherical-shock-injection/<name> --progress 148 \
  > output/research/spherical-shock-injection/logs/<name>.log 2>&1
```

- **Stop rules.** A run is stopped if a conservation contract fails (Gauss law > 10⁻¹⁰, div B > 10⁻¹²), if
  `fill` reaches 0.95, or if the dump `hot_e` is on course to exceed 60 GB (judged after its first four
  outputs). A stopped run is recorded as stopped, with the reason.
- **Records.** `runs/<date>-<commit>-<name>.json` and `runs/index.json`: commit, config, seed, GPUs, wall time,
  the place and size of the outputs, the commit of the analysis, the measures and the checks. Outputs stay
  out of git.
- **Order of the analysis** (from the review): 1 the floor; 2 the front and the valid domain; 3 the upstream
  spectrum; 4 the census against the baseline; 5 orbits of the residual.

## 10. Novelty

The direct predecessors on curved shocks are being read in parallel with the runs (Krauss-Varban & Burgess
1991; Krauss-Varban 1994; Vandas 1995, 2001; Savoini & Lembège 2001 and later; Anderson et al. 1979). Each gets
a note in `research/literature/`. A work that cannot be obtained is named as a limit.

## 11. Reviews

- **R1** (before the first run): this protocol, Sects. 3–8. Changes of numbers are committed before run 1.
- **R2** (after the analysis): results and the hypothesis, if any.

## 12. Checks of the analysis code

`analyze.py selftest` (synthetic data; all pass at the commit of this file):
the reflection map is its own inverse; the cold-electron energy; the reflected flux against the analytic
non-relativistic limit (0.01 %) and against a Monte Carlo sum at V_HT = 0.47, b = 3 (0.5 %); a fitted noisy
sphere returns V_HT = V/cos ϑ (2 %) and θ_Bn = ϑ (0.5°); back-tracing recovers synthetic reflections (time
within 0.1/ω_pe, θ_Bn within 0.1°) and returns nothing for an electron that never met the front; the
Maxwell–Jüttner tail against the loaded particles (0.03 %); the count model converges in its steps (10⁻⁵).
The front detection and the orbit measures get their synthetic checks before they are used on a run.
