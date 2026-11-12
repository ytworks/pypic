# Exploratory runs: protocol fixed before the runs

Topic `spherical-shock-injection`. Written 2026-11-13 (`ec6470e`), revised the same day after the review R1
and the 600-step cadence check, before any production run of L, M or H.
**Status: fixed. No production run exists at the commit of this revision. From here on a number changes only
with a dated notebook entry that names the new fact, and never after the run it would judge.**

This file answers the review of the exploratory conditions (issue #48, `reviews/2026-11-13-exploratory-conditions.md`)
and the review R1 of this protocol (`reviews/2026-11-13-protocol-r1.md`; the list of changes is in Sect. 11).
It replaces Sect. 6 (checks 2 and 6) and Sect. 8 of `exploratory-conditions.md`. The physical conditions of
L, M and H are unchanged. Code: `shockkin.py` (kinematics, count model), `analyze.py` (commands).
Units: c = 1, m_e = 1, lengths in d_e, times in 1/ω_pe.

## 1. The question of this exercise and its three possible answers

- **Question (from the owner, 2026-11-13).** Do runs L, M, H, analysed, give a research hypothesis about the
  electron acceleration mechanism that is novel and certain enough to build a study on?
- **A hypothesis is built only from a residual:** what is left of the energetic electrons after (i) the
  numerical floor (Sect. 5), (ii) the known classes S0, S1 and S2 (Sect. 3), (iii) the electrons that cannot be
  decided ("undetermined"), and (iv) the regions that hold no shock (Sect. 4) are taken away.
- **The answer is one of three** (criteria in Sect. 8):
  - **obtained:** all of C1–C6 hold;
  - **candidate:** a residual with a positive signature exists, and the missing evidence is named;
  - **not obtained from these runs:** with the conditions that a run would need.
- **Honest expectation.** One mass ratio (16), 2 particles per cell and V/c ≈ 0.2: the best outcome is a
  hypothesis with its test designed, not a result. Two limits are known before the runs, and each one alone
  holds the answer at "candidate" at most: two direct predecessors are unread (Sect. 10), and the baseline is
  not computed with guiding-centre orbits (S0′, Sect. 3.2).

## 2. Facts measured in the two check runs

### 2.1 The 16-step check of L (run 0)
`configs/check-L.yaml`, 8 GPUs, commit `3e97c6c`; outputs in `output/research/spherical-shock-injection/check-L`
(2.9 GB, not in git).

- All nine diagnostics wrote files. Gauss law 4 × 10⁻¹⁵, div B ≤ 1.3 × 10⁻¹⁶, fill 0.667.
- **Selection of the dump works:** with `gamma_min: 1.05` the smallest γ in the dump is 1.05, and the dump
  holds 17,782,723 ambient electrons at step 0 (a fraction 0.019628).
- **The loaded electrons are Maxwell–Jüttner** with T = v_te² = 0.01 m_e c², f ∝ exp(−(γ−1)/T): the measured
  tail fractions agree with it to 0.03 % at γ ≥ 1.05 and to 3 % at γ ≥ 1.15, while a Gaussian in u gives
  0.0166 at γ ≥ 1.05. The baseline uses the Maxwell–Jüttner form.
- **Ids.** Ambient electrons have ids below 2N (N = 768³ cells); `id // 2` is the linear index of the birth
  cell (checked: every position lies in its cell). Electrons injected at the open faces have ids ≥ 20N.
  Debris electrons are the species `de`; their weights are small and vary, so they are counted by weight.
- **Tracks:** the tracked set is `id % stride == 0`; momentum is half a step before position and fields.
- Its 0.90–1.06 s per step was the cost of writing every output every 4 steps, not the cost of a step (2.2).

### 2.2 The 600-step cadence check of L (review R1, change 5)
`configs/check-L-cadence.yaml`: L with tracks every step (stride 4733) and the production outputs, 8 GPUs,
commit `2d309c8`, 600 steps (t = 150) in 175 s. Outputs in `…/check-L-cadence` (4.9 GB); numbers in
`<run>/analysis/cadence.json` (`analyze.py cadence`).

- **Cost:** 290 ms per step (1.24 ns per particle and step). A run of 2960 steps takes about 15 minutes plus
  the cost of the denser tracks, which the null run of L measures. This corrects the "45–50 minutes" of the
  first version of this file.
- Total energy at t = 148: +1.7 × 10⁻² (the explicit scheme monitors it; no bound applies).
- **The dump `hot_e`** (γ ≥ 1.25, every 74 steps): 0.07, 0.95, 6.1, 18, 44, 73, 89, 104, 107 MB in its first
  nine outputs.
  - It is dominated by debris electrons: 3.3 × 10⁶ at t = 148 (3.4 × 10⁵ of them with γ ≥ 1.5), against
    3.5 × 10⁴ ambient electrons (329 with γ ≥ 1.5; at most 807, at t = 92).
  - These early ambient electrons are older than any shock (Sect. 4 asks for t ≥ 112). They belong to the
    start of the driver and are reported as such.
  - A whole run is expected at about 5 GB, far below the stop rule of Sect. 9.
- **Closure of the work** (184,349 tracked electrons stay in the box; the largest |B| on a track is 3.75 B₀).
  For each orbit with |Δγ| ≥ 0.1, the relative error of Σ q E·v Δt against Δγ:

  | tracks read every | orbits | median error | 90th percentile | orbits within 10 % |
  |---|---|---|---|---|
  | 1 step | 210 | 0.4 % | 0.7 % | 100 % |
  | 2 steps | 204 | 0.6 % | 3.3 % | 98 % |
  | 3 steps | 195 | 2.5 % | 10.5 % | 88 % |
  | 4 steps | 190 | 5.4 % | 21 % | 73 % |

  - Read every step, the recorded fields account for the energy change. Noise between the samples is not the
    limit; the sampling of the gyration is (88 samples per gyro-period in B₀ at every step, 23 in 3.75 B₀).
  - The orbits that met |B| ≥ 2 B₀ (192 of the 210) give the same figures within 3 points.
  - The work along B is about 1 % of the gain of these orbits: the early gain is across the field.
  - Read every 2 steps, the parallel work differs from its every-step value by 0.4 % of Δγ (median) and
    1.7 % (90th percentile); every 4 steps by 2 % and 8 %.
  - Only 10–14 orbits reach |Δγ| ≥ 0.2 in this early phase.
- **Decisions.**
  - `tracks_e` is recorded **every 2 steps**: the coarsest cadence at which at least 95 % of the orbits close
    within 10 %.
  - **Stride 947** in L and M (9.6 × 10⁵ electrons, 1481 samples, 55 GB on disk at the measured 38.6 bytes per
    sample and electron). R1 asked for 463; at every 2 steps that is 110 GB per run, and the disk has
    291 GB free. 947 is the fallback R1 names.
  - **Stride 4735 = 5 × 947** in the null run of L and in H (1.9 × 10⁵ electrons, 11 GB): every fifth tracked
    electron of L, with the same birth cells. R1 asked for one stride in L and its null. The null holds no
    reflection, and its uses (N3, N5) need thermal electrons only; H supports no hypothesis.

## 3. The known signals: reflection, leakage, trapping

### 3.1 One-step reflection (S0): the kinematics
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
- For a given V_HT the energy after a reflection does not depend on b and φ; they decide only how many
  electrons are reflected.
- **V_HT is not one number per reflection.** The foot point slows while the electron is in the layer: for a
  sphere d ln V_HT/dt = −V_n tan²θ_Bn / R. The measure is **ε = V_n tan²θ_Bn τ_dwell / R**. With
  V_n = 0.2175, R = 100 and the dwell of 3.2, ε is 0.20 at 50°, 0.43 at 60° and 1.1 at 70° in L. So in the
  band where S0 reaches γ ≥ 1.5 the energy is known only as a range. The mean gain on a curved front is below
  the plane value in adiabatic theory (`vandas1994acceleration`).

### 3.2 The per-electron test B1′ (the primary test)
For an energetic electron found ahead of the front:

1. take its guiding centre (r + u × B / (q B²)) and its parallel velocity;
2. follow it back along its field line as a free particle to the last time t_r it was at the foot point;
3. take the range of V_HT that the foot point of this field line had during [t_r − τ_dwell, t_r]: from
   V_HT(t_r) up to V_HT(t_r − τ_dwell), or up to 1 if the line touched the front later than t_r − τ_dwell;
4. undo the reflection (the map is its own inverse) for 16 values of V_HT in this range.

The electron **is explained by S0** if, for at least one of these values, its kinetic energy before the
reflection is at most E_th and it was outside the loss cone for b_max and φ = 0.

- **τ_dwell** is fixed here as 3 electron gyro-periods in B₀: 66 (L), 113 (M), 164 (H). A mirrored electron
  stays about 2.7 gyro-periods in the layer in `savoini2010origin`. If the tracked reflections of a run give a
  larger 90th percentile, that value replaces it; both results are reported.
- **The reflecting surface** is the mid-ramp of the field, not the density front: the outermost radius, within
  3 d_i of the density front, at which |B|/B₀ averaged over d_i reaches (1 + b)/2. It is fitted separately
  (`front_B`), over the bins with b ≥ 1.5. The distance between the two surfaces is reported.
- **E_th** is the energy above which the Maxwell–Jüttner distribution at the measured upstream temperature holds
  fewer than 10 macro-electrons in the whole box: 0.20 m_e c² at T = 0.010, 0.245 at T = 0.012.
- **b_max** is the largest measured compression over the detected front (bin averages). This makes the
  loss-cone test nearly empty, which errs on the side of "explained". B1′ is also reported with the local
  b(ϑ, t), and with the largest single-cell value as a test of a rippled front.
- **Two results are reported:** explained at the exit value V_HT(t_r) alone, and explained within the range.
  The residual is defined against the second.
- **The energies before the reflection are a criterion.** For the explained electrons of a band, their mean
  has to be within 30 % of the mean of the flux-weighted Maxwell–Jüttner distribution over the same selection,
  and the share above 10 T within a factor 2. An excess is reported as a pre-heated seed population, and the
  band's residual is not read until its source is known.
- **Undetermined** (neither explained nor residual):
  - the crossing lies outside the time range of the fit, or where `front_B` does not exist;
  - |v_∥ − V_HT(t_r)| < 0.1;
  - the flight time exceeds 0.1/ν of the electron's speed bin (Sect. 5, N3 and N3b).
- **Not done here: S0′,** guiding-centre orbits in the fitted, time-dependent front with the measured ramp. The
  range of V_HT is the minimum R1 accepts in its place. Consequence: a residual in a band with ε > 0.3 is at
  most "candidate".

### 3.3 Two further known classes
`savoini2010origin`, a full-PIC run of a curved front near the parameters of L, finds 70 % of the
backstreaming electrons mirrored, 20 % trapped in the parallel potential well of the overshoot (magnetic moment
not conserved), and 10 % leaked from downstream (entry near 88°, exit below about 60°). Both further classes
would fail B1′ and carry a positive signature, so both are named and tested before a residual is read.

- **S1, leakage from downstream.** An upstream electron that B1′ does not explain, and for which all hold:
  - it is inside the loss cone at its exit (local b);
  - it was in the downstream class in an earlier dump (id history; this is why `hot_e` starts at γ = 1.25 and
    is written every 74 steps);
  - its energy in the foot-point frame, γ′ − 1, is at most the thermal bound of the measured downstream T_e
    (the definition of E_th, with φ = 0).
- **S2, trapping in the overshoot.** An electron that stays in the front layer for at least 10 gyro-periods in
  B₀ with at least 2 bounces along the field (tracks), or is in the layer class in the corresponding number of
  successive dumps (id history). A dwell between 3 and 10 gyro-periods is "undetermined".
- **The residual** is what S0 (by B1′), S1, S2 and "undetermined" leave.

### 3.4 Counts (B2, a consistency check)
The number of reflected electrons per band of foot-point θ_Bn (10°) and per energy bin (γ 1.5–2, 2–3, 3–5, > 5)
that are still in the box at each dump, from the flux of a Maxwell–Jüttner upstream on the fitted front.

- Upper bound: b(ϑ, t) as measured, φ = 0, T_e as measured upstream.
- Lower bound: the same with the initial T_e and b one standard error lower; with a measured φ if its error
  is under 50 %, otherwise the lower bound carries the note "φ not measured".
- B2 is computed for the whole detected front and, separately, for the valid shock domain (Sect. 4). A moving
  compression reflects electrons whether or not it is a shock.
- If the measured count of S0-explained electrons falls outside [lower/2, 2 × upper] in the bands that hold at
  least 100 electrons, the baseline or the front measurement is wrong, and no residual is discussed until
  this is understood.

### 3.5 Prediction before the runs
A nominal sphere R = 21.12 + 0.2175 t (the shell radius, and 0.87 V_d from the example's measured front), a
Maxwell–Jüttner upstream at T = 0.01, 16 macro-electrons per d_e³, box half-width 192 d_e, φ = 0, one V_HT per
reflection (`analyze.py predict`, numbers in `data/predict-nominal.json`). The luminal angle is 77.4°.
The nominal prediction is the same for L, M and H: they share V_d, v_te, the shell and the box.

Macro-electrons with γ ≥ 1.5, both hemispheres, at the last dump (t = 740):

| b | reflected in the run | still in the box | of these γ 1.5–2 | γ 2–3 | γ ≥ 3 | bands that hold 99 % |
|---|---|---|---|---|---|---|
| 2 | 1.0 × 10³ | 5.3 × 10² | 5.3 × 10² | 0 | 0 | 50–70° |
| 3 | 6.7 × 10⁴ | 3.3 × 10⁴ | 3.2 × 10⁴ | 69 | 0 | 50–70° |
| 4 | 4.0 × 10⁵ | 1.9 × 10⁵ | 1.9 × 10⁵ | 1.8 × 10³ | 0 | 50–70° |

- **S0 gives γ ≥ 1.5 almost only for foot-point θ_Bn between 50° and 70°, and nothing at γ ≥ 3.** The cold
  limit alone would give 61°; the electrons below it are those with a large perpendicular momentum. With the
  range of V_HT of 3.1 the energies in this band spread, and the band widens to lower exit angles.
- **The count changes by a factor of about 400 between b = 2 and b = 4.** A count that exceeds S0 is therefore
  not evidence by itself; the per-electron test is the primary one.
- The real front slows and is not a sphere; the prediction is recomputed with the measured front. The nominal
  table stays as the record of what was expected.
- **The null run:** a Maxwell–Jüttner plasma at T = 0.01 holds 2 × 10⁻¹² macro-electrons with γ ≥ 1.5 in the
  whole box, and under one with γ ≥ 1.25. Every such electron in the null run is numerical or comes from the
  open faces.

### 3.6 Limits of the baseline
- The back-tracing assumes a uniform upstream field and free flight. It is valid only where the upstream is
  undisturbed (check 4 of the conditions record) and for flight times below 0.1/ν (Sect. 5).
- The front is treated as axisymmetric and smooth on the scale of the bins. Reflection by ripples smaller than
  a bin counts as "beyond S0" here; the single-cell pass of B1′ bounds this.
- The potential φ is hard to measure at 2 particles per cell. φ = 0 gives the largest reflected count, so the
  upper bound does not need it. A potential that changes in time (reformation, above about 72° in
  `savoini2010origin`) changes the energies; it is not modelled, and a residual above 70° names it as a
  possible cause.
- **Relaxation of the reflected beam.** The beam can drive Langmuir waves upstream, which spread its energies.
  Test: the spread of the energies of S0 electrons against their flight distance. A spread that grows with
  distance is this known effect, and widens E_th for the electrons concerned.

## 4. What counts as a shock: the front and its valid domain

From `slice_y` and `slice_x` (the two meridional planes through the centre; both contain B₀), every 74 steps.
Polar bins of 2 d_e in radius and 5° in polar angle about the centre. Each plane gives four quadrants; they
are compared before they are folded. d_i = 4 d_e in all runs. Ω_ci⁻¹ = 56 (L), 96 (M), 139 (H).

- **Front position R(ϑ, t).** The outermost radius at which the density of the ambient ions `n:i`, averaged over
  d_i in radius, is at least 1.5 n₀ in two neighbouring bins; the crossing is interpolated. This also works
  where B is not compressed (the quasi-parallel side). At the equator it is compared with the example's
  definition (B_z ≥ 1.25 B₀).
  - It is repeated with 1.3 n₀ and 1.7 n₀. A band is used only if R moves by less than 1 d_i and V_HT by less
    than 10 % between them.
  - The two planes have to agree within 1 d_i in a band, or the band is not used.
- **Front model.** A least-squares fit of the positions with {1, cos 2ϑ, cos 4ϑ} × {1, τ, τ², τ³}; V_HT, θ_Bn
  and the normal speed V_n follow analytically. The same for `front_B` (Sect. 3.2).
  - The rms residual has to be below d_i.
  - A bound on R does not bound R_ϑ. A band is used for B1′ only if V_HT from the fit and from differences of
    neighbouring bins agree within 15 %, and the four quadrants agree within 15 %.
- **Compression b(ϑ, t).** The largest bin average of |B|/B₀ between 3 d_i behind and 1 d_i ahead of the front.
- **Debris edge.** The outermost radius at which `n:di` or `n:de` is at least 0.1 n₀; reported also for 0.03
  and 0.3 n₀.
- **Valid shock domain.** A bin (ϑ, t) belongs to it if all of these hold:
  - (a) a front position exists;
  - (b) the debris edge lies behind the position of B_max, and the front leads it by at least the foot length
    0.7 V_n/Ω_ci with the measured V_n (at V_n = 0.2: 2 d_i in L, 3.4 d_i in M, 4.9 d_i in H). The coefficient
    is the reviewer's, quoted from memory (0.68 for specular reflection at a perpendicular shock); it was not
    checked against a paper;
  - (c) the normal speed from five neighbouring outputs is at least 1.5 times the fast speed v_f for the local
    θ_Bn, and V_n − 2σ_V > 1.1 v_f with σ_V the error of that slope. v_f uses the measured upstream
    temperatures (with the initial ones: L 0.071–0.085, M 0.046–0.062, H 0.046–0.054);
  - (d) (a)–(c) hold over 2 Ω_ci⁻¹: 6 (L), 11 (M), 15 (H) outputs in a row;
  - (e) the shock is at least 2 Ω_ci⁻¹ old (t ≥ 112, 192, 278), and R ≤ 192 − 8 − V_n/Ω_ci (at V_n = 0.2:
    173, 165, 156 d_e).
- **Consistency check, reported per band:** the normal velocity of the ambient ions behind the front (from
  `ekin:i` and the density) is within 30 % of V_n (1 − n₀/n).
- **Outside the valid domain** the electron data are reported as "fast compression, no shock". They are not
  read as shock acceleration, in any run.
- **Luminal angle.** From the measured front: the foot points with V_HT ≥ 1. The band is marked in every
  table and figure.
- **θ_Bn of an electron** is the θ_Bn at the foot point of its field line: at the time of its reflection for
  B1′, and at the time of the dump for the census.

## 5. Numerical floor

`configs/explore-L-null.yaml`: the same plasma, field, grid, seed, length and outputs as L, without the debris.

- **N1, counts.** Electrons with γ ≥ 1.5 per dump, by origin (ambient, injected at a face) and by position.
  The expectation is zero. A count in a run is a **signal** only if it is at least 100 and at least 10 times
  max(count of the null run in the same class, 3).
- **N2, heating.** T_e(x, z, t) from `spec_meridional` and T_e(t) from `spec_up` (fit to the bulk). The old
  bound of check 6 (1.25) is dropped: the baseline is computed with the initial and with the measured
  temperature, and a band in which the two differ by more than the residual is "explained by heating".
- **N3, numerical scattering of thermal electrons.** The pitch-angle scattering rate from the tracks of the
  null run, ν = ⟨(Δμ)²/(1 − μ²)⟩/Δt at lags of 4, 16 and 64 time units, per speed bin (0–0.1, 0.1–0.2,
  0.2–0.3, > 0.3). Uses:
  - the reflected count B2: if ν τ > 0.1 for thermal electrons, with τ the crossing time of the layer (its
    measured width over V_n), the loss cone is being filled numerically, and B2 carries that note;
  - a mechanism that needs scattering (stochastic shock drift acceleration) is not assigned unless ν τ < 0.1.
- **N3b, scattering of energetic electrons, in situ.** The null run holds no energetic electron. So: every
  electron of `hot_e` in L that is in the upstream class in two successive dumps gives Δμ, Δγ and Δρ_P over
  18.5 time units. From these, the pitch-angle scattering rate, the energy diffusion and the cross-field
  step, per speed bin. This rate limits the flight time of the back-tracing (Sect. 3.2). It holds physical and
  numerical scattering together; the pair of N6 separates them.
- **N3c, fluctuation level.** The power of the field fluctuations upstream of the front in L against the null
  run at the same place and time (from the slices). A ratio above 2 means that the floor of the null run does
  not stand for the upstream of L; N3 and N5 then carry that note.
- **N4, closure of the work, per orbit.** On the tracks, Σ q E·v Δt against Δγ (the momentum is recorded half a
  step before the fields, so the velocity is centred first). The split into parallel and perpendicular work is
  used for an orbit only if its own error is below 10 %. The check of 2.2 found 98 % of the orbits within it
  at the production cadence. An orbit that fails is classified by the magnetic moment and the number of
  encounters only; the share of failing orbits is reported per class.
- **N5, numerical transport across the field.** In axisymmetric fields the canonical angular momentum of an
  electron is conserved. In the uniform upstream field this is the **canonical radius**
  ρ_P = (ρ² + 2ρ u_φ/(q B₀))^(1/2), with ρ_P² = ρ_gc² − ρ_L². For a thermal electron ρ_P differs from the
  radius of its birth cell by at most about a thermal gyroradius (0.35, 0.6, 0.87 d_e in L, M, H) and a cell.
  - N5 is the distribution of Δρ = ρ_P − ρ_birth against time for the tracked electrons of the null run.
  - Limit: these electrons are thermal and never crossed a layer. N3b gives the step of energetic electrons
    upstream; the transport inside the layer has no null, only the pair of N6.
- **N6, particle number.** A pair of smaller runs of L (384³ cells, `rm: 8.0`, the same cell size) at 2 and at
  8 particles per cell. It is made only if L shows a residual or a K1 signal; its configs are written and
  plan-checked before it starts. The residual fraction (residual over S0-explained electrons, same bands)
  has to agree between the two within 30 % or within 2σ (Poisson), whichever is larger.
- A null run of M is made only if M has a valid domain.

## 6. The census of the energetic electrons

One energy threshold for all runs: **γ ≥ 1.5** (50 kT_e; the ram energy of an ion is 0.5 m_e c² at V_d).
The dump `hot_e` holds every electron with γ ≥ 1.25, every 74 steps; those between 1.25 and 1.5 serve the
histories, N3b and the statistics of S0 at lower angles. Every electron of the dump gets:

- **origin**, from the id: debris (species `de`, counted by weight); injected at a face (id ≥ 20N); ambient
  born inside the shell (birth radius ≤ 27.5 d_e); ambient swept (the front has passed its birth cell);
  ambient not yet reached;
- **position class**, from its guiding centre and the fitted front: upstream (more than 1 d_i ahead of the
  front, or on a field line that the front has not reached), front layer (within 1 d_i ahead and 3 d_i behind),
  downstream; and whether the foot point is in the valid domain;
- **foot-point θ_Bn** (10° bands) and the luminal flag;
- **canonical radius** ρ_P and Δρ = ρ_P − ρ_birth (ambient electrons);
- for upstream electrons: the result of B1′ (explained at the exit value, explained within the range, S1,
  undetermined, residual);
- **history across dumps**, by matching ids: the dump in which it first appears, and its class there.

Tables: counts per dump × origin × position class × θ_Bn band × energy bin, with the null-run count beside
each, and split by hemisphere and quadrant. A band with fewer than 100 electrons is "not measurable" and is
not used for a hypothesis.

Downstream, the reference is the thermal tail: the count above γ = 1.5 expected from a Maxwell–Jüttner fit to
the bulk of the local spectrum (`spec_meridional`), with the fit temperature varied by ±10 %. Downstream
electrons enter a hypothesis only with tracked orbits or id histories that show how they gained their energy.

## 7. Candidate signatures, looked for in the residual only

- **S0, S1, S2:** known (Sect. 3).
- **K1, transport across the field.** The foot-point θ_Bn of a field line changes without any transport
  (sin ϑ = ρ/R; published, Sect. 10), so K1 is not defined by it. A K1 electron is an ambient upstream electron
  with γ ≥ 1.5 whose |Δρ| exceeds both 2 d_e and the 99.9th percentile of N5 at the same time.
  - **Signal:** at least 100 such electrons in the valid domain, at least 10 times the count expected from the
    null percentile, and at least 1 % of the upstream energetic electrons of the same bands.
  - The sign is reported: negative is towards the axis.
  - Δρ needs no tracks and no fitted front. It is exact where the fields are axisymmetric, so a K1 signal
    means fields that are not: waves and ripples of the front (physical), or noise (numerical). N3b and the
    pair N6 decide between them; without N6 a K1 signal is "candidate".
- **K2, several encounters.** Two or more crossings of the front, and an energy above one reflection. On a
  convex front with straight upstream field lines a freely flying electron cannot return, so K2 needs
  scattering upstream. It is read only with N3b: scattering at the rate of the null floor is numerical.
- **K3, gain in the superluminal band,** where reflection is impossible, with the magnetic moment broken.
  **A missing K3 is not a result:** the Buneman wavelength is about 4 cells and the current filter damps it.
- **K4, anything else,** only with a positive signature and a control run that switches it off.

Measures on tracked orbits (the tracked set is fixed at the start; only ambient electrons are tracked):
magnetic moment, number of front crossings (with a hysteresis of 1 d_i), dwell in the layer, foot-point θ_Bn
and ρ_P against time, and the parallel and perpendicular work where N4 holds for the orbit.

**Tracked sets** are fixed by 2.2: stride 947 every 2 steps in L, M and L-s2; 4735 in the null runs and H. With
947, the nominal prediction of 3.5 puts about 1, 70 and 420 S0 reflections on tracks for b = 2, 3, 4. If L has
fewer than 30 tracked reflections, C2 rests on id histories, and the false-negative rate of C1 on its bound.

## 8. Criteria for "a hypothesis was obtained"

- **C1, the residual is real.** All of:
  - at least 100 electrons, and at least 10 times the floor (N1);
  - upstream: more than 3 × f_FN × (number of S0-explained electrons of the same bands). f_FN is the share of
    tracked single reflections with the magnetic moment conserved within 20 % that B1′ fails to explain. With
    fewer than 100 such tracks, f_FN is the 95 % upper bound (3/n for no failure);
  - downstream: at least 3 times the thermal tail at the upper fit temperature;
  - present in both hemispheres and in all four quadrants in azimuth;
  - present in a second run with another seed (L-s2). The definition of the residual (bands, energy bin,
    class, expected count) is frozen in a commit before L-s2 starts.
- **C2, the mechanism is shown on orbits:** tracked orbits or id histories show the signature (magnetic moment,
  encounters, dwell, Δρ).
- **C3, it is not numerical and not the driver.** All of:
  - N1–N5 hold, and the pair N6 agrees;
  - the origin is ambient plasma loaded at the start; electrons injected at a face (id ≥ 20N) do not count;
  - the electrons are in the valid domain, and first appear in `hot_e` at a time and place inside it;
  - the field line does not thread the debris: ρ_P > ρ_debris,max(t), the largest distance from the axis at
    which `n:de` or `n:di` reaches 0.03 n₀;
  - the late phase (Sect. 4, (e)) is excluded.
- **C4, it is testable:** a control parameter, a prediction that can fail, and the design of the run that tests
  it. The hypothesis is stated in dimensionless numbers (ε or R/L, V_HT/v_te, V_sh/c, m_i/m_e), with their
  values for a remnant. At V_sh/c ≈ 0.01 one reflection in the 50–70° band gives no such energies; a
  hypothesis has to say what survives that.
- **C5, it is new to the field:** stated against the direct predecessors on curved shocks, which have been read
  (Sect. 10). If a predecessor could not be obtained, the answer is at most "candidate".
- **C6, the research reviewer** does not return "needs rethinking" (review R2).

## 9. Runs

All on 8 GPUs, from a clean committed tree (`git status --porcelain` empty), with the GPUs free.

| Order | Name | Config | Purpose | Condition |
|---|---|---|---|---|
| 0 | check-L | `check-L.yaml` | all output paths | done (2.1) |
| 0b | check-L-cadence | `check-L-cadence.yaml` | track cadence, cost, size of `hot_e` | done (2.2) |
| 1 | L-null | `explore-L-null.yaml` | numerical floor N1–N3, N5; cost | always |
| 2 | L | `explore-L.yaml` | main run | always |
| 3 | M | `explore-M.yaml` | second condition | after steps 1–4 of the analysis of L |
| 4 | M-null | (M without debris) | floor of M | only if M has a valid domain |
| 5 | H | `explore-H.yaml` | contrast: weak coupling | last; its electrons support no hypothesis |
| 6 | L-s2 | L with another seed | C1 and orbit statistics | only if L has a residual or a K1 signal |
| 7 | L-half-p2, L-half-p8 | 384³, `rm: 8.0`, 2 and 8 per cell | N6 | as run 6 |

```bash
.venv/bin/python -m pypic.launch -n 8 examples/blast_wave/run.py \
  --config research/topics/spherical-shock-injection/configs/<config> \
  --out output/research/spherical-shock-injection/<name> --progress 148 \
  > output/research/spherical-shock-injection/logs/<name>.log 2>&1
```

- **Cost.** About 15 minutes per run at 290 ms per step, plus the tracks; measured on L-null and on L.
- **Disk.** 291 GB free. Planned: 66 GB for L, M and L-s2 each (tracks 55, `hot_e` 5, the rest 6); 19 GB for a
  null run and for H; 8 GB for the pair. The whole set is about 270 GB, so L-s2 cannot have denser tracks.
- **Stop rules.** A run is stopped, or not started, if:
  - a conservation contract fails (Gauss law > 10⁻¹⁰, div B > 10⁻¹²), or `fill` reaches 0.95;
  - the free space is below 1.5 times the planned size of the run (checked before the start);
  - `hot_e` is on course to exceed 60 GB. Judged at t = 203 (12 outputs): the ambient part scaled with R³
    (R = 21.12 + 0.2175 t), the debris part kept at its last value.

  A stopped run is recorded as stopped, with the reason.
- **Records.** `runs/<date>-<commit>-<name>.json` and `runs/index.json`: commit, config, seed, GPUs, wall time,
  the place and size of the outputs, the commit of the analysis, the measures and the checks. Outputs stay
  out of git.
- **Order of the analysis** (from the review of the conditions): 1 the floor; 2 the front and the valid domain;
  3 the upstream spectrum; 4 the census against the known classes; 5 orbits of the residual.

## 10. Novelty

The direct predecessors on curved shocks were read on 2026-11-13, before the first run; `literature.md`,
"Curved fronts", places them. What it means for the signatures of Sect. 7:

- **S0** and the change of the foot-point θ_Bn in time are published. No novelty is claimed.
- **An energy that differs from the plane-shock value on a curved front** is published in adiabatic theory
  (`vandas1994acceleration`). An S0 energy off the plane value is therefore not a finding.
- **Leakage and trapping** (S1, S2) are published for a curved front in 2-D full PIC (`savoini2010origin`),
  with electrons that enter near θ_Bn = 88° and leave near 48°. Seeing them here is a reproduction.
- **K1 as transport across the field (Δρ)** was not found in the works that could be reached. It is the
  narrow place where a new statement is possible.
- **K2 is published where a bent field line or a non-convex front closes the path** (`kong2016electron`,
  `liu2026efficient`, `pulupa2008structure`). On a straight field and a smooth convex front it was not found.
- **Limit.** `kraussvarban1991electron` and `giacalone2017acceleration` could not be obtained in full; eleven
  further works are known by abstract only. Until the two are read, criterion C5 cannot be met, and the verdict
  is "candidate" at most.

## 11. Reviews

- **R1** (2026-11-13, before the first run; `reviews/2026-11-13-protocol-r1.md`). Its verdict: with the first
  version "obtained" should not be reachable; the realistic best outcome is "candidate". Its eight required
  changes, and where they are:
  1. track stride: 947 instead of the 463 asked for, and 4735 in the null run (reasons in 2.2);
  2. `hot_e` from γ = 1.25, every 74 steps (configs, `2d309c8`); stop rule with R³ (Sect. 9);
  3. the second meridional plane `slice_x` (configs; Sect. 4);
  4. `ekin:i` in the slices (configs; Sect. 4);
  5. the cadence check (2.2);
  6. the range of V_HT, the mid-ramp surface, S1, S2, "undetermined" (Sect. 3). S0′ is not done;
  7. the numbers of the valid domain (Sect. 4); N3 per speed bin, N3b, N3c, N5 (Sect. 5);
  8. K1 by Δρ, C1 with the false-negative rate, C3 with the connection to the debris, the first appearance
     and the pair N6, C4 in dimensionless form, the frozen definition before L-s2 (Sects. 7 and 8).
- **R2** (after the analysis): results and the hypothesis, if any.

## 12. Checks of the analysis code

`analyze.py selftest` (synthetic data; 18 of 18 pass at the commit of this revision):
- the reflection map is its own inverse; the energy after a reflection; the cold-electron energy;
- the reflected flux against the analytic non-relativistic limit (0.01 %) and against a Monte Carlo sum at
  V_HT = 0.47, b = 3 (0.5 %); the count model converges in its steps (10⁻⁵);
- a fitted noisy sphere returns V_HT = V/cos ϑ (2 %) and θ_Bn = ϑ (0.5°); back-tracing recovers synthetic
  reflections (time within 0.1/ω_pe, θ_Bn within 0.1°) and returns nothing for an electron that never met the
  front; the Maxwell–Jüttner tail against the loaded particles (0.03 %);
- the front in synthetic slices: position (bias +0.7 d_e from the mean over d_i), speed and normal (2.6 %,
  2.6°), valid domain, compression and debris edge; no detection in noise;
- orbits in given fields: closure of the work (10⁻⁴), conservation of the magnetic moment, zero scattering
  rate in a uniform field;
- the temperature fit and the thermal bound; the S0 test explains every once-reflected thermal electron and no
  twice-reflected one.

The changes of this revision (the surface `front_B`, the thresholds 1.3 and 1.7, the second plane, the range
of V_HT, S1, S2, Δρ, N3b) get their synthetic checks before they are used on a run. None of them is needed to
start a run: they read the outputs that the configs already write.
