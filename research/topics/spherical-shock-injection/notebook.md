# Notebook — spherical-shock-injection

Dated lab notes, newest first. Write down what was tried, what was seen, and what was decided, failed attempts included.
Never rewrite a past entry. To correct one, add a new entry that refers to it.

## 2026-11-13 — Run L: no signal at γ ≥ 1.5; the returns to the upstream are not adiabatic mirror reflections
- **Run L.** `configs/explore-L.yaml`, 8 GPUs, commit `91f79ff`, seed 12345, 2960 steps in 1298 s (438 ms per
  step). Output `output/research/spherical-shock-injection/L`, 62 GB, not in git. Gauss law ≤ 5.2 × 10⁻¹⁵,
  div B ≤ 1.9 × 10⁻¹⁴, fill ≤ 0.674. Energy budget +10.4 % at t = 740 (monitored only; the null run has +12.7 %,
  so the growth is the numerical heating of the electrons). Record: `runs/2026-11-13-91f79ff1d-L.json`.
- **Front (`analyze.py front`).** 443 of 738 (time, θ) bins show a front; 241 are valid, 166 usable for B1′.
  - No front within 10° of the pole, and no valid front below θ = 30°: the fronts there do not pass criteria
    b–d (ahead of the debris, speed, duration). The valid domain is θ = 30–90°, t = 185–666.
  - In the valid domain V_n = 0.139–0.157, fast Mach 1.75–2.06, compression b = 1.6–2.1 per bin (largest cell
    3.7). The shock is weak: the piston moves at 3.5 v_A, the front at 2.1 v_A.
  - Foot points are superluminal (V_HT ≥ c) for θ_Bn ≥ 77–83°. The front leads the debris by 15–35 d_e; the ion
    flow agrees with V_n (1 − n₀/n) within 30 % in 240 of 241 valid bins; fit rms 0.63 d_e.
- **Census (`census`, protocol Sect. 6).** Ambient electrons with γ ≥ 1.5 per dump: 486 at t = 74, 330 at
  t = 148, 40 at t = 222, at most 4 later. The early ones are in the layer or behind it and precede every valid
  front (the start of the driver, as the 16-step check saw). Upstream of the front: 3 at t = 222 (2 on the
  detected front, both undetermined), none at any other dump. The null run has none.
  - **Result at the registered threshold: no band reaches 100 electrons, so there is no signal and no residual
    at γ ≥ 1.5 in L.** K1: one electron (t = 222); K2 and K3: none.
  - `hot_e` is dominated by debris electrons (up to 600 ambient macro-electron weights upstream at t = 296).
    They are the driver's, separated by species, and support nothing here.
  - Deviation from Sect. 6, conservative: an electron within d_i of the front is in the layer class even when
    its field line has not reached the front (the protocol reads "upstream" for it).
  - N3b: ν = 9.4 × 10⁻⁴ for γ 1.25–1.5 (2373 pairs of dumps). Only 5 pairs above γ = 1.5, so the rate of the
    null run (1.8 × 10⁻³) is used there.
- **Tracks (`tracks`, protocol Sect. 7).** 956 674 tracked electrons: 675 871 never left the upstream, 147 081
  transmitted, 96 447 in the layer at the end, 6416 back upstream after an encounter, 391 returned from
  downstream, 10 142 start behind the front, 20 326 outside the range of the fit.
  - γ_max ≥ 1.5 on 11 orbits, all of them starting behind the front. Returns with γ_out ≥ 1.25: 7. This is
    below the 30 the protocol asks for, so C2 would rest on id histories and the false-negative rate of C1 on
    its bound (upper limit 1.0): L cannot test C1 or C2.
  - N4 passes: 1289 of the 1341 orbits with |Δγ| ≥ 0.1 close the work within 10 % (49 of the 52 failures
    start behind the front). Parallel share of the work: median 0.19.
  - The returns at all energies: 6295 with one encounter, 3224 of them reached the density front. μ is kept
    within 20 % by 381 of the 3224 (12 %); μ_out/μ_in 0.25 / 1.40 / 9.0 (10 / 50 / 90 %). Gain: median 0.039,
    99 % 0.164.
  - Dwell in the layer: median 250, 90 % 418, against the 66 the protocol assumed. Caveat: the layer of this
    measure is 4 d_i = 16 d_e thick, so part of the dwell is its definition (width / V_n ≈ 107) and it is
    not by itself a sign of trapping.
- **N3 with the measured dwell.** ν τ = 0.062 at the protocol's 66 (γ 1.25–1.5), but 0.23–0.39 at the measured
  250–418. By the rule of N3 (ν τ < 0.1) the loss of μ cannot be assigned to a physical scattering mechanism
  in L at 2 particles per cell.
- **N3c (`n3c`).** Power of the fluctuations ahead of the front, L over its null: at most 1.67 in the valid
  domain (no note), up to 3.02 on the detected front outside it (the note holds there).
  - Three measures, since the protocol fixes no scale: below 1 d_i; about the mean along the front; about the
    initial field. A fourth, "below 4 d_i", was written first and dropped before the result was used: a box
    mean 4 d_i wide reaches the ramp from the 1–4 d_i shell and measured the ramp, not fluctuations (ratio in
    the tens). Reported here because it was dropped after its numbers were seen.
- **Potential (`phi`).** φ = −Σ (E·B)/B_pol ds along the in-plane field lines from 3 d_i ahead to 6 d_i behind
  the density front, by foot angle (5°) and quadrant; the null run gives the noise.
  - Chosen after seeing the data: the mean over 5 outputs, and the mean over all outputs of a bin where the
    windowed value is not measured (error above 50 %). Neither changes a registered criterion.
  - Mean over the outputs: 0.015–0.048 for θ = 10–70° (1.5–4.7 T_e; error below 50 % in every bin), largest at
    15–25°. At θ > 70° not measured (0 ± 0.005: the field lines graze the front). Windowed: measured in 157
    of 443 bins, 63 of the 241 valid ones.
- **Baseline B2 against the census (`baseline`)**, lower bound with the measured φ, upper bound with φ = 0.
  - γ ≥ 1.5: no cell with 100 measured electrons in either domain, so the rule of Sect. 3.4 (written for such
    cells) does not apply. But the model predicts a beam that is not there: on the detected front 808–2820
    electrons at t = 185, 1343–3890 at t = 296, 740–1529 at t = 407 (lower–upper bound), against 2, 0 and 0
    measured upstream.
  - γ 1.25–1.5 — **an extra bin, not registered, exploratory only.** Cells with ≥ 100 electrons: detected front
    86 (36 inside [lower/2, 2 × upper], 27 below, 23 above); valid domain 35 (11 inside, 19 below, 5 above).
    At t = 351.5 by band (model lower–upper | measured): 10–20° 0.2–3 | 201; 20–30° 19–138 | 525; 30–40°
    239–2391 | 557; 40–50° 5790–20 998 | 352; 50–60° 20 397–86 337 | 102; 60–70° 8651–40 240 | 2. The
    measured electrons sit at low θ_Bn, where a mirror gives almost nothing, and are missing by a factor
    200 or more at 50–70°, where it should give the most. By the rule of Sect. 3.4 **no residual is
    discussed**: the baseline does not describe the run.
  - b enters with its scatter over the cells of a bin as its standard error.
- **Share of the encounters that end in a return, all energies** (new; model against tracks, by foot θ_Bn at
  entry, detected front):

  | θ_Bn | model, φ = 0 | model, measured φ | tracks (± binomial) | tracks, only returns from behind the front |
  |---|---|---|---|---|
  | 10–20° | 0.131 | 0.0007 | 0.089 ± 0.012 | 0.057 |
  | 20–30° | 0.132 | 0.0027 | 0.089 ± 0.004 | 0.050 |
  | 30–40° | 0.115 | 0.0037 | 0.084 | 0.044 |
  | 40–50° | 0.093 | 0.0052 | 0.073 | 0.036 |
  | 50–60° | 0.053 | 0.0053 | 0.052 | 0.025 |
  | 60–70° | 0.013 | 0.0020 | 0.033 | 0.016 |
  | 70–80° | 0.0002 | 0.0001 | 0.021 | 0.013 |
  | 80–90° | 0 | 0 | 0.010 | 0.006 |

  - With the measured potential the adiabatic mirror returns 0.07–0.5 % of the encounters; the tracks return
    1–9 %, and still 0.6–6 % when only orbits with a sample behind the density front count. The valid domain
    gives the same picture (0.8–9.3 % against 0.01–0.6 %). The φ = 0 model is within a factor 1.5 at 10–60°,
    but φ = 0 contradicts the measurement.
  - **Reading: most returns in L are not adiabatic mirror reflections in a static potential**; μ says the
    same (12 % keep it), and so does the energy: at 50–60° the φ = 0 model has the right share (0.053 against
    0.052) but its beam at γ 1.25–1.5 is absent (above). The returning electrons do not carry 2 V_HT.
  - What this does not show: a physical mechanism. Open explanations, none excluded by L: numerical scattering
    (ν τ = 0.23–0.39, above); time dependence of the front and of φ (the model is static per foot); returns
    that are kinematic only — half of the "back upstream" orbits never reach δ ≤ 0, and a decelerating,
    expanding front lowers V_HT and θ_Bn at a field line's foot, so an overtaken electron can outrun the front
    again; the model counts only feet slower than light.
  - Two more, both known physics and both to be tested on M before anything is claimed. (1) Leakage (class S1):
    a field line first touches the sphere at θ_Bn = 90°, and its foot then moves along it at V_HT = V_n / cos θ,
    which falls towards V_n = 0.15. The heated electrons behind the front (T ≈ 0.015, thermal speed 0.12)
    can then overtake the foot. This fits the excess at low θ_Bn and the small gain. (2) Marginal
    adiabaticity by design: an upstream electron crosses a ramp one d_i wide in d_i / V_n = 27, which is 1.2
    gyro-periods (Ω_ce d_i / V_n = (m_i/m_e) v_A / V_n = 7.6). The ramp width was not measured; this is an
    estimate.
- **State of the question after L.** At the registered threshold L holds no energetic upstream electron of
  the ambient plasma, in agreement with a fast Mach number of 1.8 and b ≈ 2. The one robust finding is negative
  for the baseline: S0 with the measured potential does not describe the returning electrons at any energy.
  Whether that is physics is undecided at 2 particles per cell.
- Next: M (fast Mach 4.05 by design) from commit `ac33a41`, then the same analyses.

## 2026-11-13 — Null run of L: the numerical floor; front measurement revised to protocol Sect. 4
- **Run L-null.** `configs/explore-L-null.yaml` (no debris), 8 GPUs, commit `0114ef9`, seed 12345, 2960 steps in
  778 s (262 ms per step, 1.16 ns per particle and step). Output
  `output/research/spherical-shock-injection/L-null`, 15 GB, not in git (tracks 11 GB, reduced 4.3 GB).
  Gauss law ≤ 2.2 × 10⁻¹⁵ and div B ≤ 1.6 × 10⁻¹⁵ at the 10 steps that record them; fill ≤ 0.667.
- **Energy.** The budget grows linearly to +12.7 % at t = 740: numerical heating at 2 particles per cell, all
  of it in the electrons (electron kinetic energy +10.1 %, ion +0.2 %).
- **N1 (counts).** No electron with γ ≥ 1.5 in any of the 41 dumps, whatever its origin; largest γ 1.374; at
  most 15 per dump above γ = 1.25. With a floor of 0 the count criterion of a signal in L is the fixed one:
  100 or more.
- **N2 (temperature).** T_e in the corner column goes 0.01000 → 0.01095 at t = 740, nearly linearly (0.01004,
  0.01029, 0.01054, 0.01075 at t = 148, 296, 444, 592). The meridional map is uniform at t = 740: median
  0.0110–0.0111 in every ring, 5–95 % range 0.0104–0.0117. Above E = 0.1 the corner spectrum holds 1.37 times
  the Maxwell–Jüttner count at the fitted T (214 against 156 in weight); almost nothing above 0.2, none above
  0.5. T_i stays 0.0099. The sound speed from the measured temperatures goes 0.0455 → 0.0466.
- **N3 (scattering).** 191 335 tracked electrons, 167 026 of them in the box to the end. Pitch-angle rate ν
  at lags 4 / 16 / 64, by speed: v < 0.1: 2.6e-2 / 1.6e-2 / 8.4e-3; 0.1–0.2 (the bulk, 158 731 orbits):
  1.2e-2 / 7.8e-3 / 5.0e-3; 0.2–0.3: 5.0e-3 / 3.5e-3 / 2.6e-3; v > 0.3 (484 orbits): 2.5e-3 / 1.8e-3 / 1.3e-3.
  ν falls with the lag, so the motion is not a pure diffusion in pitch angle. |γ(end) − γ(0)|: median 0.010,
  99 % 0.058, largest 0.143; 63 orbits change by 0.1 or more.
  - Consequence, before L is seen: for the bulk electrons ν τ exceeds 0.1 once a layer takes 20 or more to
    cross (ν τ = 0.16–0.32 at τ = 20), and the dwell of 66 gives 0.33–0.55. For v > 0.3 it is 0.09–0.17 at 66.
    Unless L's measured layer is thinner than that, B2 will carry the note "loss cone filled numerically" and
    no mechanism that needs scattering can be assigned from these runs.
- **N5 (birth radius and ρ_P).** The birth cell decoded from the id lies within 0.25 d_e of the position at
  t = 0 for every tracked electron: the decoding is right. |ρ_P − ρ_birth| (median / 99 % / 99.9 % / max):
  0.26 / 0.99 / 1.26 / 1.92 d_e at t = 0 (the gyroradius), 0.76 / 2.90 / 3.71 / 5.29 at t = 740. The drift of
  ρ_P itself grows as √t: 99.9 % within 1.07 (t = 74), 2.20 (t = 296), 3.51 d_e (t = 740). A Δρ of one
  d_i = 4 d_e in L is therefore above the 99.9 % level of the null at every time.
- **Front measurement, revised to the protocol's Sect. 4** (`analyze.py front`, `shockkin.py`). Run on L-null
  as its null test: 0 of 738 (time, θ) bins show a front in both planes; nothing valid, nothing usable.
  23 of 23 self-tests pass; the new ones are a clean spherical front (flags all true, speed error against its
  own σ_V 1.3), and three fronts that must fail (planes 15 % apart, a ramp 20 d_e wide, debris 5 d_e behind).
- **Choices made while writing the code, where the protocol leaves the reading open.** None changes a number
  of the protocol.
  1. Ion-flow check (Sect. 4): from `flux:i` / `n:i`, not from `ekin:i`. The energy density cannot separate
     the bulk flow from the thermal motion. The criterion (30 %) is unchanged.
  2. Agreement of V_HT: bins in which every compared value is superluminal (V_HT ≥ 1) count as agreeing.
     A ratio of two infinite speeds carries no information.
  3. "The four quadrants agree within 15 %": the spread among the four single-quadrant fits and the folded
     fit, per plane.
  4. R is the mean of the two planes (`slice_y`, `slice_x`) where both show a front, and missing otherwise.
  5. Outputs in a row for criterion d: ceil(2/Ω_ci/Δt − 0.1) = 6 / 11 / 15 for L / M / H, as the protocol says.
  6. The mid-ramp B surface is searched only up to 3 d_i ahead of the density front, and is missing when it
     lies more than 3 d_i behind it.
  7. The comparison of V_HT between thresholds uses the fits, not the differences.
  8. E_th and v_f use the measured upstream temperature at each output (N2), since T_e rises by 9.5 %.
- **Units confirmed on the null** (needed by the checks above): the ion spectrum of `spec_up` is in m_e c²
  (T_i = 0.0099); `n:i` is in n₀ (noise 0.36 per cell at 2 per cell); `flux:i` is n v in c (noise 0.0096).
- Next: run L from the commit that holds this entry.

## 2026-11-13 — Review R1 of the protocol; cadence check; protocol revised before the first run
- **Review R1** (`reviews/2026-11-13-protocol-r1.md`, the research reviewer, on Sects. 3–8 at `ec6470e`). Its
  verdict: the structure is right and the algebra of S0 correct, but three things would let a false "hypothesis
  obtained" through: S0 with one V_HT per reflection on a curved front, two published classes missing from the
  known side (leakage, trapping), and no control of the particle number. All eight required changes are taken;
  two with a deviation, named below. The revised `exploratory-protocol.md` lists them in its Sect. 11.
- **Cadence check (R1, change 5).** `configs/check-L-cadence.yaml`, 8 GPUs, commit `2d309c8`, 600 steps in
  175 s; output `output/research/spherical-shock-injection/check-L-cadence` (4.9 GB, not in git).
  - Work of the recorded field against the energy change, per orbit with |Δγ| ≥ 0.1: read every step the
    median error is 0.4 % (all 210 orbits within 10 %); every 2 steps 0.6 % (98 % within 10 %); every 3 steps
    2.5 % (88 %); every 4 steps 5.4 % (73 %).
  - Decision: tracks every 2 steps.
- **Deviations from R1, change 1 (track stride).**
  - 947 in L and M, not 463: at every 2 steps 463 is 110 GB per run, and 291 GB are free. R1 names 947 as its
    fallback. With it the nominal prediction puts about 1 / 70 / 420 reflections on tracks (b = 2 / 3 / 4).
  - 4735 = 5 × 947 in the null run of L and in H, not the stride of L: the null holds no reflection and H
    supports no hypothesis. Its tracked electrons are every fifth one of L's.
- **Not done, against R1, change 6:** S0′ (guiding-centre orbits in the fitted front). The range of V_HT over
  the dwell is used, which R1 gives as the minimum. A residual in a band with ε > 0.3 stays "candidate".
- **Correction to the entry "16-step check of L" below:** a run does not take 45–50 minutes. A step costs
  290 ms (measured over 600 steps); the 0.90–1.06 s of the 16-step check was the cost of writing every output
  every 4 steps. A run is about 15 minutes plus the denser tracks; the null run of L will measure it.
- **Seen in the check, to keep in mind for the census:**
  - `hot_e` (γ ≥ 1.25) is dominated by debris electrons: 3.3 × 10⁶ at t = 148 against 3.5 × 10⁴ ambient ones.
    The dump levels off near 105 MB, so a run holds about 5 GB of it;
  - ambient electrons with γ ≥ 1.5 number 300–800 before any shock can be valid (t < 112). They belong to the
    start of the driver, not to a shock;
  - the work along B is about 1 % of the gain of the 210 orbits: the early gain is across the field.
- **Code state.** `analyze.py` has `front` (density front, valid domain with the numbers of the first version)
  and `cadence`; 18 of 18 self-tests pass. The R1 changes to the front measurement, and the commands for the
  floor, the census, the tracks and the records, are not written yet. They read outputs that the configs
  already write, so the null run of L starts now and the code follows while it runs.
- **Still limiting the verdict to "candidate":** the two unread predecessors (entry below), S0′, and the pair of
  runs at 2 and 8 particles per cell (made only if L shows a residual or a K1 signal).

## 2026-11-13 — Curved-shock predecessors read before the first run
- **Done:** 27 notes added to `research/literature/` (full text: 5; skimmed: 9; abstract only: 13).
  `literature.md` has a new part, "Curved fronts"; the protocol's Sect. 10 says what it means for S0, K1, K2.
- **Not obtained in full:** Krauss-Varban & Burgess 1991 (JGR 96, 143) and Giacalone 2017 (ApJ 848, 123), the
  two that matter most. The publishers' hosts refused the requests and no open copy was found. Leroy &
  Mangeney 1984 and Savoini & Lembège 1999 were not obtained at all. ADS refused requests, so the
  forward-citation search used OpenAlex only.
- **Found, against what the survey said** ("no kinetic simulation of a curved front was found"):
  - two full-PIC studies of electrons at a curved front exist, in 2-D at bow-shock parameters
    (`savoini2001twodimensional`, `savoini2010origin`; m_i/m_e = 42, M_A = 3);
  - the second shows electrons that enter near θ_Bn = 88° and leave near 48°. This is signature K1 of the
    protocol, so K1 as stated is not new;
  - the sweep of the foot-point θ_Bn in time, the source just behind the tangent line, and repeated encounters
    in closed geometries are all published.
- **Not found in what could be reached:** electrons in a 3-D spherical, slowing front in full PIC; the energy
  of an escaping electron against the one-reflection energy at its exit angle; a path across latitudes by
  drift in the ramp.
- **Consequence:** with the two key works unread, criterion C5 cannot be met; the best verdict of this
  exercise is "candidate". The owner can lift this limit with a library copy of the two papers.
- **Corrections to the survey's reference list** are kept in `literature.md` ("Vandas 1995" is two papers;
  "Vandas 2001" is two papers; Schmidt & Cairns 2014 is marked retracted in OpenAlex).
- The reading was done by a reader agent; five notes are from full texts. I checked the format of the notes,
  not each statement against the papers.

## 2026-11-13 — 16-step check of L; protocol and baseline written before the first run
- **Check run (run 0).** `configs/check-L.yaml`, 8 GPUs, commit `3e97c6c`, 16 steps in 18.8 s. All nine
  diagnostics wrote files; Gauss law 4 × 10⁻¹⁵, div B ≤ 1.3 × 10⁻¹⁶. Output:
  `output/research/spherical-shock-injection/check-L` (2.9 GB, not in git).
- **Correction to the entry of 2026-11-12 and to the plan:** a run takes about 45–50 minutes (0.90–1.06 s per
  step, measured), not 15–20. The whole run set is about 3–5 GPU-hours.
- **Changed against the plan:** the check uses `gamma_min: 1.05`, not 1.02. At 1.02 the dump would hold a
  quarter of the ambient electrons.
- **Found in the files:**
  - the loaded electrons follow a Maxwell–Jüttner distribution (T = 0.01), not a Gaussian in u: tail fraction
    at γ ≥ 1.05 measured 0.019628, Maxwell–Jüttner 0.019623, Gaussian 0.0166. The baseline uses the first;
  - `id // 2` is the birth cell of an ambient electron (checked on every electron of the dump); injected
    electrons have ids ≥ 20N; debris electrons are the species `de`, with small, varying weights;
  - reading a dump needs `import hdf5plugin` first (zstd).
- **Written:** `shockkin.py` (reflection map, foot-point formulas for an axisymmetric front, count model,
  back-tracing), `analyze.py` (`selftest`: 10 of 10 pass; `predict`), `exploratory-protocol.md`.
- **Prediction before the runs** (nominal sphere at 0.87 V_d, `data/predict-nominal.json`): one-step reflection
  gives γ ≥ 1.5 almost only at foot-point θ_Bn = 50–70°, nothing at γ ≥ 3, and a count that goes from
  5 × 10² (b = 2) to 2 × 10⁵ (b = 4) macro-electrons in the box at the end.
  - Consequence: counts cannot carry a claim. The per-electron test (does the electron trace back to a thermal
    electron reflected once?) is the primary one.
- **Not done yet:** the commands for the front, the census and the tracks; the review R1 of the protocol; no
  production run.

## 2026-11-13 — Review of the exploratory conditions (issue #48): record corrected, run set changed
- Asked by the owner: run the exploratory conditions, analyse them, and test whether a novel and certain
  hypothesis about the acceleration mechanism can be obtained; take the findings of issue #48 into the trials.
- The review is kept in `reviews/2026-11-13-exploratory-conditions.md`. Answers, point by point:
  - **A (record).** Corrected in `exploratory-conditions.md` and its Japanese version, Sects. 2.2 and 2.3:
    - each remnant value names its speed and electron temperature;
    - M_A/√β_e against the threshold of `amano2010critical` is √2 times the flow speed over the electron thermal
      speed. A remnant has about 18 at 5000 km/s and 10⁴ K and 3.2 at 3000 km/s and 10 eV; the runs have 3.5. It
      is no longer listed as matched;
    - "θ_Bn runs from 0° to 90°" now describes the initial geometry only;
    - only β of run L is listed as matched.
  - **B.1–B.5 (design).** Accepted. They go into a protocol that is committed before the first run
    (`exploratory-protocol.md`): the reflection baseline, the definition of a shock and its valid domain, the
    numerical floor from a null run, the split of energetic electrons by origin, the foot-point θ_Bn, the
    measured luminal angle, the closure test of the work, one energy threshold.
  - **C (decisions; the owner approved the plan of 2026-11-13).**
    - The record is corrected.
    - The null run of L is added as a further run; H is not replaced. The order becomes: a 16-step check, the
      null of L, L, M, then H as a contrast case only.
    - The baseline and the shock definition are written before the first run.
- The physical conditions of L, M and H are unchanged.
- Sect. 6 (checks 2 and 6) and Sect. 8 (order of the runs) of the conditions record are replaced by the protocol;
  the record keeps them as designed.

## 2026-11-12 — Three exploratory conditions designed (not run)
- Asked by the owner: two or three conditions as close to a supernova remnant as possible, for trials that
  build hypotheses: a large box, few particles, a short run. Record them in markdown; do not run yet.
- The owner chose: only what `examples/blast_wave` can express; vary the field strength only; 768³ cells with
  2 particles per cell.
- Done:
  - read the example's planner, the engine's limits and the output system;
  - designed L, M, H (nominal M_A = 3.5, 6.0, 8.7; m_i/m_e = 16; V_d = 0.25 c; v_te = 0.1 c; R_B ≈ 148 d_e);
  - wrote `configs/explore-{L,M,H}.yaml`; the planner accepts all three (`--plan-only`, CPU, plan for 8 GPUs);
  - an independent check recomputed the numbers from the planner's formulas and reviewed the output block;
  - wrote `exploratory-conditions.md` and its Japanese version.
- Found by the check, and now in the record:
  - the front runs at about 0.87 V_d and slows, so H is expected at or below the surfing threshold;
  - a track stride that shares a factor with 768 would line the tracked electrons up in planes (4733 is used);
  - the energy selection of the particle dump has no test in the repository.
- Corrected against the plan: the run time is about 15–20 minutes per run, not 7–10. The first figure counted
  the particles only, and the grid costs more here.
- Not done: no run; the electron analysis (`analyze.py`) and the run-recording tools are not written.

## 2026-11-12 — Citation check of the survey: 16 corrections
- An independent check compared about 290 located statements of Sects. 1–8 of the survey with the paper texts.
  About 280 were confirmed. Three were wrong and 13 imprecise or mis-located; all 16 are corrected in the
  English and Japanese surveys, and one wording in the summaries.
- The three wrong statements:
  - "only one paper gives Mach numbers for SNR shocks": four others quote values (M_A ≈ 20–1000);
  - "K_ep is degenerate with the magnetic field": the Tycho paper says the X-ray flux fixes K_ep independently
    of the field (the note `morlino2012strong` is corrected too);
  - "the 2-D runs last 6–52 Ω_ci⁻¹": two 2-D studies of parallel shocks run longer.
- Not checked: statements that rest on `hoshino2002nonthermal`, `amano2010critical` and `amano2022theory`
  (no text at hand), and three values that are garbled or only in figures; the survey marks them.

## 2026-11-12 — One-page summary of the survey, in English and Japanese
- Asked by the owner. `survey-snr-electron-injection-summary.md` and its `.ja.md` give the whole picture on one
  page: the problem, what observations fix, the mechanisms by θ_Bn, what decides the answer in a simulation,
  what is supported, disputed and unknown, and what it means for this topic.
- Every statement is taken from the survey, with its citekeys.

## 2026-11-12 — Japanese version of the survey
- Asked by the owner. `survey-snr-electron-injection.ja.md` has the same sections as the English survey, which
  stays the original.
- Its reference lists (Sect. 10) are taken from the English file by script; only headings and remarks are
  translated.
- Checked by script, section by section: the citekeys are identical, and no number of the English version is
  missing.

## 2026-11-12 — Survey: electron injection at supernova remnant shocks
- Asked by the owner: survey the prior work on the electron injection problem at SNR shocks, in markdown, with
  the citations kept.
- Done:
  - arXiv searches by keyword (12 queries, 238 new results), then by author and title, then reference lists;
  - 54 works selected, full texts from arXiv, read by nine readers in groups (observations; surfing line;
    perpendicular SNR series; foreshock and reviews; quasi-parallel; thresholds 2021–2026; shock drift theory;
    geometry);
  - one note per work in `research/literature/`;
  - the headline numbers of 20 works checked against the paper texts;
  - the survey written: `survey-snr-electron-injection.md`.
- Main findings:
  - mechanisms are sorted by obliquity, and each is known from planar runs at one or a few angles;
  - no study covers the range of θ_Bn with one parameter set;
  - global models of remnants assume the efficiency law ε(θ_Bn);
  - SN 1006 points to quasi-parallel acceleration and the radio population to quasi-perpendicular;
  - no kinetic simulation of a curved front at SNR-like parameters was found.
- Not done:
  - the bow-shock works on curved shocks (not on arXiv) are not read, so the last finding is not certain;
  - relativistic shocks and wakefield acceleration were not covered;
  - the survey has not been reviewed by the research reviewer.
- Corrected along the way: the reference of `amano2022theory` (ApJ 927, 132).
- No decision on the research question; it stays with the owner.

## 2026-11-12 — First review (purpose stage)
- Review: `reviews/2026-11-12-purpose.md`. Verdict:
  - "every mechanism over all θ_Bn" needs rethinking (a programme, not one question);
  - "curvature-only channels" needs the controlling number R/L_s and a target system first;
  - "the switch-on angle θ_c and its scaling with m_i/m_e and β_e" is answerable.
- Fixed in the README: the Background said that known mechanisms are framed only for planar
  shocks. amano2022theory Sect. 4.5 (Eq. 22) already has a curvature limit.
- Added to `literature.md`: the curved-shock works the review found by search (not read).
- Points to keep whatever the question becomes:
  - the mass ratio moves every threshold, with different powers;
  - θ_Bn is to be taken from the measured front normal;
  - the superluminal band;
  - Buneman waves against resolution and the current filter;
  - tail statistics per θ_Bn band;
  - tracked sets fixed at the start.
- Open: the owner's decision on the research question.

## 2026-11-12 — Topic opened
- The owner set the purpose: the electron injection problem, studied with PIC of a spherical
  shock, considering every acceleration mechanism, unknown ones included.
  - Mechanisms named: surfing, drift, wakefield, first- and second-order Fermi, combinations,
    geometric effects of the curved surface, and the shape changes as θ_Bn goes from
    perpendicular to parallel.
- Branch cut from research/trunk at 6f00c11.
- README: purpose, background and a draft research question. Prior work lists the four notes
  read so far; the survey list is in `literature.md`.
- Open: the survey, a narrowed question, the idea, and the driver and its parameters.
  `examples/blast_wave` (M_A = 3) is likely below the Mach numbers that surfing needs.
