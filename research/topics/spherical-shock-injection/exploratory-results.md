# Exploratory runs: results and verdict

Topic `spherical-shock-injection`. Japanese version: `exploratory-results.ja.md` (this English file is the
original). Written 2026-11-13, after the runs L-null, L, M, M-null and H and their analyses. Revised the same
day after review R2 (`reviews/2026-11-13-results-r2.md`); Sect. 10 lists what changed.
Judged by `exploratory-protocol.md`, which was fixed before the runs. No registered criterion was changed. Two
registered items were not applied as registered in the first version of this document: the test N3 and the
replacement rule of τ_dwell. Both are applied here (Sect. 9). The verdict is the same in both versions.
Detail and failed attempts: `notebook.md` (entries of 2026-11-13; the entry on review R2 corrects the earlier
ones). Records: `runs/` (`index.json`); the checks added after R2 are the blocks `r2`, `census-tau` (L, M) and
`r2null` (null runs) of the records.
Units: c = 1, m_e = 1, lengths in d_e, times in 1/ω_pe; d_i = 4 d_e; Ω_ci⁻¹ = 56 (L), 96 (M), 139 (H).

## 1. Verdict

**Question (protocol Sect. 1).** Do runs L, M, H, analysed, give a research hypothesis about the electron
acceleration mechanism that is novel and certain enough to build a study on?

**Answer: not obtained from these runs.**

- **No signal at the registered threshold.** No band of any run holds 100 ambient electrons with γ ≥ 1.5
  upstream of a front. Such electrons exist, but they are few: 22 electron-dumps in L (at most 3 per dump, in
  eleven dumps; 3 of the 22 in the valid domain) and 22 in M (at most 5 per dump, all with feet in bins without
  a detected front). The null runs hold none. The residual is zero because these electrons are few and all are
  classified as explained by S0 or undetermined, not because there are none (Sect. 5). There is nothing to
  build a hypothesis on.
- **The threshold lies far from the realised shocks.** γ − 1 = 0.5 is 45 (L) and 89 (M) times ½ m_e V_n² at the
  median front speed, and 48 T_e. One reflection brings a cold electron there only for V_HT ≥ 0.447. At the
  median V_n that is θ_Bn ≥ 70° (L) and 76° (M), next to the superluminal band (81°, 84°). The registered test
  therefore probed a narrow band of these fronts.
- **The runs did not sample the designed shocks.** The design had fast Mach numbers 2.95 (L) and 4.05 (M).
  Where the registered definition finds a shock, the fast Mach number is 1.5–2.5 (median 1.8) and the
  compression is b ≈ 2 in both runs. The lower end 1.5 is the cut of the definition, not a property of the
  fronts. The debris stalls at 1.3–1.5 R_m, the front satisfies the definition only from about that time on,
  and it decelerates as R ∝ t^0.55 (L) and t^0.52 (M): a blast wave that is no longer driven. H forms no shock
  in 5.3 Ω_ci⁻¹.
- **What the absence of energetic electrons says.** The first version said that the absence "is what such weak
  shocks give". That is withdrawn. The registered baseline (S0, an adiabatic mirror in a static potential)
  predicts a beam of 10²–10³ electrons at γ ≥ 1.5 for these same weak fronts, and the beam is absent. But in the
  valid domain 98 % (L) and 100 % (M) of the predicted beam comes from the bands θ_Bn = 60–80°, where the
  protocol's own measure ε is above 0.3 and the static baseline is not defined (Sect. 6). The numerical
  scattering of thermal electrons is also above the registered limit (Sect. 3). So the absence constrains
  neither the designed nor the realised Mach numbers. It is a statement about the baseline and the numerical
  floor, not about the plasma.
- **Returns beyond the baseline (not registered; exploratory).** The mirror model returns 0–6 % of the
  encounters where the tracked orbits return 1–15 %, at all energies. By the registered test N3, no loss of
  magnetic moment can be assigned to physics in these runs. The share of the returns that keep μ within 20 %
  (12 %) equals the share of electrons that keep it over the same time in the null runs (11 %; Sect. 3). This
  is recorded as an open observation, not as a candidate.
- **Why not "candidate".** The protocol asks for "a residual with a positive signature". None exists: K1 has
  one electron in L and two in M against the 100 required, none of them in the valid domain; K2 and K3 have
  none.
- **What the answer does not say.** It does not say that a curved shock has no mechanism of its own. Two
  particles per cell, m_i/m_e = 16 and a front of fast Mach 1.8 cannot test that.

What a further study would need is in Sect. 8.

## 2. Runs

All on 8 GPUs from a clean committed tree, seed 12345, 768³ cells, dx = 0.5, 2960 steps (t = 740), 2 particles
per cell, m_i/m_e = 16, V_d = 0.25. Conservation contracts hold in every run (Gauss law ≤ 5.2 × 10⁻¹⁵, div B ≤
1.9 × 10⁻¹⁴, fill ≤ 0.674).

| Run | Config | Commit | Wall time | Output | Energy budget | Record |
|---|---|---|---|---|---|---|
| L-null | `explore-L-null.yaml` | `0114ef9` | 778 s | 15.5 GB | +12.7 % | `2026-11-13-0114ef9ea-L-null.json` |
| L | `explore-L.yaml` | `91f79ff` | 1298 s | 61.8 GB | +10.4 % | `2026-11-13-91f79ff1d-L.json` |
| M | `explore-M.yaml` | `ac33a41` | 1321 s | 60.2 GB | +11.8 % | `2026-11-13-ac33a41eb-M.json` |
| M-null | `explore-M-null.yaml` | `3f5c437` | 783 s | 15.8 GB | +12.8 % | `2026-11-13-3f5c43715-M-null.json` |
| H | `explore-H.yaml` | `3f5c437` | 814 s | 16.9 GB | +12.3 % | `2026-11-13-3f5c43715-H.json` |

- Design: L has M_A = 3.5, R_m = 64 d_e = 16 d_i. M has M_A = 6.0, R_m = 44.8 d_e = 11.2 d_i. H has M_A = 8.7,
  R_m = 34.8 d_e = 8.7 d_i. R_B ≈ 148 d_e in all three.
- The gyro-radius of a debris ion in B₀ is ρ_d = M_A d_i = 14 (L), 24 (M), 34.8 (H) d_e, so R_m/ρ_d = 4.6, 1.9
  and 1.0.
- Outputs are on the NVMe, `/mnt/nvme/pypic/output/research/spherical-shock-injection/<name>`, linked from
  `output/research/spherical-shock-injection/<name>`; not in git.
- Not run: the second seed L-s2 and the 384³ pair N6. The protocol (Sect. 9) runs them only if L has a
  residual or a K1 signal.

## 3. Numerical floor (null runs)

| Measure | L-null | M-null |
|---|---|---|
| N1: electrons with γ ≥ 1.5, any dump | 0 (largest γ 1.374) | 0 |
| N2: T_e, start → end | 0.0100 → 0.0110 | 0.0100 → 0.0109 |
| N3: ν of thermal electrons (v 0.1–0.2), lags 4 / 16 / 64 | 1.2 × 10⁻² / 7.8 × 10⁻³ / 5.0 × 10⁻³ | 1.3 × 10⁻² / 7.9 × 10⁻³ / 5.0 × 10⁻³ |
| N3: ν for v > 0.3 | 2.5 × 10⁻³ / 1.8 × 10⁻³ / 1.3 × 10⁻³ | 3.3 × 10⁻³ / 2.0 × 10⁻³ / 1.4 × 10⁻³ |
| N5: 99.9th percentile of \|Δρ\|, start → end | 1.26 → 3.71 d_e | 2.07 → 5.96 d_e |
| Fronts detected | 0 of 738 bins | 0 of 738 bins |

- With a floor of 0 the count criterion of a signal is the fixed one: 100 electrons or more.
- The energy budget grows by 10–13 % in every run, with or without debris: numerical heating of the electrons
  at 2 particles per cell.
- N4 (closure of the work on tracked orbits) passes in L (1289 of 1341 orbits within 10 %) and in M (851 of
  853).

**N3 as registered (added after R2).** The protocol (Sect. 5) tests ν τ < 0.1 for thermal electrons, with τ the
crossing time of the layer, its measured width over V_n. The width of the ramp of |B| (10–90 % of the rise) was
measured in every valid bin: 2.6 / 3.6 / 4.8 d_e in L and 4.8 / 5.4 / 6.0 d_e in M (10th / 50th / 90th
percentile). These are upper bounds: the fit error of the surface inside a bin and the average over eight
half-planes widen the ramp. The ranges below span the three lags of ν.

| ν τ of thermal electrons, with τ = | τ (L / M) | L | M |
|---|---|---|---|
| crossing time of the ramp, median of the valid bins | 23.9 / 49.0 | 0.12–0.28 | 0.25–0.64 |
| the same, medians of the θ_Bn bands | | 0.10–0.37 | 0.24–0.66 |
| crossing time of the census layer (8 d_e) | 53.5 / 75.4 | 0.27–0.62 | 0.38–0.98 |
| τ_dwell assumed by the protocol for B1′ (3 gyro-periods) | 66 / 113 | 0.33–0.76 | 0.57–1.47 |
| median dwell of the tracked returns | 250 / 252 | 1.3–2.9 | 1.3–3.3 |
| 90th percentile of that dwell | 417.5 / 414.5 | 2.1–4.8 | 2.1–5.4 |

- **N3 fails in every reading.** Only the smallest value of L (the slowest rate in the thinnest bands) touches
  the limit 0.1. The registered consequences apply: the loss cone is being filled numerically, the baseline B2
  carries that note (Sect. 6), and no mechanism that needs scattering is assigned.
- The first version of this document applied N3 to the electrons of γ 1.25–1.5 and to the τ_dwell of B1′, and
  called the result at 66 and 113 a pass. That was not the registered test (Sect. 9).

**Retention of the magnetic moment in a uniform field (added after R2; not registered).** Share of the tracked
electrons of the null runs whose μ, averaged over one gyro-period, is within 20 % of its earlier value after a
gap:

| Gap | L-null (167 026 electrons) | Gap | M-null (165 498 electrons) |
|---|---|---|---|
| 0 (adjacent windows) | 27 % | 0 | 23 % |
| 66 (3 gyro-periods) | 14 % | 113 (3 gyro-periods) | 12.5 % |
| 150 | 12 % | 150 | 12 % |
| 250 | 11 % | 250 | 11 % |
| 420 | 11 % | 420 | 11 % |

- Thermal electrons alone (v 0.1–0.2) give the same shares within one point; faster ones keep μ slightly more
  often (13 % for v > 0.3 at a gap of 250 in L-null).
- The ratio of μ after to μ before spreads from 0.15 to 7.5 (10th to 90th percentile) at a gap of 250 in
  L-null. At two particles per cell the magnetic moment of a thermal electron is not conserved over the dwell
  times of these runs, with or without a shock. The 20 % test of μ therefore cannot show a physical loss of μ.

## 4. The front and its valid domain

| | L | M | H |
|---|---|---|---|
| Bins (time × θ) with a front | 443 of 738 | 326 of 738 | 0 of 738 |
| Valid bins (protocol Sect. 4) | 241 | 180 | 0 |
| Valid bins usable for B1′ | 166 | 121 | 0 |
| Valid domain | θ = 30–90°, t = 185–666 | θ = 45–90°, t = 259–703 | none |
| The same in Ω_ci⁻¹ | 3.3–11.9 | 2.7–7.3 | – |
| Fast Mach, valid domain (min / median / max) | 1.50 / 1.82 / 2.51 | 1.50 / 1.75 / 2.38 | – |
| Fast Mach, design (piston) | 2.95 | 4.05 | 4.64 |
| V_n, valid domain | 0.116–0.204 | 0.089–0.147 | – |
| Compression b, valid domain (median) | 1.58–2.78 (2.0) | 1.62–2.50 (1.94) | – |
| Debris edge (0.1 n₀) stalls at | 83 d_e = 1.3 R_m | 62–66 d_e = 1.4–1.5 R_m | ≈ 50 d_e |
| Superluminal feet | θ_Bn ≥ 77–83° | θ_Bn ≥ 78–85° | – |

The minimum 1.50 of the fast Mach number is the cut of criterion (c) of the definition. Next to the equator
the front leaves the valid domain at the end because it falls below that cut (L from t = 666, M from t = 684.5).

**Regime of the valid domain (added after R2).** Percentiles are 10th / 50th / 90th over the valid bins.

| | L | M |
|---|---|---|
| V_n | 0.127 / 0.150 / 0.185 | 0.095 / 0.106 / 0.130 |
| Alfvén Mach number | 1.78 / 2.09 / 2.59 | 2.27 / 2.55 / 3.11 |
| Fast Mach number | 1.55 / 1.82 / 2.24 | 1.57 / 1.75 / 2.13 |
| θ_Bn (min / median / max) | 20.9° / 63.6° / 87.1° | 42.8° / 69.0° / 87.6° |
| Upstream T_e, T_i | 0.0105, 0.0099 | 0.0105, 0.0099 |
| β_e, β | 0.26, 0.50 | 0.76, 1.47 |
| V_n / v_te | 1.23 / 1.46 / 1.83 | 0.91 / 1.04 / 1.28 |
| (m_i/m_e) / M_A (median M_A) | 7.6 | 6.3 |
| γ − 1 = 0.5 in units of ½ m_e V_n², ½ m_i V_n², T_e | 45, 2.8, 48 | 89, 5.6, 48 |
| Compression b | 1.71 / 2.00 / 2.58 | 1.74 / 1.93 / 2.27 |
| b at θ_Bn ≥ 80° (median) and its Alfvén Mach number | 2.13 at 2.15 (52 bins) | 2.06 at 2.60 (48 bins) |
| Perpendicular MHD jump there, adiabatic index 5/3, 2 | 1.98, 1.82 | 1.94, 1.73 |
| Width of the ramp, 10–90 % (upper bound) | 2.6 / 3.6 / 4.8 d_e | 4.8 / 5.4 / 6.0 d_e |
| Exponent m of R ∝ t^m | 0.55 | 0.52 |

- **Both runs hold the same fast Mach number, not the same shock.** The Mach scan of the design (L against M)
  did not happen: the medians are 1.82 and 1.75. The two valid domains still differ in β (0.50, 1.47), in
  V_n / v_te (1.5, 1.0), in the Alfvén Mach number (2.1, 2.6) and in θ_Bn (M has no valid front below 43°).
- **Compression.** The measured b is the largest mean |B| behind the front, so it includes any overshoot. It
  lies at or slightly above the perpendicular MHD jump for the measured Mach number and β.
- **Critical Mach number.** The first critical Mach number at these β and θ_Bn was not computed. Whether the
  valid fronts are subcritical or supercritical is therefore not known.
- **No shock towards the pole.** No valid front exists below θ = 30° (L) and 45° (M). The quasi-parallel part
  of the sphere is not sampled.
- **H:** the ambient ion density reaches 1.5 n₀ only in 3 equatorial bins at the end of the run.

**Where the front first becomes valid.** Next to the equator (θ = 87.5°):

| First valid front | L | M |
|---|---|---|
| Time | 185 = 3.3 Ω_ci⁻¹ | 259 = 2.7 Ω_ci⁻¹ |
| Radius | 73.5 d_e = 18.4 d_i | 71.6 d_e = 17.9 d_i |
| Radius in R_B, in ρ_d, in R_m | 0.50, 5.3, 1.15 | 0.48, 3.0, 1.60 |
| All polar angles: time, radius | 185–277.5, 73–86 d_e | 259–333, 71–80 d_e |

The first valid front has the same radius in d_i, the same fraction of R_B and a similar time in Ω_ci⁻¹ in
both runs. In ρ_d and in R_m it differs by a factor 1.4–1.8. Two runs cannot choose between these readings
(Sect. 8, item 4).

**Bins with a front before the valid domain (corrected after R2).** Counted are the bins that show a front
before the first valid time of their polar angle, at the polar angles that later become valid: 74 in L and 73
in M. (The first version gave 62 and 54 and said that all fail criterion (b) only. Neither number is reproduced
by the records, and the statement was wrong.)

| | L (74 bins) | M (73 bins) |
|---|---|---|
| Fail (b), lead over the debris | 74 | 73 |
| Fail (c), fast Mach ≥ 1.5 | 24 | 18 |
| Fail (e) | 28 | 35 |
| Fail (b) alone | 46 | 38 |
| Fast Mach of those that pass (c): min / median / max | 1.92 / 2.50 / 2.78 (50 bins) | 2.26 / 2.85 / 3.28 (55 bins) |
| Lead required by (b), 0.7 V_n/Ω_ci: min / median / max | 1.5 / 2.0 / 2.3 d_i | 2.3 / 2.8 / 3.3 d_i |
| Lead measured: 10th / 50th / 90th percentile | −0.8 / 0.9 / 1.7 d_i | −2.1 / −1.4 / 1.4 d_i |

- Criterion (b) asks for a lead of 0.7 V_n/Ω_ci over the debris edge, not for a fixed 2 d_i. In M that is
  2.3–3.3 d_i.
- A phase of higher fast Mach number exists (medians 2.5 and 2.85), but in it the compression has not
  separated from the piston. In M the front is mostly still inside the debris edge (negative lead).
- Next to the equator: in L the front is behind the debris edge until t = 92.5, fails (b) alone at
  t = 129.5–166.5 and is valid from t = 185. In M the front is behind the debris edge until t = 203.5, fails
  (b) alone at t = 203.5–240.5 and is valid from t = 259.
- **Criterion (b) late in M.** The debris edge exists at all three density thresholds in every output of M; at
  t = 444 it is at 63 d_e, the front is at 95 d_e, and the lead is 32 d_e against 7.5 d_e required. The dashes
  for t ≥ 444 in the equator table of the notebook were a gap of that table, not of the criterion.

**Reports asked for by protocol Sect. 4 (added after R2).**

| | L | M |
|---|---|---|
| Density front minus mid-ramp surface of \|B\| (10th / 50th / 90th) | −0.8 / 0.1 / 1.0 d_e | −0.8 / −0.1 / 0.6 d_e |
| Debris edge at 0.3 n₀ minus at 0.03 n₀ | −5.8 / −3.7 / −1.9 d_e | −14.4 / −11.2 / −6.5 d_e |
| Debris edge at 0.1 n₀ minus at 0.03 n₀ | −2.1 / −1.5 / −0.6 d_e | −6.4 / −3.8 / −2.3 d_e |
| Valid bins without a debris edge (three thresholds) | 0, 0, 0 | 0, 0, 0 |
| Front minus debris edge (0.1 n₀) | 2.7 / 6.3 / 12.8 d_i | 4.0 / 7.8 / 12.4 d_i |
| Equatorial plane: density front minus B_z ≥ 1.25 B₀ (median) | −1.6 d_e (37 outputs) | −1.9 d_e (36 outputs) |
| Rms of the fit: density front, mid-ramp surface | 0.63, 0.51 d_e | 0.29, 0.31 d_e |

The debris edge of M is diffuse: its position moves by 11 d_e between the thresholds 0.03 and 0.3 n₀. The
valid domain of M depends on that threshold more than the one of L.

**The front decays.** At the equator the fast Mach number falls from 2.2 (t = 222) to 1.6 (t = 592) in L and
from 2.2 (t = 296) to 1.6 (t = 592) in M (table in `notebook.md`). V_n falls from 0.202 to 0.130 (L) and from
0.147 to 0.099 (M) inside the valid domain. The exponents 0.55 and 0.52 lie between a driven front (1) and
the Sedov value (0.4); no model was fitted.

## 5. Census at the registered threshold (γ ≥ 1.5)

The L column of this table was wrong in the first version (486 at t = 74 as the largest dump, "≤ 40, then ≤ 4"
in the valid domain, 3 upstream at one dump). The values below are from `census.json` of the records.

| | L | M | H |
|---|---|---|---|
| Ambient electrons, all classes, largest dump | 807 (t = 92.5) | 33 (t = 92.5) | 1 (t = 111) |
| The same in Ω_ci⁻¹ | 1.65 | 0.96 | 0.80 |
| At the dumps from the first valid time on | 110, 59, 40, 21, 13, 10, then ≤ 6 | 1 in each of 5 dumps | – |
| Electron-dumps there / different electrons | 279 / 199 | 5 / 3 | – |
| Upstream of the fitted surface, electron-dumps | 22, in eleven dumps | 22, in ten dumps | – |
| …largest dump | 3 (five dumps, t = 166.5–277.5) | 5 (t = 111) | – |
| …with a foot on a detected front | 16 | 0 | – |
| …in the valid domain | 3 | 0 | – |
| …explained by S0 / undetermined / residual | 5 / 17 / 0 | 0 / 22 / 0 | – |
| Bands with ≥ 100 upstream electrons | 0 | 0 | – |
| K1: electron-dumps / different electrons (needed: 100) | 7 / 1 | 3 / 2 | – |
| Tracked orbits with γ_max ≥ 1.5 | 11, all starting behind the front | 0 | – |
| K2, K3 on tracks | 0 | 0 | – |

H has no census; its number is from the `hot_e` output.

- **Classes, all dumps.** L: 3035 electron-dumps in the layer, 338 downstream, 208 before any front is fitted,
  22 upstream. M: 54 in the layer, 52 before any front, 22 upstream.
- **The energetic electrons are made by the piston, before any valid front.** The largest dump of L (807) is at
  t = 92.5, where the front is still behind the debris edge at the equator; 805 of the 807 are in the layer
  class. This is an effect of the piston stage, not a result about a shock. Electron pre-energisation by a
  laser-driven piston is reported in full PIC by `lezhnin2021kinetic` (skimmed, not read in full).
- **Their number falls as the field falls, in every normalisation.** At their largest dumps: 807 (L), 33 (M),
  1 (H), for B₀ = 0.286, 0.167, 0.115. At equal Ω_ci t ≈ 1.3: 486 (L, t = 74) and 18 (M, t = 129.5); H has no
  census at the matching time. The first version compared the runs at equal t = 74 and described the trend
  with the wrong sign ("falls with the magnetisation").
- **Why the population of L decays.** The 807 electrons of t = 92.5 were followed by their ids: 220 are still
  above γ = 1.5 at t = 111, 44 at t = 185, none at t = 296. Meanwhile their median distance behind the front
  grows from 6 d_e to 27 d_e, and the median γ of those still above 1.25 falls from 1.55 to 1.31. They are left
  behind the front and lose energy; the cause of the loss was not analysed. They are not replaced: of the 199
  electrons above 1.5 in the dumps from t = 185 on, 107 had been above 1.5 earlier and 163 above 1.25.
- **K1.** The seven records of L are one electron (id 391104754) at t = 166.5–277.5, with γ falling from 1.84
  to 1.50. None of the seven lies in the valid domain and none outside the debris edge, so it cannot be
  separated from the piston. Whether its field line threads the debris was not traced. The three records of M
  are two electrons at t = 111–166.5 (γ 1.55–1.57), before the valid domain.
- A missing K3 is not a result (protocol Sect. 7): the current filter damps the Buneman scale.
- **Tracked returns with γ_out ≥ 1.25:** 7 (L), 3 (M). Both are below the 30 the protocol needs, so neither run
  can test C1 or C2 on orbits. The sample for the false-negative rate of B1′ is 0 (L) and 1 (M): that rate is
  not measurable.

**The replacement rule of τ_dwell (protocol Sect. 3.2; applied after R2).** The 90th percentile of the dwell of
the tracked reflections with γ_out ≥ 1.25 is 198 (L, 7 reflections) and 286.5 (M, 3 reflections), above the
assumed 66 and 113. The rule then replaces the assumed value and both results are reported. Upstream
electron-dumps by result:

| | L, τ = 66 | L, τ = 198 | M, τ = 113 | M, τ = 286.5 |
|---|---|---|---|---|
| γ ≥ 1.5 (registered): number | 22 | 22 | 22 | 22 |
| …explained by S0 | 5 | 8 | 0 | 0 |
| …undetermined | 17 | 14 | 22 | 22 |
| …residual | 0 | 0 | 0 | 0 |
| γ ≥ 1.25 (extra bin): number | 35 756 | 35 756 | 13 433 | 13 433 |
| …explained by S0 | 19 465 | 20 917 | 4277 | 4479 |
| …S1 leak | 27 | 24 | 0 | 0 |
| …undetermined | 13 820 | 13 728 | 8133 | 8127 |
| …residual with an encounter | 2103 | 746 | 316 | 120 |
| …residual without an encounter | 341 | 341 | 707 | 707 |
| …of these, valid domain: residual with an encounter | 702 of 10 155 | 150 | 111 of 2779 | 55 |

- At the registered threshold the rule changes nothing: no residual with either τ.
- The bin γ 1.25–1.5 is not registered. Its "residual" shrinks by a factor 2–5 when the measured dwell replaces
  the assumed one. What remains sits in bands with ε > 0.3 and under a failed N3, so it supports no claim.

**The downstream branch of C1 (added after R2).** Electrons above the energy that corresponds to γ = 1.5 behind
the front, in the slab |y| < 8 d_e, against 3 times the thermal tail at the upper fit temperature. L: one
macro-electron (weight 0.0625) in 16 outputs, against 6.2 in the same units for three times the tail; no
output holds 100 electrons. M: none in 14 outputs. Outside the debris edge there is none in either run. The
downstream branch fails the count as well.

## 6. The baseline against the runs

The baseline B2 (protocol Sect. 3.4) is the adiabatic mirror with the measured front, between an upper bound
(φ = 0) and a lower bound (measured φ). By the registered N3 it carries the note: **the loss cone is being
filled numerically** (thermal ν τ > 0.1, Sect. 3).

**ε per band (protocol Sect. 3.1; added after R2).** ε = V_n tan²θ_Bn τ_dwell / R is the relative change of
V_HT while the electron dwells. Medians of the valid bins per band of θ_Bn, for four values of τ_dwell: the
assumed one, the one set by the replacement rule, and the median and 90th percentile of all tracked returns.

| L: θ_Bn | Bins | V_HT median / max | ε at τ = 66 | 198 | 250 | 417.5 |
|---|---|---|---|---|---|---|
| 20–30° | 12 | 0.17 / 0.20 | 0.02 | 0.07 | 0.09 | 0.16 |
| 30–40° | 28 | 0.17 / 0.25 | 0.05 | 0.14 | 0.18 | 0.29 |
| 40–50° | 28 | 0.20 / 0.28 | 0.09 | 0.28 | 0.36 | 0.60 |
| 50–60° | 42 | 0.26 / 0.35 | 0.19 | 0.58 | 0.74 | 1.2 |
| 60–70° | 40 | 0.37 / 0.49 | 0.45 | 1.4 | 1.7 | 2.9 |
| 70–80° | 39 | 0.57 / 0.85 | 1.1 | 3.4 | 4.3 | 7.2 |
| 80–90° | 52 | 1.7 / 3.7 | 11.5 | 35 | 44 | 73 |

| M: θ_Bn | Bins | V_HT median / max | ε at τ = 113 | 286.5 | 252 | 414.5 |
|---|---|---|---|---|---|---|
| 40–50° | 20 | 0.14 / 0.20 | 0.13 | 0.34 | 0.30 | 0.49 |
| 50–60° | 40 | 0.18 / 0.24 | 0.27 | 0.69 | 0.60 | 0.99 |
| 60–70° | 40 | 0.26 / 0.35 | 0.62 | 1.6 | 1.4 | 2.3 |
| 70–80° | 32 | 0.43 / 0.60 | 1.8 | 4.6 | 4.1 | 6.7 |
| 80–90° | 48 | 0.99 / 3.4 | 10.6 | 27 | 24 | 39 |

- Valid bins with ε ≤ 0.3: 111 of 241 (L) and 47 of 180 (M) at the assumed τ_dwell; 56 and 7 at the one set by
  the rule; 26 and 0 at the 90th percentile of the returns.
- Valid bins where one reflection brings a cold electron to γ ≥ 1.5 (V_HT ≥ 0.447): 96 in L (θ_Bn 65–87°) and
  59 in M (θ_Bn 72–88°). In all of them ε is above 0.3 at the assumed τ_dwell: at least 0.69 (median 2.6) in L
  and at least 2.1 (median 6.2) in M.
- **The predicted beam lies where the baseline is not defined.** In the valid domain the bands 60–80° give
  98–99 % (L) and 100 % (M) of the predicted electrons with γ ≥ 1.5; the rest of L (0.5–1.5 %) is the band
  50–60°. The protocol reads a residual in a band with ε > 0.3 as at most "candidate" because V_HT changes by
  that fraction during one dwell. The same holds for a deficit. In the bands 60–80° ε is 0.45–1.8 at the
  assumed dwell and 1.4–7 at the measured ones; there a static V_HT per foot predicts neither the energy nor
  the number of the reflected electrons.

**Potential.** φ is measured (error below 50 %) at polar angles up to 67.5°: 0.015–0.048 (1.5–4.7 T_e) in L and
0.014–0.023 (1.4–2.3 T_e) in M. In M the bins at 42.5° and 57.5° are not measured (the first version quoted
the value of the bin at 57.5°, 0.009, as the lower end). Above 70° φ is not measured in either run.

**γ ≥ 1.5.** The model predicts on the detected front 740–3890 electrons in L (t = 185–407) and 107–1046 in M
(t = 130–407). The measured upstream electrons with a foot band on the detected front number at most 3 per
dump in L (16 electron-dumps in nine dumps between t = 166.5 and 333; 5 of them explained by S0) and 0 at
every dump in M. (The "2, 0, 0" of the first version were the values of L at the three display times
t = 185, 296 and 407.) No cell holds 100 measured electrons, so the rule of Sect. 3.4 has nothing to judge.
The sensitivity of the predicted count to the errors of V_n and θ_Bn was not computed. The count comes from
the far tail of the distribution in V_HT and depends on these errors exponentially, so the prediction is known
to no better than its order of magnitude.

**γ 1.25–1.5 (extra bin, not registered, exploratory only).** In the cells with ≥ 100 electrons, the measured
count is outside [lower/2, 2 × upper] in 50 of 86 (L) and 9 of 13 (M). The measured electrons sit at low θ_Bn,
where a mirror gives almost nothing; at 50–80°, where it should give the most, they are missing by a factor of
200 (L) to 500 (M) or more. By the rule of Sect. 3.4 no residual is discussed.

**Share of the encounters that end in a return, all energies (not registered, exploratory only)**; tracked
orbits, detected front:

| | L (θ_Bn 10–90°) | M (θ_Bn 30–90°) |
|---|---|---|
| Model, φ = 0 | 0–13 % | 0–13.5 % |
| Model, measured φ | 0.07–0.5 % | 0–6 % |
| Tracks | 1–9 % | 7–15 % |
| Tracks, only returns from behind the density front | 0.6–6 % | 4–7 % |
| Returns that keep μ within 20 % | 11.8 % (381 of 3224) | 12.3 % (488 of 3967) |
| Null run: electrons that keep μ within 20 % over 250 | 11 % | 11 % |
| Upstream Maxwellian: share faster than V_n along the normal | 7 % | 15 % |
| Foot θ_Bn at the deepest point of the returns (10th / 50th / 90th) | 20° / 37° / 54° | 22° / 45° / 61° |

- At θ_Bn = 70–90° the mirror returns almost nothing with or without a potential; the tracks return 1–2 % (L)
  and 7–8 % (M).
- The median dwell of the returns is 250 (L) and 252 (M): 11 and 7 gyro-periods. For comparison,
  `savoini2010origin` reports about 2.7 gyro-periods for mirrored electrons and about 40 for trapped ones.
- That work finds 4.7 % of the interacting electrons backstreaming, at M_A = 3 and m_i/m_e = 42. The shares
  here are of the same order. The definitions differ (a share per band of tracked encounters here), so this is
  not a test.
- The θ_Bn at the exit of a return is not in the records; only the foot angle at the deepest point is. The
  comparison of exit and entry angles asked for by R2 was not made.

**Reading (changed after R2).** The first version read: "most returns are not adiabatic mirror reflections in a
static potential". The data do not support a statement about the plasma. The returns keep μ no more often
than electrons in a uniform field do over the same time, so the test of μ shows the numerical floor. What
stands: the static mirror baseline does not describe the returns of these runs, and these runs cannot say why.

**Explanations that these runs cannot separate:**

1. *Numerical scattering.* Thermal ν τ is 0.10–0.66 over the crossing time of the measured ramp and 1.3–5.4
   over the measured dwell (Sect. 3), above the limit 0.1 of N3 in every reading. By the protocol the loss cone
   is filled numerically and no loss of μ may be assigned to a physical mechanism. This explanation alone is
   sufficient for the loss of μ; it is not shown to be sufficient for the number of returns.
2. *Marginal adiabaticity by design.* Ω_ce d_i / V_n = (m_i/m_e) v_A / V_n = 7.6 (L) and 6.3 (M). The measured
   ramp (an upper bound) is crossed in 24 (L) and 49 (M), about 1.1 and 1.3 gyro-periods.
3. *Time dependence and curvature.* ε is above 0.3 in most of the valid domain at the measured dwell; the
   model is static per foot.
4. *Thermal escape and leakage (S1), known.* V_n / v_te is 1.5 (L) and 1.0 (M): the front does not outrun the
   thermal electrons. The foot of a field line slows from V_HT = ∞ at the tangent point towards V_n, and the
   heated electrons behind the front can overtake it. This fits the excess at low θ_Bn and the small gain
   (median 0.02–0.04). Published for a curved front (`savoini2010origin`). A count model per band from the
   measured downstream distribution was not made; it is the next analysis without GPU time (Sect. 8).
5. *Kinematic returns.* Half of the "back upstream" orbits never reach the density front.
6. *Reduced gain at a curved front, known* (`vandas1994acceleration`: about 1/4 of the plane-shock gain).

## 7. Criteria C1–C6

| Criterion | Energetic residual (γ ≥ 1.5) | K1 (Δρ) | Returns beyond S0 (not registered) |
|---|---|---|---|
| C1 real: ≥ 100, ≥ 10 × floor, second seed | fails: at most 3 (L), 5 (M) upstream per dump; downstream branch fails | fails: 1 electron (L), 2 (M) | not a residual; one seed |
| C2 mechanism on orbits | not testable: 7 and 3 tracked returns | – | no: μ is lost at the rate of the null runs |
| C3 not numerical, not the driver | – | fails: not in the valid domain, inside the debris edge | fails: thermal ν τ > 0.1 in every reading; N6 not run |
| C4 testable in dimensionless numbers | – | – | not stated |
| C5 new to the field | – | – | no: items 4 and 6 are published |
| C6 review R2 | verdict stands (Sect. 10) | stands | "a statement about the baseline and the floor" |

K2, K3 and K4 have no electron and no signature. Independently of the above, C5 could not be met: the two
predecessors `kraussvarban1991electron` and `giacalone2017acceleration` are still unread (protocol Sect. 10).

## 8. What a further study would need

Rewritten after R2. The arithmetic in this section is ours and is not a result of the runs.

**What is missing is not only computing size.** The first version closed with "the difficulty is computing
size, not a missing idea". That overreached. Three things are missing before any size matters:

- **A narrowed question.** The exploratory question ("is there a mechanism of its own at a curved front?")
  has no controlling number. It is the owner's decision (open since the review of 2026-11-12).
- **A controlling number with a threshold.** The runs sit at V_n / v_te = 1.0–1.5, where thermal electrons
  keep up with the front. A statement about injection needs the regime it is meant for, fixed before the run
  as targets of V_n / v_te, β, (m_i/m_e)/M_A, R_m/ρ_d and thermal ν τ, with the energy thresholds in units of
  ½ m_e V_n² and T_e instead of a fixed γ.
- **A baseline that is valid at ε ~ 1.** Where one reflection reaches the threshold, the static mirror is not
  defined (Sect. 6). Without a baseline that includes the motion of the foot, neither a residual nor a deficit
  can be read there.

**Steps, in the order in which they decide whether a run is worth its cost.**

1. **Read the two predecessors** (`kraussvarban1991electron`, `giacalone2017acceleration`). No computing.
   Without them no hypothesis can pass C5.
2. **Narrow the question** (the owner), with the controlling number and its threshold.
3. **One analysis of the data on disk: a count model of thermal escape (S1) per band,** from the measured
   downstream distribution, compared with the return shares by the registered factor-of-two rule. With the
   null share of μ (done, Sect. 3) and the thermal ν τ (done), it decides whether "returns beyond S0" is the
   numerical floor plus thermal escape. If it is, the observation is closed and is not carried forward.
4. **A shock requirement in Mach number and duration.** The first version asked for "a shock that is valid
   while the piston drives" and R_m well above 20 d_i. Both are replaced.
   - The coupling length of the debris is its gyro-radius ρ_d = M_A d_i, not d_i. R_m/ρ_d is 4.6 (L), 1.9 (M)
     and 1.0 (H); H, with one gyro-radius inside R_m, forms no shock.
   - The first valid front admits three readings (Sect. 4): a radius of 18 d_i, half of R_B, or about
     3 Ω_ci⁻¹. Two runs cannot choose. In ρ_d and in R_m the two runs disagree.
   - What a design should register instead: a fast Mach number of at least X for at least N Ω_ci⁻¹ inside the
     valid domain, with X above the first critical Mach number at the β and θ_Bn of the run. That number was
     not computed here. In these runs a fast Mach number of 2 or more is held in 72 of the 241 valid bins of L
     and 36 of the 180 of M, for at most 2.6 (L) and 1.3 (M) Ω_ci⁻¹ at one polar angle (counted from the
     recorded front arrays).
5. **A numerical floor for thermal electrons.** Thermal ν τ < 0.1 needs ν lower by a factor 1.2–2.8 (L) and
   2.5–6.4 (M) for the crossing of the measured ramp, 2.7–6.2 and 3.8–9.8 for the census layer, and 12–33 for
   the measured dwell of the returns. If ν ∝ 1/ppc, the last is 25–65 particles per cell. That scaling is not
   measured; the pair N6 (384³) would measure it.
6. **Denser tracks.** 30 tracked returns above γ = 1.25 need a stride of about 220 instead of 947 (about
   240 GB of tracks per run; the NVMe has 2.4 TB free).
7. **A higher mass ratio.** (m_i/m_e)/M_A is 4.6 and 2.7 by design at m_i/m_e = 16 (measured 7.6 and 6.3,
   because the realised M_A is 2.1 and 2.6); at m_i/m_e = 64 it is 18 and 11.

**The physical conflict, with sizes.** A stronger shock needs a larger M_A, and so a larger ρ_d = M_A d_i; a
higher mass ratio enlarges d_i. Both enlarge the box at fixed dx = 0.5 d_e. The table holds R_m = 5 ρ_d (the
value of L, which did not give the designed shock, so the sizes are lower bounds) and a box half-width of
2–3 R_m (L has 3.0).

| M_A, m_i/m_e | ρ_d | R_m | Half-width | Cells per side | Cells relative to 768³ |
|---|---|---|---|---|---|
| 3.5, 16 | 14 d_e | 70 d_e | 140–210 d_e | 560–840 | 0.4–1.3 |
| 6, 16 | 24 | 120 | 240–360 | 960–1440 | 2.0–6.6 |
| 3.5, 64 | 28 | 140 | 280–420 | 1120–1680 | 3.1–10.5 |
| 6, 64 | 48 | 240 | 480–720 | 1920–2880 | 16–53 |
| 6, 100 | 60 | 300 | 600–900 | 2400–3600 | 31–103 |

- The 768³ box at 2 particles per cell fills 67 % of the 8 GPUs. Step 5 multiplies the particles by 4 at
  8 per cell and by 12–33 at the numbers it gives.
- In 3-D on this machine only the first row fits, and it is the run already made. A stronger shock, a quieter
  plasma and a higher mass ratio cannot be had together in 3-D here. This is a limit of computing size, and it
  is named as such; it comes after the three missing items above, not instead of them.

**The 2.5-D alternative.** The particles of one 768³ run at 2 per cell fill 8192² cells at 13.5 per cell. In
2.5-D the baseline, the numerical floor and the Mach-and-duration requirement can be tested at a higher mass
ratio. Its limit: with an ignorable coordinate a particle stays tied to its field line (Jokipii, Kóta &
Giacalone 1993, Geophys. Res. Lett. 20, 1759; checked by the reviewer, not read by us), so 2.5-D cannot test
K1, the path across field lines. The direct predecessors are 2.5-D (`savoini2001twodimensional`,
`savoini2010origin`), so a 2.5-D run must also say what it adds to them. Whether the engine can run a 2.5-D
box, or a plane shock into a curved field (`kong2016electron`), was not checked.

**The test-particle baseline is not planned as it stands.** The first version proposed test electrons in the
stored meridional slices as the next step. As written it is unsound (R2, finding 19): the slices are 18.5
apart, about one gyro-period in L, and the front moves 2.8 d_e (L) and 2.0 d_e (M) between them, against a
ramp of 3.6 and 5.4 d_e, so the fields would have to be interpolated in the frame of the front; it needs a
stated smoothing kernel, two null tests and a decision rule fixed in advance. It comes after steps 1–3, if at
all.

## 9. Deviations from the protocol and choices made after seeing data

No registered criterion was changed. The first version said in addition that no deviation touched a
registered item; that was inaccurate. All of the following are in `notebook.md`.

**Registered items first applied wrongly or late:**

- **N3.** The first version evaluated ν τ for the electrons of γ 1.25–1.5 with the τ_dwell of B1′ (66, 113) and
  with the measured dwell, and wrote that it "would pass" at 66 and 113. The registered test is for thermal
  electrons over the crossing time of the layer. Applied as registered (Sect. 3), it fails in every reading.
  The conclusion "N3 fails" is unchanged; the statement that it passes at the assumed dwell is withdrawn.
- **The replacement rule of τ_dwell** (protocol Sect. 3.2) was not applied in the first version. It is applied
  in Sect. 5, with both results.
- **The reports of protocol Sect. 4** (surface differences, sensitivity to the debris threshold, usable bins),
  ε per band (Sect. 3.1) and the downstream branch of C1 were missing in the first version. They are in
  Sects. 4, 5 and 6.
- The ramp width that the registered N3 needs was measured only after R2, from the stored slices, as an upper
  bound.

**Choices of the analysis:**

- Census: an electron within d_i of the front is in the layer class even when its field line has not reached
  the front (conservative).
- N3c: three measures of the fluctuations; a fourth ("below 4 d_i") was dropped after its numbers were seen,
  because it measured the ramp.
- Potential: the mean over 5 outputs, and the mean over all outputs of a bin where the windowed value is not
  measured.
- The tracks of M were analysed with the census made before M-null existed. The final census differs in 9
  records, none of a tracked electron.

**Not registered; made after seeing data; support no claim:**

- the bin γ 1.25–1.5 and the share of returns at all energies (Sects. 5, 6);
- the retention of μ in the null runs, the ramp width, the regime table, the comparison with the MHD jump,
  the exponent of R(t), the id histories and the count of the bins before the valid domain (the checks made
  for R2; commands `r2`, `r2null`, `census --tau` of `analyze.py`);
- ρ_d = M_A d_i as the coupling length and the three readings of the first valid front (Sect. 8);
- the withdrawn sentence "the absence … is what such weak shocks give" was itself a post hoc reading.

## 10. Review R2

The research reviewer read the first version (commit `69459cd`) against the protocol, the notebook and the
records. The report is in `reviews/2026-11-13-results-r2.md`, unchanged. It is advice, not a decision of the
owner. The revised text was not reviewed again.

**Outcome.** "Not obtained from these runs" does not need rethinking, and nothing in the review moves the
answer towards "candidate". The review found errors of fact in the first version (the census of L, the bins
before the valid domain), a registered test applied to the wrong population (N3), one sentence to withdraw
and a design section to rewrite. Each was checked against the records before it was accepted; one point of
the first version was confirmed instead (the 5 upstream electrons of M at t = 111 do have feet in bins without
a detected front).

| Required change | Where |
|---|---|
| 1. Correct the census of L | Sects. 1, 5; notebook entry on R2 |
| 2. Redo N3 for thermal electrons with the registered τ; note on B2 | Sects. 3, 6, 7, 8 (step 5), 9 |
| 3. Withdraw "what such weak shocks give"; ε per band; the beam lies where the baseline is not defined | Sects. 1, 6 |
| 4. Rule of τ_dwell; reports of protocol Sect. 4; usable bins; downstream branch of C1 | Sects. 4, 5 |
| 5. Bins before the valid domain, reproducible; the "2 d_i" wording | Sect. 4 |
| 6. Early electrons: wording and normalisation; the 1.50 cut; regime table; "same shock" | Sects. 4, 5 |
| 7. Return shares marked unregistered; no loss of μ assignable to physics | Sects. 1, 6 |
| 8. Rewrite of Sect. 8 | Sect. 8 |
| 9. README: "Next" and the Prior work paragraph | `README.md` |

**Answers to the reviewer's questions.**

1. *Share that keeps μ over 250 in the null tracks:* 11 % in both null runs (Sect. 3), equal to the share of
   the returns.
2. *Ramp width per band and thermal ν τ:* medians 2.7–4.9 d_e (L) and 4.7–5.6 d_e (M) per band, upper bounds;
   ν τ = 0.10–0.37 (L) and 0.24–0.66 (M) (Sect. 3).
3. *The K1 electron of L:* one electron in seven dumps; not in the valid domain, inside the debris edge. Its
   field line was not traced.
4. *The decay of the energetic population of L:* the electrons made by the piston are left behind the front
   and fall below the cut; they are not replaced (Sect. 5). The cause of the energy loss was not analysed.
5. *Sensitivity of B2 to the errors of V_n and θ_Bn:* not computed (Sect. 6).
6. *Critical Mach number; b against the jump conditions:* the critical Mach number was not computed. b at
   θ_Bn ≥ 80° is at or slightly above the perpendicular MHD jump (Sect. 4).
7. *Selection of the 62 and 54 bins:* not reproducible. The records give 74 and 73 with the selection stated
   in Sect. 4.

**Not done, and why.**

- The count model of thermal escape (S1): the next analysis, Sect. 8 step 3. It needs a model fixed before it
  is compared, and was not squeezed into this revision.
- The exit θ_Bn of the returns: not in the records; needs a new pass over the tracks.
- The critical Mach number and the sensitivity of B2: not computed.
- L-s2 and N6: the protocol runs them only on a residual or a K1 signal; there is none.
- References that the reviewer gave from memory are not cited. Works that the reviewer verified and we have
  not read are named as such where they are used.
