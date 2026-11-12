# Exploratory runs: results and verdict

Topic `spherical-shock-injection`. Written 2026-11-13, after the runs L-null, L, M, M-null and H and their
analyses. Judged by `exploratory-protocol.md`, which was fixed before the runs; no criterion was changed.
Detail and failed attempts: `notebook.md` (entries of 2026-11-13). Records: `runs/` (`index.json`).
Units: c = 1, m_e = 1, lengths in d_e, times in 1/ω_pe; d_i = 4 d_e.

## 1. Verdict

**Question (protocol Sect. 1).** Do runs L, M, H, analysed, give a research hypothesis about the electron
acceleration mechanism that is novel and certain enough to build a study on?

**Answer: not obtained from these runs.**

- **No signal at the registered threshold.** No band of any run holds 100 ambient electrons with γ ≥ 1.5
  upstream of a front. The largest counts upstream are 3 (L, one dump) and 5 (M, one dump, on rays without a
  detected front). The null runs hold none. There is no residual, so there is nothing to build a hypothesis on.
- **The runs did not sample the designed shocks.** The design had fast Mach numbers 2.95 (L) and 4.05 (M).
  Where the registered definition finds a shock, the fast Mach number is 1.5–2.5 (median 1.8) and the
  compression is b ≈ 2 in both runs. The debris stalls at 1.3–1.5 R_m, and the front satisfies the
  definition only from about that time on: what is measured is a free, decaying blast wave. H forms no shock
  at all in 5.3 Ω_ci⁻¹. The absence of energetic electrons is what such weak shocks give; it says nothing
  about the designed Mach numbers.
- **The known baseline does not describe the runs.** The adiabatic mirror in a static potential (S0) predicts
  a beam of 10²–10³ electrons at γ ≥ 1.5 that is absent, and returns 0–6 % of the encounters where the tracked
  orbits return 1–15 %. This is a deficit against the baseline plus an excess of low-energy returns. It is not
  a residual of energetic electrons, it has known explanations that these runs cannot separate (Sect. 6), and
  one of them is numerical. It is recorded as an open observation, not as a candidate.
- **Why not "candidate".** The protocol asks for "a residual with a positive signature". None exists: K1 has 1
  electron in L and 2 in M against the 100 required; K2 and K3 have none.
- **What the answer does not say.** It does not say that a curved shock has no mechanism of its own. Two
  particles per cell, m_i/m_e = 16 and a shock of fast Mach 1.8 cannot test that.

The conditions that a run would need are in Sect. 8.

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

- Design: L has M_A = 3.5, R_m = 64 d_e = 16 d_i, Ω_ci⁻¹ = 56. M has M_A = 6.0, R_m = 44.8 d_e = 11.2 d_i,
  Ω_ci⁻¹ = 96. H has M_A = 8.7, R_m = 34.8 d_e = 8.7 d_i, Ω_ci⁻¹ = 139. R_B ≈ 148 d_e in all three.
- Outputs are on the NVMe, `/mnt/nvme/pypic/output/research/spherical-shock-injection/<name>`, linked from
  `output/research/spherical-shock-injection/<name>`; not in git.
- Not run: the second seed L-s2 and the 384³ pair N6. The protocol (Sect. 9) runs them only if L has a
  residual or a K1 signal.

## 3. Numerical floor (null runs)

| Measure | L-null | M-null |
|---|---|---|
| N1: electrons with γ ≥ 1.5, any dump | 0 (largest γ 1.374) | 0 |
| N2: T_e, start → end | 0.0100 → 0.0110 | 0.0100 → 0.0109 |
| N3: ν of the bulk (v 0.1–0.2), lags 4 / 16 / 64 | 1.2 × 10⁻² / 7.8 × 10⁻³ / 5.0 × 10⁻³ | 1.3 × 10⁻² / 7.9 × 10⁻³ / 5.0 × 10⁻³ |
| N3: ν for v > 0.3 | 2.5 × 10⁻³ / 1.8 × 10⁻³ / 1.3 × 10⁻³ | 3.3 × 10⁻³ / 2.0 × 10⁻³ / 1.4 × 10⁻³ |
| N5: 99.9th percentile of \|Δρ\|, start → end | 1.26 → 3.71 d_e | 2.07 → 5.96 d_e |
| Fronts detected | 0 of 738 bins | 0 of 738 bins |

- With a floor of 0 the count criterion of a signal is the fixed one: 100 electrons or more.
- The energy budget grows by 10–13 % in every run, with or without debris: numerical heating of the electrons
  at 2 particles per cell.
- N4 (closure of the work on tracked orbits) passes in L (1289 of 1341 orbits within 10 %) and in M (851 of
  853).

## 4. The front and its valid domain

| | L | M | H |
|---|---|---|---|
| Bins (time × θ) with a front | 443 of 738 | 326 of 738 | 0 of 738 |
| Valid bins (protocol Sect. 4) | 241 | 180 | 0 |
| Valid domain | θ = 30–90°, t = 185–666 | θ = 45–90°, t = 259–703 | none |
| Fast Mach, valid domain (min / median / max) | 1.50 / 1.82 / 2.51 | 1.50 / 1.75 / 2.38 | – |
| Fast Mach, design (piston) | 2.95 | 4.05 | 4.64 |
| V_n, valid domain | 0.116–0.204 | 0.089–0.147 | – |
| Compression b, valid domain (median) | 1.58–2.78 (2.0) | 1.62–2.50 (1.94) | – |
| Debris edge stalls at | 83 d_e = 1.3 R_m | 62–66 d_e = 1.4–1.5 R_m | ≈ 50 d_e |
| Superluminal feet | θ_Bn ≥ 77–83° | θ_Bn ≥ 78–85° | – |

- **No shock towards the pole.** No valid front exists below θ = 30° (L) and 45° (M). The quasi-parallel part
  of the sphere is not sampled.
- **The front becomes valid as the piston stalls.** The first valid outputs are at R = 73–86 d_e (L) and
  71–80 d_e (M). Before that, 62 (L) and 54 (M) bins show a front with fast Mach 1.9–2.8 (L) and 2.5–3.3 (M),
  and every one fails only criterion b: the front is not yet 2 d_i ahead of the debris. A stronger phase
  exists, but in it the compression has not separated from the piston.
- **The front decays.** At the equator the fast Mach number falls from 2.2 (t = 222) to 1.6 (t = 592) in L and
  from 2.2 (t = 296) to 1.6 (t = 592) in M (table in `notebook.md`).
- **H:** the ambient ion density reaches 1.5 n₀ only in 3 equatorial bins at the end of the run.
- The Mach scan of the design (L against M) therefore did not happen: both valid domains hold the same weak
  shock. The runs differ in the magnetisation of the upstream, not in the strength of the shock.

## 5. Census at the registered threshold (γ ≥ 1.5)

| | L | M | H |
|---|---|---|---|
| Ambient electrons, all classes, largest dump | 486 (t = 74) | 33 (t = 92.5) | 1 (t = 111) |
| …in the dumps of the valid domain | ≤ 40 (t = 222), then ≤ 4 | ≤ 1 | – |
| Upstream of a front, largest dump | 3 (t = 222) | 5 (t = 111, no detected front) | – |
| Bands with ≥ 100 upstream electrons | 0 | 0 | – |
| K1 electrons (needed: 100) | 1 | 2 | – |
| Tracked orbits with γ_max ≥ 1.5 | 11, all starting behind the front | 0 | – |
| K2, K3 on tracks | 0 | 0 | – |

- The early energetic ambient electrons are made at the start of the driver, before any valid front, and
  their number falls with the magnetisation: 486 (L), 14 (M), 0 (H) at t = 74.
- A missing K3 is not a result (protocol Sect. 7): the current filter damps the Buneman scale.
- Tracked returns with γ_out ≥ 1.25: 7 (L), 3 (M). Both are below the 30 the protocol needs, so neither run
  can test C1 or C2 on orbits.

## 6. The baseline against the runs

The baseline B2 (protocol Sect. 3.4) is the adiabatic mirror with the measured front, between an upper bound
(φ = 0) and a lower bound (measured φ).

- **Potential.** φ = 0.015–0.048 (1.5–4.7 T_e) in L and 0.009–0.023 (1–2.3 T_e) in M, for foot angles below
  70°; not measured above.
- **γ ≥ 1.5.** The model predicts on the detected front 740–3890 electrons in L (t = 185–407) and 107–1046 in M
  (t = 130–407). Measured upstream: 2, 0, 0 in L and 0 at every dump in M. No cell holds 100 measured
  electrons, so the rule of Sect. 3.4 has nothing to judge.
- **γ 1.25–1.5 (extra bin, not registered, exploratory only).** In the cells with ≥ 100 electrons, the measured
  count is outside [lower/2, 2 × upper] in 50 of 86 (L) and 9 of 13 (M). The measured electrons sit at low
  θ_Bn, where a mirror gives almost nothing; at 50–80°, where it should give the most, they are missing by a
  factor of 200 (L) to 500 (M) or more. By the rule of Sect. 3.4 no residual is discussed.
- **Share of the encounters that end in a return, all energies** (tracked orbits, detected front):

  | | L (θ_Bn 10–90°) | M (θ_Bn 30–90°) |
  |---|---|---|
  | Model, φ = 0 | 0–13 % | 0–13.5 % |
  | Model, measured φ | 0.07–0.5 % | 0–6 % |
  | Tracks | 1–9 % | 7–15 % |
  | Tracks, only returns from behind the density front | 0.6–6 % | 4–7 % |
  | Returns that keep μ within 20 % | 12 % | 12 % |

  At θ_Bn = 70–90° the mirror returns almost nothing with or without a potential; the tracks return 1–2 % (L)
  and 7–8 % (M).
- **Reading.** Most returns are not adiabatic mirror reflections in a static potential, in L and in M alike.
- **Explanations that these runs cannot separate:**
  1. *Numerical scattering.* ν τ at the measured dwell in the layer is 0.23–0.39 (L) and 0.11–0.19 (M), above
     the limit 0.1 of N3. By the protocol, no loss of μ may then be assigned to a physical mechanism.
  2. *Marginal adiabaticity by design.* Ω_ce d_i / V_n = (m_i/m_e) v_A / V_n ≈ 7.6 (L) and 6.2 (M): an electron
     crosses a ramp one d_i wide in about one gyro-period. The ramp width was not measured.
  3. *Time dependence.* The front decelerates and φ changes; the model is static per foot.
  4. *Leakage (S1), known.* The foot of a field line slows from V_HT = ∞ at the tangent point towards V_n,
     and the heated electrons behind the front (thermal speed 0.12) can overtake it. This fits the excess at
     low θ_Bn and the small gain (median 0.02–0.04). Published for a curved front (`savoini2010origin`).
  5. *Kinematic returns.* Half of the "back upstream" orbits never reach the density front.
  6. *Reduced gain at a curved front, known* (`vandas1994acceleration`: about 1/4 of the plane-shock gain).

## 7. Criteria C1–C6

| Criterion | Energetic residual (γ ≥ 1.5) | K1 (Δρ) | Returns beyond S0 (not registered) |
|---|---|---|---|
| C1 real: ≥ 100, ≥ 10 × floor, second seed | fails: 0–5 electrons | fails: 1 (L), 2 (M) | not a residual; one seed |
| C2 mechanism on orbits | not testable: 7 and 3 tracked returns | – | μ is lost, cause undecided |
| C3 not numerical, not the driver | – | – | fails: ν τ > 0.1; N6 not run |
| C4 testable in dimensionless numbers | – | – | not stated |
| C5 new to the field | – | – | no: items 4 and 6 are published |
| C6 review R2 | Sect. 10 | Sect. 10 | Sect. 10 |

K2, K3 and K4 have no electron and no signature. Independently of the above, C5 could not be met: the two
predecessors `kraussvarban1991electron` and `giacalone2017acceleration` are still unread (protocol Sect. 10).

## 8. What a run would need

In the order in which they decide whether a further run is worth its cost.

1. **Read the two predecessors** (`kraussvarban1991electron`, `giacalone2017acceleration`). No computing.
   Without them no hypothesis can pass C5.
2. **A baseline that can fail for the right reason.** Test electrons in the measured, time-dependent fields
   of the meridional slices of L and M (E and B at full resolution, taken as axisymmetric): guiding-centre
   orbits (S0′ of the protocol, Sect. 3.2) and full orbits. Test particles in smoothed fields have no
   numerical collisions, so the full orbits separate explanation 1 from 2–5, and the guiding-centre orbits
   separate 2 from 3–5, with the data already on disk. Limit: the slices are 74 steps apart (18.5, against a
   gyro-period of 22–38), so waves faster than that are not in them.
3. **A shock that is valid while the piston drives.** In both runs the front is 2 d_i ahead of the debris only
   at R ≈ 71–86 d_e ≈ 18–21 d_i, and the debris stalls at 1.3–1.5 R_m. A driven, separated shock therefore
   needs R_m well above 20 d_i (two runs; not a law). At m_i/m_e = 16 that is R_m > 80 d_e and R_B > 185–260
   d_e for M_A = 3.5–6, against a box half-width of 192 d_e.
4. **More particles per cell.** ν τ < 0.1 at the measured dwell needs ν lower by 2.3–3.9 (L) and 1.1–1.9 (M):
   about 8 per cell if ν ∝ 1/ppc. The pair N6 would measure that scaling at 384³.
5. **Denser tracks.** 30 tracked returns above γ = 1.25 need a stride of about 220 instead of 947 (about
   240 GB of tracks per run; the NVMe has 2.4 TB free).
6. **A higher mass ratio** for adiabatic electrons at the ramp.

**Items 3, 4 and 6 pull against each other at fixed memory.** The 768³ box at 2 particles per cell fills 67 %
of the 8 GPUs. A larger R_m, more particles per cell and a higher mass ratio each need more cells or more
particles. In 3-D on this machine they cannot all be met. The difficulty is computing size, not a missing idea,
and it should be named as such (research stance). Items 1 and 2 cost no GPU time and come first: if S0′
reproduces the returns in smooth fields, the observation of Sect. 6 is known kinematics and no further run
is justified by it.

## 9. Deviations from the protocol and choices made after seeing data

None changes a registered criterion. All are in `notebook.md`.

- Census: an electron within d_i of the front is in the layer class even when its field line has not reached
  the front (conservative).
- N3c: three measures of the fluctuations; a fourth ("below 4 d_i") was dropped after its numbers were seen,
  because it measured the ramp.
- Potential: the mean over 5 outputs, and the mean over all outputs of a bin where the windowed value is not
  measured.
- The bin γ 1.25–1.5 and the share of returns at all energies are not registered. They are reported as
  exploratory and support no claim.
- The tracks of M were analysed with the census made before M-null existed. The final census differs in 9
  records, none of a tracked electron.
- N3 is evaluated at the measured dwell (250–418) as well as at the protocol's assumed 66 (L) and 113 (M). At
  the assumed dwell ν τ is 0.06 (L) and 0.05 (M) and would pass.

## 10. Review R2

Open: to be added after the review of this document by the research reviewer.
