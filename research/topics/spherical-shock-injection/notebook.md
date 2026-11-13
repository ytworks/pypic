# Notebook — spherical-shock-injection

Dated lab notes, newest first. Write down what was tried, what was seen, and what was decided, failed attempts included.
Never rewrite a past entry. To correct one, add a new entry that refers to it.

## 2026-11-13 — M6 routes: rules fixed on P5 (pilot), before any campaign run is classified

- `routes.py` (`analyze.py routes`), rules in its docstring. Revised three times on P5 only (P5 is the base run,
  not a campaign verdict run): select by the drift-frame level, not the lab γ_max (gyration phase in the moving
  downstream); R2 and R3b conditioned on where the gain is, E_φ only a diagnostic; step gains from the change of
  the smoothed γ, because the work recorded every 2 steps closes on only part of the energetic orbits (block work
  over Δγ, median 0.66); "unplaced" for gains where the sites are not defined (no magnetic wall, near the
  poles). Frozen after these.
- **P5, top 1 % of the electrons born outside the shell** (783; drift-frame level ≥ 1.41): R2 (heated, then the
  gain on the wall: foot, steep wall, cavity) 59 %; R3b (long stay in the layer) 15 %; unplaced 21 %; R3a 0.8 %;
  R5 0.1 %; R4 0.8 %; R1 and R3c 0; other 3.6 %. In P5 the electrons that come through the shock reach their
  energy mostly in the downstream wall after the shock heated them, not by reflections. 21 % gain where the sites
  do not reach (high latitude): the site definition has to be extended to the polar caps (global plan, open).

## 2026-11-13 — Correction: "Fermi" and "DSA" were used against their definitions

- **Owner:** what is DSA; is it really Fermi; does the definition agree with the prior work? It did not.
- **DSA** (first-order Fermi at a shock): scattering on both sides and repeated crossings from upstream to
  downstream and back; Δp/p ≈ (4/3)(U₁ − U₂)/v per cycle; it needs a diffusion length κ cos²θ / U_s longer than
  the shock thickness (`amano2022theory`; the founding papers have no note yet). **SSDA**: confinement in the
  transition layer by pitch-angle scattering with gradient drift along the motional field (`katou2019theory`,
  `matsumoto2017electron`).
- **Withdrawn:** (1) the label "Fermi at the polar cap (DSA-type)" for repeated reflections (`global-plan.md`
  R3, first version): the repeated returns seen in P5 are repeated SDA reflections, with no downstream
  scattering; R3 is now split into R3a (repeated SDA), R3b (SSDA), R3c (DSA). (2) The criterion "constant
  relative gain = Fermi-like" (entries "several heating sites" and "the P5 shock reflects electrons upstream"):
  for v < c the DSA gain per cycle falls as 1/v, so a relative gain that falls with γ does not exclude it.

## 2026-11-13 — Global plan: the whole sphere first (owner: the forest is missing)

- **Owner:** the parts are important, but there is no hypothesis for the acceleration over the whole spherical
  shock; plan the whole. Written: `global-plan.md` with `figures/global-map.svg`.
- **Question:** where on a spherical shock in a uniform field are electrons injected and accelerated, by which
  processes, and is the whole the area-sum of plane shocks at the local θ_Bn (the assumption of remnant models)
  or does the sphere couple its latitudes? Hypotheses G0 (local), G1 (coupling along field lines), G2 (closed
  azimuthal drift), G3 (focusing); D (the driver's wall) kept apart.
- **Campaign A** is a building block (the local reflection law) and, with the histories now saved, the first data
  for M1–M4 over the sphere. The baseline G0 needs plane oblique electron–ion shocks: an example on main (PR).
- **Changed on the way:** the campaign pass now saves compact histories of electrons with γ_max ≥ 1.25 and a
  2 % sample; the post-processing was stopped before A5-T20 and A5-T28 were pruned and restarted with it.

## 2026-11-13 — Item 2 begins: the P5 shock reflects electrons upstream at high latitude (SDA level)

- **Literature against P5** (notes; shock M_A ≈ V_sh / v_A ≈ 3.8, β_e ≈ 0.25, m_i/m_e = 16): unlimited SSA
  trapping M_A ≥ (1+α)(m_i/m_e)^(2/3) ≈ 7.6 not met (`matsumoto2012electron`, `bohdan2019kinetic`); Buneman onset
  ≈ 1.2 met; SDA-like oblique regime M_A ≲ (m_i/m_e)^(1/2) = 4 at its edge and escape for θ_Bn ≲ 84°
  (`riquelme2011electron`); potential effect for θ_Bn ≲ 75° (`amano2007electron`); efficient SDA 7 ≤ M_A ≤ 33,
  β_e ≈ 1 not met (`lezhnin2021kinetic`).
- **Missed before:** the upstream γ ≥ 1.5 electrons at polar angles 0–50° (counts) are in the upstream rest
  frame, where no frame correction applies. The upstream "thermal" verdict of the local-fit test was wrong for
  beams: a beam inside a cell raises the fitted temperature of that cell.
- **Reflections** (script `P5/analysis/figures/reflect.py`; electrons born outside the shell that enter the
  transition layer from upstream and return upstream): 57,951 events by 46,801 electrons. Returning with
  γ ≥ 1.25: 641 events (419 electrons), ≥ 1.5: 16 (10), max 1.65, none ≥ 2. Share returning with γ ≥ 1.25 by
  polar angle: 3.9 % at 10–20°, 2.7 % at 20–30°, 2.1 % at 30–40°, 1.0 % at 40–50°, about 0 above 60°; local θ_Bn
  of those 16–47°, V_HT 0.27–0.37 c. γ_out against the head-on bound γ_HT²(1 + V_HT²) from rest: 0.85–1.26
  (median 1.05). μ_out / μ_in median 1.18, wide.
- **Repeated reflections:** 4792 electrons reflect twice or more (polar ≈ 25°). Between reflections the energy
  does not change. The relative gain per reflection falls with γ: +2.0 % at γ 1.0–1.05, +0.5 % at 1.05–1.15,
  −1.3 % at 1.15–1.3 (not Fermi-like); repeated reflections settle near the SDA level ≈ 1.24. Many repeats come
  9 time units apart and may be one encounter flickering across the boundary r_f + d_i.
- **Sphere-specific idea tested:** on a field line of fixed ρ the expanding sphere is met at θ = asin(ρ/R(t)),
  so later encounters lie closer to the pole. 89 % of successive reflections move pole-ward, by 0.6° (median):
  seen, small. γ_out of electrons entering with γ < 1.1 stays at 1.04–1.06 for V_HT from 0.25 to 0.7 c, against a
  bound rising from 1.17 to 1.94: the gain does not follow V_HT; the idea does not explain the falling gain.
- **Reading.** The P5 shock does accelerate electrons by reflection at high latitude (θ_Bn 15–50°), to about the
  SDA level of its V_HT (γ ≈ 1.25, at most 1.65). At more perpendicular latitudes (larger V_HT, larger bound)
  hardly any electron is reflected. Open: why the gain stays flat in V_HT; what returns electrons for separate
  encounters.

## 2026-11-13 — Path (a): electrons born inside the initial debris shell

Script and output: `P5/analysis/figures/patha.py`, `patha.json`. Drift frame and smoothing as for path (b).

- **Birth.** 33,697 tracked electrons start inside the initial debris edge (78 d_e). Share reaching
  γ_max ≥ 1.5 (≥ 2) by birth radius: 10 % (0.4 %) at 40–50 d_e, 20 % (2.4 %) at 50–60, 15 % (1.5 %) at 60–70,
  6 % (0.3 %) at 70–80, about 0 below 30 (shell peak at 46 d_e, width 14). By polar angle: 0.35 % at 0–15°,
  12–14 % at 45–90°. Born outside the shell: 0.2 % (0.002 %).
- **γ_max ≥ 2 (336).** Main gain window (largest rise of the gyration-averaged γ' over 20 units) starts at
  t ≈ 50–150 (median 123; 21 of 336 before t = 50). Place: foot 44 %, cavity 22 %, steep wall 22 %, transition
  layer 0.2 %; 14 d_e inside the debris edge (median); |B| 0.9 B₀; polar angle 77° (median); drift speed 0.30 c.
  Work over the windows: E_φ +87, E_par +16, E_θ +15, E_r +7 (smoothed), all +90 (unsmoothed); Δγ (lab) +104,
  Δγ' (drift frame) +153. Smoothed p⊥²/B grows 2.3 times while γ' grows 1.25 times: not adiabatic. The
  gyration-averaged γ' reaches 1.92 (median): energy also in the plasma frame.
- **γ_max 1.5–2 (3420):** the same pattern, weaker (γ' level 1.42; foot 36 %, steep 34 %).
- **Control (inside-born, γ_max < 1.3; 3000):** born nearer the centre or the edge (10 % below 30 d_e), more time
  in the cavity, main windows spread to late times; γ' level 1.11.
- **Reading, not decided.** The gain is a drift along E_φ on the inner side of the magnetic wall (foot), not
  adiabatic, real in the plasma frame. What separates (a) from (b) is the side of the wall: (a) sits on the
  inner side, (b) on the outer side and is only heated. The electrons that reach the inner side are almost only
  ambient electrons placed inside the debris shell at t = 0: the mechanism may be physical, but its population
  depends on the initial overlap of ambient plasma and debris in P5.

## 2026-11-13 — Path (b) followed along its history: shock heating, and in the plasma frame no more

Owner: find the mechanism of path (b) (heated at the shock, then more) from the place on the sphere and the
fields along the particle history. Scripts in `P5/analysis/figures/`: `pathb.py`, `pathb2.py`,
`pathb_orbits.py` (figure `p5_pathb_orbits.png`), `frame.py`. Sites of `sites.md`.

- **Selection.** All tracked electrons whose first entry from upstream into the transition layer leads into
  the steep wall and which later reach γ_max ≥ 1.5 (lab frame): 751. Control: 5000 of the 161,826 that did the
  same and stayed below 1.3.
- **Place on the sphere.** Polar angle at entry: median 50° (control 65°); 10–30° is 4–5 times more frequent
  than in the control, 40–90° 0.6–0.7 times. Local θ_Bn at entry: median 50° (control 65°).
- **History.** One entry each (no repeated encounters; 13 % of the control enter twice). γ before entry 1.044
  (control 1.011), gained 4–30 d_e ahead of the front (+0.022, median). First passage of the layer: Δγ = +0.19
  (control +0.05). The largest 20-unit gain lies 83 % in the transition layer, 17 % on the steep wall; the work
  is done by E perpendicular to B (E_φ the largest part; E_par under 3 %); p⊥²/B grows 3.4 times while |B|
  grows 1.4 times (not adiabatic over that window); pitch nearly perpendicular. After entry the electrons ride
  5–10 d_e behind the front for about 1.6 / Ω_ci; γ swings over each gyration from about 1.15 to 1.8.
- **Frame.** In the local E × B frame (fields smoothed over 12.5 along the orbit; drift 0.22 c at γ_max), path
  (b) has γ' 1.34 at most (median; 3 % reach 1.5) and 1.27 on average around γ_max (p90 1.35). The lab-frame
  γ ≥ 1.5 of path (b) is mostly gyration in the outflowing downstream seen from the lab: heating, not
  acceleration beyond it. The γ_max ≥ 2 electrons born inside the shell keep γ' 1.74 on average (p90 2.11).
  The maxima γ'_max are not reliable where B is weak (the random sample shows γ'_max above its lab γ_max).
- **Reading.** In P5 the shock heats the electrons that cross it once, more at high latitude (polar 10–30°,
  where they also arrive slightly preheated), and does not accelerate them beyond the heated distribution in the
  plasma frame. The energetic electrons in the plasma frame are those born inside the debris shell (path a).
  Limits: drift frame from smoothed fields; tracks every 2 steps; one run.

## 2026-11-13 — Unselected electrons and crossings of the transition layer (sites of sites.md)

- **Owner:** "the transition layer gives a net loss" has the wrong premise: to lose, an electron must have
  gained before; the layer may just hold electrons before and during their energisation. Agreed: that statement
  came from electrons selected by γ_max ≥ 1.5 or 2, windowed up to their γ_max. **Withdrawn** as a property of
  the layer. Script: `P5/analysis/figures/crossing.py` (100,000 random tracked electrons, no selection).
- **(1) Net change of γ per site, whole run:** cavity +192 (1 % of the time), foot +834 (1 %), steep wall
  +1732 (11 %), transition layer −1204 (9 %), upstream +96 (78 %). A sum per site mixes populations (electrons
  coming in from upstream and electrons coming out from inside); its sign is not a property of the place.
- **(2) First entries from upstream into the layer (23,518):** into the steep wall 18,558, of which 93 % gain
  (Δγ p10/50/90 = 0.005 / 0.054 / 0.147, mean +0.067): shock heating. Back upstream 4891, Δγ ≈ 0.
- **(3) After the layer:** of those that gained more than 0.05 (9750), 0.7 % later reach γ ≥ 1.5; of those that
  gained nothing (1254), none. A heating-then-wall path exists and is rare.
- **(4) Origin:** γ_max ≥ 2 in the sample: 32, none ever entered the layer from upstream; all born at
  r = 51–69 d_e (inside the initial debris shell). γ_max 1.5–2: 574, of which 114 came from upstream.
- **Readings, not decided:** (a) the highest energies belong to ambient electrons placed inside the debris
  shell at t = 0 and pushed from the start, a population that may depend on the initial set-up; (b) shock
  heating followed by gain at the wall is a second path, reaching γ 1.5–2 in P5.
- `analyze.py sitegain` (run before this) groups by γ_max with the same window and has the same selection; its
  per-site sums are not used.

## 2026-11-13 — Correction: the gain at the wall is not shown to be "not the shock"

- **Owner:** the wall foot is a place the shock passed earlier; the earlier claims were one-sided.
- **Facts.** The front crossed r = 91–103 d_e at t = 72–108 (equator), so the foot at t = 144 (91–99 d_e) is
  plasma the front swept. r_Bmax moves with the front, 4–7 d_e behind it. The compressed field and the motional
  E_φ there come from the shocked downstream flow together with the field the piston pushes out.
- **Withdrawn:** the reading of the entry "does the gain zone follow the front or the debris edge" that the
  energisation is "a mechanism of the piston, not of the shock", and the same words in the replies to the owner.
  That the gain zone keeps its distance to the debris edge does not assign it to the piston: the downstream of a
  driven shock is bounded by the piston. The data do not separate the two.
- **Kept, as observations:** the gain is by drift along E_φ on the inner side of the wall; the front layer gives
  a net loss to these electrons; the gaining electrons were born at r = 52–71 d_e, inside the initial debris
  shell. What these mean for the mechanism is open.

## 2026-11-13 — P5: does the gain zone follow the front or the debris edge? The debris edge

- **Owner's questions:** how the debris edge and the front differ, and whether a gain near the front is not an
  acceleration mechanism. Definitions (`analyze.py`, protocol Sect. 4): the front is the outermost radius where
  the ambient ion density (averaged over d_i) reaches 1.5 n₀; the debris edge is the outermost radius where the
  debris ion density reaches 0.1 n₀. Two species, measured independently.
- **Test** (script and figure `P5/analysis/figures/host.py`, `p5_host.png`): net change of γ of the tracked
  electrons against r − r_front and against r − r_debris, in two windows while the two separate.
  - γ_max ≥ 2: against the front the peak moves from −17 to −27 d_e (gain-weighted mean −19 to −28) while the
    front leads the debris edge by 4 then 13 d_e; against the debris edge it stays at −11 then −9 (mean −17,
    −12).
  - γ_max 1.5–2: front −15 to −27 (mean −19 to −25); debris edge −13 to −9 (mean −15 to −12).
- **Reading.** The gain zone moves with the debris edge, 9–13 d_e inside it, and falls behind the front as
  the shock separates. It is an acceleration mechanism of the piston (the magnetic wall at the debris edge),
  not of the shock transition, whose net work on these electrons is negative. Limits: two windows, one run.

## 2026-11-13 — P5: several heating sites? By the fixed rule "one process"; the fixed sites split one structure

Scripts and figure: `P5/analysis/figures/sites.py`, `sitedef.py`, `p5_sites.png`.

- **Result by the rule fixed in the entry below.** γ_max ≥ 2 (350): net gain S1 cavity +395, S2 inner edge
  −76, S3 −1, S4 0. Criterion 1 not met; S1 gives 100 % of the positive net gain: **"one process"**.
  Criterion 2 not met (median γ_max 1.57–1.62 whatever the number of site changes). Criterion 3 not met (S1:
  mean Δγ/γ per visit 0.18 at γ 1–1.5, 0.076 at 1.5–3, ratio 0.41). γ_max 1.5–2 (5391): S1 +1774, S2 +801,
  S3 −68, S4 +18. Reference (20,000 with γ_max < 1.25): 84 % of the time upstream.
- **What the figure shows.** At the equator (t = 144) |B| is 0.5 B₀ in the cavity, rises over about 10 d_e
  to 4.5 B₀ at r_Bmax = 109 d_e, and the front is at 114 d_e; the debris edge (107 d_e) is just inside r_Bmax.
  The γ_max ≥ 2 electrons gain at r − r_Bmax = −27 to −9 d_e (|B| 0.5–0.9 B₀ at the electron: the foot of
  the wall, inside the dense debris shell) and lose at −8 to 0 (the steep wall, |B| 1–3 B₀). They spend their
  time at −25 to −3. The fixed boundary at r_Bmax − 10 cut this one wall into a gain part (S1) and a loss part
  (S2): the fixed sites were not two heating mechanisms.
- **Correction to the entry "the fields where the electrons gain energy":** its "net +520 at 0.5–1.5 B₀" was
  binned by |B| at the electron, not by place; the gain is at 0.5–0.9 B₀, 10–27 d_e inside r_Bmax, in the foot
  of the wall, not on its steep part.
- **Reading.** With a Larmor radius of about 10 d_e the electrons sample the foot and the steep part of the
  wall within one gyration; the net is a gain. Sites drawn from the structure (foot, steep wall, front layer,
  upstream) would be a definition made after seeing the data: any test with them is post hoc.

## 2026-11-13 — P5: several heating sites visited in turn? Criteria fixed before the analysis

- **Owner:** do not hold on to one process. If two heating mechanisms sit at different places and some
  electrons enter both selectively, those electrons are accelerated: the different heating sites act as the
  scatterers of a kind of second-order Fermi process.
- **Sites** (radius r of the electron against, at its polar angle and time, the radius r_Bmax of the
  B maximum and the front r_f from `front.npz`; d_i = 4 d_e):
  S1 cavity, r < r_Bmax − 10 d_e; S2 inner edge of the piston, r_Bmax − 10 ≤ r < r_Bmax; S3 outer side of the
  pile-up and the front layer, r_Bmax ≤ r ≤ r_f + d_i; S4 upstream, r > r_f + d_i. A visit is a run of
  consecutive samples (every 2 steps) at one site; its gain is the change of γ over it.
- **Data:** the tracked electrons with γ_max ≥ 1.5 (5741) and, as the reference, a random sample of 20,000
  tracked electrons with γ_max < 1.25; each one up to the time of its γ_max.
- **Fixed now. "Gains from several sites" (the owner's picture) holds if all three hold:**
  1. For γ_max ≥ 2, at least two sites each give more than 10 % of the summed positive net gain.
  2. At fixed total time at the gaining sites, γ_max rises with the number of site changes (median in the
     upper third of the number of changes above the median in the lower third, in each of three bins of
     time).
  3. At each gaining site, the mean relative gain per visit, Δγ / γ at entry, differs by less than a factor 2
     between γ at entry 1.0–1.5 and 1.5–3 (a gain proportional to the energy, as for Fermi).
- **"One process" holds** if one site gives 90 % or more of the summed net gain of the γ_max ≥ 2 electrons.
  Between the two: undecided, and the numbers are reported.
- **Second-order signature, reported in any case:** per site, the mean and the rms of the gain per visit; the
  rms against the number of visits.

## 2026-11-13 — P5: the fields where the electrons gain energy (owner: observe the acceleration itself)

Scripts and figures: `P5/analysis/figures/fields.py`, `work.py`, `orbitfield.py` (NVMe output, not in git).

- **Fields** (slices y and z = centre, t = 108–180). The debris interior is a magnetic cavity (|B| = 0–0.5 B₀)
  bounded by a pile-up of 3.5–4 B₀ between the debris edge and the front. A smooth, positive E_φ of 0.1–0.3
  (motional, V_r B_z) lies on the pile-up. Inside the cavity E fluctuates at the grid scale with ±0.3 (rms at
  the electrons 0.12–0.16, against 0.09–0.10 for the smooth part). The mean electron energy (0.1–0.2 m_e c²)
  is highest in the band just inside the debris edge, where the debris ions are dense.
- **Work on the tracked electrons up to their γ_max** (E and v averaged over 60 for the smooth part):
  - γ_max ≥ 2 (350): Δγ summed 430; work of all fields 295 (the rest not closed: tracks every 2 steps miss
    part of the fast fields). E_φ: +195, of which smooth +589 and fast −394; E_r +78 (fast); E_θ +22.
  - Gains by |B| at the electron: net +520 at 0.5–1.5 B₀ (73 % of the time), +56 below 0.5, −136 at
    1.5–3 and −146 above 3. Smoothed v_φ: median −0.09 c (10 %: −0.29 c), against E_φ = 0.04–0.12; the
    product −E_φ v_φ is positive 87 % of the time.
  - γ_max 1.5–2 (5391): the same pattern (smooth E_φ +3405, fast E_φ −1419; smooth E_r −811).
- **Orbits.** The γ ≥ 3 electrons move on large loops (Larmor radius about 10 d_e at γ ≈ 3 in B₀) across the
  inner edge of the pile-up: half a turn in the cavity, half in the pile-up. The loops drift in −φ along the
  ring. γ rises and falls within each loop. Over many loops the net gain comes from the smooth E_φ.
- **Reading (an observation, not a claim).** The electrons gain energy by drifting in −φ along the motional
  E_φ at the inner edge of the magnetic piston, where |B| grows outward. There the ∇B (and meandering) drift
  of an electron points in −φ. On the outer side of the pile-up (the shock ramp) the gradient is reversed and
  electrons lose energy in this frame: this fits the negative net work in the front layer. The grid-scale
  fluctuations take energy on balance; they are not the source.
- **Limits.** Work closure 69 % for γ_max ≥ 2. One run, mass ratio 16. The cavity and its edge are a property
  of the debris model (a dense shell expanding at 0.357 c). No literature checked yet for electron
  energisation at the edge of a diamagnetic cavity.

## 2026-11-13 — P5: were electrons accelerated at the shock and then moved away? Not in the tracks

- **Owner:** the snapshot test cannot tell an electron accelerated at the shock that then moved inward or
  escaped. Tested on the tracked electrons (956,674 ambient, every 2 steps; script
  `P5/analysis/figures/where.py`): where each one first crosses a γ threshold, where its energy changes (sum of
  dγ per region of the step), and where it ends (or whether it leaves the box). Regions: upstream (beyond the
  front + d_i), front layer (± d_i), shock downstream (debris edge to front − d_i), debris interior.
- **γ ≥ 2 (350 electrons):** 349 cross inside the debris edge (t ≈ 150, polar angle ≈ 75°), one in the shock
  downstream. Gains: +3221 inside the debris, +2.8 in the shock downstream, +0.6 in the front layer, 0 upstream.
  None leaves the box; 347 end inside the debris. Born at r = 52–71 d_e (inside the initial debris shell).
- **γ ≥ 3 (11):** all cross, gain and end inside the debris.
- **γ ≥ 1.5 (5741):** first crossing inside the debris 5037, shock downstream 577, front layer 124, upstream 3.
  Net change of γ by region: debris +1504, front layer −285, shock downstream −56, upstream +12. 13 leave the
  box, at γ ≈ 1.3.
- **Reading.** No tracked electron gains its energy at the front and then leaves it: the energetic electrons
  are energised in the debris interior, and the front layer takes energy from electrons on balance. Limits:
  tracks sample 1 in 947 of the initial ambient electrons; electrons that enter at the faces are not tracked.

## 2026-11-13 — P5: is there a nonthermal component at all? Not at the shock; yes inside the debris

- **Owner:** without a nonthermal component the question does not arise. So the first test is an excess over
  the local thermal distribution, not the counts above a fixed γ (a sum of regions of different temperature
  looks like a tail).
- **Method** (functions of `census.py`; script and output `P5/analysis/figures/nonthermal.py`, `.txt`): the
  meridional slab spectra (cells of 2 d_e), a Maxwell–Jüttner fit of the bulk of each cell, and the measured
  weight above γ − 1 = 0.56 and 1.0 against the fit (0.9, 1.0, 1.1 T). Lab-frame energies (bulk flow
  included). Regions: A, shock downstream (from 3 d_i behind the front to the debris edge); B, inside the debris
  edge; C, upstream.
- **A, at the shock:** thermal. Above 0.56: 0.06–0.69 against 0.06–2.3 expected (t = 252–360). Above 1.0: one
  or two macro-electrons (weight 0.0625 each). The region is thin (12–18 d_e) and exists in the slab only from
  t = 252: 160–1350 cells, so the test is weak.
- **B, inside the debris edge:** a tail above the local fit. Above 1.0: 32–126 against 1.9–3.0 expected
  (15–40 times, t = 108–288); above 0.56: 2–3 times. This is where the γ ≥ 3 electrons are.
- **C, upstream:** thermal within the fits (above 1.0: at most 0.8 against 0.16).
- **Reading.** In P5 no nonthermal electrons are seen at the shock; the only nonthermal component is in the
  piston region. Limits: one slab (y near 0), the cell-mixing of temperatures and the anisotropy are not
  separated, and B's tail is not yet shown to be physical rather than numerical.

## 2026-11-13 — P5 electrons, first look: maps first, highest energies first (owner's instruction)

- **Owner:** look at the distribution maps first for the whole picture, then from the highest energies down.
  Figures and their scripts (not in git; observation only): `P5/analysis/figures/` on the NVMe output.
- **Why the census was replaced by counts.** P5 dumps 1–7 million ambient electrons with γ ≥ 1.25 per output
  (S25: 9–74 thousand); the census, which traces each electron along its field line, did not finish in 30 min.
  `analyze.py counts` gives O8 as fixed (counts by class, band and energy; no result, origin or mechanism).
- **Whole box.** γ ≥ 1.5: 0.6 million at most (t = 180); γ ≥ 3: 1469 at most (t = 144, 2.6 / Ω_ci); largest
  γ 4.8. The tail of the spectrum reaches γ − 1 ≈ 3–4 at t = 144–216 and recedes after. L-null: no electron
  with γ ≥ 1.5 in any dump.
- **The highest energies are not at the shock.** γ ≥ 3 electrons sit at polar angles 60°–88°, 5–20 d_e inside
  the debris edge and 10–30 d_e behind the front, with pitch cosines near 0 (median |μ| 0.17–0.28). The six
  tracked electrons with the largest γ (4.15 at most; 11 tracked reach γ ≥ 3, 350 reach γ ≥ 2) are ambient
  electrons born inside the debris shell (r = 55–68 d_e). They stay inside the debris edge and never meet the
  front. Their γ rises and falls by 1–2 within each gyration (period 45–50, Larmor radius about 6 d_e); the
  work is done by E perpendicular to B; the envelope peaks near t = 150 and decays.
- **Upstream.** γ ≥ 1.5 ahead of the front: up to 3069 per dump (t = 360), all at foot θ_Bn and polar angles
  of 0°–50° (outside the valid domain), mostly γ 1.5–2; none at 50°–90°. At the last dump their distance to the
  z faces is 29 d_e (median; 10 % within 6.5 d_e), and 90 of 2719 were injected at a face.
- **Not decided.** Whether the piston-region energization is physical (debris ions streaming through ambient
  electrons) or numerical (dense shell, ω_pe dt up to 1.05) is open: L-null has no debris and cannot tell.

## 2026-11-13 — P5 run: a strong shock on a sphere was made (verdict fixed before the run)

- **Run.** `P5` from `e6f1890`: 1440 steps, 755 s, 53 GB; Gauss 1.3e-14, div B 1.8e-14, `fill` 0.71, energy
  budget +4.4 %. Record `runs/2026-11-13-e6f189046-P5.json`. Account: `strong-shock-run.md` Sect. 7.
- **Verdict: made.** Valid bins at 7 polar angles (52.5°–87.5°, θ_Bn 48°–87°), all driven and strong; valid
  and high (fast Mach ≥ 3) up to 2.09 / Ω_ci; median fast Mach of the valid bins 3.25. All end by
  criterion (e), the box, with the fast Mach number still 3.0–3.1.
- **Predictions.** E1 time missed (3.86 against 2.8–3.4 / Ω_ci), radius met; E2 missed (3.32 against 3.6–3.9);
  E3 met; E4 below (0.83 and 0.74 of V_d at R_m and 1.2 R_m): the front speed did not follow R / R_m as in
  L and S25; E5 missed (0.50) and was built wrongly (m from the centre of a front that starts at 0.65 R_m).
- **Seen.** The fast Mach number was highest (3.6–4.0) before the front was valid; the lead criterion (b) sets
  the start. 40 of the 96 valid bins are beyond the luminal angle (71–75°).
- **Slip.** The config dumps `hot_e` every 72 steps and the slices every 36; `census` needs them at the same
  times and refused. O8 waits for an adapted census.
- **Literature.** 37 + 48 notes committed (spot-checked against OpenAlex). The piston-shock survey: formation
  in about 1–2.5 / Ω_ci at 1–2.5 ρ (Schaeffer 2020), strong shocks in the literature all have a sustained
  source, finite clouds gave M_A ≈ 2. P5 is a finite shell that gave fast Mach 3–3.6 while it was driven.

## 2026-11-13 — Strong shock first: run P5 designed (before the run); S25 read in part

- **Instructions of the owner (two).** (1) "Why are the shocks still weak? Start from making a strong shock on
  a sphere." (2) Works that cannot be read in full are read by their abstract, so that their outline is known.
- **The fact that changes the plan** is instruction (1). The next run is designed for the strength of the
  shock: `strong-shock-run.md`, config `configs/strong-P5.yaml`, runner `run_piston.py`. Sections 1–6 of that
  document (reasons, design, observables, predictions E1–E5, verdict rule, limits) are committed before the run.
- **Why weak** (Sect. 2 there). The planner ties the box to the stopping radius, so with the box fixed
  R_m / ρ_d falls as M_A^(−5/3): a high M_A gave a cloud that does not couple (M, H), a low one a weak shock by
  construction (L, S25).
- **The lever.** A heavy piston: R_m = 140 d_e = 7 ρ_d in the box of L, M_A = 5 (V_d = 0.357 c), the debris
  does not stop inside the box, and the run ends when the front reaches the edge of the usable box.
- **S25, read so far** (run from `66a3629`: 2960 steps, 1310 s, nine checks passed; `front` and `drive` run
  with the code of `c122e4d`; census done, tracks still running, nothing recorded yet).
  - First valid front at 3.24 / Ω_ci, R = 71.0 d_e = 0.88 R_m; debris edge at 0.79 R_m and at 0.66 of the front
    speed. Valid and driven bins last 2.78 / Ω_ci at the longest: by the rule fixed before the run,
    **undecided** (between 2.5 and 4). Fast Mach in the valid bins 1.54 / 1.74 / 2.03: weak, as expected.
  - Reading P (the stalling of the debris) is rejected. Readings T and R are still not separated.
  - Used for P5 (Sect. 4 there): in L and S25 the front speed over V_d follows R / R_m, and the front becomes
    valid at 3.2–3.3 / Ω_ci, 1.9–2.3 ρ_d outside the initial debris edge.
  - The full account goes into `shock-forming-run.md` Sect. 7 with the record.
- **A slip in the order.** `drive.py` and `record.py` were extended (the count of bins with fast Mach ≥ 3,
  commit `f09d19a`) before S25 was recorded. The fixed observables of S25 do not change; `drive` is run again
  on S25 with the committed code before its record is written.
- **Literature, in progress.** Three searches run: abstract-level notes for the works listed as cited and
  not read, for those found and not read, and a survey of piston-driven strong shocks. Their notes are checked
  and committed when they arrive. If one of them gives a fact against the design of P5, that is a new entry.

## 2026-11-13 — The owner narrows the theme: the acceleration mechanism at a shock on a sphere

Written while S25 runs (started 12:44 from commit `66a3629`); none of its outputs has been read yet.

- **Decision of the owner:** the theme is narrowed to the acceleration mechanism at a shock on a sphere.
- **What it settles** (README, Research question): the formation of the shock and the coupling of the debris
  are conditions of the experiment, not subjects. Our reading, not yet confirmed: the subject is what the
  sphere adds (transport along the front between obliquities, drifts along it, mirroring in the draped field,
  expansion and deceleration); planar mechanisms that act locally are the baseline.
- **What changes, and the fact that changes it.** The fact is this decision.
  - The second run named in `shock-forming-run.md` Sect. 6 (cells of 1 d_e, electron scales not resolved) is
    withdrawn as the next run: it is an observation of ion scales and cannot show how electrons gain energy.
    The next design starts from the narrowed question, not from the drive.
  - S25 keeps its role and its fixed observables. It answers a condition of the theme: whether this box can
    hold a shock that is still driven. Its strength is expected below the first critical Mach number, so it
    is not expected to show acceleration.
  - **Added to the analysis of S25, before any output was read:** O7, the registered measures on the tracked
    electrons (`analyze.py tracks`: reflections, encounters, the change of the foot θ_Bn between gain and
    escape, K1–K3 of the protocol), as counts beside the census O6. No residual is claimed, as for O6.
- **What does not change:** the two predecessors are read in full before any hypothesis, and the question
  still needs a controlling number and a threshold. Both full texts are still not obtained.

## 2026-11-13 — Run S25 designed: a longer drive (before the run)

- **Instruction of the owner:** first run a condition in which a shock forms, and look at it. This changes the
  order recommended in `exploratory-results.md` Sect. 8 (read the two predecessors, narrow the question, then
  design). The readings and the narrowing stay conditions before any hypothesis.
- **New command `analyze.py drive`** (`drive.py`), run on the recorded outputs of L and M before the design:
  - L: first valid front at 3.30 / Ω_ci with the debris edge at 1.00 R_m and at 0.57 of the front speed;
    35 of 241 valid bins are driven (edge at half the front speed or faster), longest 1.32 / Ω_ci.
  - M: first valid front at 2.70 / Ω_ci with the debris edge moving inwards; no driven valid bin.
  - So both shocks became valid as free, decaying fronts. Criterion (b), the lead over the debris edge, is the
    last to hold; the speed criterion (c) holds from the first output that has a speed.
- **Design.** In this box only M_A moves R_m / ρ_d = (R_B / d_i) M_A^(−5/3). S25: M_A = 2.5, `rm` 20.1,
  R_m / ρ_d = 8.0 (L: 4.6), nominal fast Mach 2.27 (L: 2.95), β = 0.25. Everything else as L.
- **Fixed before the run** (`shock-forming-run.md`, Sect. 5): observables O1–O6, the predictions of three
  readings of where the front becomes valid (a time, a length, the stalling of the debris), and the verdict
  rule (driven valid bins for 4 / Ω_ci or longer at one polar angle: a driven shock formed; 2.5 or less: no
  better than L).
- **Not done:** no review before the run (an observation with the registered front definition; nothing is
  claimed), no null run for this condition.
- **Note on a number.** `exploratory-results.md` gives "at most 2.6 / Ω_ci" for the time L holds a fast Mach
  number of 2 or more at one polar angle. `drive` gives 2.31: it counts k outputs in a row as k − 1 intervals.

## 2026-11-13 — Review R2 of the results: verdict stands; corrections to the entries on L and on M, M-null and H

- **What was reviewed.** The first version of `exploratory-results.md` (commit `69459cd`), by the research
  reviewer, against the protocol, this notebook and the records. Report, unchanged:
  `reviews/2026-11-13-results-r2.md`. It is advice, not a decision of the owner.
- **Outcome.** "Not obtained from these runs" does not need rethinking. Nothing moves the answer towards
  "candidate". Nine changes were required before a merge; all are taken (table in `exploratory-results.md`,
  Sect. 10). The revised text was not reviewed again.
- **How the points were handled.** Each finding was checked against the records before it was accepted. New
  checks were added as commands (`r2`, `r2null`, `census --tau` of `analyze.py`, commit `aac5d4a`) and recorded
  (blocks `r2`, `census-tau`, `r2null` of the records, commit `f22fb3b`). They are not registered tests: they
  were written after the data were seen.

**Corrections to the entry "Run L" (below; not rewritten).**

- **Census series.** The entry gives "486 at t = 74, 330 at t = 148, 40 at t = 222, at most 4 later". The dumps
  between were left out, and the largest is not t = 74. From `census.json`: 486 (t = 74), **807 (t = 92.5)**,
  753, 533, 330 (t = 148), 207, 110 (t = 185), 59, 40 (t = 222), 21, 13, 10, then at most 6.
- **Upstream counts.** The entry gives "3 at t = 222 …, none at any other dump". Wrong: 22 electron-dumps
  upstream of the fitted surface, in eleven dumps (3 each at t = 166.5, 185, 222, 259, 277.5; 2 at 203.5; 1
  each at 240.5, 314.5, 333, 481, 610.5). 16 have a foot on a detected front, 3 are in the valid domain. By
  result: 5 explained by S0, 17 undetermined, no residual. The verdict is unchanged: no band holds 100.
- **"2, 0 and 0 measured upstream".** These are the baseline's counts at the three display times t = 185, 296
  and 407 only. Over all dumps the baseline counts the same 16 electron-dumps.
- **K1.** "One electron (t = 222)" is one electron in seven dumps: id 391104754, t = 166.5–277.5, γ falling from
  1.84 to 1.50. Not in the valid domain; inside the debris edge. Its field line was not traced.
- **N3.** The entry evaluates ν τ for the electrons of γ 1.25–1.5 with the τ_dwell of B1′ ("0.062 at the
  protocol's 66") and with the measured dwell. The registered N3 is for thermal electrons and the crossing
  time of the layer. The numbers are in the block below. The reading "N3 would pass at 66" is withdrawn.
- **"Most returns in L are not adiabatic mirror reflections …; μ says the same (12 % keep it)".** Withdrawn as a
  statement about the plasma. In L-null 11 % of the electrons keep μ within 20 % over a gap of 250, with no
  shock at all (block below).
- **"The ramp width was not measured".** It is measured now, as an upper bound (block below).

**Corrections to the entry "Runs M, M-null and H" (below; not rewritten).**

- **Bins before the valid domain.** "62 in L, 54 in M, all fail criterion b only: the front is not yet 2 d_i
  ahead" is not reproduced by the records and is wrong. Counting the bins with a front before the first valid
  time of their polar angle, at the polar angles that become valid: 74 (L) and 73 (M). All fail (b); 24 and 18
  also fail (c); 28 and 35 also fail (e); (b) alone: 46 and 38. Those that pass (c) have fast Mach
  1.92 / 2.50 / 2.78 (L, 50 bins) and 2.26 / 2.85 / 3.28 (M, 55 bins). The selection that gave 62 and 54 could
  not be recovered.
- **Criterion (b)** asks for a lead of 0.7 V_n/Ω_ci over the debris edge: 1.5–2.3 d_i in L, 2.3–3.3 d_i in M.
  It is not a fixed 2 d_i.
- **Equator table, M, t ≥ 444.** The dashes for the debris radius were a gap of that table. The debris edge
  exists in every output of M (62–64 d_e at 0.1 n₀); at t = 444 the lead is 32 d_e against 7.5 required.
- **"At t = 74, 486 (L), 14 (M), 0 (H) … their number falls with the magnetisation".** Two errors. The
  comparison at equal t mixes different phases, and the wording has the wrong sign: the number falls as the
  field falls. At the largest dumps: 807 (L, t = 92.5 = 1.65 Ω_ci⁻¹), 33 (M, t = 92.5 = 0.96 Ω_ci⁻¹), 1 (H,
  t = 111). At equal Ω_ci t ≈ 1.3: 486 (L) and 18 (M).
- **"The absence of energetic electrons is what a fast Mach number of 1.5–2.5 and b ≈ 2 give".** Withdrawn. The
  registered baseline predicts a beam of 10²–10³ at these same fronts. The beam is absent, but in the valid
  domain 98 % or more of it is predicted in the bands θ_Bn = 60–80°, where ε > 0.3 and the baseline is not
  defined, and the thermal N3 fails. The absence is a statement about the baseline and the floor.
- **N3 of M** ("0.051 at the protocol's 113"): the same error as in L; see the block below.
- **Small numbers.** Ω_ce d_i / V_n of M is 6.28, not 6.2. φ of M is 0.014–0.023 where measured; the lower end
  0.009 quoted before is the bin at 57.5°, which is flagged not measured.
- **Confirmed, not corrected.** The 5 upstream electrons of M at t = 111 "on rays without a detected front":
  all 22 upstream electron-dumps of M have feet in bins without a detected front and are undetermined.

**New numbers (records, blocks `r2`, `r2null`, `census-tau`).**

- **Ramp width of |B|** (10–90 % of the rise, upper bound; 10th / 50th / 90th percentile of the valid bins):
  2.6 / 3.6 / 4.8 d_e (L), 4.8 / 5.4 / 6.0 d_e (M). Crossing time: 24 (L), 49 (M) at the median.
- **N3 as registered** (thermal ν at lags 64 … 4 times τ): crossing of the ramp 0.12–0.28 (L), 0.25–0.64 (M);
  per θ_Bn band 0.10–0.37, 0.24–0.66; census layer 0.27–0.62, 0.38–0.98; τ_dwell of B1′ 0.33–0.76, 0.57–1.47;
  median dwell of the returns 1.3–2.9, 1.3–3.3. Above 0.1 in every reading (the lowest value of L is at the limit).
- **μ in the null runs.** Share of tracked electrons with μ within 20 % after a gap of 0 / 3 gyro-periods /
  150 / 250 / 420: 27 / 14 / 12 / 11 / 11 % (L-null), 23 / 12.5 / 12 / 11 / 11 % (M-null). The returns of the
  shock runs: 11.8 % (L), 12.3 % (M). No difference.
- **The rule of τ_dwell** (protocol Sect. 3.2), not applied before: the 90th percentile of the dwell of the
  reflections with γ_out ≥ 1.25 is 198 (L, n = 7) and 286.5 (M, n = 3), above 66 and 113. Census with the new
  τ: at γ ≥ 1.5 no residual in either (L: explained by S0 5 → 8, undetermined 17 → 14; M: 22 undetermined in
  both). In the extra bin γ ≥ 1.25 the residual with an encounter falls from 2103 to 746 (L) and from 316 to
  120 (M).
- **Regime of the valid domain.** V_n / v_te = 1.46 (L), 1.04 (M) at the median; β = 0.50, 1.47; γ − 1 = 0.5 is
  45 and 89 times ½ m_e V_n². R ∝ t^0.55 (L), t^0.52 (M). b at θ_Bn ≥ 80°: 2.13 (L), 2.06 (M) against a
  perpendicular MHD jump of 1.98, 1.94 (adiabatic index 5/3).
- **ε per band.** At the assumed τ_dwell, ε ≤ 0.3 in 111 of 241 (L) and 47 of 180 (M) valid bins. In every
  valid bin where one reflection reaches γ = 1.5 (V_HT ≥ 0.447; 96 bins in L, 59 in M) ε is above 0.3.
- **Why the energetic population of L decays.** The 807 electrons of t = 92.5, followed by id: 220 above 1.5 at
  t = 111, 44 at t = 185, 0 at t = 296, while their median distance behind the front grows from 6 to 27 d_e.
  They are made at the piston stage, left behind and not replaced. The cause of the energy loss was not
  analysed.
- **Downstream branch of C1.** One macro-electron above the cut behind the front in L (16 outputs), none in M.
  Fails the count.
- **Reports of protocol Sect. 4.** In `exploratory-results.md`, Sect. 4. The debris edge of M moves by 11 d_e
  between the thresholds 0.03 and 0.3 n₀ (L: 3.7 d_e).

**What failed or was not done.**

- An edit of `README.md` was made while a background job still had a `record` step pending; `record` needs a
  clean tree. Caught before the step ran; the edit was parked outside the tree and restored afterwards.
- Not done: the count model of thermal escape (S1) per band; the exit θ_Bn of the returns (not in the
  records); the first critical Mach number; the sensitivity of B2 to the errors of V_n and θ_Bn; the field
  line of the K1 electron; a second review. L-s2 and N6 were not run (no residual, no K1 signal).
- References that the reviewer gave from memory are not cited anywhere. The reviewer's estimate of V_sh / v_te
  for young remnants is the reviewer's and is not used as ours.

**Reading.** The runs answer the registered question with "not obtained". What they add beyond that is a list
of what a design must fix first: a narrowed question with a controlling number, a baseline valid where ε ~ 1,
a thermal numerical floor, and a shock requirement in Mach number and duration (`exploratory-results.md`,
Sect. 8).

## 2026-11-13 — Runs M, M-null and H: no signal at γ ≥ 1.5; the measured shocks are weak and decaying, the Mach scan did not happen
- **Runs.** All 8 GPUs, seed 12345, 2960 steps, clean tree.
  - M: `configs/explore-M.yaml`, commit `ac33a41`, 1321 s (446 ms per step), 60.2 GB. Gauss law ≤ 4.5 × 10⁻¹⁵,
    div B ≤ 1.1 × 10⁻¹⁴, fill ≤ 0.669, energy budget +11.8 % (monitored). Record `runs/2026-11-13-ac33a41eb-M.json`.
  - M-null: `configs/explore-M-null.yaml`, commit `3f5c437`, 783 s, 15.8 GB. Gauss law ≤ 2.1 × 10⁻¹⁵, div B ≤
    1.8 × 10⁻¹⁵, fill ≤ 0.667, budget +12.8 %. Record `runs/2026-11-13-3f5c43715-M-null.json`.
  - H: `configs/explore-H.yaml`, commit `3f5c437`, 814 s (274 ms per step), 16.9 GB. Gauss law ≤ 4.7 × 10⁻¹⁵,
    div B ≤ 7.0 × 10⁻¹⁵, fill ≤ 0.668, budget +12.3 %. Record `runs/2026-11-13-3f5c43715-H.json`.
  - M has a valid domain (θ ≥ 45°), so M-null was run, as planned: the census and N3c of M need its floor.
- **Storage.** All outputs of the topic now sit on the NVMe (owner's instruction), under
  `/mnt/nvme/pypic/output/research/spherical-shock-injection/`. `output/research/spherical-shock-injection/<name>` is a link to it, so the paths in the
  records stay valid. Nothing was deleted.
- **Null floor of M (`null`).** No front in any of the 738 bins. T_e 0.0100 → 0.0109. No electron at γ ≥ 1.5 at
  any dump (N1 = 0). N3 at lags 4 / 16 / 64: ν = 1.30 × 10⁻² / 7.9 × 10⁻³ / 5.0 × 10⁻³ (v 0.1–0.2), 3.3 × 10⁻³ /
  2.0 × 10⁻³ / 1.4 × 10⁻³ (v 0.3–1.0). N5: the 99.9th percentile of |Δρ| rises from 2.07 to 5.96 d_e.
- **Front of M (`front`).** 326 of 738 bins show a front; 180 are valid, 121 usable for B1′. No front is
  detected below θ = 35°; the valid domain is θ = 45–90°, t = 259–703.
  - In the valid domain (min / 10 % / median / 90 % / max): V_n 0.089 / 0.094 / 0.106 / 0.130 / 0.147; fast
    Mach 1.50 / 1.57 / 1.75 / 2.13 / 2.38; V_n / v_A 2.14 / 2.27 / 2.55 / 3.11 / 3.53; b 1.62 / 1.74 / 1.94 /
    2.27 / 2.50 (largest cell 2.92). Lead over the debris 10–59 d_e (median 31). Ion flow within 30 % in 180 of
    180; fit rms 0.29 d_e. Superluminal feet for θ_Bn ≥ 78–85°.
  - L on the same statistic: V_n 0.116–0.204 (median 0.150); fast Mach 1.50 / 1.55 / 1.82 / 2.24 / 2.51;
    b 1.58–2.78 (median 2.0). The entry of L quoted a narrower range (1.75–2.06); this one is the full range.
  - **The designed Mach scan did not happen.** The design had fast Mach 2.95 (L) and 4.05 (M). The valid
    domains of both hold fast Mach 1.5–2.5 with a median of 1.8, and b ≈ 2.
- **Front history, L and M (new; from `front.npz`).** The debris edge stalls far inside R_B = 148 d_e: at 83 d_e
  = 1.3 R_m in L (about 80 d_e from t ≈ 370), at 62–66 d_e = 1.4–1.5 R_m in M (from t ≈ 200). After that the
  front is a free blast wave that decays. Equator:

  | t | L: R front | R debris | V_n | fast Mach | valid | M: R front | R debris | V_n | fast Mach | valid |
  |---|---|---|---|---|---|---|---|---|---|---|
  | 148 | 65.8 | 59.5 | 0.215 | 2.54 | no | 52.9 | 58.9 | 0.197 | 3.20 | no |
  | 222 | 80.7 | 67.5 | 0.188 | 2.22 | yes | 66.0 | 64.8 | 0.160 | 2.59 | no |
  | 296 | 94.0 | 74.2 | 0.175 | 2.06 | yes | 76.9 | 60.9 | 0.136 | 2.20 | yes |
  | 370 | 106.6 | 79.2 | 0.161 | 1.89 | yes | 86.3 | 62.9 | 0.118 | 1.90 | yes |
  | 444 | 118.0 | 81.8 | 0.149 | 1.76 | yes | 94.9 | – | 0.112 | 1.80 | yes |
  | 518 | 128.8 | 82.9 | 0.142 | 1.67 | yes | 102.8 | – | 0.106 | 1.70 | yes |
  | 592 | 139.0 | 83.2 | 0.134 | 1.57 | yes | 110.4 | – | 0.099 | 1.59 | yes |
  | 666 | 148.6 | 82.6 | 0.125 | 1.47 | no | 117.9 | – | 0.099 | 1.59 | yes |

  - The first valid outputs are t = 185–278 in L (R = 73–86 d_e, fast Mach 1.9–2.5) and t = 259–333 in M
    (R = 71–80 d_e, fast Mach 2.1–2.4). The bins that show a front before t = 300 but are not valid (62 in L,
    54 in M) all fail criterion b only: the front is not yet 2 d_i ahead of the debris. They have fast Mach
    1.9–2.8 (L) and 2.5–3.3 (M).
  - So the stronger phase exists, but it is the phase in which the front has not separated from the piston, and
    the registered definition does not read it as a shock. The front becomes valid just as the piston stalls.
    With R_m = 16 d_i (L) and 11.2 d_i (M) there is no room between the two.
- **Census of M (`census --null`).** Ambient electrons with γ ≥ 1.5 per dump (all classes | of them upstream):
  5 | 0 (t = 55.5), 14 | 0 (74), 33 | 0 (92.5), 32 | 5 (111), 18 | 4 (129.5), 9 | 3 (148), 6 | 3 (166.5), 3 | 2
  (185), 2 | 1 (203.5), then at most 1 per dump up to t = 425.5 and none later. Every upstream one is on a ray
  without a detected front (all undetermined). The null run has none.
  - **Result at the registered threshold: no band reaches 100 electrons, so there is no signal and no residual
    at γ ≥ 1.5 in M.** K1: 3 records of 2 electrons (t = 111; t = 148 and 166.5), before any valid front.
  - b_max 3.39 (cell 4.33), τ_dwell 113. N3b: ν = 4.5 × 10⁻⁴ for γ 1.25–1.5 (1229 pairs); 3 pairs above γ = 1.5,
    so the rate of the null run (2.0 × 10⁻³) is used there.
  - The tracks were made with a preliminary census (before M-null existed). The final census differs from it in
    9 records (`k1`, its threshold, `und`), none of a tracked electron, and the join reads no other changed
    array. The tracks were therefore not re-run.
- **Tracks of M (`tracks`).** 956 674 tracked: 791 967 never left the upstream, 61 517 transmitted, 57 709 in the
  layer at the end, 6964 back upstream after an encounter, 579 returned from downstream, 6510 start behind the
  front, 31 428 outside the range of the fit.
  - γ_max ≥ 1.25 on 188 orbits, **γ_max ≥ 1.5 on none**. Returns with γ_out ≥ 1.25: 3 (below the 30 of the
    protocol; f_FN upper bound 1.0). M cannot test C1 or C2 either. K2, K3, returns from downstream at γ ≥ 1.5:
    none.
  - N4 passes: 851 of 853 orbits close the work within 10 % (both failures start behind the front). Parallel
    share of the work: median 0.23.
  - Returns at all energies: 6850 with one encounter, 3967 reached the density front; μ within 20 % in 488
    (12 %, as in L); μ_out/μ_in 0.20 / 1.15 / 6.0 (10 / 50 / 90 %). Gain: median 0.022, 99 % 0.131.
  - Dwell in the layer: median 252, 90 % 414.5 (protocol: 113). The caveat of L holds (the layer is 4 d_i thick).
- **N3 with the measured dwell, M.** ν τ = 0.051 at the protocol's 113, 0.11–0.19 at the measured 252–414.5:
  above 0.1, less than in L (0.23–0.39). The rule of N3 still forbids assigning the loss of μ to a physical
  scattering mechanism.
- **N3c, M.** Largest ratio 1.25 in the valid domain (no note), 3.77 on the detected front outside it.
- **Potential, M (`phi --null`).** Mean over the outputs 0.009–0.023 (1–2.3 T_e) at θ = 37–72°, measured (error
  below 50 %) in 5 of the 5° bins there; not measured above 70°. Windowed: measured in 49 of 323 detected bins,
  23 of the 180 valid ones. About half the potential of L.
- **Baseline B2 against the census, M (`baseline`).**
  - γ ≥ 1.5: no cell with 100 measured electrons, so the rule of Sect. 3.4 does not apply. The model again
    predicts a beam that is not there: on the detected front 107–223 electrons at t = 129.5, 560–1017 at t = 296,
    587–1046 at t = 351.5, falling to 21–50 at t = 740 (lower–upper bound); in the valid domain 76–100 (t = 296)
    to 187–264 (t = 462.5). **Measured upstream of a detected front: 0 at every dump.**
  - γ 1.25–1.5 (extra bin, not registered): 13 cells with ≥ 100 electrons on the detected front (4 inside
    [lower/2, 2 × upper], 4 below, 5 above), none in the valid domain. At t = 499.5 by band (model lower–upper |
    measured): 30–40° 1.2–1.3 | 33; 40–50° 58–164 | 166; 50–60° 1264–2250 | 93; 60–70° 5629–18 773 | 11;
    70–80° 6650–9347 | 0. The pattern of L: measured at low θ_Bn, missing by a factor 500 or more at 60–80°.
    **No residual is discussed.**
  - Share of the encounters that end in a return, all energies, detected front:

  | θ_Bn | model, φ = 0 | model, measured φ | tracks (± binomial) | tracks, only returns from behind the front |
  |---|---|---|---|---|
  | 30–40° | 0.135 | 0.064 | 0.146 ± 0.011 | 0.065 |
  | 40–50° | 0.130 | 0.036 | 0.142 ± 0.005 | 0.073 |
  | 50–60° | 0.104 | 0.047 | 0.121 ± 0.004 | 0.067 |
  | 60–70° | 0.049 | 0.009 | 0.094 ± 0.004 | 0.056 |
  | 70–80° | 0.004 | 0.003 | 0.077 ± 0.004 | 0.045 |
  | 80–90° | 0 | 0 | 0.067 ± 0.007 | 0.041 |

  - Tracks return 7–15 % of the encounters, the model with the measured potential 0–6 %. At 70–90°, where
    the mirror returns almost nothing with or without a potential, the tracks return 7–8 % (4–5 % with a sample
    behind the front). The valid domain agrees (6–13 % against 0–4.5 %).
- **H (contrast only).** No front in any of the 738 bins. The ambient ion density reaches 1.5 n₀ in 3
  equatorial bins at the end (r ≈ 102–107 d_e); the debris edge is near 50 d_e. In 5.3 Ω_ci⁻¹ the compression
  does not become a shock by the registered definition. Ambient electrons at γ ≥ 1.5: at most 1 per dump
  (one at t = 111).
  - The early energetic ambient electrons scale with the field: at t = 74, 486 (L), 14 (M), 0 (H). They are
    made at the start of the driver, not at a shock, and their number falls with the magnetisation.
- **Readings.**
  - L and M agree in everything that was measured: no ambient electron population at γ ≥ 1.5 upstream; a
    predicted mirror beam that is absent; returns at all energies that exceed the mirror with the measured
    potential and that keep μ in 12 % of the cases; returns concentrated at low θ_Bn.
  - The absence of energetic electrons is what a fast Mach number of 1.5–2.5 and b ≈ 2 give. It says nothing
    about the designed Mach numbers, which no valid front reached.
  - The failure of S0 is a property of the baseline as much as of the runs: marginal adiabaticity
    (Ω_ce d_i / V_n = (m_i/m_e) v_A / V_n ≈ 7.6 in L, 6.2 in M), a decaying front, and ν τ above 0.1. The
    candidate readings of the entry of L (leakage, kinematic returns from a decelerating front, numerical
    scattering) all remain; M removes none and adds none.
- **State of the question after L, M and H.** No signal at the registered threshold in any run. K1: one
  electron in L, two in M (before its first valid front). K2, K3: none. No residual with a positive signature. By
  Sect. 9 of the protocol the second seed (L-s2) and the 384³ pair are not run.
- Next: the results document with the criteria C1–C6 and the verdict, then review R2.

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
