# Inputs of the four cases

Collected on 2026-11-14 (plan step 3). The values the runs and the transport model take from observations, each with
its source and its kind. The read range of each source is in its note (`research/literature/<citekey>.md`); (A) and
(S) mark notes on `research/shock-alfvenic-turbulence` and `research/spherical-shock-injection`.

Kind of a value:
- **measured**: from data of the event, by the authors of the source.
- **derived**: computed by the authors from measured values.
- **model**: from a model the authors name (solar wind propagation, shock surface, coronal reconstruction).
- **assumed**: a typical or adopted value the authors use because the event has no measurement.
- **own**: our calculation from the values above. Unverified; not a published number.

A quantity without a source is listed under "Open inputs" and is not guessed. No value here is a decision: the run
conditions and the criteria are fixed by the review of step 3 and the owner's approval.

## PSP, interplanetary shock of 2023-03-13 (0.24 au)

| Quantity | Value | Kind | Source |
|---|---|---|---|
| θ_Bn | 8° ± 4° | derived | `jebaraj2024acceleration` (A) |
| θ_Bn | 19.6° ± 9.1° | derived | `chiappetta2026evolution` |
| Upstream speed in the shock frame u_1 | 2400 ± 300 km/s | derived | `jebaraj2026wave`, `kouloumvakos2026wave` |
| Shock speed, spacecraft frame | 2800 ± 300 km/s | derived | `jebaraj2024acceleration` (A) |
| Shock speed | ~1350 km/s | model (hybrid) | `giacalone2025source` (abstract) |
| Shock speed from the 3-D reconstruction | ~2500 km/s at the nose (04:38 UT); 1800-1000 km/s (04:53 UT) | model | `chiappetta2026evolution` |
| M_A | 7-7.5 (v_A 320-360 km/s) | derived | `jebaraj2026wave`, `kouloumvakos2026wave` |
| M_A | 9.1 ± 1.35 | derived ("coarser estimate" in `kouloumvakos2026wave`) | `jebaraj2024acceleration` (A) |
| M_A | 4.68 ± 3.45 (M_ms 4.23 ± 3.13) | derived | `chiappetta2026evolution` |
| M_A | ~4 | model (hybrid) | `giacalone2025source` (abstract) |
| B_0 | ~80 nT; 90 → 140 nT across the shock | measured | `jebaraj2026wave`; `jebaraj2024acceleration` (A) |
| Density | ~30 cm⁻³ (31 → 38 cm⁻³ across the precursor) | measured | `jebaraj2026wave`, `jebaraj2026diffusive` |
| T_i | 90-100 eV | measured | `jebaraj2026wave` |
| T_e | "≲ 40 eV, unresolved over this interval" | not measured | `jebaraj2026wave` |
| T_e | 20 ± 5 eV | assumed (earlier PSP measurements) | `jebaraj2024acceleration` (A) |
| Ion β | ~0.2; 0.16 ± 0.03 | derived | `jebaraj2026wave`; `jebaraj2024acceleration` (A) |
| β (one value) | 0.28 ± 0.02 | derived | `chiappetta2026evolution` |
| Electron β for T_e 20-40 eV, 30 cm⁻³, 80 nT | 0.04-0.08 | own | — |
| Density compression | 4 | derived | `jebaraj2024acceleration` (A), `kouloumvakos2026wave` |
| Density compression | 1.39 ± 0.70 (SPAN-i); ~2.54 with the Langmuir-wave density | derived | `chiappetta2026evolution` |
| Magnetic compression | 1.78 ± 0.55 | derived | `chiappetta2026evolution` |
| Age at the crossing | eruption 03:13-03:15 UT, crossing 07:13 UT: ~4 h | measured | `kouloumvakos2026wave`, `chiappetta2026evolution` |
| ω_ci t for 4 h at 80 nT | 1.1 × 10⁵ | own (the field was stronger nearer the Sun) | — |
| Precursor (developed part) | ~58 min, ~12 R_s | measured | `jebaraj2026wave` |
| Wave amplitude δB/B | ~10 % total (0.05-1 Hz), 5-6 % resonant with 1-2.4 MeV protons | measured | `jebaraj2026wave` |
| Proton transport | L_EP ≈ 5.48 R_s, κ_eff ≈ 9.2 × 10¹⁵ m²/s | derived | `jebaraj2026wave` |
| Electrons | 50 keV-2.5 MeV with index ~2.7; isotropic to 6.7 MeV (near the instrument's limit) | measured | `jebaraj2024acceleration` (A) |
| Low-energy instrument | EPI-Lo saturated ~07:16 to ~08:55 UT | measured | `kouloumvakos2026wave` |
| v_sh/c | 0.0080 (u_1 2400 km/s) | own | — |

- Reason the authors give for the spread: `chiappetta2026evolution` attributes its low compression and M_A to "an
  underestimation of the gas density upstream of the shock from SPAN-i"; `kouloumvakos2026wave` calls the 9.1 of
  `jebaraj2024acceleration` (A) a coarser estimate. `giacalone2025source` was read as an abstract only.
- Along the heliocentric distance the same shock is quasi-perpendicular at Solar Orbiter (0.60 au: θ_Bn 57.6° ± 24.1°,
  compression 3.60 ± 0.43, M_A 4.81 ± 0.58) and at Wind (0.99 au: 71.5° ± 10.3°, 2.98 ± 0.17, 3.04 ± 0.20)
  (`chiappetta2026evolution`).

## Saturn, bow shock crossing of 2007-02-03 (Cassini)

| Quantity | Value | Kind | Source |
|---|---|---|---|
| B upstream | ~0.1 nT | measured | `masters2013electron`, `masters2017situ` |
| M_A | ~100 | model (Cassini data with a propagated solar wind; `masters2011electron`) | `masters2013electron` |
| M_s, M_f | ~25 | assumed (typical temperatures; `achilleos2006orientation`) | `masters2013electron` |
| θ_Bn | ~20° (~60° before 00:05 UT) | model normal (`went2011new`) with the measured field | `masters2013electron` |
| Upstream flow speed | ~450 km/s; ~500 km/s | assumed (typical at Saturn's orbit) | `masters2017situ` |
| Solar wind speed at the bow shock | 300-600 km/s | range of the parameter survey | `oka2025maximum` |
| Shock motion toward the spacecraft | 100 km/s | assumed ("large relative uncertainty") | `masters2017situ` |
| Density | not given for the event | — | — |
| Density implied by M_A 100, 0.1 nT, 450 km/s | 0.24 cm⁻³ | own | — |
| β (total) from M_A 100 and M_s 25, γ 5/3 | ~19 | own | — |
| T_e upstream | ≲ 10 eV; below 100 eV mixed with photoelectrons | measured | `masters2016suprathermal`, `masters2017situ` |
| Locally quasi-parallel connection | from ~00:05 UT until after the crossing | measured | `masters2017situ` |
| Thermal plasma transition | ~01:05 UT | measured | `masters2017situ` |
| Electron signature resolved | from ~00:35 UT; two power laws ~00:30 to ~01:20 UT | measured | `masters2017situ` |
| Connection at the spacecraft before the crossing | ≥ 60 min | own (00:05 to 01:05 UT) | — |
| ω_ci t for 60 min at 0.1 nT | ≥ 34 | own | — |
| ω_ci t from size / speed, 50-60 R_S and 300-600 km/s at 0.1 nT | 46-112 | own | — |
| Size of the shock | 50-60 R_S = (2.9-3.5) × 10⁹ m | survey value | `oka2025maximum` |
| Hillas limit | 0.2-2 MeV (for B 0.2-1.0 nT) | derived | `oka2025maximum` |
| Electrons | to ~1 MeV (top channel); indices −4 ± 1 below and −2 ± 1 above 100 keV | measured | `masters2013electron`, `masters2017situ` |
| Electrons, maximum as tabulated elsewhere | 700 keV | survey value | `oka2025maximum` |
| Seeds | electrons to ~100 keV upstream and downstream from 00:05 UT | measured | `masters2013electron` |
| Statistics | 31 of 508 crossings with a signal; 3 to ~1 MeV | measured | `masters2016suprathermal` |
| v_sh/c | 0.0015 (450 km/s) | own | — |

- The event papers give no measured solar wind speed, density or β. `masters2011electron` (abstract) states why: ion
  measurements upstream are scarce and there is no upstream monitor, so the speed comes from a model propagated from
  near-Earth data.
- The connection time at the spacecraft is a lower bound on how long the field line had touched the shock; the time a
  field line stays connected while it is swept along the shock is not given in any source read.
- The Hillas range of `oka2025maximum` uses B 0.2-1.0 nT; the event has ~0.1 nT.

## Earth, quasi-parallel bow shock without foreshock transients

| Quantity | Value | Kind | Source |
|---|---|---|---|
| M_A of bow shock crossings | most crossings at 5-15 (2797 MMS crossings, OMNI upstream) | derived | `lalti2022database` |
| M_A, one quasi-parallel example | 9.4 (θ_Bn 33°) | derived | `lalti2022database` |
| M_A as quoted in a scaling over planets | ~4 | stated without derivation | `raptis2026relativistic` |
| M (events of the foreshock study) | < 20 | derived | `wilson2016relativistic` |
| Share of quasi-parallel crossings | 45.5 % (quasi-perpendicular 51.4 %) | measured | `lalti2022database` |
| Solar wind speed | 400-600 km/s | survey value | `oka2025maximum` |
| B | 5-8 nT ("fast" wind) | survey value | `oka2025maximum` |
| One event | 4.3 nT, 4 cm⁻³, flow ~520 → 190 km/s | measured | `liu2026efficient` (A) |
| A simulated event | 5 nT, 6 cm⁻³, v_A ≈ 44.5 km/s | set-up values | `liu2019relativistic` |
| T_e | 12.06 eV | assumed (typical; used for M_f) | `lalti2022database` |
| Electron β for 12 eV, 5 cm⁻³, 5 nT | ~1 | own | — |
| Compression | approaches 4 at large M_A | measured | `lalti2022database` |
| Electrons outside the disturbances | "only noise > 12 keV" | measured (instrument limit) | `wilson2016relativistic` |
| Electrons inside disturbances | ≥ 30 keV for tens of seconds to a few minutes; to ~300 keV in 10 of 30 | measured | `wilson2016relativistic` |
| Model without transients | 10-20 keV | model | `shi2025compound` |
| Size of the shock | 15-25 R_E = (9.6-16) × 10⁷ m | survey value | `oka2025maximum` |
| Hillas limit; observed maximum | 190-760 keV; 300 keV | derived; survey value | `oka2025maximum` |
| Connection time | no source found | — | — |
| ω_ci t from size / speed, 15-25 R_E, 400-600 km/s, 5-8 nT | 76-305 | own | — |
| v_sh/c | 0.0013-0.0020 (400-600 km/s) | own | — |

- The "no" of this case is a detection limit near 12 keV, not a measured absence at lower intensity.
- θ_Bn of the planned run E (~25°) is a choice inside the quasi-parallel range, not a measured input.

## Interplanetary shocks near 1 au

| Quantity | Value | Kind | Source |
|---|---|---|---|
| 51 shocks of the electron search (mean distance 0.79 au), averages ± standard errors | θ_Bn 54° ± 3°; speed 512 ± 24 km/s (spacecraft frame); M_A 3.3 ± 0.3; M_fast 2.8 ± 0.3; gas compression 2.0 ± 0.1; magnetic compression 1.9 ± 0.1; upstream δB/B 0.25 ± 0.02; density 19.6 ± 3.4 cm⁻³; β 0.84 ± 0.12 | measured | `kartavykh2025statistical` |
| The 4 shocks with an electron response | θ_Bn 68° ± 9° (88°, 42°, 80°, 61°); 545 ± 81 km/s; M_A 6.1 ± 1.0 (7.4, 3.2, 5.3, 8.4); M_fast 4.6 ± 0.7; gas compression 1.9 ± 0.3; β 1.05 ± 0.34 | measured | `kartavykh2025statistical` |
| Uncertainty of single values | up to 20 % in θ_Bn or speed | stated by the authors | `kartavykh2025statistical` |
| 475 shocks | mean magnetosonic Mach number 1.6; five electron events at 1.25-2.85 | measured | `dresing2016efficiency` |
| Nine MeV-electron events of 587 shocks | M_ms 1.48-2.85 (mean 2.04); transit speeds > 900 km/s; θ_Bn ~45-88° | measured | `talebpoursheshvan2023role` |
| Quasi-parallel subset | 15 shocks; increase at 40 keV in ~7 % (quasi-perpendicular 59 shocks, ~17 %) | measured | `yang2019electron` (abstract) |
| Mach numbers of that subset | not available (abstract only) | — | — |
| Threshold for an effect at ≥ 2 keV | speed along the upstream field ≳ 250 km/s | measured | `tsurutani1985acceleration` (abstract) |
| ω_ci t for 1-3 days at 5 nT | 4 × 10⁴ - 1.2 × 10⁵ | own | — |

- The averages of `kartavykh2025statistical` are over all geometries; a table for the quasi-parallel shocks alone was
  not found. Each entry of its table has a second value besides the average; a check on the four listed M_A did not
  reproduce it as a median, so only averages and standard errors are used here.
- The plan's "M_A ~1.5-3" for the 15 quasi-parallel shocks of `yang2019electron` has no source in what was read. The
  measured values above are magnetosonic Mach numbers of 1.25-2.85 for other samples and M_A 3.3 ± 0.3 for the 51
  Solar Orbiter shocks.

## Remote context: coronal shocks
`dresing2022role` (abstract): in-situ peak intensities of ~1 MeV electrons correlate with the Mach number of the
coronal shock at the connected field line (coefficient ~0.7), electrons below 100 keV less (0.47); the correlations
vanish when the two spacecraft are exchanged. Not an input of the four cases; it bears on H7 (the history of the PSP
shock near the Sun).

## Open inputs
- **Instrument backgrounds** as intensities for each case (THEMIS, Cassini LEMMS, PSP IS⊙IS, Wind/STEREO/Solar
  Orbiter). Needed for the confrontation by intensity (H8). Not collected.
- **Earth:** a published connection time of a field line to the quasi-parallel shock; the 76-305 above is our size
  over speed.
- **Saturn:** density and the split of β into β_e and β_i; the solar wind speed and the shock motion of the event.
- **PSP:** T_e (unresolved in the data); the history of M_A, θ_Bn and B between the Sun and 0.24 au.
- **1 au:** Mach numbers and compressions of the quasi-parallel shocks of `yang2019electron`.
- **Convention of D_Bohm** and of the acceleration-time integral: `jebaraj2026diffusive` uses t_acc = ∫ (20 κ_∥/u_1²)
  d ln p for compression 4; `gupta2025speed` (S) uses D_B = u r_L. To be fixed in run V1.

## The plan's run conditions against the collected values
Facts only; no condition is changed here.

| Run | Planned | Collected |
|---|---|---|
| P0 (PSP) | M_A 7.5, β_i 0.2, T_e 20-40 eV, θ 8° | M_A 7-7.5 in two works, 9.1 ± 1.35, 4.68 ± 3.45 and ~4 in others; β_i ~0.2 and 0.16 ± 0.03; T_e not measured; θ_Bn 8° ± 4° and 19.6° ± 9.1° |
| P1, P3 (one change from P0) | M_A 5; θ 30° | M_A 5 lies inside the spread of the PSP estimates; θ 30° lies above both published ranges |
| E (Earth) | M_A 8, β ~1, θ ~25°, ω_ci t 240 | M_A 5-15 (most crossings), ~4 in one scaling; β_e ~1 from typical values (own); ω_ci t 76-305 (own) |
| S1-S3 (Saturn) | M_A 100, θ 20°, β ~20, ω_ci t 50 | M_A ~100 (model), θ_Bn ~20° (model normal), β ~19 (own, from assumed M_s); ω_ci t ≥ 34 at the spacecraft, 46-112 from size over speed (own) |
| W (1 au) | M_A 3, β ~1, θ 30° | M_A 3.3 ± 0.3, β 0.84 ± 0.12, gas compression 2.0 ± 0.1 over all geometries; no values for quasi-parallel shocks alone |
