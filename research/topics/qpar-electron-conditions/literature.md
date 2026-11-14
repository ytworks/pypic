# Prior work — qpar-electron-conditions

This file places the topic in its field. Each work has its own file in `research/literature/<citekey>.md`; here
they are related to the question. The search behind it is recorded in `notebook.md` (entry of 2026-11-14).

Notes marked (S) are on `research/spherical-shock-injection`, notes marked (A) on
`research/shock-alfvenic-turbulence`. Neither branch is merged into trunk yet, so those files are not on this branch;
they are not copied, to keep one file per work.

## Map

- **Events with relativistic electrons at a quasi-parallel shock:**
  - Saturn (Cassini): `masters2013electron` (2007-02-03, M_A ~100, to ~1 MeV), `masters2016suprathermal` (31 of 508
    crossings, 3 to ~1 MeV), `masters2017situ` (a quasi-perpendicular crossing with the same index above 100 keV).
    `bebesi2019observations`: SLAMS at Saturn's bow shock (conference abstract only).
    Where the shock parameters of the Saturn event come from: `masters2011electron` (M_A from a propagated solar wind
    model), `went2011new` (the model normal behind θ_Bn), `achilleos2006orientation` (typical temperatures; shock
    speeds up to ~400 km/s). All three read as abstracts.
  - PSP 2023-03-13: `jebaraj2024acceleration` (A) (θ_Bn 8°, electrons to 6.7 MeV), `jebaraj2024direct` (compared with
    the quasi-perpendicular 2022-09-05 shock), `dresing2025reason` (flank crossing, saturated instrument).
  - Jupiter (Juno): `raptis2026relativistic` (≥ 1 MeV inside a foreshock transient).
- **Events without them, or only with help:**
  - Earth: `wilson2016relativistic` (30 foreshock disturbances, 10 with electrons to ~300 keV; noise above 12 keV
    outside), `liu2019relativistic` (100-400 keV inside a foreshock bubble), `raptis2025revealing` (S) (> 500 keV only
    with 1-5 keV seeds), `raptis2025multi` (electrons confined in transmitted hot flow anomalies),
    `shi2025compound`, `liu2026efficient` (A), `qin2026electron` (abstract), `bai2023electron` (betatron gain inside
    a SLAMS; abstract). The parameter range of the bow shock itself: `lalti2022database` (2797 MMS crossings, most
    at M_A 5-15).
  - Interplanetary shocks, statistics: `tsurutani1985acceleration`, `dresing2016efficiency`, `yang2019electron`,
    `talebpoursheshvan2023role`, `kartavykh2025statistical`. Electron effects are rare, of low energy, and mostly at
    oblique or quasi-perpendicular shocks.
  - A planar interplanetary shock without mature SLAMS: `raptis2026compressive`.
  - MeV electrons inferred at a weak, oblique shock: `wilson2025large` (abstract). A near-relativistic electron
    peak at a shock at PSP, 0.81 au: `mitchell2021energetic` (abstract; geometry not given there).
- **The PSP shock as a transport laboratory (protons):** `jebaraj2026wave` (precursor ~12 R_s, κ_eff, wave amplitudes),
  `kouloumvakos2026wave` (confinement boundary ~2.46 MeV), `jebaraj2026diffusive` (wave growth, not age or size, limits
  confinement; preprint), `giacalone2025source` (hybrid model of the event; abstract), `chiappetta2026evolution`
  (other shock parameters for the same crossing).
- **Fully kinetic simulations of quasi-parallel shocks with electrons:**
  - The direct predecessors: `gupta2024electron` (S), `gupta2025speed` (S); proceedings `gupta2023regulates`.
  - The opposite dependence on M_A in 1-D parallel shocks: `shalaby2022mechanism` (S); mass ratio 100 against 1836:
    `shalaby2024energy`.
  - Very high M_A, SLAMS: `zekovic2025slams` (S), `hemler2025shock` (test particles in prescribed SLAMS; abstract).
  - Single conditions with reduced mass ratio: `kato2015particle` (S), `park2015simultaneous` (S), `fang2019early` (S),
    `fang2022kinetic`, `yu2024nonthermal`.
  - The transition of Earth's shock in 2-D, with reconnection: `bessho2023electron` (abstract).
  - Long-term self-regulation: `galishnikova2026self` (preprint; abstract).
- **Measuring a diffusion coefficient in a kinetic run:** ions from the upstream profile, `caprioli2014simulationsb`;
  electrons from the mean square displacement in unmagnetized shocks, `grassi2023electron`; test electrons in field
  snapshots, `gupta2025speed` (S); test electrons in hybrid fields, `hanusch2020electron`.
- **Models that reproduce an event with prescribed fields:** `guo2015acceleration` (S) (Saturn), `shi2025compound`
  (Earth, quasi-perpendicular with transients).
- **Remote evidence for the shock near the Sun:** `dresing2022role` (peak intensities of ~1 MeV electrons against the
  Mach number of the connected coronal shock; abstract).
- **Limits by size:** `oka2025maximum` (Hillas limit over heliospheric shocks), `raptis2026relativistic` (Bohm
  scaling over three planets).
- **Remote:** `zhou2023magnetic` (S) (SN 1006, the field behind the shock is nearly radial).

## Open

- The electron injection rate and D̃_e at the conditions of an observed event. `gupta2025speed` (S) measured D_e at
  θ_Bn 30°, M_A 5 and 20, and names weak and more oblique shocks as future work.
- The dependence of D̃_e on v_sh/c. `gupta2024electron` (S) varied the piston speed for spectra and efficiency only;
  `jikei2026dependence` (S) cites that work for a possible efficiency minimum near 0.1 c.
- Whether the low-M_A stall of `gupta2024electron` (S) is a matter of run time ("at least" within 300 ω_ci⁻¹). The PSP
  shock has M_A 7-7.5 (or ~4) and accelerates electrons to MeV.
- Which of the two orders in M_A holds with the real mass ratio at event conditions (`gupta2024electron` (S) against
  `shalaby2022mechanism` (S), `shalaby2024energy`).
- Electron transport at the PSP event: the works above treat protons; `jebaraj2024acceleration` (A) leaves open how
  sub-relativistic electrons become relativistic.
- Saturn: what produces the turbulence that `guo2015acceleration` (S) has to assume (variance 10 B_0²), and the role of
  the ~100 keV magnetospheric electrons upstream (`masters2013electron`).
- Earth without transients: the absence is an instrument limit (`wilson2016relativistic`); the geometry of the
  efficient events is not separated from position in `qin2026electron` (abstract).
- Whether the size limit with Bohm scattering and "sufficient time" (`raptis2026relativistic`) holds once D_e and the
  age are measured rather than assumed.

## This topic

- Adds: the injection rate and D̃_e(x, p̃) of electrons at event conditions with their dependence on v_sh/c, from
  self-consistent kinetic runs; a transport model fed with them, confronted with four cases.
- Reproduces first, and compares against:
  - `gupta2024electron` (S): slope and p_max(t) of run B1; the stall of run B2_II near p ≈ m_i v_pt.
  - `gupta2025speed` (S): D_e ≈ 0.6-1 Bohm just upstream and a tenth of it downstream at M_A 20; the p_max formula.
  - `caprioli2014simulationsb`: the ion D(p) from the upstream profile, as a check of the method at p̃ ≳ 1.
- Numbers of the events:
  - PSP: M_A 7-7.5, u_1 ~2400 km/s, B_0 ~80 nT, ion β ~0.2, T_e ≲ 40 eV, wave amplitude ~10 % (5-6 % resonant),
    κ_eff ~9.2 × 10¹⁵ m²/s, precursor ~12 R_s (`jebaraj2026wave`); electrons to 6.7 MeV with index ~2.7 between 50 keV
    and 2.5 MeV (`jebaraj2024acceleration` (A)); the other parameter sets (`dresing2025reason`,
    `chiappetta2026evolution`, `giacalone2025source`).
  - Saturn: B ~0.1 nT, M_A ~100, θ_Bn ~20°, ~1 MeV, indices −4 ± 1 below and −2 ± 1 above 100 keV
    (`masters2013electron`); 3 of 508 (`masters2016suprathermal`).
  - Earth: noise above 12 keV without disturbances (`wilson2016relativistic`); 10-20 keV in the model without
    transients (`shi2025compound`).
  - 1 au: ~7 % of 15 quasi-parallel shocks at 40 keV (`yang2019electron`); 5 of 475 at 65-75 keV
    (`dresing2016efficiency`); effects at 2-10 keV (`tsurutani1985acceleration`).
- All of these, with source and kind (measured, model, assumed, own), are in `inputs.md`.
- Conflicts to settle before the runs, as `inputs.md` now states them:
  - PSP: M_A 7-7.5, 9.1 ± 1.35, 4.68 ± 3.45 or ~4, and speeds of 2400-2800 or ~1350 km/s, depending on the work;
    T_e is not measured.
  - Saturn: M_A, θ_Bn and the speed are model or typical values; density and β_e are not given; the connection is
    bounded from below by the 60 min seen at the spacecraft.
  - Earth: no published connection time was found; M_A 5-15 against "~4" in `raptis2026relativistic`.
  - 1 au: no Mach numbers or compressions for the quasi-parallel shocks alone.
  - Conventions for D_Bohm and the acceleration-time integral (`gupta2025speed` (S), `jebaraj2026diffusive`).
