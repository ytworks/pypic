# Prior work — shock-alfvenic-turbulence

This file places the topic in its field. Each work has its own file in `research/literature/<citekey>.md`; here
they are related to the question. The search behind it: `survey-shock-alfvenic-turbulence.md`.

## Map

- **Fully kinetic shocks in pre-existing turbulence, nonrelativistic electron–ion** (the direct predecessors):
  - `fulat2023kinetic`: θ_Bn = 90°, B0 in the plane. Compressive turbulence, δn/n ≤ 10 %. No effect on electrons.
  - `fulat2025kinetic`: ICRC proceedings of the oblique case.
  - `fulat2026pic`: θ_Bn = 60°, B0 out of the plane, δn/n 15 %, δB/B ≈ 10 %. The non-thermal electron fraction rises
    from 0.14 to 0.23 %.
  - Both papers leave Alfvénic turbulence for later. The amplitude is capped by heating.
- **Fully kinetic shocks in turbulence, relativistic pair plasma:**
  - `bresci2023particle`: magnetic turbulence driven by an antenna; controls without turbulence and without a
    shock; per-particle work split. This is the design template.
  - `demidem2023relativistic`: density waves corrugate the front.
- **Field-line meandering, test electrons or test particles** (the picture to be tested):
  - Test electrons in fields from a hybrid code: `guo2010effect`, `guo2012particle`.
  - Test electrons in kinematic fields: `guo2015acceleration`.
  - Protons: `giacalone2005particle` (test particles in prescribed fields) and `giacalone2005efficient` (self-consistent).
  - Review: `guo2021shock`. Latest work of the group: `lo2026effect`.
  - The 2-D limit: `jokipii1993perpendicular`.
- **Hybrid shocks in turbulence (ions only, electrons a fluid):** `trotta2021phase`, `trotta2024study` (2-D,
  δB/B0 up to 2.1); `trotta2023three` (3-D, δB/B0 0.5 and 1); `nakanotani2022turbulence`; `nakanotani2021interaction`
  (current sheets).
- **Electron acceleration at shocks without imposed turbulence** (the baseline mechanisms):
  - Ripples and test electrons: `burgess2006simulations`, `trotta2019electron`.
  - SSDA: `katou2019theory`, `amano2020observational`, `lindberg2024electron`.
  - Whistlers: `riquelme2011electron`, `oka2006whistler`.
  - Surfing and Weibel: `matsumoto2017electron`, `bohdan2017electron` (2-D in-plane vs out-of-plane B0).
  - Low-Mach ICM shocks: `guo2014nonthermal`.
- **Front corrugation by turbulence (theory):** `jebaraj2025turbulence`.
- **Relativistic shocks in clumps (MHD plus test particles):** `morikawa2024particle`.
- **Turbulence in PIC without a shock** (null test and set-up):
  - Large-amplitude Alfvén waves: `barnes1974large` (exact |B| = const waves).
  - Parametric decay: `matteini2010kinetics`, `gonzalez2023particle`, `shankar2026particle`.
  - Driving by an antenna: `tenbarge2014oscillating`.
  - Acceleration by turbulence alone: `comisso2018particle`, `zhdankin2017kinetic`, `comisso2022ion`
    (nonrelativistic e-i, δB/B0 = 1).
  - Dependence on β: `parashar2018dependence`.
  - Turbulence fed to an obstacle without an inflow boundary: `behar2024impact`.
- **Observations and parameters:**
  - Interplanetary shocks: `trotta2025overview`, `kilpua2015properties`, `jebaraj2023relativistic`,
    `jebaraj2024acceleration`.
  - Bow shock: `lindberg2024electron`, `liu2026efficient`, `raptis2025revealing`.
  - Solar wind: `wilson2018statistical` (β), `cuesta2022isotropization` (correlation length), `chen2020evolution`
    (δB/B near the Sun).

## Still to read
- Not reached:
  - Matthaeus et al. 2005 PRL 95, 231101 (correlation scale): only a summary was reachable.
  - Bruno & Carbone 2013 (Living Rev. Sol. Phys.): blocked.
  - Pitňa et al. 2021: abstract only, no δB/B values.
  - Ha 2025 (Ap&SS 370, 27, review of shocks in pre-existing turbulence): no abstract reachable.
- Linear theory of parametric decay for our β and amplitude: Derby 1978; Goldstein 1978.
- Statistics of δB/B0 upstream of interplanetary shocks at 1 au: the catalogue of `trotta2025overview` has the
  column, but no summary was read.
- Behar 2022; Behar & Henri 2023 (the method behind `behar2024impact`).
- Abstract-only notes to be read in full when the text is reachable: `giacalone2005particle`,
  `giacalone2005efficient`, `jokipii1993perpendicular`, `matteini2010kinetics`, `kilpua2015properties`, `oka2006whistler`,
  `nakanotani2021interaction`, `lo2026effect`, `morikawa2024particle`, `jebaraj2025turbulence`, `shankar2026particle`.
  `barnes1974large` rests on citing works only.

## Open
- Whether a fully kinetic (kinetic-electron) run of a nonrelativistic shock in pre-existing magnetic or Alfvénic
  turbulence exists. None was found in the search of 2026-11-13; see `survey-shock-alfvenic-turbulence.md`, Sect. 3.
- Whether electrons gain at all from meandering when they start from the thermal pool. The test-electron works
  start suprathermal: 100 eV in `guo2010effect`.

## This topic
The question asks for the fully kinetic test that the meandering picture lacks. Its parts relate to the works as
follows:
- **Electrons:** `guo2010effect` and `guo2015acceleration` follow test electrons in fields without electron
  scales. The topic follows kinetic electrons from the thermal pool, and test electrons in the same fields.
- **Kind of turbulence:** `fulat2023kinetic` and `fulat2026pic` impose compressive turbulence at kinetic scales.
  The topic imposes Alfvénic turbulence with δB/B0 0.3–1 and L_c ≥ 10 d_i.
- **Design:** the controls and the per-particle work split follow `bresci2023particle`.
- **Geometry:** `jokipii1993perpendicular` decides what a 2-D run can test.
- **Null test without a shock:** the acceleration by turbulence alone in `comisso2022ion` makes it necessary.
