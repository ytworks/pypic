# Electron injection at a nonrelativistic shock in large-amplitude Alfvénic turbulence

Branch: `research/shock-alfvenic-turbulence`. Started: 2026-11-13. Notebook: `notebook.md`. Survey: `literature.md`,
`survey-shock-alfvenic-turbulence.md`.
Sections follow `.claude/rules/research/stance.md`. Mark a section as open rather than guessing.

## Background
- **Electron injection problem.** Electrons must be pre-accelerated from the thermal pool before diffusive shock
  acceleration can take them. The known pre-acceleration mechanisms at quasi-perpendicular shocks are:
  - shock drift acceleration and its stochastic form, SSDA (`katou2019theory`, `amano2020observational`);
  - whistler-mediated injection (`riquelme2011electron`);
  - surfing and Weibel-mediated processes (`matsumoto2017electron`, `bohdan2017electron`).

  All of them were found with shocks that move into a uniform upstream.
- **Upstream turbulence is the rule.** Heliospheric shocks move into solar-wind turbulence:
  - δB/B0 ∼ 1 in the Alfvénic wind near the Sun (`chen2020evolution`);
  - a correlation length of 1.5–2 × 10⁶ km at 1 au, about 10⁴ d_i (`cuesta2022isotropization`);
  - β_e ≈ β_p ≈ 1 (`wilson2018statistical`).

  Electrons are seen accelerated at interplanetary shocks and at the bow shock: `jebaraj2023relativistic`,
  `jebaraj2024acceleration`, `lindberg2024electron`, and `liu2026efficient`, where an upstream discontinuity sets up the
  geometry that traps the electrons.
- **Field-line meandering.** In test-particle work, large-scale magnetic fluctuations let field lines cross a
  perpendicular shock several times. Electrons then return to the front and gain energy repeatedly
  (`guo2010effect`, `guo2015acceleration`, `guo2021shock`).

## Research question
Draft. It is to be fixed after the review and before any run.

> When a nonrelativistic supercritical shock moves into pre-existing Alfvénic turbulence with δB/B0 = 0.3–1 and a
> correlation length L_c ≥ 10 d_i:
> - how do the electron injection fraction, the spectrum and the maximum energy depend on δB/B0, θ_Bn and L_c, for
>   electrons followed self-consistently from the thermal pool;
> - is the gain due to field-line meandering across the front (the test-particle picture), or to other channels?
>
> The other channels are SSDA in the turbulent front, the corrugation of the front, and acceleration by the
> turbulence itself.

## Prior work
- **`guo2010effect`, `guo2012particle`, `guo2015acceleration`:** meandering field lines and ripples give fast
  electron energization at perpendicular shocks.
  - The fields come from a hybrid code (massless electron fluid) or are kinematic, and the electrons are test
    particles injected suprathermal (100 eV).
  - The shares of meandering and of ripples are not separated.
  - The authors list the stated limits: not self-consistent, no electron-scale shock structure, 2-D (`guo2010effect`, Sect. 2–3).
- **`fulat2023kinetic`, `fulat2025kinetic`, `fulat2026pic`:** the only fully kinetic nonrelativistic runs with
  pre-existing turbulence.
  - The turbulence is compressive and kinetic-scale: δn/n ≤ 15 %, δB/B ≈ 10 %.
  - At 90° there is no effect on electrons. At 60° the non-thermal fraction goes from 0.14 to 0.23 %.
  - Both papers state: "Future studies will explore the effect of Alfvénic turbulence."
  - Why it was left open, as they state it: heating caps the amplitude of compressive turbulence (M_s ≳ 30 must be
    kept), and their set-up imposes no magnetic fluctuation.
- **`trotta2021phase`, `trotta2023three`, `trotta2024study`:** hybrid shocks in MHD turbulence with δB/B0 up to 2.1
  (2-D) and 1 (3-D). They show corrugation and a laminar-to-turbulent transition, for ions only.
- **`bresci2023particle`:** the fully kinetic test with magnetic turbulence exists for relativistic pair shocks
  only. It has controls without turbulence and without a shock, and a per-particle split of the work.

## Idea
- **What makes the question answerable now:** a fully kinetic run that holds a thermal electron pool, an
  electron-scale shock structure and Alfvénic turbulence at δB/B0 ∼ 1 together. Then add controls that switch the
  meandering on and off while the turbulence keeps the same amplitude and spectrum:
  1. **Polarization control (2-D, B0 in the plane along the front, θ_Bn = 90°):**
     - fluctuations δB_x(y) carry field lines across the front;
     - fluctuations δB_z(y) of the same spectrum and amplitude move them out of the plane, along the front.
     - By `jokipii1993perpendicular`, the second cannot return electrons to the front by meandering.
  2. **Test electrons in the same fields:** test species of the engine, PIC-D9, pushed in the self-consistent fields
     of the same run. This tests the test-particle picture directly.
  3. **3-D:** no ignorable coordinate.
  4. **Per-particle work split** (`pypic.analysis.particles.work_by_component`): the motional field, the fluctuating
     field, and E∥.
- **What the difficulty is:** much of it is engine work and computing time, not an idea.
  - The current engine cannot load or inject an Alfvénic mode with δB_x in a piston box (Approach, E1–E4).
  - L_c is limited to 10–20 d_i, about 10³ times below the solar-wind value.

## Related work and novelty
- **Search:** 2026-11-13, `survey-shock-alfvenic-turbulence.md`.
  - Forward citations of eight seed works: 318 works, 95 screened.
  - About 12 arXiv queries, about 200 hits.
  - Web searches; observations and parameters.
  - Not done: asking experts; theses and 2026 meeting abstracts.
- **Novelty claim, as narrow as it is:**
  - No fully kinetic run of a nonrelativistic electron–ion shock in pre-existing magnetic or Alfvénic turbulence was found.
  - The claim would be a fully kinetic measurement of electron injection from the thermal pool in such turbulence,
    with the meandering channel isolated by polarization, dimensionality and test-electron controls.
- **Objections an expert would raise:**
  - L_c ≪ the solar-wind value, so the result holds for the inertial-range part only.
  - m_i/m_e = 100 and ω_pe/Ω_ce ∼ 10 (solar wind ≈ 150) distort electron-scale waves.
  - 2-D restricts cross-field transport (`jokipii1993perpendicular`).
  - Turbulence alone accelerates electrons (`comisso2022ion`), so a gain may not be due to the shock.
  - The turbulence decays and heats the plasma during the run (`parashar2018dependence`).
  - Typical interplanetary shocks have M_A ≲ 3, below the whistler critical Mach number (`trotta2025overview`).

## Approach
Draft. The numbers are fixed by the pilot runs; the criteria and null tests are fixed before the production runs.

**Geometry:** x is the shock normal, with a reflecting wall at x = 0 (piston). y lies along the front and is
periodic. z is out of the plane; in 3-D it is periodic.

**Engine prerequisites.** Each one is done on main: a branch from main, a pull request, a review by Codex (sol 6.1
ultra) repeated until it reports no findings, the gates of `pypic/`, then the merge, and `research-sync-main` into
trunk and this topic. Runs that need a change start after its merge (owner, 2026-11-13).

| ID | Change | Needed for |
|---|---|---|
| E1 | Non-uniform initial E tangential to a reflecting face (extends PIC-D18; today `initial_fields._check_tangential_E` refuses it) | loading modes with δB_x(y) in the whole upstream: E_z(y) = (v0 − δv_x) B0; pilot |
| E2 | Turbulent inflow face: time-dependent δB and E = −v × B at the injection face, and δv(y, z, t) on the injected particles (extends PIC-D13, a boundary drive like PIC-D20) | turbulence entering at a fixed age; the set-up of `guo2010effect` |
| E3 | Restart checkpoints (roadmap R2) | extending long runs; saving the state of a periodic turbulence box |
| E4 | A recorded turbulence slab fed through the inflow face (on E2 and E3) | relaxed turbulence (family T2) |

- An alternative to E2/E4 to be weighed in the plan of the engine work: turbulence decayed in a periodic box, then
  a reflecting wall moving through it in the plasma frame (`behar2024impact`). It needs a moving reflecting wall and E3.
- **Turbulence families:**
  - **T1:** a slab spectrum of parallel-propagating Alfvén waves, as in `guo2010effect`, linear or circular. Its
    lifetime is limited by parametric decay (onset about 150 Ω_ci⁻¹ at δB/B0 = 0.5, β = 0.25 in
    `gonzalez2023particle`; to be computed for our parameters).
  - **T2:** turbulence relaxed in a periodic box.
- **Runs (2-D production; to be fixed after the pilot):**
  - m_i/m_e = 100, β_e = β_i = 1.
  - M_A: 3–10 (heliospheric) and about 30, to connect to `fulat2026pic`.
  - δB/B0 = 0, 0.3, 0.6, 1.
  - θ_Bn: 90° and 60°.
  - Polarization: δB_x, δB_z, circular.
  - L_c: two values, 10 and 20 d_i.
  - B0 out of the plane, for the geometry of `fulat2026pic`.
  - Checks: m_i/m_e = 400, and two values of ω_pe/Ω_ce.
  - 3-D: a few runs at m_i/m_e = 25–50, with the largest box the engine holds (about 600³ cells).
- **Cost estimate:**
  - Box: 300 × 50 d_i at d_e = 10 cells, i.e. 1.5 × 10⁸ cells; 20 particles per cell in all, about 3 × 10⁹
    particles when full.
  - 20 Ω_ci⁻¹ at ω_pe/Ω_ce = 10 is about 3 × 10⁵ steps.
  - At the measured 1.1 × 10¹⁰ particle-steps/s and half fill on average, that is about 11 h per run.
  - The campaign occupies the node for weeks.
- **Controls and null tests** (draft, to be fixed before the runs):
  - N1: δB = 0 (planar baseline).
  - N2: the same turbulence without a shock, in a periodic box. It gives the electron gain from the turbulence alone
    (`comisso2022ion`).
  - N3: test electrons in the same run.
  - N4: a second seed.
  - N5: twice the particles per cell, and twice the resolution.
  - N6: the upstream spectrum, δB/B0 and temperatures recorded as functions of time and position.
- **Measured quantities** (thresholds to be fixed before the runs):
  - the downstream non-thermal electron fraction and its energy share, with the definition of `fulat2026pic`
    (γ − 1 ≥ 3) and one scaled to the shock speed;
  - the spectral index and the largest energy;
  - their spread over seeds;
  - the δB_x / δB_z difference;
  - the ratio of test to self-consistent electrons;
  - the share of each work term on the tracked electrons.
- **The claim would be shown wrong if:**
  - with δB_x and δB_z at equal amplitude, the injection fractions agree within the seed spread (meandering not the
    channel);
  - or the gain at the shock does not exceed that of N2 by more than the seed spread.

## Claims and evidence
None yet. No run has been made.

## Publication
Open. Intended: a peer-reviewed journal paper (ApJ or JGR).

## Status
2026-11-13:
- **Done:**
  - topic chosen by the owner (theme D);
  - branch cut;
  - survey with 51 literature notes (39 new, 12 copied verbatim from `research/spherical-shock-injection`);
  - this draft;
  - review (`reviews/2026-11-13-plan.md`).
- **Review verdict:** the novelty claim is plausible and needs evidence. The question, the polarization control, the
  N3 control, the 3-D part and the falsification criteria need rethinking (Majors 1–6, 9, 10).
  - The Approach above is the reviewed draft; it is not yet revised.
- **Missing before the question is fixed** (review, "Must change"):
  - one question, with a controlling number Ξ = ξ_x / L_sh and the access ratio v_te cos θ_loc / U1;
  - a minimal encounter/escape model, or a published number to falsify (`guo2012particle`);
  - a local-θ_Bn superposition null, an L_c (Ξ) scan and A_z field-line labelling, in place of the polarization
    control as the main evidence;
  - test electrons in ion-scale-filtered fields, in place of N3;
  - the Hall rotation of the polarization;
  - full reading of `jokipii1993perpendicular`, `giacalone2005*`, `barnes1974large`;
  - notes on Jokipii & Giacalone 2007, Qin et al. 2018, Xu et al. 2025;
  - an expert check;
  - the owner's decision on the question.
- **Engine work E1–E4:** each on main, after the question is fixed; their need may change with it.
