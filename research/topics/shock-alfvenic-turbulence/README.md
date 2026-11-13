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
- **Field-line meandering.** Large-scale magnetic fluctuations let a field line cross a perpendicular shock several
  times. Fast electrons then return to the front along it and gain energy repeatedly. The mechanism is from
  `jokipii2007adiabatic`; test-particle work followed (`guo2010effect`, `guo2015acceleration`, `guo2021shock`), with
  a conflicting result in `qin2018effects`.

## Research question
Proposed 2026-11-13 after the review (`reviews/2026-11-13-plan.md`, Major 6). It is the owner's decision and is not
fixed yet.

> At a supercritical perpendicular shock (M_A ≈ 4.5, β ≈ 1) moving into large-scale incompressible Alfvénic
> fluctuations, with electrons followed self-consistently from the thermal pool:
> - does the multiple connection of convected, meandering field lines to the front (`jokipii2007adiabatic`,
>   `guo2010effect`) inject electrons beyond what the distribution of local θ_Bn of a laminar shock gives;
> - and how does that excess scale with Ξ = ξ_x / L_sh and with the access ratio set by the mass ratio?
>
> Here Ξ is the peak-to-peak excursion of a field line across the front over the shock thickness.

The first draft asked two things: a scan over δB/B0, θ_Bn and L_c, and an attribution of the gain to mechanisms.
The review found it was two questions; the draft is kept in the history of this file.

## Prior work
- **The mechanism: `jokipii2007adiabatic`.**
  - Fast, low-rigidity particles follow large-scale meandering field lines back and forth across a shock that moves
    normal to the mean field.
  - They gain energy from the difference of the flow speeds, without pitch-angle scattering.
  - The paper proposes this as a pre-injection stage for suprathermal electrons. The note is from the abstract only.
- **Test electrons in fields from a hybrid code (massless electron fluid) or in kinematic fields:
  `guo2010effect`, `guo2012particle`, `guo2015acceleration`.**
  - Meandering field lines and ripples give fast electron energization at perpendicular shocks.
  - In `guo2010effect` the electrons are injected suprathermal (100 eV). In `guo2012particle` they are thermal: about
    10 % exceed 15 keV at M_A ≈ 2, δB²/B0² ≈ 0.3.
  - The shares of meandering and of ripples are not separated.
  - The authors state their limits: not self-consistent, no electron-scale shock structure, 2-D (`guo2010effect`,
    Sect. 2–3).
- **A conflicting test-particle result: `qin2018effects`.**
  - Static slab plus 2-D field, planar tanh shock, M_A = 10, electrons at 1 keV.
  - At a perpendicular shock the efficiency falls as (b/B0)² rises from 0.01 to 10.
  - The authors attribute the difference from Guo & Giacalone to the turbulence model and do not test it.
  - Reading recorded in the note (not the authors'): their fluctuations are static in the shock frame and not frozen
    in, so the multiple connection of convected lines is absent.
  - The test-particle picture is therefore not one prediction.
- **Fully kinetic, nonrelativistic, pre-existing turbulence: `fulat2023kinetic`, `fulat2025kinetic`, `fulat2026pic`
  (ApJ 1005, 107).**
  - The turbulence is compressive and kinetic-scale: δn/n ≤ 15 %, δB/B ≈ 10 %.
  - At 90° there is no effect on electrons. At 60° the non-thermal fraction goes from 0.14 to 0.23 %.
  - Both papers state: "Future studies will explore the effect of Alfvénic turbulence."
  - Why it was left open, as they state it: heating caps the amplitude of compressive turbulence (M_s ≳ 30 must be
    kept), and their set-up imposes no magnetic fluctuation.
- **Laminar kinetic baseline: `xu2020electron`.** 1-D PIC at 63°, M_A ≈ 63: repeated SDA cycles between
  self-generated waves and the shock, then DSA. A turbulence effect is measured against the laminar shock at the same
  parameters.
- **Hybrid shocks in MHD turbulence: `trotta2021phase`, `trotta2023three`, `trotta2024study`.**
  - δB/B0 up to 2.1 (2-D) and 1 (3-D).
  - They show corrugation and a laminar-to-turbulent transition, for ions only.
- **Fully kinetic test with magnetic turbulence: `bresci2023particle`.** It exists for relativistic pair shocks only.
  It has controls without turbulence and without a shock, and a per-particle split of the work.

## Idea
- **What makes the question answerable now:** a fully kinetic run holds a thermal electron pool, an electron-scale
  shock structure and convected Alfvénic fluctuations together. In 2-D with B0 in the plane, every field line is a
  contour of A_z (`jokipii1993perpendicular`). Its intersections with the front, and an electron's encounters with
  the front along it, can then be counted directly.
- **How meandering is separated from local obliquity** (review, Major 1):
  - **(a) Local-θ_Bn null.** Laminar runs at uniform θ_Bn = 90, 80, 70, 60 and 45°, weighted by the distribution of
    θ_loc measured at the front of the turbulent run. Multiple connection shows as the excess over this weighted sum.
  - **(b) Ξ scan.** At fixed δB_x/B0 the distribution of θ_loc stays the same, while Ξ grows with the box width L_y
    (`design/estimates.md`, Sect. 2): Ξ = 3.4, 6.9, 13.6, 27 for L_y = 50, 100, 200, 400 d_i at δB_x/B0 = 0.6.
  - **(c) Field-line labelling.** For every tracked electron: the number of front intersections of its A_z contour,
    the encounters made, and the gain per encounter, over the whole energized population.
- **The access ratio is a second control number** (`design/estimates.md`, Sect. 2–3).
  - t_conn / t_along: the time a line stays multiply connected, ξ_pp / U1, over the thermal transit between two
    intersections, (L_y / 2) / v_te.
  - It does not depend on L_y. It is about 0.66 at m_i/m_e = 100 and about 2.8 at 1836 (δB_x/B0 = 0.6, M_A = 4.5).
  - The fraction of thermal electrons that can stream upstream along the local field is 7 % at m_i/m_e = 100 and 29 %
    at 1836.
  - So the reduced mass ratio suppresses the channel. The mass ratio is an explicit variable: 100 and 400.
- **What the test-particle picture adds, separated from kinetic electron physics** (review, Major 5): test species
  (PIC-D9) of two kinds.
  - (i) A suprathermal shell as in `guo2010effect`, in the self-consistent fields.
  - (ii) Thermal and shell electrons pushed offline in the stored fields low-pass filtered to ion scales (k d_e ≪ 1).
  - Kinetic electrons minus (ii) is what electron-scale physics changes.
- **What the difficulty is:**
  - engine work: E1–E3 (Approach);
  - computing time, about a month of the node (Approach);
  - a minimal analytic model that gives the number to be tested. It is not written yet (open; review, Major 4).

## Related work and novelty
- **Search:** 2026-11-13, `survey-shock-alfvenic-turbulence.md`, plus the works the review asked for
  (`jokipii2007adiabatic`, `qin2018effects`, `xu2025electron`, `xu2020electron`, `zank2015diffusive`).
  - Forward citations of eight seed works: 318 works, 95 screened.
  - About 12 arXiv queries, about 200 hits.
  - Web searches; observations and parameters.
  - Not done: asking experts (the authors of `fulat2026pic`); theses and 2026 meeting abstracts.
- **Novelty claim, as narrow as it is:**
  - The mechanism is known (`jokipii2007adiabatic`). Its test-particle tests disagree (`guo2015acceleration` against
    `qin2018effects`).
  - No fully kinetic run of a nonrelativistic electron–ion shock in large-scale incompressible magnetic or Alfvénic
    fluctuations was found. The kinetic runs with pre-existing turbulence are compressive (`fulat2023kinetic`,
    `fulat2026pic`).
  - The claim would be the first self-consistent kinetic test of multiple connection as an injection channel: its
    excess over the local-obliquity null, as a function of Ξ and of the mass ratio.
- **Objections an expert would raise:**
  - L_c and L_y ≪ the solar-wind correlation length (about 10⁴ d_i), so the result holds for Ξ ≤ 27 only. Guo &
    Giacalone had Ξ ≈ 38.
  - m_i/m_e = 100 suppresses access (above). The step to 400 shows the trend only.
  - 2-D restricts cross-field transport (`jokipii1993perpendicular`). A slab of waves is not turbulence.
  - Turbulence alone accelerates electrons (`comisso2022ion`). Null N2 covers that.
  - Typical interplanetary shocks have M_A ≲ 3 (`trotta2025overview`). M_A = 4.5 is above the ripple threshold
    (3.5, `trotta2019electron`), so ripples are part of the baseline.

## Approach
Proposed with the question. The numbers are fixed by the pilot runs, and the criteria are fixed before the
production runs.

**Geometry:**
- x is the shock normal, with a reflecting wall at x = 0 (piston).
- y lies along the front and is periodic. B0 = B0 ŷ, so θ_Bn = 90°.
- Slab fluctuations δB_x(y), δv_x(y): linear Alfvén waves along B0, spectrum as in `guo2010effect`, L_c = L_y / 4.
  Circular polarization is a check.
- The laminar null runs tilt B0 in the plane: B0 = B0 (cos θ, sin θ, 0).

**Engine prerequisites.** Each one is done on main: a branch from main, a pull request, a review by Codex (sol 6.1
ultra) repeated until it reports no findings, the gates of `pypic/`, then the merge, and `research-sync-main` into
trunk and this topic. Runs that need a change start after its merge (owner, 2026-11-13).

| ID | Change | Needed for |
|---|---|---|
| E1 | Non-uniform initial E tangential to a reflecting face (extends PIC-D18; today `initial_fields._check_tangential_E` refuses it) | δB_x(y) loaded over the whole upstream: E_z(y) = (v0 − δv_x) B0; pilot |
| E2 | Turbulent inflow face: time-dependent δB and E = −v × B at the injection face, and δv(y, t) on the injected particles (extends PIC-D13, a boundary drive like PIC-D20) | waves that meet the shock at a fixed age; needed because a linear polarization rotates within 7 Ω_ci⁻¹ at L_y = 50 d_i (`design/estimates.md`, Sect. 1) |
| E3 | Restart checkpoints (roadmap R2) | runs of 20–160 h |
| Analysis | A_z labelling, encounter counting, offline pusher in stored filtered fields (topic code; output cadence to be planned) | Idea (a)–(c) and test electrons (ii) |

- E4 of the first draft (a relaxed turbulence slab) is not needed for this question and is dropped.
- Alternative to E2, to be weighed in the plan of the engine work: a reflecting wall moving through a periodic
  plasma (`behar2024impact`).

**Resolution (from `design/estimates.md`, Sect. 5).**
- At ω_pe/Ω_ce = 10, β_e = 1 and 5 cells per d_e, λ_De / dx = 0.35. The explicit scheme may then heat the grid.
- The pilot decides between three options:
  - ω_pe/Ω_ce = 4, which gives λ_De / dx ≈ 0.9;
  - more cells;
  - the energy-conserving implicit scheme (M5), which is slower per step.

**Runs** (m_i/m_e = 100 unless stated; M_A = 4.5; β_e = β_i = 1; δB_x/B0 = 0.6; 40 Ω_ci⁻¹; cost at 5 cells per d_e,
20 particles per cell, ω_pe/Ω_ce = 10):

| Set | Runs | Cost |
|---|---|---|
| Laminar null | θ_Bn = 90, 80, 70, 60, 45° × 2 seeds, 200 × 50 d_i | 10 × 5 h |
| Ξ scan | L_y = 50, 100, 200, 400 d_i × 4 realizations, L_x = 200 d_i | 4 × (5 + 10 + 20 + 40) h |
| N2, no shock | the same fluctuations in a periodic box, one per L_y | about 20 h |
| Mass ratio | m_i/m_e = 400 at L_y = 100 d_i, turbulent and laminar 90°, 1 each | 2 × 160 h |
| Total | | about 700 h, about 29 days of the node |

**Controls and null tests** (to be fixed before the runs):
- N1: the laminar runs (the local-θ_Bn null).
- N2: no shock.
- N3: test electrons (i) and (ii) of the Idea.
- N4: 4 realizations per turbulent case.
- N5: twice the particles per cell, and twice the resolution, on one case.
- N6: upstream δB_x/B0, δn/n, |B| variance, polarization angle and spectrum at arrival, vs time. A tolerance for
  "equal upstream" is to be fixed.

**Definitions** (review, Major 7; to be fixed before the runs):
- Thermal fits in patches of a few d_i, in the local bulk frame.
- A physical injection threshold: a gyroradius of the order of the thermal ion gyroradius, or the SSDA → DSA
  threshold of `amano2022theory`.
- The downstream region measured from the local, corrugated front.
- A second definition, to show that the conclusion does not depend on it.

**The claim would be shown wrong if:**
- the excess of the turbulent runs over the local-θ_Bn null does not exceed the spread of the 4 realizations at the
  largest Ξ;
- or the excess does not grow with Ξ;
- or electrons with one front encounter carry the same gain as those with several.
- The size of the effect expected from the minimal model (open) sets whether 4 realizations can see it.

**3-D:** reframed as a check of one 2-D result at matched Ξ (review, Major 9). Its box in d_i, L_c and Ξ is open.

## Claims and evidence
None yet. No run has been made.

## Publication
Open. Intended: a peer-reviewed journal paper (ApJ or JGR).

## Status
2026-11-13:
- **Done:**
  - topic chosen by the owner (theme D);
  - branch cut;
  - survey with 56 literature notes;
  - first draft;
  - review (`reviews/2026-11-13-plan.md`);
  - design estimates (`design/estimates.py` and `estimates.md`);
  - notes on the works the review asked for;
  - this revision, which proposes one question and a new Approach.
- **The owner's decision on the question is open.**
- **Missing before the question is fixed:**
  - the minimal model and the expected size of the effect (review, Major 4);
  - an expert check (the authors of `fulat2026pic`);
  - full texts of `jokipii1993perpendicular`, `giacalone2005*`, `jokipii2007adiabatic` and `barnes1974large`, not
    reachable from here (IOP, Wiley and ADS refuse).
- **Then:** the plan of E1 on main.
