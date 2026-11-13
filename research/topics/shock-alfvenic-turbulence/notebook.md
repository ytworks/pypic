# Notebook — shock-alfvenic-turbulence

Dated lab notes, newest first. Write down what was tried, what was seen, and what was decided, failed attempts included.
Never rewrite a past entry. To correct one, add a new entry that refers to it.

## 2026-11-13 — Revision after the review: one question, design estimates

- **Paper work** (review, "Shortest path"): `design/estimates.py` (output `design/estimates.md`).
  - Assumptions are in its docstring. L_sh is taken as the foot width 0.68 M_A d_i (Woods / Livesey, not read).
  - **Polarization rotation (kinetic, parallel Maxwellian):** about 1.5 × the Hall-MHD rate. The x-polarized share
    of a linear slab falls to 0.75 after 6.9 Ω_ci⁻¹ at L_y = 50 d_i, after 87 at L_y = 200, and after 307 at L_y = 400.
  - **Ξ at δB_x/B0 = 0.6, M_A = 4.5:** 3.4 / 6.9 / 13.6 / 27 for L_y = 50 / 100 / 200 / 400 d_i. The set-up of
    `guo2010effect` gives 38.
  - **t_conn / t_along** does not depend on L_y: about 0.66 at m_i/m_e = 100 and about 2.8 at 1836.
  - **Thermal electrons able to stream upstream along the local field:** 7 % at m_i/m_e = 100, 18 % at 400, 29 % at
    1836 (δB_x/B0 = 0.6, M_A = 4.5).
  - **Grid heating:** λ_De / dx = 0.35 at ω_pe/Ω_ce = 10, β_e = 1, 5 cells per d_e. This is a new risk, to be settled
    by the pilot.
- **Notes the review asked for:**
  - `jokipii2007adiabatic` (abstract), `qin2018effects` (full), `xu2025electron` (abstract), `zank2015diffusive`
    (skimmed);
  - `xu2020electron`, copied verbatim from `research/spherical-shock-injection` and checked against a full read.
  - The full texts of the abstract-only design notes were not reachable (IOP bot check, Wiley 403, ADS 403/504).
- **Corrections:**
  - `barnes1974large`: the abstract, found through OpenAlex, describes a transverse mode with constant pressure,
    density and |B|. It names a linearly polarized finite-amplitude wave for a small range of directions only. It
    does not mention arc polarization, which the first README draft and the note cite.
  - `fulat2026pic`: the published reference is ApJ 1005(1), 107 (2026-06-29), DOI 10.3847/1538-4357/ae71b6 (Crossref).
    The note is kept verbatim (shared with another branch), so this entry is the record.
- **README revised:**
  - One proposed question: multiple connection against the local-θ_Bn null, scaled by Ξ and by the mass ratio.
  - The Approach is rebuilt on the review's controls.
  - E4 is dropped.
  - About 700 h of the node.
- **Not done:** the minimal model, the expert check, and the owner's decision on the question.

## 2026-11-13 — Review of the plan

- The `research-reviewer` agent reviewed README, survey, literature map and notebook. Its report is kept in
  `reviews/2026-11-13-plan.md`.
- **Main points:**
  - The δB_x/δB_z polarization control does not isolate meandering. δB_x also lowers the local θ_Bn, which switches
    on SDA, whistlers and SSDA; it changes the 2-D wave content, the front forcing and the downstream turbulence.
  - Hall dispersion rotates a linear polarization within about 8–30 Ω_ci⁻¹ at L_c = 10–20 d_i (reviewer's estimate).
  - There is no controlling number: Ξ = ξ_x / L_sh and v_te cos θ_loc / U1 are proposed.
  - There is no analytic model.
  - N3 cannot differ from the kinetic electrons.
  - The question is two questions; the matrix is about 288 runs.
  - The 3-D box is too narrow.
  - Survey gaps: Jokipii & Giacalone 2007, Qin et al. 2018 (opposite trend), Xu et al. 2025.
- **Not yet acted on.** The revision of the question and of the Approach follows.

## 2026-11-13 — Topic chosen and opened; engine prerequisites found

- **Request (owner):** a research theme on shock acceleration that uses this machine fully, found from a survey of prior
  work and new to the field; first cut a branch and write the plan.
  - Machine: 8 × RTX PRO 6000 (96 GB each), one node.
  - Throughput measured on `blast_wave/full`: 3.1e9 particles × 2800 steps in 788 s, i.e. about 1.1e10 particle-steps/s.
  - Particle capacity: about 3–4e9 particles.
- **Candidates surveyed** (agent searches of arXiv, ADS, journal pages):
  - A: a shock crossing an upstream current sheet.
  - B: ε_e(θ_Bn) at high Mach number in 2-D.
  - C: positron selection at magnetised e±p shocks.
  - D: a shock in large-amplitude Alfvénic turbulence.
  - The owner asked for the comparison of A and D as facts only. The facts for D are in the README (Prior work).
- **Facts that decided the closest gap for D:**
  - Fulat et al. (2023, 2026) are the only fully kinetic runs of shocks in pre-existing turbulence found. Their
    turbulence is compressive, δn/n ≤ 15 % and δB/B ≈ 10 %, with the amplitude capped by heating.
  - Both papers state: "Future studies will explore the effect of Alfvénic turbulence."
  - The field-line-meandering picture (Guo & Giacalone 2010, 2015) rests on test electrons in hybrid or kinematic fields.
- **Owner's decision (2026-11-13):** theme D. Research quality comes first, and the plan must not assume an unchanged
  engine.
- **Engine check** (source read; `pypic/initial_fields.py`, `pypic/simulation.py`, `pypic/config.py`):
  - Usable as they are:
    - velocity perturbations, any number of sinusoidal modes along a periodic axis, any component (PIC-D9);
    - B profiles along another axis, added together (`_add`);
    - a uniform `E0` tangential to reflecting faces (PIC-D18);
    - test species (`SpeciesConfig.test`).
  - Refused: an E profile whose component is tangential to a non-periodic face (`_check_tangential_E`).
    - With B0 = B0 ŷ in the plane and flow −v0 x̂, an Alfvénic mode with δB_x(y) needs E_z(y) = (v0 − δv_x) B0, and
      that profile is refused in a piston box (reflecting x faces).
    - The δB_x part is what carries field lines across the shock.
    - The Guo & Giacalone set-up (δB_x, δB_z varying along the shock surface, injected continuously at the inflow face)
      therefore needs engine work: E1 to E4 in the README (Approach).
  - **Procedure for the engine changes (owner, 2026-11-13):**
    - Each change E1–E4 goes on a branch cut from main and is opened as a pull request into main.
    - It is merged only when a review by Codex (model sol 6.1 ultra) reports no findings; fixes and re-reviews repeat
      until then.
    - The runs that depend on a change start after its merge and its sync into trunk and this topic.
    - The repository's own gate rules for `pypic/` still apply.
  - **Correction of a statement made during the comparison (chat, same day):** candidate A was described as buildable
    with an E_z(x) profile. The same check refuses that profile with reflecting x faces.

- **Opened:** branch `research/shock-alfvenic-turbulence`, cut from `origin/research/trunk` at 3d6b50e.
  - Literature notes shared with `research/spherical-shock-injection` were copied verbatim from that branch, so a later
    merge of both topics adds identical files: fulat2023kinetic, fulat2026pic, guo2010effect, guo2015acceleration,
    trotta2019electron, amano2020observational, katou2019theory, bohdan2017electron, matsumoto2017electron,
    riquelme2011electron, guo2014nonthermal.
