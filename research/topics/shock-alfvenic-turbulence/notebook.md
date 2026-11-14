# Notebook — shock-alfvenic-turbulence

Dated lab notes, newest first. Write down what was tried, what was seen, and what was decided, failed attempts included.
Never rewrite a past entry. To correct one, add a new entry that refers to it.

## 2026-11-14 — Set-up of the runs; purpose corrected; pilot P1a fixed

- **Purpose corrected by the owner.** While writing the set-up I treated the differences from `guo2010effect`
  (wave mix, variance, L_c) as things to match. The owner: the purpose is not to replicate prior work but to do the
  next problem that the prior kinetic work left.
  - My reading: the future work named in `fulat2023kinetic` and `fulat2026pic`. Stated to the owner twice; not yet
    confirmed or rejected.
  - Consequence: `guo2010effect` is the source of the mechanism only. `design/setup.md` records the differences as
    facts (Sect. 2) and sets the set-up against the statements of Fulat et al. (Sect. 1).
- **Corrections of earlier text:**
  - Entry "Topic chosen and opened" (2026-11-13) and the README said both Fulat papers state "Future studies will
    explore the effect of Alfvénic turbulence." That is the 2026 wording. 2023: "Further studies should explore the
    effect of the Alfvénic turbulence, where the magnetic field fluctuations dominate." (arXiv texts, read
    2026-11-14.) 2023 also names larger scales with a wider box, and lower M_A.
  - README, E2 row: "waves that meet the shock at a fixed age" was wrong reasoning. The modes have k_x = 0, so the
    upstream state depends on (y, t) only and has no age along x. The face supplies the wave for the plasma that
    enters after t = 0.
- **`guo2010effect` prescription verified** against the TeX source (arXiv:1003.5946): equal mixture of R and L,
  forward and backward waves; P(k) ∝ 1 / (1 + (k L_c)^{5/3}), L_c = L_z; σ = δB²/B_1² = 0.3; δv = v_A δB/B_1;
  loaded at t = 0 and injected at the inflow boundary; piston with a conducting wall.
- **Upstream state chosen: one forward circular R wave**, rms δB_x/B0 = 0.6 (`design/setup-tables.md`, Table 1).
  - Random phases: std / mean of |B|² = 0.21–0.37 at t = 0. The question says incompressible.
  - Several modes with uniform |B| at t = 0 (Blaschke product): the spread reaches 0.142 at t = 40 for L_y = 50 d_i,
    0.036 for 100, 0.010 for 200.
  - One wave: 0 at all times (linear phases); exact cold two-fluid solution.
  - Not read: any source on the parametric or modulational stability of the finite-amplitude wave at β = 2. The
    pilot measures it; R and L are both run.
- **Numbers that changed with the state** (new facts for the README amendments):
  - Ξ = 4.4 / 8.8 / 17.7 / 35.3 at L_y = 50 / 100 / 200 / 400 d_i (was 3.4 / 6.9 / 13.6 / 27, random phases).
  - Access ratio 0.85 at m_i/m_e = 100 (was 0.66).
  - rms δB_x/B0 = 0.15 at L_y = 50 d_i gives Ξ = 1.10, not 0.85: the Ξ < 1 set needs a lower amplitude (0.1 gives
    0.74), fixed after L_sh is measured.
  - With one wave a "realization" is a shift along y: N4 becomes 4 particle seeds.
- **How the wave meets the front** (`design/setup.md`, Sect. 5): the front holds all phases at once; local θ_Bn from
  49.7° to 130.3°; local normal inflow 5.30 to 3.70 v_A; the pattern slides along the front at v_ph ≈ v_A; a field
  line stays doubly connected for ξ_pp / U1 = 3.0–24 Ω_ci⁻¹.
  - The laminar null must therefore sample the full circle of phases, with B_z and the bulk velocity of the wave.
  - That needs a uniform B normal to the reflecting wall. `B0` cannot hold it (`docs/design/engine.md`, PIC-D12:
    the wall zeroes the normal B of the Yee fields). External fields (PIC-D12) may be normal to a face; untested for
    a piston shock. Called E4 from here on (not the dropped E4 of the first draft).
- **Box and face:** nz = 2 (engine minimum), 2D3V, ppc 5 per species. For L_x = 200 d_i the plasma that entered
  through the face reaches the shock at 44.4 Ω_ci⁻¹, a fast-mode signal from the face at 32.6 (uniform-state speeds).
- **Smoke runs** (scratch, 8 GPUs, `configs/smoke.yaml`, on 458bfd9 with the untracked topic code): Gauss max
  7.3e-13, div B max 5.3e-14, dE/E0 = +3.49e-3, 13.89 ms per step.
  - Failed first: the periodic branch of `plan.py` passed `injection=None`; the engine needs `Injection()`. Fixed.
- **Grid heating measured** (scratch, laminar periodic 4 × 10 d_i, 1 Ω_ci⁻¹, 5 cells per d_e, ppc 5):
  - ω_pe/Ω_ce = 10 (λ_De/dx = 0.35): kinetic energy of the electrons +4.85 % per Ω_ci⁻¹ (first half +2.1 %); ions
    −1.1 %; budget drift +7.5e-3.
  - ω_pe/Ω_ce = 4 (λ_De/dx = 0.88): electrons −0.41 % (first half −0.34 %); ions −0.18 %; budget drift +3.4e-4.
  - With `analyze.py` (thermal energy in the local bulk frame, 1 d_i blocks): +4.88 % and −0.31 %.
  - ω_pe/Ω_ce = 4 is the candidate; one Ω_ci⁻¹ does not show the rate over 40.
- **`analyze.py pilot`** written and tried on the scratch outputs.
  - Bug found: `gauss` and `divB` scalars hold NaN between samples; `max` gave NaN. Now `nanmax`.
- **Pilot fixed before the runs** (`design/setup.md`, Sect. 8): P1a = eight periodic boxes 10 × 50 d_i
  (`configs/p1a-*.yaml`); criteria C1 (laminar |ΔT_e/T_e| ≤ 0.05 over 40 Ω_ci⁻¹), C2 (one wave: δn_e/n and |B|²
  spread ≤ 0.05), C3 (rms δB_x/B0 ≥ 0.54 at t = 40), several modes within a factor 2 of the linear spread. Then
  P1b, P2, P3.
- **Review:** the owner asked for a cross-review of the set-up by Codex (sol 6.1 ultra) on a pull request, on
  whether the way the shock meets the large-amplitude fluctuations is what the prior work intended.

## 2026-11-14 — Engine prerequisites merged

- PR #50 merged into main (fa0972664): E1 (PIC-D24), E2 (PIC-D26), E3 (PIC-D25, OUT-D6). The Codex review (sol 6.1
  ultra) reported no findings in its 13th round; gates all passed on 9666c49 in fp64, mixed and compact.
- Taken: main → `research/trunk` (abc1099, baseline recorded in 244a05f) → this branch (c70c322).
- What the example `examples/alfven_inflow` measured on ef4fc4f (records
  `bench/results/2026-11-14-ef4fc4f08-example-alfven_inflow-piston-run-*.json`):
  - piston shock speed 0.0320 c, against 0.0336 c from the jump conditions; B_z behind the shock 4.25 B0;
  - a restart continued to step 8200 agreed with the run without restart to 1.2e-16;
  - TF/SF face: about 2 ms per step of overhead, and a field error of 5–15 % next to the face; its cause is not
    established.
- Differences from the plan of the README (first draft of E1): the motional E is loaded by a charge displacement,
  not as an E profile, and only away from non-periodic faces: the pilot loads the wave under an envelope h(x).
- Wiki page of the example published (Alfven-Inflow, EN and JA).

## 2026-11-13 — Expert check skipped

- The owner decided that the expert check (asking the authors of `fulat2026pic`) is skipped, because it cannot be
  done.
- The novelty claim rests on the literature search alone (`survey-shock-alfvenic-turbulence.md` and the notes the
  review asked for). The README says so under "Related work and novelty".

## 2026-11-13 — Minimal model; amendment of the falsification criteria

- **`design/model.py`** (output `design/model.md`): an encounter chain on a multiply connected field line.
  - Adiabatic reflection in the HTF of each intersection, with the loss cone 1/r.
  - Arrival speed raised by 2 U1/cos θ at the next intersection.
  - Energies counted in the downstream frame, with the magnetic moment conserved.
  - Neglected: the cross-shock potential, scattering, ripples, and the time of each encounter.
- **First version wrong:** it counted HTF energies. Those contain the frame speed U1/cos θ, so 70 % of the null was
  above 10 T_e. It was corrected before any use.
- **Results** (M_A = 4.5, β = 1, r = 2.78):
  - At δB_x/B0 = 0.6, m_i/m_e = 100: reflected at the first encounter 0.355; N = 2 0.325; above 10 T_e 1.0e-2 (null)
    against 7.7e-2 (multi); above 30 T_e 1e-5 against 1.5e-3.
  - The chain ends after a few encounters (loss cone). No power law without scattering.
- **New fact for the plan:** in the periodic slab, every field line that touches the front returns to it
  (A_z = −B0 x + a(y), a periodic).
  - The model then has no Ξ dependence above the threshold Ξ ≳ 1, where the arc leaves the shock layer.
  - So the accepted criterion "the excess does not grow with Ξ → claim wrong" would reject a correct claim.
- **Amendment of the Approach (README):**
  - criteria: excess at Ξ > 1, none at Ξ < 1, gain carried by multiply-encountering electrons;
  - model numbers to compare;
  - a set of Ξ < 1 runs (δB_x/B0 = 0.15, L_y = 50 d_i, Ξ ≈ 0.85; 4 × 5 h);
  - limitation: the runs test what multiple connection does, not how often field lines return in nature.
  - Reported to the owner.

## 2026-11-13 — Question fixed by the owner

- The owner accepted the question and the Approach of the revision (b793c45): "これでok".
- Still to come before the production runs:
  - the minimal model and the expected size of the effect;
  - the expert check;
  - the engine work E1–E3 on main;
  - the pilot runs.

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
