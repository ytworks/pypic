# Survey: electron acceleration at nonrelativistic shocks in pre-existing Alfvénic turbulence

Date: 2026-11-13. Notes on every work cited here: `research/literature/<citekey>.md`; map: `literature.md`.
The searches were done by agents with web access, and each note states its read level (abstract, skimmed, full).

## 1. Scope
- **In scope:** what is known about electron injection and acceleration at a nonrelativistic collisionless shock that
  moves into pre-existing magnetic fluctuations with δB/B0 of order 0.1–1, and about how such fluctuations behave in
  kinetic simulations.
- **Out of scope:** turbulence that the shock generates itself (foreshock waves, ripples), except as the baseline.

## 2. Search routes
1. **Backward citations** of the direct predecessors (`fulat2023kinetic`, `fulat2026pic`, `guo2010effect`,
   `guo2015acceleration`), read from their reference lists.
2. **Forward citations** through OpenAlex (`filter=cites:`). Semantic Scholar returned HTTP 429.
   - Seeds: Guo & Giacalone 2010 and 2015, Giacalone 2005, Nakanotani 2022, Trotta 2021 and 2023, Fulat 2023 and 2026.
   - 318 citing works. A title filter (shock with turbulence, kinetic, electron or similar) left 95, which were screened
     by title and abstract.
3. **arXiv API, about 12 queries,** about 200 hits screened:
   - "pre-existing turbulence" shock;
   - shock + turbulence + PIC + electrons;
   - "field line" meandering shock;
   - "Alfvenic turbulence" shock;
   - "turbulent upstream"; "upstream fluctuations";
   - "shock-turbulence interaction" kinetic;
   - "parametric decay instability Alfvén wave" particle-in-cell;
   - nonrelativistic electron–ion PIC turbulence electron heating;
   - turbulent inflow boundary PIC/hybrid shock.
4. **Web searches:**
   - the Pohl group's 2025–2026 output: found ICRC 2025 (`fulat2025kinetic`) and an HPC project summary, which is not a
     publication;
   - observations named in the brief (Lindberg 2024, Jebaraj 2023 and 2024, Trotta 2025, the 2026 Nature
     Communications paper);
   - solar-wind parameters (correlation length, β, δB/B).
5. **Earlier scoping searches of the same day,** for the choice of topic: four candidate themes, two rounds each.
   Their results are in `notebook.md` (2026-11-13).

**Corrections found during the search:**
- arXiv:2404.07404 is not Lindberg et al. 2024 but Amano et al. 2024 (Phys. Plasmas); no arXiv id was found for
  `lindberg2024electron`.
- Giacalone 2005 (ApJ 624, 765) uses test particles in prescribed fields, not a hybrid code.
- OpenAlex gives `fulat2026pic` the DOI 10.3847/1538-4357/ae71b6 (ApJ 1005, 107); the note, copied verbatim from
  another topic, says it has none.

## 3. What the search found on the question
- **Fully kinetic, nonrelativistic, pre-existing magnetic or Alfvénic turbulence: none found.**
  - The fully kinetic nonrelativistic works with pre-existing turbulence are the Fulat line only: `fulat2023kinetic`,
    `fulat2025kinetic`, `fulat2026pic`, and the ICRC 2023 proceedings (PoS 444, 286), which repeat the 2023 paper.
  - Their turbulence is compressive at kinetic scales: tiles of λ_si/4, 50 velocity modes up to 6 × 2π/L_y.
  - Amplitude: δn/n ≤ 15 %, δB/B ≈ 10 %, capped by heating.
  - Their result: no effect on electrons at 90° (B0 in the plane), and 0.14 → 0.23 % non-thermal at 60° (B0 out of the
    plane, Ly ≈ 12 λ_si, about 20 Ω_ci⁻¹).
- **Fully kinetic with magnetic turbulence: relativistic pair plasma only.**
  - `bresci2023particle`: l_c ≈ 200 c/ω_p, tails s ≈ 2.5–3.5. Controls without turbulence and without a shock show no
    tail. The work done by the mean motional field and by δE is split per particle.
  - `demidem2023relativistic`: density waves.
- **The meandering picture rests on test particles.**
  - `guo2010effect`: 2-D hybrid fields, δB²/B² = 0.3, L_c = 1024 d_i, electrons from 100 eV; 200–300 × in 50 Ω_ci⁻¹
    at θ_Bn = 90°.
  - `guo2012particle`: about 10 % of the electrons above 15 keV at M_A ≈ 2.
  - `guo2015acceleration`: kinematic 3-D fields; at σ² = 10, about equal acceleration at every θ_Bn.
  - `guo2021shock` says the effect of turbulence "has not been understood in a quantitative way".
  - `lo2026effect`, the newest work of the group, is still hybrid and test-particle.
  - The shares of meandering and of ripples are not separated in `guo2010effect`. `burgess2006simulations` gets only
    2–3 × without ripples.
- **Hybrid with ions only:**
  - `trotta2021phase` and `trotta2024study` (2-D, θ_Bn 45°, B0 in the plane, δB/B0 up to 2.1): the shock goes from
    laminar to turbulent, and v_sh rises from 1.4 to 2.0 v_A.
  - `trotta2023three` (3-D, 128³ d_i³, δB/B0 0.5 and 1, M_A ≈ 6, 20 Ω_ci⁻¹): corrugation at the turbulent scales;
    ripples survive at 0.5.
  - `nakanotani2022turbulence`: no mean field.
- **2-D geometry** (`jokipii1993perpendicular`): with an ignorable coordinate, a particle stays within a gyroradius
  of one contour of A along the ignorable axis. A 2-D run follows meandering along in-plane field lines, but not the
  transport that needs field lines to separate.
- **Observations that tie upstream structure to electrons:**
  - `liu2026efficient` (MMS): a rotational discontinuity bends field lines into a V; electrons between two converging
    field-line–shock intersections go from about 17 eV to above 200 k_B T_e.
  - `raptis2025revealing`: an injection threshold of 1–5 keV in a foreshock transient.
  - `jebaraj2023relativistic`: 30 keV–1 MeV beams ahead of an M_A ≈ 3.3 shock.
  - `jebaraj2024acceleration`: about 6 MeV at a near-Sun shock; its δB/B ≥ 0.5 was generated by the shock.

## 4. Parameters of heliospheric shocks and the solar wind
| Parameter | 1 au | Near the Sun | Source |
|---|---|---|---|
| M_A, interplanetary shocks | mostly ≲ 3 | one strong case 9.1 at 0.24 au | `trotta2025overview`, `jebaraj2024acceleration` |
| M_A, bow shock | 5–11 (7 crossings) | – | `lindberg2024electron` |
| θ_Bn | quasi-perpendicular more frequent | more quasi-parallel inside 0.8 au | `trotta2025overview` |
| β_e, β_p (medians) | 1.09, 1.05; in ICME flux ropes 0.40, 0.13 | β_i 0.16 (one event) | `wilson2018statistical`, `jebaraj2024acceleration` |
| δB/B0 | no statistic read | ∼ 1 in Alfvénic wind | `chen2020evolution` |
| Correlation length | 1.5–2 × 10⁶ km, ≈ 1.5–2 × 10⁴ d_i | ≈ 5 × 10⁴ km at 0.1 au | `cuesta2022isotropization` |
| ω_pe/Ω_ce | ≈ 150 (B = 4.3 nT, n = 4 cm⁻³) | ≈ 24 (90 nT, 46.5 cm⁻³) | arithmetic from `liu2026efficient`, `jebaraj2024acceleration` |
| Whistler critical Mach number 0.5 √(m_i/m_e) cos θ_Bn | 10.7 at 60°, 3.7 at 80° | same | form of `lindberg2024electron` |

A kinetic box cannot hold the correlation length at 1 au (about 10⁴ d_i). It holds a part of the inertial range, as
`trotta2023three` (128 d_i) and `comisso2022ion` (ℓ = 20 d_i) do.

## 5. Numbers for the design
- **Turbulence alone accelerates** (`comisso2022ion`, 3-D e-i PIC, δB/B0 = 1, v_A = 0.14c, m_i/m_e = 50,
  ℓ = 20 d_i):
  - about 2/3 of the magnetic energy goes to particles by about 8 ℓ/v_A, shared equally by the two species;
  - tails p_i ≈ 3.5, p_e ≈ 4.2 at β = 0.04.
  - A null run without a shock is therefore needed.
- **Decay of the turbulence** (`parashar2018dependence`, 2.5-D PIC, δb = δv = 0.5, box 150 d_i): relaxation over
  Ω_ci t < 40, then decay over more than 300 Ω_ci⁻¹, faster at higher β.
- **Parametric decay of a coherent |B| = const pump** (`gonzalez2023particle`; δB/B0 = 0.5, β = 0.25, k0 d_i = 0.245,
  m_i/m_e = 400, v_A/c = 0.05):
  - onset at about 150 Ω_ci⁻¹, saturation at about 225 Ω_ci⁻¹;
  - electrons get about 10 % of the pump energy.
  - The values for our β and amplitude need the linear theory (Derby 1978; Goldstein 1978; not read).
- **Set-up methods in the literature:**
  - turbulent slabs prepared in periodic boxes and matched at the upstream end (`fulat2023kinetic`);
  - a Langevin antenna (`tenbarge2014oscillating`, `bresci2023particle`);
  - continuous injection of a slab spectrum at the inflow boundary (`guo2010effect`);
  - a snapshot of an MHD run imposed upstream (`trotta2023three`);
  - turbulence decayed in a periodic box, then a moving obstacle in the plasma frame (`behar2024impact`).

## 6. Limits of this survey
- **Not reached:**
  - Matthaeus 2005; Bruno & Carbone 2013 (full text);
  - Pitňa 2021 (numbers);
  - Ha 2025 (review);
  - theses and AGU/EGU 2026 abstracts;
  - Derby 1978 and Goldstein 1978.
- **Several notes are abstract-only** (list in `literature.md`). `barnes1974large` rests on citing works only.
- **No expert was asked.** The stance lists that route.
