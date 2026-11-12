# Electron injection at a spherical collisionless shock

Branch: `research/spherical-shock-injection`. Started: 2026-11-12. Notebook: `notebook.md`. Survey: `literature.md`.
Sections follow `.claude/rules/research/stance.md`. Mark a section as open rather than guessing.

**Purpose (owner, 2026-11-12):** the electron injection problem. We use PIC simulations of a
spherical shock to find out which mechanisms accelerate electrons. Every mechanism is considered,
including unknown ones:
- known mechanisms: surfing, drift, wakefield, first- and second-order Fermi acceleration, and
  their combinations;
- mechanisms that come from the geometry of the curved front;
- changes in the shape of the front as the shock normal goes continuously from perpendicular to
  parallel to the field.

## Background
- **The injection problem.** Nonthermal electrons are inferred at astrophysical shocks such as
  supernova remnants and observed at heliospheric shocks. Diffusive (first-order Fermi) shock
  acceleration works only for particles whose gyroradius exceeds the shock thickness. Thermal
  electrons are far below that and must be pre-accelerated first; this is the electron injection
  problem (hoshino2002nonthermal Sect. 1; amano2010critical; amano2022theory Sect. 1).
- **Known injection mechanisms are framed for planar shocks with a fixed angle θ_Bn** between the
  shock normal and the upstream field:
  - surfing in electrostatic solitary waves of the Buneman instability at the foot
    (hoshino2002nonthermal);
  - a critical Mach number for whistler-mediated injection (amano2010critical);
  - shock drift acceleration and injection at oblique shocks of finite thickness (amano2022theory);
  - second-order Fermi acceleration in multiple magnetic islands, studied for reconnection
    (hoshino2012stochastic).

  This is not the whole picture: amano2022theory (Sect. 4.5, Eq. 22) already limits stochastic shock drift
  acceleration by the curvature of the bow shock (correction from the review of 2026-11-12). The survey of
  2026-11-12 found curved fronts treated for the Earth's bow shock only (`survey-snr-electron-injection.md`,
  Sect. 6.5); those works were read in part on 2026-11-13 (Prior work).
- **Real shocks are curved:** supernova blast waves, planetary bow shocks, shocks driven by coronal
  mass ejections, laser-driven shocks. On a spherical shock in a uniform field, θ_Bn runs
  continuously from 90° at the equator to 0° at the poles. The front is curved, and a blast wave
  slows down as it expands. The engine already drives such a shock in 3-D full PIC: a debris shell
  expanding into a magnetised plasma (`examples/blast_wave`, PIC-D17, gated at M_A = 3).

## Research question
Draft (the owner's purpose; to be narrowed after the survey and the review):

> At a spherical collisionless shock expanding into a uniform magnetic field:
> - which mechanisms inject thermal electrons into the nonthermal population — known ones alone
>   or combined, and any not yet described;
> - where on the front each one acts, as a function of the local θ_Bn, Mach number and curvature;
> - which channels exist only because the front is curved, so that a planar shock with the same
>   local parameters would not have them?

**Mechanisms in scope** (the owner's list, with our grouping):
1. Shock surfing acceleration (electrostatic solitary waves at the foot).
2. Shock drift acceleration, including its stochastic form.
3. Wakefield acceleration.
4. First-order Fermi (diffusive) and second-order Fermi (stochastic, in turbulence or islands) acceleration.
5. Combinations and sequences: one mechanism pre-accelerates electrons for the next.
6. Geometric mechanisms of a curved front:
   - drifts along the front;
   - transport along the surface from quasi-perpendicular to quasi-parallel regions;
   - converging or diverging field lines and mirroring;
   - the expansion and deceleration of the front.
7. The shape of the front as θ_Bn changes from perpendicular to parallel: quasi-perpendicular
   against quasi-parallel structure, rippling, departures from a sphere.
8. Unknown mechanisms: energy gain that none of the above accounts for.

## Prior work
Surveyed for supernova remnant shocks on 2026-11-12: `survey-snr-electron-injection.md` (54 works read, with
references; the map by line of work is in `literature.md`). What it establishes for this topic:
- **Mechanisms are known by obliquity**, each from planar simulations at one or a few angles: surfing, shock
  drift with scattering, whistlers and Weibel turbulence at quasi-perpendicular shocks; an electron foreshock
  at oblique shocks; reflection and proton-driven waves at quasi-parallel shocks (survey Sect. 4).
- **No study covers the range of θ_Bn with one set of parameters** (survey Sect. 6.1).
- **Models of whole remnants assume the efficiency law ε(θ_Bn)**; SN 1006 and the radio population of remnants
  point in opposite directions (survey Sect. 3.3).
- **Curved fronts:** treated for the bow shock only. These are the direct predecessors for curvature
  (`literature.md`, "Curved fronts", read 2026-11-13):
  - read in full: `savoini2001twodimensional` and `savoini2010origin` (2.5-D full PIC of a cylindrical piston,
    m_i/m_e = 42, M_A = 3), `vandas1994acceleration`, `vandas2001shockb`, `kong2016electron`;
  - known by abstract only, full texts not obtained: `kraussvarban1991electron` (the first treatment of
    electrons at a curved shock) and `giacalone2017acceleration` (protons at a spherical blast wave), the two
    that matter most for novelty, and the other works listed in `literature.md` under "Still to read".

The works read before the survey:
- **hoshino2002nonthermal:** electron surfing at high-Mach shocks. Threshold M_A > 43–75 at the
  real mass ratio (Eq. 8). 1-D, m_i/m_e = 20.
- **amano2010critical:** critical Mach number M_A/√β_e ≳ 30 for injection (Eq. 3). Linear
  theory, compared with spacecraft statistics.
- **amano2022theory:** injection at oblique shocks of finite thickness, M_A/cos θ_Bn ≳ 55 (Eq. 14).
  A small mass ratio makes injection artificially easy (Eq. 21, Sect. 4.4).
- **hoshino2012stochastic:** stochastic (second-order Fermi type) acceleration in multiple islands,
  for reconnection.

Works still to read are listed in `literature.md`.

## Idea
Open; one candidate, a hypothesis:
- A single spherical shock samples every θ_Bn and curvature at once, and electrons can move
  along the front between regions. A combination of mechanisms can then act along one trajectory,
  for example pre-acceleration in a quasi-perpendicular foot followed by further gain where θ_Bn
  is oblique. Planar boxes at a fixed θ_Bn exclude this by construction.
- What would make the question answerable is attributing the gain of each injected electron to
  its channels. This is done from tracked trajectories, by splitting the work of the field into
  parallel and drift terms, in coordinates of the front (radius, θ_Bn).
- The main difficulty may be computing time rather than an idea: a 3-D front large against ion
  scales, at a Mach number high enough for the planar thresholds above, with a reduced mass ratio.
  This is to be assessed.

## Related work and novelty
A first systematic search was done on 2026-11-12 (arXiv, by keyword, author and reference lists; log in the
survey, Sect. 9). Still missing:
- the part of the bow-shock literature on curved shocks whose full texts were not obtained (J. Geophys. Res.;
  list in `literature.md`, "Still to read");
- forward citations of the key works through a citation index;
- the heliospheric literature on electrons at interplanetary and coronal shocks.

The gaps that the survey finds relative to this topic are in its Sect. 8.3.

**Candidate novelty (unverified):** attributing electron injection to its mechanisms over the whole
range of θ_Bn on one curved 3-D full-PIC shock, and finding the channels that exist only on a
curved front.

**Strongest expected objection:** a reduced mass ratio and the Mach numbers reachable in 3-D could
make injection artificially easy (amano2022theory Sect. 4.4) or stay below the surfing threshold
(hoshino2002nonthermal Eq. 8). The front could also be too small against ion scales for the planar
mechanisms to develop.

## Approach
Draft; the criteria and null tests are to be fixed before any production run.

**Exploratory runs (designed 2026-11-12, run 2026-11-13):** three trials to build hypotheses, recorded in
`exploratory-conditions.md` with their configs in `configs/`.
- They use the debris shell of `examples/blast_wave` in a 768³ box with 2 particles per cell, and differ in the
  field strength only: nominal M_A = 3.5, 6.0 and 8.7 at m_i/m_e = 16.
- They are not evidence for a claim. Their criteria, floors and null tests were fixed before the runs in
  `exploratory-protocol.md`; the results and the verdict are in `exploratory-results.md`.
- Their main limit, expected and now measured: the front is slower than the debris. Where the protocol finds
  a shock, its fast Mach number is 1.5–2.5 in L and in M, and H forms none. The debris cloud is small against
  the debris gyroradius ρ_d = M_A d_i: R_m/ρ_d = 4.6 (L), 1.9 (M), 1.0 (H).
- Their second limit, measured: numerical scattering of thermal electrons at 2 particles per cell is above the
  registered limit in every reading (`exploratory-results.md`, Sect. 3), so no loss of magnetic moment in these
  runs can be assigned to physics.

The plan for the later, claim-bearing runs:
- **Driver:** a spherical shock driven outwards into a uniform field. The existing debris shell
  (`examples/blast_wave`, M_A = 3) is one option, and is likely too low in Mach number for
  surfing. The driver and the parameters (M_A, m_i/m_e, β_e, ω_pe/Ω_ce, radius in ion scales) are
  to be chosen against the thresholds in the Prior work section.
- **Diagnostics:**
  - energy spectra of electrons in bands of θ_Bn and radius;
  - tracked electrons (particle ids, track sets) with the fields along their orbits;
  - the work split into parallel-field and drift terms, and each injected electron classified by
    its dominant channel;
  - the gain that no channel explains, kept as a quantity of its own.
- **Controls:**
  - planar shocks at matched local θ_Bn and Mach number, to separate what curvature adds;
  - test particles in frozen fields;
  - a second seed;
  - scans of particles per cell, resolution and front radius.

## Claims and evidence
No claim.

Observations of the exploratory runs (2026-11-13). They are not claims: one seed, 2 particles per cell,
m_i/m_e = 16. Limits and detail: `exploratory-results.md`.
- No ambient electron population at γ ≥ 1.5 upstream of a front in L, M or H (at most 3 electrons per dump in
  L and 5 in M, against the 100 required); the null runs hold none.
  Records: `runs/2026-11-13-91f79ff1d-L.json`, `-ac33a41eb-M.json`, `-3f5c43715-H.json`, `-0114ef9ea-L-null.json`,
  `-3f5c43715-M-null.json`.
- The valid shocks of L and M are a decaying blast wave of fast Mach 1.5–2.5 and compression about 2; the
  debris stalls at 1.3–1.5 R_m (records of L and M, `front`).
- The adiabatic mirror in a static potential does not describe the electrons that return upstream, in number,
  energy or magnetic moment (records of L and M, `baseline`, `tracks`). This is a statement about the baseline
  and the numerical floor, not about the plasma (review R2): in the valid domain 98 % or more of the predicted
  beam comes from bands where the baseline is not defined (ε > 0.3), and the share of returns that keep the
  magnetic moment equals the share in the null runs (records of L, M, L-null, M-null; blocks `r2`, `r2null`).

## Publication
Open.

## Status
2026-11-12: the purpose is stated and the branch exists.

Missing:
- the curved-shock part of the survey (Prior work, Related work and novelty);
- a narrowed research question;
- a confirmed idea;
- the driver and parameters, and the criteria and null tests.

Reviewed 2026-11-12 (`reviews/2026-11-12-purpose.md`): the purpose is a programme, not yet a question.
The reviewer proposes to narrow it to the θ_Bn at which injection switches on, and to its scaling with
m_i/m_e and β_e. That decision is the owner's and is open.

Surveyed 2026-11-12: electron injection at supernova remnant shocks (`survey-snr-electron-injection.md`).

Designed 2026-11-12: three exploratory runs (`exploratory-conditions.md`).

Run and analysed 2026-11-13: L-null, L, M, M-null, H, against the protocol fixed before them
(`exploratory-protocol.md`; reviews in `reviews/`). **Verdict: no hypothesis was obtained from these runs**
(`exploratory-results.md`): no signal at the registered threshold, and the runs sampled a weak, decaying shock
instead of the designed Mach numbers. Review R2 (`reviews/2026-11-13-results-r2.md`) confirms the verdict; its
corrections are in the results document and in the notebook.

Next (order of `exploratory-results.md`, Sect. 8):
- reading `kraussvarban1991electron` and `giacalone2017acceleration` in full (access beyond arXiv);
- the owner's decision on the question, with a controlling number and a threshold;
- one analysis of the existing data without GPU time: a count model of thermal escape per band, to close the
  observation "returns beyond the mirror baseline" as numerical floor plus known leakage, or to keep it open;
- only then a design. A 2.5-D design comes first for the baseline and the numerical floor; a 3-D run is for
  what only 3-D can show. A test-particle baseline in the stored slices is not planned as it stands (the slices
  are one gyro-period apart).
