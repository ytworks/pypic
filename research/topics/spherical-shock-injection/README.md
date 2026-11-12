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
  acceleration by the curvature of the bow shock (correction from the review of 2026-11-12). How far the
  literature treats curvature is for the survey to establish.
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
Open: the survey has not been done. The works read so far (notes in `research/literature/`):
- **hoshino2002nonthermal:** electron surfing at high-Mach shocks. Threshold M_A > 43–75 at the
  real mass ratio (Eq. 8). 1-D, m_i/m_e = 20.
- **amano2010critical:** critical Mach number M_A/√β_e ≳ 30 for injection (Eq. 3). Linear
  theory, compared with spacecraft statistics.
- **amano2022theory:** injection at oblique shocks of finite thickness, M_A/cos θ_Bn ≳ 55 (Eq. 14).
  A small mass ratio makes injection artificially easy (Eq. 21, Sect. 4.4).
- **hoshino2012stochastic:** stochastic (second-order Fermi type) acceleration in multiple islands,
  for reconnection.

Further candidates are listed in `literature.md`. They are not yet read and are not cited here
until they are.

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
Open: no systematic survey yet. Planned searches:
- citations of the notes above, in both directions;
- keywords: curved or spherical shock electron acceleration, blast wave PIC electron injection,
  θ_Bn dependence of electron acceleration, shock rippling and electron acceleration;
- 2015–2026 of ApJ, PRL, JGR Space Physics and Phys. Plasmas;
- later papers of the same groups.

**Candidate novelty (unverified):** attributing electron injection to its mechanisms over the whole
range of θ_Bn on one curved 3-D full-PIC shock, and finding the channels that exist only on a
curved front.

**Strongest expected objection:** a reduced mass ratio and the Mach numbers reachable in 3-D could
make injection artificially easy (amano2022theory Sect. 4.4) or stay below the surfing threshold
(hoshino2002nonthermal Eq. 8). The front could also be too small against ion scales for the planar
mechanisms to develop.

## Approach
Draft; the criteria and null tests are to be fixed before any production run.
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
None yet.

## Publication
Open.

## Status
2026-11-12: the purpose is stated and the branch exists.

Missing:
- the survey (Prior work, Related work and novelty);
- a narrowed research question;
- a confirmed idea;
- the driver and parameters, and the criteria and null tests.

Reviewed 2026-11-12 (`reviews/2026-11-12-purpose.md`): the purpose is a programme, not yet a question.
The reviewer proposes to narrow it to the θ_Bn at which injection switches on, and to its scaling with
m_i/m_e and β_e. That decision is the owner's and is open.
Next: the owner's decision on the question; the survey, starting with the curved-shock works the review names.
