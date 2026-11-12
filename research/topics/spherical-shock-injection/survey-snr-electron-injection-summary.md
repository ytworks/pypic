# Electron injection at supernova remnant shocks: the picture on one page

Summary of `survey-snr-electron-injection.md` (2026-11-12). Section numbers refer to that survey, and citekeys to
the notes in `research/literature/`. Japanese version: `survey-snr-electron-injection-summary.ja.md`.

## The problem
- Diffusive shock acceleration (DSA) works only for particles that can cross the shock. At a real supernova
  remnant (SNR) shock, thermal electrons are short of that by a factor 10⁵–10⁶ in energy
  (`bohdan2023electron` Eq. 1). Something must pre-accelerate them. This is the injection problem (Sect. 2).
- Even the threshold is disputed: a gyroradius, a momentum, or simply a speed above the shock speed
  (`gupta2025speed`; Sect. 2).

## What observations fix
- Young remnants hold 10–100 TeV electrons at their shocks (`vink2012supernova` Sect. 6.3.1).
- Accelerated electrons are about 10⁻³ of the protons: K_ep = 1.6 × 10⁻³ in Tycho (`morlino2012strong` Sect. 3)
  and 9.3 × 10⁻⁴ in SN 1006 (`mcginness2026sn` Table 1).
- Where on the shell electrons are accelerated is disputed. SN 1006 points to the quasi-parallel caps
  (`bocchino2011constraints`; `winner2020evolution`). The radio population of bilateral remnants points to
  quasi-perpendicular shocks (`west2017connection`; Sect. 3.3).
- No observation measures the injected fraction or identifies the mechanism.

## The mechanisms, by the angle θ_Bn between shock normal and field (Sect. 4)

| θ_Bn | First stage | Second stage | How far simulations got |
|---|---|---|---|
| Quasi-perpendicular, 75–90° | surfing on Buneman waves, or acceleration by whistler waves | shock drift with scattering (SSDA); reconnection and stochastic scattering in Weibel turbulence | about 0.5 % nonthermal electrons, reaching at most 10–100 kT, in 2-D at 90° (`bohdan2023electron` Sect. 3.3); a power law of index 3.5 in the one 3-D run, at 74.3° (`matsumoto2017electron`); DSA in 1-D only, at 63° and 75° (`xu2020electron`; `kumar2021nonthermal`) |
| Oblique, 50–75° | mirror reflection; it falls from 5 % at 30° to none at 74.3° (`morris2022preacceleration`) | confinement in whistler cavities of the foreshock, energy gain at the shock | 0.8 % of the electrons with 7 % of the energy at 60°; no DSA (`morris2023preacceleration`) |
| Quasi-parallel, below 50° | mirror reflection and shock drift | trapping between the shock and proton-driven waves | DSA with p⁻⁴ and K_ep = 3.8 × 10⁻³ in 1-D at 30° (`park2015simultaneous`); at most 0.5 % nonthermal by number (`gupta2024electron`); other mechanisms compete (`shalaby2022mechanism`; `zekovic2025slams`) |

Every result that reaches high energy is a chain: a first stage lifts electrons out of the thermal pool, and a
second takes over (Sect. 4.4).

## What decides the answer in a simulation (Sect. 7)
- **Mass ratio.** The thresholds scale differently. At the real mass ratio, surfing needs M_A above about 180
  (`wieland2016nonrelativistic` Sect. 2.2, from `matsumoto2012electron` Eq. 8). Whistlers need M_A below about 43
  (our arithmetic from `riquelme2011electron` Sect. 3.1). The two cannot hold at the same Mach number.
- **Shock speed.** Most runs use more than 0.1 c; real shocks have 0.007–0.02 c. At 0.26 c, 26 % of a sphere is
  superluminal, so that no electron can outrun the shock there; at 5000 km/s it is under 2 % (our arithmetic).
- **Dimension.** Each 2-D field orientation holds only part of the instabilities. Only 3-D holds all, and one
  3-D run exists (`matsumoto2017electron`).

## What is supported, disputed and unknown

**Supported by several works**
- Electrons need a pre-acceleration stage, and the mechanism changes with θ_Bn.
- As θ_Bn grows, fewer electrons are reflected and each gains more (`amano2007electron` Sect. 1;
  `morris2022preacceleration`).
- In hybrid simulations, protons are accelerated efficiently and drive strong upstream waves only below about
  45°, and not at all above 60° (`caprioli2014simulations` Fig. 3).
- In 1-D, electrons reach DSA at 30° with the help of protons and at 63° without it (`park2015simultaneous`;
  `xu2020electron`).

**Disputed** (Sect. 8.1)
- Whether surfing matters at real parameters.
- Whether quasi-parallel or quasi-perpendicular regions inject more electrons: 1-D simulations give more at 63°
  than at 30°, and the models of SN 1006 need the opposite.
- The mechanism at parallel shocks, and whether a reduced mass ratio can be trusted.

**Unknown, by the works' own words** (Sect. 8.2)
- The efficiency law ε(θ_Bn) for electrons. Every model of a whole remnant assumes it.
- Which mechanism dominates at the real mass ratio and shock speed.
- How pre-acceleration passes into DSA in more than one dimension.

## What it means for this topic (reading of the survey, Sect. 8.3)
- No kinetic simulation of electron injection on a front curved on the large scale was found for SNR-like
  parameters. The works on the curved bow shock of the Earth are not read, so this is not yet certain.
- No study covers the range of θ_Bn with one set of parameters.
- On a sphere in a uniform field, the three regimes of the table lie on one front: 29 % of the area has
  θ_Bn below 45°, 45 % lies between 45° and 75°, and 26 % above 75° (our arithmetic).
- Which works a spherical run can be compared with depends on its Mach number: SNR shocks have M_A of tens to
  hundreds.
- The research question of the topic is not decided by the survey; it is the owner's decision.

## Limits of the survey (Sect. 9)
- Only arXiv could be reached. 54 works were read.
- Not read: the bow-shock literature on curved shocks, the theory of the 1980s and 1990s, and several
  observational papers.
- Not covered: relativistic shocks and wakefield acceleration.
- The survey has not yet been reviewed.
