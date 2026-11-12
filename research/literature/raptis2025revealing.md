# Raptis et al. (2025): Revealing an unexpectedly low electron injection threshold via reinforced shock acceleration

- Ref: Nature Communications 16(1), 488 (2025); DOI 10.1038/s41467-024-55641-9; arXiv:2502.10643 (v1, read here);
  peer-reviewed publication
- Read: 2026-11-12, introduction, results, discussion and methods (tables skimmed); topics: spherical-shock-injection

## Claim
In the foreshock upstream of Earth's bow shock, a foreshock transient accelerated electrons to above 500 keV when
the solar wind already carried an enhanced 1–5 keV electron population. The authors propose a multi-step
"reinforced shock acceleration": a suprathermal seed, shock acceleration at the transient's own compressive
shock, resonant energization by intense whistlers there, and confinement between the transient's turbulent core
and the bow shock. They conclude that the seed energy needed at Earth is only ~1–10 keV and that the maximum
energy grows with the size of the system.

## Method
- Observation, no simulation and no transport model. MMS (FPI, FEEPS, FGM, SCM; burst data) just upstream of
  the bow shock at (10.8, 11.1, 5.1) RE GSE, and ARTEMIS P2 near the Moon at (65.2, −5.3, 5.1) RE as far-upstream
  monitor; OMNI for solar-wind context. Main event 2017-12-17, 17:52:45–17:53:25 UT.
- Comparison set: six events with a seed and six without (Fig. 3, Tables II–IV); solar-wind origin from the
  four-class scheme of Xu & Borovsky with a Gaussian-process classifier.
- Shock at the transient's sunward compressive edge: quasi-perpendicular, θBn above 55° by coplanarity, timing
  and minimum variance; speed ~300 km/s in the spacecraft frame (Methods). Alfvén Mach number and plasma beta of
  the transient's shock or of the bow shock: not stated.
- Wavelet and singular-value wave analysis; Bohm diffusion length L ~ (1/3)(v/V_sh) r_ce (Methods, Eq. 4).

## Numbers to compare against
- Electron flux enhancement to more than 500 keV; highest energy 523 keV (Results; Methods).
- Fitted index of the phase space density versus energy reaches p ~ −2 for 80–200 keV (Fig. 1h, Results).
- Seed: 1–5 keV flux at ARTEMIS 2–5 times the background in every seeded event, all in fast coronal-hole wind;
  no relativistic electrons in the unseeded slow-wind events (Fig. 3, Results).
- Efficiency: electron energy density above 100 eV over the solar-wind ion ram energy density ≈ 5 % (n = 10 cm⁻³,
  600 km/s); below 0.5 % if counted above 1 keV (Methods, Eq. 3).
- Whistlers: narrowband at 0.1–1 fce, right-handed, parallel-propagating (Fig. 2); δB/B = 0.015–0.13 over the
  twelve events, 0.060 for the main event (Table III); ω ~ 0.3 Ωce, wavelength ~25 km; resonant pitch angles
  ~79° at 20 keV and ~83° at 50 keV; pitch-angle distribution of 40–200 keV electrons peaks at 90° (Methods, Fig. 2).
- fpe/fce ~ 30 in the main event (Methods). Core field ~1 nT, up to 10 nT in extreme cases (Methods).
- Diffusion length for 0.5 MeV: 10–100 RE, the size of foreshock transients at Earth; 100–1000 RE and 0.5–1 nT
  (Jupiter-like) give above 2 MeV; fields of 0.01–1 G near ultra-hot Jupiters give GeV–TeV (Results; Methods).
- Field lines map unambiguously to the bow shock in 2 of 6 seeded and 1 of 6 unseeded events (Fig. 9, Methods).

## Assessment
- Shown (observation): one event with electrons above 500 keV inside a foreshock transient, together with a
  quasi-perpendicular shock at its edge, intense whistlers, and a 90°-peaked distribution; a twelve-event
  comparison in which the far-upstream suprathermal seed separates energetic from non-energetic transients.
- Not shown: which process supplies the energy; the authors say the shares of shock acceleration, betatron,
  Fermi and wave resonance cannot be separated observationally. "Injection threshold" here is the seed energy
  that precedes a successful event, not a measured resonance threshold. p ~ −2 is likened to DSA; no model is fitted.
- Stated/conjectured: generalization to other planets and stars; a sizeable share of cosmic-ray electrons from
  planetary quasi-parallel shocks in stellar winds; support for quasi-parallel shocks as the source. These rest
  on a gyroradius-versus-size estimate with Bohm scattering. SNR shocks enter only as "MA > 10".
- Stated limits: small sample set by spacecraft conjunctions; fast wind necessary but not sufficient; seed origin
  (solar or electron foreshock) not settled; quantification left to global kinetic simulations.
- Geometry: confinement relies on field lines of the discontinuity that meet the bow shock at a point moving
  slowly across it as they convect (Methods, Fig. 9A).
- Relation to other work: reports a harder spectrum than expected for quasi-perpendicular SDA, citing Amano et
  al. (2020); that paper expects no relativistic electrons at ordinary heliospheric shocks, which this event does
  not contradict because it needs a transient plus a seed. Re-reads the events of Wilson et al. (2016) and Liu et
  al. (2019) as fast-wind cases.
- Doubts: 20 keV electron speed printed as 8.4e7 km/s (~8e4 km/s expected); diffusion length 6e4–6e5 km in the
  Results but 5e4–5e5 km in the Methods; the frequency statement after "ω ~ 0.3 Ωce" is unclear.
