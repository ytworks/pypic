# Winner, Pfrommer, Girichidis, Werhahn & Pais (2020): Evolution and observational signatures of the cosmic ray electron spectrum in SN 1006

- Ref: MNRAS 499, 2785 (2020); DOI 10.1093/mnras/staa2989; arXiv:2006.06683v2; peer-reviewed.
- Read: 2026-11-12, introduction, Sects. 2-4, 7; Sects. 5-6 skimmed; Appendix not read; topics: spherical-shock-injection

## Claim
In a 3-D MHD model of SN 1006 with sub-grid acceleration, only an electron efficiency that is highest at
quasi-parallel shocks reproduces the radio and X-ray maps; an efficiency that is constant or highest at
quasi-perpendicular shocks does not. The fit also needs a volume-filling field of about 35 microgauss and a
Bell-amplified field that decays right behind the shock.

## Method
- 3-D ideal MHD on a moving mesh (AREPO) with cosmic-ray protons as a second fluid (adiabatic index 4/3), injected at
  shocks found on the fly; electron spectra evolved on tracer particles in post-processing (CREST): adiabatic changes,
  Coulomb, synchrotron, inverse-Compton (Sect. 2.1).
- Not a kinetic simulation: m_i/m_e, particles per cell do not apply; Mach numbers and beta not stated.
- Set-up (Sect. 2.3): 1e51 erg thermal energy in the central cell of a periodic 22 pc box, no ejecta phase; n = 0.12
  cm^-3, mu = 1.25, T = 5.1e3 K; uniform B = 1 microgauss along the diagonal of the plane of the sky; 1e6 cells
  (parameter study) or 5e6 cells (maps). Shock speed about 3000 km/s after 1000 yr (Sect. 2.4). Distance 1660 pc from
  matching the angular size (Sect. 3.1).
- Protons (Eq. 1): zeta_p = (zeta_max/2) [tanh((theta_crit - theta)/delta) + 1], zeta_max = 0.15, theta_crit = pi/4,
  delta = pi/18; theta is the upstream obliquity; form from Pais et al. (2018), after the hybrid simulations of
  Caprioli & Spitkovsky (2014a).
- Electrons (Eq. 4): zeta_e = ((zeta_par - zeta_perp)/2) [tanh((theta_crit - theta)/delta) + 1] + zeta_perp. Three
  cases (Fig. 1): zeta_par = 10 zeta_perp; constant; zeta_perp = 10 zeta_par. The larger of the two, zeta_e,max, is
  set by the radio flux. Efficiencies are fractions of the dissipated energy. The paper states that the form of Eq. 4
  is not known from first principles (Sect. 2.2).
- Field (Sect. 2.4): compression of the perpendicular component; post-shock field multiplied by a factor to mimic a
  turbulent dynamo; Bell amplification at the shock by f_Bell with the obliquity dependence of Eq. 1 (Eq. 5).

## Numbers to compare against
- Motivation from particle-in-cell work (Sect. 2.2): electron energy fraction about 1e-2 at quasi-perpendicular strong
  shocks (0.1 x 0.1: Xu et al. 2020 and a private communication) against at most about 1e-3 at quasi-parallel shocks
  (Park et al. 2015; Xu et al. 2020).
- Best fit (Sect. 3.1): zeta_e,max = 5e-4 (quasi-parallel case), spectral indices 2.1 (electrons) and 1.9 (protons),
  electron cutoff momentum 3.5e7 m_e c, f_Bell = 20, equivalent field 35 microgauss. The constant and the
  quasi-perpendicular case need zeta_e,max = 2e-4 (Sects. 4.2, 4.3).
- Field: 35 microgauss at parallel, 140 at perpendicular, 76 at theta = 30 deg (Sect. 7); about 4 microgauss at
  quasi-perpendicular shocks from compression alone (Sect. 2.4); dynamo ceiling about 250 microgauss (Eq. 6).
- Radio and X-ray brightness peak about 35 deg from the field axis, not at 0 deg (Sect. 3.2, Fig. 5).
- Result: efficiency for GeV electrons at quasi-perpendicular shocks lower "at least by a factor ten" (Sect. 7).
- Bell field must decay within about one cell, about 100 gyroradii of TeV particles (Sects. 5.2, 7).
- Pre-shock densities quoted from observations: 0.04, 0.05, 0.045 cm^-3; NW 0.35; SW 0.3-0.5 cm^-3 (Sect. 6.3).

## Assessment
- Shown by forward modelling: with this field model, the constant and the quasi-perpendicular case give bright
  south-east and north-west limbs in radio and X-rays, unlike the data (Figs. 7, 8).
- Only three cases with a fixed contrast of ten are tried; the factor is assumed, not fitted.
- Field direction rests on the gamma-ray map being hadronic above 500 GeV: a 90 deg turn would spoil it (Sect. 4.2).
- Reader's remark: the imposed post-shock field is four times stronger at perpendicular than at parallel shocks, which
  favours synchrotron emission there; the conclusion on electrons depends on this field model.
- Stated against kinetic work: the result "precludes extrapolation" of 1-D particle-in-cell runs favouring
  quasi-perpendicular electron acceleration (Riquelme & Spitkovsky 2011; Bohdan et al. 2017; Xu et al. 2020) (Sect. 7).
- Stated limits: acceleration, Bell growth and dynamo are sub-grid; uniform initial field, no upstream turbulence;
  density, explosion energy and indices are degenerate (Sect. 6.2).
- Model against data: shock 3000 km/s (about 5000 km/s observed, see the other notes); density 0.12 vs 0.04-0.05 quoted.
- Inconsistencies: zeta_e,max = 5.4e-3 for index 2.1 in Sect. 6.1 against 5e-4 in Sect. 3.1; f_Bell 10 in the Fig. 15
  caption against 20 in the text; GeV emission in equal parts (abstract, Sect. 7) against leptonic-dominated (Sect.
  3.1); post-shock densities of 1.44-1.99 cm^-3 next to a pre-shock 0.04 cm^-3 in Sect. 6.3.
