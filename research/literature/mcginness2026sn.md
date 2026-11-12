# McGinness, Diesing, Caprioli & Acero (2026): SN 1006: A Cosmic Laboratory for Investigating Shock Acceleration Physics

- Ref: DOI 10.3847/1538-4357/aea3f3 (an Astrophysical Journal DOI; volume and page not in the metadata: journal
  reference not checked); arXiv:2606.12111v2, draft of 30 Sep 2026; the acknowledgements thank a referee.
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
A semi-analytic non-linear shock acceleration model, run separately for four quadrants of SN 1006, reproduces the
spectra and the X-ray profiles with a cosmic-ray efficiency of 20 % in the quasi-parallel polar caps and below 1 % in
the quasi-perpendicular quadrants. Most gamma rays come from electrons; the north-west, where the shock hit a cloud,
is hadronic.

## Method
- Model, no kinetic simulation (m_i/m_e, beta, particles per cell do not apply; Mach number not stated).
- Hydrodynamics: 1-D, spherical, per quadrant. E_SN = 1e51 erg, M_ej = 1.2 M_sun, age 1020 yr, distance 2.2 kpc, B_0 =
  3 microgauss, T = 1e4 K. Ejecta-dominated stage from Tang & Chevalier (2017), then a thin-shell Sedov stage (Sect.
  2.1, Eq. 1). NW: density step to n_cl at R_cl = 5 pc.
- Acceleration: multizone kinetic model of non-linear diffusive shock acceleration at a quasi-parallel shock (Caprioli
  2012; Diesing & Caprioli 2019, 2021): proton injection momentum from hybrid simulations, Bell instability driven by
  escaping particles (Eq. 3), drift of the scattering centres at the Alfven speed up- and downstream (Eq. 4),
  free-escape boundary at 0.02 R_sh (Sects. 2.2, 2.3). This quasi-parallel model is applied to all quadrants.
- Obliquity enters only through per-quadrant choices: proton injection is turned down outside the polar caps (Sect.
  2.2), and diffusion is slowed by a factor eta = D/D_Bohm (Sect. 2.4).
- Electrons: spectrum after Zirakashvili & Aharonian (2007), tied to the protons by one ratio K_ep, the same in all
  quadrants (Sect. 2.2, Table 1). Radiation with naima. Filling factors f.

## Numbers to compare against
- Table 1, NE/SW | SE | NW: v_sh 5300 | 5300 | 2700 km/s; R_sh 9.7 | 9.7 | 7.8 pc; n_ISM 1.8e-2 cm^-3 (NW: step to
  0.18); xi_CR 20 % | 0.4 % | 0.8 %; B_2 40 | 6 | 10 microgauss; R 5.0 | 4.0 | 4.0; eta 1.0 | 8.1 | 8.1; f 0.1 | 0.4 |
  0.4; K_ep 9.3e-4 in all.
- Observational inputs (Sect. 2.1, Fig. 2): v_sh = 5000 +/- 400 km/s (NE, SW, SE) and 3000 +/- 400 km/s (NW); R_sh =
  9.1 +/- 0.9 pc and 7.6 +/- 0.8 pc; these allow n_ISM = 0.018-0.034 cm^-3 and n_cl = 0.12-0.18 cm^-3.
- Proton and uncooled electron spectra E^-2.18 from the radio slope (Sect. 3.1); E_max about 1e5 GeV (Sect. 2.3).
- K_ep about 9.3e-4, against about 1.1e-3 from the formula of Gupta et al. (2025) for these parameters (Sect. 3.1).
  Typical range quoted for remnants: 1e-2 to 1e-4 (Sect. 2.2).
- Inside a cap (Sect. 3.2): centre with xi_CR = 20 % and edges with xi_CR = 0.5 %, deltaB/B_0 about 0.55, each with f
  = 0.05; then K_ep = 1.9e-3. The spectrum is nearly unchanged.
- eta: X-ray cutoff energy falls by about 1.8 over about 20 deg (Miceli et al. 2016); extrapolated over about 90 deg
  this gives eta from 1 to about 8.1 (Sect. 2.4).
- Radio rim (Sect. 3.4, Fig. 5): damping at 2e-2 per yr fits the profile but then needs K_ep = 1.9e-2 for the radio
  and 2.8e-4 for the gamma rays, so damping is rejected.
- Shock/contact-discontinuity radius ratio: model 1.12 (caps), 1.16 (SE); observed about 1.00, 1.04 (Sect. 3.4).

## Assessment
- Shown by fitting: one parameter set per quadrant reproduces the spectra (Fig. 3) and the X-ray rim profiles of the
  caps (Fig. 4). B_2 about 40 microgauss accounts for slope, compression ratio, rim width and synchrotron
  normalisation together (Sect. 3.3).
- Not tested: an obliquity dependence of electron injection separate from that of protons. K_ep is held constant, so
  the electron efficiency follows the proton efficiency by construction. The paper says electron injection "is not
  completely understood" (Sect. 2.2).
- Stated limits: the modelled emission outside the caps is an upper limit (Fig. 3 caption) and unconstrained outside
  the Fermi band; diffusion there may not be a rescaled Bohm coefficient (Sect. 2.4); radio and X-ray emission of the
  forward shock outside the caps is assumed negligible; the radio profile is not reproduced and is attributed to 3-D
  effects at the contact discontinuity (Sect. 3.4).
- Stated, not shown: particles in the quasi-perpendicular NW may come from shock drift acceleration (Orusa & Caprioli
  2023; Orusa et al. 2026); re-acceleration of Galactic cosmic rays "cannot be ruled out" (Sect. 3.2).
- Reader's comparison (the paper cites Winner et al. 2020 without discussing it): Winner et al. need a damped Bell
  field and hadronic TeV emission at density 0.12 cm^-3; here no damping, leptonic TeV emission, 0.018 cm^-3. Giuffrida
  et al. (2022) infer a compression ratio up to about 7; here R = 5.0, which the paper calls consistent (Sect. 3.2).
