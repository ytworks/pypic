# West, Safi-Harb & Ferrand (2017): The connection between supernova remnants and the Galactic magnetic field: An analysis of quasi-parallel and quasi-perpendicular cosmic ray acceleration for the axisymmetric sample

- Ref: A&A 597, A121 (2017); DOI 10.1051/0004-6361/201628079; arXiv:1609.02887v1; peer-reviewed (accepted 8 Sep 2016).
- Read: 2026-11-12, full (Appendix A image strips not readable in the text file); topics: spherical-shock-injection

## Claim
For 32 radio SNRs with bilateral or one-sided shells, synthetic synchrotron images in a Galactic field model match the
observed orientation when the electron density is weighted towards quasi-perpendicular parts of the shock, and mostly
fail when it is weighted towards quasi-parallel parts. SN 1006 is the exception: its radio polarisation is better
matched by the quasi-parallel weighting.

## Method
- Model, no kinetic simulation (m_i/m_e, Mach numbers, beta, particles per cell: do not apply).
- A bubble-like Sedov-phase remnant is inserted by a coordinate transformation into a Galactic magnetic field model
  (Jansson & Farrar 2012, "JF12"; Sun et al. 2008 as a check); the frozen-in field is dragged into its post-shock
  configuration; the Hammurabi code gives Stokes I, Q, U (Sect. 2).
- Electron weighting (Sect. 2): density scaled by cos^2(phi_Bn2) (quasi-parallel) or sin^2(phi_Bn2)
  (quasi-perpendicular); phi_Bn2 is the angle between shock normal and POST-shock field. Source of the scaling:
  Leckband et al. (1989), Fulbright & Reynolds (1990). No threshold angle, no proton model.
- Each SNR modelled at 11 distances (0.5, 1, 2, ... 10 kpc); the angle psi of the bilateral axis to the Galactic plane
  is compared with the radio image by eye (Sect. 3; Appendix A).
- SN 1006 in addition: uniform field (Bx line of sight, By horizontal, Bz vertical) compared with total intensity,
  polarised intensity and B-vectors of Reynoso et al. (2013) (Sect. 4.1, Fig. 3).

## Numbers to compare against
- Sample: 33 SNRs minus G1.9+0.3 = 32 (Sect. 2).
- Quasi-perpendicular: about 75 % of the models match well, most at published distances (Sect. 5).
- Quasi-parallel: no reasonable match in at least 10 of 32 cases (Sect. 3; 10 named in Appendix A); 7 cases match in
  shape but at a distance inconsistent with the published one; 3 match in both (G65.1+0.6, G156.2+5.7, G332.0+0.2),
  and only 2 of them where quasi-perpendicular fails (Sect. 3).
- Measured field directions: 13 of 15 cases agree better with the quasi-perpendicular model; the exceptions are
  G46.8-0.3 and SN 1006 (Sect. 3, Fig. 2).
- SN 1006: diameter 30 arcmin, distance 1.6-2.2 kpc (upper limit 2.1 kpc), bilateral axis at 83 +/- 5 deg to the
  Galactic plane (Sect. 4); age 1010 yr, Type Ia (Table 1).
- SN 1006 in JF12: quasi-perpendicular fit at 1 +/- 1 kpc; quasi-parallel has no fit at any distance, 0.5 kpc is
  closest in orientation (Sect. 4).
- Uniform-field models (Sect. 4.1, Fig. 3): quasi-parallel By = 1.0, Bz = 0.12 microgauss, Bx < By (0.03 or 0.55)
  gives radial B-vectors in the limbs as observed; quasi-perpendicular Bz = 1.0, By = -0.12, Bx = 1.28 microgauss
  (close to JF12 at SN 1006) gives tangential vectors; with Bx = -12.5 the vectors turn radial but the intensity
  becomes ring-like. tan(psi) = Bz/By (perp), -By/Bz (par).
- Limb-to-centre ratio from Rothenflug et al. (2004): 0.3 X-ray, 0.7 radio; an equatorial belt needs >= 0.5 (Sect. 3.1).

## Assessment
- Shown (model against images, qualitative): the orientation statistics above. The comparison is by eye; brightness
  ratios and azimuthal profiles are deliberately not used (Sect. 3.1).
- Stated, not shown: SN 1006 may combine both geometries: quasi-perpendicular electron acceleration in the compressed
  Galactic field sets the radio shape, turbulence in the young remnant makes a radial field, and quasi-parallel ion
  acceleration then gives the gamma rays at the limbs (Sect. 4.2, 5).
- Stated limits: no field amplification or turbulence (Sect. 4.2); Sedov phase assumed, doubtful for SN 1006 and
  G1.9+0.3 (Sect. 2, 4.3); JF12 has no local structure (Sect. 4.1); psi cannot be defined for many quasi-parallel
  models with filled centres (Sect. 3).
- The weighting uses the post-shock angle, not the upstream theta_Bn; the paper gives no conversion.
- Neither weighting reproduces the observed B-vector pattern of SN 1006 convincingly: the quasi-parallel vectors
  converge at the poles, which is not observed (Sect. 4.2).
- The introduction cites Matsumoto et al. (2012, 2013) and Caprioli & Spitkovsky (2014) for electrons being
  accelerated best at quasi-perpendicular and ions at quasi-parallel shocks; this is cited, not tested here.
- Disagreements it lists for SN 1006: quasi-perpendicular (Petruk et al. 2009; Schneiter et al. 2010) against
  quasi-parallel (Rothenflug et al. 2004; Bocchino et al. 2011; Schneiter et al. 2015).
- Internal tension: Sect. 3 says quasi-perpendicular matches "for all cases" at some distance; Sect. 5 says about 75 %.
