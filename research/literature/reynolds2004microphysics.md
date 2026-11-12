# Reynolds (2004): Microphysics of shock acceleration from observations of X-ray synchrotron emission from supernova remnants

- Ref: Advances in Space Research (Proc. 34th COSPAR Assembly, Symposium E1.4; "to be published" in
  the arXiv comment; volume and pages not in the metadata - journal reference not checked); DOI
  10.1016/j.asr.2003.02.013; arXiv:astro-ph/0308217v1. Conference proceedings paper.
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
The synchrotron X-ray cut-off of a remnant measures the acceleration rate and so the diffusion
coefficient kappa ~ E^beta. Model images of a spherical remnant in a uniform field look unlike real
remnants for beta < 1 (Kolmogorov, Kraichnan), so the scattering waves are probably self-generated
(k^-1 spectrum, beta = 1).

## Method
Theory and synthetic images; no kinetic simulation (m_i/m_e, plasma beta, particles per cell do not
apply). Spherical Sedov blast wave in a uniform medium with a uniform upstream field B1 inclined by phi
to the sky plane, so theta_Bn runs from 0 (poles) to 90 deg (equator). Each shocked fluid element gets
N(E) = K E^-s exp(-E/Emax), with Emax set by shock speed, age, kappa and theta_Bn, then evolved with
radiative and adiabatic losses (code of Reynolds 1998, generalised to kappa ~ E^beta). Standard picture
assumed valid for strong shocks, M_A and M_s > 10. Shock speed of SN 1006: "known", value not given;
RCW 86: 600-900 km/s. The normalisation K (injection) is given no obliquity dependence in this paper.

## Numbers to compare against
- kappa_par = eta r_g c/3, eta = gyrofactor, eta = 1 is the Bohm limit (Eq. 1); kappa_perp =
  kappa_par/(1+eta^2) (Eq. 2); kappa = kappa_par [cos^2 theta_Bn + sin^2 theta_Bn/(1+eta^2)] (Eq. 3).
- Wave spectrum I(k) ~ k^-n gives beta = 2-n: Kolmogorov 1/3, Kraichnan 1/2, constant gyrofactor 1.
- Acceleration time: parallel shock tau ~ E^beta (Eq. 6); perpendicular shock tau ~ E^(2-beta)
  (Eq. 7); tau(90)/tau(0) ~ (E/E_h)^(2-2 beta), E_h = energy at which mean free path = gyroradius.
- SN 1006: roll-off frequency 7e17 Hz; with mean B = 10 muG (from TeV inverse Compton) Emax = 53 TeV;
  synchrotron half-life at Emax > 1200 yr, acceleration time < 800 yr: neither losses nor age set the
  cut-off, so escape above a wavelength of 2e17 cm (about 0.1 pc) is inferred.
- SN 1006 halo: the model predicts upstream emission at 5-10 % of the post-shock peak (3 % for an
  edge-on plane shock); the Chandra limit is 1 %. A thin halo (< 10 arcsec) would need an upstream
  field of at least 20 muG, judged unreasonable 570 pc above the Galactic plane.
- RCW 86: nu_c = (7-10)e16 Hz; loss-limited nu_c = 5e16 eta (B2/4B1)(u_sh/1000 km/s)^2 Hz (Eq. 4),
  so eta > 1 and/or compression r > 4 is needed even for a perpendicular shock.
- All shell remnants examined (14 Galactic, 11 in the LMC) are fainter in X-rays than their
  extrapolated radio spectrum: the electron spectrum cuts off below X-ray-emitting energies.
- Fig. 2: four model images for different beta, field at phi = 30 deg to the sky plane (the beta
  values of the panels are not legible in the text version).

## Assessment
- Shown by theory and model: with kappa_perp << kappa_par, acceleration to a given energy is faster
  at quasi-perpendicular shocks (Jokipii 1987), so Emax is highest on the equatorial belt; this
  preference becomes extreme for beta < 1.
- Stated, not shown here: in SN 1006 "the shock is almost certainly perpendicular at the bright
  non-thermal rims" (Conclusions). This is the equatorial-belt reading; Rothenflug et al. (2004),
  bocchino2011constraints and reynoso2013radio conclude the opposite (polar caps).
- Conjecture: the missing halo requires a field jump larger than the fluid compression, i.e.
  post-shock turbulent amplification (cites Lucek & Bell 2000).
- Limits: test-particle picture; obliquity enters only through the acceleration rate, not through the
  number of injected electrons; a well-defined theta_Bn needs weak turbulence (large eta), while for
  eta near 1 the obliquity dependence is small anyway.
- Differs from later work: mean field 10 muG here against 100 muG or more from rim widths
  (vink2012supernova Sect. 11.2; bocchino2011constraints Sect. 1).
