# Vink (2012): Supernova remnants: the X-ray perspective

- Ref: Astronomy and Astrophysics Review 20, 1 (January 2012; as given in the arXiv metadata); DOI
  10.1007/s00159-011-0049-1; arXiv:1112.0576v2. Published invited review (78 pages).
- Read: 2026-11-12, sections 1, 5.1-5.3, 5.5, 6.3, 11.1-11.4 and 12 (not read: ejecta, abundances,
  supernova types, mature remnants); topics: spherical-shock-injection

## Claim
X-ray data show that fast SNR shocks leave electrons much cooler than ions, and that thin X-ray
synchrotron rims imply 10-100 TeV electrons in fields amplified to 100-600 muG, hence efficient
cosmic-ray acceleration. How particles are injected is called one of the outstanding uncertainties.

## Method
Review of observations (Chandra, XMM-Newton, Suzaku, optical/UV line widths, gamma rays) with tutorial
derivations (jump conditions with cosmic-ray pressure, test-particle DSA, synchrotron loss scales). No
simulation of its own (m_i/m_e, M_A, beta, theta_Bn, box, particles per cell do not apply); hybrid and PIC
runs are only cited, without parameters (Cargill & Papadopoulos 1988; Shimada & Hoshino 2000; others).

## Numbers to compare against
- Mean post-shock temperature kT = (3/16) mu m_p Vs^2 = 1.2 (Vs/1000 km/s)^2 keV, mu = 0.6 (Eq. 20);
  with cosmic-ray pressure fraction w: kT2 = (1-w)(1/chi)(1-1/chi) mu m_p Vs^2 (Eq. 15).
- Shock thickness 10-100 c/omega_pe, about 1e7 n_e^(-1/2) cm; Coulomb mean free path 1e20/n_p cm (Sect. 5.2).
- Expected m_e/m_p << Te/Tp < 1; Shimada & Hoshino (2000) quoted with Te/Tp = 17 % (Sect. 5.3).
- Balmer shocks: equilibrated below Vs = 500 km/s, not above; Te/Tp ~ 1/Vs^2 (Ghavamian et al. 2007)
  is not confirmed by van Adelsberg et al. (2008) and Helder et al. (2011) (Sects. 5.3, 11.1).
- No young SNR shows kTe > 5 keV, although Cas A (Vs = 5100 km/s) and Tycho (4000 km/s) should have
  kT > 19 keV by Eq. 20 (Sect. 11.1).
- SN 1006: NE knot kT(O VII) = 528 +- 150 keV against kTe = 1.5 keV, Vs about 4000 km/s; NW UV lines
  of H, He, C, N, O all 2300 km/s FWHM, so ions are not equilibrated (Sect. 11.1). Radius 9.3 pc at
  2.0 kpc; Tycho 3.8 pc at 3 kpc (Table 2).
- Coulomb: tau_ep = 3.1e11/n_p (kT/1 keV)^1.5 s (Eq. 25); full equilibration for n_e t > 1e12 cm^-3 s (Sect. 5.3).
- Bykov (2004): electron thermal speed exceeds the shock speed for M < sqrt(m_p/m_e) = 43, i.e.
  Vs < 400 km/s for a sound speed of 10 km/s (Sect. 11.1).
- DSA: Delta E/E = (4/3)((chi-1)/chi) Vs/c (Eq. 33); q = (chi+2)/(chi-1) (Eq. 37); D = eta E c/(3eB)
  (Eq. 39); tau_acc = 124 eta (B/100 muG)^-1 (Vs/5000 km/s)^-2 (E/100 TeV) yr at chi = 4 (Eq. 40).
- Synchrotron: h nu = 13.9 (B_perp/100 muG)(E/100 TeV)^2 keV (Eq. 47); tau_syn = 12.5 (E/100 TeV)^-1
  (B/100 muG)^-2 yr (Eq. 49); loss-limited cut-off h nu = 1.4/eta (Vs/5000 km/s)^2 keV at chi = 4,
  independent of B (Eq. 52). X-ray synchrotron needs Vs > 2000 km/s and eta < 10 (Sect. 6.3.1).
- Rims 1-2 arcsec wide in Cas A, Tycho, Kepler; B2 = 26 (l_adv/1e18 cm)^(-2/3) eta^(1/3) muG at
  chi = 4 (Eq. 62); B = 100-600 muG against B_ISM = 5 muG (Sects. 11.2, 11.3).
- Radio index alpha = 0.6, electron index q = 2.2; X-ray photon index 2-3.5 (Sect. 6.3.1).
- Cosmic-ray pressure > 50 % of the post-shock pressure in RCW 86 NE (Helder et al. 2009), > 20 % in
  B0509-67.5; Cas A cosmic-ray energy < 4 % of the explosion energy (Fermi) (Sects. 11.3, 11.4).
- Suprathermal electron tails thermalise by Coulomb collisions within n_e t = 1e10 cm^-3 s (Sect. 6.3.2).

## Assessment
- Shown by observation (as reviewed): Te << Ti behind fast shocks; narrow synchrotron rims and the
  field strengths they imply; in SN 1006 synchrotron X-rays from two caps only.
- Stated or conjectured: for SN 1006 the "inferred polar geometry" means X-ray synchrotron emission is
  favoured where B is parallel to the shock normal; the need for turbulence is reconciled by assuming
  that acceleration starts at a parallel field and then makes its own turbulence (Voelk et al. 2003)
  (Sect. 11.3). A link between equilibration and cosmic-ray acceleration is called suggestive only.
- On injection only a few sentences: thermal-tail injection is "often assumed"; simulations show a non-thermal
  seed population behind the shock (Bykov & Uvarov 1999; Riquelme & Spitkovsky 2011 for electrons; Sect. 5.5);
  electron injection is "poorly known" (Sect. 6.3.2). No K_ep is given; cosmic rays at Earth are 99 % nuclei.
- Stated limits: Balmer-line results are model dependent and ignore precursors; X-ray ion temperatures
  so far concern reverse-shocked ejecta; the scaling of B^2 with Vs (Vs^2 or Vs^3) is uncertain; rim
  widths could reflect field damping (Pohl et al. 2005), which the radio profiles disfavour.
- Nothing on shock-front curvature, rippling or non-stationarity.
- Slip: Sect. 5.3 prints "m_p/m_e = 0.01-1" for PIC mass ratios (m_e/m_p is meant).
