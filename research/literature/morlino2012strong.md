# Morlino & Caprioli (2012): Strong evidences of hadron acceleration in Tycho's supernova remnant

- Ref: Astronomy & Astrophysics (accepted, per arXiv comment; volume and article number not in the
  metadata - journal reference not checked); DOI 10.1051/0004-6361/201117855; arXiv:1105.6342v2.
  Peer-reviewed article.
- Read: 2026-11-12, introduction, set-up (Sects. 2-3), electron treatment (Sects. 2.5, 4.1-4.4) and
  conclusions; topics: spherical-shock-injection

## Claim
A semi-analytic non-linear DSA model with self-generated field amplification fits Tycho from radio to
TeV: protons take about 12 % of the kinetic energy flux and reach 500 TeV, and the gamma rays are
hadronic. The electron normalisation is not predicted: K_ep = 1.6e-3 is fitted to the synchrotron data.

## Method
Model, not a kinetic simulation (m_i/m_e, particles per cell, box do not apply). Stationary NLDSA for
quasi-parallel, non-relativistic shocks (Caprioli, Amato & Blasi 2010) applied in time steps on the
Truelove & McKee (1999) evolution of a Type Ia remnant (1e51 erg, 1 Msun of ejecta, uniform medium).
Spherical symmetry; no theta_Bn dependence. Proton injection by thermal leakage (Blasi, Gabici &
Vannoni 2005): particles with p > xi_inj times the downstream thermal momentum are injected. Bohm
diffusion in the amplified field; free-escape boundary at 0.1 R_sh; resonant streaming instability
saturating at P_w = P_cr/2M_A; Alfvenic drift in the amplified field; non-linear Landau damping
downstream. Electrons: proton spectrum times K_ep with a synchrotron-loss cut-off (Eq. 12), evolved
downstream with adiabatic and radiative losses (Eqs. 13-14). Te = Tp m_e/m_p at the shock, then
Coulomb heating only. Mach numbers are V_sh/c_s and V_sh/v_A, V_sh = shock speed in the ambient medium.

## Numbers to compare against
- Environment: n0 = 0.3 cm^-3, T0 = 1e4 K (Sect. 2.1); far-upstream field about 5 muG (Sect. 5).
  Literature limits n0 < 0.6 and < 0.4 cm^-3 (Sect. 3).
- At 439 yr: V_sh = 4990 km/s, R_sh = 3.94 pc (Eqs. 1-2), distance 3.3 kpc; X-ray proper motion gives
  4600 +- 400 (d/2.3 kpc) km/s; Halpha widths give 1500-2800 km/s (Sect. 3).
- Sonic Mach number M_0 = 443; Alfvenic Mach number in the amplified field M_A,1 = 15 at xi_cr = 0.06
  (Eq. 7, Sect. 2.3). Field about 90 muG just upstream, B2 = sqrt(11) B1 = 300 muG downstream (Sects. 3, 5).
- xi_inj = 3.7, i.e. a fraction of about 1e-4 of the particles crossing the shock (Sect. 2.2).
- Cosmic-ray pressure 0.06 rho0 V_sh^2, about 12 % of (1/2) rho0 V_sh^3 (Sect. 4.4); R_tot = 4.2,
  R_sub = 3.9 (Sect. 3); proton spectrum p^-4.2 (E^-2.2), p_max = 470 TeV/c (Fig. 4, Sect. 4.4).
- K_ep = 1.6e-3 from the radio and Suzaku X-ray data together (Sect. 3); "a factor a few" below the
  ratio in Galactic cosmic rays at Earth.
- Electron spectrum p^-4.2 below p_roll = 200 GeV/c, p^-5.2 up to p_cut = 20 TeV/c, cut-off shape
  exp[-(p/p_e,max)^2] (Sect. 2.5).
- Post-shock T_p = 5e8 K; Coulomb-heated Te reaches 2e7 K near the contact discontinuity (Fig. 3).
- Radio: 40.5 Jy at 1.4 GHz, index 0.6-0.65 (Sect. 1); model index 0.606 (Sect. 4.1). Fitting radio
  and X-rays together needs B2 > 200 muG (Fig. 10, Sect. 4.3).
- nu_syn = 24 (B/100 muG)(gamma/100)^2 MHz (Eq. 15): electrons near injection (gamma = 10-100) emit at
  1-100 MHz in 100-300 muG.
- Contact discontinuity at 0.87 R_sh in the model against 0.93 +- 2 % observed (Warren et al. 2005) (Sect. 3).

## Assessment
- Shown by the model fit: a single K_ep reproduces radio and X-ray synchrotron fluxes; inverse Compton
  on the CMB is 20-30 times below pion decay (Fig. 11), so a leptonic origin fails in flux and slope.
- Stated: K_ep "very likely" reflects different injection mechanisms for leptons and hadrons; electron
  injection is "far from being understood from first principles" (cites Amano & Hoshino 2010;
  Riquelme & Spitkovsky 2010; Morlino 2009). K_ep and xi_inj are free parameters for that reason.
- Conjecture: the low-frequency radio curvature (Reynolds & Ellison 1992) may probe electron injection
  rather than a precursor (Sect. 4.1).
- Limits: stationary NLDSA; quasi-parallel theory applied to the whole sphere; neutrals ignored;
  K_ep taken constant in time although "in principle" time dependent; the steep spectrum rests on
  taking the Alfven speed in the amplified field.
- K_ep and B are separated in the model: when synchrotron losses dominate, the X-ray flux in the cut-off region
  depends only on the electron density and fixes K_ep independently of B; the radio flux then needs
  B_2 of at least about 200 μG (Sect. 4.3). (Corrected 2026-11-12: this line said K_ep is degenerate with B.)
- Slips: Sect. 4.4 prints "A Kep smaller by a factor 20-30" where larger is meant; the last paragraph
  prints "500 GeV" for 500 TeV; 1/(gamma M_0^2) with M_0 = 443 is 3.1e-6, the text prints 5.6e-6.
