# Pohl, Hoshino & Niemiec (2020): PIC simulation methods for cosmic radiation and plasma instabilities

- Ref: Prog. Part. Nucl. Phys., article 103751 (from the DOI; volume not checked, the arXiv metadata give
  none); DOI 10.1016/j.ppnp.2019.103751; arXiv:1912.02673; peer-reviewed review (arXiv comment: accepted).
- Read: 2026-11-12, sections 1.2.3-1.3, end of 2.2 (limits of PIC), 2.4, 5 intro, 5.2-5.3 (remnant
  shocks), 6 (low Mach number, high beta), 7-7.2 (reformation, rippling), 8; topics: spherical-shock-injection

## Claim
Review. For remnant shocks: PIC runs show electron pre-acceleration by several processes inside the
shock transition, but electrons must reach a momentum far above m_i V_sh before diffusive acceleration
can take over, and the published runs are restricted in dimension, shock speed and mass ratio.

## Method
Review of PIC results by object class; hybrid, MHD-PIC and Vlasov-Fokker-Planck work is quoted for ions.

## Numbers to compare against
- Injection threshold: electron mean free path ~ gyroradius must exceed the ion gyroradius (proxy for the
  shock thickness), so p_e >> m_i V_sh: about 10 MeV/c for V_sh ~ 0.01 c, i.e. gamma_e >> 20 (Eq. 34).
- Reflected ion fraction about 0.2 (Sect. 5.3); 20-30 % for M_A > 3 (Sect. 7.1).
- Buneman waves grow if M_s >~ [(1+a)/2] sqrt(m_i/m_e) sqrt(T_e/T_i) (Eq. 35); surfing works if
  M_A >~ (1+a) (m_i/m_e)^(2/3) (Eq. 36); a = reflected fraction (after Matsumoto et al. 2012).
- Drift (mirror) reflection gives a speed of about 2 v_sh / cos(theta_Bn) upstream; reflected electrons
  can drive upstream whistlers if M_A >~ (cos(theta_Bn)/2)(m_i/m_e) (Eq. 31, Amano & Hoshino 2010); for
  quasi-parallel remnant shocks this needs v_A below a few km/s (Sect. 5.2).
- Ions (Caprioli & Spitkovsky 2014): about 10 % of post-shock energy above 10 times the flow energy;
  per cent level or less for theta_Bn >~ 50 deg (Fig. 13, Sect. 5.3).
- Quasi-parallel PIC: no significant electron acceleration when strictly parallel (Niemiec et al. 2012);
  theta_Bn = 5.7 deg at v_sh = 0.5 c gives short large-amplitude magnetic structures that scatter
  electrons (Dieckmann et al. 2010); 1-D run of 310 Omega_ci^-1 at 30 deg, M_A = 20 (Park et al. 2015).
- Reformation for beta_i < 0.5 at M_A ~ 5; critical Mach number 2.76 for a perpendicular shock (Sect. 7).
  Superluminal shock: v_sh tan(theta_Bn) >= c (Eq. 7).

## Assessment
- Numerical constraints stated: nearly all runs are 2D3V; essentially all use v_sh > 0.1 c, ten times a
  real remnant shock; reduced mass ratio (Sect. 5). Runs last a few Omega_ci^-1, so ions do not isotropize
  (Sect. 5.2). In 2-D only wave vectors in the plane exist: Buneman waves are complete only for an
  out-of-plane field, other instabilities need other orientations (Sect. 5.3). The ranking of processes
  can change with dimension and mass ratio; scans in mass ratio allow extrapolation to the real one;
  full 3-D needs exascale machines (Sect. 2.2). A reflecting wall may cause unphysical initial
  reflection at oblique and perpendicular shocks (Sect. 2.4).
- Listed as open or unsettled for electrons: surfing contributes little at the real mass ratio (Bohdan
  et al. 2019); reconnection in ion-Weibel current sheets as a second channel (Matsumoto et al. 2015);
  the first 3-D run (Matsumoto et al. 2017) accelerates faster than either 2-D orientation, by surfing
  followed by pitch-angle diffusion in Weibel turbulence; "only 3D simulations can provide a complete
  picture" (Sect. 5.3). Fully kinetic runs of oblique shocks are called for, because electrons may fill
  a foreshock where ions do not (Sect. 5.3).
- Disagreements reported: hybrid runs find no ion acceleration above about 50 deg, MHD-PIC runs (van
  Marle et al. 2018) find it with large-scale corrugation, and a Vlasov-Fokker-Planck study (Reville &
  Bell 2013) finds a disordered upstream field whatever the far-upstream direction (Sect. 5.3).
- Rippling: ion-scale, Alfvenic ripples (Lowe & Burgess 2003), origin disputed (ion temperature
  anisotropy or Kelvin-Helmholtz), seen by MMS (Johlander et al. 2016) (Sect. 7.2). Two opposite effects
  on electrons are quoted: more energization (Trotta & Burgess 2019; Sect. 5.3), and, at beta = 3 and low
  Mach number, less drift reflection because electrons pass through weak-field patches (Matsukiyo &
  Matsumoto 2015; Sect. 6).
- No simulation of a globally curved or spherical shock is discussed. Burgess (2006) and Guo & Giacalone
  (2010, 2015) are not in the reference list.
- Doubts: Sect. 1.2.3 ties the de Hoffmann-Teller frame to "if Eq. 7 applies" (reads inverted); Trotta &
  Burgess (2019) is called a study of "superluminal" shocks (that paper: supercritical); Eq. 31 is
  broken across lines in the text file and is given as printed.
