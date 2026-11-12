# Matsumoto, Amano & Hoshino (2013): Electron acceleration in a nonrelativistic shock with very high Alfvén Mach number

- Ref: Phys. Rev. Lett. 111, 215003 (2013) (volume and page from the reference list of Matsumoto et al. 2017; the
  arXiv metadata gives only the DOI); DOI 10.1103/PhysRevLett.111.215003; arXiv:1310.8367 (v1); peer-reviewed (PRL)
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
A 2-D perpendicular shock at M_A ~ 45 with M/m = 225, chosen to satisfy M_A >~ (M/m)^(2/3), shows three
instabilities at once in different places: Buneman waves at the leading edge of the foot, ion-beam Weibel
filaments throughout the foot, and ion-acoustic-like drift waves at the front. Electrons are heated and
accelerated mainly by surfing on the Buneman waves; the resulting tail exceeds the upstream ion kinetic energy.

## Method
- 2-D PIC with second-order shape function and Esirkepov charge-conserving deposition (code not named);
  injection method, downstream rest frame, periodic in y.
- Purely perpendicular (theta_Bn = 90 deg), B_0 along z (out of plane).
- M/m = 225; M_A = 44.8 in the shock rest frame; beta_e = 0.5; omega_pe/Omega_ge = 10; T_i/T_e = 1. Sonic Mach
  number not stated. V_0/c not stated (0.2 follows from the quoted "5 B_0 (25 E_0)").
- L_y = 6.8 c/omega_pi, L_x grows with the shock; 24000 x 2048 cells of one Debye length; omega_pe dt = 0.025;
  40 particles per cell per species; ~1e10 particles late in the run.
- Results at Omega_gi T = 4; the spectrum is said not to change much afterwards.

## Numbers to compare against
- Condition for effective SSA: M_A >~ (M/m)^(2/3) (Eq. 1, from Matsumoto et al. 2012); 225^(2/3) = 37 (my arithmetic).
- Buneman wave amplitude about twice the local field (~B_0) at the leading edge; oblique, coherent, electron
  scale (Fig. 1b, 2a); this meets the "unlimited SSA" condition of Hoshino & Shimada (2002).
- SSA: electrons reach gamma ~ 4 within about an electron gyro-period (Fig. 2a); a nonthermal distribution exists
  in the foot (50.0-53.3 c/omega_pi) before the front (Fig. 2b).
- Downstream (43.3-45.0 c/omega_pi): maximum gamma ~ 12, i.e. 2.5 times the upstream ion kinetic energy; fit by two
  relativistic Maxwellians, T_1/mc^2 = 0.27 (main) and T_2/mc^2 = 0.58 ("nonthermal") (Fig. 2b).
- Weibel mode: wavelength about one ion inertial length (|k| c/omega_pi ~ 2 pi), wave vector tilted from x; B_z
  varies between 5 and 10 B_0 within an ion inertial length near x = 47, forming current sheets (Fig. 1c, 3a).
- Density overshoot 25 N_0; downstream ~4 N_0; an m = 1 modulation along the front (Fig. 1a).
- Drift waves behind the overshoot: |E| ~ 5 B_0 = 25 E_0 at k_y lambda_e ~ 1; electron drift 0.1 c; electrons there
  already have T_e ~ mc^2 (Fig. 1d, 4).

## Assessment
- Shown by simulation: one run; coexistence and location of the three instabilities; SSA on traced electrons;
  after SSA the further gain is described as "basically adiabatic".
- The tail is characterised by a second Maxwellian, not a power law; no injected fraction is given.
- Conjectured: magnetic reconnection in the Weibel current sheets as an extra dissipation channel in 3-D at
  higher M_A; ion-cyclotron turbulence at the overshoot (excluded by the out-of-plane field) would break the
  coherent front and scatter the pre-accelerated electrons; SSA is a "robust" seed for DSA at young SNRs.
- Stated limit: the out-of-plane configuration restricts the possible instabilities.
- Limits I note: perpendicular only; M/m = 225; 4 Omega_gi^-1; V_0/c ~ 0.2; a single run, so Eq. 1 is confirmed at
  one point above the threshold and not bracketed from below at this mass ratio.
- Context: the paper calls the efficiency of SSA in several dimensions "controversial" (Ohira & Takahara 2007;
  Amano & Hoshino 2009; Riquelme & Spitkovsky 2011) and claims the first confirmation of Eq. 1 at a mass ratio
  large enough to separate ion and electron scales.
- The ion-acoustic heating at the front is called complementary only: electrons are heated chiefly at the leading
  edge, by the Buneman instability, rather than at the ramp.
