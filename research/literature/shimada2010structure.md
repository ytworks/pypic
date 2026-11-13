# Shimada, Hoshino & Amano (2010): Structure of a strong supernova shock wave and rapid electron acceleration confined in its transition region

- Ref: Phys. Plasmas 17, 032902 (2010); DOI 10.1063/1.3322828; arXiv:1004.2771 (v1); peer-reviewed.
- Read: 2026-11-13, skimmed (abstract, set-up, conclusions); topics: spherical-shock-injection
- Source: the arXiv full text (HTML of v1). Read: abstract, Sect. II (set-up), Sect. III (results, by figure)
  and the summary with its Eq. 2. The derivations of the appendix were not followed.

## Claim
Inside the shock transition region (STR) of a perpendicular shock with a very high Mach number, ion-driven
instabilities grow into two kinds of large field structures of the ion inertial scale: a leading wave packet
and a series of magnetic solitary humps. Electrons scatter many times between these structures and gain
energy directly. Electrons take almost half of the incoming ion energy within the STR. The hot electrons are
not Maxwellian and "could be the seed of further non-thermal population" (abstract).

## Method
- 1-D relativistic electromagnetic PIC. Plasma is injected from one side and reflected at a wall; the shock
  is exactly perpendicular (Sect. II).
- u_0/c = 0.096, beta_e = beta_i = 0.1, M_A = 174, M_s = 159, omega_pe/Omega_ce = 120, M/m = 100;
  150 particles per cell for each species; the electron plasma period is 115 time steps (Sect. II).
- Shock speed u_s = 1.5 u_0. Unit of length L_STR = u_0/Omega_ci = 115.2 c/omega_pi (Sect. III).

## Numbers to compare against
- Electron heating: "Almost half of the energy initially carried by the incoming ion" goes to the electrons
  within < 0.1 L_STR. The field energies stay near 1 % (B_z) and 0.1 % (E_x) of the total (Fig. 5).
- Time-averaged cross-shock potential: 0.3-0.4 M u_s^2 / 2e (Sect. III, text on Fig. 4).
- Spectrum: "we do not find a clear non-thermal tail, some excess is confirmed in higher energy region
  compared to the Maxwellian" (Fig. 8).
- Trajectories: gamma ~ 3 by bouncing between the leading packet and the humps, which the authors call a
  "1st-order Fermi type process", and gamma ~ 4 by stochastic bouncing among humps (Figs. 10, 11).
- Maximum energy from STR width = electron Larmor radius (Eq. 2): gamma_max ~ 5 against ~4 measured; with
  M/m = 1836, gamma_max up to 21, about 10 MeV.
- The Buneman instability "dose not have strong effect on the electron energization here" (sic, Sect. III).

## Assessment
- For this topic: a reference for the perpendicular end (the equator) at very high Mach number. The
  energization is by nonlinear ion-scale structures in the transition region, not by the Buneman instability.
  The energetic electrons stay inside the STR because the field is perpendicular. On an expanding sphere a
  field line that touches the front at the equator is cut at oblique angles later (our remark); the paper
  cannot say what follows from that.
- Limits the authors state: 1-D, reduced mass ratio, and "1D simulation tends to overheat electrons compared
  to the two-dimensional simulations". The result T_e >= T_i is against observations of remnants. See
  `ohira2008oblique` for the two-dimensional argument against such heating.
- M_A = 174 is far above what a debris-driven front in a 3-D box reaches. Use it for the list of mechanisms
  to test for, not for numbers.
