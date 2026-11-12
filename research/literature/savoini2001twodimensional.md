# Savoini & Lembège (2001): Two-dimensional simulations of a curved shock: Self-consistent formation of the electron foreshock

- Ref: J. Geophys. Res. 106(A7), 12975–12992 (2001); DOI 10.1029/2001JA900007; HAL hal-01646986 (open
  copy read); peer-reviewed, the venue's core scope (bow-shock physics).
- Read: 2026-11-13, full; topics: spherical-shock-injection

## Claim
A 2-D full-particle simulation of a shock curved on the large scale forms the electron foreshock
self-consistently. Most backstreaming electrons are mirror-reflected (loss-cone signature); about a
quarter have a high parallel energy and no loss cone, which mirror reflection alone does not explain.
- Question (Sect. 1): can the electron foreshock be reproduced with a self-consistent curved shock, without
  the assumptions of the fast-Fermi theory (stationary, plane, no scattering, prescribed fields)?
- Idea: drive a curved shock with a cylindrical magnetic piston in a 2-D full PIC code, so that
  theta_Bn varies continuously along one front.
- Prior work named: Leroy & Mangeney 1984 and Wu 1984 (theory); Krauss-Varban et al. 1989, Krauss-Varban &
  Burgess 1991 (test particles); Vandas 1995a,b (analytic, curved). The only earlier self-consistent curved
  shock named is Thomas & Winske 1990 (2-D hybrid, ion foreshock).

## Method
- 2-1/2-D electromagnetic relativistic full PIC; kinetic ions and electrons (Sect. 2).
- Box 768 x 1024 cells = 40 x 53 c/omega_pi; 8,388,608 particles; run 1.1 tau_ci.
- Curved shock from a cylindrical magnetic piston. B_0 partly out of the plane (theta_0 = 65 deg,
  phi_0 = 45 deg), so theta_Bn runs from 90 deg to 65 deg along the front (Eq. 2).
- m_i/m_e = 42; c = 3, v_the = 0.3 (code units), T_e/T_i = 1.58; beta_e = 0.24, beta_i = 0.15;
  v_A = 0.23; M_A = 3 at theta_Bn = 90 deg. Shock speed / c: about M_A v_A / c = 0.23 (derived, not printed).
- Radius of curvature R_c >= 34 rho_ci after t <= 0.2 tau_ci (Sect. 2; the text adds "about 200" in a
  second unit whose symbol the scanned copy does not render).
- Selection of backstreaming electrons (Sect. 3.2): upstream at the end of the run, inside the area swept by
  the shock, kinetic energy >= 6 x the initial thermal energy (E_k >= 0.8).

## Numbers to compare against
- Backstreaming electrons: about 1 % of the incident ones (Sect. 3.3, Fig. 4a). Their density falls
  steeply toward theta_Bn = 90 deg; none are found at the leading edge of the foreshock (the authors
  say this may be a limit of particle statistics).
- Time-dependent geometry (Sect. 1 and Sect. 3.4, Fig. 5b): "In the solar wind frame (i.e., frame of the
  present simulation) ... a given magnetic field line connected from the expanding curved shock to the
  observer point S will scan different angles theta_Bn in time"; Sect. 1: reflected electrons "will
  experience a time-dependent angle theta_Bn which can modify the effectiveness of the acceleration process".
- Theory quoted (Eq. 1 and text): mean parallel energy of reflected electrons about 2 m_e V_up^2 / cos^2
  theta_Bn; reflected density falls as exp[-1/cos^2 theta_Bn].
- Local distributions (Sect. 4.1): bump-in-tail; two peaks, P1 at zero drift and P2 with a drift up to
  <p_par> about 1.5 at theta_Bn about 71–65 deg. Density, thermal spread and drift of the backstreaming
  electrons all increase as theta_Bn moves away from 90 deg.
- Two populations (Sect. 4.2): low parallel energy, about 74 %, loss cone / ring (Fig. 9) = mirror
  reflection; high parallel energy (E_k,par > 0.8), about 26 %, no loss cone, Gaussian in v_perp.
- Leakage from downstream is argued against: T_perp,d / T_perp,u = 0.9 while B_d / B_u = 2.7 (Sect. 4.2).
- Self-reformation of the front is seen near 90 deg and disappears below about 85 deg.

## Assessment
- Shown: a self-consistent electron foreshock at a curved front, with the angular dependence of the
  backstreaming population between 90 and 65 deg, at m_i/m_e = 42 and M_A = 3.
- Stated: the first 2-D electromagnetic full-particle simulation of a curved shock.
- Stated limits: no particle trajectories (left for further study); too few particles per cell to resolve
  the leading edge; unrealistic mass ratio.
- Open (Sect. 4.2.1, Sect. 5): the origin of the high-parallel-energy population. Three hypotheses are
  listed (non-mirror interaction with the nonstationary front; mirror reflection plus scattering; leakage);
  the share of acceleration by the macroscopic fields against wave scattering is "not yet determined".
- For this topic: the statement that a field line connected to an expanding curved shock scans theta_Bn
  in time is printed here (2001). The paper does not relate the electron energy to the theta_Bn history
  of one electron; that is taken up in Savoini et al. 2010.
- Doubts: with rho_ci = 28 cells (Table 1 of the 2010 paper, same parameters) 34 rho_ci is about 950
  cells, more than the 768 x 1024 box holds as a radius early in the run; the value of R_c in ion scales
  should be re-derived from the figures before it is used. Low M_A and a hot, slow-light plasma
  (c / v_the = 10).
