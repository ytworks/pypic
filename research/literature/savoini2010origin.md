# Savoini, Lembège & Stienlet (2010): Origin of backstreaming electrons within the quasi-perpendicular foreshock region: Two-dimensional self-consistent PIC simulation

- Ref: J. Geophys. Res. 115, A09104 (2010); DOI 10.1029/2010JA015263; HAL hal-00475416 (open copy read);
  peer-reviewed, the venue's core scope.
- Read: 2026-11-13, full; topics: spherical-shock-injection

## Claim
Trajectories of self-consistent electrons at a curved full-PIC shock show three classes of backstreaming
electrons: mirrored (about 70 %), trapped in the parallel electrostatic well of the overshoot (about
20 %), and leaked from downstream (about 10 %). A leaked electron can enter the front near
theta_Bn = 88 deg and leave it near 48 deg.
- Question (Sect. 1): what produces the high-parallel-energy population of Savoini & Lembège 2001, and do
  leaked electrons exist?
- Idea: follow the time trajectories of electrons of the PIC run itself (not test particles) across a
  front whose theta_Bn runs from 90 to 45 deg.

## Method
- Same 2-1/2-D electromagnetic full PIC and cylindrical magnetic piston as the 2001 paper (Sect. 2, Table 1).
- Box L_x = 95, L_y = 110 c/omega_pi; theta_Bn from 90 to 45 deg; run 2.2 tau_ci.
- m_i/m_e = 42; c = 3; T_e/T_i = 1.58; v_A = 0.23; M_A = 3 at 90 deg. Table 1 (electrons / ions):
  v_th 0.3 / 0.04, rho_c 0.84 / 28, c/omega_p 3 / 19, omega_c 0.5 / 0.012, omega_p 1 / 0.155,
  tau_c 13 / 524, beta 0.16 / 0.10. Shock speed / c about 0.23 (derived).
- R_c >= 34 rho_ci after t <= 0.2 tau_ci (paragraph [6]); the text also prints "660 rho_ci" for the end
  of the run (paragraph [8]).
- Selection: E_end / E_init >= 2. 1500 sampled trajectories.

## Numbers to compare against
- Backstreaming electrons: 4.7 % of the electrons that interact with the shock (Sect. 3.2, Fig. 2).
- Mirrored (Sect. 3.2.1, Fig. 3): interaction time about 2.7 tau_ce; parallel energy up to about 7 (gain
  >= 30); mu conserved; some short multibounces. Penetration about rho_ce.
- Trapped (Sect. 3.2.2, Fig. 4): stays in the overshoot about 40 tau_ce; bounces in the parallel
  electrostatic well (period 12.5 predicted, about 11 measured); parallel energy about 4; mu not conserved
  (mu/mu_0 about 4). Penetration >= rho_ci.
- Leaked (Sect. 3.2.3, Fig. 5): stays downstream about 2 tau_ci; "the electron does not see the same shock
  orientation theta_Bn as it goes into and out from the shock front. The present particle 3 originally hits
  the shock front in a nearly perpendicular region (theta_Bn about 88 deg), where no reflected electrons are
  observed and finally succeed to escape from the downstream region in a much more oblique
  quasi-perpendicular domain (theta_Bn about 48 deg), where the electrostatic barrier has a much lower
  amplitude". p_par reaches about 2.1. The gain is attributed to multiple magnetic bounces in the
  downstream fluctuations and is shared between parallel and perpendicular.
- Shares (Sect. 4): mirrored about 70 %, trapped about 20 %, leaked about 10 %. "All 'leaked' electrons
  finally escape into the upstream region within an angular range below theta_Bn about 60 deg"; other
  leaked electrons "may originally penetrate the front ... within a wide range of oblique shock normal
  directions".
- Front structure: reformation disappears for theta_Bn <= 72 deg; whistler precursor for <= 60 deg; the
  shock amplitude and the parallel field E_l,par (maximum about 0.22) drop by about 2 near 53 deg, which
  the authors call a "funnel" for the escape. All three example electrons end near theta_Bn about 55 deg.
- The time-dependent theta_Bn of a field line connected to the expanding shock is repeated (paragraph [22]).

## Assessment
- Shown: by trajectories, an electron path that enters a curved front close to 90 deg and exits at a much
  lower theta_Bn, in a full-PIC run with kinetic electrons.
- Stated limits: unrealistic m_i/m_e and omega_ce/omega_pe; the results are called preliminary; the
  statistical study is left for a later paper (not found among the later papers of these authors, which
  treat the ion foreshock).
- Open: how the shares depend on M_A, mass ratio and curvature; the energy of leaked electrons against
  the entry and exit angles.
- For this topic: the cross-latitude path is published here, but only for the leaked class (10 % of 4.7 %),
  with a downstream mechanism. The paper does not compare the final energy with the single-reflection
  energy V_sh / cos theta_Bn at the exit point, and does not treat drift along the front in the ramp as the
  transport. Shock deceleration is not discussed for electrons.
- Doubts: "660 rho_ci" for the final R_c does not fit a box of about 100 c/omega_pi; read it as a
  misprint until checked. Three example trajectories carry the argument; the 70/20/10 split comes from
  1500 sampled electrons with class criteria stated only by penetration depth.
