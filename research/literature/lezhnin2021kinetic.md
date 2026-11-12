# Lezhnin, Fox, Schaeffer, Spitkovsky, Matteucci, Bhattacharjee & Germaschewski (2021): Kinetic simulations of electron pre-energization by magnetized collisionless shocks in expanding laboratory plasmas

- Ref: Astrophys. J. Lett. 908, L52 (2021); DOI 10.3847/2041-8213/abe407; arXiv:2001.04945 (read);
  peer-reviewed.
- Read: 2026-11-13, skimmed (abstract, set-up, the passage on curvature, conclusions);
  topics: spherical-shock-injection

## Claim
Full-PIC simulations of a laser-ablation piston driving a magnetized shock show electron pre-energization
by shock drift acceleration under laboratory conditions; the curvature of the front is argued to be
unimportant.

## Method
- 2-D full PIC with Coulomb collisions; box 40000 x 40 cells = 30000 x 30 d_e,ab, i.e. quasi-1-D with a
  planar front (z along the shock normal, y out of the plane). B_0 = B_0 (sin theta_Bn e_y + cos theta_Bn
  e_z), with theta_Bn from 50 to 90 deg as a run parameter: one angle per run.
- m_i/m_e = 100 (1-D check with 400); T_e,ab/m_e c^2 = 0.04; beta_e,up = 0.5–2; n_up/n_ab = 0.05; M_A about
  15 in the reference run.

## Numbers to compare against
- Regime stated for efficient SDA: 33 >= M_A >= 7 and beta_e,up about 1; "only quasi-perpendicular shocks
  with theta_Bn < 70 deg show significant acceleration" (as printed, summary of requirements).
- Out-of-plane runs reflect at least an order of magnitude more electrons than in-plane runs (auxiliary
  runs; 3-D left open).
- Curvature passage (verbatim): "we do not expect a significant influence of the shock curvature on the SDA
  acceleration of electrons. In principle, the shock front curvature does affect the efficiency of the SDA,
  since it requires a significant transverse motion of the pre-accelerated particle. But in our case, at
  later stages of shock evolution, the radius of curvature scales as M_A Omega_i t and, after a few ion
  gyrotimes, will be significantly larger than transverse acceleration scales within the shock. Thus, the
  SDA mechanism will not be affected."

## Assessment
- This is not a curved-shock simulation: each run has one theta_Bn and a planar front.
- For this topic: it states the opposite expectation for a laser-driven blast (curvature negligible after a
  few ion gyrotimes) as an argument, not as a result. The argument compares R_c with the drift length in
  the ramp; it does not consider the change of theta_Bn along a field line in time.
