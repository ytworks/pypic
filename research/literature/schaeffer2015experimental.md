# Schaeffer, Everson, Bondarenko, Clark, Constantin, Winske, Gekelman & Niemann (2015): Experimental study of subcritical laboratory magnetized collisionless shocks using a laser-driven magnetic piston

- Ref: Phys. Plasmas 22(11), 113101 (2015); DOI 10.1063/1.4934983; read from OSTI (purl 1493543);
  peer-reviewed.
- Read: 2026-11-13, full; topics: spherical-shock-injection

## Claim
A laser-driven magnetic piston in the LAPD made a subcritical quasi-perpendicular shock in one of three
configurations. The three outcomes (no shock, pre-shock, shock) agree with a set of formation criteria:
a super-Alfvenic piston, enough coupling (epsilon), and enough magnetised room (mu).

## Method
- Experiment: Raptor laser on a carbon target, debris into the magnetised H or He plasma of the LAPD;
  magnetic flux probes along the blow-off axis. Three runs (Table I): B_0 = 275, 200, 300 G; ambient
  ion inertial length 20, 10, 7 cm; debris gyroradius rho_d = 35, 32, 31 cm.
- 2-D hybrid simulation of Run 3 (Sect. V): a C+5 debris cloud expands conically at M_A = 2.65 from a
  disk 0.25 c/omega_pi wide into a Gaussian ambient core.

## Numbers to compare against
- Conditions named in Sect. I: the piston must drive a super-Alfvenic pulse "long enough for
  instabilities to grow in the pulse front"; "the system size must be large enough for a shocked ambient
  ion moving at the shock speed to complete one gyro-orbit"; hybrid work predicts "the piston width must
  be larger than an ion inertial length and the piston radius must be larger than a piston ion directed
  gyroradius" (their Refs. 25, 26).
- Sect. II: critical Mach number M_c <= 2.76. Coupling parameter epsilon: "originally derived by
  Bashurin et al. using the equal charge radius"; the paper uses "an alternate, empirically motivated
  version using the equal mass radius", Eq. (1): epsilon = R_M^2/(rho_d rho_a)
  = (R_M/rho_d)^2 (m_d Z_a/(m_a Z_d)), with R_M = (3 N_d m_d/(4 pi m_a n_a))^(1/3) and rho_d, rho_a the
  directed gyroradii of debris and ambient ions. It "was derived strictly for large-Mach number
  shocks". "Hybrid simulations [Clark 2013] suggest that epsilon >~ 0.5 is necessary for shock
  formation."
- Magnetisation parameter mu = D_0/rho_a,s > 1 (Sect. II, Table II). rho_a,s is the gyroradius of a
  shocked ambient ion. "The system size D_0 is defined as the distance from the target to the edge of
  the (suitable) ambient plasma." The same condition as a time: D_0/v_s > 1/Omega_ci,s.
- Table II, Runs 1/2/3: M_A = 2.8/1.3/1.9; epsilon = 0.04/0.9/0.6; mu = 0.8/0.7/1.8; "only Run 3
  satisfies all formation conditions".
- Run 3 (Sect. IV): at 1.9/Omega_ci and x = 29 cm (0.9 rho_d) the pulse has B_z/B_0 = 2.0, M_A = 2.1,
  ramp 0.5 c/omega_pi. The front moves at M_A 2.1 and the back at 1.6: it begins to leave the piston.
- Run 2: B_z/B_0 = 1.7, a pre-shock that dissipates. Run 1: B_z/B_0 = 1.3, no shock.
- "the most difficult criterion to satisfy was sufficient magnetization mu. This in turn appears to be
  most sensitive ... to the amount of space available with suitable ambient conditions."
- Fig. 10 regimes: weak (epsilon < 0.5, mu < 1), moderate (epsilon > 0.5, mu < 1), strong (both above).
  "continuing to increase epsilon or mu will eventually cause M_A < 1."
- Summary: the shock "was beginning to separate from the piston before dissipating in a low-density
  ambient plasma".

## Assessment
- For this topic: this is the source of the two-parameter picture. Coupling (epsilon > 0.5) decides
  whether the debris hands its energy over. Room (mu > 1) decides whether the pulse steepens to a shock.
- Derived here, for equal mass and charge of debris and ambient ions: epsilon = (R_m/rho_d)^2. A design
  with R_m/rho_d = 6-8 has epsilon = 36-64. Coupling is then not the limit; room and the speed that is
  left at formation are.
- The experiment reached only M_A about 2 with pistons of M_A about 2-3. It does not show a strong
  shock from a finite debris cloud. The criteria come from hybrid runs and say nothing about electrons.
- Left open: the Mach number of the shock as a function of epsilon and of the debris M_A.
