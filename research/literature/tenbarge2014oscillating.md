# TenBarge, Howes, Dorland & Hammett (2014): An oscillating Langevin antenna for driving plasma turbulence simulations

- Ref: Comput. Phys. Commun. 185(2), 578-589 (2014); DOI 10.1016/j.cpc.2013.10.022; arXiv:1305.2212;
  peer-reviewed.
- Read: 2026-11-13, skimmed (ar5iv HTML, Sects. 2-6); topics: shock-alfvenic-turbulence

## Claim
A driven, damped oscillator with random forcing (Langevin antenna), coupled to the plasma as an external
parallel current at the box-scale wavevectors, mimics the energy input from scales larger than the box:
it reaches steady Alfvenic turbulence within about one domain-scale turnaround time and can drive
strong or weak turbulence.

## Method
- Antenna amplitude per mode: a_{n+1} = a_n exp(-i omega_a dt) + F_a dt, F_a = sigma u_n (Eqs. 10-11),
  omega_a = omega0 + i gamma0 with gamma0 < 0; continuous form da/dt = -i omega0 a + gamma0 a + F_a (Eq. 15).
- u_n: complex, real and imaginary parts uniform on [-1/2, 1/2] (<|u|^2> = 1/6);
  sigma = A0 sqrt(12 |gamma0| / dt) gives <|a|^2> = A0^2 (Eqs. 12-13). Gaussian u needs another factor.
- Coupling: parallel vector potential A_za = a exp(i k.r) (Eq. 9), i.e. a body current
  J_z = -(c/4 pi) nabla^2 A_z, delta B_perp = -z x grad A_z (Eq. 8).
- Driven modes: four lowest, (k_x, k_y, k_z) = 2 pi (1,0,+-1)/(L_perp,L_perp,L_par) and (0,1,+-1)
  (Sect. 5.2); amplitudes independent (imbalance possible).
- Implemented in AstroGK (gyrokinetic, triply periodic). PIC is not used; PIC/hybrid appear only in
  the comparison of driving methods (Sect. 6), incl. an impedance mismatch at high driving frequency
  noted by Parashar et al. (2011).

## Numbers to compare against
- Recommended omega0 ~ 0.9 omega_l, gamma0 ~ -0.8 omega_l (omega_l = k_par v_A); Sect. 5.3 run:
  omega0 = 1.14 k_zD v_A, gamma0 = -0.9 k_zD v_A.
- Amplitude from critical balance: A0 = omega_l B0 k_par sqrt(delta) / (C2 k_perp^2 sqrt(N)) (Eq. 22),
  delta = 2, C2 = 1 typical, N = number of driven modes.
- Steady state across the domain within ~1 tau_D = 2 pi / omega_D (Fig. 8).

## Assessment
- The standard driver used later in PIC turbulence (zhdankin2017kinetic uses it with J_ext). For our
  null test it offers driven, statistically steady turbulence, whose heating rate is set by the antenna,
  in place of decaying turbulence whose amplitude falls during the run.
- For a shock run it is less natural: the antenna is a volume source in a periodic box; driving inside the
  upstream of a shock run would inject energy that a real upstream does not have.
- Limit: designed for delta B << B0 (gyrokinetic ordering); delta B ~ B0 needs care, as in the PIC
  papers that adopt it.
