# Matsumoto, Amano & Hoshino (2012): Electron accelerations at high Mach number shocks: two-dimensional particle-in-cell simulations in various parameter regimes

- Ref: Astrophys. J. 755, 109 (2012) (volume and page from the reference lists of Matsumoto et al. 2013, 2017; the
  arXiv metadata gives only the DOI); DOI 10.1088/0004-637X/755/2/109; arXiv:1204.6312 (v1); peer-reviewed (ApJ)
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
Whether shock surfing acceleration (SSA) makes relativistic electrons in 2-D is decided by two conditions from
Buneman theory, the instability threshold (Eq. 5) and a stricter trapping condition M_A >= (1+alpha)(M/m)^(2/3)
(Eq. 8), together with a low upstream beta_e. Runs meeting all of them (A, C) show nonthermal tails; runs failing one
(B, D) do not. The authors present this as the reason why earlier 2-D studies reached opposite conclusions.

## Method
- 2-D electromagnetic PIC (MPI + OpenMP; code not named), injection method, downstream rest frame.
- Purely perpendicular shock (theta_Bn = 90 deg), B_0 along z (out of plane).
- omega_pe/Omega_ge = 10; T_i/T_e = 1; 40 particles per cell per species; cell = Debye length; omega_pe dt = 0.025.
- Table 1 (M_A in the shock rest frame): A: M/m = 25, M_A = 14.4, beta_e = 0.5 (same as Amano & Hoshino 2009);
  B: 100, 16.2, 0.5; C: 100, 30.0, 0.5; D: 100, 30.0, 4.5. V_0/c = 0.2 in run C (Sect. 4.1.1). Sonic Mach
  numbers not stated.
- Box L_x = 12 V_0/Omega_gi, L_y = 5 c/omega_pi; results at Omega_gi T = 6. Run C: 24001 x 1024 cells, 5e9 particles.

## Numbers to compare against
- Eq. 5: M_A >= [(1+alpha)/2] sqrt(beta_e M/m), alpha = n_r/n_0 ~ 0.2 (0.22 measured, Sect. 4.1.1).
- Eq. 6-8: trapping needs E_BI/B_0 >= 1; with saturation energy C m n_0 dV^2/2, C ~ (m/M)^(1/3), cut to ~25 % in
  2-D (Amano & Hoshino 2009b), this gives M_A >= (1+alpha)(M/m)^(2/3): 10.3 for M/m = 25, 25.8 for 100 (Table 1).
  For 1836 it gives ~180 (my arithmetic; the paper shows it only as the line in Fig. 10).
- Electrostatic field at the leading edge, |E_est|/B_0: ~2.4 (A), ~0.5 (B), ~2 (C), ~0.5 (D) (Sect. 4.1-4.3, Fig. 8).
  Run C has coherent oblique wave trains, run A isolated peaks. Run C: dV ~ 0.4 c against v_the = 0.07 c.
- Run C, stage 1 (SSA): gamma ~ 1 to ~3 in ~20 Omega_ge^-1 = 0.2 Omega_gi^-1, typically three encounters with the
  Buneman layer, each made possible by the larger gyro-radius (Fig. 4-5). Foot spectrum reaches gamma ~ 3.5.
- Run C, stage 2 near the overshoot: adiabatic gain by gradient-B drift (SDA) plus betatron acceleration from the
  inductive field of reformation (~0.3 B_0), gamma ~ 3 to ~7 (Sect. 4.1.2, 4.2); downstream tail to gamma ~ 9 above
  a thermal part with T_e/mc^2 = 0.13 (Fig. 3). Overshoot field 20 B_0 (C) against 8 B_0 (A) (Fig. 7).
- Downstream T_e/mc^2: 0.062 (A), 0.033 (B), 0.13 (C), 0.11 (D) (Sect. 4.4).
- Scalings: B_ov ~ 0.4 B_0 M_A^(7/6); gamma_ov = sqrt(0.2) M_A^(7/12) gamma_ft (Eq. 12-13).

## Assessment
- Shown by simulation: four runs whose outcomes follow Eq. 5 and Eq. 8 at beta_e = 0.5; multi-encounter SSA on
  one traced electron; no spectral index and no nonthermal fraction is quoted.
- Not explained: run D (beta_e = 4.5) satisfies both conditions yet has weak waves and almost no tail; the authors
  say the cold-plasma estimates fail at high beta_e and place Kato & Takabe (2010) in the same class.
- Extrapolated, not simulated: SSA works at SNR shocks "provided" the upstream electrons are cold enough (Fig. 10);
  the adiabatic stage could add an order of magnitude at SNR Mach numbers (Eq. 13), but that energy is given
  back downstream unless electrons are scattered at the front (Sect. 5).
- Stated limits: strictly perpendicular; out-of-plane field, so no waves along B, no ripples; drift waves along the
  front were seen but not linked to acceleration; a 2-D run holds either 2-D surfing (out of plane) or rippling
  (in plane, where surfing would be 1-D-like and electrons could not leave the trap); 3-D needed (Sect. 5).
- Limits I note: M/m <= 100; V_0/c = 0.2; one run per parameter set; trapping condition tested at two mass ratios.
- Disagreements: Riquelme & Spitkovsky (2011) concluded SSA fails at realistic mass ratio from M_A = 7 runs that
  do not meet Eq. 8; Kato & Takabe (2010) blamed the in-plane field, here the high beta_e is blamed. The paper
  quotes beta_e = 26 for Kato & Takabe, whose own text gives beta ~ 26 (total, as far as I can tell).
