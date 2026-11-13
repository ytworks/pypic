# Fox, Matteucci, Moissard, Schaeffer, Bhattacharjee, Germaschewski & Hu (2018): Kinetic simulation of magnetic field generation and collisionless shock formation in expanding laboratory plasmas

- Ref: Phys. Plasmas 25(10), 102106 (2018); DOI 10.1063/1.5050813; arXiv:1712.00152 (read); peer-reviewed.
- Read: 2026-11-13, skimmed (abstract, Sect. IV on the shock runs, Sect. V); topics: spherical-shock-injection

## Claim
Full-PIC simulations with a laser-ablation model reproduce laboratory plasma expansions. For a plume that
expands into a magnetised background, a curved magnetised shock forms in about one upstream ion
gyro-time, leaves the piston, and carries Biermann fields in its front.

## Method
- 2-D full PIC (PSC) with collisions and a heating operator that keeps ablating the target.
- Shock runs (Table 4): R_H/d_i0 = 30 (heating radius), B_up/B_0 = 0.04, n_bg/n_ab = 0.01,
  M_i/m_e = 100, one ion species. B_up/B_0 was "slightly increased ... to obtain faster shock
  formation in a given simulation wall time".
- The front is curved: the plume expands from a spot into a background with the field in the plane.

## Numbers to compare against
- Table 4 results: t_shock omega_ci,up ~ 1; V_shock/C_s0 = 4.5.
- Fig. 6: formation at t omega_ci,up = 0.9; the shock "completely separates from the piston" at 1.9.
  The field is compressed to about 6 B_up by 1.9.
- "these ablation-driven shocks continue to expand through the whole simulation" (Sect. IV.2).
- Biermann fields in the shock reach about 50 % of the compressed shock field (Sect. V).
- Derived here, not stated in the parts read: v_A,up/C_s0 = 0.04/sqrt(0.01) = 0.4, so V_shock/v_A,up
  is about 11.
- The box size in upstream d_i and the Mach number of the shock are not stated in the parts read.

## Assessment
- For this topic: the closest full-PIC predecessor with a curved, strong, piston-driven shock. It is 2-D,
  mass ratio 100, and the piston is fed by ablation during the whole run.
- It confirms in 2-D the times of Schaeffer 2020: onset near 1 and separation near 2 gyro-times.
- No electron acceleration is analysed in the parts read. The variation of theta_Bn along the curved
  front is not used.
- Left open: a finite debris mass, three dimensions, lower mass ratio.
