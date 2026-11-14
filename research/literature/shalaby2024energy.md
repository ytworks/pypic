# Shalaby (2024): Energy Dissipation in Strong Collisionless Shocks: The Crucial Role of Ion-to-electron Scale Separation in Particle-in-cell Simulations

- Ref: ApJL; DOI 10.3847/2041-8213/ad99d8; arXiv:2412.03530 (volume and page not recorded). Peer reviewed.
  Read in the arXiv version.
- Read: 2026-11-14, abstract and the passages on the mass-ratio comparison (arXiv PDF text); topics:
  qpar-electron-conditions

## Claim
In 1-D parallel shocks a reduced ion-to-electron mass ratio changes electron acceleration in opposite
directions at high and low Alfven Mach number, and distorts the heating of both species.

## Method
- 1D3V PIC (SHARP code), parallel shocks; mass ratio 100 against 1836 at M_A 21.3 and at M_A 5.3.

## Numbers to compare against
- M_A 21.3: mass ratio 100 gives more efficient electron acceleration than 1836, and an unrealistically
  high ion flux at high momentum.
- M_A 5.3: mass ratio 100 suppresses electron acceleration completely; 1836 accelerates efficiently.
- Reduced mass ratio: too much electron heating, too little ion heating. With 1836 the ion-to-electron
  temperature ratio does not depend on the upstream field.

## Assessment
- For this topic this is the main warning on reduced mass ratio (null test N5 and runs V1, P6, S1-S3,
  which use mass ratio 100). The sign of the error depends on M_A, so a mass-ratio-100 result cannot be
  corrected by a fixed factor.
- Follows shalaby2022mechanism (on `research/spherical-shock-injection`): efficient acceleration at
  M_A 5.3 with the real mass ratio, the opposite of the low-M_A stall in gupta2024electron (theta_Bn
  30 deg, 2-D-like quasi-1-D). The two differ in geometry (0 against 30 deg), sonic Mach number and code.
- 1-D and exactly parallel; whether the statements hold at 8-30 deg is not shown.
