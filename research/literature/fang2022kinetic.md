# Fang, Xia, Tian, Zhou & Yu (2022): Kinetic simulation of electron, proton and helium acceleration in a non-relativistic quasi-parallel shock

- Ref: MNRAS; DOI 10.1093/mnras/stac886; arXiv:2202.02689 (volume and page not recorded). Peer reviewed.
  Read in the arXiv version.
- Read: 2026-11-14, abstract, set-up and the main results (arXiv PDF text); topics:
  qpar-electron-conditions

## Claim
A 1-D PIC simulation of a quasi-parallel shock with protons, helium ions and electrons accelerates all
three species; the electrons form a power law over part of the momentum range.

## Method
- 1D3V PIC (Smilei); helium-to-proton number ratio 0.1.
- Shock speed 0.27c, theta_Bn 30 deg, mass ratio 30, M_s 14, M_A 19.5.

## Numbers to compare against
- Electron power law over ~3-20 m_e c at t = 2.16e5 omega_pe^-1.
- ~7 % and ~5.4 % of the bulk kinetic energy in non-thermal protons and helium.

## Assessment
- One condition with mass ratio 30 and v_sh 0.27c; no diffusion coefficient, no injection rate.
- For this topic: with fang2019early (on `research/spherical-shock-injection`) and yu2024nonthermal,
  the 1-D quasi-parallel runs of this group at M_A ~17-20. Consistent with tails at M_A >~ 10
  (gupta2024electron); not a prior work for the measurement proposed here.
