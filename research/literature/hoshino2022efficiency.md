# Hoshino (2022): Efficiency of Nonthermal Particle Acceleration in Magnetic Reconnection

- Ref: Phys. Plasmas 29, 042902 (2022); doi:10.1063/5.0086316; arXiv:2203.15169
- Read: 2026-11-12, full (arXiv v2); topics: (none yet)

## Claim
In 2D pair-plasma reconnection, the share of exhaust energy held by nonthermal particles rises with
plasma temperature (magnetization) in the nonrelativistic regime and saturates in the relativistic
regime. There the exhaust spectrum is a Maxwellian plus a kappa distribution with kappa ~ 3
(power-law index ~ 2), and the nonthermal part carries about 95 % of the exhaust internal energy.

## Method
- 2D PIC; Harris pair-plasma sheet, uniform T, background Nb = 0.05 N0, no guide field; periodic x,
  conducting y walls. Fiducial box 430 x 430 lambda on 5376^2 cells with 8.4e10 particles; seed
  perturbation max(By/B0) = 1e-2 at the centre.
- Scan: T/mc^2 = 10^(n/2-2), n = 0..7 (0.01 to ~32), so sigma = B0^2/(8 pi N0 mc^2) = 2 T/mc^2;
  gyroradius-to-thickness ratio rg/lambda = 0.45 held fixed.
- Spectra integrated over the exhaust inside the separatrix (the last reconnected Az contour, at least
  lambda wide) at the late, nearly steady stage (one big plasmoid).
- Fit: N = N_M + N_kappa with an exponential cutoff above gamma_cut (Eqs. 4-6); residual (Eq. 7).
- Measures: fractions of number and energy in the kappa part, xi_den and xi_ene. The nonthermal
  fractions eps_den and eps_ene subtract the kappa -> inf (Maxwellian) limit of the fitted kappa
  (Sec. IV, Eq. 8).
- Exhaust cut into equal-dAz flux tubes (Fig. 5): spectra per tube (Fig. 6) and for one tube in
  time (Fig. 7). App. A: eps of a pure kappa distribution against kappa and T (Fig. 8).

## Numbers to compare against
- Residual < 1e-1 over wide energy ranges (Fig. 2, Eq. 7); fit parameters vary by about 10 % between
  late snapshots (Sec. III).
- kappa large (near Maxwellian) at T/mc^2 = 0.01; kappa ~ 3 (s = kappa - 1 ~ 2) for T/mc^2 > 1
  (Fig. 3b).
- gamma_cut / (T_kappa/mc^2) ~ 50-100, except at T/mc^2 = 0.01 (Fig. 3c).
- xi_ene > 95 %, xi_den > 70 % (Fig. 4a).
- eps_ene ~ 95 % and eps_den saturating at ~ 60 % for T/mc^2 > 1; both rise with T below
  T/mc^2 = 1 (Fig. 4b). The summary gives eps_den ~ 40-60 % (Sec. VI).
- Fixed-tube final fits (T_M/mc^2, T_kappa/mc^2, kappa, gamma_cut) = (0.013, 0.078, 6.1, 620),
  (0.79, 0.44, 3.0, 22), (15, 9.0, 2.6, 470) for T/mc^2 = 0.01, 1/sqrt(10), 10 (Fig. 7, Sec. V).
- Highest energies several 1e2-1e3 x the initial thermal energy, gyroradius ~ 100-200 lambda, about
  the island size (Sec. III).
- Box test at T/mc^2 = 0.1 (boxes 215^2, 860 x 215, 1720 x 430 lambda^2): T_M, T_kappa, kappa
  change by < 15 %; gamma_cut ~ 6.7 in the smallest box, 9-12 in the others (Sec. VI).
- Background Nb/N0 = 1/80 (preliminary): kappa = 2.2, 2.6, 2.7, 4.3 for T/mc^2 = 10, 1, 0.1, 0.01,
  about 30 % below the reference runs (Sec. VI).
- Pure kappa, T -> inf: kappa = 3 gives eps_den -> 7/9, eps_ene -> 1; kappa = 4 gives 5/8 and
  29/32 (App. A, Fig. 8).

## Assessment
- Shown by simulation and fitting: the trend of the nonthermal fraction with T/mc^2 over eight runs.
  It also shows where energy is gained: a jump when tubes reach the X-line and separatrix, then a
  milder gain downstream (Figs. 6-7).
- The fractions depend on the model. They come from the Maxwellian + kappa decomposition with the
  thermal part defined as the kappa -> inf limit. T_M stays near the initial T (unheated inflow
  inside the separatrix). The kappa form is "not necessarily a good approximation", only better than
  a Maxwellian. Another cutoff shape fits better at T/mc^2 = 10.
- Inconsistency: Sec. V quotes kappa < 2 for relativistic runs, while Fig. 3b and the summary give
  kappa ~ 3 (s ~ 2), probably s confused with kappa. The individual mechanisms are only reviewed.
- Stated, with separate papers announced: a cold background (Tb/mc^2 = 1e-4) gives almost the same
  efficiency, and a lower background density hardens kappa.
- Limits (stated): 2D, no guide field, pair plasma (ion-electron and ion acceleration need longer
  runs and larger boxes). The finite box and run time set gamma_cut.
