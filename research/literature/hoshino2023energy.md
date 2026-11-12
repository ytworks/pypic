# Hoshino (2023): Energy Partition of Thermal and Nonthermal Particles in Magnetic Reconnection

- Ref: Astrophys. J. 946, 77 (2023); doi:10.3847/1538-4357/acbfb5; arXiv:2302.13517
- Read: 2026-11-12, full (arXiv v1); topics: (none yet)

## Claim
This paper adds a guide field BG to hoshino2022efficiency. In 2D pair-plasma reconnection the
nonthermal energy share exceeds 90 % for relativistic, weak-BG cases and falls as BG grows.
Nonrelativistic anti-parallel reconnection mostly heats. A moderate weak guide field
(BG/B0 ~ 1/4 at T0/mc^2 = 0.01) raises the nonthermal output. The trends follow the downstream
(non-ideal) reconnection electric field.

## Method
- Same 2D set-up as hoshino2022efficiency (430 x 430 lambda, 5376^2 cells, 8.4e10 particles,
  ~2.9e3 per cell; periodic x, conducting y; Nb = 0.05 N0; rg/lambda = 0.45), plus a uniform BG.
- Run grid: BG/B0 = 0, 1/8, 1/4, 1/2, 1 times T0/mc^2 = 0.01, 0.1, 1, 10 (sigma = 0.02-20)
  (Eqs. 3-4).
- Exhaust spectra (inside the separatrix) fitted by chi-square to N_M + N_kappa with an exponential
  cutoff (Eqs. 5-9); residual (Eq. 10).
- Measures: mean thermal temperature T_th (Eq. 11). eps_ene (Eq. 12) has the preheated plasma in its
  denominator. eps'_ene is relative to the heating increment plus the nonthermal energy (Eqs. 13-14),
  with the adiabatic index switched from 5/3 to 4/3 between regimes.
- Theory: diffusion-region acceleration against escape gives N(eps) (Eqs. 15-16); collisionless
  tearing growth rate with a guide field (Eq. 17).
- Field diagnostic: time maximum of the downstream integrals of |E_z| and |E'_z|, with
  E' = E + v x B/c (Fig. 5).

## Numbers to compare against
- Residual < 1e-1 over wide ranges for all 20 runs (Fig. 2, Eq. 10).
- T_M/T0 = 0.8-1.6 (Fig. 3b); T_kappa/T0 = 7.6 at T0/mc^2 = 0.01, BG = 0 (Fig. 3c); T_kappa is
  several times T_M for T0/mc^2 < 0.1 and BG/B0 < 0.5 (Sec. 4).
- kappa falls with T0 and rises with BG, except an anomalously low kappa at T0/mc^2 = 0.01,
  BG/B0 = 1/4 (Fig. 4a).
- eps_ene > 90 % for relativistic BG = 0 and ~30 % at T0/mc^2 = 0.01, BG = 0 (Fig. 4b); eps_ene
  and eps'_ene nearly equal for T0/mc^2 > 1, BG ~ 0 (Fig. 4b-c).
- Model: N ~ exp(-(Bn/E)(v/c)) near eps ~ mc^2 and N ~ (eps/mc^2)^(-Bn/E) for eps >> mc^2
  (Eqs. 15-16). In the relativistic case E ~ (VA/c) Bn, so the index ~ c/VA ~ 1 (Sec. 5).
- Tearing growth gamma tau_A ~ (k lambda)(rg/lambda)^2 B0/BG for a strong guide field (Eq. 17).
- |E_z|_S and |E'_z|_S rise with T0 and fall with BG; |E'_z|_S has a bump near T0/mc^2 = 0.01,
  BG/B0 = 0.25 (Fig. 5c-d).
- Observed flare and magnetotail power-law indices of about 4-6 are cited for context (Sec. 5).

## Assessment
- Shown by simulation: the (T0, BG) maps of kappa, eps_ene and eps'_ene, and their resemblance to
  the maps of downstream |E_z| and |E'_z|.
- The theory is heuristic. Eqs. 15-16 treat anti-parallel Speiser acceleration only, and their index
  (~1) is not compared quantitatively with the fitted kappa - 1. The weak-BG enhancement is put
  down to partial magnetisation trapping particles ("may explain"). The paper leaves open whether
  E_z or E'_z drives the acceleration.
- Correlation, not causation: the electric-field argument compares maps, not the work done on
  particles.
- Agreement is noted with French et al. (2022): stronger guide fields suppress acceleration.
- Limits (stated): 2D and pair plasma. Ion-electron partition is deferred, and ion acceleration
  needs larger runs. Background density and temperature are fixed (they set the inflow vA), and BG
  only goes up to B0.
