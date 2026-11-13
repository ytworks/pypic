# Qin, Kong & Zhang (2018): Effects of shock and turbulence properties on electron acceleration

- Ref: Astrophys. J. 860, 3 (2018) (Crossref); DOI 10.3847/1538-4357/aac26f; arXiv:1802.08367 (v2 read);
  peer-reviewed
- Read: 2026-11-13, full; topics: shock-alfvenic-turbulence

## Claim
In test-particle runs with a prescribed planar shock and static slab + 2-D turbulence, electron acceleration at a
perpendicular shock gets weaker as the turbulence level (b/B0)^2 rises, and strongest at 0.01, where electrons stay
near the front for long shock drift acceleration (SDA). At a parallel shock the trend is the reverse. Oblique shocks
(30-60 deg) accelerate weakly at every level.

## Method
- Test electrons, Lorentz equation in the shock frame, RK4 with adaptive step (accuracy 1e-9). Shock: planar,
  tanh flow profile U(z) of thickness L_diff = 9.28e-6 AU (Eq. 1); E = -U x B. Nonrelativistic or relativistic
  pusher not stated (momentum form, Eq. 4).
- Table 1: U1 = 500 km/s; compression s = 4 (also 2.6); B01 = 3 nT; M_A1 = 10; theta_Bn = 0-90 deg in 15 deg
  steps; (b/B0)^2 = 0.01, 0.1, 1, 10; electrons 1 keV (upstream plasma frame), isotropic, released at
  z0 = -5.8e-5 AU; run time 23.2 min (the acceleration time of Guo & Giacalone 2015).
- Field: B = B0(z) + b, b time-independent (magnetostatic) and transverse to B0; slab : 2-D energy = 20 : 80.
  Slab correlation scale lambda = 0.02 AU, 2-D scale lambda/2.6. Slab spectrum has an inertial range (5/3) and a
  dissipation range (index 2.7, also 3.4) from k_b = 3e4/lambda; the 2-D part has no dissipation range (Eq. 3).
  Turbulence boxes 10 lambda (2-D) and 50 lambda (slab), Nz = 2^22, Nx = Ny = 4096. Fully 3-D particle motion.
- Efficiency measures: R% = fraction above 10 keV and mean energy E_aver at the end (Fig. 4); residence time
  within one gyroradius of the front t_rs and number of crossings n_cross (Fig. 5).

## Numbers to compare against
- Perpendicular: spectra alike at 1-10 keV for all levels; above 10 keV, 0.01 and 0.1 exceed 1 and 10; above
  100 keV, 0.01 is "much larger" than 0.1 (Fig. 2). Parallel: (b/B0)^2 = 10 reaches ~2 MeV, 1 reaches ~1 MeV.
- R% and E_aver fall from 0 to 15 deg and rise from 30 to 90 deg; they rise with turbulence for theta_Bn <= 60 deg
  and fall with it for theta_Bn > 75 deg (Fig. 4).
- t_rs and n_cross: perpendicular and parallel equal at (b/B0)^2 ~ 0.3; E_aver(perp) > E_aver(par) below ~0.4 and
  < above (Fig. 5).
- Shock thickness: a bend-over thickness L_diff,b with R_gyro,b ~ L_diff,b; efficiency flat below, falling fast
  above (Figs. 6-9); stated to be of order an ion inertial length for these parameters.

## Assessment
- Shown (test particles): with a static field and a planar, structureless shock, perpendicular-shock electron
  energization is SDA-dominated and is reduced when turbulence adds a field component along the normal.
- Opposite trend to Guo & Giacalone (2010, 2015). The authors' only stated reason: "the difference in the
  turbulence models we adopt" (Sect. 4): slab+2-D with a dissipation range vs Guo & Giacalone's isotropic field
  from the induction equation. They do not test it.
- My reading (not stated by the authors): b is static in the shock frame and is not advected or compressed by the
  flow, while E = -U x B; such a field is not frozen in (curl E != 0 is not balanced by dB/dt), so field lines do
  not cross the shock with the flow, the downstream field is not amplified, and there are no front ripples. The
  meandering-line return of Jokipii & Giacalone (2007) needs convected lines; this model may not contain it. The
  2-D component (80 %) also has lambda_x = 0.0077 AU, far above an electron gyroradius.
- The level is varied at fixed M_A = 10 and fixed spectrum shape; no seed statistics; one initial energy.
- The abstract's "SDA is more efficient than first-order Fermi" is an assumption offered as a reason, not shown.
