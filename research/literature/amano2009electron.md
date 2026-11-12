# Amano & Hoshino (2009): Electron shock surfing acceleration in multidimensions: two-dimensional particle-in-cell simulation of collisionless perpendicular shock

- Ref: Astrophys. J. 690, 244-251 (2009); DOI 10.1088/0004-637X/690/1/244; arXiv:0805.1098 (v2); peer-reviewed (ApJ)
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
In a self-consistent 2-D perpendicular shock with the field out of the plane, Buneman waves at the leading edge of
the foot still accelerate electrons (power law of index 2.0-2.5), but not by trapping as in 1-D. Electrons resonate
stochastically with oblique electrostatic waves, are carried back upstream by the turbulence ("reflected" by
electron-scale fields, not by the cross-shock potential) and then gain energy from the upstream motional field,
in analogy with ion shock surfing. Multidimensionality and the inhomogeneous shock structure are both essential.

## Method
- 2-D electromagnetic PIC (code not named), injection method, downstream rest frame, periodic in y.
- B_0 along z (out of plane), theta_Bn = 90 deg.
- m_i/m_e = 25; omega_pe/Omega_ce = 10; beta_e = beta_i = 0.5; v_A/c = 0.02; injection four-velocity U_0/c = 0.2;
  M_A ~ 14 in the shock rest frame (shock moves at ~0.4 U_0 in the box); sonic Mach number not stated.
- 4096 x 256 cells of one Debye length: L_x ~ 204 c/omega_pe, L_y ~ 12.8 c/omega_pe (~2.6 c/omega_pi, my conversion);
  40 particles per cell per species.
- Analysis at omega_pe t = 1000-1150, with Omega_ci t = omega_pe t/250 (so ~4-4.6); total run time not stated.
- A 1-D run with the same parameters and a 2-D periodic (homogeneous) model are described in words only (Sect. 3).

## Numbers to compare against
- Electron spectral index 2.0-2.5 in the transition region, steepening slightly inward; downstream spectrum equal
  to that at the overshoot; high-energy hump just ahead of the front from reflected electrons (Fig. 4, Sect. 2.3).
- Sample electron: energy x ~20 when ejected upstream, x ~40 on return; mu also x ~40 (Fig. 5, Sect. 2.4).
- Estimated gains in units of the upstream electron bulk energy: ~14 in the turbulent layer (Eq. 3, path
  L_1 ~ 5 c/omega_pe, reflected-ion drift V_r/V_0 = 2 sqrt(2)); ~20 during the upstream half-gyration (Eq. 4,
  L_2 ~ 20 c/omega_pe). The remainder comes from the wave field at the edge of the turbulent layer.
- Overshoot B/B_0 ~ 8 in 2-D against ~13 in 1-D; reformation weaker in 2-D than the 1-2 Omega_ci^-1 cycle of 1-D
  (Sect. 2.2, 3).
- Width of the Buneman layer ~ (m_i/m_e)^(1/3) V_0/omega_pe (Papadopoulos 1988): the real mass ratio widens it ~4
  times in units of the wavelength (Sect. 3).
- For comparison they quote index 3-4 from the 1-D oblique runs of Amano & Hoshino (2007) (Sect. 3).

## Assessment
- Shown by simulation: Buneman wave fronts oblique to the normal and roughly perpendicular to the gyrating
  reflected-ion beam; power-law spectrum; the two-stage history of one electron. Only single reflections were found.
- Stated without figures: the 1-D run with identical parameters grows no strong waves at the leading edge, so the
  1-D/2-D efficiency could not be compared (they suggest the Buneman threshold is lower in 2-D because the
  y-drift of the reflected ions also drives it); the periodic model accelerates less, "basically consistent" with
  Ohira & Takahara (2007).
- Conjectured: reflection persists at the real mass ratio; multiple reflections at higher M_A; maximum energy may
  scale with the wave phase speed (independent of mass ratio), which would lower it relative to the ion energy.
- Stated limits: m_i/m_e = 25; U_0/c = 0.2 (argued harmless because the dominant mode stays electrostatic);
  one run; strictly perpendicular; no direction along B in the plane (no whistlers, no ripples); mass-ratio
  dependence of the maximum energy "not yet clear"; mechanism "not yet fully understood".
- Disagreements it discusses: Ohira & Takahara (2007) (homogeneous model has no edge and no gradient); Dieckmann
  et al. (2008) (planar waves; stronger field omega_pe/Omega_ce = 5 or short run suspected); Umeda et al. (2008)
  (in-plane field gives less acceleration; attributed to the plane in which the reflected-ion beam rotates).
- For SSA+SDA injection at oblique shocks, SSA must reach about the upstream ion bulk energy (their reading of
  Amano & Hoshino 2007); this run does, at m_i/m_e = 25.
