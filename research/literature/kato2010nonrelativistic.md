# Kato & Takabe (2010): Nonrelativistic collisionless shocks in weakly magnetized electron-ion plasmas: two-dimensional particle-in-cell simulation of perpendicular shock

- Ref: Astrophys. J. 721, 828-842 (2010); DOI 10.1088/0004-637X/721/1/828; arXiv:1008.0265 (v1); peer-reviewed (ApJ)
- Read: 2026-11-12, introduction, set-up, Sect. 3.1, 3.2.1, 3.3, discussion and conclusion (rest skimmed);
  topics: spherical-shock-injection

## Claim
At M_A ~ 130 (sigma ~ 6e-5, SNR-like magnetisation) a 2-D perpendicular shock with in-plane field is shaped by the
ion-beam Weibel (filamentation) instability in the foot, which generates fields comparable to or stronger than
the compressed upstream field. Electrons are heated in the foot by a slowed Buneman instability and then by the
ion-acoustic instability, but end up Maxwellian with T_e/T_i ~ 0.3-0.4: no efficient electron acceleration.

## Method
- 2-D relativistic electromagnetic PIC (standard scheme after Birdsall & Langdon; code not named), injection
  method, downstream rest frame, periodic in y.
- B_0 along y (in the simulation plane), theta_Bn = 90 deg.
- m_i/m_e = 30; bulk speed V = 0.25 c (box frame); V_sh ~ 0.33 c in the upstream frame; sigma~ = 1e-4 (defined with
  V), sigma = 5.9e-5 and M_A ~ 130 with the shock speed (V/v_A = 100); omega_ce/omega_pe ~ 1.4e-2;
  k_B T/m_e c^2 = 1.25e-3 for both species; beta ~ 26 (total; beta_e ~ 13 by my arithmetic). M_s not stated.
- 16384 x 1024 cells, 3200 x 200 c/omega_pe (cell ~0.2 c/omega_pe; L_y ~ 36.5 c/omega_pi, my conversion);
  ~40 particles per cell per species.
- Run to omega_pe t = 8000 (~4 Omega_ci^-1 with their Omega_i^-1 ~ 2000 omega_pe^-1). The shock starts as an
  unmagnetised Weibel-mediated shock and becomes magnetised at omega_pe t ~ 3000-4000 (Sect. 3.1, 3.6).

## Numbers to compare against
- Magnetic energy ~15 % of the upstream bulk kinetic energy density in the transition region and ~2 % downstream
  (Sect. 3.1); |B|/B_0 ~ 40 in the transition region, ~15 downstream (Fig. 17).
- Downstream Maxwellian fits: k_B T_e/m_e c^2 = 0.14, k_B T_i/m_e c^2 = 0.42, T_e/T_i ~ 0.33 (Fig. 14); 0.38 from
  the temperature profile (Fig. 7). Full equipartition would give 0.32 each; no exchange 0.021 and 0.63 (Sect. 3.3).
- Electrons: no suprathermal tail. Ions: tail at E_kin/m_e c^2 ~ 3-6 from ions reflected at the ramp (Sect. 3.3-3.4).
- Jump conditions: V_1/V_2 ~ 3.9, N_2/N_1 ~ 4.1, in agreement with MHD Rankine-Hugoniot (Sect. 3.3).
- Buneman growth rate at x = 2300: 0.0057 omega_pe against ~0.087 omega_pe for cold beams, because the reflected
  ions are hot (Sect. 3.2.1, Fig. 9). Filament width about one ion inertial length (Sect. 3.1).
- Whistler precursor needs roughly M_A < 2 (m_i/m_e)^(1/2) (Eq. 23): <= 11 for m_i/m_e = 30, <= 86 for 1836; the
  filaments here are therefore not whistlers (Sect. 4).

## Assessment
- Shown by simulation: Weibel filaments and field generation in a magnetised perpendicular shock; the two-step
  electron heating of Papadopoulos (1988) checked against local linear theory; a Maxwellian electron spectrum.
- Stated, not tested: the absence of electron acceleration is put down to the field lying in the plane, against
  out of the plane in Amano & Hoshino (2009), "which used similar parameters" (Sect. 4). My remark: that run has
  M_A ~ 14 and beta_e = 0.5, so the two differ in more than geometry.
- Other readings of this result: Matsumoto et al. (2012) blame the high beta_e; Riquelme & Spitkovsky (2011) take
  it as evidence that whistler injection fails at M_A >> (m_i/m_e)^(1/2) for any theta_Bn.
- Stated limits: m_i/m_e = 30; V_sh ~ 0.3 c against ~0.01 c in SNRs; shock speeds fitted by eye and perhaps not
  yet steady; filaments along x cannot merge in 2-D, and 3-D deflection would differ (Sect. 3.5, 3.6).
- Limits I note: one run; high beta; ~4 Omega_ci^-1 of which about half is the unmagnetised start-up phase.
- Conjectured: the Weibel field may relate to the amplified fields inferred at SNR shocks (Sect. 1).
