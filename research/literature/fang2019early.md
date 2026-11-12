# Fang, Lu, Yan & Yu (2019): Early acceleration of electrons and protons at the nonrelativistic quasiparallel shocks with different obliquity angles

- Ref: Res. Astron. Astrophys. 19, 182 (2019) (volume and article number as printed in the reference list of
  Kumar & Reville 2021; the arXiv metadata gives only the DOI); DOI 10.1088/1674-4527/19/12/182;
  arXiv:1908.08170v1; peer-reviewed (comment: accepted by RAA)
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
1-D PIC runs at 15, 30 and 45 deg with otherwise equal parameters all give non-thermal electron and proton tails.
Shock drift acceleration is more effective at the larger angle (hardest upstream electron spectrum at 45 deg),
while the downstream tails reach the highest energies at 15 deg. DSA is not yet developed when the runs end.

## Method
- PIC, 1D3V, code of Derouillat et al. (2018) (name not given in the text); plasma at v_0 = 0.2 c against a
  reflecting wall, downstream frame.
- m_p/m_e = 30. T_e = T_p = 0.2e-3 m_e c^2 as printed. Stated M_s = 13.4 and M_A = 16.5 with v_sh = v_0 r/(r-1),
  r = 4, c_s = sqrt(2 gamma_ad k T_p/m_p), v_A = B_0/sqrt(4 pi n_0 m_p). The paper calls these "in the downstream
  shock frame"; see the doubt under Assessment.
- beta not stated. Proton gyro-time 338/omega_pe and gyroradius m_p v_0 c/(e B_0) = 67.6 c/omega_pe (so
  omega_pe/Omega_ce = 11.3, my arithmetic; the same magnetisation as Kato 2015).
- theta_Bn = 15, 30, 45 deg; B_0 in the x-y plane.
- Box 6.4e4 c/omega_pe, dx = 0.1 c/omega_pe, dt = 0.045/omega_pe, 64 particles per cell per species.
- Run to 1.89e5/omega_pe (= 559 proton gyro-times of 338/omega_pe, my arithmetic).

## Numbers to compare against
- 15 deg, downstream (x = 1.0-1.2e4 c/omega_pe): electron spectrum leaves the Maxwellian at ~5 m_e c and cuts off
  near 50 m_e c; proton tail from ~15 to ~100 m_e c; thermal fit k T = 0.2 m_e c^2 for both (Sect. 3.1, Fig. 5).
- Precursor wave power F(k) at 15 deg is largest for 2e-3 <= k <= 2e-2 omega_pe/c, steeper than k^-1; from the
  polarisation the authors judge the non-resonant instability to dominate at the end (Sect. 3.1, Fig. 4).
- Compression ratio 4 far downstream at 15 and 45 deg (Figs. 1, 6).
- Downstream electron tails have index "about 4" in momentum space (Sect. 3.2, Fig. 9); no fit is described.
- Obliquity trends (Sect. 3.2, Figs. 9-10): upstream (x = 1.3-1.5e4 c/omega_pe) electron spectrum hardest at
  45 deg; upstream proton spectra similar at 15 and 45 deg; downstream tails of both species extend furthest
  at 15 deg; at 45 deg more energetic particles stay confined upstream.
- K_ep, injected fractions, energy fractions, injection momenta: not given.

## Assessment
- Shown by simulation: tails at all three angles; reflected protons upstream and wave growth at each angle; a
  qualitative ordering of the spectra with angle at one time.
- Interpretation, not measured: the split into SDA and DSA contributions. No particle orbits are shown; the
  mechanism (mirror reflection, SDA while trapped by upstream waves, then DSA) is taken from Park et al. (2015).
  The paper itself says that at the end most accelerated particles are still in the SDA stage.
- The 30 deg run appears only in the spectra of Figs. 9-10; no figure values are quoted in the text for it.
- Limits stated: early stage only; longer runs with a larger box are planned. Not discussed: 1-D geometry, the
  small mass ratio, or the particle number (64 per cell at dx = 0.1 is ~640 per skin depth, close to the value
  at which Kato 2015 still saw numerical energy loss; my remark).
- Doubts: (1) from the printed T, m_p/m_e and sound-speed formula I get c_s = 4.7e-3 c and M_s ~ 57 with
  v_sh = 0.267 c, not 13.4; 13.4 follows from v_0/c_s with T ten times larger (2e-3 m_e c^2, the value used by
  Kato 2015). M_A = 16.5 is consistent with v_sh and the stated gyro-time. (2) "n_0 = 0.1 cm^-1" and a
  duplicated definition of the Stokes parameter I (Eq. 2 text) are misprints. (3) The 45 deg run is called
  quasi-perpendicular in the captions of Figs. 7-8 and quasi-parallel in the title.
- Relation to others: the introduction cites Wieland et al. (2016), 2-D PIC at 45 deg, as finding a supra-thermal
  ion tail but no electron acceleration; this paper finds electron tails at 45 deg in 1-D and does not discuss
  the difference.
