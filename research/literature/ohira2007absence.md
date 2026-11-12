# Ohira & Takahara (2007): Absence of electron surfing acceleration in a two-dimensional simulation

- Ref: Astrophys. J. 661, L171 (2007) (volume and page from the reference lists of Amano & Hoshino 2009 and Kato &
  Takabe 2010; the arXiv metadata gives no journal reference, comment "accepted for publication in ApJL");
  DOI 10.1086/518888; arXiv:0705.2061 (v2); peer-reviewed (ApJ Letters)
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
In two dimensions the Buneman instability in the foot of a perpendicular shock makes potential structures that are
coherent over only about one wavelength along the direction of the motional electric field. Electrons leave the
potential sideways within a gyro-period, so the surfing acceleration of 1-D simulations does not occur; the outcome
is heating to k T_e ~ m_e V_d^2, far below injection energies.

## Method
- 2-D electrostatic PIC (Poisson equation, no electromagnetic waves), doubly periodic box that models the foot only
  (no shock), upstream rest frame; B along z (out of plane); geometry of a perpendicular shock (theta_Bn = 90 deg).
- Three uniform Maxwellian populations: incoming electrons, incoming protons, reflected protons drifting at
  V_d = -0.04 c = 2 V_sh; n_e = 1.25 n_p = 5 n_r = 1 cm^-3; all at T = 7 eV (V_d ~ 7 electron thermal speeds).
- m_i/m_e = 1836 (real).
- B scanned from 0 to 72 muG (omega_ce/omega_pe = 0-0.08); the case shown is 27 muG (0.03), chosen because it gives
  the largest electron energy in 1-D.
- M_A and beta not stated (from the stated n, B, T, V_sh: v_A ~ 59 km/s, M_A ~ 100, beta_e ~ 0.4; my arithmetic).
- Box 16 x 16 lambda_Bun (most unstable Buneman wavelength), 256 x 256 cells, cell = 3 Debye lengths; 80 electrons
  per cell in 2-D and 10240 in the 1-D comparison run; dt = 0.01 omega_pe^-1; run to 1000 omega_pe^-1.
- Further runs at T = 1.75, 28 and 112 eV, described in words (Sect. 3).

## Numbers to compare against
- Field energy peaks at ~100 omega_pe^-1 with e-folding time ~10 omega_pe^-1, as linear theory predicts; E_y is
  comparable to E_x; the 2-D field is a few times weaker than in 1-D (Fig. 1, Sect. 2.2).
- Coherence length along y of order lambda_Bun (Fig. 2); only a few times longer at T = 28-112 eV (Sect. 3).
- Cold dispersion relation (Eq. 1) does not contain k_y: all k_y are unstable once k_x is (Sect. 2.2).
- 1-D: jet-like accelerated electrons and a high-energy component; 2-D: none (Fig. 3 at 690 omega_pe^-1, Fig. 4 at
  1000 omega_pe^-1). The 2-D run with B is slightly hotter than the one without.
- Final electron temperature ~ m_e V_d^2 (Sect. 2.4), of order keV, "4 orders of magnitude" short of what DSA
  injection needs (Sect. 3).
- Energy change over one coherence length: dE = 2 pi alpha m_e V_d^2 (omega_ce/omega_pe), alpha = O(1) (Eq. 3-4):
  negligible at omega_ce/omega_pe = 0.03.

## Assessment
- Shown by simulation: in a homogeneous, periodic, electrostatic model with the real mass ratio, no nonthermal tail
  forms in 2-D for the field strengths and temperatures tried; the 1-D counterpart does show surfing.
- Stated beyond the evidence: that surfing fails "in real cases", and that it cannot serve as injection; the 3-D
  case is argued from the equivalence of y and z at weak field, not simulated.
- Stated limits: electrostatic only (electromagnetic modes may matter after saturation); one parameter set shown;
  no 3-D.
- Limits I note: no shock structure, hence no finite-width turbulent layer, no gradient and no escape upstream;
  run time 1000 omega_pe^-1 is far below an ion gyro-period.
- Later works on this result: Amano & Hoshino (2009) find the same weak effect in a periodic model but efficient
  acceleration in a self-consistent shock; Matsumoto et al. (2012, 2013) find multi-encounter surfing on coherent
  oblique wave trains when M_A >~ (m_i/m_e)^(2/3). The paper itself disagrees with the 1-D results of Shimada &
  Hoshino (2000), Hoshino & Shimada (2002) and McClements et al. (2001) as a guide to real shocks.
- The paper prints lambda_Bun = 2 pi omega_pe/V_d (Sect. 2.1); dimensionally it must be 2 pi V_d/omega_pe.
