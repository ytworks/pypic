# Riquelme & Spitkovsky (2011): Electron injection by whistler waves in non-relativistic shocks

- Ref: Astrophys. J. 733, 63 (2011) (volume and page from the reference list of Matsumoto et al. 2012; the arXiv
  metadata gives only the DOI); DOI 10.1088/0004-637X/733/1/63; arXiv:1009.3319 (v2); peer-reviewed (ApJ)
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
Oblique whistler waves in the foot of quasi-perpendicular shocks energise electrons through their electric field
component along B. Those electrons can outrun the shock along the field, are scattered again, and later gain
energy from the motional field in an SDA-like way. The process needs theta_Bn != 90 deg and M_A <~ (m_i/m_e)^(1/2),
and it hardens with mass ratio, whereas Buneman surfing fades at realistic mass ratio. For SNRs it requires
M_A <~ 20, hence strong amplification of the upstream field.

## Method
- PIC code TRISTAN-MP, reflecting wall, downstream rest frame, jumping injector; spectra at t omega_ci = 10.
- 2-D with B_0 in the plane (21 runs, Table 1); 2-D with B_0 15 deg from the plane normal (7 runs, Table 3); 3-D
  with only 3 c/omega_pe along z (3 runs, Table 2).
- M_A = 3.5, 7, 14 (21 in two out-of-plane runs), with v_sh measured in the upstream frame: v_sh/c = 0.14
  (v_in = 0.1 c) or 0.042; m_i/m_e = 25, 100, 400, 1600; beta_e = beta_i from 0.005 to 0.5 (4.5 at M_A = 21).
- theta_Bn = 90, 80, 75, 70, 60, 45 deg at M_A = 7, m_i/m_e = 400; mass-ratio scans at 75 deg (M_A = 3.5, 7) and
  60 deg (M_A = 14); 3-D and out-of-plane runs at 75 deg.
- omega_pe/Omega_ce not stated; the times in the Fig. 9 caption imply 10, 5, 2.5 for m_i/m_e = 25, 100, 400 at
  M_A = 7, v_sh = 0.14 c (my arithmetic), so it varies along the mass-ratio scan.
- c/omega_pe = 10-15 cells (5 in 3-D); L_y = 34-102 c/omega_pe; 10-100 particles per cell in total; L_x not stated.

## Numbers to compare against
- Run 2D-3 (M_A = 7, 75 deg, 400): downstream index alpha ~ 3.6; Gamma_max ~ 50 ~ (m_i/m_e) v_sh/c, i.e. electron
  and ion gyro-radii equal; no ion tail (Fig. 2, 13). Whistlers ~10 c/omega_pe, ~45 deg to the normal (Sect. 3.1).
- theta_Bn scan (Fig. 7): hardest at 70 deg (alpha ~ 3.6); Maxwellian with a small bump at 90 deg; softer below 70 deg.
- Mass ratio at M_A = 7, 75 deg: alpha ~ 5.4, 4.7, 3.6 for 25, 100, 400 (Sect. 4.2). At 1600: alpha = 2.7 at
  M_A = 3.5 rising to 4 at M_A = 14 (Sect. 4.4, Fig. 12). Hard and soft regions split along M_A = (m_i/m_e)^(1/2).
- Nonthermal electrons: ~2-10 % by number, ~10-30 % by energy (Sect. 4.4, footnote 4).
- v_sh 0.14 c to 0.042 c: same indices, tail normalisation lower by ~3 (Sect. 4.3). beta_e 0.5 to 0.05: no change.
- Whistler conditions quoted: MTSI needs cos(theta) >~ 4(1-r) M_A/(m_i/m_e)^(1/2) or 4 r M_A/(m_i/m_e)^(1/2),
  r ~ 0.2; M_w = |cos theta_Bn| (m_i/m_e)^(1/2)/2 (Sect. 3.1). Escape ahead of the shock:
  100 cos^2(theta_Bn) >~ M_A/(m_i/m_e)^(1/2) (Sect. 5).
- SNR: alpha ~ 1 - log(eta_inj)/log(c/v_sh) (Eq. 1) gives alpha ~ 4 for eta_inj = 1e-6 and 3000 km/s, reached only
  for M_A <= 14-20 at m_i/m_e = 1600; unamplified M_A ~ 300, so the upstream field must grow ~25 times (Sect. 5).
- Surfing (Appendix, out-of-plane 2-D): Gamma-1 ~ 1 at v_sh = 0.14 c and ~0.1 at 0.042 c, the same for all mass
  ratios; Buneman field down ~4 times from m_i/m_e = 25 to 400. 2-D and 3-D spectra differ at 25, agree at 400.

## Assessment
- Shown by simulation: the scans above; work done by E_parallel dominates the first energisation of a traced
  electron (Fig. 3-4); whistler acceleration vanishes when the direction along B_0 is not in the plane.
- Assumed: the best angle does not depend on m_i/m_e or v_sh. Extrapolated: bow shock and SNR speeds.
- Called speculative by the authors: cosmic-ray-driven turbulence amplifies and reorients a quasi-parallel
  upstream field, so quasi-perpendicular patches cross the shock intermittently and inject electrons locally.
- Stated limits: partial parameter coverage; whistler scale not resolved along z in 3-D; Gamma_max limited by
  particle statistics; no upstream turbulence, so no DSA.
- Limits I note: M_A <= 14 with in-plane field; many runs with 10 particles per cell; mass ratio and
  omega_pe/Omega_ce change together.
- Disagreements: surfing judged negligible at 1836, against Amano & Hoshino (2009) and later Matsumoto et al.
  (2012), who object that M_A = 7 fails their trapping condition; Levinson (1992, 1994) predicted whistler
  injection at M_A >~ 43, not seen here nor by Kato & Takabe (2010).
