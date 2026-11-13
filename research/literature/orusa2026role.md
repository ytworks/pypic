# Orusa, Caprioli, Sironi & Spitkovsky (2026): The role of three-dimensional effects on ion injection and acceleration in perpendicular shocks

- Ref: ApJ 1001, 158 (2026); DOI 10.3847/1538-4357/ae563e; arXiv:2507.13436 (v1 2025; v2 17 Apr 2026,
  "Matches version published by ApJ"); peer-reviewed. The topic's survey lists it under 2025, the year of the
  first arXiv version.
- Read: 2026-11-13, skimmed (abstract, set-up, conclusions); topics: spherical-shock-injection
- Source: the arXiv full text (HTML of v2). Read: abstract, Sect. II (set-up, Table 1), Sect. III (results,
  by figure) and Sect. IV (conclusions). "Paper I" in it is `orusa2023fast`.

## Claim
At exactly perpendicular shocks, ions are injected and accelerated by shock drift acceleration (SDA) only in
3-D. The cause is the "porosity" of the downstream magnetic turbulence near the shock: "how easily the
post-shock region allows particles to traverse it and return upstream without being trapped" (abstract). In
2-D the amplified field forms "walls" that span the invariant direction and no ion returns. Acceleration
"remains an on–off process governed by the dimensionality of the simulation" (Sect. IV).

## Method
- Hybrid code dHybridR: kinetic ions, fluid electrons. A flow against a reflecting wall; theta = 90 deg;
  thermal speed v_th = v_A, T_e = T_i, plasma beta = 2; M_s = M_A/sqrt(5/3) (Sect. II).
- Transverse box 20 d_i in y and z; 8 ions per cell. The 2-D runs have B_0 in the plane (Sect. II).
- Runs (Table 1): A M_A = 30, B M_A = 60, D M_A = 100 at dx = 0.1 d_i; C M_A = 100 at dx = 0.4 d_i.
- Resolution chosen by agreement of a 2-D hybrid run with a 2-D full PIC run (Tristan-MP, m_i/m_e = 49; 256
  also tested) at M_A = 30 (Fig. 1). Finer grids bring non-physical whistler modes.

## Numbers to compare against
- Table 1, at t = 32/omega_c: energy fraction above 10 E_sh and spectral index alpha: A 0.3 %, 6; B 0.6 %,
  5; D 0.9 %, 4; C (coarse) 10 %, 2.5. E_sh = m v_sh^2 / 2.
- Spectra are the same in 2-D and 3-D up to 6-7 E_sh; ions above that exist only in 3-D (Fig. 2).
- Test particles at 7.5 E_sh, with r_L* = (M_A/4) sqrt(E/E_sh) d_i: in 2-D the probability of return drops
  to zero at a "return horizon" of 1.31 r_L* for M_A = 30 and 100; in 3-D "a significant fraction of
  particles starting well beyond 1.31 r_L*" return (Fig. 5).
- "shocks with M_A <~ 10 do not efficiently accelerate particles"; they stay laminar (Sect. III.5, citing
  Paper I).

## Assessment
- For this topic: it is about ions, and the hybrid model has no kinetic electrons. Its weight here is on the
  method: at the perpendicular end the return of particles upstream is a 3-D effect, and a 2-D result there
  can be qualitatively wrong. That supports doing this topic in 3-D; whether the same holds for electrons is
  not shown.
- The threshold M_A <~ 10 for no ion acceleration at perpendicular shocks is a number to place our front
  against. By it, a front below that value should show no SDA-accelerated ions at the equator (our reading).
- The authors state that oblique shocks "are not covered in this work"; there "particles could escape
  upstream by streaming along the magnetic field lines". A sphere has all obliquities side by side. The
  flat-shock result does not say what the neighbouring oblique regions do to the perpendicular one.
