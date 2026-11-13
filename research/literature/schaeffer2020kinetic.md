# Schaeffer, Fox, Matteucci, Lezhnin, Bhattacharjee & Germaschewski (2020): Kinetic simulations of piston-driven collisionless shock formation in magnetized laboratory plasmas

- Ref: Phys. Plasmas 27(4), 042901 (2020); DOI 10.1063/1.5123229; accepted manuscript read from OSTI
  (purl 1643946); peer-reviewed.
- Read: 2026-11-13, full; topics: spherical-shock-injection

## Claim
A piston-driven perpendicular shock forms in three steps that scale with the upstream ion gyro-time and
with the gyroradius of an ion at the piston speed: onset, separation from the piston, and a developed
downstream. The times and distances depend little on the Mach number. Compression and ion heating
appear before any shock exists, so they do not prove a shock.

## Method
- Full PIC (PSC), "quasi-1D": 30000 x 12 cells, uniform drive in the transverse direction, so the
  ablation is planar. Strictly perpendicular field. m_i/m_e = 100 (one check at 400). Binary collisions on.
- The piston is made by a heating operator on a dense target (2.5 n_e,ab), which models laser ablation.
  New particles are added to keep the target density, and "the heating operator was applied for the
  lifetime of the simulation" (Sect. II).
- Runs cover M_A = v_sh/v_A = 3 to 57 and M_ms = 3 to 27 (Sect. II.2).
- Reference run (Table I): L_z = 900 d_i,ab; d_i0 = 11.2 d_i,ab; v_p = 730 km/s; v_sh = 980 km/s;
  M_A = 14; M_ms = 13; beta_0 = 0.2. Derived here: box about 80 d_i0, run about 6.5/omega_ci0.

## Numbers to compare against
- Piston speed v_p: speed of the peak field at the cavity edge (Sect. III A).
- Shock speed against piston speed (Sect. III A, Eq. 2, Fig. 3b): "The shock speeds asymptote to about
  4/3 times the piston speed at low background fields and then increase with increasing B_0, and are
  consistent with the RH relations." In the limit of large M_A, v_sh/v_p -> 4/3.
- Formation times and places (Sect. III B, Fig. 6), with rho_i0 = v_p/omega_ci0:
  - onset: t*_1 ~ 1/omega_ci0 at z*_1 ~ 1 rho_i0 (highest rate of reflected ambient ions);
  - separation from the piston: t*_2 ~ 2.5/omega_ci0 at z*_2 ~ 2.5 rho_i0;
  - developed downstream: t*_3 ~ 5/omega_ci0 at z*_3 ~ 5-6 rho_i0.
  They are "nearly constant for a wide range of conditions, with little dependence on the Alfvenic Mach
  number M_A or ablation beta".
- Drive duration (Sect. III B): "The formation time is also not sensitive to the duration of the heating
  operator, as long as heating is applied for at least ~1 omega_ci0^-1. This is because the shock is
  driven by the interaction of the piston and ambient ions, which primarily occurs near the leading edge
  of the piston plasma."
- Seven criteria for a formed shock (Sect. III B): super-magnetosonic; collisionless; n_e/n_e0 > 2;
  B/B_0 > 2; ramp width about d_i0; reflected ambient ions; separation from the piston. The first six
  hold at t*_1, all seven at t*_2.
- Reference run (Fig. 5): at t*_2 the shock is near 29 d_i0 and the piston near 25 d_i0; at t*_3 they
  are near 56 and 46 d_i0.
- Fig. 9 text: at t omega_ci0 ~ 1.5 "the shock-piston separation is most pronounced at lower Mach numbers".
- Later evolution (Sect. III C): the shock leaves the piston because v_sh > v_p; it reforms about every
  1.5/omega_ci0; the Rankine-Hugoniot values are reached by about 10/omega_ci0 (runs at M_A = 7.1, 3.2).
- Mass ratio (Sect. III F): 100 against 400 gives t*_1 = 1.1/omega_ci0 in both, z*_1 = 1.1 and 1.2 rho_i0.
- Each ion species forms its shock in about one of its own gyro-times (Sect. III D).
- Sect. IV: "for a given Z/m slower shocks (lower M_A) are favorable as the shock will form closer to the
  target".
- Sect. I, on earlier hybrid work at low Mach number: those runs were "in a regime in which the piston
  lifetime is comparable to the shock formation time".
- Electrons: "the electrons do not play a significant role at this stage in shock formation" (Sect. III).
  No electron acceleration is analysed.

## Assessment
- For this topic: the only quantitative formation table found. It gives the room a driven front needs:
  1 rho for the onset, 2.5 rho for a front that has left the piston, 5-6 rho for a developed downstream.
  rho here uses the piston speed, not the first debris speed.
- A design with 3.2-4.5 rho_d of room covers the onset and the separation, but not the developed
  downstream. This is our reading; the paper has no finite piston and no sphere.
- The Mach number of the shock is set by the piston speed at the time of formation (4/3 v_p at high
  M_A). A debris shell that has slowed down before 2.5 rho gives a weaker shock than its first M_A.
- The piston here is driven for the whole run and is planar. The paper states the limit itself:
  "Laser-driven plasmas expand approximately hemispherically, and could thus drive shocks over a range
  of non-perpendicular magnetic geometries" (Sect. IV). A finite debris mass and the geometric weakening
  of a spherical front are left open.
- The lowest mass ratio tested is 100. Nothing is said about 16.
