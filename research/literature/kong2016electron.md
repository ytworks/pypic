# Kong, Chen, Guo, Feng, Du & Li (2016): Electron acceleration at a coronal shock propagating through a large-scale streamer-like magnetic field

- Ref: Astrophys. J. 821, 32 (2016); DOI 10.3847/0004-637X/821/1/32; arXiv:1602.08170 (the arXiv version
  was read); peer-reviewed, the venue's core scope.
- Read: 2026-11-13, full; topics: spherical-shock-injection

## Claim
Test electrons at a circular shock front that expands through a streamer-like coronal field are
accelerated mainly where a closed field line crosses the front at two points: the two mirror points
approach each other as the shock moves out (a collapsing trap), and the efficient region moves from the
flank to the nose of the shock in time.
- Question (Sect. 1): how does the large-scale shape of the coronal field act on electron acceleration at
  an outward-moving curved coronal shock (source of type II bursts)?
- Idea: replace the planar shock of Kong et al. 2015 by a circular front, so that the tangent points of
  front and field lines move in time.

## Method
- Test-particle electrons, full Lorentz orbits (Bulirsch–Stoer), m_e = m_p/1836, in prescribed fields
  (Sect. 2). Not self-consistent; no m_i/m_e, M_A or beta of a kinetic shock.
- Field: axisymmetric streamer model (Low 1986), polar field 10 G, domain r = 1–3 R_sun, cusp at 2.5 R_sun,
  B_y = 0. Closed loops below the cusp, open field above and at the sides.
- Shock: circular front centred at (x = 1 R_sun, z = 0); constant speed U_sh about 1000 km/s
  (U_sh/c about 0.0033); radius R_sh from 0.3 to 1.8 R_sun; compression X_sh = 2.5; tanh flow profile of
  thickness delta_sh = 0.01 U_1/Omega_ci (B_0 about 0.2 G); MHD jump with E = -U x B / c. The radius of
  curvature is thus of order 1e7–1e8 shock thicknesses (our estimate from these numbers, not printed).
- 1e7 electrons of E_0 = 300 eV injected at 10 U_1/Omega_ci ahead of the front at eight times t1–t8.
- Ad hoc pitch-angle scattering: random redirection every tau = 1e4 Omega_ci^-1 (about 5 s).

## Numbers to compare against
- "Because of the variation of shock radius and streamer field topology, the tangent points of the shock
  front and field lines (where theta_Bn about 90 deg) change with time" (Sect. 3, Fig. 2).
- Maximum energy 440 E_0 overall; per injection time t1–t8: 52, 44, 202, 195, 114, 151, 62, 20 E_0.
- Sample electrons: flank on closed field, about 35 E_0 after about 8 reflections (Fig. 3); loop-top trap,
  about 10 to about 70 E_0 within about 2 s (Fig. 4); flank on open field, about 20 E_0 (Fig. 5).
- Fig. 4 text: "the shock front can intersect with the same field line at two different points. With the
  outward propagation of the shock, the two mirroring points approach each other giving rise to a
  collapsing magnetic trap geometry".
- Where the energetic electrons were injected (Fig. 7): "At earlier time (t1 and t2), energetic electrons
  are injected near (but not at) the tangent points"; at t1, 74 deg (81 deg) <~ theta_Bn <~ 87 deg for final
  energy > 10 (20) E_0. Near the streamer axis, theta_Bn >~ 86.5 deg, no efficient acceleration: "particles
  can not be reflected by the shock when theta_Bn is near 90 deg" (citing Holman & Pesses 1983; Ball &
  Melrose 2001).
- The efficient region shifts from the flank to the nose once the radius of curvature of the front exceeds
  that of the field lines it meets (Fig. 9). The flank is less efficient: "only one mirroring point".
- Spectra close to power laws: index about -3 at t3–t4, softening to about -6 at t8 (Fig. 10).

## Assessment
- Shown: in a test-particle model, multiple encounters of one electron with an expanding curved front
  along one field line, and the motion in time of the tangent points and of the injection region.
- Stated limits: analytic field, MHD jump, prescribed scattering, no feedback, constant shock speed.
- Open (Sect. 4): "Further studies should be conducted to investigate electron acceleration in a
  self-consistently-solved coronal shock environment."
- For this topic: the trap needs field lines that are curved more strongly than the front (concave toward
  the shock). With a straight uniform field and a convex front, the segment between the two contact points
  lies downstream, so this trap does not carry over. The time dependence of theta_Bn at a contact point
  and "injection near but not at the tangent point" are published here for an expanding front; the shock
  does not decelerate, and the electrons are test particles.
- Doubts: the result depends on the scattering time chosen; a shock 0.01 U_1/Omega_ci thick makes every
  encounter a clean mirror reflection.
