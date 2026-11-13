# A spherical shock through a planet's own dipole: design

Owner, 2026-11-13: the aim is the acceleration mechanism with the planet as a magnet: a shock that spreads outwards
from a magnetised planet through its dipole magnetosphere, with the dipole dominant over the whole box. The runs in
a uniform field (P5, campaign A) are discarded and their outputs deleted; they are background only. Runner:
`run_dipole.py`. Configs: `configs/dipole_configs.py`.

## 1. Geometry and numbers

- **Field.** A dipole along +z at the centre of the box, inside an absorbing planet of radius 32 d_e. The
  equatorial field is B_eq(r) = 0.333 (70 / r)³ (in ω_ce / ω_pe); at the pole of a sphere it is twice that. The
  uniform field B0 of the example is kept negligible: the dipole equals it 576 d_e out (three half-widths).
- **θ_Bn on a sphere around the planet:** tan θ_Bn = tan θ / 2 (θ the polar angle): 90° at the equator, 0° at
  the poles, 27° at θ = 45°. Most of the sphere is quasi-parallel compared with a uniform field.
- **Plasma at rest:** n0, v_te = 0.1, T_i = T_e, m_i/m_e = 16, everywhere. Local values change with r:

| r [d_e] | B_eq | M_A at the equator (pole) for V = 0.25 c | β_e (equator) | ion gyroradius at V / r (equator) |
|---|---|---|---|---|
| 56 | 0.65 | 1.5 (0.8) | 0.05 | 0.11 |
| 70 | 0.33 | 3.0 (1.5) | 0.18 | 0.17 |
| 88 | 0.17 | 5.9 (2.9) | 0.71 | 0.27 |
| 105 | 0.10 | 10 (5.1) | 2.0 | 0.39 |
| 140 | 0.042 | 24 (12) | 11 | 0.69 |

- **What follows:** the front can be a magnetised, supercritical shock only in a window of about 70–105 d_e at
  the equator (M_A from 3 to 10, ion gyroradius below 0.4 r); further out the ions are not magnetised on the scale
  of the front, and the plasma is of high β (the sound speed, 0.046 c, then sets a fast Mach number of at most
  about 5). At the poles M_A is half that at the equator: the shock forms later and further out there.
- **Debris:** a shell at 48 d_e (width 5 d_e) with a radial drift V_d = 0.25 c and the mass of the ambient
  plasma inside R_m = 90 d_e; its kinetic energy is 3.2 times the dipole field energy outside 48 d_e, so the shell
  pushes the field outwards and is not stopped by the field alone.
- **Grid:** 768³ cells of 0.5 d_e (box 384 d_e), dt = 0.25, 1760 steps (t = 440); Ω_ce Δt ≤ 1.74 at the pole of
  the planet's surface; ω_pe Δt ≤ 0.90 at the peak of the shell; 58.8 GB per GPU.

## 2. Hypotheses on the route to high energy in a dipole

- **D-A, converging mirrors on closed field lines.** A sphere of radius R cuts a dipole field line of apex L > R
  at two points, north and south (sin² θ = R / L). Electrons reflected at one point travel along the field line,
  outside the shock, to the other; as R grows the two points approach and close at the apex. A first-order Fermi
  process between two converging shock points, specific to closed field lines. *Signature:* upstream electrons that
  change hemisphere along the field and are reflected again with a gain.
- **D-B, drift around the planet.** Gradient and curvature drifts of the dipole are azimuthal and close around
  the planet; the expanding flow gives a motional E_φ. *Signature:* gain done by E_φ with a long azimuthal path.
- **D-C, compression of trapped flux tubes.** Electrons trapped on closed field lines are compressed as the shock
  and the flow sweep the flux tube (betatron and Fermi-type compression). *Signature:* gain with p⊥²/B kept
  (betatron) or with the bounce invariant kept (compression of the tube length).
- **D-D, local shock processes:** reflection (SDA) and SSDA at the front, with V_HT and the Mach number set by the
  latitude and the radius. *Signature:* as in `global-plan.md`.
- **D-E, the driver:** gain at the magnetic wall of the debris (as seen in P5): kept apart.

## 3. First runs: D1 and D1-null (fixed before they run)

- **D1:** the design above. **D1-null:** the same planet, dipole and plasma without the debris. A uniform,
  isotropic Maxwellian in a static field is a Vlasov equilibrium; the planet absorbs particles and opens loss
  cones. The null run measures what that does (and the numerical floor) over the same time.
- **Observables of D1:** where and when an ambient compression front (ambient ions ≥ 1.5 n0) appears, by polar
  angle; its speed against V_d and against the local fast speed of the dipole field; the local M_A, β and θ_Bn
  at the front; the radius up to which the ion gyroradius stays below 0.4 r; the energetic electrons (counts,
  local spectra against Maxwell–Jüttner fits, places and fields along tracks) with the tools of P5, adapted to a
  local field instead of B0 before D1 is read.
- **Predictions:**
  - The front first appears at the equator, between r = 65 and 85 d_e (M_A 2.5–5), at t = 70–160 (1–3
    local ion gyro-times at 70 d_e, 1 / Ω_ci = 48 there); near the poles it appears later and further out, or
    not inside the window.
  - Near r = R_m = 90 d_e the equatorial front moves at 0.8–0.9 V_d (the regularity of L, S25 and P5 in units of
    the debris; it may not hold in a field that falls as r⁻³).
  - No electron completes a bounce between the hemispheres on a closed field line within the run: at mid
    latitudes the field line is about 250 d_e long, a round trip takes 500 time units at 0.5 c, and field lines
    with an apex beyond the box leave it. **D-A is not tested by D1**; D1 gives the numbers to design a run that
    tests it.
- **Not shown by D1:** any mechanism as a claim. D1 is an observation run, like S25 and P5 were in a uniform field.

## 4. Limits

- m_i/m_e = 16, 2 particles per cell, one seed. The planner's printed M_A, β and R_B refer to the negligible B0.
- The window of a magnetised shock is a factor 1.5 in radius; the box holds it, not the free expansion beyond.
- The analysis tools of the uniform-field runs take θ_Bn and the fast speed from B0; they are to be adapted to
  the local dipole field before D1 is read.
