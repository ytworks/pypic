# Schaeffer, Fox, Haberberger, Fiksel, Bhattacharjee, Barnak, Hu & Germaschewski (2017): Generation and Evolution of High-Mach-Number Laser-Driven Magnetized Collisionless Shocks in the Laboratory

- Ref: Phys. Rev. Lett. 119, 025001 (2017); DOI 10.1103/PhysRevLett.119.025001; arXiv:1610.06533 (read);
  peer-reviewed.
- Read: 2026-11-13, full (arXiv version; it has no methods section, so the PIC set-up was not read);
  topics: spherical-shock-injection

## Claim
First laboratory high-Mach-number magnetised collisionless shock. A laser-driven piston expands into a
pre-formed magnetised ambient plasma. The shock forms in about one ion gyro-time and later leaves the
piston.

## Method
- Omega EP; CH piston into a pre-formed CH ambient plasma in B_0 = 8 T. Time-resolved 2-D images of
  density and field.
- 2-D PIC (PSC) with a model of the ablation. Its box, mass ratio and resolution are not in the text read.

## Numbers to compare against
- Measured: v_s = 700 +- 30 km/s; M_ms = 12 +- 5; density compression at least 3.2-3.7; B/B_0 >~ 3.
- PIC: a shock in the ambient hydrogen by 2.35 ns at M_ms ~ 15. The times 2.35 and 2.85 ns are
  t omega_ci,H ~ 1.4 and t omega_ci,C ~ 0.9, "consistent with the expected shock formation timescale
  t omega_ci ~ 1 for each species".
- At 3.85 ns the profile has two bumps: the shock has left the piston pile-up. The piston ions are then
  trapped behind and "their role in coupling energy and momentum to the ambient ions has effectively
  ended".
- Abstract: "the shocks form on time scales as fast as one gyroperiod, aided by the efficient coupling of
  energy, and the generation of a magnetic barrier between the piston and ambient ions".
- No shock forms without the field or without the ambient plasma.

## Assessment
- For this topic: proof that a piston can make a shock far above the critical Mach number. The piston is
  a laser ablation flow far faster than v_A, not a debris shell of known finite mass.
- The front in the experiment is curved, but the analysis uses 1-D lineouts. Curvature is not studied.
- Electrons are not analysed.
- Left open for us: whether a finite shell keeps its speed up to the separation.
