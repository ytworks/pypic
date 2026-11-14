# Simon, Caprioli, Haggerty & Reville (2026): A "Faux-Shock" Method for Hybrid Simulations of Astrophysical Shocks

- Ref: ApJ 997, 144; DOI 10.3847/1538-4357/ae2611; arXiv:2507.14282. Peer reviewed. Read in the arXiv
  version.
- Read: 2026-11-14, arXiv PDF text: abstract, method, the comparison with the reflecting-wall runs,
  discussion; topics: qpar-electron-conditions

## Claim
A hybrid shock run can be replaced by a box fixed in the shock frame whose downstream side is a
semi-permeable boundary (the "faux shock"): thermal particles pass, energetic particles are returned with a
prescribed probability and energy kick. The upstream length then stays constant and the precursor can be
followed for longer than in a reflecting-wall run.

## Method
- Hybrid (kinetic ions, fluid electrons). Return probability after Peacock (1981), with isotropy assumed
  downstream. A separate cosmic-ray population is injected by prescription.
- 2-D, 10000 × 300 d_i (faux shock) against 20000 × 300 d_i (reflecting wall); M_A 20 in the downstream
  frame (26.67 in the shock frame); n_cr = 0.01 n_g; ~450 ω_c⁻¹.
- One faux-shock run at θ 80°, M_A 100 in 2-D is compared with a 3-D reflecting-wall run.

## Numbers to compare against
- The two methods agree after a time offset of 150 ω_c⁻¹ (the initial transient of the faux shock).

## Assessment
- On the reflecting-wall set-up (ours): the shock propagates through the box, "gradually shrinking the
  upstream length in time".
- Stated price of the method: "giving up the ability to self-consistently account for all the effects that
  occur at the shock, such as particle injection, as well as downstream particle scattering on the
  post-shock magnetic turbulence"; it is "not fully self-consistent in the sense that the shock is
  represented by a sharp discontinuity … and that CR injection is prescribed rather than self-consistently
  injected from the thermal pool". It is to be "used in conjunction with traditional RW simulations that
  establish the appropriate energy gain and escape probability per acceleration cycle, as well as the
  normalization and the fraction of injected particles".
- Ions only. Boundary effects at the far upstream edge are noted.
- For this topic: a boundary-area method that takes exactly the quantities we plan to measure (injected
  fraction, return probability, energy gain per cycle) as inputs. It states the division of labour: a
  self-consistent run fixes injection, a reduced model carries the long times. It also names a limit of our
  wall set-up: the upstream shrinks while the precursor grows (null test N6).
