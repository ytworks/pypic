# van Marle, Bohdan, Morris, Pohl & Marcowith (2022): Diffusive Shock Acceleration at Oblique High Mach Number Shocks

- Ref: ApJ 929, 7; DOI 10.3847/1538-4357/ac5962; arXiv:2203.00353. Peer reviewed. Read in the arXiv
  version.
- Read: 2026-11-14, arXiv PDF text: abstract, set-up of the PIC and PIC-MHD runs, the injection rates,
  discussion and conclusions; topics: qpar-electron-conditions

## Claim
A two-step chain: 2-D PIC runs give the ion injection rate at oblique shocks; a 2-D PIC-MHD run with that
rate as input then follows the acceleration on larger scales. At M_A 20, diffusive shock acceleration is
unlikely at θ ≥ 60°.

## Method
- PIC: modified TRISTAN, 2-D, m_i/m_e = 50, 40 particles per cell and species, electron skin depth 8
  cells, v_sh 0.133 c and 0.264 c, M_A 20-30 and more, θ 55-65°, field in and out of the plane.
- PIC-MHD: 2-D, 240 × 30 r_g, with "a fixed fraction of the gas crossing the shock … converted into
  non-thermal particles".

## Numbers to compare against
- Injection rate against shock speed: ε_inj ≈ 12 ε_sh at 0.133 c and ≈ 23 ε_sh at 0.264 c. "the injection
  rates are not strongly dependent on the shock speed at non-relativistic speeds. A strong dependence does
  develop at relativistic speeds". Ions, oblique shocks.

## Assessment
- Limits of PIC, as stated: "The simulation time of PIC simulations is limited by the computational
  expense, especially at high Alfvénic Mach numbers and low shock velocities. Also, if a simulation is too
  long, artificial heating can substantially change the upstream conditions."
- On dimension: "the 2D results are not fully representative of the 3D reality. Therefore, we choose an
  injection rate in between"; "injection is, fundamentally, a 3D process".
- On the box: "an insufficient length of the simulation box in the direction of flow may change the
  outcome"; "Box-size limitations can be an issue if one want to investigate the raise of the maximum
  particle momentum as function of time".
- On electrons: "the exact fraction of injected electrons still should be defined with further PIC
  numerical experiments".
- On PIC-MHD: it "cannot model the process by which particles start to deviate from thermal equilibrium".
- SN 1006: with M_A 100-200 ion injection would be possible up to ~65°.
- For this topic: the closest published form of our chain "measure injection in PIC, feed a larger model",
  done for ions at oblique shocks. It shows what such a chain had to concede: the injection rate fed in was
  chosen between two 2-D values. It names long-run artificial heating (null test N1) and box length (N6) as
  limits, and gives a measured dependence of an injection rate on shock speed over a factor 2.
