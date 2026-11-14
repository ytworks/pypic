# Caprioli & Spitkovsky (2014): Simulations of Ion Acceleration at Non-relativistic Shocks. III. Particle Diffusion

- Ref: ApJ 794, 47; DOI 10.1088/0004-637X/794/1/47; arXiv:1407.2261. Peer reviewed. Read in the arXiv
  version. Citekey letter "b": `caprioli2014simulations` (paper I of the series) is on
  `research/spherical-shock-injection`.
- Read: 2026-11-14, abstract and the passages on how the diffusion coefficient is obtained (arXiv HTML);
  topics: qpar-electron-conditions

## Claim
In hybrid simulations of parallel shocks the ions diffuse in their self-generated turbulence with a mean
free path comparable to their gyroradius. Moderately strong shocks follow quasi-linear theory with the
self-generated spectrum; very strong shocks show Bohm-like diffusion in the amplified field.

## Method
- 2-D and 3-D hybrid simulations (kinetic ions, fluid electrons).
- D(p) from the steady upstream profile of the accelerated ions, D(p) ~ u F_up(p) / f_sh(p), and from
  particle tracking.

## Numbers to compare against
- Moderately strong shocks: D from the energy in modes resonant with the gyroradius (quasi-linear).
- Very strong shocks: D close to Bohm in the amplified field.
- The enhanced diffusion sets the maximum energy reached in a given time.

## Assessment
- Ions only. This is the source of the D(p) method and the "Bohm in the amplified field" result that
  gupta2025speed extends to electrons.
- For this topic: at p >~ m_i v_sh an electron has the same rigidity as an ion and, in the magnetostatic
  limit, the same D. What is new for electrons is the range below that, which this paper cannot reach.
- The upstream-profile method needs a steady precursor of the species; for electrons in a short run the
  mean-square-displacement method is used instead.
