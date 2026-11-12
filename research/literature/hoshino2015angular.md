# M. Hoshino (2015): Angular Momentum Transport and Particle Acceleration during Magnetorotational Instability in a Kinetic Accretion Disk

- Ref: Phys. Rev. Lett. 114, 061101 (2015); DOI 10.1103/PhysRevLett.114.061101; arXiv:1502.02452
- Read: 2026-11-12, full; topics: (none yet)

## Claim
The paper describes itself as the first 3D PIC run of collisionless MRI (pair plasma, shearing box). Kinetic MRI both
accelerates particles and raises angular-momentum transport above MHD values: α ~ 1 during active reconnection and
~ 0.1 afterwards. The proposed cause is that reconnection inside the channel flows leaves p∥ > p⊥ there. That
anisotropy suppresses the next reconnection, so the field saturates higher and the Maxwell stress grows. Outside the
channels, p⊥ > p∥ from the MRI dynamo still favours fast reconnection and acceleration.

## Method
- Same scheme as hoshino2013particle (Eq. 1: Coriolis + tidal force, Ω0 r0 ≪ c, open shearing box).
- 3D, 300³ cells, box 1.91 λ in each direction (λ = 2πVA/Ω0), β = 1536, Ωc/Ω0 = 10, VA/c = 6.25e-3,
  40 particles/cell, pair plasma, uniform vertical B0, isotropic drifting Maxwellian.
- Diagnostics: field lines and channel-flow surfaces (density > mean + 2σ) in Fig. 1; energy spectra (Fig. 1f);
  energy histories and α split into Reynolds, Maxwell and anisotropy stress (Eq. 2, Fig. 2); anisotropy maps
  and histograms (Fig. 3).
- Control experiment: an artificial isotropization. White-noise δB with |δB|/|B| = 2.5 enters the pusher only where
  |B/B0| < 50 (the channel flows), after Torbit = 6.89 (Eq. 5). Fig. 4 compares the two runs.

## Numbers to compare against
- Channel flow at Torbit ≈ 6.89 and first reconnection at ≈ 7.18 (Fig. 1b–d).
- Spectrum (Fig. 1f): nonthermal above ε/mc² ≈ 0.2 at 7.18, against a Maxwellian fit T/mc² = 0.121. It is ∝ ε^-3/2 at
  8.84 and harder over 10–100 mc² after Torbit 9–10. The box-limited εmax is ~ 10^2 mc².
- Magnetic energy rises steeply at Torbit ~ 6, peaks at ~ 8 and stays roughly flat after ~ 9. The instantaneous β is ~ 1 (Fig. 2 top).
- α with the instantaneous pressure (Eq. 2, Fig. 2 bottom): O(1) for 7 < Torbit < 9 with wM > wR ~ wA, and O(0.1) after 9
  with wM ~ wA > wR. MHD shearing boxes quoted for comparison give 1e-3–1e-2.
- Estimate α ~ 1/β (Eq. 3).
- Tearing estimate (Eq. 4, Chen & Palmadesso): B/B0 ~ 230 and δ/Δ ≥ 10 give (rg/δ)^{3/2} ≤ 0.089, while
  1 − p⊥/p∥ < 0.5 in the channel (Fig. 3d).
- With isotropization the volume-mean p⊥/p∥ is higher, α and the Maxwell stress are lower, and wA is unchanged (Fig. 4).

## Assessment
- Shown by simulation: in 3D the channel flows survive the first reconnection and reconnection keeps occurring at several
  sites; in 2D a few large islands formed instead. Anisotropy is p⊥/p∥ ≤ 1 at channel centres and > 1 elsewhere (Fig. 3).
  The control run lowers α.
- Argued, not demonstrated: that p∥ > p⊥ is what suppresses the next reconnection. This rests on a linear-tearing
  order-of-magnitude estimate (Eq. 4), not on a measured onset delay. The late-phase role of p⊥ > p∥ in acceleration
  is carried over from the 2D paper. The link to high-energy particles near massive black holes is stated.
- The spectral hardening may come from multiple reconnection or from particles piling up at the box-limited
  εmax; the paper names both and does not decide.
- Limits: one physical run (plus the control) at one β; 40 particles/cell; a small box (1.91 λ) that caps εmax;
  pair plasma; Ωc/Ω0 = 10; ~14 orbits; the isotropization model is ad hoc (noise amplitude and |B| threshold chosen).
  No resolution or box-size convergence is reported, as expected for the PRL format.
- Left open: the dependence of α on box size, β, Ωc/Ω0 and mass ratio, and the long-time saturated state.
