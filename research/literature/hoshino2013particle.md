# M. Hoshino (2013): Particle Acceleration during Magnetorotational Instability in a Collisionless Accretion Disk

- Ref: ApJ 773, 118 (2013); DOI 10.1088/0004-637X/773/2/118; arXiv:1306.6720
- Read: 2026-11-12, full; topics: (none yet)

## Claim
In a 2D (radial–vertical) kinetic shearing box, the field stretching of the MRI drives p⊥ > p∥ by
double-adiabatic (CGL) behaviour. The anisotropy left in the channel-flow current sheets
(p⊥/p∥ ≈ 1.3) makes tearing/reconnection grow quickly (tearing coupled to the mirror mode). Repeated
formation and disruption of current sheets then produces a hard nonthermal tail, N(ε) ∝ ε^-1, at
every β studied.

## Method
- Relativistic PIC in a local rotating frame: Coriolis and tidal forces in the pusher (Eq. 8, q = 3/2).
  The non-inertial terms in Maxwell's equations are dropped (v0 = 0, Sec. 2.1); test runs are said to show little change.
- Pair plasma to save cost. The linear pair-plasma MRI dispersion reduces to the MHD one (App. A, Eq. A12).
- Open shearing-box boundary in x by Lorentz-transforming fields, density and current (Eqs. 11–14; Hawley et al. 1995);
  periodic in y, z. The code is 3D, but only the x–z plane is run.
- Initial state: vertical B with a small correction for charge separation (Eq. 10, App. B), an isotropic drifting
  Maxwellian with Keplerian shear, β = 96, 384, 1536, 6144 (Runs A–D, Table 1), Ωc/Ω0 = 10, box 2.2–2.7 λ
  (λ = 2πVA/Ω0), 400–6400 particles/cell.
- Semi-implicit, damped field solver (App. C) with charge-conserving deposit plus a Poisson correction every step for the boundary.
- Linear theory of MRI/mirror/firehose with anisotropy in a Landau-fluid model (Figs. 2, 9).

## Numbers to compare against
- Ideal MRI maximum growth 0.75 Ω0 (Eq. 18); anisotropy shifts the most unstable kz to longer wavelength (Eq. 17).
- Measured γ/Ω0 for β = 1536 (Sec. 3.1, Fig. 3): m = 1, 2, 3 → 0.18, 0.48, 0.30 before Torbit = 2.5; m = 1 → 0.55 after.
  (The text says m = 2 grows fastest early; the Fig. 3 caption says m = 1.)
- Mirror vs ion-cyclotron threshold crossover at (β∥, p⊥/p∥ − 1) ≈ (187, 0.039) (Sec. 3.1).
- p⊥/p∥ up to 3.16 at Torbit = 4.82 (Fig. 1) and 18.7 at 5.39 (Sec. 3.2). The x-averaged value in the current
  sheet before onset is ≈ 1.3 for all four β (Fig. 5, Sec. 4.1). After reconnection it is 0.3–0.4 between islands (Sec. 3.2).
- Reconnection onset at Torbit ≈ 5.4–5.5 for β = 1536 (Fig. 4). Late fields at Torbit 7.95: Bpol/B0 ~ 330, Btor/B0 ~ 50 (Sec. 3.2).
- Field amplification Bf/Bi ~ 10^2 for all β (Fig. 10, Sec. 4.2). The final β, defined with the initial pressure, is ~ 0.01–1 (Sec. 4.2).
- α = wxy/p (Eqs. 20–21): 0.1–1 in the active phase and 1e-3–1e-2 in the quiescent phase (Fig. 11).
  With the initial pressure as the norm, α is 0.1–100 (Sec. 5).
- Spectrum ∝ ε^-1 with εmax ~ 10^2 mc² (Figs. 7, 12). The box-gyroradius estimate gives εmax/mc² ≈ 174, 150, 100, 75
  for β = 96…6144 (Eq. 22, Sec. 4.2).

## Assessment
- Shown by simulation: anisotropy growth consistent with CGL (Eqs. 15–16); mirror-mode signatures (density
  anticorrelated with |B|, spots of gyroradius size); channel flows that thin into current sheets, then reconnect; power-law tails across a 4-point β scan.
- Shown by theory: the linear MRI growth maps with anisotropy (Figs. 2, 9) and the pair-plasma reduction to MHD (App. A).
- Stated or argued, not tested: that the anisotropy causes the fast onset. The support is an order-of-magnitude use of
  the Chen–Palmadesso tearing–mirror growth rate (Eq. 19); there is no control run without anisotropy. That stochastic
  "Fermi-reconnection" acceleration works in large disks, inferred from energetic-particle voids inside islands
  (Fig. 13, Sec. 5). The relevance to Sgr A*. The measured growth rates fall below 0.75 Ω0, which the paper attributes to
  time-dependent anisotropy without checking.
- Stated limits: 2D only, with 3D transport and saturation expected to differ (Sec. 5); small box, so εmax is set by the box
  and the final island size may be too; Ωc/Ω0 = 10 exaggerates the anisotropy (Sec. 3.1, Sec. 5); pair plasma;
  non-inertial-frame terms and vertical gravity neglected; ~8–9 orbits.
- Left open: 3D runs, realistic Ωc/Ω0, ion–electron plasmas, larger boxes, and whether α stays large in 3D.
