# M. Hoshino (2021): Nonlinear Explosive Magnetic Reconnection in a Collisionless System

- Ref: Phys. Plasmas 28, 062106 (2021); DOI 10.1063/5.0050389; arXiv:2105.13554
- Read: 2026-11-12, full; topics: (none yet)

## Claim
For thick pair-plasma current sheets (gyroradius rg well below the thickness λ), linear tearing is slow
(γτA ∝ (rg/λ)^{3/2}). It is followed by an explosive nonlinear phase whose growth rate scales with b1 = By/B0 and does
not depend on the initial thickness. Two kinetic effects raise the inertia resistivity at the X-line: the sheet shrinks
to about two gyroradii, and the meandering particles become magnetized by the reconnected field. This extends Galeev
et al. (1979). The explosive phase lasts tens of Alfvén times.

## Method
- 2D PIC, periodic in x, conducting walls at y = ±8λ; box 108λ × 16λ, 5376 × 800 cells, 1.2e11 particles. The
  particle count is chosen so that slow tearing can rise above thermal noise ∝ 1/√np (Sec. II).
- Pair Harris sheet with T = mc²/10, vA = vth = 0.447c, background 5 %. The drift β = 0.05, 0.1, 0.2 gives
  rg/λ = 0.11, 0.22, 0.45 (Sec. II).
- Seed: a localized X-type δAz with peak By/B0 = 1e-3 at the box centre. A run without the seed is said to be nearly
  the same. A filter outside the sheet damps waves for t/τA < 16.
- Diagnostics: Fourier spectra of Az near y = 0 (Fig. 2); linear rates fitted to exponentials (Fig. 3); reconnected flux
  Φ (Fig. 4a); instantaneous growth rate against b1 (Fig. 4b); maximum and normalized nonlinear rates against rg/λ
  (Fig. 4c); sheet thickness from a cosh⁻² fit at the X-line (Figs. 5–6).
- Theory: an energy principle that balances the adiabatic field energy against Ohmic loss in a diffusion region
  dx × dy (Eq. 1, App. A), with the collisionless inertia conductivity (App. B, Eq. B6).

## Numbers to compare against
- Linear theory Eq. 2 / A7, with a non-constant-A1 correction (Eq. A9). It is scaled by 0.73 to fit (Sec. III, Fig. 3).
  The fit windows are t/τA = [10,150], [10,70], [10,55], and the maximum growth is at kλ ≈ 0.2–0.5.
- The rapid phase starts at t/τA ≈ 520, 230, 110 for rg/λ = 0.11, 0.22, 0.45 (Fig. 2, Sec. II).
- Φ saturates at ≈ 0.1, which the paper ties to a reconnection rate of ~ 0.1 (Fig. 4a); kmin λ = 0.0584.
- Maximum nonlinear γτA ≈ 0.139, 0.151 and "1.85" as printed, almost certainly 0.185 since the Fig. 4b axis ends at 0.2.
  These are roughly thickness-independent (Fig. 4c).
- Minimum thickness λ(t)/λ ≈ 0.26, 0.45, 0.83 for rg/λ = 0.11, 0.22, 0.45, about 2 rg (Fig. 5). The shrinkage grows at the linear rate (Fig. 6).
- The Galeev rate ∝ (rg/λ)^{1/2} b1 (Eq. 3) disagrees with the PIC. The new γnl τA = (√π/4)(1 − k²λ²) b1/(kλ) (Eq. 4)
  does not. With the kλ + 2rg/λ correction it gives γnl τA/b1 = 2.33, 1.41, 0.77, "roughly consistent" with Fig. 4c (Sec. IV B).
- Onset condition b1 > rg k; critical b1 = 6.4e-3, 1.3e-2, 2.6e-2 against the seed 1e-3 (Sec. V).
- Finite-time singularity at Δt/τA = 4kλ/(b1(0)√π) (Eq. 5).

## Assessment
- Shown by simulation: two-stage growth in three thickness runs, a nonlinear rate nearly independent of rg/λ, shrinkage
  to ~ 2 rg, and the time of minimum thickness matching the time of steepest growth (Figs. 4–5).
- Theory: an order-of-magnitude model (straight meandering orbits, dy = 2 rg assumed, constant A1). Its stated range is
  2rg/λ < kλ < 1, while the measured mode kλ = 0.0584 lies outside it (Sec. IV B).
- Stated or conjectured: the mechanism carries over to ion–electron plasma (first at electron, then at ion scale); it
  does not work with a strong guide field (Sec. V); it links the linear phase to the late plasmoid phase.
- Limits: 2D, pair plasma, antiparallel field only, a Harris sheet with 5 % background, one T and vA, three thicknesses
  only. Thicker sheets are out of reach because linear growth is too slow against PIC noise. Saturation is left outside the scope.
