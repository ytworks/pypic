# M. Hoshino (2020): Stabilization of Magnetic Reconnection in Relativistic Current Sheet

- Ref: ApJ 900, 66 (2020); DOI 10.3847/1538-4357/aba59d; arXiv:2006.15501
- Read: 2026-11-12, full; topics: (none yet)

## Claim
In a relativistically hot pair Harris sheet, the collisionless tearing growth rate does not keep rising with the
drift speed βc of the current carriers. It follows the standard (Zelenyi–Krasnoselskikh) β^{3/2} law up to β ≈ 0.5 and
then falls once Γβ > 1. The paper derives an extended linear theory: the resonant (Landau) conductivity grows by Γβ,
through relativistic inertia and a retarded-time effective collision frequency kc/Γβ. The result is
γτc ∝ kλ(1 − k²λ²)β^{3/2}/Γβ. Reconnection in a current sheet with ud → c (e.g. a pulsar striped wind) would then be suppressed.

## Method
- 2D PIC, periodic in x, conducting walls at y = ±16λ; box 64λ × 32λ, 1600 × 800 cells, 4.5e8 particles.
- Pair Harris sheet with a drifting relativistic Maxwellian (Eq. 6): T/mc² = 10, σ = 20, λ = 25 cells,
  background density 5 % of the peak. β is scanned up to Γβ ~ 10³ (Fig. 2).
- Growth rates come from exponential fits to Fourier modes of By averaged over ±4 cells about the sheet. The 3–4
  largest rates around kλ ≈ 0.39 are averaged, with 10 % error bars (Sec. 2).
- Theory: linearized Vlasov equation along unperturbed orbits, meandering orbits taken as straight lines, vz² ≈ c²,
  and an energy principle that balances the adiabatic (∝ γ) and resonant (∝ γ²) energy changes (Eqs. 7–19).
- Physical motivation: an estimate of the drift speed in a pulsar wind current sheet (Eqs. 1–2).

## Numbers to compare against
- Standard theory: γτc ≈ (2√2/π) kλ(1 − k²λ²) β^{3/2} (Eq. 3). New theory: the same divided by Γβ (Eq. 19).
- Relativistic effective conductivity ∝ N e² Γβ mc² / (m k c T) (Eq. 13), against the nonrelativistic N e²/(m k vth) (Eq. 21).
- In the simulations the fastest modes are at kλ ≈ 0.3–0.5 (mode numbers 3–5) and kmin λ ≈ 0.098 (Sec. 2).
- PIC agrees with Eq. 3 for β < 0.5 and departs above it, with γ decreasing in β for Γβ > 1 (Fig. 1). Eq. 19 matches PIC
  up to Γβ ≈ 10³ (Fig. 2). The theory curve is scaled by an O(1) constant (Sec. 2).
- The growth rate does not depend on T for T/mc² > 1 (Sec. 2; agrees with Zenitani & Hoshino 2007).
- Drift-kink (from earlier work): γτc ≈ 1/(16 Γβ β) and most unstable kλ ≈ 1/(8Γβ β²) (Eq. 20).
- Pulsar wind: ud/c = r/(4κεΓwind rlc) (Eq. 2), so ud reaches c before the termination shock if εκΓwind < 10^9.
  With a growth time > 10 λ/c, ud/c > (5/2)Γwind/κ (Sec. 4).

## Assessment
- Shown: PIC growth rates against drift speed over three decades of Γβ − 1, and a closed-form theory that fits them
  after one O(1) rescaling. The theory relies on straight meandering orbits, a constant-A1 profile (corrected by hand for
  kλ ≪ 1) and vz ≈ c, so it is an order-of-magnitude derivation.
- Stated, not simulated: Eq. 19 also holds for a cold sheet with relativistic drift (argued in Sec. 3). The drift-kink
  mode is also suppressed (cited from earlier work). Tearing and drift-kink stay decoupled until the early nonlinear
  phase for Γβ ≫ 1. Reconnection ceases in pulsar winds (an order-of-magnitude argument depending on Γwind vs κ).
- Limits: 2D tearing plane only, with 3D competition against the drift-kink left open (Sec. 4); pair plasma; linear
  stage only; the β scan at one T and σ; how a β ~ 1 sheet forms and survives is open, and displacement current is
  neglected in Eq. 1 (Sec. 4).
- The review guo2023magnetic (Sec. 3.3) uses this result to say that plasmoid formation in collisionless pair
  plasma need not follow resistive-MHD plasmoid-instability theory.
