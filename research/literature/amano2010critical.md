# Amano & Hoshino (2010): A Critical Mach Number for Electron Injection in Collisionless Shocks

- Ref: Phys. Rev. Lett. 104, 181102 (2010); DOI 10.1103/PhysRevLett.104.181102; arXiv:1004.3838
- Read: 2026-11-12, full; topics: (none yet)

## Claim
Shock drift acceleration mirror-reflects electrons upstream as a beam with a loss cone. Cyclotron resonance at
the loss-cone edge lets this beam drive whistler waves that scatter the energetic electrons themselves, which
injects them into DSA. Waves grow against thermal cyclotron damping above a critical Alfvén Mach number that
scales as cosθBn (mi/me)^(1/2) βe^(1/2). The controlling parameter MA/√βe does not depend on density or field
strength. The authors read this as why electron injection is poor at heliospheric shocks and efficient at young SNRs.

## Method
- Analytic linear theory only, no new simulation. Reflected beam speed Vb ≃ 2Vs/cosθBn (in the HTF picture).
- Rough conditions: self-generated Alfvén waves need Eq. 1 (from the authors' 2007 ApJ paper); whistlers at the
  loss-cone resonance need |Vb| ≳ ve, which gives Eq. 2.
- Refined: cold ring-beam electrons (density η = nb/ni, ring speed Vr from mirror reflection with normalized HTF
  potential φ̃ = 0.2 and a loss cone from Bmax that includes upstream δB = B0, compression r = 6). Cold-plasma
  R-mode dispersion, perturbative growth rate (∝ η^(1/3) Vr^(2/3)) minus thermal cyclotron damping with βi = βe.
  The critical MA is the root of γ_tot = 0 for η = 1e-3, 1e-4, 1e-5 (Fig. 2).
- Evidence check: Geotail bow-shock statistics (Oka et al. 2006) and Gosling et al. (1989).

## Numbers to compare against
- Eq. 1: MA ≳ (cosθBn/2)(mi/me). Eq. 2: MA_inj = (cosθBn/2)(mi/me)^(1/2) βe^(1/2), lower by ~43/√βe (p. 2).
- Numerical/analytic ratio MA_inj*/MA_inj ≈ 2 or slightly less for βe ~ 1, θBn ≳ 60°, 1e-5 ≲ η ≲ 1e-3;
  ≈ 1–3 over a wider survey (Fig. 2, pp. 3–4). Weak dependence on φ̃ and loss-cone choice at θBn ≳ 45° (p. 3).
- Vr/Vb ≈ 4.7 at θBn = 30°, ≲ 1 at quasi-perpendicular shocks (p. 3).
- γ_tot/Ωci ≈ 21 at MA = 10, βe = 1, θBn = 80°, η = 1e-4, so growth length ≪ Vs/Ωci (p. 4).
- Oka et al. whistler critical Mach number (α/2) cosθBn (mi/me)^(1/2), α = 27/16; spectral index flattens slightly above it (p. 3).
- SNR condition: MA/√βe ≃ 68 (Vs/3000 km/s)(Te/10 eV)^(-1/2) ≳ 30 (Eq. 3).
- Assumed inputs: φ̃ = 0.2 (observed 0.1–0.2), r = 6, δB = B0, injection rate η ~ 1e-4 from the authors' 2007 1D-PIC-based estimate (p. 4).

## Assessment
- Shown (theory): growth-vs-damping threshold curves for a model ring-beam (Fig. 2); a closed-form scaling (Eq. 2)
  that the numerical root follows within a factor of 1–3.
- Supported (observations, after the fact): consistent with the Geotail result that electron spectra harden
  above roughly twice the whistler critical Mach number, and with an enhanced-suprathermal high-MA (≈ 41)
  quasi-parallel event. Not tested by kinetic simulation here.
- Stated/conjectured: the injection is called "robust"; the paper says nothing prevents it at astrophysical
  shocks and that SNRs satisfy Eq. 3 (pp. 3–4). These go beyond a linear necessary condition.
- Stated limits: linear analysis only (whistler cascades may change efficiency); cold ring-beam idealization;
  η a free parameter; quasi-parallel estimate uncertain by a factor of order unity and "can be studied only"
  with self-consistent simulations; loss cones seen upstream imply no whistler scattering there, so electron
  data need re-analysis vs MA; CR back-reaction (amplified field, pre-heating, upstream deceleration) left open.
- Open: self-consistent kinetic test of the threshold and of the scattering it implies.
