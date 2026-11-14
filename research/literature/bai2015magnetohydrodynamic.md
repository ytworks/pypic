# Bai, Caprioli, Sironi & Spitkovsky (2015): Magnetohydrodynamic-particle-in-cell method for coupling cosmic rays with a thermal plasma: application to non-relativistic shocks

- Ref: ApJ 809, 55; DOI 10.1088/0004-637X/809/1/55; arXiv:1412.1087. Peer reviewed. Read in the arXiv
  version.
- Read: 2026-11-14, arXiv PDF text: abstract, the formulation and its range of validity, the injection
  recipe, discussion; topics: qpar-electron-conditions

## Claim
The thermal plasma can be treated as an MHD fluid and the cosmic rays as particles, coupled through the
cosmic-ray current and a cosmic-ray Hall term. The method reaches scales that full PIC and hybrid PIC
cannot; particle injection has to be prescribed.

## Method
- MHD-PIC, valid on scales much larger than the ion skin depth with n_CR ≪ n_i.
- Injection: a fixed fraction η of the fluid crossing the shock becomes particles; typical value 2 × 10⁻³.

## Numbers to compare against
- η = 2 × 10⁻³ (ions). `caprioli2018diffusive` states that this recipe "does not apply to oblique shocks".

## Assessment
- On full PIC: it "has to resolve the microscopic plasma scale … Even with artificially reduced
  ion-to-electron mass ratio, properly capturing the physics at both scales still demands tremendous
  computational cost".
- On its own limit: "the main disadvantage of the MHD-PIC approach … is that injection of CR particles has
  to be handled with ad hoc prescriptions. However, this may be remedied by properly designing controlled
  experiments to calibrate the injection recipes with hybrid-PIC simulations". η is "our only parameter
  for particle injection, which encapsulates our ignorance of the physics of particle injection".
- Electron-scale streaming instabilities "are not captured in either MHD or hybrid-PIC methods".
- The formulation is "equally applicable to sufficiently energetic CR electrons (e.g., whose energy is
  higher than energy of thermal ions)".
- For this topic: the neighbouring method that would carry electrons to long times, and the statement of
  what it needs from a kinetic run: a calibrated injection recipe. For electrons such a recipe does not
  exist in what we read; an electron injection rate at event conditions would be its input.
  `galishnikova2026self` uses this class of method with prescribed electron injection.
