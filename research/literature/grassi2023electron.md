# Grassi, Rinderknecht, Swadling, Higginson, Park, Spitkovsky & Fiuza (2023): Electron Injection via Modified Diffusive Shock Acceleration in High-Mach-number Collisionless Shocks

- Ref: ApJL 958, L32; DOI 10.3847/2041-8213/ad0cf9 (report number LLNL-JRNL-859893). Peer reviewed. Read
  in the accepted manuscript (OSTI 2293646).
- Read: 2026-11-14, abstract and parts of the body (text extracted from the PDF; many digits lost, so only
  statements checked in context are quoted); topics: qpar-electron-conditions

## Claim
In high-Mach-number shocks electrons are injected by scattering in the ion-Weibel turbulence of the shock
transition: a "modified DSA" that works before the electrons can cross the whole shock. The model agrees
with PIC simulations and with laser experiments at the National Ignition Facility.

## Method
- PIC of initially unmagnetized counterstreaming plasmas, mass ratio 128-512.
- Electron diffusion coefficient from the mean square displacement of tracked electrons,
  D = <dx^2> / (2 dt); energy gain rate against energy.

## Numbers to compare against
- Energy gain rate dW/dt proportional to energy^alpha with alpha = 1.2 +- 0.1.

## Assessment
- This is the "LLNL-JRNL-859893" paper of the plan: an electron diffusion coefficient taken from PIC.
  It is measured in Weibel turbulence of an unmagnetized shock, not in the precursor of a magnetized
  quasi-parallel shock, and not at the conditions of any heliospheric event.
- For this topic: prior use of the mean-square-displacement method for electrons in PIC, next to
  gupta2025speed. The factor in D (1/2 here) must be fixed against the convention of gupta2025speed
  before values are compared (null test N7).
- The digits lost in extraction mean any further number needs the published version.
