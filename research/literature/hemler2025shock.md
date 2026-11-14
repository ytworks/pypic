# Hemler, Spitkovsky & Zekovic (2025): Shock Acceleration via Idealized Short Large-amplitude Magnetic Structures

- Ref: ApJ 988, 144; DOI 10.3847/1538-4357/adb7e3. Peer reviewed.
- Read: 2026-11-14, abstract only (full text not obtained); topics: qpar-electron-conditions

## Claim
Test particles in analytically prescribed SLAMS, built to resemble those of PIC simulations, reproduce the
electron spectra of PIC. Electrons with gyroradius well below the SLAMS size are accelerated by shock
drift and quasi-periodic shock acceleration; larger-gyroradius electrons by DSA.

## Method
- Prescribed magnetic structures moving with the flow; test electrons and ions.

## Numbers to compare against
- Persistent SLAMS for M_A >~ 20, driven by the Bell instability.
- Broadband Bohm diffusion for both magnetized (r_L << lambda) and unmagnetized electrons.
- Quasi-periodic shock acceleration is "vital" when the SLAMS are superluminal.

## Assessment
- Abstract only.
- The authors state that "the test-particle method can be used to extrapolate PIC results to longer
  spatial and temporal scales". This is the same strategy as step 2 of this topic's idea (measure in PIC,
  extend by a transport model), for prescribed structures at M_A >~ 20.
- Differences: the structures are prescribed, not taken from a run at event conditions; no confrontation
  with a heliospheric event in the abstract. Follows zekovic2025slams (on
  `research/spherical-shock-injection`).
- M_A >~ 20 for persistent SLAMS agrees with the absence of mature SLAMS at the PSP shock
  (jebaraj2026wave, M_A 7-7.5) and at a planar interplanetary shock (raptis2026compressive, M_A 5.4).
