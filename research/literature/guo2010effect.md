# Guo & Giacalone (2010): The effect of large scale magnetic turbulence on the acceleration of electrons by perpendicular collisionless shocks

- Ref: Astrophys. J. 715, 406 (2010) (volume and page from the DOI and from the entry in Guo & Giacalone
  2015); DOI 10.1088/0004-637X/715/1/406; arXiv:1003.5946; peer-reviewed (arXiv comment: accepted by ApJ).
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
A shock that moves through large-scale magnetic fluctuations accelerates suprathermal electrons quickly:
meandering field lines bring the fast electrons back to the front several times, and ripples of the front
mirror them. The effect is strongest when the mean field is perpendicular to the shock normal.

## Method
- 2-D hybrid run (kinetic ions, massless electron fluid) with turbulence imposed upstream, then test
  electrons in the stored fields (non-relativistic, Bulirsch-Stoer, step 5e-5 to 0.1 Omega_ce^-1).
- Plane x-z, three vector components; piston method (plasma injected at x = 0, reflected at x = L_x);
  periodic in z. Box 400 x 1024 c/omega_pi, cells 0.5, time step 0.01 Omega_ci^-1, 40 ions per cell,
  resistivity 1e-5 (4 pi / omega_pi).
- Fluctuations: 1-D, B(z,t), equal mix of circularly polarized Alfven waves along the mean field,
  spectrum P(k) ~ 1/[1 + (k L_c)^(5/3)], L_c = L_z = 1024 c/omega_pi; variance dB^2/B_1^2 = 0.3
  (also 0.1 and 0.5); injected continuously at x = 0. Mean field in the simulation plane.
- Inflow Mach number 4.0 in the simulation frame, about 5.6 in the shock frame (average); shock speed
  about 1.6 v_A1 in the simulation frame. beta_e = beta_i = 0.5. Sonic Mach number: not stated.
- Mean theta_Bn = 90 deg (main case), 75 and 60 deg. c/v_A1 = 8696; Omega_ce/Omega_ci = 1836.
- Electrons: 1.6e6, shell of 100 eV (30.7 U_1 = 5.7 v_the), released at 70 Omega_ci^-1; run ends at 120.

## Numbers to compare against
- Mean theta_Bn = 90 deg, variance 0.3: energies up to 200-300 times the initial 100 eV within
  50 Omega_ci^-1; flat energy-flux spectrum from 0.1 to 0.7 keV; slope about -3 above 1 keV (Fig. 4).
- Number fraction and highest energy both fall as the mean angle goes 90 -> 75 -> 60 deg (Fig. 4; the
  text gives no values).
- Variance 0.1, 0.3, 0.5: spectra alike at low energy, different above 2 keV, hardest for 0.5 (Fig. 5).
- Energetic electrons reach about 100 c/omega_pi upstream at z = 200 c/omega_pi, but there is no
  significant upstream flux at z = 800 c/omega_pi (Fig. 6, at 100 Omega_ci^-1).
- One trajectory: several reflections between ripples (labelled a-b, c, d, f-g, h, j), each with a jump
  in energy, while the shock hardly moves (Fig. 3).

## Assessment
- Shown by simulation: fast energization at a perpendicular mean field with turbulence; patchy energetic
  electron density along the front, which the authors tie to the local theta_Bn (Fig. 2; Sect. 3).
- Shown for single particles only: mirroring between ripples (Fig. 3). The shares of field-line meandering
  and of ripples in the total gain are not separated; there is no run without imposed turbulence.
- The front is irregular from ion-scale ripples up to the scale of the imposed fluctuations (Sect. 3); the
  wavelengths are not quantified. Electron gyroradius: not stated.
- Stated limits: 2-D with an ignorable coordinate, so electrons stay tied to their field line and stop
  gaining once the line is downstream (Sect. 3); no electron-scale structure (hybrid); the full range from
  coherence scale to electron scale is out of reach (Sect. 2).
- Agrees with Burgess (2006) on the plateau and on lower efficiency for higher initial energy (Sect. 3).
- Stated, not shown: relevance to the injection problem (Sect. 4). Initial energy is suprathermal (halo).
- The shock is planar on average; no globally curved front. Bale et al. (1999) and Pulupa & Bale (2008)
  are cited for large-scale irregular shock surfaces at interplanetary shocks.
- Doubt: the in-plane direction of the mean field for the oblique cases is not spelled out.
