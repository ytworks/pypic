# Trotta & Burgess (2019): Electron acceleration at quasi-perpendicular shocks in sub- and supercritical regimes: 2D and 3D simulations

- Ref: Mon. Not. R. Astron. Soc. 482, 1154 (2019) (volume and page as printed in Pohl et al. 2020,
  ref. 307; the arXiv metadata give none); DOI 10.1093/mnras/sty2756; arXiv:1808.00812; peer-reviewed
  according to the DOI, the arXiv version is the submitted preprint.
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
Ion-scale ripples of the front exist only above a critical Mach number (about 3.5), and only then are
suprathermal electrons accelerated well beyond adiabatic mirror reflection, also at theta_Bn down to
80 deg. What counts is the local theta_Bn, not the mean one. With the front resolved in 3-D the gain is
slightly larger than in 2-D, while 2-D shows partly artificial trapping.

## Method
- 2-D and 3-D hybrid runs (code HYPSI, CAM-CL scheme; massless adiabatic electron fluid), fields stored at
  every step; then test electrons (4th-order scheme, adaptive step 0.1 to 6e-6 Omega_ce^-1, tricubic
  interpolation in space, linear in time).
- Shock by injection against a reflecting wall; upstream flow along the normal; periodic in y and z.
  Upstream field in the x-y plane (in-plane for the 2-D runs).
- Box 60 x 60 d_i (2-D) and 60 x 60 x 20 d_i (3-D); cells 0.25 d_i; ion step 0.01 Omega_ci^-1; more than
  100 ions per cell; small resistivity. Shocks chosen not to reform.
- M_A = 2.9, 3.5, 4.6, 5.1, 5.6, 6.6 at theta_Bn = 87 deg; theta_Bn = 80, 85, 87 deg at M_A = 2.9 and 6.6
  (Sect. 3.3). Frame of the Mach number, beta, sonic Mach number, run time: not stated.
- c/v_A = 5000 (v_A = 60 km/s). Electrons: about 1e6, kappa distribution (kappa = 4, T = 1e5 K) built from
  mono-energetic shells released 5 d_i upstream; collected 7 d_i upstream and 5 d_i downstream.

## Numbers to compare against
- Spread of the local field-to-normal angle in a 1 d_i thick layer at the front (Gaussian fit): 1.4 deg
  at M = 2.9, 15.2 deg at M = 6.6, both at mean theta_Bn = 87 deg (Fig. 2, Sect. 3.1).
- M = 2.9: only moderate gain, as in adiabatic theory, in 2-D and 3-D (Fig. 3). M = 3.5: longer tail,
  maximum about 10 keV (Fig. 5). M = 6.6: plateau-like spectrum over two decades in energy (Fig. 4).
- M = 2.9, theta_Bn = 80 deg: maximum below 1 keV (Fig. 6). M = 6.6: flattened spectra at all three
  angles, also at 80 deg (Fig. 7).
- 3-D against 2-D at M = 6.6: about the same maximum energy, shallower fall-off in 3-D (Fig. 4).
- Trajectories: 100 eV -> 2e3 eV at a 2-D M = 2.9 shock by a long drift along z, called an artefact of
  2-D (Fig. 10); 100 eV -> 3.4e4 eV at a 3-D M = 6.6 shock in several drift phases (Fig. 11).
- Among reflected electrons that end above 1 keV (100 eV shell, 3-D, M = 6.6), the more energetic ones
  reached deeper into the shock before returning (Fig. 12).
- Critical Mach number for surface fluctuations: about 3.5 (Sect. 5).

## Assessment
- Shown by simulation: the step from sub- to supercritical in the electron spectra; the broadening of the
  useful theta_Bn range; the measured spread of the local angle; first hybrid plus test-particle study
  with a 3-D front.
- Interpretation supported by single trajectories and Fig. 12, not by a statistical measure: ripples act
  as scatterers and traps, an electron meets several ripples in one encounter, and it may arrive where
  the front is locally perpendicular and leave where it is less so (Sect. 5).
- Stated limits: no reformation or other large-scale non-stationarity; no upstream turbulence (its share
  against ripples is left for later work); hybrid fields without electron scales; electrons start
  suprathermal, so thermal injection is not addressed; confirmation by full PIC is asked for (Sect. 5).
- Against Burgess (2006): the 2-D plateau is accentuated by trapping in ripples that are artificially
  coherent along z. Compared with Amano & Hoshino (2009): similar reflect-and-scatter picture, but there
  the scatterers are electron-scale waves.
- The mean front is planar. Curvature appears only in the introduction (Krauss-Varban & Burgess 1991:
  electrons travel far during reflection, so the curvature of the shock matters).
- Maximum energies above M = 3.5 are not given as numbers. Pohl et al. (2020) describe this paper as
  being about "superluminal" shocks; the paper itself says supercritical.
