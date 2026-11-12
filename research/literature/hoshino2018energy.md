# Hoshino (2018): Energy Partition between Ion and Electron of Collisionless Magnetic Reconnection

- Ref: Astrophys. J. Lett. 868, L18 (2018); doi:10.3847/2041-8213/aaef3a; arXiv:1811.03226
- Read: 2026-11-12, full (arXiv v1); topics: (none yet)

## Claim
Following the plasma in one reconnecting flux tube (bounded by Az contours) separates heating into
two stages. A short nonadiabatic stage happens while the tube crosses the diffusion region. An
adiabatic stage (T V^(2/3) ~ const) follows as the merged tube shrinks. The ion/electron ratio of
temperature increments is set in the nonadiabatic stage and scales as (mi/me)^(1/4). An Ohmic-heating
estimate over species-dependent diffusion regions reproduces this scaling.

## Method
- 2D PIC; anti-parallel Harris sheet (no guide field), Ti0 = Te0, background 5 % of peak density;
  periodic in x, conducting walls at y = +-40 lambda; box 160 x 80 lambda on 6400 x 3200 cells
  (lambda = 40 cells, ion inertial length 56.5 cells, electron 2.8 cells); 2.2e10 particles
  (538 x 2 per cell). mi/me = 400, plus runs at 6.25, 25 and 100.
- Lagrangian diagnostic: in 2D, frozen-in flux means Az is advected, so a tube between two Az
  contours follows the flow. The paper tracks the particle count in the tube, its volume V (the area
  between the lines) and a tube-averaged "temperature". That temperature is the second moment in the
  simulation frame, so it includes the bulk-flow energy.
- Theory (Sec. 3): Q_j = E J_j S_j over each species' diffusion region. It uses a meandering width
  d_j and a Speiser length Delta_j, with an inertial collision frequency nu_j = v_out/Delta_j.

## Numbers to compare against
- Particle number in the tube conserved to about 1 %; the largest slippage is at
  t Omega_ci ~ 27-33 (Fig. 2a).
- Nonadiabatic stage t Omega_ci ~ 27-33: Te about doubles and Ti grows about 4x at nearly fixed V
  (Fig. 2b-c).
- Adiabatic stage t Omega_ci > 35: T_j proportional to V^-(gamma-1) with gamma = 5/3 (Fig. 2b-c,
  dashed).
- dTi/dTe: slope 4.31 for the main tube (Fig. 3a); 4.06 +- 0.26 over several tubes (Sec. 2);
  mi/me = 400.
- Fit over mi/me = 6.25-400: dTi/dTe = (mi/me)^(0.237 +- 0.005) (Eq. 1, Fig. 3b); a ratio of 1
  is assumed for pairs.
- Model: dTi/dTe = (mi Ti0 / me Te0)^(1/4) (Eq. 2) = 4.47 for mi/me = 400, Ti0 = Te0 (Sec. 2).

## Assessment
- Shown by simulation: two-stage heating of one tracked tube and the match to adiabatic compression
  afterwards. The mass-ratio scaling is a power-law fit to four mass ratios plus the assumed pair
  point.
- The theory is a scaling estimate. It assumes the inflow speeds equal the thermal speeds (valid here
  because Ti0 = Te0 and the background is as hot as the sheet) and equal J*Delta for ions and
  electrons. Separatrix heating enters only implicitly, through beams from the diffusion region.
- Stated or conjectured: the final partition is set only in the nonadiabatic stage, because Ti/Te
  stays fixed in the adiabatic stage. The model is said to suit plasmoid-dominated reconnection, and
  the paper suggests it applies to heavy ions.
- Definitions matter. The temperature includes bulk flow, and the method follows flux tubes, while
  prior work compares upstream with exhaust (Shay et al. 2014; Haggerty et al. 2015).
- Limits (stated): 2D and anti-parallel only. A guide field changes the dissipation (Drake & Lee
  1977). 3D effects (lower-hybrid drift, flux ropes) are excluded. The box is only several
  diffusion-region sizes, and the periodic x boundary makes the late-time plasmoid oscillate.
