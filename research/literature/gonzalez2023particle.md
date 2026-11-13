# Gonzalez, Innocenti & Tenerani (2023): Particle-in-cell simulations of Alfven wave parametric decay in a low-beta plasma

- Ref: J. Plasma Phys. 89(2), 905890208 (2023) (Crossref); DOI 10.1017/S0022377823000120;
  arXiv:2301.07646; peer-reviewed.
- Read: 2026-11-13, skimmed (ar5iv HTML, set-up, figures and conclusions); topics:
  shock-alfvenic-turbulence

## Claim
1-D fully kinetic runs of the parametric decay of a parallel, circularly polarized Alfven wave at low
beta: the pump decays into an ion-acoustic mode and a backward Alfven wave; protons gain about 50 % of
the pump energy as internal energy and electrons about 10 %. A fast-type compressive front at the
steepened edge of the backward wave drives a backward field-aligned proton beam at v_A.

## Method
- Codes: ECsim (semi-implicit PIC), VPIC (explicit PIC), CAMELIA (hybrid, reference) (Table 1).
- 1-D, L_x = 102.4 d_i; dx = 0.025 d_i (ECsim), 0.0125 d_i (VPIC), 0.1 d_i (hybrid); 500 ppc (PIC),
  10^4 (hybrid). m_i/m_e = 400, v_A/c = 0.05.
- Pump: circularly polarized, delta B/B0 = 0.5, mode m0 = 4, k0 d_i = 0.2454 (weakly dispersive).
- Runs: TiTe1 (beta_i = beta_e = 0.125), TiTe1B (beta_i = beta_e = 0.078), TiTe4 (T_i/T_e = 4).

## Numbers to compare against
- Decay sets in at t Omega_ci ~ 150; saturation with delta n_rms/<n> ~ 0.2 near t Omega_ci ~ 225;
  then nearly complete reflection (Fig. 1).
- Ion-acoustic mode damps over t Omega_ci = 200-400; steepened backward front at t Omega_ci ~ 400
  (Fig. 5).
- Energy: Alfvenic energy drops by about 60 % of the pump; protons +~50 %, electrons +~10 % of the pump
  energy (Fig. 4). Electron heating nearly isotropic, at a roughly constant rate; electron internal
  energy still doubles after the decay (t Omega_ci > 500). No electron beam.
- Growth rate: set by total beta "in agreement with fluid theory" (Fig. 2); the value appears only as a
  figure legend, not in the text.
- Hybrid and PIC agree on the density growth (Fig. 1). VPIC gives slightly higher final T_e than ECsim;
  energy error < 1 % in both.

## Assessment
- Directly usable for the design: at delta B/B0 = 0.5, beta ~ 0.25, k0 d_i ~ 0.25 the decay takes
  ~150 Omega_ci^-1 to start. A shock run that crosses a 10-100 d_i upstream in less than that is
  little affected by parametric decay of a coherent pump; a longer run is not.
- Electrons gain ~10 % of the wave energy with no shock: a turbulence-only null run must subtract this.
- Limits: 1-D, k parallel only (no perpendicular cascade, no oblique decay); m_i/m_e = 400; low beta
  only; a single monochromatic pump. Multi-D is announced, not shown.
