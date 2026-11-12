# Zeković, Spitkovsky & Hemler (2025): SLAMS-propelled electron acceleration at high-Mach number astrophysical shocks

- Ref: ApJ 988, 40 (2025); DOI 10.3847/1538-4357/ade05f; arXiv:2408.02084v2 (accepted version).
  Peer-reviewed publication.
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
At parallel shocks with M_A >~ 80, cosmic-ray-driven waves steepen into SLAMS that amplify the
transverse field >~ 10 times. At fast shocks this makes the field intermittently superluminal for
magnetized electrons and raises their DSA threshold to the momentum where their gyroradius reaches the
SLAMS scale. Below it, electrons gain energy by bouncing between the shock and the next approaching
SLAMS maximum ("quasi-periodic shock acceleration", QSA), with a spectrum steeper than DSA:
f ~ p^-4.7 ... p^-5.7.

## Method
- PIC, TRISTAN-MP, upstream frame, moving piston with a graded-field wall (App. A). B0 along the flow
  (theta_Bn = 0; in-plane in 2-D). M_s = M_A in all runs, Ti = Te (beta not stated).
- Table 1, 1-D: run 1: m_i/m_e = 32, M_A = 80, v_sh = 0.133 c, 217 /omega_ci; run 2: 32, 40, 0.067 c,
  158; run 3: 32, 20, 0.033 c, 221; run 4: 100, 80, 0.133 c, 82; run 5: 32, 200, 0.33 c, 36;
  run 6: 32, 300, 0.267 c, 234. 10 cells per c/omega_pe, 1024 particles per cell, up to 8.1e5 cells.
- 2-D: run 7: 32, 80, 0.267 c, 61 /omega_ci, 8.9e4 x 11200 cells (width ~5 ion Larmor radii in B0);
  run 8: same, 204 /omega_ci, 1.1e5 x 2800 cells (~1.2 Larmor radii). 5 cells per c/omega_pe,
  32 particles per cell. sigma_e = 2.5e-5 ... 3.6e-4; 16 current-filter passes.
- Analytic QSA model (Sect. 4.3, App. C); test-particle follow-up cited as Hemler et al. (2024, subm.).

## Numbers to compare against
- B_perp/B0 near the shock: 3-5 (M_A = 20), ~7 (40), >~10 (80), >~20 (200), >~30 (300) (Table 1).
- SLAMS scales 500, 1000, 2000 c/omega_pe ~ 1, 2, 4 r_L,sh, with r_L,sh = M_A sqrt(m_i/m_e) c/omega_pe;
  initial SLAMS move upstream at v_sh/4 - v_sh/2; n_CR/n0 ~ 0.01-0.05 (Sect. 3).
- Escape along the field needs v_par > v_g ~ (B_perp/B0) v_sh; superluminal if v_sh B_perp/B0 > c
  (Sect. 4.2), which holds in all the high-M_A runs.
- DSA threshold: p_res,e/(m_e c) = [lambda_SLAMS/(2 pi r_Lsh,i)] (v_sh/c) (m_i/m_e) (B/B0) (Eq. 1);
  ~30 in run 1. Real mass ratio, v_sh ~ 1e4 km/s, M_A >~ 300: p_e/(m_e c) >~ 1000; for (1-4)e4 km/s:
  0.6-3.5 GeV (resonant growth) or 0.6-16.5 GeV (non-resonant) (Sects. 4.4, 5 iii).
- Slope: delta_e ~ 6 with the initial SLAMS, 5.0-5.5 later (Sect. 4.1); 5.6 quoted for runs 1 and 8
  (Sect. 4.3). Model (Eq. 7) gives 4.7-5.7 for v_sh = 0.033-0.133 c and P_ref = 0.7-0.9; radio index
  alpha = (delta_e - 3)/2 = 0.85-1.35.
- QSA time ~ one SLAMS advection time, lambda_SLAMS/v_sh ~ 0.01 T_end (Eq. C8, Sect. 4.1).
- No injected fraction or K_ep is given; ions reach only p_i/(m_i c) ~ 2 (Sect. 4.4).

## Assessment
- Shown by simulation: SLAMS, shocklets and steep electron tails in 1-D and 2-D at m_i/m_e = 32;
  bouncing trajectories whose p(t) follows the QSA curve, not SDA (Fig. 5). In 2-D the SLAMS-driven
  cavities corrugate the front, yet trajectories and spectra stay close to 1-D (Sects. 3, 4.1).
- Model: P_ref is taken from another work and v_SLAMS is chosen to recover delta_e ~ 5.6, so the
  agreement of slopes is a consistency check, not a prediction.
- Conjectured: flattening to p^-4 above p_res and a GeV break at SNRs; no run reaches the transition
  (Sect. 5 v). The m_i/m_e = 100 run only begins to form a tail at its end (App. B).
- Stated limits: v_sh >= 0.033 c, reduced mass ratio, 2-D width of 1.2-5 gyroradii; the cavities
  reach the box width by the end of run 8. 3-D electron-proton runs and wider boxes are planned
  (Sect. 5 v).
- Geometry: B_perp/B0 >~ 10 is a local obliquity >~ 84 deg; for magnetized electrons the parallel
  shock resembles a quasi-perpendicular one (Sects. 4.4, 5 iv), linked to Bell et al. (2011).
  Time variability: strong reformations and a variable SLAMS envelope (Sect. 3).
- Disagreements: runs with M_A < 40 (e.g. Park et al. 2015) saturate at B_perp ~ 2-3 B0 and miss this
  regime; the steepening needs no downstream "postcursor" (Caprioli et al. 2020). Against
  gupta2025speed: here the DSA threshold at high M_A is a gyroradius condition, far above u = v_sh.
