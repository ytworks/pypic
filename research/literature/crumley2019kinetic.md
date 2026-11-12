# Crumley, Caprioli, Markoff & Spitkovsky (2019): Kinetic simulations of mildly relativistic shocks I: particle acceleration in high Mach number shocks

- Ref: Mon. Not. R. Astron. Soc. 485, 5105 (2019) (volume and page as printed in the reference list of Xu et al.
  2020; the arXiv metadata gives only the DOI); DOI 10.1093/mnras/stz232; arXiv:1809.10809v1. Peer-reviewed
  journal article; the text read is the submitted version (v1), which may differ from the published one.
- Read: 2026-11-12, introduction, set-up, Sects. 4-6, conclusions and appendices A-B;
  topics: spherical-shock-injection

## Claim
Wide 2-D PIC runs of a mildly relativistic (v ~ 0.75 c), M_A ~ 15 shock show that a quasi-parallel shock turns from
Weibel- to Bell-mediated and then accelerates ions and electrons with E_max growing linearly in time, while a
superluminal quasi-perpendicular shock accelerates neither. The non-thermal electron energy fraction is ~5e-4 and
is measured correctly only if the box is wider than the gyroradius of the current-carrying ions.

## Method
- PIC, TRISTAN-MP, 2-D (x-y) with 3 field/momentum components; upstream rest frame, moving reflecting wall with
  Lorentz factor 1.5 (v_0 ~ 0.75 c); shock speed beta_s ~ 0.83.
- B_0 in the simulation plane in the main runs; out-of-plane control runs (Table 1, App. B).
- m_i/m_e = 64 (fiducial), 16 and 160 (App. A). sigma_0 = 0.007 (downstream frame); M_A ~ 15 (shock speed over
  upstream Alfven speed, Eq. 2); M_s ~ 15 (definition not spelled out); beta_p = 1; T_i = T_e.
- theta_B ~ 10 deg (15 deg in the downstream frame) and 55 deg (65 deg); superluminal above ~34 deg.
- Fiducial box width 16000 cells = 4000 c/omega_pe = 500 c/omega_pi; narrower runs 25-200 c/omega_pi.
  4 cells per c/omega_pe; dt = 0.1125/omega_pe; 4 particles per cell in total; current filtered 32 times per step.
- Run to omega_pi t = 4148; not given in ion gyro-times (my estimate from M_A and beta_s: ~230/Omega_ci).

## Numbers to compare against
- Bell growth time ~60/omega_pi for a non-thermal ion fraction ~0.01 (Eq. 6); waves appear at omega_pi t ~ 560;
  Bell-mediated from ~1500, about 25 growth times (Sect. 3). Upstream waves right-handed (Fig. 4).
- Field: y-averaged ~2 B_0 upstream, > ~5 B_0 at cavity walls near the shock (Sect. 3, Fig. 6; conclusion i).
- Electrons: f(p) ~ p^-4.2 downstream; eps_e ~ 5e-4 of the incoming ion kinetic energy, integrated above
  p_inj = 2 gamma_0 beta_0 m_i c in the downstream frame (Eq. 5, Sect. 4, Fig. 7); number ratio of non-thermal
  to all electrons 5e-5 (conclusion viii).
- Ions: eps_p ~ 0.1; index softer than p^-4 but poorly constrained (Sect. 4).
- T_e/T_i downstream: ~0.23 (Fig. 7 and Sect. 4) or ~0.3 (Sect. 4, conclusion v). Upstream pre-heating once
  Bell-mediated: T_e ~ 3 T_0, T_p ~ 4 T_0 one ion gyroradius ahead of the shock (conclusion v).
- Minimum energy of electrons escaping upstream ~ gamma_0^2 m_i c^2 (Sect. 4, Figs. 3 and 8).
- E_max equal for both species and proportional to t after omega_pi t ~ 1500 (Fig. 9).
- Box width: eps_e drops when L_y = 25 c/omega_pi (below the reflected-ion gyroradius ~30 c/omega_pi); no
  difference between 200 and 500 c/omega_pi (Sect. 5, Fig. 10).
- Mass ratio: m_i/m_e = 16 gives eps_e more than 10 times larger and no Bell stage by omega_pi t ~ 2000;
  64 and 160 agree, eps_e slightly larger at 64 (App. A).
- Out-of-plane B_0 at 10 deg: eps_e about twice as large early on, the same once turbulence saturates (App. B).
- theta_B = 55 deg: no particles escape upstream, no upstream turbulence, no power law, E_max constant; ion SDA
  bump only; downstream T_p ~ 6 T_e (Sect. 6, Figs. 11-12).

## Assessment
- Shown by simulation: the Weibel-to-Bell transition; filaments and cavities that make patches of the front
  locally sub- and superluminal; suppression of eps_e in narrow boxes; no acceleration at the superluminal angle.
- Not shown: the electron injection path. Electron orbits were not followed; SDA followed by DSA "as in Park et
  al. (2015)" is an inference from the minimum escape energy (Sect. 4, conclusion vi). Whistlers seen in the ion
  foot at late times are only suggested as a pre-heating source (Sect. 3).
- Conjectured: for v << c the filamentary mode may matter less for electrons because even the amplified field
  stays subluminal; this would explain the agreement of eps_e with the 1-D result of Park et al. (Sect. 7).
- Limits stated: 2-D (structures may differ in 3-D; a 3-D run must be wider than the ion gyroradius); one speed,
  one magnetisation; reduced mass ratio; dependence of eps_e on inclination, sigma_0 and speed left open.
- Caution: mildly relativistic, so the 55 deg result is a superluminal case and does not transfer to
  non-relativistic oblique shocks, where 55 deg is subluminal.
- The paper quotes Park et al. (2015) as eps_e ~ 1e-4 (Sect. 1) and "a few 1e-4" (Sect. 7); Park et al. give
  K_ep only, not eps_e. The speed of that run is quoted as 0.1 c in one place and 0.05 c in another.
