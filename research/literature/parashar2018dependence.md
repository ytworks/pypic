# Parashar, Matthaeus & Shay (2018): Dependence of kinetic plasma turbulence on plasma beta

- Ref: Astrophys. J. Lett. 864(1), L21 (2018) (Crossref); DOI 10.3847/2041-8213/aadb8b;
  arXiv:1807.11371; peer-reviewed.
- Read: 2026-11-13, skimmed (ar5iv HTML, set-up and figure text); topics: shock-alfvenic-turbulence

## Claim
In fully kinetic PIC of decaying turbulence, beta changes spectra, intermittency, decay rates and the
partition of energy exchange: higher beta gives more total heating and a larger proton share.

## Method
- P3D, 2.5-D (x-y plane, three components), L = 149.6 d_i, 4096^2 cells, 3200 ppc per species
  (~1.07e11 particles), m_i/m_e = 25.
- beta = 0.3, 0.6, 1.2 (beta_i = beta_e); initial delta b/B0 = 0.5 and delta v/v_A = 0.5 (rms),
  modes 2 <= |k| <= 4 (in units of 2 pi/L), random phases. Run > 300 Omega_ci^-1. v_A/c not read.

## Numbers to compare against
- Early relaxation (Omega_ci t < 40) nearly beta-independent; called "a few nonlinear times".
- Fluctuations decay more at larger beta (Fig. 1); proton heating rises with beta and electron heating
  change becomes smaller, so Q_p/Q_e grows with beta (Fig. 2). Values are only in figures.

## Assessment
- Parameter regime closest to our upstream (nonrelativistic e-i, beta ~ 0.3-1, delta b/B0 = 0.5,
  energy-containing scale ~ L/3 to L/2 ~ 40-75 d_i). Time scale: decay over a few tens to a few hundred
  Omega_ci^-1; the null run must cover the same time as the shock run.
- Initial state includes velocity fluctuations (Alfvenic delta v = delta b), unlike the magnetic-only
  start of comisso2018particle; a closer analogue of an Alfvenic upstream.
- Limits: 2.5-D, m_i/m_e = 25; heating fractions only in figures (not quoted here); electron
  acceleration (nonthermal tails) not addressed.
