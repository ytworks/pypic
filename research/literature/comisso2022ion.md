# Comisso & Sironi (2022): Ion and electron acceleration in fully kinetic plasma turbulence

- Ref: Astrophys. J. Lett. 936, L27 (2022); DOI 10.3847/2041-8213/ac8422; arXiv:2209.04475;
  peer-reviewed.
- Read: 2026-11-13, skimmed (ar5iv HTML, set-up, figures, text); topics: shock-alfvenic-turbulence

## Claim
In 3-D electron-ion PIC of decaying turbulence at low beta, both ions and electrons form nonthermal
power laws; ions are harder than electrons and both steepen with beta. Electrons are injected by E_par
at reconnecting current sheets and become field-aligned; ions gain mainly from E_perp and stay isotropic.

## Method
- TRISTAN-MP, 3-D periodic, 1400^3 cells, d_e = 3.3 cells, 64 ppc on average.
- m_i/m_e = 50, L/d_i = 60, omega_pe/Omega_ce = 1, v_A = 0.14 c (trans-relativistic electrons,
  nonrelativistic ions).
- beta0 = 0.32, 0.16, 0.08 (fiducial), 0.04; beta_i = beta_e; k_B T0 = 0.08-0.01 m_e c^2 (fiducial 0.02).
- Initial transverse magnetic fluctuations, delta B_rms0/B0 = 1, energy-containing scale l = L/3
  (= 20 d_i). Decaying (no driving).

## Numbers to compare against
- About 2/3 of the turbulent magnetic energy ends in particles, split equally ions/electrons
  (kappa = 1/3 each), by the end of the run; mean gain Delta eps_i,e ~ 0.5 m_e c^2.
- beta0 = 0.04: p_i ~ 3.5, p_e ~ 4.2 (Fig. 2); other beta only in insets.
- Snapshots at t = 1.25 and 8 l/v_A (Figs. 1, 4); pitch-angle statistics over t v_A/l = 2-6 (Fig. 3).
- Injection: electrons where |J|/J_rms ~ 2.5, ions ~ 1; criterion d eps/dt >= 0.2 R_rec v_A m_e c omega_pe,
  R_rec = 0.1, averaged over 150 omega_pe^-1 (Fig. 4).
- Characteristic energy eps* = (3/2 + 2 kappa/beta0) k_B T0.

## Assessment
- Closest published analogue of our "turbulence without shock" null run: e-i PIC, delta B/B0 = 1, but
  beta <= 0.32, v_A = 0.14 c and l = 20 d_i. It shows that such turbulence alone heats electrons to
  ~ eps* and makes a steep electron tail (p_e ~ 4.2) within ~ 8 l/v_A.
- Time scale for us: the magnetic energy is mostly converted within a few l/v_A; a shock run that lasts
  several l/v_A in the upstream frame will see the upstream decay substantially unless driven.
- Limits: m_i/m_e = 50; high v_A/c (electron energies trans-relativistic); no decay-rate value in the
  text; electron-only heating rate not separated from acceleration; low beta only.
