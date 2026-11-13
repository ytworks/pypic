# Comisso & Sironi (2018): Particle acceleration in relativistic plasma turbulence

- Ref: Phys. Rev. Lett. 121, 255101 (2018); DOI 10.1103/PhysRevLett.121.255101; arXiv:1809.01168;
  peer-reviewed.
- Read: 2026-11-13, skimmed (ar5iv HTML, set-up, figures, Eqs. 1-3); topics: shock-alfvenic-turbulence

## Claim
Decaying turbulence in magnetically dominated pair plasma makes a power law that hardens with sigma
and with delta B/B0. Particles are first injected at reconnecting current sheets by E_par and then
gain energy stochastically from the turbulent fluctuations.

## Method
- TRISTAN-MP, 2-D (to L/d_e0 = 6560, 10 cells per d_e0, 16 ppc) and 3-D (to L/d_e0 = 820, 3 cells
  per d_e0, 4 ppc; l = L/4). Convergence checked to 256 (2-D) and 16 (3-D) ppc.
- Initial state: only delta B_x, delta B_y, modes m, n = 1..N (N ~ 8), equal power, random phases;
  in 3-D modulated in z with wavelengths L and L/2. No initial velocity; not in pressure balance, so it
  relaxes quickly to turbulence. l = 2 pi / k_N = L/N.
- theta0 = kT/mc^2 = 0.3; sigma0 = 2.5-160; delta B_rms0/B0 = 0.5-4.

## Numbers to compare against
- sigma0 = 10, delta B/B0 = 1: p ~ 2.9 (Fig. 2); hardest cases p ~ 2.
- ~70 % of the initial turbulent energy into particles by the end (Fig. 1 inset gives the decay; no
  decay time is quoted in the text).
- Cutoff gamma_c grows linearly with l (not L); slope size-independent (Fig. 2); same in 3-D (Fig. 3).
- Injection: 97 % of tail particles (gamma > 30 at ct/l = 12) show a sudden gain; 95 % of injected
  sit where |J_z| > 2 J_rms, vs 9 % of all particles (Fig. 4). E_par ~ beta_R delta B_rms, beta_R ~ 0.05
  (Eq. 1); stochastic stage fitted by a Fokker-Planck model (Eqs. 2-3).

## Assessment
- Gives the standard decaying-turbulence initial condition for PIC (random-phase delta B only), and
  its known drawback: a fast non-equilibrium transient. For a shock upstream this transient would be
  confused with shock physics, so an upstream must be pre-relaxed or start from an MHD-consistent state.
- Injection at current sheets with E_par is the turbulence-only injection channel our null test must
  measure, but sigma >> 1 here; the efficiency at sigma << 1 is not given.
