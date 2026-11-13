# Bresci, Lemoine & Gremillet (2023): Particle acceleration at magnetized, relativistic, turbulent shock fronts

- Ref: Phys. Rev. Research 5, 023194 (2023) (Crossref); DOI 10.1103/PhysRevResearch.5.023194; arXiv:2303.11394;
  peer-reviewed
- Read: 2026-11-13, full (arXiv v2; appendices skimmed); topics: shock-alfvenic-turbulence

## Claim
First PIC simulations of a relativistic shock propagating into well-developed magnetized turbulence: strong
upstream turbulence (dB > B0) revives acceleration in a superluminal configuration that forbids it in a laminar
flow, with power laws s about 2.5-3.5; the mechanism (diffusive-type vs stochastic shock drift) depends on r_g/l_c;
at larger turbulent magnetization stochastic acceleration in the pre-shock turbulence takes over.

## Method
- PIC (CALDER), 2D3V, pair plasma; inflow gamma = 2 (v = -0.87c), reflecting wall; guide field B0 along z (out of
  plane); turbulence driven by a Langevin antenna (in-plane dB_x, dB_y) in a band 4200-4800 c/omega_p, then decays
  before reaching the shock.
- l_c about 200 c/omega_p; L_y = 600 c/omega_p; dx = 0.1 c/omega_p; box 48000 x 6000 cells; dt = 0.99 dx/c.
- Runs {sigma0, sigma_dB, kT/mc^2}: S1 {2e-4, 2e-2, 0.1}, S2 {6e-4, 1e-2, 4}, S3 {6e-4, 0.1, 4}; controls S2a (no
  turbulence), S2b (no shock), S4.

## Numbers to compare against
- Shock speed about 0.4c in the downstream frame, gamma_sh about 2.7; compression about 3.5.
- Tail indices s about 2.5 (S1, diffusive-type) and 3.5 (S2, stochastic shock drift); r_g/l_c about 0.1 (S1) and
  about 1 (S2) at the onset of the tail (Figs. 2, 4).
- No tail in S2a (no turbulence) nor S2b (no shock) - the paper's two null tests (Fig. 2).
- Mechanism separated per particle by the work of E0 (mean motional field) vs of dE_z (Fig. 4).

## Assessment
- Relativistic pair plasma, no electron-ion physics; turbulence fluctuations in the plane with B0 out of plane.
- Methodologically closest to what our topic needs: controls without turbulence and without shock, and a
  per-particle work decomposition that separates shock drift from scattering-based acceleration. A meandering test
  could reuse that design (work by mean motional field during repeated front encounters).
- 2D3V only, chosen for stated cost constraints (large l_c, long interaction time, rest-frame driving).
