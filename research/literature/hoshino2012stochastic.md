# Hoshino (2012): Stochastic Particle Acceleration in Multiple Magnetic Islands during Reconnection

- Ref: Phys. Rev. Lett. 108, 135003 (2012); doi:10.1103/PhysRevLett.108.135003; arXiv:1201.0837
- Read: 2026-11-12, full (arXiv v1); topics: (none yet)

## Claim
When a current layer breaks up into many merging islands, the most energetic particles stay outside
the islands in the strong-field region. They gain energy each time they meet the Alfvénic outflow jet
of an island-merging reconnection site. This makes Fermi's stochastic scheme first order: the gain
per encounter scales as VA/c, not (Vc/c)^2.

## Method
- 2D relativistic PIC with Maxwell's equations solved in Fourier space (to limit numerical
  dispersion); doubly periodic, 512 x 512 cells, 200 particles per cell.
- Electron-positron plasma; four Harris sheets (peak n0) on a background nb = 0.05 n0;
  T0/mc^2 = 5, drift +-0.5c, (Omega_c/Omega_p)^2 = 20, sheet thickness lambda = 5 cells
  (box about 100 lambda); uniform guide field Bz = B0. Run to t = 320 lambda/c.
- Diagnostics: density and field lines (Fig. 1a-c); energy spectrum with a Maxwellian fit
  (Fig. 1d); density maps of momentum-selected particles, p/mc in (50, 200) and (200, inf)
  (Fig. 1e-f); one representative particle's energy against the local in-plane |B| and its orbit
  over field maps (Fig. 3).
- Theory: a cartoon (Fig. 2c) and a de Hoffmann-Teller frame argument; no quantitative model.

## Numbers to compare against
Values as stated in the text; plots not digitised.
- Maxwellian fit to the final spectrum T/mc^2 = 7.3; tail above gamma ~ 40; maximum gamma ~ 10^3,
  gyroradius about half the box (Fig. 1d, t = 320 lambda/c).
- Nonthermal particles carry about 66 % of the total particle energy (text with Fig. 1d; the
  separation method is not given).
- p/mc = 200 has a gyroradius of about 10 lambda, the island size in the merging stage (Fig. 1e-f).
- In-plane |B| <= 0.4 B0 is roughly inside islands, >= 0.6 B0 outside; the traced particle gains
  energy at |B| ~ 0.8 B0 between t = 180 and 200 lambda/c (Fig. 3).
- gamma_max = eBL/(2mc^2) ~ 10^3 (L / 100 lambda); a 4x larger box gave about 4x the maximum energy
  (text, p. 8; run not shown).
- Gain per jet reflection delta_eps/eps proportional to VA/c, against (Vc/c)^2 for classical Fermi
  (text with Fig. 2c).

## Assessment
- Shown by simulation: during merging, high-energy particles avoid island interiors (Fig. 1e-f).
  One traced particle is first accelerated near an X-line (Speiser-type), then gains energy at a
  merging-island jet in the strong-field region (Fig. 3).
- Argued but not demonstrated: the first-order VA/c scaling (no measured gain per encounter, no
  statistics over many particles). Losses in head-tail encounters are taken as negligible from an
  assumed ordering Vc < Vs < VA for beta < 1. The island random velocity is called an open issue.
- Stated but not shown: the same behaviour in nonrelativistic runs, the 4x box test, and that
  escape will produce a power law. The index is not predicted; the paper says it may depend on the
  island dynamics and the losses.
- Limits: 2D, pair plasma (ion-electron left for future work), small box (~100 lambda), a single
  guide-field value and parameter set. The injection problem (how particles leave the islands) is
  left open.
- The evidence for the mechanism rests on one representative trajectory.
