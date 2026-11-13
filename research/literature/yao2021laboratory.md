# Yao et al. (2021): Laboratory evidence for proton energization by collisionless shock surfing

- Ref: Nature Physics 17(10), 1177-1182 (2021); DOI 10.1038/s41567-021-01325-w; arXiv:2011.00135 (read);
  peer-reviewed.
- Read: 2026-11-13, skimmed (abstract, main text, Table 1, PIC part of the Methods);
  topics: spherical-shock-injection

## Claim
A laser-driven piston in a hydrogen gas with a 20 T field drives a supercritical quasi-perpendicular
collisionless shock. Protons are energised by shock surfing while the shock is fast; the shock then slows
down and becomes subcritical.

## Method
- JLF/Titan and LULI2000; a laser-ablated piston expands into ambient hydrogen (n_e about 1e18 cm^-3)
  in a homogeneous 20 T field. Interferometry, Thomson scattering, proton spectra.
- 3-D MHD (FLASH) for the piston and the shock formation.
- PIC (SMILEI, 1D3V): L_x = 2048 d_e = 11 mm, dx = 0.2 d_e, m_p/m_e = 1836 (tests at 400-1836 and in
  2-D). It does not model the piston: a shocked flow at 1500 or 500 km/s is put next to the ambient
  plasma, with a shock width of one d_i.

## Numbers to compare against
- Supercritical threshold used: M_ms about 2.7. Table 1: the experiment has M_ms = 3.1; the extended
  table gives M_A = 3.4.
- Both the piston front and the shock front are seen at 4 ns.
- The shock speed is about 1500 km/s early; "after 4-5 ns, the shock velocity decreases rapidly to about
  500 km/s, thus becoming sub-critical". Proton acceleration takes place in the first 2-3 ns.
- PIC, high-speed case: forward and reverse shocks at 0.5 ns, reformation at 1.5 ns. The proton gyro-period
  is 3 ns and the gyroradius 0.8 mm at 1500 km/s (Methods).
- PIC, low-speed case: M_ms <~ 1.0 < 2.7; few reflected protons.
- The FLASH run gives a steady shock speed and not the fast decay.

## Assessment
- For this topic: an experiment in which the shock of a finite piston is strong only for a short time
  (about one gyro-period, derived here from 2-3 ns against 3 ns) and then decays below the critical
  Mach number. Our weak shocks may be the same decay seen late.
- The PIC starts from a shock that has already left the piston. It says nothing about formation or about
  the room needed.
- Protons only. No electron acceleration.
