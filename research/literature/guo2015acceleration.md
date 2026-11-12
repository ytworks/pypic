# Guo & Giacalone (2015): The acceleration of electrons at collisionless shocks moving through a turbulent magnetic field

- Ref: Astrophys. J. 802, 97 (2015) (volume and article from the DOI; Trotta & Burgess 2019 print
  "Astrophysical Journal, 802"); DOI 10.1088/0004-637X/802/2/97; arXiv:1409.5854; peer-reviewed
  (arXiv comment: accepted by ApJ).
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
When the upstream magnetic fluctuations are much stronger than the mean field (variance about 10 B_0^2),
electrons are accelerated about equally well for every angle between the mean field and the shock normal,
because the local field at the front is often nearly perpendicular and the upstream field mirrors
electrons back to the front.

## Method
- Test-particle electrons (relativistic equations, Bulirsch-Stoer) in prescribed, kinematic fields; no
  self-consistent shock. Planar shock at x = 0 with a tanh flow profile, U_1 = 500 km/s, compression 4,
  width 0.2 U_1/Omega_i; the front stays planar by assumption.
- Field: mean B_0 in the x-z plane plus isotropic 3-D fluctuations (random wave modes, Kolmogorov-like),
  correlation length 50 U_1/Omega_i, wavelengths 5 to 500 U_1/Omega_i, frozen into the flow; transverse
  components compressed at the shock; E = -U x B/c. No overshoot, no cross-shock potential.
- Mach numbers: only "much larger than 1" (assumed). beta, m_i/m_e, U_1/v_A: not stated.
- Electrons: 1 keV, isotropic, injected at x = -5 U_1/Omega_i at random times; initial gyroradius about
  0.02 U_1/Omega_i; about 5e6 per run; followed to 400 Omega_i^-1; no spatial boundary.
- 21 runs: variance 0.1, 1, 10 B_0^2 and mean theta_Bn = 0, 15, ..., 90 deg (Table 1).

## Numbers to compare against
- Fraction of electrons above 10 keV at the end (Table 1, Fig. 2):
  variance 0.1: 0.0 % for 0-45 deg, 0.2 % (60), 0.97 % (75), 9.3 % (90);
  variance 1: 1.8 % (0), 2.1 % (15), 2.58 % (30), 3.78 % (45), 6.04 % (60), 9.5 % (75), 13.5 % (90);
  variance 10: 7.3 % (0), 7.5 % (15), 7.9 % (30), 8.6 % (45), 9.7 % (60), 10.9 % (75), 11.4 % (90).
- At 90 deg the fraction is higher for variance 1 (13.5 %) than for 10 (11.4 %): stronger turbulence
  lowers the local angle there (Sect. 3).
- Maximum energy: about 1 MeV for variance 10 at all angles; for variance 0.1 only 30 keV at 90 deg and
  7 keV at 0 deg (Sect. 3, Fig. 5).
- Scales: smallest wavelength = 250 initial gyroradii = 25 shock widths (from Sect. 2 values). Upstream
  mirror distance near the shock about 20 U_1/Omega_i (Sect. 3).
- One electron at mean angle 15 deg, variance 10: energy x 7 in one encounter while |B| rises (drift
  signature), plus many small gains; more than 100 E_0 in total (Fig. 3, 4).

## Assessment
- Shown (test particles in a model field): the dependence on the mean angle weakens as the variance
  grows; accelerated electrons pile up at the front (Fig. 1).
- Not self-consistent: front shape, turbulence level and spectrum are inputs. The origin of fluctuations
  this strong (variance 10 B_0^2) is left open (Sect. 4). No fluctuations at the electron gyroradius, so
  no pitch-angle scattering at injection; returns are by mirroring (Sect. 3).
- Stated, not shown: this injects electrons into diffusive shock acceleration; the link to supernova
  remnants rests on the inferred field amplification (Sect. 4). The spectral index is not given in the text.
- Contrast it draws: heliospheric shocks accelerate electrons mostly where quasi-perpendicular; remnants
  seem to do so at any angle (Reynolds 2008); Saturn's quasi-parallel bow shock gave MeV electrons with
  dB/B_0 ~ 10 (Masters et al. 2013). Earlier runs with variance <= B_0^2 (Guo & Giacalone 2010) and PIC
  runs (Kato & Takabe 2010; Riquelme & Spitkovsky 2011; Niemiec et al. 2012) found no such acceleration.
- It argues that the upstream-reflection gain seen in 1-D PIC (Kato 2014; Park et al. 2014; U_1 =
  0.1-0.3 c, gain ~ 2 U_1 E/c per reflection) is much smaller at real, slower shocks (Sect. 3).
- Doubts: captions of Fig. 3 and 4 are identical and give times 120-230, the text 300-380 and 300-325;
  line styles in the Fig. 5 caption disagree with the text.
