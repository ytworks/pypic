# Yao et al. (2022): Detailed characterization of a laboratory magnetized supercritical collisionless shock and of the associated proton energization

- Ref: Matter Radiat. Extremes 7(1), 014402 (2022; online 2021-12-08); DOI 10.1063/5.0055071;
  arXiv:2104.12170 (read); peer-reviewed.
- Read: 2026-11-13, skimmed (abstract, Sect. I, II.1-II.3, III, IV, V of the arXiv version);
  topics: spherical-shock-injection

## Claim
Companion paper of Yao et al. 2021 (note `yao2021laboratory`). A laser-driven piston in a magnetised
hydrogen gas makes a "globally mildly super-critical" quasi-perpendicular shock that leaves the piston.
1-D and 2-D PIC runs show that the protons are energised by shock surfing.

## Method
- JLF/Titan: a laser (1 ns, 70 J, 1.6e13 W/cm2) on a Teflon (CF2) target; H2 gas of 1e18 cm^-3; 20 T
  (Sect. II.1). Interferometry, x-ray spectroscopy, proton radiography, Thomson scattering.
- FLASH 3-D MHD of the piston and the shock (Sect. III).
- PIC (SMILEI, 1D3V; Sect. IV): L_x = 2048 d_e = 11 mm, m_p/m_e = 1836. It does not model the piston: a
  shocked plasma that flows at v_1 = 1500 km/s fills the left half, the ambient plasma the right half.
  The shock width is initialised to one d_i = 200 micron. The run lasts about 2.5 ns.
- 2-D check: L_x = 8 mm, L_y = 0.8 mm, to 0.7 ns, "due to the limitation of the computational resources".

## Numbers to compare against
- Sect. I: "Contrary to Schaeffer et al. 2017a, in our setup, the expanding plasma and the magnetic field
  were decoupled as the higher Z piston evacuates the magnetic field and was thus unmagnetized." The
  result is a highly magnetised ambient plasma and a high-beta piston.
- Sect. I: with 20 T against the 8 T of Schaeffer et al. 2017 "the shock was able to fully separate from
  the piston, which is crucial for its characterization".
- Sect. II.1, Fig. 2: at 4 ns the piston front and the shock front "are well detached from each other".
  No shock forms without the field (ambient gas only).
- Sect. III, FLASH at 2 ns: forward shock near x = 1.4 mm, reverse shock in the piston near 0.8 mm.
  The MHD run has no foot and too low temperatures: "the shock evolution is dominated by kinetic effects".
- Sect. IV: the box is "more than 10 r_Li", with r_Li = v_1/omega_ci about 0.8 mm. The compression is
  n_i,DS/n_i,US about 4.
- Sect. IV, on the curved front: "the PIC simulation represents the tip of the semi-sphere shaped
  expanding shock front at a single slice of z-direction, where the B-field is strictly perpendicular to
  the plasma flow and the shock is the strongest". The radiography integrates over "all other plasma flow
  directions in the xz-plane, which are not perpendicular to the B-field and the corresponding shocks
  are weaker".
- Supercritical threshold used: M_ms >~ 2.7 (Sect. I).

## Assessment
- For this topic: the experiment has a shock on a half sphere in a uniform field, and the authors say
  that its strength varies along the front. They model only the perpendicular tip, in 1-D, at the real
  mass ratio, in a box of about 10 ion gyroradii. The variation along the front is not simulated.
- A stronger ambient field made the shock leave the piston. In our units this is a lower M_A at a fixed
  debris speed; it agrees with Schaeffer et al. 2020 (separation is clearest at lower Mach number).
  This reading is ours.
- The piston is high-beta and unmagnetised, a second way to drive besides the Larmor-coupled piston.
- "Mildly super-critical" is the best that this finite piston gives (M_ms = 3.1 in the first paper).
- Protons only; electrons are not analysed. The PIC starts from a formed shock, so it gives no formation
  length.
