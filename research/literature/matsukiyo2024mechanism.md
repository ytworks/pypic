# Matsukiyo et al. (2024): Mechanism of generating collisionless shock in magnetized gas plasma driven by laser-ablated target plasma

- Ref: Phys. Plasmas 31(11), 112902 (2024); DOI 10.1063/5.0230232; arXiv:2411.09083 (read); peer-reviewed.
- Read: 2026-11-13, full (arXiv version); topics: spherical-shock-injection

## Claim
In 1-D full PIC, a target plasma injected for a finite time into a magnetised gas plasma makes a
"seed-shock". A longer injection keeps the magnetic piston strong and gives a faster shock with a higher
Mach number. A short injection disperses and the piston weakens fast.

## Method
- 1-D full PIC; system about 130 c/omega_Gi long (G: gas ions). m_G/m_e = 100, with checks at 50 and 200.
- Target plasma: charge 6, m_T/m_G = 2, density 7 times the gas (one run with 14), field 6 times the
  ambient field. It is injected at v_in = 21.7 v_A (4.4 T) or 15.9 v_A (6 T) for a time T_in with
  Omega_Gi T_in = 0.019 to 0.21.

## Numbers to compare against
- Sect. IV.1, Fig. 5: "As T_in increases, both v_sh and M_A of the seed-shock increase significantly.
  This is attributed to the suppression of target plasma deceleration by injecting target plasma for an
  extended duration, thereby maintaining a strong magnetic piston effect."
- "in all runs, M_A > 3". The values of Fig. 5 are not in the text.
- The formation mechanism changes near Omega_Gi T_in = 4-5e-2. The seed-shock needs a precursor of
  gyrating gas ions.
- The structure varies with a period of about a quarter of the gas-ion gyro-period.
- A denser target or a stronger target field gives a stronger piston and a faster seed-shock.
- A lower mass ratio lets more target ions enter the gas.
- The authors call the 1-D assumption "quite strong" (discussion).

## Assessment
- For this topic: the clearest statement that the Mach number of the shock follows the duration and the
  density of the piston, not only its first speed. This supports a heavier piston at a fixed speed.
- The injection speed is 16-22 v_A, far above our M_A = 5; the shocks are above M_A = 3, not at the
  injection speed.
- The remark on the mass ratio matters at 16: more debris ions may cross into the ambient plasma than
  at 100. The paper gives no number for this.
- 1-D and planar: no geometric weakening, no curvature. Electron acceleration is not the subject.
