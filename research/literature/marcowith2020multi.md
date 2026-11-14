# Marcowith, Ferrand, Grech, Meliani, Plotnikov & Walder (2020): Multi-scale simulations of particle acceleration in astrophysical systems

- Ref: Living Rev. Comput. Astrophys. 6, 1; DOI 10.1007/s41115-020-0007-6; arXiv:2002.09411. Peer-reviewed
  review. Read in the arXiv version.
- Read: 2026-11-14, arXiv PDF text: the passages on PIC noise, dimension, mass ratio, non-relativistic
  shocks (table 2), hybrid and PIC-MHD methods, the open questions and the conclusion; topics:
  qpar-electron-conditions

## Claim
A review of the numerical methods for particle acceleration (PIC, Vlasov, hybrid, PIC-MHD, kinetic-fluid)
with their ranges of validity. Its stated role of PIC: "an efficient tool to probe the efficiency of
particle acceleration for a given parameter set", limited in duration and box size; "The goal and common
approach is to provide robust scalings which can be included in fluid simulations as sub-grid
prescriptions."

## Method
Review. Table 2 lists the non-relativistic PIC shock runs known then, among them `park2015simultaneous` (S)
(1-D, mass ratio 100, 30°, M_A 20, 0.1 c) and `kato2015particle` (S) (1-D, mass ratio 30, 30°, 0.37 c).

## Numbers to compare against
- Duration: "a couple of ion gyro-periods for full-PIC and 10^3 ion gyro-periods for hybrid-PIC"; boxes
  "typically less than thousands of ion skin depths". Elsewhere: typical full-PIC runs "not longer than a
  10 ω_ci⁻¹". Our runs: 50-300 ω_ci⁻¹.
- Hybrid is about two orders of magnitude cheaper than full PIC.
- Whistler criterion M_A / (m_i/m_e)^(1/2) < 1 (Matsukiyo & Scholer 2003).
- Non-thermal electron-to-proton ratio at high-M_A quasi-parallel shocks ~10⁻³, "roughly proportional to
  u_sh/c" (from Park 2015 and Kato 2015).
- Energy gain per Fermi cycle ΔE/E ≃ u_sh/c: the reason the full tail is hard to capture at non-relativistic
  speed.

## Assessment
- Noise: "due to the introduction of a finite number of macro-particles, PIC simulation suffers from the
  highly exaggerated level of noise. This well known short-coming of the PIC method makes it less adapted
  to treating problems for which regions of phase-space where the distribution function assumes small
  values (e.g., in its tail) can impact the physics at play". Vlasov codes are "virtually noise-free" but
  costly.
- Dimension: Weibel filamentation is suppressed in 1-D; "This restriction is less severe in
  non-relativistic case for quasi-parallel shocks, where non-resonant streaming or Bell instability can be
  triggered in 1D and sustains particle scattering (or mirroring) on both sides of the shock front.
  Additional step from 2D to 3D is important to correctly deal with particle scattering properties as the
  topology of the turbulent magnetic field is different."
- Mass ratio: the conditions for shock surfing and the Buneman instability "depend on the ion to electron
  mass ratio adopted in the PIC simulations"; "a small ion to electron mass ratio suppresses the
  propagation of oblique whistler waves (Scholer and Matsukiyo, 2004), whereas these waves can become
  over-dominant to heat/energize electrons in the foot".
- Diffusion: the Bohm limit "has been widely favored in the literature … may be in part due to the fact
  that it is the most favorable case, and also may just stem from habit given the lack of clear theoretical
  alternatives"; the diffusion coefficient, "although it should in principle be computed from f itself, is
  commonly prescribed as a function of p".
- Box width: a statement that a too narrow box suppresses electron acceleration because it suppresses Bell
  modes is made there for mildly relativistic shocks (Crumley et al. 2019), not for our regime.
- PIC-MHD: "promising when one is mainly interested in dynamics of supra-thermal particles … Yet, this
  requires to robustly prescribe how some part of thermal particles … are promoted to non-thermal status".
- Open question 2 of the review: "which configuration, quasi-parallel or quasi-perpendicular, in
  non-relativistic shocks is more efficient for ion/electron acceleration?".
- "The most precise methods, PIC simulations, are restricted to tiny scales and narrow dynamical ranges."
- For this topic: the review places what a PIC shock run can state: an efficiency or a scaling at a given
  parameter set, to be carried further by another model. That is the form of our plan (measure the
  injection rate and D̃_e, then a transport model). It also names the three limits our results must carry:
  noise in the tail, the 1-D to 2-D to 3-D steps, and the mass ratio. The ratio "roughly proportional to
  u_sh/c" is a published dependence on speed that bears on null test N4.
