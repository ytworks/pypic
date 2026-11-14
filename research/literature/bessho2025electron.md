# Bessho, Chen, Hesse, Ng, Wilson & Stawarz (2025): Electron Acceleration in Thinning Nonreconnecting Current Sheets in a Quasi-parallel Shock

- Ref: ApJ 992, 183; DOI 10.3847/1538-4357/ae0327. Peer reviewed.
- Read: 2026-11-14, abstract; topics: qpar-electron-conditions

## Claim
In the turbulence of a quasi-parallel shock at Earth's bow shock parameters, thinning current sheets that do
not reconnect energize electrons through the motional electric field; the electrons "can reach the injection
energy for diffusive shock acceleration".

## Method
- 2-D PIC of a quasi-parallel shock "Using parameters in the Earth's bow shock"; electron trajectory
  tracing; an analytic rate for unmagnetized electrons.

## Numbers to compare against
- The number of electrons energized in non-reconnecting current sheets "is about a quarter of that in
  reconnection regions".

## Assessment
- Mass ratio, shock speed and duration are not in the abstract; `bessho2023electron` (same group and
  set-up family) ran ~20 ω_ci⁻¹.
- For this topic: with `bessho2023electron` the closest published fully kinetic work at one of our cases
  (Earth). It answers "by which process do electrons gain energy in the transition region" and stops at
  the injection energy; it does not follow acceleration beyond it and does not confront an event. The
  current sheets live in the plane of a 2-D run; a quasi-1-D run cannot hold them. That is a named limit of
  our run E and one more reason for the 2-D run X.
