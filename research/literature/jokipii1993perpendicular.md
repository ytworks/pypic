# Jokipii, Kota & Giacalone (1993): Perpendicular transport in 1- and 2-dimensional shock simulations

- Ref: Geophys. Res. Lett. 20(17), 1759-1761 (1993) (Crossref); DOI 10.1029/93GL01973; peer-reviewed
- Read: 2026-11-13, abstract; topics: shock-alfvenic-turbulence

## Claim
In 1D and 2D simulations in which all quantities are independent of a coordinate not parallel to the magnetic
field, charged particles are tied to the convected field lines (to within a gyroradius) because the canonical
momentum along the ignorable coordinate is conserved. Perpendicular diffusion is therefore suppressed, and
acceleration at quasi-perpendicular shocks cannot be studied properly in such simulations.

## Method
- Analytic argument (conservation of canonical momentum along the ignorable coordinate).

## Numbers to compare against
- None (analytic).

## Assessment
- Design constraint for the topic (my reading of the argument): with ignorable coordinate y, p_y + q A_y / c is
  conserved, so a particle stays within about a gyroradius of one A_y contour, i.e. of the in-plane projection of
  its field line. Meandering of that line within the plane is followed (electrons can return to a rippled or
  corrugated front along it), but true cross-field transport between lines is absent. A 2D run can test
  "meandering along in-plane field lines"; only 3D can test transport that needs field-line separation.
- Guo & Giacalone (2010) state this limit for their 2D hybrid runs (electrons stop gaining once their line is
  downstream). A fully kinetic test of meandering in 2D inherits it; 3D removes it.
