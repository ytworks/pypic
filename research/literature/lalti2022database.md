# Lalti, Khotyaintsev, Dimmock, Johlander, Graham & Olshevsky (2022): A database of MMS bow shock crossings compiled using machine learning

- Ref: JGR Space Physics (2022); DOI 10.1029/2022JA030454; arXiv:2203.04680 (volume and article number not
  recorded). Peer reviewed. Read in the arXiv version.
- Read: 2026-11-14, abstract, the description of the parameters and the statistical results (arXiv PDF
  text; some digits may be lost in the extraction, only values checked in context are quoted); topics:
  qpar-electron-conditions

## Claim
A machine-learning search of the MMS data gives a public list of 2797 crossings of Earth's bow shock with
their upstream parameters; the ion acceleration efficiency is higher at quasi-parallel crossings.

## Method
- Crossings from October 2015 to December 2020, found by a neural network on ion distributions and checked
  by eye.
- Upstream parameters from time-shifted OMNI data (349 events lack at least one); M_A in the normal
  incidence frame; the fast Mach number uses a typical T_e of 12.06 eV (Newbury et al. 1996).

## Numbers to compare against
- 45.5 % quasi-parallel, 51.4 % quasi-perpendicular, 3.1 % undetermined.
- Most crossings at M_A 5-15, as in the solar wind statistics; some values up to ~150 at very low field
  with large uncertainty.
- The compression approaches 4 at large M_A.
- Ion acceleration efficiency 14 +- 11 % (quasi-parallel) against 8 +- 8 % (quasi-perpendicular); no clear
  dependence on M_A.
- Example quasi-parallel crossing: theta_Bn 33 deg, M_A 9.4.

## Assessment
- Ions only; no electron acceleration, no connection times.
- For this topic: the source of the M_A range of Earth's bow shock in `inputs.md`. The value "~8" of the
  case table lies inside 5-15; raptis2026relativistic writes "M_A ~ 4" for Earth.
