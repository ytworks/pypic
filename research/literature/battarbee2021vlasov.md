# Battarbee, Brito, Alho, Pfau-Kempf, Grandin, Ganse, Papadakis, Johlander et al. (2021): Vlasov simulation of electrons in the context of hybrid global models: an eVlasiator approach

- Ref: Ann. Geophys. 39, 85-103; DOI 10.5194/angeo-39-85-2021. Peer reviewed. 11 authors.
- Read: 2026-11-14, abstract; topics: qpar-electron-conditions

## Claim
Electron distributions can be evolved with a Vlasov solver inside snapshots of a global hybrid-Vlasov
simulation; the scheme is stable and gives Langmuir waves and Bernstein modes.

## Method
- eVlasiator: Vlasov electrons in the fixed ion-scale fields of a Vlasiator run. Test case in the near-Earth
  magnetotail plasma sheet.

## Numbers to compare against
- None in the abstract.

## Assessment
- The electrons do not act back on the ion-scale fields; the test case is not a shock.
- For this topic: the noise-free relative of "test electrons in saved fields" (`gupta2025speed` (S),
  `hanusch2020electron`). It shows the class of method "kinetic electrons in snapshots of an ion-kinetic
  global run" exists for Earth's geometry. No electron acceleration at the bow shock is reported in the
  abstract.
