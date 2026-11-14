# Nyberg, Vuorinen, Afanasiev, Annie John & Vainio (2026): Electron Acceleration Simulations (EASI): An open-source tool to simulate electron acceleration in shocks

- Ref: A&A 708, A67; DOI 10.1051/0004-6361/202557951. Peer reviewed.
- Read: 2026-11-14, abstract (the journal page could not be read); topics: qpar-electron-conditions

## Claim
An open-source 1-D Monte Carlo model of electron acceleration at shocks, with a shock of finite thickness and
an adjustable profile of the mean free path, accelerates electrons from 1 keV to above 100 keV at coronal
shock parameters.

## Method
- Monte Carlo solution of a transport equation with magnetic focusing and pitch-angle scattering.
- Injection mono-energetic or Maxwellian; the mean free path is an input profile.

## Numbers to compare against
- ">100 keV from a 1 keV monoenergetic injection" at coronal shock parameters; electron beams at high
  obliquity.

## Assessment
- The scattering (mean free path in space and energy) and the injected population are inputs.
- For this topic: the transport-model side of our chain exists as an open tool, for coronal shocks. Its
  inputs are what our runs are to measure: D̃_e(x, p̃) is a mean free path profile, the injection rate a
  source. Two readings for the plan: our transport model can be checked against this code; or this code can
  serve as the transport model, if its equation covers our cases. Not decided here; the equation was not
  read.
