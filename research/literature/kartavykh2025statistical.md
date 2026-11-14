# Kartavykh et al. (2025): A statistical study of energetic particle events associated with interplanetary shocks observed by Solar Orbiter in solar cycle 25

- Ref: A&A 699, A24; DOI 10.1051/0004-6361/202553775 (18 authors: Kartavykh, Rodriguez-Garcia, Heber,
  Wimmer-Schweingruber, Trotta, Gieseler, Droge, Kollhoff, ...). Peer reviewed; open access.
- Read: 2026-11-14, abstract, section 2.1 (electrons), conclusions, table 1, table A.1 header (text
  extracted from the journal PDF; some symbols lost); topics: qpar-electron-conditions

## Claim
Electron acceleration at interplanetary shocks between ~0.3 and 1 au is rare and belongs to oblique or
quasi-perpendicular shocks; the mechanism "is unlikely to be DSA", shock drift acceleration remains a
candidate. Protons near 100 keV respond at most shocks, mainly parallel and quasi-parallel ones.

## Method
- Solar Orbiter EPD (lowest electron channel used 35.6-63.5 keV; residuals after a correction for proton
  contamination) at the forward shocks of the mission's list; shock parameters from in situ data.

## Numbers to compare against
- 58 forward shocks; 51 used for the electron search.
- Electron response at <= 100 keV in 4 of 51 (~8 %); theta_Bn of the four: 88, 42, 80, 61 deg.
- Above 100 keV in one case (~2 %): 2022-07-25, theta_Bn 61 deg, M_A 8.4, 820 km/s, electrons to ~250 keV.
- About 83 % of the shocks show a response in ~100 keV protons.
- Earlier 1-au values quoted there: 17 % at 45 keV (Lario et al. 2003), 1 % (dresing2016efficiency).
- Table 1, averages +- standard errors of the 51 shocks (mean distance 0.79 au): theta_Bn 54 +- 3 deg,
  shock speed 512 +- 24 km/s (spacecraft frame), M_A 3.3 +- 0.3, M_fast 2.8 +- 0.3, gas compression
  2.0 +- 0.1, magnetic compression 1.9 +- 0.1, upstream dB/B 0.25 +- 0.02, density 19.6 +- 3.4 cm^-3,
  beta 0.84 +- 0.12.
- The same for the 4 shocks with an electron response: theta_Bn 68 +- 9 deg, 545 +- 81 km/s, M_A 6.1 +- 1.0
  (7.4, 3.2, 5.3, 8.4), M_fast 4.6 +- 0.7, gas compression 1.9 +- 0.3, beta 1.05 +- 0.34. The text repeats
  "6.1 +- 1.0 compared to 3.3 +- 0.3".
- Each table entry carries a second value; for the four M_A it is not their median in our check, so it is
  not used. The authors state uncertainties up to 20 % for theta_Bn or the speed of a single shock.

## Assessment
- The four electron cases have theta_Bn 42-88 deg; the authors call them oblique or quasi-perpendicular
  (42 deg lies just below the usual 45 deg limit). No case has theta_Bn below 40 deg, while protons respond
  mainly at parallel and quasi-parallel shocks.
- The only > 100 keV case has M_A 8.4, close to the M_A of the PSP shock (7-7.5, jebaraj2026wave), but
  theta_Bn 61 deg against ~8 deg; Mach number alone does not separate the two.
- Proton contamination removes candidates; the 8 % is a lower bound for that reason.
