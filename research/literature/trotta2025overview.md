# Trotta, Dimmock, Hietala, Blanco-Cano, Horbury, Vainio, et al. (2025): An Overview of Solar Orbiter Observations of Interplanetary Shocks in Solar Cycle 25

- Ref: Astrophys. J. Suppl. Ser. 277, 2 (2025); DOI 10.3847/1538-4365/ada4a7; arXiv:2410.24007; peer-reviewed.
- Read: 2026-11-13, skimmed (arXiv HTML text up to ~100k characters, via extraction; last ~9k characters not
  read); topics: shock-alfvenic-turbulence

## Claim
The first 100 IP shocks seen by Solar Orbiter (2020 – end 2023, 0.3–1 au) are mostly weak (M_A, M_fms ≲ 3,
compression < 4) and favour quasi-perpendicular geometry. Compared with 1 au surveys, inner-heliosphere shocks
are slightly faster, have higher M_fms and are more often quasi-parallel. About half show shock-induced upstream
fluctuations. IP shocks rarely accelerate electrons efficiently in situ.

## Method
- Detection: TRUFLS algorithm (jumps B ≥ 1.2, n ≥ 1.2, ΔV ≥ 20 km/s), then visual check.
- Normals: MX3 (mixed), magnetic coplanarity when plasma data missing; speed from mass flux.
- SerPyShock: up/downstream averaging windows varied 30 s – 8 min for uncertainties.
- Upstream δB/B0: increments with 1 min lag, averaged over 8 min upstream (catalogue column "Mean upstream δB/B0").
- Waves: wavelet power/ellipticity, visual classification (foreshock 0.01–1 Hz; precursors ~0.5–10 Hz).
- Energetic particles: EPD/STEP, EPT, HET; responses classified spike / plateau / irregular / none.

## Numbers to compare against
- Mach numbers mostly ≤ 3; gas and magnetic compression < 4 (Sect. IV.1, Fig. 2). Medians not given in text.
- Comparison (Fig. 4): 30 fast-forward shocks below 0.8 au vs 586 fast-forward shocks at 1 au (STEREO-A + Wind,
  1995–2009; Kilpua et al. 2015).
- Wave foreshocks in 46 % of shocks; low-frequency precursors in ~62 % of 72 usable shocks; overall ~50 %
  (Sect. IV.2); foreshocks favour low θBn and high Mach number.
- Ions: at ~15 MeV, 88.5 % no response; at ~1 MeV, 72 % no protons (Sect. IV.3).
- Upstream β and δB/B0 are catalogued (Table 1 columns) but no summary values are given in the text read.

## Assessment
- Shown: the parameter space of real heliospheric shocks is low-Mach (≲ 3) and mostly quasi-perpendicular;
  any PIC set-up meant to represent "typical" IP shocks sits there, not at M_A ~ 10.
- Upstream fluctuation detection is visual and refers to shock-induced waves, not to the pre-existing solar
  wind turbulence level; the catalogue's δB/B0 column (1 min lag) would give the latter but was not inspected.
- Internal inconsistencies noted by extraction: ~70 % with 50 keV increases (conclusions) vs 36 % without
  enhancement (Sect. IV.3); 10 % vs 9 % irregular at 15 MeV.
- Electron statistics are deferred to the SERPENTINE catalogue; no electron–θBn or electron–M relation given.
- Open for us: the catalogue (SERPENTINE web) is the place to extract δB/B0, β and M_A distributions.
