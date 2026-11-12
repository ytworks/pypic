# Bohdan, Pohl, Niemiec, Morris, Matsumoto, Amano, Hoshino & Sulaiman (2021): Magnetic field amplification by the Weibel instability at planetary and astrophysical high-Mach-number shocks

- Ref: Phys. Rev. Lett. 126, 095101 (2021) (volume/article from the DOI; meta.txt has no journal_ref);
  DOI 10.1103/PhysRevLett.126.095101; arXiv:2102.04328 (v1, "Accepted to PRL", read here); peer-reviewed
- Read: 2026-11-12, full (a letter; figures not readable in the text version); topics: spherical-shock-injection

## Claim
At perpendicular high-Mach-number shocks the field in the shock transition is amplified beyond compression by the
ion-ion Weibel instability, and the amplification depends on MA alone: |Bsh| ≈ 2 sqrt(MA) B0. Upstream β, mi/me
and shock speed do not matter. Mock spacecraft crossings of the simulated shocks reproduce the peak fields that
Cassini measured at Saturn's bow shock, so the scaling is proposed for SNR shocks as well.

## Method
- PIC, 2D3V, TRISTAN-based code; flow-flow set-up, beams at ±0.2c; θBn = 90°; B0 in the simulation plane (φ = 0°).
- Table I, runs (mi/me; MA; Ms cold/warm): A (50; 22.6; 1104/35), B (100; 31.8; 1550/49), C (100; 46; 2242/71),
  D (200; 32; 1550/49), E (200; 44.9; 2191/69), F (400; 68.7; 3353/106), G (50; 68.7; 3353/106),
  H (50; 100; 4870/154), I (50; 150; 7336/232). βe = 5·10⁻⁴ (runs *1) and 0.5 (*2). H1 and I1 are not used
  (numerical noise upstream). Mach numbers in the upstream frame.
- Shock speed 0.263c in the upstream frame; r ≈ 4; ωpe/Ωe = 8.5-80; λse = 20Δ; time step 1/40 ωpe⁻¹; 20 particles
  per cell per species; Ly = (8-24) λsi; run time ≈ 8 Ωi⁻¹.
- "Shock region" = slab of width rgi,up/3 = MA λsi/3 centred on the overshoot; fields averaged over reformation.
- Linear Weibel growth rates (method of Paper III) for five mass ratios and for vsh = 0.26c, 0.026c, 0.0026c with
  the foot parameters of run F2 rescaled (Fig. 3).
- Data: 422 quasi-perpendicular (θBn ≥ 45°) Cassini crossings of Saturn's bow shock (Sulaiman et al. 2016); virtual
  spacecraft on straight paths through the simulations, averaged over crossing point, speed and direction (Fig. 4).

## Numbers to compare against
- Earlier model (Leroy 1983 with hybrid runs): Bover ≈ 0.4 B0 MA^(7/6) (Eq. 1).
- This work: |Bsh| ≈ 2 sqrt(MA) B0 in the shock region (Eq. 2, Fig. 2a); magnetic over upstream ion energy density
  ≈ 4/MA (Eq. 3, Fig. 2b); peak along a crossing Bmax/B0 = 5.5 (sqrt(MA) − 2), said to describe the Cassini data
  for MA ≳ 10 (Eq. 5, Fig. 4).
- Density overshoot Nover/N0 ≈ 7 in all runs, with no MA dependence; the field rises about twice as much (Fig. 1b).
- Floor: unmagnetised-shock value UB = 0.006 Ui (Kato & Takabe 2008) is met at MA ≈ 670, |Bsh| ≈ 50 B0.
- Weibel growth rate ≈ 0.1 ωpi; constant within ~10 % over five mass ratios; Γmax ∝ vsh ωpi, i.e. ∝ MA Ωi (Eq. 4).
- Far downstream field 4 B0; amplified field lasts a few ion gyroradii behind the shock.
- Reforming shocks at Saturn (16 of 54 crossings, Sulaiman et al. 2015) have Bmax/B0 1.42 times higher; simulated
  peaks are divided by 1.26 before the comparison. Reformation period ≈ 1.5 Ωi⁻¹ in runs and in the Cassini data.
- Injection momentum for DSA pinj ∝ rsh Bsh ∝ sqrt(MA) (Eq. 6). Synchrotron emissivity up by (Bsh/B0)², "may
  reach a thousand".

## Assessment
- Shown by simulation: the sqrt(MA) scaling over MA = 22.6-150 at fixed shock speed; Bx and Bz (Weibel folding)
  grow beyond compression while By follows the density; no visible effect of βe or mi/me (runs with equal MA and
  different mass ratio agree).
- Shown only by linear theory: independence of shock speed (two orders of magnitude between runs and Saturn).
- Agreement with Cassini is stated from Fig. 4 without a quantitative measure; Eq. 5 is a fit to the simulation
  points; the factor 1.26 assumes that all simulated shocks reform and 16/54 of the observed ones do.
- Stated or conjectured: validity for SNR shocks; the stronger field shortens the electron mean free path and
  raises the gyrofrequency, so the cut-off of stochastic shock drift acceleration at quasi-perpendicular shocks
  "may depend on MA"; injection into DSA needs a larger momentum (Eq. 6). Nothing on electron spectra is computed.
- Limits stated: runs too short for the downstream relaxation at high MA; in-plane 2-D is taken as a proxy for 3-D
  (Matsumoto et al. 2017; Bohdan et al. 2017); the exponential stage is not seen, the field "defies an analytical
  derivation"; below MA ≈ 20 the fit is not evidence that the Weibel mode sets the field.
- Disagreement: Matsumoto et al. (2012) confirmed Eq. 1 in out-of-plane 2-D, a geometry that misses the Weibel
  mode; agreement of Eq. 1 with the data below MA ≈ 60 is called possibly a coincidence. The constant Nover is
  "not in line" with Leroy's model.
- Doubt: the text says the number of Weibel growth cycles is proportional to MA and, further on, that it "scales
  inversely" with MA; Eq. 4 with a time ~Ωi⁻¹ supports the first.
- Perpendicular only in the runs; the observations cover θBn ≥ 45° without an obliquity analysis.
