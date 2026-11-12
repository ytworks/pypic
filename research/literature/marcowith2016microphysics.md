# Marcowith, Bret, Bykov et al. (2016): The microphysics of collisionless shock waves

- Ref: Rep. Prog. Phys. 79, 046901 (2016); DOI 10.1088/0034-4885/79/4/046901; arXiv:1604.00318;
  peer-reviewed review (the arXiv text is the accepted manuscript).
- Read: 2026-11-12, sections 2.6 (acceleration processes), 3.1.2-3.1.3 (non-stationarity, rippling,
  particle energization), 3.2.1 and 3.2.3 (radio and X-ray observations of remnants), 3.4.1 (injection in
  PIC and hybrid simulations), start of 3.5, summary; topics: spherical-shock-injection

## Claim
Review. On electron injection at non-relativistic shocks it states that diffusive shock acceleration is
the accepted final mechanism, while the injection process, and how its efficiency depends on the shock
conditions, is open (Sect. 3.4.1).

## Method
Review of theory, simulations (PIC, hybrid, Vlasov-Fokker-Planck, PIC-MHD) and observations. No own runs.
Mechanisms defined in Sect. 2.6: diffusive acceleration, shock drift (drift along the convective electric
field caused by the field gradient "or the curvature of the shock front"), shock surfing.

## Numbers to compare against
- Low M_A (~ 3), high beta_e (> ~ 20), clusters and flares: electron acceleration by shock drift
  (Park et al. 2013; Guo, Sironi & Narayan 2014a,b) (Sect. 3.4.1).
- High M_A: surfing with Buneman waves proposed for M_A > (m_i/m_e)^(2/3) ~ 150 at quasi-perpendicular
  shocks (Matsumoto et al. 2012); it yields non-thermal electrons but has not been shown to reach an
  electron gyroradius comparable to the ion one (Sect. 3.4.1).
- Whistler waves from the modified two-stream instability in the foot (Riquelme & Spitkovsky 2011) do
  reach that gyroradius; power-law tail of index about -3; most efficient for M_A below about 20,
  the instability needing M_A < (m_i/m_e)^(1/2) (Sect. 3.4.1).
- Ions (hybrid, Caprioli & Spitkovsky 2014): up to about 15 % of the energy, best for theta < 45 deg and
  high M_A; amplification factor scales as M_A^(1/2) (Sect. 3.4.1).
- Injection criterion used: a sizeable fraction of electrons must reach a gyroradius comparable to the
  ion gyroradius, which sets the shock width (Sect. 3.4.1).
- Radio spectral index of remnants: mean 0.5, dispersion about 0.2 (Sect. 3.2.1). X-ray rims: electrons
  above about 10 TeV at the forward shock; SN 1006 edges, 100 TeV electrons (Sect. 3.2.3).
- Reformation period below one upstream ion gyroperiod; strongest for beta_i much less than 1; emitted
  whistlers have wavelengths of one to a few ion inertial lengths along the front (Sect. 3.1.2).

## Assessment
- Listed as open (Sect. 3.4.1): (1) the dependence of electron injection on shock conditions; (2) for the
  whistler route, the passage to diffusive acceleration with a p^-4 spectrum is not shown; (3) whether
  high-M_A quasi-parallel shocks accelerate electrons (suggested by Masters et al. 2013, proposed by
  Levinson 1992, 1994); (4) a conjecture that ion-driven field amplification turns parts of a quasi-parallel
  shock into low-Mach, quasi-perpendicular patches where whistler injection could work.
- Stated tension: electrons favour quasi-perpendicular shocks of moderate M_A (< ~ 20), ions favour
  high-M_A quasi-parallel shocks (Sect. 3.4.1).
- Rippling and non-stationarity (Sect. 3.1.2-3.1.3): small-scale rippling from lower-hybrid drift type
  instabilities in the ramp (Lembege & Savoini 1992; also Lowe & Burgess 2003, Moullard et al. 2006);
  large-scale rippling from ion cyclotron and mirror instabilities (Winske & Quest 1988). In PIC runs
  reformation gives bursts of reflected electrons in time and rippling gives "packs" in space (Lembege &
  Savoini 2002); whether these survive foot micro-turbulence is open. For ions, drift acceleration stays
  dominant with rippling (Yang et al. 2012).
- Curved shock: only the foreshock of a curved front in 2-D PIC (Savoini & Lembege 2001; Savoini et al.
  2013), whose outer edge is the field line tangent to the shock (Sect. 3.1.2). No electron acceleration
  result for a curved shock is reviewed. The reference list has neither Burgess (2006) nor Guo &
  Giacalone (2010, 2015): the hybrid plus test-electron line of work is absent from this review.
- Numerical constraints stated: ion injection in PIC is hard at even moderately realistic mass ratios
  (Sect. 3.4.1); the thresholds above scale with m_i/m_e; reformation against whistler emission depends on
  dimension and field orientation, with one 3-D run at the real mass ratio (Shinohara et al. 2011).
- Doubts: says the highest M_A in Caprioli & Spitkovsky (2014) was 100, that paper has 50; its ref. 219
  is printed as "Astrophys. Journ. 749 46" (Pohl et al. 2020 print ApJ 794, 46 for the same paper).
