# Jebaraj et al. (2026): The Wave-Regulated Precursor of a Near-Parallel Interplanetary Shock Observed by Parker Solar Probe

- Ref: arXiv:2608.12606v1 (10 authors: Jebaraj, Colomban, Agapitov, Gedalin, Malkov, Kouloumvakos, Husidic,
  Mondal, ...). Preprint; the arXiv metadata gives no journal reference.
- Read: 2026-11-14, abstract, the event and input-value sections, the transport passages and the appendix
  passages on wave amplitude (arXiv HTML); the wave-family classification was skimmed;
  topics: qpar-electron-conditions

## Claim
At the near-parallel shock of 2023-03-13 (0.24 au) the upstream wave field separates into four families.
Three are cyclotron-resonant with the suprathermal-to-MeV protons streaming from the shock, which drive
them; a weak oblique compressive (fast-mode) component lies outside this loop. The precursor is
self-regulated and stays linear.

## Method
- PSP FIELDS and plasma data upstream of the ramp (crossed 07:13:13 UT); shock-normal distance from
  x = u_1 (t_shock - t); polarization and wavenumber analysis; energetic-particle pressure profile from
  IS(sun)IS protons.

## Numbers to compare against
- Inputs: shock speed ~2800 km/s; V_sw ~390 km/s; upstream speed in the shock frame u_1 ~2400 km/s;
  v_A ~320-360 km/s; c_s ~120-140 km/s; M_A ~7-7.5; T_i ~90-100 eV measured; T_e <~ 40 eV "unresolved";
  n ~30 cm^-3; B_0 ~80 nT; ion beta ~0.2. (The M_A of 9.1 in jebaraj2024acceleration is superseded.)
- Precursor: first wave activity ~06:00 UT; fully developed ("stationary") inner foreshock from ~06:15 UT
  to the ramp: ~58 min, ~12 solar radii (~0.06 au).
- Fluctuations stay linear; no SLAMS-like jumps. Total transverse amplitude ~10 %; 5-6 % at the 1-2.4 MeV
  proton resonances (App. B).
- Energetic-particle pressure e-folding length L_EP ~5.48 solar radii; kappa_eff = u_1 L_EP ~9.2e15 m^2/s
  (an average over the energies that dominate the pressure); lambda_par = 3 kappa_eff / v ~1.9-2.9 solar
  radii; from the resonant wave amplitude instead, lambda_par ~0.7-1.7 solar radii.
- Compressive fluctuations shift the resonance energies by up to 13 % along the precursor.

## Assessment
- Protons only; electrons are not treated. The measured precursor is the one electrons of the same event
  have to cross, so its amplitude and scale are inputs any electron model must meet.
- kappa_eff is one number for a range of energies; an energy-resolved mean free path is given in
  jebaraj2026diffusive.
- The spacecraft crossed the flank of the shock (dresing2025reason, chiappetta2026evolution), and other
  works give different parameters for this crossing (giacalone2025source: M_A ~4).
- Companions: kouloumvakos2026wave, jebaraj2026diffusive.
