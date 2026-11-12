# Giuffrida, Miceli, Caprioli et al. (2022): The supernova remnant SN 1006 as a Galactic particle accelerator

- Ref: Nature Communications 13, 5098 (2022) (volume and article number as printed in the reference lists of Zhou et
  al. 2023 and McGinness et al. 2026); DOI 10.1038/s41467-022-32781-4; arXiv:2208.14491v1; peer-reviewed.
- Read: 2026-11-12, full (data points of Fig. 3 not read from the plot); topics: spherical-shock-injection

## Claim
The density of the shocked interstellar gas behind the forward shock of SN 1006 is lowest where the shock is
perpendicular to the ambient field and rises towards the quasi-parallel limb. The authors read this as a compression
ratio growing from 4 to about 7 through efficient hadron acceleration.

## Method
- Observation plus an analytic model; no kinetic simulation is run (hybrid results are quoted).
- Chandra (669.85 ks in 2012 plus 68.87 ks in 2008) and XMM-Newton (obs. 0555630201) spectra of nine narrow regions
  between contact discontinuity and forward shock, theta = 0-122 deg, with theta = 0 at the centre of the
  north-eastern radio limb (Methods; Fig. 1).
- Thermal non-equilibrium-ionisation component with kT fixed at 1.35 keV and N_H at 7e20 cm^-2, plus loss-limited
  synchrotron in regions +2 to +5. Density from the emission measure and a numerically computed volume (Methods).
- Compression ratio r_t = 4 imposed at theta about 90 deg (region 0); other regions scaled by density.
- Model: r_t from post-shock cosmic-ray and magnetic pressures including the "postcursor" (Haggerty & Caprioli 2020),
  each pressure modulated with obliquity (Methods, Eq. 1), dth = 20 deg:
  xi(theta) = (xi_i/2) [2 + tanh((th_i - theta)/dth) - tanh((pi - th_i - theta)/dth)].
  Protons from the thermal pool: xi_p = 12 %, th_p = 45 deg; re-accelerated Galactic seeds: xi_s = 6 %, th_s = 70 deg;
  magnetic pressure: xi_B = 5 %, th_B = 70 deg. The values are not fitted; they come from Caprioli & Spitkovsky
  (2014a,b), Caprioli, Pop & Spitkovsky (2015) and Caprioli, Zhang & Spitkovsky (2018).

## Numbers to compare against
- Post-shock density (Table 1, Chandra): 0.164 (+0.014, -0.016) cm^-3 in region 0 (theta about 90); 0.209 in +1; 0.20
  in +2; 0.21 in +3; 0.34 (+0.17, -0.13) in +4; 0.29 (+0.10, -0.07) in +5 (theta about 0). XMM-Newton: 0.206 in +3;
  0.26 (+0.07, -0.06) in +4-5. Errors are 68 %.
- Revised XMM values for theta = 53-120 deg: 0.150 to 0.206 cm^-3, minimum at 80-96 deg (Table 3).
- r_t "up to approximately 7" in quasi-parallel regions (Results; Fig. 3b).
- Model at theta = 0: xi_p + xi_s = 18 %, xi_B = 5 % give r_t = 6.34 and a spectrum E^-2.19; observed radio index 0.6
  corresponds to E^-2.2 (Discussion).
- Quoted from simulations: thermal-proton injection 10-15 % at quasi-parallel shocks, suppressed above about 45 deg;
  re-acceleration of seeds 2-6 %, up to about 60 deg (Discussion).
- SN 1006: shock speed 5000 km/s in region 5; distance 2.2 kpc; shock radius 14.4 arcmin (region +5) to 14.55 arcmin
  (region 0), variation below 0.15 arcmin; about 600 pc above the plane; time since shock passage in the regions
  100-200 yr (Results; Methods; Fig. 4).
- Ambient field: south-west to north-east; tilt 38 +/- 4 deg to the line of sight and a gradient of 1.5 B over 10 pc,
  both from Bocchino et al. (2011) (Introduction; Discussion).
- A density origin of the modulation would need a speed difference of about 1200 km/s and a radius difference of about
  1e18 cm (0.5 arcmin) in 250 yr, which is not seen (Results).

## Assessment
- Shown by observation: the post-shock density rises by a factor of about 1.6-2 from region 0 to +4/+5, in two
  instruments; errors are large where synchrotron dominates. Checks: region size, ejecta, thermal model, free kT, N_H.
- Inferred, not measured: that the rise is a compression ratio. This needs a uniform upstream density, argued from the
  round shock, faint uniform H-alpha and absent H I structure. The upstream density itself is not quoted.
- The model curve is not fitted; three parameter sets are overplotted (Fig. 3). A field-strength gradient is said to
  give a narrow minimum "remarkably similar" to the observed one (Discussion, Fig. 6).
- Electrons: no parametrisation and no statement on electron injection. The only electron evidence is that TeV
  electrons (X-ray synchrotron) sit in the quasi-parallel limbs while quasi-perpendicular regions show radio (GeV
  electrons) but no X-ray synchrotron (Discussion).
- Stated limits: theta is measured in the plane of the sky; a tilted field stacks different obliquities along the line
  of sight and lowers the contrast; regions with theta < 0 and > 122 deg are excluded (several shock fronts, ejecta
  knots, a cloud); refined field geometry is deferred.
- Noted alternative: the small shock-to-ejecta gap can also come from ejecta clumping (Orlando et al. 2012) (Intro.).
- Eq. 1 is garbled in the text file; the form above is rebuilt and reproduces the stated maxima at theta = 0 and pi.
