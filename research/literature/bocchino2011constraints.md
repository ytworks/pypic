# Bocchino, Orlando, Miceli & Petruk (2011): Constraints on local interstellar magnetic field from non-thermal emission of SN1006

- Ref: Astronomy & Astrophysics ("in press" in the arXiv comment; volume and article number not in
  the metadata - journal reference not checked); DOI 10.1051/0004-6361/201016341; arXiv:1105.2689v1.
  Peer-reviewed article.
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
The 1 GHz morphology of SN 1006 is best reproduced when electron injection is most efficient at
quasi-parallel shocks, provided the ambient field strength has a gradient and the field makes an angle
of 38 +- 4 deg with the line of sight. The bright limbs are then polar caps tilted out of the sky plane.

## Method
3-D MHD simulation (FLASH; model of Orlando et al. 2007) of a point explosion (1.4 Msun, 1.3e51 erg)
in a uniform-density medium with a non-uniform field, stopped at 1000 yr; radio maps synthesised with
REMLIGHT. Not kinetic: m_i/m_e and particles per cell do not apply; Mach numbers, plasma beta, upstream
temperature and grid resolution: not stated. theta_Bn covers 0-90 deg over the sphere.
Emissivity i(nu) ~ K B_perp^(alpha+1) nu^-alpha, alpha = (s-1)/2, for N(E) = K E^-s, with
K ~ eps V_sh^-b and b = -1.5 (results unchanged for b between -1.5 and +2). Injection recipes (Sect. 2):
- quasi-parallel eps ~ cos^2(Theta_Bn2); quasi-perpendicular eps ~ sin^2(Theta_Bn2); isotropic:
  eps independent of obliquity (footnote 1, citing Ellison et al. 1995 and Voelk et al. 2003).
- Theta_Bn2 is the angle between the shock normal and the POST-shock field. Origin given: Fulbright &
  Reynolds (1990); alternative forms in Voelk et al. (2003) and Petruk et al. (2009a).
Five gradients (Table 1): |B| changes by a factor 1.0, 1.2, 1.4, 2.0 or 4.7 over 10 pc; the gradient
is across B for quasi-parallel and along B for quasi-perpendicular injection. 315 maps per model
(aspect angle phi_B 0-90 deg in 2 deg steps, gradient angle 0-90 deg in 15 deg steps) are compared
through four numbers measured in circles of 45 arcsec (5 % of the radius): A = rim max/rim min,
C = centre/rim max, R_max = ratio of the two limbs, theta_D = angle between the limb maxima.

## Numbers to compare against
- Set-up: n0 = 0.05 cm^-3; B = 30 muG at the explosion site, chosen high to mimic amplification (a few
  hundred muG post-shock); at 1000 yr V_sh = 4600 km/s and R = 9.6 pc (Sect. 2).
- Quoted for SN 1006: B about 100 muG, V_s about 4000 km/s, eta = 1-10 (Sect. 1); distance 2.2 kpc,
  Galactic latitude 14.6 deg, about 550 pc from the plane (Sect. 4).
- Acceleration time of electrons radiating at nu: 3e-3 eta (nu/GHz)^0.5 (B/100 muG)^-1.5
  (V_s/3000 km/s)^-2 yr (Sect. 1, from Uchiyama et al. 2007): under a year at 1 GHz, so the radio
  brightness as a function of obliquity is taken to trace injection alone.
- Observed: A = 2.7 +- 0.2, C = 0.36 +- 0.03, R_max = 1.2 +- 0.1, theta_D = 135 +- 10 deg (Sect. 3.2).
- Uniform field: quasi-parallel cannot match (A, C); quasi-perpendicular and isotropic can, near
  phi_B = 70 deg (Fig. 3, Sect. 3.3); but any uniform field gives R_max = 1 and theta_D = 180 deg.
- Best fit: quasi-parallel, gradient 1.4 per 10 pc (allowed 1.2-2.0), phi_B = 38 +- 4 deg, gradient
  15 +- 15 deg from the sky plane (Fig. 4, Sect. 3.4.1). Polar caps in the sky plane (90 deg) or along
  the line of sight (0 deg) are excluded.
- Quasi-perpendicular with the same gradient: acceptable (A, C) near phi_B = 70 deg, worse in
  theta_D-C, and R_max = 1 always (Fig. 6, Sect. 3.4.2).
- Rothenflug ratio R_pi/3: model 0.16 at phi_B = 90 deg, 1.00 at 38 deg; observed in the radio 1.01 (Sect. 4).
- The fitted gradient is much steeper than the large-scale Galactic field (factor 1.4 over more than
  100 pc); the random field component (scales of a few pc) is suggested instead (Sect. 4).

## Assessment
- Shown by forward modelling: a quasi-parallel recipe fits the radio limbs once a field gradient and a
  tilted aspect are allowed. This removes the radio argument against polar caps (Fulbright & Reynolds
  1990; Petruk et al. 2009b; Schneiter et al. 2010), which assumed a uniform field.
- Not shown: any mechanism. The cos^2 and sin^2 laws are assumed shapes; measuring in selected regions is
  said to reduce the sensitivity to their exact form. Isotropic injection is tested only in the uniform field.
- The quasi-perpendicular fit is "worse", not excluded; no goodness-of-fit statistic is given.
- Limits: amplification is mimicked by a 30 muG ambient field, the same at all obliquities; MHD
  without cosmic-ray back-reaction; uniform density (the NW cloud is left out).
- Footnote 5: most bilateral remnants have their axis parallel to the Galactic plane (Gaensler 1998);
  SN 1006 has it perpendicular and is called "rather peculiar".
