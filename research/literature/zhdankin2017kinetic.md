# Zhdankin, Werner, Uzdensky & Begelman (2017): Kinetic turbulence in relativistic plasma: from thermal bath to nonthermal continuum

- Ref: Phys. Rev. Lett. 118, 055103 (2017) (Crossref); DOI 10.1103/PhysRevLett.118.055103;
  arXiv:1609.04851; peer-reviewed.
- Read: 2026-11-13, skimmed (ar5iv HTML, text and methods); topics: shock-alfvenic-turbulence

## Claim
Driven kinetic turbulence in a closed box of relativistically hot pair plasma turns an initial thermal
bath into a nonthermal power law; the index hardens with magnetization and steepens with system size.

## Method
- Zeltron (explicit EM PIC), 3-D periodic cube, B0 along z; N^3 = 256^3 to 1024^3 (fiducial 768^3,
  ~1e11 particles), rho_e0 >= 1.5 dx, >= 128 ppc for converged spectra.
- Pair plasma, Maxwell-Juttner with mean gamma ~ 300 (ultra-relativistic).
- Driving: Langevin antenna (tenbarge2014oscillating) as J_ext; J_z at 8 modes, J_x, J_y at 4 each;
  omega0 = 0.6 * 2 pi v_A0/(sqrt3 L), Gamma0 = 0.5 * 2 pi v_A0/(sqrt3 L); amplitude tuned to delta B_rms ~ B0.
- sigma0 = 0.25-4; L/(2 pi rho_e0) = 27-109 (fiducial 61.1); duration >= 10 L/c (fiducial ~22 L/c).

## Numbers to compare against
- Index alpha = 3.44, 3.05, 2.43, 2.04, 1.55 for sigma0 = 1/4, 1/2, 1, 2, 4 (Fig. 4, 512^3).
- Fit alpha ~ 1 + C0 xi0^(-1/2), C0 ~ 0.075, xi0 = sigma0 rho_e0 / L (Fig. 4); steeper for larger L/rho_e0.
- Energy injection per volume ~ B0^2 c/(8 pi L) per unit time; internal energy grows linearly (Fig. 2).
- Magnetic spectrum ~ k_perp^-5/3 inertial, ~ k_perp^-4 below the Larmor scale (Fig. 3).

## Assessment
- Relevance to us is indirect: relativistic pairs, no ions, v_A ~ c. It shows that delta B ~ B0
  turbulence alone makes a power law in a closed box when sigma >~ 1, and that the tail softens fast
  as sigma falls (alpha ~ 3.4 at sigma = 1/4). Our nonrelativistic e-i upstream has sigma ~ (v_A/c)^2
  << 1, so the expected turbulence-only electron tail is far steeper; that is a hypothesis to test
  with the null run, not a result of this paper.
- Driven heating is linear in time: a null run with driving has a known, controllable heating rate.
- Limit: no nonthermal fraction is quoted; spectra keep steepening with box size, so the asymptote
  is open.
