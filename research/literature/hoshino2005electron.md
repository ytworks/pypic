# Hoshino (2005): Electron Surfing Acceleration in Magnetic Reconnection

- Ref: J. Geophys. Res. Space Physics (2005), volume/article number not checked; DOI 10.1029/2005JA011229; arXiv:astro-ph/0507528
- Read: 2026-11-12, full; topics: (none yet)

## Claim
In externally driven reconnection of an ion-scale-thin plasma sheet, Hall-related polarization (ambipolar)
fields at the lobe/plasma-sheet boundary become much stronger than the reconnection field. Lobe electrons
trapped in this potential well are surfed along the inductive field Ey on their way to the X-line. This
"first-stage" pre-acceleration feeds the Speiser (second) and ∇B/curvature pile-up (third) stages of the
author's earlier work, and the author argues that the three stages together can quickly produce MeV electrons.

## Method
- 2D PIC, driven reconnection. Doubly periodic Lx = 38.4λ, Lz = 32λ (λ = initial sheet thickness),
  plasma injection layers at z = ±12λ with E_dri/E0 = α1/cosh²(x/2λ) + α2. The strongly driven case
  α1 = α2 = 0.2 is shown; the weakly driven case (0.1) evolves about 2× slower with similar structure (Sect. 2).
- 1536×1280 cells, ~4e8 particles, mi/me = 64, c/VA = 16.8, c/ωpi = 1.15λ (46 cells), c/ωpe = 0.14λ
  (6 cells), Ti/Te = 4 in sheet, lobe density 5 %, VA/vth,i = 1.15 (Sect. 2).
- Diagnostics: Ex, Ez maps (Fig. 1); pitch-angle-integrated spectra at t/τA = 7.5, 44.7, 52.2, 59.6 (Figs. 2, 7);
  one tracked electron, orbit and energy (Fig. 3), (Py, Pz) hodogram (Fig. 6); cartoons (Figs. 4, 5).
- Analytic model: sheet thinning under driving (Eqs. 1–2); surfing force balance v_max ~ c Ez'/Bx' (Eq. 6);
  polarization field from eφ ~ Ti over the hybrid scale c/(ωpe ωpi)^(1/2) (Eqs. 7–12); surfing maximum
  energy from trapping time along the boundary (Eqs. 13–15).

## Numbers to compare against
- Polarization field ≈ 10–20 E_dri, width below c/ωpi (Sect. 2); E/E0 ≈ 5–10 driven vs ≈ 1 spontaneous (Sect. 3).
- Nonthermal/thermal electron energy density ≈ 10 % driven vs ≈ 1 % in the spontaneous run of the earlier
  paper (Sect. 3, Fig. 2).
- Tracked electron: ε/me c² ≈ 0.3 after surfing, ≈ 0.6 after Speiser, ≈ 0.7 in pile-up (Sect. 3, Fig. 3);
  (Py, Pz) ≈ (0.8, −0.2) after stage 1 (Fig. 6).
- Magnetotail estimate Ep ≈ 13 mV/m (l = c/ωpi) to ≈ 540 mV/m (l = c/ωpe) (Eqs. 10–11); Ep/(B0VA/c) ~ c/(ωpi l) (Eq. 12).
- Surfing maximum energy ~ O(1e-2) me c² ("several keV") for tail parameters (Eq. 15); ≈ 0.8 me c² for the
  run (θ = 4, (ωce/ωpe)² = 0.05, Lωpi/c ≈ 2, ve0 ≈ 0.5 vth,e) vs spectral break ≈ 0.7 me c² (Sect. 4, Fig. 7).
- Double power law: index ≈ 3 below and ≈ 7 above the break, t/τA = 59.6 (Fig. 7).

## Assessment
- Shown (simulation): strong V-shaped boundary fields and outward-moving bipolar wave trains under driving;
  more nonthermal electrons than in the spontaneous case; one electron whose first gain is along Py,
  perpendicular to both B and Ep (Figs. 3, 6).
- Shown (theory vs simulation): the Eq. 15 estimate is close to the spectral break. The low-energy segment is
  attributed to surfing on this coincidence plus one orbit: no tail statistics, no switch-off test. Eq. 15
  scales as ve0^(-2) (from Eqs. 13–14), and ve0 is assumed, so the match is order-of-magnitude. My
  re-evaluation with the squared factor gives ≈ 1.6 rather than 0.8; the text extraction of Eq. 15 is ambiguous,
  so check this against the published version.
- Stated/conjectured: MeV electrons in the tail and in flares via three stages (not simulated at real parameters);
  links to RHESSI double power laws and kappa distributions; in-situ test via spiky boundary fields (Sect. 4).
- Stated limits: mi/me = 64 "to save" CPU time, so validity at general parameters needs care (Sect. 4); the
  trapping condition along the path needs closer analysis (Sect. 4); boundary placement at ±12λ was chosen to
  separate the driven front from the X-line (Sect. 2). 2D geometry and the periodic box are not discussed as limits.
- Open: trapping along a varying Bx', observational confirmation of the polarization-field structure.
