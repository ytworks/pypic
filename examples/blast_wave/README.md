# Blast wave around a magnetised body

A shell of dense electron–ion plasma (the debris) expands from the centre of an open 3-D box into a magnetised plasma at rest. It pushes the field and the ambient plasma outwards and leaves a magnetic cavity behind; the compression ahead of it steepens into a shock when the debris is stopped within a distance larger than the ion scales. An absorbing sphere with a magnetic dipole sits at the centre, inside the shell; the uniform field B₀ and the dipole moment are along +z.

- **Prior work the set-up follows**:
  - Le et al. 2021 (JGR Space Physics 126, e2021JA029125; arXiv:2109.00583) and Keenan et al. (arXiv:2112.12114), hybrid simulations: a debris cloud with a radial drift V_d in a uniform field. The parameters are the Alfvén Mach number M_A = V_d / v_A, the equal-mass radius R_m (the ambient mass inside R_m equals the debris mass), the directed gyroradius of the debris ions ρ_d = V_d / Ω_ci and the magnetic stopping radius R_B = R_m M_A^(2/3) (debris energy = field energy inside R_B). The shock leaves the debris when ρ_d < R_m.
  - Dieckmann et al. 2018 (Phys. Plasmas 25), 2-D full PIC: a dense shell expanding into a magnetised ambient plasma.
  - Nersisyan & Osipyan (arXiv:0905.0051): an explosion in a dipole field.
- **Differences to them**: full PIC instead of hybrid; a reduced mass ratio; a shell around the body instead of a Gaussian cloud at the centre; R_m = 8 d_i; the body and the dipole at the centre are additions of this example. No PIC or hybrid work was found that drives a flow onto a magnetised sphere from all directions, so the shock is made the way the literature makes spherical shocks: outwards from the centre.
- **Theory** (`theory.py`, `pypic/setup/mhd.py`):
  - R_m, R_B and the peak density of the shell that holds the debris mass.
  - The fast magnetosonic speed across B₀ and the compression of a perpendicular MHD shock (γ = 5/3).
  - The largest extent of the debris ions across the field, Eq. 4 of Le et al. 2021: Δ_max / R_m ∼ (1 + c₀ M_A) / (1 + c₁ / M_A²)^(1/D) with c₀ = 0.03 and c₁ = 2. The coefficients are a fit to 2-D hybrid runs (D = 2) with a Gaussian cloud and a debris-to-ambient mass ratio of 3; for D = 3 the paper gives the form and no fit. The comparison here is against that extrapolation.

| config | debris | grid | use |
|---|---|---|---|
| `gate.yaml` (gate) | m_i/m_e = 16, V_d = 0.1875 c, R_m = 32 d_e, R_B = 66.6 d_e | 384³ (half-width 1.44 R_B), 1760 steps | front, symmetry, debris extent, conservation |
| `full.yaml` | m_i/m_e = 25, V_d = 0.15 c, R_m = 40 d_e, R_B = 83.2 d_e | 576³ (half-width 1.73 R_B), 2800 steps | production |

Both: M_A = 3 (fast Mach number 2.4), R_m = 8 d_i, ρ_d = 0.375 R_m, ω_pe/Ω_ce = 4, v_te = 0.1 c (Debye length 0.2 Δx), T_i = T_e (β = 0.64), 8 particles per cell and beam, run length 1.25 R_B / V_d ≈ 7 / Ω_ci.

```bash
python -m pypic.launch -n 8 examples/blast_wave/run.py --config examples/blast_wave/configs/gate.yaml
python examples/blast_wave/analyze.py output/blast_wave/gate
```

- **Set-up** (explicit scheme; engine features PIC-D12 to PIC-D17 in `docs/design/engine.md`):
  - Four species. The ambient electrons and ions fill the box at rest and enter through the six open faces. The debris electrons and ions are loaded as a shell (Gaussian in the radius: centre 0.33 R_m, width 0.1 R_m, peak density 16.5 n₀) with the radial drift V_d and the ambient temperature; they are not injected. Each debris species splits its density over four beams, which lowers the weight of a debris particle to 4.1 ambient particles.
  - B₀ and the dipole are external fields; the grid holds the difference to them, and the outputs hold the total fields. The equatorial field of the dipole is opposite to B₀: the total field vanishes at r_dip = 0.2 R_m and the sphere (radius r_dip / 1.5) lies inside that radius. The energy of the dipole field outside the sphere is about 0.6 % of the kinetic energy of the debris: the dipole acts near the body, the far front is set by B₀.
  - The sphere removes the particles that enter it and keeps their charge (an insulating body at floating potential).
  - A compensated binomial current filter (4 passes) and the Marder zone at the open faces.
  - The planner refuses: Ω_ce Δt > 2 at the surface of the body; ω_pe Δt > 1.2 at the peak of the shell; a shell with more than 1 % of its particles inside the body; a body radius below 6 cells; a half-width below 1.3 R_B; more particle slots per species and GPU than the kernels index; particle ids beyond 2⁵³.
  - `--set debris=false` and `--set body=false` switch the shell and the body off (the null tests of the plan).
- **Gate** (definitions in `docs/plans/2026-11-11-blast-wave.md`, fixed before the first run), measured in the plane through the centre across B₀:
  - Gauss's law ≤ 1e-10, div B ≤ 1e-12;
  - the front (outermost radius with B_z > 1.25 B₀ on eight rays) moves faster than the fast magnetosonic speed while it is between 0.3 and 0.7 R_B, and the largest ring mean of B_z there exceeds 1.5 B₀;
  - the front radii on the axes and on the diagonals agree within 5 %;
  - the largest extent of the debris ions (98 % of the line integral of their density, Eq. 3 of Le et al. 2021) lies within 25 % of Eq. 4 with D = 3 and is reached before the end of the run;
  - the ambient densities and |B| ahead of the front stay within 5 % of their initial values.
- **Reported**: the compression against the MHD shock; the distance between the front and the debris; the extent of the debris along and across B₀; the radius of the magnetic cavity; the charge and the potential of the body; the energy budget; the largest fill ratio of the particle arrays.
