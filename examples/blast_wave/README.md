# Blast wave around a magnetised body

A shell of dense electron–ion plasma (the debris) expands from the centre of an open 3-D box into a magnetised plasma at rest. It pushes the field and the ambient plasma outwards and leaves a magnetic cavity behind; the compression ahead of it steepens into a shock when the debris is stopped within a distance larger than the ion scales. An absorbing sphere with a magnetic dipole sits at the centre, inside the shell; the uniform field B₀ and the dipole moment are along +z.

- **Prior work the set-up follows**:
  - Le et al. 2021 (JGR Space Physics 126, e2021JA029125; arXiv:2109.00583; 2-D hybrid simulations) and Keenan et al. 2022 (Phys. Plasmas 29, 012107; arXiv:2112.12114; 3-D hybrid simulations): a debris cloud with a radial drift V_d in a uniform field. The parameters are the Alfvén Mach number M_A = V_d / v_A, the equal-mass radius R_m (the ambient mass inside R_m equals the debris mass), the directed gyroradius of the debris ions ρ_d = V_d / Ω_ci and the magnetic stopping radius R_B = R_m M_A^(2/3) (debris energy = field energy inside R_B). The debris ions couple to the ambient plasma through the fields when ρ_d < R_m and decouple from it when ρ_d > R_m (Le et al. 2021); ρ_d = 0.375 R_m here, as in Keenan et al. 2022.
  - Dieckmann et al. 2018 (Phys. Plasmas 25), 2-D full PIC: a dense shell expanding into a magnetised ambient plasma.
  - Nersisyan & Osipyan (arXiv:0905.0051): an explosion in a dipole field.
- **Differences to them**: full PIC instead of hybrid; a reduced mass ratio; debris ions with the mass of the ambient ions (the hybrid runs: three times that mass); a shell around the body instead of a Gaussian cloud at the centre (of width 0.1 R_m in Le et al. 2021, 4 d_i in Keenan et al. 2022); R_m = 8 d_i (10 to 150 d_i in Le et al. 2021, 40 d_i in Keenan et al. 2022); the body and the dipole at the centre are additions of this example. No PIC or hybrid work was found that drives a flow onto a magnetised sphere from all directions, so the shock is made the way the literature makes spherical shocks: outwards from the centre.
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

## Result of the gate configuration (8 GPUs, FP64, two runs with the same seed)

| quantity | measured | reference |
|---|---|---|
| speed of the front between 0.3 and 0.7 R_B | 0.1627 c = 2.10 fast speeds = 0.87 V_d | gate: above the fast speed, 0.0774 c |
| largest ring mean of B_z / B₀ there | 2.69 | gate: above 1.5; perpendicular MHD shock at the speed of the front: 2.26 |
| front radii, axes against diagonals | within 0.21 % | gate 5 % |
| largest extent of the debris ions Δ98 | 39.5 d_e = 1.24 R_m, at step 616 of 1760 | 32.6 d_e = 1.02 R_m (Eq. 4 of Le et al. 2021 with D = 3); gate ± 25 % |
| ambient n_e, n_i and \|B\| ahead of the front | within 0.05 % of step 0 | gate 5 % |
| Gauss's law, div B | 7e-15 and 8e-15, 9e-15 | gate 1e-10, 1e-12 |

- **The debris extent is 21 % above the scaling**, 4 points inside the band of the gate. The front speed, the ring mean, the symmetry and Δ98 are the same in the two runs in all digits; the other values differ from the fifth digit on (the plasma ahead of the front: 0.05003 % in both) and the round-off residual of Gauss's law in the first, so the 21 % are not a difference between runs. The scaling is a 2-D fit used with D = 3, the cloud here is a shell with the ambient ion mass, and R_m is 8 d_i; these runs do not separate the three. The two runs are the first run of the configuration and the gate record (`bench/results/`); the plan lists both.
- **The front leaves the debris.** It is 10.5 d_e ahead of Δ98 when it leaves the window (step 660) and at 77.4 d_e = 1.16 R_B at the end, where Δ98 has fallen to 34.3 d_e. The ring mean of B_z falls from 2.6 B₀ (step 352) to 1.5 B₀ at the end: the value of the table is the largest one of the window, not a mean behind a plane shock.
- **Magnetic cavity** (B_z < 0.5 B₀): its largest radius is 33.3 d_e = 0.50 R_B = 1.04 R_m, constant over the last 350 steps; the static field of the dipole alone gives 8.1 d_e.
- **Along B₀ the debris spreads further**: at the end the rms distance of the debris ions from the centre is 36.8 d_e along B₀ and 14.8 d_e per direction across it.
- **Energy** is monitored only (explicit scheme). The debris ions keep 20 % of their kinetic energy and the ambient ions gain 41 % of the kinetic energy of the debris. The budget with the outflow and the absorbed particles changes by +1.95 % of the energy at step 0, which is 10.9 % of the kinetic energy of the debris; the ambient plasma alone (null test N1 of the plan) gains 2.2 % of its energy over the same number of steps.
- **The body** absorbs 106,480 particles and holds the charge +22.5 at the end (largest value 32.9; units of e n₀ d_e³). Its potential against the plasma at two dipole radii is 4.9 T_e / e in the plane across B₀ and 1.3 T_e / e in the plane along it. These are reported without a limit: the null test with the body alone did not show a level of the charge.
- **Not shown**: whether the 21 % hold in another configuration (another R_m / d_i, mass ratio, shell shape, seed); a level of the body charge.
- **Cost**: 93 ms per step (0.80 ns per particle and step, times 8 GPUs), 164 s for the 1760 steps. The planner's memory figure is 21.5 GB per GPU for full particle arrays (66 GB for `full.yaml`); the arrays reach 67 % of their capacity (capacity factor 1.5). The device memory was not measured.
