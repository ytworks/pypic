**English** · [日本語](Blast-Wave-ja) · source: [`examples/blast_wave`](https://github.com/ytworks/pypic/tree/main/examples/blast_wave)

# Blast wave around a magnetised body

A shell of dense electron–ion plasma (the debris) expands from the centre of an open 3-D box into a
magnetised plasma at rest (Alfvén Mach number {{ c.gate.params.mach_alfven }}, m_i/m_e =
{{ c.gate.params.mass_ratio }}). It pushes the field and the ambient plasma outwards and leaves a
magnetic cavity behind. An absorbing sphere with a magnetic dipole sits at the centre, inside the shell.
In this run the front of the compressed field moves at {{ r.gate.front_mach_fast | .2f }} times the fast
magnetosonic speed, and the debris ions reach {{ r.gate.debris_extent | .1f }} d_e across the field,
{{ r.gate.debris_extent_ratio | .2f }} times the scaling of Le et al. (2021).

## Physics

**Scales of the expansion** (Le et al. 2021). A debris cloud of mass M_d leaves the centre at the
speed V_d into a plasma of density n₀ in a uniform field B₀. The equal-mass radius R_m is where the
ambient mass inside equals the debris mass, and the magnetic stopping radius R_B is where the field
energy inside equals the kinetic energy of the debris:

```math
\frac{4\pi}{3} R_m^3\, n_0 m_i = M_d , \qquad
\frac{4\pi}{3} R_B^3\, \frac{B_0^2}{2} = \frac{1}{2} M_d V_d^2
\;\;\Rightarrow\;\; R_B = R_m M_A^{2/3} , \quad M_A = \frac{V_d}{v_A} .
```

The debris ions couple to the ambient plasma through the fields when their directed gyroradius
ρ_d = V_d/Ω_ci is smaller than R_m; for ρ_d > R_m they decouple from it (Le et al. 2021).

**Extent of the debris.** Le et al. (2021) measure the radius Δ98 across the field that holds 98 % of
the line-integrated debris density and fit its largest value by

```math
\frac{\Delta_{max}}{R_m} \sim \frac{1 + c_0 M_A}{\left(1 + c_1 / M_A^2\right)^{1/D}} ,
\qquad c_0 = 0.03 , \quad c_1 = 2 .
```

The coefficients are a fit to 2-D hybrid runs (D = 2) with a Gaussian cloud and a debris-to-ambient
ion mass ratio of 3. For D = 3 the paper gives the form and no fit; the comparison below is against
that extrapolation.

**Front.** The compression ahead of the debris moves across B₀ and is compared with the fast
magnetosonic speed and with the compression X of a perpendicular MHD shock:

```math
v_f^2 = v_A^2 + \frac{\gamma (T_e + T_i)}{m_i} , \qquad
(2-\gamma) X^2 + \left[\gamma(\beta+1) + (\gamma-1) M^2\right] X - (\gamma+1) M^2 = 0 , \qquad \gamma = 5/3 ,
```

with M the Alfvén Mach number of the front.

| quantity | value |
|---|---|
| debris speed V_d | {{ c.gate.planner.v_d | .4f }} c |
| Alfvén speed, fast magnetosonic speed | {{ c.gate.planner.v_a | .4f }} c, {{ c.gate.planner.v_fast | .4f }} c |
| plasma β (electrons and ions) | {{ c.gate.planner.beta | .2f }} |
| R_m, R_B | {{ c.gate.planner.r_m | .1f }} d_e, {{ c.gate.planner.r_b | .1f }} d_e |
| directed gyroradius of the debris ions ρ_d | {{ c.gate.planner.debris_gyroradius | .1f }} d_e |
| Δ_max of the scaling with D = 3 | {{ c.gate.planner.debris_extent_le | .1f }} d_e |

## Experimental conditions

Normalised units: c = ε₀ = 1, lengths in d_e = c/ω_pe, times in 1/ω_pe, the ambient density n₀ = 1.
B₀ and the dipole moment are along +z.

- **Four species.** The ambient electrons and ions fill the box at rest and enter through the six open
  faces. The debris electrons and ions are loaded as a shell, Gaussian in the radius (centre
  {{ c.gate.planner.shell_radius | .1f }} d_e, width {{ c.gate.planner.shell_width | .1f }} d_e, peak
  density {{ c.gate.planner.shell_peak_density | .1f }} n₀), with the radial drift V_d and the ambient
  temperature. The debris ions have the mass of the ambient ions.
- **External fields**: B₀ and the dipole are external fields; the grid holds the difference to them and
  the outputs hold the total fields. The equatorial field of the dipole is opposite to B₀, and the total
  field vanishes at {{ c.gate.planner.dipole_radius | .1f }} d_e from the centre.
- **Sphere**: radius {{ c.gate.planner.body_radius | .2f }} d_e. It removes the particles that enter it
  and keeps their charge, as an insulating body at floating potential.
- **Run length**: 1.25 R_B/V_d, which is Ω_ci t = {{ c.gate.planner.omega_ci_t_end | .1f }}.

{{ table:conditions }}

<details><summary>Planner parameters</summary>

{{ table:params }}

</details>

## Animation

![blast wave](images/blast_wave/gate.gif)

The sphere is blanked.

- **Left, B_z / B₀ in the plane z = centre** (the plane across B₀): the compressed field at the front
  and the cavity behind it.
- **Middle, ambient ion density in the plane z = centre.**
- **Right, debris ion density in the plane y = centre** (the plane that holds B₀), logarithmic colour
  scale.

## Comparison with the scalings

![front, debris and cavity](images/blast_wave/gate-blast_wave.png)

B_z in the plane across B₀ with the front on eight rays, the debris ions in a plane along B₀, and the
radii of the front, the debris and the cavity against time. The front is the outermost radius with
B_z > 1.25 B₀. The front criteria use the samples with the front between 0.3 and 0.7 R_B.

| quantity | PIC | reference | gate |
|---|---|---|---|
| front speed in the window | {{ r.gate.front_speed | .4f }} c | {{ r.gate.fast_speed | .4f }} c (fast magnetosonic speed) | faster |
| largest ring mean of B_z / B₀ in the window | {{ r.gate.compression | .2f }} | — | > 1.5 |
| front radii, axes against diagonals | {{ r.gate.asymmetry | .2% }} | — | within 5 % |
| largest Δ98 of the debris ions | {{ r.gate.debris_extent | .1f }} d_e | {{ r.gate.debris_extent_le | .1f }} d_e (Le et al. 2021, D = 3) | within 25 %, before the end |
| ambient n_e, n_i and \|B\| ahead of the front | within {{ r.gate.ambient_deviation | .2% }} of step 0 | — | 5 % |
| Gauss's law, max\|∇·E − ρ\| | {{ r.gate.gauss_max | .0e }} | — | ≤ 1e-10 |
| max\|∇·B\| | {{ r.gate.divB_max | .0e }} | — | ≤ 1e-12 |
| compression of a perpendicular MHD shock at the measured front speed | — | {{ r.gate.compression_mhd_front | .2f }} | reported |
| largest radius of the magnetic cavity (B_z < 0.5 B₀) | {{ r.gate.cavity_max | .2f }} d_e | — | reported |

The gate criteria were fixed before the first run.

## Notes

- **The debris extent lies above the scaling.** The scaling is a 2-D fit used with D = 3, the cloud
  here is a shell with the ambient ion mass, and R_m is 8 d_i. The run does not separate these.
- **The front leaves the debris.** At the end of the run the front is at {{ r.gate.front_end | .1f }} d_e,
  beyond R_B, while the debris ions stay inside their largest extent across the field. Along B₀ the
  debris spreads further: its rms distance from the centre is {{ r.gate.debris_rms_along | .1f }} d_e
  along B₀ and {{ r.gate.debris_rms_across | .1f }} d_e across.
- **The compressed field weakens as the front moves out**; the value in the table is the largest one in
  the window.
- **Energy.** At the end the debris ions keep {{ r.gate.debris_ion_energy_left | .0% }} of their kinetic
  energy. The ambient ions have gained {{ r.gate.ambient_ion_energy_gain | .0% }} of the kinetic energy
  that the debris (electrons and ions) had at step 0. The total budget is monitored only (explicit
  scheme); it changes by {{ r.gate.energy_budget_change | +.1% }} of that debris energy.
- **The body** has absorbed {{ r.gate.n_absorbed }} particles and holds the charge
  {{ r.gate.body_charge | +.1f }} (units of e n₀ d_e³) at the end. A static plasma test without debris
  did not show a level of the charge; the potential of the body is reported, not gated.
- **Particle number**: the particle arrays reach {{ r.gate.fill_max | .0% }} of their capacity on the
  fullest rank.
- **A larger R_m / d_i** (`rm12`, a recorded run): R_m = {{ c.rm12.params.rm }} d_i against
  {{ c.gate.params.rm }} d_i in the run above, {{ c.rm12.particles | , }} particles. The debris ions reach
  {{ r.rm12.debris_extent | .1f }} d_e across the field, {{ r.rm12.debris_extent_ratio | .2f }} times the
  scaling of Le et al. (2021), against {{ r.gate.debris_extent_ratio | .2f }} times in the run above. The
  front moves at {{ r.rm12.front_mach_fast | .2f }} times the fast magnetosonic speed and the largest ring
  mean of B_z is {{ r.rm12.compression | .2f }} B₀. The run passes the gate criteria.

  <details><summary>Front, debris and cavity of `rm12`</summary>

  ![rm12](images/blast_wave/rm12-blast_wave.png)

  </details>
- **Run-to-run differences**: this run and three earlier FP64 runs of the configuration with the same
  seed give the same front speed, ring mean and symmetry in all digits and the same Δ98 in 12 digits;
  the other values differ by 4e-5 of their value or less. A run with `precision: mixed` passes the same
  criteria; it is in the gate records of the repository.

## Reproduce

```bash
python -m pypic.launch -n 8 examples/blast_wave/run.py \
    --config examples/blast_wave/configs/gate.yaml --out output/wiki/blast_wave/gate
python examples/blast_wave/analyze.py output/wiki/blast_wave/gate --movie
```

## Data provenance

{{ table:provenance }}

## References

- A. Le, D. Winske, A. Stanier, W. Daughton, M. Cowee, B. Wetherton and F. Guo, "Astrophysical
  explosions revisited: collisionless coupling of debris to magnetized plasma", J. Geophys. Res. Space
  Physics 126, e2021JA029125 (2021).
- B. D. Keenan, A. Le, D. Winske, A. Stanier, B. Wetherton, M. Cowee and F. Guo, "Hybrid
  particle-in-cell simulations of electromagnetic coupling and waves from streaming burst debris",
  Phys. Plasmas 29, 012107 (2022).
- M. E. Dieckmann et al., "Expansion of a radially symmetric blast shell into a uniformly magnetized
  plasma", Phys. Plasmas 25, 052108 (2018).
- H. B. Nersisyan and D. A. Osipyan, "Collisionless plasma expansion in the presence of a dipole
  magnetic field", Contrib. Plasma Phys. 49, 351 (2009).
