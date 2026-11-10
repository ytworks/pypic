**English** · [日本語](Bow-Shock-ja) · source: [`examples/bow_shock`](https://github.com/ytworks/pypic/tree/main/examples/bow_shock)

# Bow shock of a magnetised planet

A magnetised electron–ion wind (Alfvén Mach number {{ c.gate.params.mach_alfven }}, m_i/m_e =
{{ c.gate.params.mass_ratio }}) flows through an open 3-D box and meets an absorbing sphere with a
magnetic dipole at its centre. The dipole field stops the wind at the magnetopause. The wind is faster
than the fast magnetosonic speed, so a bow shock stands upstream of the magnetopause, with the sheath of
compressed plasma in between. In this run the magnetopause lies at {{ r.gate.magnetopause | .2f }} of the
pressure-balance distance and the shock at {{ r.gate.shock | .2f }} of it, against
{{ r.gate.shock_spreiter | .2f }} from the gas-dynamic stand-off of Spreiter et al. (1966).

## Physics

**Magnetopause.** The wind stops where its ram pressure equals the magnetic pressure of the compressed
dipole field (Chapman–Ferraro; Spreiter et al. 1966). The field in front of a plane current sheet is
twice the dipole field B_d = M/R³, so

```math
\frac{(2 B_d)^2}{2} = n\, m_i v_{sw}^2 , \qquad B_d = \frac{M}{R_{mp}^3} .
```

The planner takes the distance R_mp and sets the moment M.

**Bow shock.** On the line from the planet towards the wind (the subsolar line) the interplanetary
field (IMF) is perpendicular to the flow. The Rankine–Hugoniot conditions of a perpendicular MHD shock
give the density compression X = n₂/n₁ as the positive root of

```math
(2-\gamma) X^2 + \left[\gamma(\beta+1) + (\gamma-1) M_A^2\right] X - (\gamma+1) M_A^2 = 0 , \qquad \gamma = 5/3 ,
```

and gas-dynamic models tie the thickness of the sheath to it (Spreiter et al. 1966):

```math
\frac{R_{bs} - R_{mp}}{R_{mp}} = \frac{1.1}{X} .
```

| quantity | value |
|---|---|
| plasma β (electrons and ions) | {{ c.gate.planner.beta | .2f }} |
| fast magnetosonic Mach number | {{ c.gate.planner.mach_fast | .2f }} |
| compression X of the MHD shock | {{ c.gate.planner.compression_mhd | .2f }} |
| sheath thickness 1.1 / X | {{ c.gate.planner.sheath_spreiter | .2f }} R_mp |

**Known results that are not analytic**: a shock needs a magnetosphere larger than the ion scales,
Dp = R_mp/d_i ≳ 1, and is Earth-like from Dp ≈ 20 (Omidi et al. 2004). The PIC model of Mercury of
Lavorenti et al. (2022) runs at Dp ≈ 8; this run uses Dp = {{ c.gate.params.dp }}.

## Experimental conditions

Normalised units: c = ε₀ = 1, lengths in d_e = c/ω_pe, times in 1/ω_pe, the wind density n₀ = 1.
The wind moves along +x at v_sw = {{ c.gate.planner.v_sw | .3f }} c. The IMF is along +z and the dipole
moment along −z, so the two fields are parallel on the day side. R_mp = {{ c.gate.planner.r_mp | .0f }} d_e;
the sphere has the radius {{ c.gate.planner.body_radius | .1f }} d_e.

- **External fields**: the IMF, the motional field E = −v × B and the dipole are external fields. The
  grid holds the difference to them, so the absorbing layers of the open faces relax the fields to the
  wind values. The outputs hold the total fields.
- **Open faces**: all six faces are open. The wind enters through the upstream face (electrons and ions
  in pairs at the same place) and through the four side faces; the downstream face only removes
  particles.
- **Sphere**: it removes the particles that enter it and keeps their charge, as an insulating body at
  floating potential. The box starts filled with the wind outside the sphere, and the dipole is there
  from t = 0.
- **Run length**: {{ c.gate.params.crossings }} crossing times L_x/v_sw, which is
  Ω_ci t = {{ c.gate.planner.omega_ci_t_end | .1f }} with the ion gyrofrequency of the IMF.

{{ table:conditions }}

<details><summary>Planner parameters</summary>

{{ table:params }}

</details>

## Animation

![bow shock](images/bow_shock/gate.gif)

The wind enters from the left; the sphere is blanked.

- **Left, ion density in the plane z = centre** (the plane perpendicular to the IMF and the dipole).
- **Middle, ion density in the plane y = centre** (the plane that holds the IMF and the dipole axis).
- **Right, B_z / B_imf in the plane z = centre**, logarithmic colour scale.

## Comparison with the pressure balance and the MHD shock

![subsolar line](images/bow_shock/gate-bow_shock.png)

Ion density in the two planes, averaged over the last half crossing time, and its profile on the
subsolar line with the magnetic pressure. The magnetopause is where the density falls back to the
wind density behind its maximum; the shock is where it first reaches twice the wind density.

![stand-off distances](images/bow_shock/standoff.png)

| quantity | PIC | reference | gate |
|---|---|---|---|
| magnetopause distance / R_mp | {{ r.gate.magnetopause | .3f }} | 1 (pressure balance) | within 25 % |
| shock distance / R_mp | {{ r.gate.shock | .3f }} | {{ r.gate.shock_spreiter | .3f }} (Spreiter et al. 1966) | upstream of the magnetopause |
| largest density on the subsolar line / n₀ | {{ r.gate.compression | .2f }} | {{ r.gate.compression_mhd | .2f }} (perpendicular MHD shock) | ≥ 2 |
| sheath thickness / magnetopause distance | {{ r.gate.sheath | .2f }} | {{ r.gate.sheath_spreiter | .2f }} (Spreiter et al. 1966) | reported |
| magnetopause by magnetic pressure / R_mp | {{ r.gate.magnetopause_magnetic | .3f }} | 1 | reported |
| upstream density and u_x, last crossing | within {{ r.gate.upstream_deviation | .2% }} of the wind | — | 5 % |
| Gauss's law, max\|∇·E − ρ\| | {{ r.gate.gauss_max | .0e }} | — | ≤ 1e-10 |

The gate criteria were fixed before the first run.

## Notes

- **The density maximum is not the mean sheath density.** It lies above the MHD compression.
  Supercritical perpendicular shocks are known to overshoot behind the ramp (Leroy et al. 1982), and the
  sheath here is about 2 d_i thick, too thin for a plateau behind an overshoot; this run does not
  separate the two.
- **The sheath is thicker than the gas-dynamic value** because the magnetopause lies inside R_mp while
  the shock is near the predicted distance.
- **The body charges positively**: the absorbed charge is {{ r.gate.body_charge | +.0f }} at the end, in
  units of e n₀ d_e³. The rate of absorption falls from {{ r.gate.absorbed_per_step.0 | .0f }} particles
  per step in the first crossing time to {{ r.gate.absorbed_per_step.2 | .0f }} in the third; it has not
  become constant by the end of the run.
- **Particle number**: the particle arrays, which hold twice the particles of the start, reach
  {{ r.gate.fill_max | .0% }} of their capacity on the fullest rank.
- **Run-to-run differences**: the GPU deposit adds currents atomically in an order that differs from
  run to run. Runs of this configuration with the same seed differ in the third digit of the
  distances above and by 2 % in the charge of the body (+730 to +749 in four runs).
- **Energy** is monitored only (explicit scheme).

## Reproduce

```bash
python -m pypic.launch -n 8 examples/bow_shock/run.py \
    --config examples/bow_shock/configs/gate.yaml --out output/wiki/bow_shock/gate
python examples/bow_shock/analyze.py output/wiki/bow_shock/gate --movie
```

## Data provenance

{{ table:provenance }}

## References

- J. R. Spreiter, A. L. Summers and A. Y. Alksne, "Hydromagnetic flow around the magnetosphere",
  Planet. Space Sci. 14, 223 (1966).
- N. Omidi, X. Blanco-Cano, C. T. Russell and H. Karimabadi, "Dipolar magnetospheres and their
  characterization as a function of magnetic moment", Adv. Space Res. 33, 1996 (2004).
- F. Lavorenti et al., "Electron dynamics in small magnetospheres: insights from global, fully kinetic
  plasma simulations of the planet Mercury", A&A 664, A133 (2022).
- M. M. Leroy, D. Winske, C. C. Goodrich, C. S. Wu and K. Papadopoulos, "The structure of perpendicular
  bow shocks", J. Geophys. Res. 87, 5081 (1982).
- E. R. Priest, *Solar Magnetohydrodynamics* (Reidel, 1982).
