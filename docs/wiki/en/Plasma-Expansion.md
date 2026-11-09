**English** · [日本語](Plasma-Expansion-ja) · source: [`examples/plasma_expansion`](https://github.com/ytworks/pypic/tree/main/examples/plasma_expansion)

# Plasma expansion into vacuum (Mora 2003)

A slab of hot electrons and cold ions is left free to expand into vacuum. The electrons run ahead and
set up a charge-separation sheath at the ion front, and its field accelerates the ions. This is the
mechanism of target-normal sheath acceleration (TNSA) of ions by intense lasers. Mora (2003) solved
the expansion exactly for isothermal electrons. With the proton mass the PIC ions follow the predicted
self-similar velocity profile to {{ r.isothermal_hydrogen.velocity_profile_median_rel_error | .1% }}, and
their maximum energy lies {{ r.isothermal_hydrogen.E_max_median_rel_error | .0% }} below the isothermal
law: the kinetic electrons of the PIC are cooler than the model's ideal ones.

## Physics

With τ = ω_pi t/√(2e), ion charge Z, electron temperature T_e and initial Debye length λ_D0, the
isothermal self-similar solution gives:

```math
E_{max} = 2 Z T_e \left[\ln\left(\tau + \sqrt{\tau^2 + 1}\right)\right]^2 ,
\qquad
x_f = 2\sqrt{2e}\,\lambda_{D0}\left[\tau \ln\left(\tau + \sqrt{\tau^2+1}\right) - \sqrt{\tau^2+1} + 1\right],
```

```math
\frac{dN}{dE} = \frac{n_{i0}\, c_s t}{\sqrt{2 E Z T_e}} \exp\left(-\sqrt{\frac{2E}{Z T_e}}\right),
\qquad v = c_s + \frac{x}{t}\ \ \text{(quasi-neutral rarefaction)} ,
```

where c_s = √(ZT_e/m_i) is the ion sound speed. The front field decays as E_f = 2E₀/√(2e + ω_pi²t²).

| config | m_i/m_e | slab | run to | use |
|---|---|---|---|---|
| `isothermal_hydrogen` | {{ c.isothermal_hydrogen.params.mass_ratio | .0f }} | {{ c.isothermal_hydrogen.params.slab_lambda_D }} λ_D | τ = {{ c.isothermal_hydrogen.params.tau_end }} | gate configuration |
| `isothermal` | {{ c.isothermal.params.mass_ratio }} | {{ c.isothermal.params.slab_lambda_D }} λ_D | τ = {{ c.isothermal.params.tau_end }} | fast; here the front speed approaches the electron thermal speed |
| `thin_foil` | {{ c.thin_foil.params.mass_ratio }} | {{ c.thin_foil.params.slab_lambda_D }} λ_D | τ = {{ c.thin_foil.params.tau_end }} | the electrons cool as the foil expands: the isothermal law is an upper bound |
| `isothermal_hydrogen_long` | {{ c.isothermal_hydrogen_long.params.mass_ratio | .0f }} | {{ c.isothermal_hydrogen_long.params.slab_lambda_D }} λ_D | τ = {{ c.isothermal_hydrogen_long.params.tau_end }} | recorded production run |

## Experimental conditions

Normalised units: c = ε₀ = 1, lengths in c/ω_pe, times in 1/ω_pe, energies in m_e c². T_e =
{{ c.isothermal_hydrogen.params.Te }} m_e c² (5.11 keV), so λ_D0 = {{ c.isothermal_hydrogen.planner.lambda_D0 | .2f }} c/ω_pe. The
grid is quasi-1-D (4 × 4 transverse cells) with Δx = λ_D0; electrons and ions share one step density
profile, so the plasma is exactly neutral at the start. Both x walls reflect particles.

{{ table:conditions }}

<details><summary>Planner parameters</summary>

{{ table:params }}

</details>

## Animation

![plasma expansion](images/plasma_expansion/isothermal_hydrogen_long.gif)

`isothermal_hydrogen_long`, the recorded production run to τ = 30.

- **Top left, ion phase space (x, u_x)**: a rarefaction wave eats into the slab from its edge, and the
  ions beyond the edge are accelerated into a straight ramp, the self-similar velocity profile v = c_s + x/t.
- **Top right, ion and electron densities**: the density falls off smoothly into the vacuum; electrons
  and ions overlap closely (quasi-neutral) except at the very front.
- **Bottom left, E_x**: the sheath field peaks at the ion front and decays as the front moves out.
- **Bottom right, maximum ion energy**: E_max(τ) measured in the PIC (dots), drawn as time advances, over
  Mora's isothermal law (dashed).

## Comparison with Mora's solution

![expansion, proton mass](images/plasma_expansion/isothermal_hydrogen-expansion.png)

`isothermal_hydrogen` (gate). **Left**: E_max(τ). **Middle**: ion-front position. **Right**: ion energy
spectrum at the end, against Mora's dN/dE.

![expansion, long run](images/plasma_expansion/isothermal_hydrogen_long-expansion.png)

`isothermal_hydrogen_long`, to τ = 30. The spectrum follows the self-similar dN/dE over three decades
of energy.

| config | E_max at τ_end, PIC / Mora [m_e c²] | E_max median error | front median error | velocity profile median error | T_e end / start | gate |
|---|---|---|---|---|---|---|
| `isothermal_hydrogen` | {{ r.isothermal_hydrogen.E_max_final | .4f }} / {{ r.isothermal_hydrogen.E_max_theory_final | .4f }} | {{ r.isothermal_hydrogen.E_max_median_rel_error | .1% }} | {{ r.isothermal_hydrogen.front_median_rel_error | .1% }} | {{ r.isothermal_hydrogen.velocity_profile_median_rel_error | .1% }} | {{ r.isothermal_hydrogen.electron_temperature_final_over_initial | .2f }} | velocity ≤ 10 %, E_max ≤ 20 % |
| `isothermal` | {{ r.isothermal.E_max_final | .4f }} / {{ r.isothermal.E_max_theory_final | .4f }} | {{ r.isothermal.E_max_median_rel_error | .1% }} | {{ r.isothermal.front_median_rel_error | .1% }} | {{ r.isothermal.velocity_profile_median_rel_error | .1% }} | {{ r.isothermal.electron_temperature_final_over_initial | .2f }} | reported |
| `isothermal_hydrogen_long` | {{ r.isothermal_hydrogen_long.E_max_final | .4f }} / {{ r.isothermal_hydrogen_long.E_max_theory_final | .4f }} | {{ r.isothermal_hydrogen_long.E_max_median_rel_error | .1% }} | {{ r.isothermal_hydrogen_long.front_median_rel_error | .1% }} | {{ r.isothermal_hydrogen_long.velocity_profile_median_rel_error | .1% }} | {{ r.isothermal_hydrogen_long.electron_temperature_final_over_initial | .2f }} | reported |
| `thin_foil` | {{ r.thin_foil.E_max_final | .4f }} / {{ r.thin_foil.E_max_theory_final | .4f }} | {{ r.thin_foil.E_max_median_rel_error | .1% }} | {{ r.thin_foil.front_median_rel_error | .1% }} | {{ r.thin_foil.velocity_profile_median_rel_error | .1% }} | {{ r.thin_foil.electron_temperature_final_over_initial | .2f }} | E_max ≤ 1.05 × isothermal law |

The errors are medians of |PIC − Mora|/Mora over 2 ≤ τ ≤ τ_end. In every thick-slab run the PIC values
lie below the isothermal model. Gauss's law holds to {{ r.isothermal_hydrogen.gauss_max | .1e }}
(`isothermal_hydrogen`).

![thin foil](images/plasma_expansion/thin_foil-expansion.png)

`thin_foil`: with a slab of only {{ c.thin_foil.params.slab_lambda_D }} λ_D the electrons lose a large part
of their energy to the ions (T_e falls to {{ r.thin_foil.electron_temperature_final_over_initial | .0% }} of its
initial value), and E_max stays well below the isothermal law, as expected for an adiabatic expansion
(Mora 2005).

<details><summary>m_i/m_e = 100</summary>

![expansion, m_i/m_e = 100](images/plasma_expansion/isothermal-expansion.png)

</details>

## Why the PIC lies below the isothermal model

- The model assumes isothermal (Boltzmann) electrons backed by an infinite heat bath. Collisionless
  PIC electrons have no such bath: the expanding region is cooler than the slab.
- All the deviations (E_max, front position, velocity profile) are consistent with an effective T_e
  about 15 % lower in the expanding region.
- The offset does not change with Δx = λ_D/2, twice the particles or a 2.5× longer vacuum (sensitivity
  runs described in the [example README](https://github.com/ytworks/pypic/tree/main/examples/plasma_expansion)),
  so it is not numerical. It stays at the same level from τ = 10 to τ = 30, and it grows at m_i/m_e = 100, where the
  front speed approaches the electron thermal speed and the model's assumption v_ion ≪ v_te fails.

## Reproduce

```bash
python -m pypic.launch -n 4 examples/plasma_expansion/run.py \
    --config examples/plasma_expansion/configs/isothermal_hydrogen.yaml --out output/wiki/plasma_expansion/isothermal_hydrogen
python examples/plasma_expansion/analyze.py output/wiki/plasma_expansion/isothermal_hydrogen --movie
```

## Data provenance

{{ table:provenance }}

## References

- P. Mora, "Plasma expansion into a vacuum", Phys. Rev. Lett. 90, 185002 (2003).
- P. Mora, "Thin-foil expansion into a vacuum", Phys. Rev. E 72, 056401 (2005).
