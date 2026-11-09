**English** · [日本語](Relativistic-Shock-ja) · source: [`examples/relativistic_shock`](https://github.com/ytworks/pypic/tree/main/examples/relativistic_shock)

# Relativistic pair shock

A cold electron–positron flow with Lorentz factor γ₀ = {{ c.pair_2d.params.gamma0 }} hits a reflecting wall
(the piston method). The incoming and the reflected plasma stream through each other, the
filamentation (Weibel) instability builds a magnetic field, and a collisionless shock forms and moves
back into the flow. Across a shock the speed and the density jump follow from conservation of particles,
momentum and energy, whatever happens inside the shock. The PIC shocks reproduce the relativistic
jump conditions within the ±10 % gate in 2-D and in 3-D, and the long run starts to build the
non-thermal tail of Fermi acceleration.

## Physics

The wall frame is the downstream rest frame. With the downstream plasma a relativistic Maxwell–Jüttner
(Synge) gas of temperature θ₂ = T₂/(mc²), energy and particle conservation give

```math
\langle\gamma\rangle(\theta_2) = \gamma_0, \qquad
\beta_{sh} = \frac{\theta_2}{\gamma_0\beta_0}, \qquad
\frac{n_2}{n_1} = 1 + \frac{\beta_0}{\beta_{sh}} ,
```

where ⟨γ⟩(θ) is the mean Lorentz factor of the gas, which depends on the number of degrees of freedom
of the particle motion. A quasi-2-D run has two, a 3-D run three:

| gas | θ₂ | β_sh | n₂/n₁ |
|---|---|---|---|
| 2-D Synge gas | {{ c.pair_2d.planner.jump_2d.theta2 | .2f }} | {{ c.pair_2d.planner.jump_2d.beta_sh | .3f }} | {{ c.pair_2d.planner.jump_2d.compression | .3f }} |
| 3-D Synge gas | {{ c.pair_2d.planner.jump_3d.theta2 | .2f }} | {{ c.pair_2d.planner.jump_3d.beta_sh | .3f }} | {{ c.pair_2d.planner.jump_3d.compression | .3f }} |
| constant adiabatic index Γ = 3/2 | — | {{ c.pair_2d.planner.jump_gamma_3_2.beta_sh | .3f }} | {{ c.pair_2d.planner.jump_gamma_3_2.compression | .3f }} |

**Known results that are not analytic**: particles that cross the shock many times gain energy (Fermi
acceleration) and form a power-law tail dN/dγ ∝ γ^−p, with p ≈ 2.4 in 2-D (Spitkovsky 2008; Sironi et
al. 2013) and p ≈ 2.2–2.3 from the theory of isotropic diffusion (Kirk et al. 2000; Achterberg et al.
2001); the maximum energy grows as γ_max ∝ t^{1/2}.

## Experimental conditions

Normalised units: c = ε₀ = 1, lengths in c/ω_pe, times in 1/ω_pe; the relativistic skin depth is
√γ₀ c/ω_pe = {{ c.pair_2d.planner.skin | .3f }}. The flow enters from the right along −x and the wall is at
x = 0. A compensated binomial current filter suppresses the numerical Cherenkov instability of the
relativistic flow. The useful time is bounded by the vacuum gap that opens at the far wall:
t < L_x/(1 + β_sh).

{{ table:conditions }}

<details><summary>Planner parameters</summary>

{{ table:params }}

</details>

## Animation

![relativistic shock](images/relativistic_shock/pair_2d.gif)

`pair_2d`, x–y plane (the wall is on the left).

- **Top left, B_z**: magnetic filaments form where the reflected particles meet the flow; behind the
  shock the field is turbulent.
- **Top right, electron density**: the compressed downstream region grows from the wall; its edge, the
  shock, moves to the right at β_sh ≈ 0.5.
- **Bottom left, electron phase space (x, u_x)**: the incoming beam at u_x = −{{ c.pair_2d.planner.u0 | .1f }}
  is thermalised at the shock; reflected particles run ahead of it at positive u_x.
- **Bottom right, n_e(x, t)** (x across, t up): the shock is the straight line from the origin. The
  distance it moves along x divided by the elapsed time is the shock speed β_sh (≈ 0.5), so the line's
  slope in this plot, Δt/Δx, is 1/β_sh.

## Comparison with the jump conditions

### Shock speed and compression of every run

![jump conditions of every run](images/relativistic_shock/jump_conditions.png)

Measured over theory for both quantities, with the gas of each run's dimensionality (2-D for the
quasi-2-D runs, 3-D for `pair_3d`). The shaded band is the ±10 % gate.

| run | gas | β_sh PIC / theory | error | n₂/n₁ PIC / theory | error | gate |
|---|---|---|---|---|---|---|
| `pair_2d` | {{ r.pair_2d.dims_eff }}-D | {{ r.pair_2d.beta_sh | .3f }} / {{ r.pair_2d.beta_sh_theory | .3f }} | {{ r.pair_2d.beta_sh_rel_error | +.1% }} | {{ r.pair_2d.compression | .2f }} / {{ r.pair_2d.compression_theory | .2f }} | {{ r.pair_2d.compression_rel_error | +.1% }} | ±10 % |
| `pair_2d_long` | {{ r.pair_2d_long.dims_eff }}-D | {{ r.pair_2d_long.beta_sh | .3f }} / {{ r.pair_2d_long.beta_sh_theory | .3f }} | {{ r.pair_2d_long.beta_sh_rel_error | +.1% }} | {{ r.pair_2d_long.compression | .2f }} / {{ r.pair_2d_long.compression_theory | .2f }} | {{ r.pair_2d_long.compression_rel_error | +.1% }} | reported |
| `pair_2d_long_run2` | {{ r.pair_2d_long_run2.dims_eff }}-D | {{ r.pair_2d_long_run2.beta_sh | .3f }} / {{ r.pair_2d_long_run2.beta_sh_theory | .3f }} | {{ r.pair_2d_long_run2.beta_sh_rel_error | +.1% }} | {{ r.pair_2d_long_run2.compression | .2f }} / {{ r.pair_2d_long_run2.compression_theory | .2f }} | {{ r.pair_2d_long_run2.compression_rel_error | +.1% }} | reported |
| `pair_3d` | {{ r.pair_3d.dims_eff }}-D | {{ r.pair_3d.beta_sh | .3f }} / {{ r.pair_3d.beta_sh_theory | .3f }} | {{ r.pair_3d.beta_sh_rel_error | +.1% }} | {{ r.pair_3d.compression | .2f }} / {{ r.pair_3d.compression_theory | .2f }} | {{ r.pair_3d.compression_rel_error | +.1% }} | reported |

The 3-D run jumps to the 3-D gas values (β_sh ≈ 1/3, n₂/n₁ ≈ 4) while the quasi-2-D runs follow the 2-D
gas (β_sh ≈ 1/2, n₂/n₁ ≈ 3): the shock knows how many directions its particles can move in.

### Shock structure and downstream spectrum

![pair_2d_long](images/relativistic_shock/pair_2d_long-shock.png)

`pair_2d_long` (t = 2000). **Left**: n_e(x, t) with the jump-condition shock trajectory (dashed).
**Middle**: the density profile at the measurement time and the n₂/n₁ level. **Right**: the downstream
spectrum, a relativistic Maxwellian with a non-thermal tail and the power-law fit.

![pair_3d](images/relativistic_shock/pair_3d-shock.png)

`pair_3d` (3-D, t = 600).

<details><summary>pair_2d (gate) and the second long realization</summary>

![pair_2d](images/relativistic_shock/pair_2d-shock.png)
![pair_2d_long_run2](images/relativistic_shock/pair_2d_long_run2-shock.png)

</details>

### Particle acceleration (reported)

| run | non-thermal fraction | tail index p (window fit) | γ_max ∝ t^α | ε_B at the end |
|---|---|---|---|---|
| `pair_2d_long` | {{ r.pair_2d_long.nonthermal_fraction | .1% }} | {{ r.pair_2d_long.tail_index_window_p | .2f }} ± {{ r.pair_2d_long.tail_index_window_err | .2f }} | α = {{ r.pair_2d_long.gamma_max_growth_exponent | .2f }} | {{ r.pair_2d_long.eps_B_end | .4f }} |
| `pair_2d_long_run2` | {{ r.pair_2d_long_run2.nonthermal_fraction | .1% }} | {{ r.pair_2d_long_run2.tail_index_window_p | .2f }} ± {{ r.pair_2d_long_run2.tail_index_window_err | .2f }} | α = {{ r.pair_2d_long_run2.gamma_max_growth_exponent | .2f }} | {{ r.pair_2d_long_run2.eps_B_end | .4f }} |
| literature (2-D, converged) | — | ≈ 2.4 | α = 0.5 | — |

The tail exists but is still young at t = 2000 ω_pe⁻¹ (≈ 520 relativistic plasma times √γ₀/ω_pe): its
index ({{ r.pair_2d_long.tail_index_window_p | .1f }}) is much steeper than the converged 2.4, which
Spitkovsky measured at about 8700 ω_pe⁻¹. A converged tail needs t ≳ 8000 ω_pe⁻¹ and about 10⁹ particles
in quasi-2-D (roadmap item R11).

## Notes

- **Two realizations of the long run**: same configuration and seed, run at two commits that do not
  change this set-up. The GPU deposit adds currents atomically in a non-deterministic order, and the
  turbulent flow amplifies the round-off differences, so the runs are independent realizations. Their
  spread is a measure of the statistical uncertainty.
- **Shock-speed fit window**: 0.3–0.7 of the valid time; later, the reflected-particle precursor reaches
  the receding plasma edge.

## Reproduce

```bash
python -m pypic.launch -n 4 examples/relativistic_shock/run.py \
    --config examples/relativistic_shock/configs/pair_2d.yaml --out output/wiki/relativistic_shock/pair_2d
python examples/relativistic_shock/analyze.py output/wiki/relativistic_shock/pair_2d --movie
```

## Data provenance

{{ table:provenance }}

## References

- A. Spitkovsky, "On the structure of relativistic collisionless shocks in electron-ion plasmas",
  ApJ 673, L39 (2008); "Particle acceleration in relativistic collisionless shocks: Fermi process at
  last?", ApJ 682, L5 (2008).
- L. Sironi, A. Spitkovsky and J. Arons, "The maximum energy of accelerated particles in relativistic
  collisionless shocks", ApJ 771, 54 (2013).
- J. G. Kirk, A. W. Guthmann, Y. A. Gallant and A. Achterberg, ApJ 542, 235 (2000).
- A. Achterberg, Y. A. Gallant, J. G. Kirk and A. W. Guthmann, MNRAS 328, 393 (2001).
- J. L. Synge, *The Relativistic Gas* (North-Holland, 1957).
