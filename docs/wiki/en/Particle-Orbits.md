**English** · [日本語](Particle-Orbits-ja) · source: [`examples/particle_orbits`](https://github.com/ytworks/pypic/tree/main/examples/particle_orbits)

# Test-particle orbits in prescribed fields

The particle pusher is the heart of a PIC code. Here test electrons move in fields with known
closed-form relativistic orbits: they feel the fields but deposit no current, so the fields stay exactly
as prescribed. Every orbit is tracked at every step and compared with the exact solution. The largest
position error of any case is {{ r.gyration.max_rel_position_error | .2% }} of the orbit size (gyration), and
halving Δt divides the error by {{ r.plane_wave.convergence_ratio | .2f }} (plane wave) to
{{ r.gyration.convergence_ratio | .2f }}: the scheme is second order.

## Physics

| config | fields | exact solution |
|---|---|---|
| `gyration` | uniform B_z, u⊥ = 0.1 … 10 | rotation at qB/(γm), radius u⊥m/(qB) |
| `exb` | E_y and B_z, E/B = 0.5 | gyration around a guiding centre that drifts at E×B/B², solved exactly in the drifting frame |
| `hyperbolic` | uniform E_x, starting at rest | u = qEt/m, x = (m/qE)(√(1 + (qEt/m)²) − 1) |
| `e_gt_b` | E_y > B_z (E/B > 1) | hyperbolic motion in the frame moving at B/E, boosted back to the lab |
| `plane_wave` | vacuum plane wave, a₀ = {{ c.plane_wave.params.a0 }}, the exact Yee eigenmode | light-front invariant γ − u_x = 1: figure-8 orbit, drift a₀²/(4 + a₀²) = {{ c.plane_wave.planner.drift_theory }} c |

The solutions are in `examples/particle_orbits/theory.py` and are themselves checked against an ODE
integration to 1e-11. The relativistic Boris push rotates the momentum in B exactly, so |u| is conserved
to round-off in a pure magnetic field.

## Experimental conditions

Normalised units: c = 1, lengths in c/ω, times in 1/ω, where ω is the characteristic frequency of each
case (cyclotron or wave frequency). Eight test electrons per case, one GPU; every case is run at Δt and
again at Δt/2 (the plane wave also halves Δx).

{{ table:conditions }}

<details><summary>Planner parameters</summary>

{{ table:params }}

</details>

## Animation

![test-particle orbits](images/particle_orbits/orbits.gif)

Each panel draws the PIC orbits (solid, coloured) as they are traced, on top of the exact orbits
(dashed black). The dashed curves are hidden under the PIC ones because the two agree to a fraction of a
percent.

- **Gyration**: three electrons with u⊥ = 2.7, 5.2 and 10; the larger momentum gives a larger radius.
- **E×B drift**: cycloids drifting along x.
- **Hyperbolic**: x(t) under a constant E_x, with the speed approaching c.
- **E > B**: the electric field wins and the electron runs away.
- **Plane wave**: in the frame that moves with the mean drift, the orbit is the classic figure 8.

## Comparison with the exact solutions

### Second-order convergence

![error against the time step](images/particle_orbits/convergence.png)

The maximum position error (relative to the orbit size) at Δt and at Δt/2 for every case. The dashed
guides have slope 2: a second-order scheme reduces the error 4-fold when Δt is halved.

| config | Δt | max error at Δt | max error at Δt/2 | ratio | gate |
|---|---|---|---|---|---|
| `gyration` | {{ r.gyration.dt | .4g }} | {{ r.gyration.max_rel_position_error | .2e }} | {{ r.gyration.refined_max_error | .2e }} | {{ r.gyration.convergence_ratio | .3f }} | error ≤ 2 %, ratio 3–5 |
| `exb` | {{ r.exb.dt | .4g }} | {{ r.exb.max_rel_position_error | .2e }} | {{ r.exb.refined_max_error | .2e }} | {{ r.exb.convergence_ratio | .3f }} | error ≤ 2 %, ratio 3–5 |
| `hyperbolic` | {{ r.hyperbolic.dt | .4g }} | {{ r.hyperbolic.max_rel_position_error | .2e }} | {{ r.hyperbolic.refined_max_error | .2e }} | {{ r.hyperbolic.convergence_ratio | .3f }} | error ≤ 2 %, ratio 3–5 |
| `e_gt_b` | {{ r.e_gt_b.dt | .4g }} | {{ r.e_gt_b.max_rel_position_error | .2e }} | {{ r.e_gt_b.refined_max_error | .2e }} | {{ r.e_gt_b.convergence_ratio | .3f }} | error ≤ 2 %, ratio 3–5 |
| `plane_wave` | {{ r.plane_wave.dt | .4g }} | {{ r.plane_wave.max_rel_position_error | .2e }} | {{ r.plane_wave.refined_max_error | .2e }} | {{ r.plane_wave.convergence_ratio | .3f }} | error ≤ 2 %, ratio 3–5 |

Plane-wave drift velocity: measured {{ r.plane_wave.drift_measured | .4f }} c, theory
{{ r.plane_wave.drift_theory_node | .4f }} c.

### Orbits and errors of every case

Left: x–y orbits, PIC (colour) over exact (dashed). Right: |x_PIC − x_exact| against time.

![gyration](images/particle_orbits/gyration-orbits.png)
![E×B](images/particle_orbits/exb-orbits.png)
![hyperbolic](images/particle_orbits/hyperbolic-orbits.png)
![E > B](images/particle_orbits/e_gt_b-orbits.png)
![plane wave](images/particle_orbits/plane_wave-orbits.png)

## Notes

- **Why the ratio is 4, not 2**: the explicit scheme starts with a backward half kick of the momentum,
  so positions and momenta are staggered by Δt/2 from the start, as the leapfrog scheme assumes. The
  global error is then second order.
- **Plane wave**: the wave is an eigenmode of the Yee grid, so this run refines Δx together with Δt. Its
  ratio ({{ r.plane_wave.convergence_ratio | .2f }}) measures the combined space–time convergence and is a
  little below 4.

## Reproduce

```bash
python -m pypic.launch -n 1 examples/particle_orbits/run.py --config examples/particle_orbits/configs/plane_wave.yaml \
    --out output/wiki/particle_orbits/plane_wave
python -m pypic.launch -n 1 examples/particle_orbits/run.py --config examples/particle_orbits/configs/plane_wave.yaml \
    --set dt_scale=0.5 --out output/wiki/particle_orbits/plane_wave_half
python examples/particle_orbits/analyze.py output/wiki/particle_orbits/plane_wave \
    --refined output/wiki/particle_orbits/plane_wave_half
```

## Data provenance

{{ table:provenance }}

## References

- J. P. Boris, "Relativistic plasma simulation — optimization of a hybrid code", Proc. 4th Conf. Num.
  Sim. Plasmas, 3 (1970).
- L. D. Landau and E. M. Lifshitz, *The Classical Theory of Fields* (motion of a charge in constant
  uniform fields and in a plane wave).
