# Test-particle orbits with exact relativistic solutions

Test electrons (PIC-D9: they feel the fields but deposit nothing) move in prescribed fields. They are tracked every step and compared with closed-form orbits (`theory.py`, checked against a 1e-11 ODE integration).

| config | fields | exact solution |
|---|---|---|
| `gyration.yaml` | uniform B_z, u⊥ = 0.1 … 10 | rotation at qB/(mγ), radius u⊥m/qB |
| `exb.yaml` | E_y, B_z, E/B = 0.5 | exact orbit in the drifting frame, guiding centre at E/B |
| `hyperbolic.yaml` | uniform E_x | u = qEt/m, x = (m/qE)(√(1 + (qEt/m)²) − 1) |
| `e_gt_b.yaml` | E_y > B_z | hyperbolic motion in the frame moving at B/E, boosted back |
| `plane_wave.yaml` | vacuum plane wave, a₀ = 2 (exact Yee eigenmode) | light-front invariant γ − u_x = 1; figure-8 orbit, drift a₀²/(4 + a₀²) = 0.5c |

```bash
python examples/particle_orbits/run.py --config examples/particle_orbits/configs/plane_wave.yaml --out output/particle_orbits/pw
python examples/particle_orbits/run.py --config examples/particle_orbits/configs/plane_wave.yaml --set dt_scale=0.5 --out output/particle_orbits/pw_half
python examples/particle_orbits/analyze.py output/particle_orbits/pw --refined output/particle_orbits/pw_half
```

**Gate**: maximum position error ≤ 2 % of the orbit size, and second-order convergence (the error ratio between Δt and Δt/2 lies in [3, 5]; for the plane wave Δx is refined too).

| case | error ratio |
|---|---|
| gyration | 4.000 |
| E×B | 4.000 |
| hyperbolic | 4.000 |
| E > B | 3.99 |
| plane wave, 128 → 256 cells | 3.76 |

The backward half kick of the explicit scheme is what makes these ratios 4 instead of 2.
