**English** · [日本語](Two-Stream-Instability-ja) · source: [`examples/two_stream`](https://github.com/ytworks/pypic/tree/main/examples/two_stream)

# Two-stream and Buneman instabilities

Two interpenetrating particle beams are unstable. A small electrostatic wave along the drift bunches
the beams, the bunching charge reinforces the wave, and the wave grows exponentially until it traps
the beam particles. pypic runs this in 3-D for five configurations and measures the growth rate of the
fastest-growing mode, which linear theory predicts exactly. All five measured rates lie within the
±10 % gate: from {{ r.buneman.relative_error | +.1% }} (Buneman) to {{ r.ei_implicit.relative_error | +.1% }}
(implicit scheme).

## Physics

The planner solves the relativistic warm-fluid dispersion relation of all beams s (density n_s, drift
velocity v_s, Lorentz factor γ_s, thermal speed v_ts) for the fastest-growing wavenumber k_max and its
growth rate γ_max:

```math
1 = \sum_s \frac{\omega_{ps}^2 / \gamma_s^3}{(\omega - k v_s)^2 - 3 k^2 v_{ts}^2},
\qquad \omega = \omega_r + i\gamma .
```

The box is four k_max wavelengths long, so the fastest mode is the m = 4 Fourier mode of E_x. Its
amplitude must grow as e^{γ_max t} during the linear phase.

| config | set-up | theory k_max [ω_pe/c] | theory γ_max [ω_pe] | theory ω_r [ω_pe] |
|---|---|---|---|---|
| `ei` | two electron beams at u = ±{{ c.ei.params.u_drift }} c through ions at rest (m_i/m_e = {{ c.ei.params.mass_ratio | .0f }}) | {{ c.ei.planner.k_max | .3f }} | {{ c.ei.planner.gamma_max | .4f }} | 0 |
| `ep` | pair plasma: electrons and positrons, each in two beams at u = ±{{ c.ep.params.u_drift }} c | {{ c.ep.planner.k_max | .3f }} | {{ c.ep.planner.gamma_max | .4f }} | 0 |
| `eip` | electron beams ±{{ c.eip.params.u_drift }} c over ions (m_i/m_e = {{ c.eip.params.mass_ratio }}) with a positron fraction {{ c.eip.params.positron_fraction }} | {{ c.eip.planner.k_max | .3f }} | {{ c.eip.planner.gamma_max | .4f }} | {{ c.eip.planner.omega_r | .3f }} |
| `ei_implicit` | `ei` with the energy-conserving implicit scheme | {{ c.ei_implicit.planner.k_max | .3f }} | {{ c.ei_implicit.planner.gamma_max | .4f }} | 0 |
| `buneman` | one electron beam (u = {{ c.buneman.params.u_drift }} c) through ions at rest | {{ c.buneman.planner.k_max | .3f }} | {{ c.buneman.planner.gamma_max | .4f }} | {{ c.buneman.planner.omega_r | .4f }} |

- **Symmetric beams** (`ei`, `ep`): the unstable root is purely growing (ω_r = 0), so the wave
  pattern stands still while it grows.
- **Buneman**: the electrons stream through ions at rest. The ion inertia makes the instability slow;
  the cold-plasma limit is γ = (√3/2)(m_e/2m_i)^{1/3} ω_pe = 0.0561 ω_pe for hydrogen. The planner uses the
  warm root above. The wave has a real frequency, so it travels.
- **Implicit scheme** (`ei_implicit`): Crank–Nicolson fields with an implicit-midpoint push conserve the
  total energy exactly. The run uses Δx ≈ 13 λ_De and cΔt/Δx = 1, beyond the explicit stability limits.

## Experimental conditions

Normalised units: c = ε₀ = 1, lengths in c/ω_pe, times in 1/ω_pe, densities in n₀, momenta as u = γv/c,
temperatures θ = T/(m c²).

{{ table:conditions }}

<details><summary>Planner parameters (the physics inputs of each config)</summary>

{{ table:params }}

</details>

## Animations

### Electron–ion two-stream (`ei`)

![ei two-stream animation](images/two_stream/ei.gif)

- **Top left, electron phase space (x, u_x)**: the two beams at u_x = ±{{ c.ei.params.u_drift }} start to
  wiggle at the k_max wavelength. They roll up into vortices, where particles are trapped in the wave,
  and finally merge into a broad plateau.
- **Top right, E_x(x, t)**, averaged over y and z: four wavelengths grow out of the thermal noise. The
  pattern stands still because ω_r = 0.
- **Bottom left, E_x in an x–y plane**: the wavefronts are planes, because the wave is electrostatic
  with k along x.
- **Bottom right, energies**: the electric energy rises as a straight line on the log axis, at 2γ, and
  saturates once the beams are trapped. Some magnetic energy grows too: the 3-D beams are also
  unstable to transverse filamentation, which the planner estimates at γ ≈ {{ c.ei.planner.filamentation_gamma_est | .2f }} ω_pe.

### Buneman instability (`buneman`)

![Buneman animation](images/two_stream/buneman.gif)

The growth is about six times slower than in `ei`. The tilted stripes of E_x(x, t) show that the wave
travels (ω_r ≠ 0). At saturation the wave traps the electron beam in large phase-space vortices, and the
ion kinetic energy rises by an order of magnitude: the instability takes energy from the electron drift
and gives part of it to the ions.

## Comparison with linear theory

### Growth of the fastest mode

![growth of the k_max mode in every configuration](images/two_stream/growth.png)

|E_x| of the k_max Fourier mode in every configuration. The dashed line has the slope of the
linear-theory growth rate, and the shaded band is the window of the exponential fit. The measured
exponential phase follows the theory slope over two to three decades, up to saturation.

### Measured against predicted growth rate

![measured over predicted growth rate](images/two_stream/growth_ratio.png)

| config | γ theory [ω_pe] | γ PIC [ω_pe] | error | gate | max \|∇·E − ρ\| | max \|∇·B\| |
|---|---|---|---|---|---|---|
| `ei` | {{ r.ei.gamma_theory | .4f }} | {{ r.ei.gamma_fit_kmax | .4f }} | {{ r.ei.relative_error | +.1% }} | ±10 % | {{ r.ei.gauss_max | .1e }} | {{ r.ei.divB_max | .1e }} |
| `ep` | {{ r.ep.gamma_theory | .4f }} | {{ r.ep.gamma_fit_kmax | .4f }} | {{ r.ep.relative_error | +.1% }} | ±10 % | {{ r.ep.gauss_max | .1e }} | {{ r.ep.divB_max | .1e }} |
| `eip` | {{ r.eip.gamma_theory | .4f }} | {{ r.eip.gamma_fit_kmax | .4f }} | {{ r.eip.relative_error | +.1% }} | ±10 % | {{ r.eip.gauss_max | .1e }} | {{ r.eip.divB_max | .1e }} |
| `ei_implicit` | {{ r.ei_implicit.gamma_theory | .4f }} | {{ r.ei_implicit.gamma_fit_kmax | .4f }} | {{ r.ei_implicit.relative_error | +.1% }} | ±10 % | {{ r.ei_implicit.gauss_max | .1e }} | {{ r.ei_implicit.divB_max | .1e }} |
| `buneman` | {{ r.buneman.gamma_theory | .4f }} | {{ r.buneman.gamma_fit_kmax | .4f }} | {{ r.buneman.relative_error | +.1% }} | ±10 % | {{ r.buneman.gauss_max | .1e }} | {{ r.buneman.divB_max | .1e }} |

Gauss's law and ∇·B hold to round-off in every run (contracts: 1e-10 and 1e-12). The charge
conservation comes from the deposition itself, never from a Poisson correction.

### Energy conservation: explicit and implicit

![total energy error, explicit and implicit](images/two_stream/energy.png)

The implicit scheme keeps the total energy (particles + fields) to {{ r.ei_implicit.energy_error_max | .1e }}
of its initial value over the whole run, well inside its 1e-8 contract. The explicit scheme drifts by
{{ r.ei.energy_error_max | .1e }}. Explicit PIC does not conserve energy exactly, so this quantity is only
monitored for it.

<details><summary>Full diagnostic figure of each configuration</summary>

Each figure shows the energy history, the conservation history, all E_x modes with the fit, the
electron phase space at mid-run and at the end, and E_x(x, t).

![ei](images/two_stream/ei-summary.png)
![ep](images/two_stream/ep-summary.png)
![eip](images/two_stream/eip-summary.png)
![ei_implicit](images/two_stream/ei_implicit-summary.png)
![buneman](images/two_stream/buneman-summary.png)

</details>

## Notes

- The fit window is chosen automatically (`pypic.analysis.spectral.fit_growth`). In `eip` the positrons
  break the symmetry, so ω_r ≠ 0, and a two-wave model is fitted.
- `buneman` sets `remove_mean_current`: in a periodic box the net current of the drifting electrons would
  otherwise drive a uniform E_x that decelerates the whole beam. Removing the mean current models a
  return current.

## Reproduce

```bash
python -m pypic.launch -n 4 examples/two_stream/run.py --config examples/two_stream/configs/ei.yaml \
    --out output/wiki/two_stream/ei
python examples/two_stream/analyze.py output/wiki/two_stream/ei --movie
```

All the runs, charts and pages of this wiki: `python scripts/build_wiki.py run`, then `collect` and
`render` (see [`scripts/build_wiki.py`](https://github.com/ytworks/pypic/blob/main/scripts/build_wiki.py)).
`run` makes only the fresh gate-size runs; `collect` also needs the output of the recorded production runs
of the other pages, which is not in git (see [Home](Home#how-the-pages-are-made)). `render` alone works from
the checked-in `docs/wiki/data.json`.

## Data provenance

{{ table:provenance }}

## References

- O. Buneman, "Dissipation of currents in ionized media", Phys. Rev. 115, 503 (1959).
- C. K. Birdsall and A. B. Langdon, *Plasma Physics via Computer Simulation* (IOP, 1991).
