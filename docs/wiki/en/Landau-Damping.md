**English** · [日本語](Landau-Damping-ja) · source: [`examples/landau_damping`](https://github.com/ytworks/pypic/tree/main/examples/landau_damping)

# Langmuir oscillations: Landau damping, Dawson oscillation, wave breaking

Electron plasma (Langmuir) oscillations are the simplest collective motion of a plasma, and several of
their properties are known exactly. This example checks five of them: the collisionless (Landau)
damping of a small wave in a warm plasma, the nonlinear damping and trapping of a large one, the exact
cold oscillation of Dawson, the moment a cold wave breaks, and the period of a relativistic
oscillation. The linear Landau damping rate comes out {{ r.landau.gamma_rel_error | +.1% }} from the kinetic
root, and the three exact cold solutions are reproduced to between
{{ r.relativistic.freq_rel_error_median | .2% }} (relativistic period) and {{ r.wavebreaking.rel_error | abs | .1% }}
(wave-breaking time).

## Physics

| config | set-up | exact or theoretical result |
|---|---|---|
| `landau` | warm electrons (θ_e = {{ c.landau.params.theta_e }}), n = n₀(1 + α cos kx), α = {{ c.landau.params.alpha }}, kλ_De = {{ c.landau.params.k_lambda_D }} | kinetic root of the Maxwellian dielectric function: ω = {{ r.landau.omega_theory | .4f }} ω_pe, γ = {{ r.landau.gamma_theory | .4f }} ω_pe |
| `nonlinear` | the same with α = {{ c.nonlinear.params.alpha }} | strong-damping benchmark (Cheng & Knorr 1976): the envelope first decays at γ₁ ≈ −0.29, then grows again by particle trapping at γ₂ ≈ +0.08 |
| `dawson` | cold electrons, velocity v = v₀ sin kx₀, kv₀ = 0.5 ω_pe | every electron sheet oscillates as x = x₀ + (v₀/ω_pe) sin(kx₀) sin(ω_pe t), exactly |
| `wavebreaking` | cold, kv₀ = 1.2 ω_pe | the sheets first cross (the wave breaks) at ω_pe t = arcsin(ω_pe/kv₀) = {{ r.wavebreaking.t_cross_theory | .3f }} |
| `relativistic` | cold, β_max = {{ r.relativistic.beta_max | .2f }} | exact relativistic sheet period from the invariant γ + ω_pe²ξ²/(2c²) = const; T ≈ (2π/ω_pe)(1 + 3β²/16) |

**Landau damping.** In a warm plasma the wave exchanges energy with the electrons moving at its phase
velocity. With a Maxwellian there are more slower electrons than faster ones, so the wave loses energy
and decays even without collisions. The complex frequency is the root of

```math
\varepsilon(k,\omega) = 1 - \frac{1}{2k^2\lambda_{De}^2}\, Z'\!\left(\frac{\omega}{\sqrt{2}\,k v_{te}}\right) = 0 ,
```

where Z is the plasma dispersion function (Fried & Conte).

**Dawson's cold oscillation.** In a cold plasma every electron sheet is a harmonic oscillator at ω_pe as
long as the sheets do not cross, whatever the amplitude. If kv₀ > ω_pe, neighbouring sheets overtake
each other: the density becomes infinite and the wave breaks.

## Experimental conditions

Normalised units: c = ε₀ = 1, lengths in c/ω_pe, times in 1/ω_pe, θ = T/(m c²). The box is one
wavelength long; ions are a fixed neutralising background (m_i/m_e = 1836).

{{ table:conditions }}

<details><summary>Planner parameters</summary>

{{ table:params }}

</details>

**Initial state.** The density perturbation is made by moving the loaded electrons; the move is
deposited as a current with the charge-conserving scheme, so the initial E_x satisfies Gauss's law to
round-off. For the explicit scheme the momenta are then pushed back by Δt/2. The cold cases start at
equilibrium with a velocity perturbation and no field.

## Animations

### Nonlinear Landau damping and trapping (`nonlinear`)

![nonlinear Landau damping](images/landau_damping/nonlinear.gif)

- **Top left, electron phase space (x, u_x)**: the large wave first damps; electrons near the phase
  velocity are caught by the wave and form vortices (holes) that rotate at the bounce frequency.
- **Top right, E_x(x, t)**: the wave decays in the first ω_pe t ≈ 15, then partially recovers as the
  trapped electrons give energy back.
- **Bottom left, energy spectrum** and **bottom right, energies**: electric energy flows into the
  electrons and partly back.

### Wave breaking (`wavebreaking`)

![wave breaking](images/landau_damping/wavebreaking.gif)

Drawn from 1,500 tracked electrons (sheets), coloured by their initial position x₀.

- **Left, phase space (x, u_x)**: the cold electrons lie on one curve. With kv₀ = 1.2 ω_pe the curve
  steepens until it becomes vertical and then folds over: faster sheets overtake slower ones.
- **Middle, Lagrangian map x₀ → x**: while the sheets keep their order the map is monotonic. It stops
  being monotonic at the moment the wave breaks; theory puts this first crossing at
  ω_pe t = {{ r.wavebreaking.t_cross_theory | .3f }}.
- **Right, electron density**: where the sheets cross the density spikes (it would be infinite for
  ideal sheets).

The sheet orbits are no longer the Dawson solution after the crossing: the oscillation is destroyed
and the electrons mix in phase space.

## Comparison with theory

### Errors against the tolerances

![error of every comparison](images/landau_damping/errors.png)

Each marker is |PIC − theory|/theory of one quantity; the vertical tick is its gate tolerance. Hollow
markers are the nonlinear benchmark, which is reported without a tolerance.

| config | quantity | PIC | theory | error | gate |
|---|---|---|---|---|---|
| `landau` | damping rate γ [ω_pe] | {{ r.landau.gamma_fit | .4f }} | {{ r.landau.gamma_theory | .4f }} | {{ r.landau.gamma_rel_error | +.1% }} | ±10 % |
| `landau` | frequency ω [ω_pe] | {{ r.landau.omega_fit | .4f }} | {{ r.landau.omega_theory | .4f }} | {{ r.landau.omega_rel_error | +.2% }} | ±2 % |
| `nonlinear` | early damping γ₁ | {{ r.nonlinear.gamma_1_fit | .3f }} | −0.29 (Cheng & Knorr) | — | reported |
| `nonlinear` | trapping growth γ₂ | {{ r.nonlinear.gamma_2_fit | .3f }} | +0.08 (Cheng & Knorr) | — | reported |
| `dawson` | sheet frequency (median over {{ r.dawson.sheets }} sheets) | — | ω_pe | {{ r.dawson.freq_rel_error_median | .2% }} | 1 % |
| `dawson` | sheet amplitude (median) | — | v₀ \|sin kx₀\|/ω_pe | {{ r.dawson.amp_rel_error_median | .2% }} | 2 % |
| `wavebreaking` | first sheet crossing [1/ω_pe] | {{ r.wavebreaking.t_cross | .3f }} | {{ r.wavebreaking.t_cross_theory | .3f }} | {{ r.wavebreaking.rel_error | +.1% }} | ±5 % |
| `relativistic` | sheet frequency (median over {{ r.relativistic.sheets }} sheets) | — | exact period | {{ r.relativistic.freq_rel_error_median | .3% }} | 0.5 % |

In the relativistic case the frequency shift itself is {{ r.relativistic.nonrelativistic_shift | .1% }}
(the median of 1 − ω/ω_pe over the sheets), so the measured period resolves the relativistic correction
with a large margin. Gauss's law holds to round-off in every run (max |∇·E − ρ|: `landau` {{ r.landau.gauss_max | .1e }}, `nonlinear` {{ r.nonlinear.gauss_max | .1e }}, `dawson` {{ r.dawson.gauss_max | .1e }}, `wavebreaking` {{ r.wavebreaking.gauss_max | .1e }}, `relativistic` {{ r.relativistic.gauss_max | .1e }}).

### Linear Landau damping

![damped mode, linear](images/landau_damping/landau-mode.png)

|E_x(k)| of the perturbed mode on a log axis. The oscillation decays along the dashed line, whose slope
is the kinetic damping rate γ = {{ r.landau.gamma_theory | .4f }} ω_pe, until the thermal noise floor is
reached.

### Nonlinear damping and trapping

![damped mode, nonlinear](images/landau_damping/nonlinear-mode.png)

The large-amplitude wave first decays faster than linear theory predicts (dashed), at a rate close to
Cheng & Knorr's γ₁ (dotted). Trapped electrons then make the envelope grow again at a rate close to
γ₂ (dash-dot).

### Cold sheets: Dawson oscillation, wave breaking, relativistic oscillation

![Dawson sheets](images/landau_damping/dawson-sheets.png)
![wave-breaking sheets](images/landau_damping/wavebreaking-sheets.png)
![relativistic sheets](images/landau_damping/relativistic-sheets.png)

Positions x(t) of tracked electrons (sheets). In `dawson` every sheet is a pure sinusoid at ω_pe; in
`wavebreaking` the trajectories cross shortly after ω_pe t ≈ 1; in `relativistic` the period is longer
than 2π/ω_pe, as the exact solution predicts.

## Reproduce

```bash
python -m pypic.launch -n 4 examples/landau_damping/run.py --config examples/landau_damping/configs/landau.yaml \
    --out output/wiki/landau_damping/landau
python examples/landau_damping/analyze.py output/wiki/landau_damping/landau --movie
```

## Data provenance

{{ table:provenance }}

## References

- L. D. Landau, "On the vibrations of the electronic plasma", J. Phys. USSR 10, 25 (1946).
- J. M. Dawson, "Nonlinear electron oscillations in a cold plasma", Phys. Rev. 113, 383 (1959).
- C. Z. Cheng and G. Knorr, "The integration of the Vlasov equation in configuration space",
  J. Comput. Phys. 22, 330 (1976).
- B. D. Fried and S. D. Conte, *The Plasma Dispersion Function* (Academic Press, 1961).
