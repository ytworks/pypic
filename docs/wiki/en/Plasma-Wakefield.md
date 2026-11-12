**English** · [日本語](Plasma-Wakefield-ja) · source: [`examples/wakefield`](https://github.com/ytworks/pypic/tree/main/examples/wakefield)

# Plasma wakefield of a laser pulse and of an electron bunch

A short driver that crosses a cold plasma leaves a plasma wave behind it: every point oscillates at
ω_pe, and the pattern follows the driver. The driver is a laser pulse sent in by an antenna or a rigid
electron bunch, and the box follows it (moving window). Two weak drivers are compared with the linear
theory: the amplitude of the wake on the axis is {{ r.laser_linear.amplitude_ratio | .3f }} of the theory
behind the pulse and {{ r.beam_linear.amplitude_ratio | .3f }} behind the bunch.
A third, strong pulse (a₀ = {{ c.laser_bubble.params.a0 }}) empties a cavity of electrons behind it and traps electrons
in it: the cavity has the radius {{ r.laser_bubble.cavity_radius | .2f }} c/ω_pe (scaling of Lu et al. 2007:
{{ r.laser_bubble.theory.radius | .0f }}), and the fastest trapped electrons reach u_x = {{ r.laser_bubble.trapped.u_max | .0f }}.

## Physics

**Linear wake.** In a cold plasma with ions at rest the electron density answers the driver point by
point as a driven oscillator at ω_pe, so behind a driver of the speed v the wake is a wave with the
wavenumber k_w = ω_pe / v. Fields are in units of m c ω_pe / e.

**Laser pulse.** A pulse with the field envelope a₀ exp(−t²/τ²) exp(−r²/w²) pushes the electrons with
its ponderomotive force. Behind it

```math
E_x = k_w A\, e^{-2r^2/w^2} \cos(t - k_w x) , \qquad
E_r = \frac{4r}{w^2} A\, e^{-2r^2/w^2} \sin(t - k_w x) , \qquad
A = a_0^2\, \frac{\tau}{4} \sqrt{\frac{\pi}{2}}\, e^{-\tau^2/8} ,
```

which is largest for τ = 2 (the resonant Gaussian pulse of Schroeder et al. 2011). a₀ and w are the
local values of the diffracting beam, and a₀² is 1/η larger inside the plasma than in vacuum
(η = √(1 − ω_pe²/ω₀²): the pulse keeps its energy flux through the density ramp). The ratio
E_r / E_x = 4 r v / w² does not depend on a₀.

**Speed of the pulse.** On the Yee grid the pulse travels with the group velocity of the grid,

```math
v_g = \frac{\sin(k\,\Delta x)/\Delta x}{\sin(\omega_0 \Delta t)/\Delta t} , \qquad
\frac{\sin^2(\omega_0 \Delta t/2)}{\Delta t^2} = \frac{\sin^2(k\,\Delta x/2)}{\Delta x^2} + \frac{1}{4} ,
```

which is {{ c.laser_linear.planner.theory.v_g | .5f }} here, against √(1 − ω_pe²/ω₀²) = 0.99499 without
the grid. The comparison uses the grid value.

**Electron bunch.** A rigid Gaussian bunch of the peak density n_b, the lengths σ_x and σ_r and the
speed β leaves

```math
E_x = n_b \sqrt{2\pi}\, \sigma_x\, e^{-k_w^2 \sigma_x^2/2}\, R(r) \cos(k_w \xi) , \qquad
R(r) = k_w^2 \int_0^\infty r'\, e^{-r'^2/2\sigma_r^2}\, I_0(k_w r_<)\, K_0(k_w r_>)\, dr' ,
```

with ξ = x − βt and k_w = 1/β (Muggli, arXiv:2007.05226, Eqs. 1 and 3, for β = 1).

| quantity | laser pulse | electron bunch |
|---|---|---|
| speed of the driver | {{ c.laser_linear.planner.speed | .5f }} c | {{ c.beam_linear.planner.speed | .5f }} c |
| E_x amplitude on the axis | {{ c.laser_linear.planner.theory.Ex_axis_focus | .3e }} (at the focus) | {{ c.beam_linear.planner.theory.Ex_axis | .4f }} |
| radial structure | E_r / E_x = {{ c.laser_linear.planner.theory.Er_over_Ex_half_waist | .4f }} at r = w/2 | R(σ_r) / R(0) = {{ c.beam_linear.planner.theory.R_ratio | .4f }} |
| Rayleigh length of the pulse | {{ c.laser_linear.planner.theory.rayleigh | .1f }} c/ω_pe | — |

**Strong pulse (blowout).** A pulse with a₀ ≳ 2 and the matched waist w = 2√a₀ expels the electrons from
a nearly spherical cavity of the radius R = 2√a₀ behind it; the ions stay. On the axis the field of the
cavity is linear, E_x = (x − x_c)/2 through its centre x_c, so its largest value is √a₀. The wake follows
the pulse at the speed 1 − 3ω_pe²/(2ω₀²); an electron trapped at the back of the cavity stays in the
accelerating half over the dephasing length (2/3)(ω₀/ω_pe)² R and gains (2/3)(ω₀/ω_pe)² a₀ m c² there
(Lu et al. 2007). These are scalings, not a solution, so the strong pulse has no gate.

| quantity (a₀ = {{ c.laser_bubble.params.a0 }}, ω₀ = {{ c.laser_bubble.params.omega0 }} ω_pe) | scaling of Lu et al. 2007 |
|---|---|
| radius of the cavity | {{ r.laser_bubble.theory.radius | .0f }} c/ω_pe |
| slope of E_x on the axis | {{ r.laser_bubble.theory.slope | .1f }} |
| speed of the wake | {{ r.laser_bubble.theory.wake_speed | .3f }} c |
| dephasing length | {{ r.laser_bubble.theory.dephasing_length | .0f }} c/ω_pe |
| energy gain there | {{ r.laser_bubble.theory.energy_gain | .0f }} m c² |

## Experimental conditions

Normalised units: c = ε₀ = 1, lengths in c/ω_pe, times in 1/ω_pe, the plasma density n₀ = 1, fields in
m c ω_pe / e. The driver moves along +x.

- **Plasma**: electrons and ions (mass 1836) at {{ c.laser_linear.ppc }} particles per cell each, loaded
  at the same places, so the plasma has no field before the driver arrives. The particles of a cell sit
  on a 2 × 2 × 2 lattice, so the plasma also starts without the density noise of random places. The
  density rises over a ramp.
- **Laser pulse** (`laser_linear`): a₀ = {{ c.laser_linear.params.a0 }}, τ = {{ c.laser_linear.params.duration }},
  waist w = {{ c.laser_linear.params.waist }}, ω₀ = {{ c.laser_linear.params.omega0 }} ω_pe, linearly
  polarised along y. The antenna is a plane of currents in the vacuum before the ramp. The focus is at
  x = {{ c.laser_linear.params.focus | .2f }}; since the Rayleigh length is not long against the path, a₀²
  on the axis changes by up to 13 % along it and the theory is evaluated at the place of the measurement.
- **Strong pulse** (`laser_bubble`, a recorded run): a₀ = {{ c.laser_bubble.params.a0 }}, τ = {{ c.laser_bubble.params.duration }},
  w = {{ c.laser_bubble.params.waist }} = 2√a₀, ω₀ = {{ c.laser_bubble.params.omega0 }} ω_pe, focus at the top of the density ramp
  (x = {{ c.laser_bubble.params.focus }}). The box is {{ c.laser_bubble.L.1 | .0f }} c/ω_pe wide with transverse cells of
  {{ c.laser_bubble.dx.1 }} c/ω_pe, and the pulse crosses 172 c/ω_pe of plasma with
  {{ c.laser_bubble.window.window_shift_end }} shifts of the window.
- **Electron bunch** (`beam_linear`): n_b = {{ c.beam_linear.params.bunch_density }} n₀,
  σ_x = σ_r = {{ c.beam_linear.params.sigma_x }}, γ = {{ c.beam_linear.params.gamma }}. The bunch has the
  mass 10⁶, so it does not slow down, and it is loaded on a partner of the opposite charge at rest, so
  the run starts neutral and without a field. The partner stays behind as a charge that the plasma
  screens.
- **Moving window**: the box starts to move when the driver is {{ c.laser_linear.params.window_lead }} c/ω_pe
  behind its upper x face and then moves one cell at a time at the speed of light. The lower x face is
  open; the upper one is a conducting wall behind a vacuum margin of 1 c/ω_pe, and the plasma is loaded
  plane by plane as the box reaches it. Positions and outputs are laboratory values. The pulse run
  shifts the window {{ c.laser_linear.window.window_shift_end }} times, the bunch run
  {{ c.beam_linear.window.window_shift_end }} times.
- **Transverse faces**: periodic for the pulse, conducting walls for the bunch (between periodic faces
  the charge of the partner would add a uniform E_x of about 2 % of the wake).

{{ table:conditions }}

<details><summary>Planner parameters</summary>

{{ table:params }}

</details>

## Animation

![wake of a laser pulse](images/wakefield/laser_linear.gif)

The plane through the axis normal to the polarisation, at its laboratory place; the frame follows the
window.

- **Top, E_x**: the wake.
- **Middle, E_z**: the radial field of the wake in this plane.
- **Bottom, E_y**: the laser pulse.

![wake of an electron bunch](images/wakefield/beam_linear.gif)

- **Top, E_x** and **middle, E_z** as above.
- **Bottom, density of the bunch.**

![wake of a strong laser pulse](images/wakefield/laser_bubble.gif)

The strong pulse (`laser_bubble`) in the same plane.

- **Top, electron density**: the cavity behind the pulse and the cavities that follow it.
- **Second, E_x**: electrons are accelerated forward where E_x is negative (blue), in the back half of a
  cavity.
- **Third, E_y**: the laser pulse.
- **Bottom, phase space (x, u_x) of the electrons**: the trapped electrons, from the back of the first
  cavity forward.

## Comparison with the linear theory

![wake of the pulse](images/wakefield/laser_linear-wakefield.png)

![wake of the bunch](images/wakefield/beam_linear-wakefield.png)

E_x in the plane through the axis at the last output, and E_x on the axis with the fit and the amplitude
of the theory. The fit is E_x = c A_th(x) cos(k x + φ) + b with the amplitude A_th of the theory and free
c, k, φ, b, behind the driver and no closer than 1 c/ω_pe to the lower face of the box. The radial
structure is measured on the two lines at r = w/2 (pulse) or r = σ_r (bunch) with the same k and φ.

| quantity | PIC | linear theory | gate |
|---|---|---|---|
| pulse: amplitude of E_x on the axis, PIC / theory | {{ r.laser_linear.amplitude_ratio | .4f }} | 1 (amplitude {{ r.laser_linear.predicted_amplitude | .3e }} at the focus) | within 5 % |
| pulse: wavenumber, k v_g | {{ r.laser_linear.k_times_speed | .4f }} | 1 | within 2 % |
| pulse: E_r / E_x at r = w/2, PIC / theory | {{ r.laser_linear.Er_over_Ex_ratio | .4f }} | 1 (ratio {{ r.laser_linear.Er_over_Ex_theory | .4f }}) | within 10 % |
| bunch: amplitude of E_x on the axis, PIC / theory | {{ r.beam_linear.amplitude_ratio | .4f }} | 1 (amplitude {{ r.beam_linear.predicted_amplitude | .4f }}) | within 5 % |
| bunch: wavenumber, k β | {{ r.beam_linear.k_times_speed | .4f }} | 1 | within 2 % |
| bunch: amplitude at r = σ_r over the one on the axis | {{ r.beam_linear.R_ratio | .4f }} | {{ r.beam_linear.R_ratio_theory | .4f }} | within 10 % |
| Gauss's law, max\|∇·E − ρ\| | {{ r.laser_linear.gauss_max | .0e }} (pulse), {{ r.beam_linear.gauss_max | .0e }} (bunch) | — | ≤ 1e-10 |

The gate criteria were fixed before the first run.

## Strong pulse: the cavity and the trapped electrons

![cavity of the strong pulse](images/wakefield/laser_bubble-wakefield.png)

The last output (t = {{ r.laser_bubble.t | .1f }}). Top: the electron density in the plane through the axis, with
E_x on the axis (orange) and the line fitted at its zero inside the first cavity (dashed). Middle: the
phase space (x, u_x) of the electrons. Bottom: the spectrum of u_x of all electrons of the box; electrons
right of the dashed line (u_x > {{ r.laser_bubble.trapped.u_min | .0f }}) count as trapped.

| quantity | PIC | scaling of Lu et al. 2007 |
|---|---|---|
| radius of the cavity | {{ r.laser_bubble.cavity_radius | .2f }} | {{ r.laser_bubble.theory.radius | .0f }} |
| slope of E_x at its zero in the cavity | {{ r.laser_bubble.Ex_slope | .3f }} | {{ r.laser_bubble.theory.slope | .1f }} |
| a₀ of the pulse | {{ r.laser_bubble.a0_measured | .2f }} (largest of the outputs {{ r.laser_bubble.history.a0.9 | .2f }}, at t = {{ r.laser_bubble.history.t.9 | .0f }}) | {{ c.laser_bubble.params.a0 }} at the focus |
| electrons with u_x > {{ r.laser_bubble.trapped.u_min | .0f }} | {{ r.laser_bubble.trapped.number | .1f }} n₀ (c/ω_pe)³; for a wavelength of 0.8 µm {{ r.laser_bubble.trapped_electrons_800nm | .2e }} electrons, {{ r.laser_bubble.trapped_charge_pC_800nm | .0f }} pC | — |
| u_x of these electrons | mean {{ r.laser_bubble.trapped.u_mean | .0f }}, largest {{ r.laser_bubble.trapped.u_max | .0f }} | {{ r.laser_bubble.theory.energy_gain | .0f }} at the dephasing length |
| Gauss's law, max\|∇·E − ρ\| over the run | {{ r.laser_bubble.gauss_max | .1e }} | ≤ 1e-10 (contract of the engine) |

- **The pulse does not keep its shape.** Its a₀ rises to {{ r.laser_bubble.history.a0.9 | .2f }} at
  t = {{ r.laser_bubble.history.t.9 | .0f }} and falls back to {{ r.laser_bubble.a0_measured | .2f }}; the radius of the cavity
  is {{ r.laser_bubble.history.cavity_radius.2 | .2f }} at t = {{ r.laser_bubble.history.t.2 | .0f }} and {{ r.laser_bubble.cavity_radius | .2f }} at the last output. The
  slope of E_x stays between 0.41 and 0.46.
- **Electrons are trapped all through the run**, from t = {{ r.laser_bubble.history.t.3 | .0f }} on, so the spectrum is
  broad. The largest u_x is above the {{ r.laser_bubble.theory.energy_gain | .0f }} of the scaling from
  t = {{ r.laser_bubble.history.t.11 | .0f }} on ({{ r.laser_bubble.history.u_max.11 | .0f }}) and still grows when the run ends,
  after 172 of the {{ r.laser_bubble.theory.dephasing_length | .0f }} c/ω_pe of the dephasing length.
- **The extremes of E_x on the axis**, {{ r.laser_bubble.Ex_axis_max | .2f }} and {{ r.laser_bubble.Ex_axis_min | .2f }}, are a spike
  0.06 c/ω_pe wide at the trapped electrons, not the field of the cavity, which the slope measures.

## Notes

- **Without the driver** (test runs of this example, recorded in
  [`docs/plans/2026-11-11-wakefield.md`](https://github.com/ytworks/pypic/blob/main/docs/plans/2026-11-11-wakefield.md),
  as are the next two notes) both boxes stay exactly at rest: E_x = 0 on the axis and the density of the
  plateau is n₀ in every sample. The lattice load leaves no noise for the window to move.
- **Window against a box at rest.** The same runs in a box as long as the plasma give the same wake: on
  the same laboratory interval of the last output E_x on the axis of the pulse run differs by at most
  1.9e-9 of its maximum, and the one of the bunch run in no stored bit (float32 output). The box at rest
  holds 2.8 (pulse) and 2.3 (bunch) times the cells.
- **Lattice load.** With 8 particles at random places in the cell both amplitudes fail the gate in the
  box at rest (−8.3 % for the pulse, −6.9 % for the bunch). The loss grows with the age of the
  oscillation and becomes smaller with more particles per cell. Its cause is not established, and
  further behind the driver the wake decays with the lattice load too.
- **The ripple at the front of the E_x curve of the pulse run** is at the place of the laser pulse; the
  fit starts behind it.
- **∇·B in the run of the strong pulse** is 1.6e-11 at most while the pulse enters and 1.7e-12 to 2.7e-12
  afterwards, above the absolute 1e-12 of the engine's contract: one rounding unit of a field of 50 over
  the cell of 0.0196 c/ω_pe is 3.6e-13. The Gauss residual rises through the run (2.0e-11, 2.9e-11,
  3.6e-11 and 3.7e-11 in its four quarters); the cause of its part on the axis is not established.

## Reproduce

```bash
python -m pypic.launch -n 8 examples/wakefield/run.py \
    --config examples/wakefield/configs/laser_linear.yaml --out output/wiki/wakefield/laser_linear
python examples/wakefield/analyze.py output/wiki/wakefield/laser_linear --movie
```

`beam_linear` and `laser_bubble` run the same way. `--set drive=false` runs a box without the driver.

## Data provenance

{{ table:provenance }}

## References

- C. B. Schroeder, C. Benedetti, E. Esarey and W. P. Leemans, "Nonlinear pulse propagation and phase
  velocity of laser-driven plasma waves" (2011), [OSTI 1013067](https://www.osti.gov/servlets/purl/1013067).
- P. Muggli, "Beam-driven systems, plasma wakefield acceleration", CERN Accelerator School 2019,
  arXiv:2007.05226.
- W. Lu, M. Tzoufras, C. Joshi, F. S. Tsung, W. B. Mori, J. Vieira, R. A. Fonseca and L. O. Silva,
  "Generating multi-GeV electron bunches using single stage laser wakefield acceleration in a 3D
  nonlinear regime", Phys. Rev. ST Accel. Beams 10, 061301 (2007).
