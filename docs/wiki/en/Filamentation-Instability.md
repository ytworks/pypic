**English** · [日本語](Filamentation-Instability-ja) · source: [`examples/filamentation`](https://github.com/ytworks/pypic/tree/main/examples/filamentation)

# Current filamentation (Weibel) instability

Two plasma beams streaming through each other along ±x are unstable to perturbations across the flow.
A small transverse ripple of the current makes a magnetic field that pinches like currents together and
pushes opposite currents apart. The ripple grows into current filaments, and the magnetic field they
carry is what later mediates collisionless shocks. Over the modes m = 2…10, the growth rates measured in
the pair plasma lie a median {{ r.pair.median_rel_error_m2_10 | .1% }} from the cold relativistic
theory (gate 10 %).

## Physics

For cold beams of drift β (Lorentz factor γ) along x and a wave vector k across the flow, the
relativistic fluid theory gives the growth rate Γ of a purely growing mode (c = 1):

```math
\Gamma^4 + \Gamma^2\left(k^2 + \sum_s \frac{\omega_{ps}^2}{\gamma_s^3}\right) - k^2 \sum_s \frac{\omega_{ps}^2\beta_s^2}{\gamma_s} = 0 .
```

The response along the drift carries the longitudinal mass γ³m. Γ increases with k and tends to
βω_p/√γ for k ≫ ω_p/c. The planner evaluates Γ for every Fourier mode that fits in the box.

| config | beams | Γ of the fastest mode in the box | Γ for k → ∞ | trapping estimate of B_sat |
|---|---|---|---|---|
| `pair` | electrons and positrons, each in two beams at u = ±{{ c.pair.planner.u0 }} | {{ c.pair.planner.gamma_max | .3f }} ω_pe | {{ c.pair.planner.gamma_limit | .3f }} ω_pe | {{ c.pair.planner.B_sat_trapping | .3f }} |
| `electron_ion` | electron beams at u = ±{{ c.electron_ion.planner.u0 }} over ions at rest | {{ c.electron_ion.planner.gamma_max | .3f }} ω_pe | {{ c.electron_ion.planner.gamma_limit | .3f }} ω_pe | {{ c.electron_ion.planner.B_sat_trapping | .3f }} |

The magnetic field saturates roughly when the filaments' field traps the beam particles, i.e. when the
bounce frequency in the field reaches the growth rate: B_sat ≈ γ₀mΓ²/(qβ₀k) (Davidson et al. 1972).

## Experimental conditions

Normalised units: c = ε₀ = 1, lengths in c/ω_pe, times in 1/ω_pe. The box is only 4 cells long along the
flow, which excludes the two-stream modes (k ∥ v), and resolves k in the (y, z) plane with 10 cells per
skin depth c/ω_p.

{{ table:conditions }}

<details><summary>Planner parameters</summary>

{{ table:params }}

</details>

## Animation

![filamentation in the pair plasma](images/filamentation/pair.gif)

`pair`, in the (y, z) plane across the flow.

- **Top left, B_z** and **top right, J_x**: from noise, a pattern of current filaments of alternating
  sign appears, each wrapped in its own magnetic field. Small filaments then merge into fewer, larger
  ones (the coalescence that follows saturation).
- **Bottom left, phase space (y, u_x)**: the two beams at u_x = ±{{ c.pair.planner.u0 }} are deflected
  as the magnetic field grows; late in the run the particles oscillate across the filaments, trapped in
  their magnetic field.
- **Bottom right, energies**: the magnetic energy grows exponentially, at twice the growth rate, and
  saturates at about {{ r.pair.eps_B_peak | .0% }} of the initial kinetic energy.

## Comparison with linear theory

![growth rate against k, pair plasma](images/filamentation/pair-filamentation.png)

`pair`. **Left**: the growth rate of every mode, measured from |B(k_⊥)| over one common linear window
(dots), on the cold-theory curve Γ(k) (dashed). **Middle**: |B(k)| of three modes against time: the
faster modes saturate first. **Right**: B_z at the magnetic-energy peak.

![growth rate against k, electron-ion](images/filamentation/electron_ion-filamentation.png)

`electron_ion` (reported, not a gate configuration), the same analysis.

| config | median \|Γ_PIC − Γ_theory\|/Γ_theory, m = 2…10 | maximum | gate | linear window [1/ω_pe] | ε_B at peak | peak time [1/ω_pe] |
|---|---|---|---|---|---|---|
| `pair` | {{ r.pair.median_rel_error_m2_10 | .1% }} | {{ r.pair.max_rel_error_m2_10 | .1% }} | 10 % | {{ r.pair.linear_window.0 | .1f }} – {{ r.pair.linear_window.1 | .1f }} | {{ r.pair.eps_B_peak | .3f }} | {{ r.pair.t_peak | .1f }} |
| `electron_ion` | {{ r.electron_ion.median_rel_error_m2_10 | .1% }} | {{ r.electron_ion.max_rel_error_m2_10 | .1% }} | reported | {{ r.electron_ion.linear_window.0 | .1f }} – {{ r.electron_ion.linear_window.1 | .1f }} | {{ r.electron_ion.eps_B_peak | .3f }} | {{ r.electron_ion.t_peak | .1f }} |

ε_B is the magnetic energy divided by the initial kinetic energy. Gauss's law holds to
{{ r.pair.gauss_max | .1e }} (`pair`) and {{ r.electron_ion.gauss_max | .1e }} (`electron_ion`).

## Notes

- **One common fit window.** Many modes grow at once, and all of them leave the linear regime when the
  fastest ones saturate. The window is therefore taken from the total magnetic energy: from two e-folds
  (in amplitude) above its initial noise to two e-folds below its peak.
- **Mode m = 1** is excluded from the gate: it grows slowest and is still buried in noise when the fast
  modes saturate.
- **`electron_ion`** is reported, not gated. Its measured rates follow the theory curve at small k but
  lie systematically below it at large k (median deviation {{ r.electron_ion.median_rel_error_m2_10 | .1% }});
  the cause has not been investigated.

## Reproduce

```bash
python -m pypic.launch -n 4 examples/filamentation/run.py --config examples/filamentation/configs/pair.yaml \
    --out output/wiki/filamentation/pair
python examples/filamentation/analyze.py output/wiki/filamentation/pair --movie
```

## Data provenance

{{ table:provenance }}

## References

- E. S. Weibel, "Spontaneously growing transverse waves in a plasma due to an anisotropic velocity
  distribution", Phys. Rev. Lett. 2, 83 (1959).
- B. D. Fried, "Mechanism for instability of transverse plasma waves", Phys. Fluids 2, 337 (1959).
- R. C. Davidson, D. A. Hammer, I. Haber and C. E. Wagner, "Nonlinear development of electromagnetic
  instabilities in anisotropic plasmas", Phys. Fluids 15, 317 (1972).
