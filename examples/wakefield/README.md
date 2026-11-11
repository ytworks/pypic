# Plasma wakefield: the wake of a laser pulse and of an electron bunch

A short driver that crosses a cold plasma leaves a plasma wave behind it: every point oscillates at ω_pe, and the pattern follows the driver with the wavenumber k_w = ω_pe / v. The driver is a laser pulse from the antenna (PIC-D20) or a rigid electron bunch loaded on a heavy partner of the opposite charge (PIC-D21).

- **Theory** (`theory.py`; E in m c ω_pe / e, lengths in c/ω_pe): linear, cold plasma, ions at rest.
  - Pulse with the field envelope a₀ exp(−(t/τ)²) exp(−r²/w²): behind it

    E_x = (A / v_g) exp(−2r²/w²) cos(t − k_w x), E_r = (4r/w²) A exp(−2r²/w²) sin(t − k_w x), A = a₀² (τ/4) √(π/2) exp(−τ²/8),

    largest for τ = 2 (Schroeder et al. 2011, OSTI 1013067: the resonant Gaussian pulse). a₀ and w are the local ones of the diffracting beam, and a₀² is 1/η larger inside the plasma than in vacuum (η = √(1 − ω_pe²/ω₀²): the pulse keeps its energy flux through the ramp).
  - The pulse travels with the group velocity of the Yee grid, v_g = 0.99142 at 32 cells per wavelength and dt = dx/2 (0.99499 without the grid).
  - Bunch with the peak density n_b: E_x = n_b √(2π) σ_x exp(−k_w²σ_x²/2) R(r) cos(k_w ξ) behind it, with R(r) = k_w² ∫ r′ exp(−r′²/2σ_r²) I₀(k_w r_<) K₀(k_w r_>) dr′ and k_w = 1/β (Muggli, arXiv:2007.05226, Eqs. 1 and 3, for β = 1).

| config | driver | box | predicted |
|---|---|---|---|
| `laser_linear.yaml` (gate) | a₀ = 0.1, τ = 2, w = 4, ω₀ = 10 ω_pe, focus at the centre of the box | 3200×128², dx = λ₀/32, transverse cells 0.1875, periodic transverse faces | E_x = 3.853e-3 on the axis at the focus; E_r/E_x = 0.4957 at r = w/2 |
| `beam_linear.yaml` (gate) | n_b = 0.05, σ_x = σ_r = 1, γ = 10, bunch and partner of the mass 10⁶ | 900×192², dx = 1/8, reflecting transverse faces | E_x = 0.03505 on the axis; R(σ_r)/R(0) = 0.7516; k_w = 1.00504 |

```bash
python -m pypic.launch -n 8 examples/wakefield/run.py --config examples/wakefield/configs/laser_linear.yaml
python examples/wakefield/analyze.py output/wakefield/laser_linear
```

- **Set-up**: electrons and ions (mass 1836) at 8 particles per cell each, loaded at the same places, so the plasma has no field before the driver arrives. Open x faces; the density rises over a ramp and falls again before the far face, so no plasma touches an open face. The run ends when the driver is at `x_end`.
  - Pulse: the antenna stands in the vacuum before the ramp. The Rayleigh length (80 c/ω_pe) is not long against the box, so a₀² on the axis changes by up to 13 % along it and the theory is evaluated at the place of the measurement.
  - Bunch: it starts inside the plateau. Its partner stays behind as a charge that the plasma screens; between periodic transverse faces that charge would add a uniform E_x of about 2 % of the wake, which the conducting walls remove.
- **Analysis**: on the plane through the axis (for the pulse the one normal to its polarisation, where the laser has no E_x and E_z is the radial field), at the last output, behind the driver.
  - On the axis, a least-squares fit of E_x = c A_th(x) cos(k x + φ) + b with the amplitude A_th of the theory; free c, k, φ, b.
  - On the two lines at r = w/2 (pulse) or r = σ_r (bunch), the amplitudes with the same k and φ: E_x and, for the pulse, E_r = −c_r A_r sin(k x + φ).
- **Gates** (also Gauss ≤ 1e-10):
  - `laser_linear`: |c − 1| ≤ 5 %, |k v_g − 1| ≤ 2 %, |c_r/c_x − 1| ≤ 10 %.
  - `beam_linear`: |c − 1| ≤ 5 %, |k β − 1| ≤ 2 %, amplitude at σ_r over the one on the axis within 10 % of R(σ_r)/R(0).
- **Reported**: the phase φ, the amplitude at the output before the last, the place of the driver.
- **Null run**: `--set drive=false` runs the same box for the same time without the driver. Criteria: max |E_x| on the axis ≤ 1 % of the predicted amplitude, mean electron density on the plateau within 1 % of n₀.
