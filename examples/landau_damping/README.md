# Langmuir oscillations: Landau damping, Dawson oscillation, wave breaking

| config | set-up | exact / theoretical result | gate |
|---|---|---|---|
| `landau.yaml` | warm electrons (θ_e = 10⁻³), n = n₀(1 + α cos kx), kλ_De = 0.5, α = 0.01 | kinetic root of the Maxwellian dielectric function: ω = 1.4157, γ = −0.1534 ω_pe | γ within 10 %, ω within 2 % |
| `nonlinear.yaml` | the same with α = 0.5 | initial damping, then trapping: field-amplitude oscillation at the bounce frequency, phase-space vortices | reported |
| `dawson.yaml` | cold, v = v₀ sin(kx₀), kv₀ = 0.5 ω_pe | every sheet: x = x₀ + (v₀/ω_pe) sin(kx₀) sin(ω_pe t), exactly | frequency 1 %, amplitude 2 % |
| `wavebreaking.yaml` | cold, kv₀ = 1.2 ω_pe | first sheet crossing at ω_pe t = arcsin(ω_pe/kv₀) = 0.985 | crossing time within 5 % |
| `relativistic.yaml` | cold, u₀ = 0.3145 (β_max = 0.3) | exact relativistic sheet period from γ + ω_pe²ξ²/2 = γ₀; T ≈ (2π/ω_pe)(1 + 3β²/16) | sheet frequency within 0.5 % |

```bash
python -m pypic.launch -n 8 examples/landau_damping/run.py --config examples/landau_damping/configs/landau.yaml
python examples/landau_damping/analyze.py output/landau_damping/landau --movie
```

- **Initial state**:
  - The density perturbation is created by moving the loaded electrons to x with x + (α/k) sin kx = x₀, which gives exactly n₀(1 + α cos kx).
  - This move is deposited as a current with the charge-conserving Esirkepov scheme, so the initial E satisfies Gauss's law to round-off (PIC-D9).
  - For the explicit scheme, the momenta are then pushed back by Δt/2.
  - The cold cases start at equilibrium with a velocity perturbation and no field.
- **Output**:
  - The modes of E_x are written at the Nyquist rate of 1.3 ω and with at least 10 samples per damping time.
  - Tracked electrons give the sheet orbits.
  - Phase-space snapshots feed the movie.
