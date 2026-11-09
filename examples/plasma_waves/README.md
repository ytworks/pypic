# Plasma waves from thermal noise (omega–k dispersion)

A uniform Maxwellian plasma at rest is never quiet. Its thermal fluctuations excite every normal
mode of the plasma. The power of a space–time Fourier transform, |F(k, ω)|², of transverse-averaged
lineouts therefore shows the dispersion branches, ω(k), which can be compared directly with linear
theory. No initial perturbation is needed.

| config | B0 | component: branch (theory, `pypic.theory.waves`) |
|---|---|---|
| `unmagnetized.yaml` | 0 | E_x: Langmuir, kinetic root of the Maxwellian ε (Bohm–Gross for kλ_D ≪ 1); E_y, E_z: EM waves, ω² = ω_pe² + c²k² with the Yee/leapfrog correction |
| `perpendicular.yaml` | along z, ω_ce = 0.5 ω_pe, θ_e = 10⁻³ | E_z: O mode; E_y: X mode (upper/lower branches, upper-hybrid resonance); E_x: electron Bernstein waves between cyclotron harmonics (bands 1–3 gated) |
| `parallel.yaml` | along x, ω_ce = 1.5 ω_pe, θ_e = 10⁻³ | E_y, E_z: whistler, R and L branches (R/L for k ≤ 2, where their windows are disjoint); E_x: Langmuir |
| `ion_acoustic.yaml` | 0 | E_x: ion-acoustic waves, kinetic root with electron and ion Landau damping (T_e/T_i = 50, m_i/m_e = 25) |

```bash
python -m pypic.launch -n 8 examples/plasma_waves/run.py --config examples/plasma_waves/configs/unmagnetized.yaml
python examples/plasma_waves/analyze.py output/plasma_waves/unmagnetized --movie
```

- **Set-up** (`plan.py`):
  - The box is quasi-1-D along k ∥ x, with 8×8 transverse cells that the lineouts average over.
  - Δx = λ_De.
  - The frequency resolution is Δω = 2π / t_end.
- **Output cadence**: the planner gives the highest frequency in the analysis window as `shortest_period`. The lineouts are then written at its Nyquist rate and time-averaged over each interval. The average acts as an anti-aliasing filter for the high-k electromagnetic noise, whose frequency reaches about 2/Δt.
- **Analysis**:
  - Each theory branch is located by a peak search within ±15 % of the theoretical ω at every k.
  - A branch passes when it is found at at least half of the points.
  - At ≥ 90 % of the found points, |ω_meas − ω_th|/ω_th ≤ max(3 %, 2Δω/ω).
- **Gates**: `unmagnetized` and `perpendicular` are GPU gates (`scripts/gpu_gates.py`).

- **Measured**:
  - unmagnetized: Langmuir 1.4 %, EM 0.3 % (median).
  - perpendicular: O / X 0.05–0.07 %, Bernstein bands 1–3 0.5–0.7 %.
  - parallel: whistler 0.4 %, R 0.08 %, L 0.01 %, Langmuir 0.2 %.
  - ion acoustic: 4.2 % (reported).
- **Why θ_e = 10⁻³ for the magnetized cases**: at θ_e = 10⁻² the relativistic mass increase (⟨γ⟩ − 1 ≈ 1.5 θ_e) shifts the cyclotron harmonics, and the Bernstein bands come out 1.5–4 % below the non-relativistic theory.
