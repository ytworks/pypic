# Current filamentation (Weibel) instability

Counter-streaming beams along ±x are unstable to transverse perturbations, k ⊥ v. Current filaments form, and so does the magnetic field that later mediates collisionless shocks.

- **Box**: thin along x (4 cells), which excludes the two-stream modes; it resolves k in the (y, z) plane.
- **Theory** (`pypic.theory.streaming.filamentation_growth`): cold relativistic fluid,

  Γ⁴ + Γ²(k² + Σω_ps²/γ³) − k² Σω_ps²β²/γ = 0 (c = 1).

  The response along the drift carries the longitudinal mass γ³m. Γ tends to βω_p/√γ for k ≫ ω_p/c.

| config | set-up | Γ_max (cold) |
|---|---|---|
| `pair.yaml` (gate) | electron and positron beams, u₀ = ±0.5 | 0.598 ω_pe |
| `electron_ion.yaml` | electron beams u₀ = ±1 over ions at rest | — (reported) |

```bash
python -m pypic.launch -n 8 examples/filamentation/run.py --config examples/filamentation/configs/pair.yaml
python examples/filamentation/analyze.py output/filamentation/pair --movie
```

- **Analysis**: B_z(k_y) and B_y(k_z) mode amplitudes, fitted over one common linear window. The window is taken from the total magnetic energy, from two e-folds above its noise to two e-folds below its peak, because all modes leave the linear regime when the fastest ones saturate.
- **Gate**: median |Γ_fit − Γ_th|/Γ_th ≤ 10 % over m = 2…10 (measured: 2.2 % at 16 particles per cell, 1.1 % at 64).
- **Saturation**: compared with the magnetic-trapping estimate B ≈ γ₀mΓ²/(qβ₀k).
