# Relativistic pair shock: jump conditions and Fermi acceleration

A cold pair flow at γ₀ = 15 hits a reflecting wall, which is the piston method. The counter-streaming of the incoming and reflected plasma is Weibel unstable. A collisionless shock forms and moves back into the flow; the wall frame is the downstream frame.

- **Theory** (`theory.py`): Synge-gas jump conditions.
  - β_sh = θ₂/(γ₀β₀), n₂/n₁ = 1 + β₀/β_sh, with ⟨γ⟩(θ₂) = γ₀.
  - 2-D: β_sh = 0.497, n₂/n₁ = 3.01. 3-D: β_sh = 0.332, n₂/n₁ = 4.01.
  - With a constant Γ = 3/2: β_sh = 0.47, n₂/n₁ = 3.13 (Spitkovsky 2008).
- **Known results** (not analytic):
  - a non-thermal tail dN/dγ ∝ γ^−p with p ≈ 2.4 in 2-D (Spitkovsky 2008b; Sironi et al. 2013);
  - the isotropic-diffusion theory gives p ≈ 2.2–2.3 (Kirk et al. 2000; Achterberg et al. 2001);
  - γ_max ∝ t^{1/2}.

| config | box | use |
|---|---|---|
| `pair_2d.yaml` (gate) | quasi-2-D (4 cells along z), t = 400 ω_pe⁻¹ | jump conditions |
| `pair_2d_long.yaml` | quasi-2-D, t = 2000 | Fermi tail and γ_max(t) |
| `pair_3d.yaml` | 3-D | 3-D jump conditions |

```bash
python -m pypic.launch -n 8 examples/relativistic_shock/run.py --config examples/relativistic_shock/configs/pair_2d.yaml
python examples/relativistic_shock/analyze.py output/relativistic_shock/pair_2d --movie
```

- **Numerics**:
  - A compensated binomial current filter (PIC-D10) suppresses the numerical Cherenkov instability of the γ₀ = 15 flow.
  - Both x walls are reflecting. The vacuum gap that opens at the far wall bounds the useful time to t < L/(1 + β_sh).
- **Gate**: β_sh and n₂/n₁ within 10 % of the jump conditions for the run's effective dimensionality.
- **Reported**: downstream temperature, tail index, γ_max(t), magnetic energy fraction. The electron tracks record E and B at the particle, for the work-by-component analysis (`pypic.analysis.particles.work_by_component`).
