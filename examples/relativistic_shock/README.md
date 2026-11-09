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

## Production results (8 GPUs, recorded in `bench/results/`)

| run | box | shock speed | n₂/n₁ | tail | γ_max(t) |
|---|---|---|---|---|---|
| `pair_2d` (gate) | quasi-2-D, t = 400 | 0.472 (−5.1 % vs 2-D Synge) | 3.06 (+1.7 %) | thermal only | — |
| `pair_2d_long` | quasi-2-D, t = 2000 | −5.9 / −7.6 % | +6.9 / +8.1 % | p = 3.91 / 3.85 on γ ≈ 40–100; 15–16 % non-thermal | ∝ t^0.36 / t^0.41 |
| `pair_3d` | 3-D, t = 600 | 0.320 (−3.7 % vs 3-D Synge 0.332) | 4.15 (+3.5 % vs 4.01) | thermal | — |

- **Two realizations**: the long run has two recorded runs, each a record of its own in `bench/results/`: one at 200af13, and run 40c0ca85 at 20b2e74. The engine changes between the two commits do not touch this configuration, and the seed is the same. The two values in each cell come from those two runs.
  - The GPU deposit adds currents atomically in a non-deterministic order, so round-off differences grow in this turbulent flow and the runs are independent realizations. The spread is a measure of the statistical uncertainty.
  - An earlier run (−4.2 %, +7.8 %, p = 3.93, t^0.42) was quoted here before but never recorded, so it is no longer quoted.
- **Dimensionality**: the 3-D run confirms the jump to the 3-D gas values (β_sh ≈ 1/3, n₂/n₁ ≈ 4), while the quasi-2-D runs follow the 2-D gas values (0.50, 3.0).
- **The tail is still young at t = 2000 ω_pe⁻¹**: that is ≈ 520 relativistic plasma times ω_p,rel⁻¹ = √γ₀/ω_pe. Spitkovsky's p ≈ 2.4 was measured at ≈ 2250 ω_p,rel⁻¹, about 8700 ω_pe⁻¹; the index flattens and the cutoff grows (∝ t^1/2) with time. A converged index needs t ≳ 8000 ω_pe⁻¹ and a box L ≳ 1.5 c t, which is about 10⁹ particles in quasi-2-D.
- **Shock-speed fit**: the speed is fitted between 0.3 and 0.7 of the validity time. Later on, the reflected-particle precursor reaches the receding plasma edge, and an earlier window that ran to 0.95 t_valid gave −20 % in 3-D.
