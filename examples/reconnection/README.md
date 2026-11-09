# Relativistic pair magnetic reconnection (double Harris sheet)

A periodic double Harris sheet in a pair plasma with magnetisation σ = B₀²/(n_b,tot mc²) = 10 tears from noise into a chain of plasmoids, which converts magnetic energy into a power-law particle spectrum.

- **Analytic equilibrium** (`theory.py`):
  - B_x = B₀[tanh((y−y₁)/δ) − tanh((y−y₂)/δ) − 1] and n = n_cs sech²((y−y_i)/δ);
  - Ampère's law gives n_cs β_d = B₀/(2δ), and pressure balance gives B₀²/2 = 2(n_cs/Γ_d)θ_cs.
- **Empirical comparisons** (not analytic):
  - reconnection rate ≈ 0.1 v_A;
  - dN/dγ ∝ γ^−p with p ≈ 1.9 + 0.7/√σ = 2.1 (Werner et al. 2016), or p ≈ 2 at σ = 10 (Sironi & Spitkovsky 2014).

| config | box | use |
|---|---|---|
| `pair_sigma10.yaml` (gate) | 100 × 100 c/ω_pe, quasi-2-D | equilibrium, conservation |
| `pair_sigma10_large.yaml` | 400 × 200, t = 2000 | power law, rate |

```bash
python -m pypic.launch -n 8 examples/reconnection/run.py --config examples/reconnection/configs/pair_sigma10.yaml
python examples/reconnection/analyze.py output/reconnection/pair_sigma10 --movie
```

- **Set-up**:
  - Each species has one background beam and two sheet beams (sech² profiles with opposite drifts, PIC-D9).
  - B_x is given as a field profile.
  - A compensated current filter is applied (PIC-D10).
- **Gate**: before the onset (t ≤ 10 ω_pe⁻¹), the B_x(y) profile stays within 5 % of B₀ of the analytic sheet; Gauss's law holds.
- **Reported**:
  - the peak reconnection rate, from the flux function along the lower sheet;
  - the magnetic energy converted;
  - the power-law index;
  - the share of the work done by E∥ on the 1 % most energetic tracked electrons.

## Production result (8 GPUs, recorded in `bench/results/`)

`pair_sigma10_large.yaml` runs a 400 × 200 c/ω_pe box to t = 2000 with 10⁸ particles. It was run twice with the same code and seed; the two runs differ only through the non-deterministic order of the atomic current deposit, which turbulence amplifies, so they are independent realizations.
- **power law**: dN/dγ ∝ γ^−p with p = 2.15 ± 0.03 and 2.10 ± 0.03 over σ/2 ≤ γ−1 ≤ 5σ, against 1.9 + 0.7/√σ = 2.12 (Werner et al. 2016), within 1.5 %;
- **reconnection rate**: 0.080 and 0.064 v_A on average over the active phase, with a peak of 0.19 from the initial tearing and a second burst from plasmoid mergers;
- **energy**: 77 % of the magnetic energy is converted;
- **acceleration channel**: the 1 % most energetic tracked electrons gained essentially nothing from E∥ (share ≈ 0); their energy came from E⊥, the motional field of the outflows and plasmoids. The tracks are sampled every 2.5 ω_pe⁻¹, so this split is indicative.

The small gate box (100 × 100, t = 400) gives p = 2.4: the power law needs the larger system and the longer time.
