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
| `pair_sigma10_xl.yaml` | 800 × 400, t = 4000 | cut-off of the power law against the box size |

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

`pair_sigma10_large.yaml` runs a 400 × 200 c/ω_pe box to t = 2000 with 10⁸ particles.

It has two recorded runs, each a record of its own in `bench/results/`: one at 200af13, and run e879cfd6 at d9637cd. The engine changes between the two commits do not touch this configuration. The seed is the same; the runs differ only through the non-deterministic order of the atomic current deposit, which turbulence amplifies, so they are independent realizations.
- **power law**: dN/dγ ∝ γ^−p with p = 2.10 ± 0.03 in both runs, over σ/2 ≤ γ−1 ≤ 5σ, against 1.9 + 0.7/√σ = 2.12 (Werner et al. 2016), within 1 %;
- **reconnection rate**: 0.064 and 0.074 v_A on average over the active phase, with a peak of 0.19 from the initial tearing and a second burst from plasmoid mergers;
- **energy**: 77 % and 76 % of the magnetic energy is converted;
- **acceleration channel**: the 1 % most energetic tracked electrons gained essentially nothing from E∥ (share +0.005 and −0.005); their energy came from E⊥, the motional field of the outflows and plasmoids. The tracks are sampled every 2.5 ω_pe⁻¹, so this split is indicative.
- An earlier run (p = 2.15, rate 0.080) was quoted here before but never recorded, so it is no longer quoted.

The small gate box (100 × 100, t = 400) gives p = 2.4: the power law needs the larger system and the longer time.

### Twice the box and twice the time (`pair_sigma10_xl.yaml`)

800 × 400 c/ω_pe (3200 × 1600 × 4 cells) to t = 4000 with 3.7 × 10⁸ particles: 32,000 steps in 22 min on 8 GPUs, FP64. One run, 0e8c5abd at 2945145, recorded in `bench/results/`. The gate criteria hold (the sheet profile before the onset, Gauss's law 2.8e-12); the total energy changes by −0.34 %.

| | `large`, t = 2000 (two runs) | `xl`, t = 2000 | `xl`, t = 4000 |
|---|---|---|---|
| index p over σ/2 ≤ γ−1 ≤ 5σ | 2.10 ± 0.03, 2.10 ± 0.03 | 2.05 ± 0.02 | 2.04 ± 0.02 |
| γ−1 above which 10⁻³ of the particles lie | 90, 87 | 90 | 98 |
| the same for 10⁻⁴ | 140, 158 | 158 | 185 |
| the same for 10⁻⁵ | 225, 286 | 264 | 286 |

- **Index**: p = 2.04 ± 0.02 at t = 4000, against 2.12 (Werner et al. 2016), −3.7 %. It has not changed since t = 1000 (2.05, 2.05, 2.04, 2.04 at t = 1000, 2000, 3000, 4000), so it is converged in time in this box. It is 0.06 flatter than in the half-size box, which is 1.8 times the combined fit error of one run of each.
- **Cut-off against the box size**: with twice the box and twice the time the energy that 10⁻⁴ of the particles exceed is 185 against 140 and 158, that is 17–32 % higher, and not twice as high. At equal time (t = 2000) the two boxes agree within the spread of the two `large` runs. Between t = 3000 and t = 4000 the 10⁻³ and 10⁻⁴ values do not change. The spectrum is stored on 160 logarithmic bins, so these energies come in steps of 4 %; no cut-off form was fitted. The values are `pypic.analysis.particles.energy_cutoff` on the `spectrum` output and are not part of the record.
- **Rate**: 0.066 v_A on average over the active phase, peak 0.17 — as in the half-size box (0.064 and 0.074, peak 0.19).
- **Energy**: 78 % of the magnetic energy is converted (77 % and 76 % in the half-size box).
- **Acceleration channel**: E∥ share of the 1 % most energetic tracked electrons −0.004, as before.
- **Not shown**: a second realization of this box; a box of a third size, which a scaling of the cut-off with the box would need.
