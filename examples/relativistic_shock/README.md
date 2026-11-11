# Relativistic pair shock: jump conditions and Fermi acceleration

A cold pair flow at γ₀ = 15 hits a reflecting wall, which is the piston method. The counter-streaming of the incoming and reflected plasma is Weibel unstable. A collisionless shock forms and moves back into the flow; the wall frame is the downstream frame.

- **Theory** (`theory.py`): Synge-gas jump conditions.
  - β_sh = θ₂/(γ₀β₀), n₂/n₁ = 1 + β₀/β_sh, with ⟨γ⟩(θ₂) = γ₀.
  - 2-D: β_sh = 0.497, n₂/n₁ = 3.01. 3-D: β_sh = 0.332, n₂/n₁ = 4.01.
  - With a constant Γ = 3/2: β_sh = 0.47, n₂/n₁ = 3.13 (Spitkovsky 2008).
  - Magnetised flow (`jump_conditions_magnetised`): a field B₁ perpendicular to the flow with σ = B₁² / (γ₀ n₁ m) and the motional field E = −v × B. Particle, momentum and energy flux with the Maxwell stress and the Poynting flux, and flux freezing B₂/B₁ = n₂/n₁, in the downstream frame. σ → 0 gives the relations above; γ₀ → ∞ in 3-D gives the form of Kennel & Coroniti (1984), which was checked against the expression as remembered and not against the printed paper (`docs/plans/2026-11-11-more-examples.md`).
  - σ = 0.3, field out of the plane (2-D gas): β_sh = 0.645, n₂/n₁ = B₂/B₁ = 2.55.
- **Known results** (not analytic):
  - a non-thermal tail dN/dγ ∝ γ^−p with p ≈ 2.4 in 2-D (Spitkovsky 2008b; Sironi et al. 2013);
  - the isotropic-diffusion theory gives p ≈ 2.2–2.3 (Kirk et al. 2000; Achterberg et al. 2001);
  - γ_max ∝ t^{1/2}.

| config | box | use |
|---|---|---|
| `pair_2d.yaml` (gate) | quasi-2-D (4 cells along z), t = 400 ω_pe⁻¹ | jump conditions |
| `pair_2d_sigma.yaml` (gate) | quasi-2-D, t = 400, σ = 0.3, B₁ along z (out of the plane) | magnetised jump conditions |
| `pair_2d_long.yaml` | quasi-2-D, t = 2000 | Fermi tail and γ_max(t) |
| `pair_2d_t4000.yaml` | quasi-2-D, t = 4000, twice the box | the same, later |
| `pair_2d_tail.yaml` | quasi-2-D, t = 8000, four times the box, 3.3 × 10⁸ particles | the same, at the time of the published tail |
| `pair_3d.yaml` | 3-D | 3-D jump conditions |

```bash
python -m pypic.launch -n 8 examples/relativistic_shock/run.py --config examples/relativistic_shock/configs/pair_2d.yaml
python examples/relativistic_shock/analyze.py output/relativistic_shock/pair_2d --movie
```

- **Numerics**:
  - A compensated binomial current filter (PIC-D10) suppresses the numerical Cherenkov instability of the γ₀ = 15 flow.
  - Both x walls are reflecting. The vacuum gap that opens at the far wall bounds the useful time to t < L/(1 + β_sh).
  - With `sigma` > 0 the flow carries E = −v × B as the initial field `E0`, which is tangential to the reflecting walls (PIC-D18). The wall zeroes the tangential field on its nodes, so it launches an electromagnetic front into the plasma at t = 0. `field: out` puts B₁ along z, `field: in` along y.
- **Gate**: β_sh and n₂/n₁ within 10 % of the jump conditions for the run's effective dimensionality. `pair_2d_sigma` also: the mean downstream field over B₁ within 10 % of the density jump of the theory. Gauss's law ≤ 1e-10.
- **Reported**: downstream temperature, tail index, γ_max(t), magnetic energy fraction. The electron tracks record E and B at the particle, for the work-by-component analysis (`pypic.analysis.particles.work_by_component`).

## Production results (8 GPUs, recorded in `bench/results/`)

| run | box | shock speed | n₂/n₁ | tail | γ_max(t) |
|---|---|---|---|---|---|
| `pair_2d` (gate) | quasi-2-D, t = 400 | 0.472 (−5.1 % vs 2-D Synge) | 3.06 (+1.7 %) | thermal only | — |
| `pair_2d_long` | quasi-2-D, t = 2000 | −5.9 / −7.6 % | +6.9 / +8.1 % | p = 3.91 / 3.85 on γ ≈ 40–100; 15–16 % non-thermal | ∝ t^0.36 / t^0.41 |
| `pair_2d_t4000` | quasi-2-D, t = 4000 | 0.460 (−7.5 %) | 3.16 (+5.1 %) | p = 4.00 ± 0.06 on γ ≈ 40–121; 16 % non-thermal | 172 at the end; ∝ t^0.21 |
| `pair_2d_tail` | quasi-2-D, t = 8000 | 0.463 (−7.0 %) | 3.11 (+3.6 %) | p = 3.86 ± 0.10 on γ ≈ 43–136; 18 % non-thermal | 195 at the end; ∝ t^0.10 |
| `pair_3d` | 3-D, t = 600 | 0.320 (−3.7 % vs 3-D Synge 0.332) | 4.15 (+3.5 % vs 4.01) | thermal | — |

- **Two realizations**: the long run has two recorded runs, each a record of its own in `bench/results/`: one at 200af13, and run 40c0ca85 at 20b2e74. The engine changes between the two commits do not touch this configuration, and the seed is the same. The two values in each cell come from those two runs.
  - The GPU deposit adds currents atomically in a non-deterministic order, so round-off differences grow in this turbulent flow and the runs are independent realizations. The spread is a measure of the statistical uncertainty.
  - An earlier run (−4.2 %, +7.8 %, p = 3.93, t^0.42) was quoted here before but never recorded, so it is no longer quoted.
- **Dimensionality**: the 3-D run confirms the jump to the 3-D gas values (β_sh ≈ 1/3, n₂/n₁ ≈ 4), while the quasi-2-D runs follow the 2-D gas values (0.50, 3.0).
- **No extended tail up to t = 8000 ω_pe⁻¹**. `pair_2d_t4000` (run 014c11db) and `pair_2d_tail` (run b012ba98), both at 2945145 with one run each, double the box and the time twice. t = 8000 ω_pe⁻¹ is 2070 relativistic plasma times ω_p,rel⁻¹ = √γ₀/ω_pe; Spitkovsky's p ≈ 2.4 was measured at ≈ 2250.
  - **Index**: the slope above the thermal peak stays at 3.9–4.0 and does not flatten: 3.91 / 3.85 at t = 2000, 4.00 ± 0.06 at t = 3900, 3.86 ± 0.10 at t = 7800, over the window from 3 × the peak of γ dN/dγ to 0.7 × the cut-off. In `pair_2d_tail` the automatic search finds no straight range, so `tail_index_p` is NaN there and only the window fit exists.
  - **Cut-off** (γ that 10⁻⁴ of the downstream particles exceed): 144–153 at t = 2000 in the four runs, 172 at t = 4000 in the two new runs, 183 from t = 4500 to 7500 and 195 at t = 8000. From t = 4000 to 8000 it rises by 13 %, where t^1/2 gives 41 %. The exponent fitted over the second half of each run falls from 0.36 / 0.41 to 0.21 and 0.10. The histogram has 120 logarithmic bins up to γ = 1500, so the values come in steps of 6.3 % and are not capped by it.
  - **Non-thermal fraction**: 15–16 %, 16 % and 18 % of the downstream particles.
  - **Magnetic energy**: over the initial kinetic energy of the box it is 0.0079 / 0.0067, 0.0040 and 0.0021 at the end of the runs. The box and its kinetic energy double from one run to the next, so the magnetic energy in the box at t = 4000 and 8000 is 1.02 and 1.07 times that of the first t = 2000 run (1.20 and 1.26 times the second), while the shocked volume grows four times.
  - **Meaning**: time alone does not bring the published tail in this set-up. An earlier version of this section expected a converged index at t ≳ 8000 with about 10⁹ particles; the time was reached with 3.3 × 10⁸ particles (L = 1.57 c t) and neither the index nor the cut-off approaches the published values.
  - **Not separated**: the transverse size (62 c/ω_pe, which is 16 relativistic skin depths, in all runs), 8 particles per cell, the 4-pass current filter, the 4 cells along z, and the region of the spectrum (0.3–0.8 of the shock position, the body of the downstream and not a slab behind the shock). None of these was varied.
  - Energy budget +0.29 % and +0.56 % (at step 20,000 of 20,656 and 40,000 of 41,312), Gauss's law ≤ 1e-14; 10 min and 41 min on 8 GPUs, FP64.
- **Shock-speed fit**: the speed is fitted between 0.3 and 0.7 of the validity time. Later on, the reflected-particle precursor reaches the receding plasma edge, and an earlier window that ran to 0.95 t_valid gave −20 % in 3-D.

## Magnetised shock (8 GPUs, FP64, one run each)

`pair_2d_sigma.yaml`: γ₀ = 15, σ = 0.3, B₁ = 3.0 along z, 1792×80×4 cells, 4.6e6 particles, 2066 steps (7 s). The limits were fixed before the first run (`docs/plans/2026-11-11-more-examples.md`), and the first run passed them. The other three rows are `--set` variants of it (`field=in`, `sigma=0.1`, `sigma=1`); the criteria are evaluated and are not gates there. These are targeted runs without a record in `bench/results/`. The gate record of `pair_2d_sigma` is in `bench/results/2026-11-12-c54e5d104-gpu-gates.json` (and `…-mixed.json`): β_sh −7.2 %, n₂/n₁ +4.2 % (mixed +4.4 %), field +4.0 %; the second digit of the density jump differs between runs with the same seed.

| run | jump used | β_sh | n₂/n₁ | ⟨B⟩/B₁ downstream | precursor δB/B₁ (rms, max) | Gauss |
|---|---|---|---|---|---|---|
| σ = 0.3, B₁ along z (gate) | 2-D gas: 0.645, 2.547 | 0.598 (−7.2 %) | 2.666 (+4.7 %) | 2.651 (+4.1 %) | 0.71, 1.9 | 1.2e-12 |
| σ = 0.3, B₁ along y | 3-D gas: 0.549, 2.817 | 0.571 (+4.0 %) | 2.800 (−0.6 %) | 2.746 (−2.5 %) | 0.59, 2.2 | 1.5e-13 |
| σ = 0.1, B₁ along z | 2-D gas: 0.561, 2.780 | 0.522 (−6.8 %) | 2.995 (+7.7 %) | 2.972 (+6.9 %) | 1.08, 3.5 | 1.1e-12 |
| σ = 1, B₁ along z | 2-D gas: 0.783, 2.274 | 0.809 (+3.3 %) | 2.262 (−0.6 %) | 2.281 (+0.3 %) | 0.34, 0.73 | 4.7e-13 |

- **The jump follows the magnetised relations within 8 %** in the four runs, the gate within 7.2 %. The deviations of β_sh and n₂/n₁ have opposite signs in every run, as mass conservation asks (β_sh = β₀ / (r − 1) in the wall frame).
- **Flux freezing**: the mean downstream field over B₁ is within 2 % of the measured density jump in the four runs.
- **A field in the plane gives the 3-D gas**: the particles gyrate out of the simulation plane, and the run follows the 3-D jump (n₂/n₁ −0.6 %), not the 2-D one (+10 %).
- **Precursor**: an electromagnetic wave train runs ahead of the shock in all four runs; its amplitude falls with σ between 0.1 and 1, and the density in it deviates by 1 % rms. Its leading edge has travelled at c since the start of the run, so it contains the front that the wall launches at t = 0 and the emission of the forming shock; the two are not separated, and no comparison with a published precursor amplitude was made.
- **No non-thermal tail** at t = 397 in any of the four runs.
- **Not shown**: the dependence on the seed, on the box width and on the run length; a 3-D magnetised run.
