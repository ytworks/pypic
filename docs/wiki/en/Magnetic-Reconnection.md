**English** · [日本語](Magnetic-Reconnection-ja) · source: [`examples/reconnection`](https://github.com/ytworks/pypic/tree/main/examples/reconnection)

# Relativistic magnetic reconnection

Two thin current sheets separate regions of opposite magnetic field in a strongly magnetised
electron–positron plasma: magnetisation σ = B₀²/(n_b,tot mc²) = {{ c.pair_sigma10.params.sigma }}, so the
magnetic energy density B₀²/2 is σ/2 = {{ c.pair_sigma10.params.sigma | half }} times the rest-mass energy
density of the background plasma. The
sheets tear into chains of magnetic islands (plasmoids), field lines break and reconnect, and the released
magnetic energy goes into particles, which end up in a power-law spectrum. In the large production box the
measured power-law index is p = {{ r.pair_sigma10_large.power_law_index | .2f }} and
{{ r.pair_sigma10_large_run2.power_law_index | .2f }} in two realizations, against
{{ r.pair_sigma10_large.werner2016_index | .2f }} from the σ-scaling of Werner et al. (2016).

## Physics

**Analytic equilibrium.** A periodic double Harris sheet: with half-thickness δ, sheets at y₁ and y₂,

```math
B_x(y) = B_0\left[\tanh\frac{y-y_1}{\delta} - \tanh\frac{y-y_2}{\delta} - 1\right],
\qquad n(y) = n_{cs}\,\mathrm{sech}^2\frac{y-y_i}{\delta} + n_b .
```

The sheet particles drift at ±β_d along z to carry the current. Ampère's law gives n_cs β_d = B₀/(2δ), and
pressure balance gives B₀²/2 = 2(n_cs/Γ_d)θ_cs. With σ = B₀²/(n_b,tot mc²) = {{ c.pair_sigma10.params.sigma }}:
B₀ = {{ c.pair_sigma10.planner.B0 | .3f }}, β_d = {{ c.pair_sigma10.planner.beta_d | .3f }},
θ_cs = {{ c.pair_sigma10.planner.theta_cs | .3f }}, and the Alfvén speed v_A = c√(σ/(1+σ)) =
{{ c.pair_sigma10.planner.v_A | .3f }} c.

**Empirical comparisons (not analytic).** Relativistic pair reconnection proceeds at a rate of about
0.1 v_A, and the particle spectrum is a power law dN/dγ ∝ γ^−p whose index depends on σ:
p ≈ 1.9 + 0.7/√σ (Werner et al. 2016); Sironi & Spitkovsky (2014) report p ≈ 4, 3, 2, 1.5 for σ = 1, 3,
10, 50.

## Experimental conditions

Normalised units: c = ε₀ = 1, lengths in c/ω_pe (of the background), times in 1/ω_pe. Quasi-2-D (4 cells
along z), 4 cells per skin depth, periodic in every direction; each species has one background beam and
two sheet beams. A compensated binomial current filter is applied.

{{ table:conditions }}

<details><summary>Planner parameters</summary>

{{ table:params }}

</details>

## Animation

![reconnection](images/reconnection/pair_sigma10_large_run2.gif)

`pair_sigma10_large` (second realization, 400 × 200 c/ω_pe, to t = 2000).

The black lines are magnetic field lines: contours of the flux function A_z (B = ∇×(A_z ẑ)), drawn at
the same levels in every frame. × marks X points, where field lines break and reconnect; ○ marks O points,
the centres of the magnetic islands (plasmoids).

- **Top left, J_z with field lines**: the straight, antiparallel field lines on both sides of each sheet
  break at a row of X points and reconnect into closed loops around O points. The islands move along the
  sheet and merge; field lines keep reconnecting at the X points between them, until each sheet has
  become one or two big islands.
- **Top right, electron density with field lines**: the plasma is swept into the islands, which become
  dense blobs bounded by closed field lines.
- **Bottom left, reconnected flux Ψ(t)** of the lower sheet: |A(O) − A(X)| between the main O point and the deepest X point, which
  grows each time field lines are reconnected.
- **Bottom right, electron spectrum**: from the initial thermal distributions (grey) a hard power-law tail
  grows to γ ≈ 50 and beyond.

## Field lines and reconnection

![field lines at six times](images/reconnection/fieldlines.png)

`pair_sigma10_large` (second realization). J_z with the field lines (the same contour levels of A_z in
every panel) at t = 0 and at the times when the reconnected flux of the lower sheet reaches 15, 40, 65 and
90 % of its final value, and at the end; Ψ(t) below with those times marked. Ψ = |A(O) − A(X)| is measured
between the marked points themselves (the main O point and the deepest X point of the lower sheet), so the
curve starts when the tearing has made X and O points (before that there is nothing to measure).

- **t = 0**: the field lines are straight and antiparallel across each sheet: the Harris equilibrium.
- **Tearing**: the sheets break into a chain of X points (×) and small islands (○); the field lines
  that crossed at the X points are now closed loops around the islands.
- **Merging**: neighbouring islands coalesce; the field lines between them reconnect at a new X point, and
  Ψ rises in steps (bottom).
- **End**: each sheet is left with a few large islands and few X points; Ψ levels off as the free
  magnetic energy is used up.

## Comparison with published results

### Power-law index against σ

![power-law index against sigma](images/reconnection/power_law.png)

The measured index of each run at σ = 10 (dots with fit errors; small horizontal offsets only separate
the points), the Werner et al. (2016) fit (dashed) and the Sironi & Spitkovsky (2014) values (squares).
The large box lands on the Werner curve; the small gate box (100 × 100, t = 400) gives a steeper
spectrum, because the power law needs a larger system and a longer time to develop.

### Rate, spectrum and current structure

![large run](images/reconnection/pair_sigma10_large-reconnection.png)

`pair_sigma10_large` (first realization). **Left**: J_z at the time of the peak reconnection rate.
**Middle**: the reconnection rate as the example's analysis reports it, d(Ψ_y₁)/dt/(B₀v_A), against the
published ≈ 0.1. Ψ_y₁ = max_x A_z − min_x A_z is taken along the fixed line y = y₁ (the initial position of
the lower sheet). This is a different measure from the X/O-point flux of the field-line chart: when the sheet
bends, the fixed line also crosses the bulging field lines and Ψ_y₁ grows although no extra flux reconnects. **Right**: the final spectrum with the power-law fit over σ/2 ≤ γ − 1 ≤ 5σ.

<details><summary>The second realization, the gate run and the run in twice the box</summary>

![large run, second realization](images/reconnection/pair_sigma10_large_run2-reconnection.png)
![gate run](images/reconnection/pair_sigma10-reconnection.png)
![xl run](images/reconnection/pair_sigma10_xl-reconnection.png)

</details>

| run | box [c/ω_pe] | p (dN/dγ ∝ γ^−p) | Werner 2016 | rate on y = y₁, peak | rate on y = y₁, mean while active | magnetic energy converted | E∥ share, top 1 % | equilibrium deviation / B₀ |
|---|---|---|---|---|---|---|---|---|
| `pair_sigma10` (gate) | 100 × 100 | {{ r.pair_sigma10.power_law_index | .2f }} ± {{ r.pair_sigma10.power_law_err | .2f }} | {{ r.pair_sigma10.werner2016_index | .2f }} | {{ r.pair_sigma10.reconnection_rate_peak | .3f }} | {{ r.pair_sigma10.reconnection_rate_mean_active | .3f }} | {{ r.pair_sigma10.magnetic_energy_converted | .0% }} | {{ r.pair_sigma10.E_parallel_share_top1pct | +.3f }} | {{ r.pair_sigma10.equilibrium_max_dev_over_B0 | .4f }} |
| `pair_sigma10_large` | 400 × 200 | {{ r.pair_sigma10_large.power_law_index | .2f }} ± {{ r.pair_sigma10_large.power_law_err | .2f }} | {{ r.pair_sigma10_large.werner2016_index | .2f }} | {{ r.pair_sigma10_large.reconnection_rate_peak | .3f }} | {{ r.pair_sigma10_large.reconnection_rate_mean_active | .3f }} | {{ r.pair_sigma10_large.magnetic_energy_converted | .0% }} | {{ r.pair_sigma10_large.E_parallel_share_top1pct | +.3f }} | {{ r.pair_sigma10_large.equilibrium_max_dev_over_B0 | .4f }} |
| `pair_sigma10_large_run2` | 400 × 200 | {{ r.pair_sigma10_large_run2.power_law_index | .2f }} ± {{ r.pair_sigma10_large_run2.power_law_err | .2f }} | {{ r.pair_sigma10_large_run2.werner2016_index | .2f }} | {{ r.pair_sigma10_large_run2.reconnection_rate_peak | .3f }} | {{ r.pair_sigma10_large_run2.reconnection_rate_mean_active | .3f }} | {{ r.pair_sigma10_large_run2.magnetic_energy_converted | .0% }} | {{ r.pair_sigma10_large_run2.E_parallel_share_top1pct | +.3f }} | {{ r.pair_sigma10_large_run2.equilibrium_max_dev_over_B0 | .4f }} |
| `pair_sigma10_xl` | 800 × 400 | {{ r.pair_sigma10_xl.power_law_index | .2f }} ± {{ r.pair_sigma10_xl.power_law_err | .2f }} | {{ r.pair_sigma10_xl.werner2016_index | .2f }} | {{ r.pair_sigma10_xl.reconnection_rate_peak | .3f }} | {{ r.pair_sigma10_xl.reconnection_rate_mean_active | .3f }} | {{ r.pair_sigma10_xl.magnetic_energy_converted | .0% }} | {{ r.pair_sigma10_xl.E_parallel_share_top1pct | +.3f }} | {{ r.pair_sigma10_xl.equilibrium_max_dev_over_B0 | .4f }} |

- **Gate** (`pair_sigma10`): before the tearing onset (t ≤ 10) the B_x(y) profile stays within 5 % of B₀
  of the analytic equilibrium, and Gauss's law holds ({{ r.pair_sigma10.gauss_max | .1e }}). Everything
  else is reported.
- **Rates**: the peak comes from the initial tearing. The mean over the active phase is close to the
  published 0.1 v_A in the gate box and lower in the large box. There the rate shows the initial peak and
  then several later bursts (middle panel above), which the example README attributes to plasmoid mergers.
- **Twice the box and twice the time** (`pair_sigma10_xl`, 800 × 400 c/ω_pe, to t = 4000 ω_pe⁻¹,
  {{ c.pair_sigma10_xl.particles | , }} particles): the index is
  {{ r.pair_sigma10_xl.power_law_index | .2f }} ± {{ r.pair_sigma10_xl.power_law_err | .2f }}, against
  {{ r.pair_sigma10_large.power_law_index | .2f }} and {{ r.pair_sigma10_large_run2.power_law_index | .2f }}
  in the half-size box and {{ r.pair_sigma10_xl.werner2016_index | .2f }} from Werner et al. (2016). The
  index does not depend on the box within 0.06. The high-energy cut-off of the spectrum is 17–32 % higher
  than in the half-size box, not twice as high, and equal to it at equal time; the example README has the
  table.
- **Acceleration channel**: the 1 % most energetic tracked electrons gained essentially nothing from the
  field component along B (E∥ share ≈ 0); their energy came from E⊥, the motional field of the outflows
  and plasmoids. The tracks are sampled every 2.5 ω_pe⁻¹, so this split is indicative.

## Notes

- **Two realizations of the large run**: same configuration and seed, at two commits that do not change
  this set-up. The atomic GPU current deposit has a non-deterministic order, and turbulence amplifies the
  round-off differences, so the runs are independent realizations. The index agrees between them to
  within the fit error ({{ r.pair_sigma10_large.power_law_index | .3f }} and
  {{ r.pair_sigma10_large_run2.power_law_index | .3f }}); the mean rate is
  {{ r.pair_sigma10_large.reconnection_rate_mean_active | .3f }} and
  {{ r.pair_sigma10_large_run2.reconnection_rate_mean_active | .3f }} v_A.

## Reproduce

```bash
python -m pypic.launch -n 4 examples/reconnection/run.py --config examples/reconnection/configs/pair_sigma10.yaml \
    --out output/wiki/reconnection/pair_sigma10
python examples/reconnection/analyze.py output/wiki/reconnection/pair_sigma10 --movie
```

## Data provenance

{{ table:provenance }}

## References

- G. R. Werner, D. A. Uzdensky, B. Cerutti, K. Nalewajko and M. C. Begelman, "The extent of power-law
  energy spectra in collisionless relativistic magnetic reconnection in pair plasmas", ApJL 816, L8 (2016).
- L. Sironi and A. Spitkovsky, "Relativistic reconnection: an efficient source of non-thermal particles",
  ApJL 783, L21 (2014).
- G. R. Werner, D. A. Uzdensky, M. C. Begelman, B. Cerutti and K. Nalewajko, "Non-thermal particle
  acceleration in collisionless relativistic electron–proton reconnection", MNRAS 473, 4840 (2018).
- E. G. Harris, "On a plasma sheath separating regions of oppositely directed magnetic field",
  Nuovo Cimento 23, 115 (1962).
