# Alfvénic wave carried by a flow into a shock

A magnetised electron–ion plasma flows at −V₀ along x across the field B₀ = B₀ ŷ. The face x+ injects it and carries in one circularly polarised wave that propagates along B₀ (k = k ŷ): the cold two-fluid R wave. The pattern is uniform along x, so the flow carries it into the box unchanged. In the open box (the gate) the wave crosses the box. With a reflecting wall at x− (the piston) the flow piles up into a perpendicular shock that runs into the plasma that carries the wave.

The example exercises three engine features of `docs/plans/2026-11-13-upstream-state-restart.md`:
- PIC-D26, the turbulent inflow face: fluctuation table, total-field / scattered-field surface on the face plane, `flux` face with the Marder zone;
- PIC-D24, the frozen-in initial state: vector potential with an envelope h(x) that vanishes at the wall, shaped velocity perturbations, motional E from a displacement of the electrons;
- PIC-D25, restart: the piston configuration writes checkpoints.

- **Theory** (`theory.py`, CPU-tested in `tests/cpu/test_alfven_inflow.py`):
  - **Parallel waves.** Every species at rest in B₀ responds with the cold velocity −iω v = (q/m)(E + v × B₀). The current σ(ω) E closes Maxwell's equations, σ E = i(ω − k²/ω) E, on a circular vector of the (z, x) plane. The two roots below the cyclotron resonances are the R wave (it rotates with the electrons, ω > k v_A) and the L wave (ω < k v_A). They satisfy the Stix relations k²/ω² = 1 − Σ ω_ps²/(ω(ω ± Ω_s)). For k d_i → 0 both tend to k v_A / (1 + v_A²)^½. To first order the Hall term splits them: (ω_R − ω_L)/(k v_A) = k d_i (m_i − m_e)/(m_i(m_i + m_e))^½ / (1 + v_A²)². These are the same solver and table as `tests/gpu/mgpu_inflow.py` (agreement 1e-14).
  - **Box frame.** k is normal to the flow, so the phase k y − ωt is the same in both frames (to order V₀²). E = E_rest − v_d × B adds the motional field −V₀ δB_z ŷ. The table holds B, E and the velocity of each species as cos and sin parts.
  - **Piston shock.** Downstream frame = wall frame. With U₁ = V₀ + V_sh, r = U₁/V_sh, so V_sh = V₀/(r − 1). r is the compression of a perpendicular MHD shock at M_A = U₁/v_A (γ = 5/3, `pypic.setup.mhd`); r → 4 for M_A → ∞. The downstream T_e + T_i follows from the momentum flux; the tests check the energy flux too. A frozen-in field keeps the normal δB_x and compresses the tangential δB_z by r.

| config | set-up | grid | use |
|---|---|---|---|
| `gate.yaml` (gate) | open flow, V₀ = 2 v_A; the wave loaded over the whole box and carried in through x+ | 256 × 128 × 2 (64 × 32 × 0.5 d_e), 3000 steps (t = 300 = 1.9 / Ω_ci) | the wave against the incident one, densities, Gauss's law |
| `piston.yaml` | reflecting wall at x−, V₀ = 3 v_A (M_A = 4.39 in the shock frame); the wave loaded with the envelope h(x) and carried in through x+; checkpoints every 4000 steps | 384 × 128 × 2 (96 × 32 × 0.5 d_e), 13000 steps (t = 1300 = 8.1 / Ω_ci) | reported: shock, the wave across it, spectra |
| `piston.yaml --set wave=false` | the piston without the wave (a table of zero amplitude) | as above | null run: the shock's own k-mode |

All: m_i/m_e = 16, ω_pe/Ω_ce = 10 (B₀ = 0.1), θ_e = 10⁻³ and T_i = T_e (β = 0.40, Debye length 0.13 Δx), Δx = 0.25 d_e, Δt = 0.1, 32 particles per cell and species. The wave: R, k d_i = 0.785 (one wavelength, 32 d_e, along y), ω = 1.385 k v_A, |δB| = 0.2 B₀.

```bash
python -m pypic.launch -n 2 examples/alfven_inflow/run.py --config examples/alfven_inflow/configs/gate.yaml
python examples/alfven_inflow/analyze.py output/alfven_inflow/gate
python -m pypic.launch -n 2 examples/alfven_inflow/run.py --config examples/alfven_inflow/configs/piston.yaml
python -m pypic.launch -n 2 examples/alfven_inflow/run.py --config examples/alfven_inflow/configs/piston.yaml \
    --set wave=false --out output/alfven_inflow/piston_null
# continue a piston run from a checkpoint (PIC-D25)
python -m pypic.launch -n 2 examples/alfven_inflow/run.py --config examples/alfven_inflow/configs/piston.yaml \
    --restart output/alfven_inflow/piston/checkpoint
```

- **Set-up** (explicit scheme):
  - **Flow and face.** Electrons and ions enter through x+ as a `flux` face, each species with its own flux. A `paired` face would tie the electrons' birth rate to the ions' flux, which differs from theirs in a wave by the wave current (plan, Results). The motional field E₀ = V₀ B₀ ẑ is uniform. Compensated binomial current filter (4 passes), Marder zone 0.125, CPML of 12 cells.
  - **The wave at t = 0 (gate)**: B_x, B_z and E_x as profiles along y; E_z from a displacement of the electrons along z (ξ = −E/(q n)); every velocity component as a shaped perturbation along y.
  - **The wave at t = 0 (piston)**: everything is multiplied by an envelope h(x), a table of sin² ramps. h is 0 below 4 d_e, 1 above 12 d_e, falls to 0 over 8 d_e and is exactly 0 over the last 1 d_e before the face. B comes from the vector potential: A_z = h a_x(y) for B_x = ∂_y A_z and A_x = −h a_z(y) for B_z = −∂_y A_x, so ∇·B = 0 bitwise. The scale uses the discrete derivative, 2 sin(k Δy/2)/Δy, so B on the Yee grid is the wave's B exactly where h = 1. The ramps add B_y = −h′ a_x, of the order of δB since λ/(2π × ramp) ≈ 0.6–1.3. E_x and E_z come from displacements of the electrons along x and z.
  - **Why h vanishes at the face too.** The engine refuses an A tangential to a non-periodic face that is not zero on its plane (PIC-D24; a wall would zero the normal B), and a displacement along the normal of such a face within one cell of it. A B_x(y) on the wall itself would leave ∇·B ≈ B_x/Δx in the first cell. So either the wall or the face plane has no B_x at t = 0; here the face does. That plasma (the last 9 d_e) is flushed in at V₀ and reaches the shock near t ≈ 850. After that the shock meets the wave that the face carried in.
  - **What is not loaded.** The box-frame E_y = −V₀ δB_z needs the charge v₀ δv/c², which the particles do not carry (PIC-D26), so it is left out and the plasma sets it up. The ramps' extra currents (−h′ δB_z along y, −h″ a_x along z) are not loaded either.
  - **Why the null run keeps a table.** PIC-D26 allows E₀ tangential to the wall only when the open face of that axis carries a fluctuation table. The null run carries one of zero amplitude, so the face holds the uniform motional field as in the wave run. The analysis normalises by the wave of the planner's table in both runs.
  - **Particle arrays.** The capacity factor is 1.2 r for the piston, because the downstream holds r times the density.
- **Gate** (fixed before the first run: the criteria of test (g) in `docs/plans/2026-11-13-upstream-state-restart.md`, amendment of 2026-11-14, unchanged). Measured at the last output (t = 300) in the plane z = 0, per x row: a(x) is the complex amplitude of e^{iky} of B_x and B_z, b the same of the incident wave of the table.
  - |a − b| / mean|b| ≤ 12 % on every row with x > 10 d_e (the outflow face x− carries no incident field);
  - the same ≤ 4 % on average over the rows 10 < x < 60 d_e;
  - the densities of electrons and ions, averaged over bins of 2 d_e from x = 10 d_e, within 10 % of 1;
  - Gauss's law ≤ 1e-10 and ∇·B ≤ 1e-12.
- **Reported (piston)**, not gated:
  - the shock position (the outermost x where the ion density, averaged over d_i/2, exceeds (1 + r)/2) and its speed from a fit over the second half of the run;
  - the downstream density and B_y between 0.2 and 0.8 x_sh, and T_e, T_i (2/3 of the mean kinetic energy per particle in the ions' frame);
  - the k-mode of B_x and B_z ahead (x_sh + 4 to x_sh + 20 d_e, the foot of the reflected ions lies there), far ahead (to lx − 4) and behind (0.2 x_sh to x_sh − 4), against the incident wave;
  - the electron energy spectra downstream and far ahead.

## Result of the gate configuration (2 GPUs, run of a502a2c)

| quantity | FP64 | `mixed` | `compact` | gate |
|---|---|---|---|---|
| B_x row error, largest / mean | 8.30 % / 2.63 % | 8.30 % / 2.63 % | 8.30 % / 2.63 % | 12 % / 4 % |
| B_z row error, largest / mean | 10.86 % / 2.90 % | 10.86 % / 2.90 % | 10.86 % / 2.90 % | 12 % / 4 % |
| densities in bins of 2 d_e | within 1.89 % | within 1.89 % | within 1.89 % | 10 % |
| Gauss's law, ∇·B | 4.7e-15, 1.1e-14 | 3.4e-15, 1.1e-14 | 7.1e-16, 1.1e-14 | 1e-10, 1e-12 |

- **All three precisions pass** and agree to 4 or more digits in every criterion. `scripts/gpu_gates.py --gpus 2 --examples alfven_inflow:gate` passes with each of `--precision fp64`, `mixed` and `compact` (single-gate runs, not recorded).
- **The largest row errors are those of the rows next to the inflow face** (x = 63.75 d_e for B_x, 63.875 d_e for B_z); without the last 5 rows (1.25 d_e) the largest is 5.9 % (B_x) and 6.2 % (B_z). The mean amplitude ratio is 0.986 (B_x) and 0.980 (B_z), the mean phase offset 0.015 rad. The test with V₀ = 0.05 (this gate: 2 v_A = 0.0485) measured 8 % and 2.6 %.
- **A fault found by this gate.** Before d341a49 (merged into this branch) the `compact` run failed with row errors of 72 % (B_x) and 74 % (B_z): with the compact storage the entering particles took the wave's velocity at their index coordinate instead of their physical place. The errors grew from the face inward at the flow speed. FP64 was not affected. The Codex review found the same fault; the fix and its GPU test are on `upstream-state-restart`.
- **Cost**: 9.1 ms per step on 2 RTX PRO 6000 (4.2e6 particles), 27 s for the run.

## Result of the piston (2 GPUs, FP64, one run each, ef4fc4f)

Records: `bench/results/2026-11-14-ef4fc4f08-example-alfven_inflow-piston-run-d6f52e61.json` (with the wave) and `…-run-b4f6699d.json` (the null run, `wave=false`), on the head that the Codex review left without findings. Outputs: `/mnt/nvme/pypic/output/final3/alfven_inflow/piston` and `…/piston_null` (1.5 GB each). The earlier records of c571302 (`…-run-9245d02a.json`, `…-run-d6549bc2.json`) predate the review fixes of the gather next to the face and agree with these to 1 % or better.

| quantity | with the wave | null run | jump conditions |
|---|---|---|---|
| shock speed (fit over t = 650–1300) | 0.0320 c | 0.0367 c | 0.0336 c |
| M_A in the shock frame from the measured speed | 4.32 | 4.51 | 4.39 |
| downstream density n₂ (incoming flow: 1; measured far ahead 0.999 / 0.996) | 2.94 | 2.71 | 3.16 |
| downstream B_y / B₀ | 2.94 | 2.72 | 3.16 |
| downstream T_e + T_i against the jump conditions (0.0280 m_e c²) | 1.015 | 0.988 | 1 |
| downstream T_e, T_i [m_e c²] | 4.6e-3, 2.38e-2 | 3.8e-3, 2.38e-2 | |
| B_x mode: ahead / far ahead / behind (mean of 59 outputs, units of the wave) | 0.97 / 1.04 / 2.40 | 0.014 / 0.012 / 0.62 | frozen-in: behind = ahead |
| B_z mode: ahead / far ahead / behind | 1.13 / 0.92 / 4.25 | 0.028 / 0.032 / 0.64 | frozen-in: behind = 3.16 × ahead |
| B_x, B_z behind at the end (t = 1300) | 2.55, 5.23 | 1.23, 2.14 | |
| Gauss's law, ∇·B | 1.2e-14, 7.4e-14 | 9.1e-15, 7.1e-14 | 1e-10, 1e-12 |

- **Shock.** With the wave the shock is 4.9 % slower and the downstream density 7 % lower than the jump conditions; without it the shock is 9.0 % faster and the density 14 % lower. B_y follows the density in both runs (frozen-in), and T_e + T_i is within 2 % of the jump conditions in both. These are single runs measured at one time; the difference between the two runs is not separated from the run-to-run spread of the shock position, which this example did not measure.
- **The wave ahead of the shock** keeps its amplitude on average (0.97 B_x, 1.13 B_z of the incident wave), but varies with time: B_z is 0.97–1.03 until t = 600, falls to 0.55–0.69 at t = 750–850, when the plasma loaded without the wave at the face (the ramp of h) reaches the shock, and rises to 1.74–1.78 at t = 1100–1200 in the wave that the face carried in (B_x: 0.62–1.19 from t = 600 on). The window ahead holds the foot of the reflected ions. The test of PIC-D26 (h) measured 0.78–1.50 ahead of a shock at M_A ≈ 3.
- **Behind the shock the k-mode is larger than a frozen-in field gives.** Over the run B_x grows by 2.48 and B_z by 3.75 across the shock (frozen-in: 1 and 3.16). The shock alone, in the null run, has a k-mode of 0.62 (B_x) and 0.64 (B_z) of the wave on average and 1.23 and 2.14 at the end, against 2.55 and 5.23 with the wave. So the downstream mode at this k is not only the transmitted wave. One run each; separating the shock's own ripples from the processed wave needs more runs (seeds, amplitudes) or a phase analysis, which this example does not do.
- **Electron spectrum.** Downstream the electrons follow a Maxwellian of the measured T_e = 4.6e-3 m_e c² (4.6 times the upstream temperature) within 3 % up to 5 T_e. Above that the spectrum lies above it: 2.0 times at 8 T_e, 3.5 at 10 T_e and 50 at 15 T_e. The null run has the same excess (2.0, 3.6 and 81 times its own Maxwellian of T_e = 3.8e-3), so this run does not show a tail that the wave adds. T_e/T_i = 0.19 downstream.
- **Restart** (PIC-D25), on ef4fc4f. A copy of the run continued from the checkpoint of step 8000 agrees with the uninterrupted run to round-off for the first 200 steps (magnetic energy ≤ 1.2e-16 relative). The difference then grows, as round-off does in a shock run: 8.8e-14 at step 9000 and 1.2e-5 at step 13000. The shock speed agrees to 8 digits, and the downstream density and the B_z mode behind the shock to 5 digits.
- **Cost**: 10.1 ms per step with the wave (9.0 ms without), 131 s for the 13000 steps on 2 GPUs (ef4fc4f). On c571302 the same runs took 9.3 and 7.4 ms: the corrections of the review (gather and outputs across the face, the conductor) evaluate the incident field about ten times per step, about 2 ms whatever the number of particles (plan, Results). The particle number doubles to 1.25e7; the arrays reach 78 % of their capacity.
