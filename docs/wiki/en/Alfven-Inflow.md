**English** · [日本語](Alfven-Inflow-ja) · source: [`examples/alfven_inflow`](https://github.com/ytworks/pypic/tree/main/examples/alfven_inflow)

# Alfvénic wave carried by a flow into a shock

A magnetised electron–ion plasma (m_i/m_e = {{ c.gate.params.mass_ratio }}) flows along −x across the
field B₀ = B₀ ŷ. The face x+ injects it, and the plasma carries in a circularly polarised wave that
propagates along B₀: the cold two-fluid R wave at k d_i = {{ c.gate.planner.k_di | .3f }} with
|δB| = {{ c.gate.params.amplitude }} B₀. In the open box (the gate) the wave crosses the box with its
k-mode within {{ r.gate.Bx_mean | .1% }} (B_x) and {{ r.gate.Bz_mean | .1% }} (B_z) of the incident wave
on average. With a reflecting wall at x− (the piston) the flow piles up into a perpendicular shock
(Alfvén Mach number {{ r.piston.mach_alfven_shock | .1f }} in its frame), which runs into the plasma
that carries the wave.

The example exercises three engine features: the turbulent inflow face (PIC-D26), the frozen-in initial
state (PIC-D24: vector potential, shaped velocities, E from a displacement of the electrons) and restart
checkpoints (PIC-D25; the piston configuration writes them).

## Physics

**Parallel waves of a cold plasma.** Electrons and ions at rest in B₀ = B₀ ŷ respond to a wave
E e^{i(ky − ωt)} with the cold velocities −iω v = (q/m)(E + v × B₀). The current closes Maxwell's
equations on a circular vector of the (z, x) plane. Its two roots below the cyclotron resonances are the
R wave (rotating with the electrons) and the L wave (with the ions). They are the textbook relations

```math
\frac{k^2 c^2}{\omega^2} = 1 - \sum_s \frac{\omega_{ps}^2}{\omega(\omega \pm \Omega_s)} ,
\qquad \frac{\omega}{k v_A} \simeq 1 \pm \frac{k d_i}{2} \quad (k d_i \ll 1) ,
```

with + for R. In the box frame the plasma moves at −V₀ x̂; the wave's phase k y − ωt is the same there
(k ⊥ the flow), and its field is E = E_rest − v_d × B, which adds the motional field −V₀ δB_z ŷ.

**The inflow face.** The table of the face holds the wave's B, E and the velocity of each species, as
cos and sin parts of k·r − ωt + φ. Each entering particle takes the velocity of the wave at its birth
point and step; the fields enter through a total-field / scattered-field surface on the face plane.

**Piston shock.** A wall at rest stops the flow; the shock moves out at V_sh and the downstream is at
rest. In the shock frame the upstream speed is U₁ = V₀ + V_sh, so r = n₂/n₁ = U₁/V_sh. r is the
compression of a perpendicular MHD shock at M_A = U₁/v_A (γ = 5/3):

```math
(2-\gamma) r^2 + \left[\gamma(\beta+1) + (\gamma-1) M_A^2\right] r - (\gamma+1) M_A^2 = 0 ,
\qquad V_{sh} = \frac{V_0}{r - 1} .
```

The downstream temperature T_e + T_i follows from the momentum flux. A field frozen into the plasma keeps
the normal δB_x and compresses the tangential δB_z by r.

| quantity | gate | piston |
|---|---|---|
| Alfvén speed v_A | {{ c.gate.planner.v_a | .4f }} c | {{ c.piston.planner.v_a | .4f }} c |
| flow speed V₀ | {{ c.gate.planner.v0 | .4f }} c | {{ c.piston.planner.v0 | .4f }} c |
| wave frequency ω / (k v_A) (R wave) | {{ c.gate.planner.omega_over_kva | .3f }} | {{ c.piston.planner.omega_over_kva | .3f }} |
| plasma β | {{ c.gate.planner.beta | .2f }} | {{ c.piston.planner.beta | .2f }} |
| jump conditions: r, V_sh, M_A | — | {{ r.piston.compression_theory | .2f }}, {{ r.piston.shock_speed_theory | .4f }} c, {{ r.piston.mach_alfven_shock | .2f }} |

## Experimental conditions

Normalised units: c = ε₀ = 1, lengths in d_e = c/ω_pe, times in 1/ω_pe, the density of each species
n₀ = 1. Quasi-2-D: two cells along z.

- **Gate (`gate`)**: x− is an open face where the plasma leaves; the wave fills the box at t = 0 (B_x,
  B_z and E_x as profiles along y, E_z from a displacement of the electrons along z, the velocities as
  shaped perturbations). The x+ face is a `flux` face with the Marder zone.
- **Piston (`piston`)**: x− is a reflecting wall. The wave is loaded times an envelope h(x) that is zero
  within 4 d_e of the wall, rises to 1 at 12 d_e and falls to 0 over the last 9 d_e before the face:
  B from the vector potential A_z = h a_x(y), A_x = −h a_z(y) (∇·B = 0 exactly), E_x and E_z from
  displacements of the electrons. The face carries the wave in from t = 0. Restart checkpoints every
  4000 steps.
- **Null run (`piston_null`)**: the piston without the wave (a table of zero amplitude), for the shock's
  own k-mode.

{{ table:conditions }}

<details><summary>Planner parameters</summary>

{{ table:params }}

</details>

## Animation

Gate: the wave carried through the open box.

![gate](images/alfven_inflow/gate.gif)

Piston: the shock runs into the wave.

![piston](images/alfven_inflow/piston.gif)

- **Left, B_z / B₀ in the plane z = 0.**
- **Middle, ion density in the plane z = 0.**
- **Right, electron and ion density averaged over y and z.**

## Gate: the wave against the incident one

![gate](images/alfven_inflow/gate-alfven_inflow.png)

B_z at the end, the k-mode of every x row against the incident wave, and the densities.

| quantity | PIC | gate |
|---|---|---|
| B_x row error, largest (x > 10 d_e) | {{ r.gate.Bx_row_max | .1% }} | ≤ 12 % |
| B_x row error, mean (10 < x < 60 d_e) | {{ r.gate.Bx_mean | .1% }} | ≤ 4 % |
| B_z row error, largest | {{ r.gate.Bz_row_max | .1% }} | ≤ 12 % |
| B_z row error, mean | {{ r.gate.Bz_mean | .1% }} | ≤ 4 % |
| densities in bins of 2 d_e | within {{ r.gate.density_deviation | .1% }} | 10 % |
| Gauss's law, max\|∇·E − ρ\| | {{ r.gate.gauss_max | .0e }} | ≤ 1e-10 |
| max\|∇·B\| | {{ r.gate.divB_max | .0e }} | ≤ 1e-12 |
| mean amplitude ratio B_x, B_z | {{ r.gate.Bx_ratio | .3f }}, {{ r.gate.Bz_ratio | .3f }} | reported |

The criteria were fixed before the first run. The largest row errors sit next to the inflow face.

## Piston: the shock and the wave

![piston](images/alfven_inflow/piston-alfven_piston.png)

B_z at the end with the shock (dashed), the ion density in x and t with the shock position, the
k-mode of B_x and B_z ahead of, far ahead of and behind the shock, and the electron spectra. These
quantities are reported, not gated.

| quantity | with the wave | null run (no wave) | jump conditions |
|---|---|---|---|
| shock speed (fit over the second half) | {{ r.piston.shock_speed | .4f }} c | {{ r.piston_null.shock_speed | .4f }} c | {{ r.piston.shock_speed_theory | .4f }} c |
| downstream density n₂ / n₁ | {{ r.piston.compression | .2f }} | {{ r.piston_null.compression | .2f }} | {{ r.piston.compression_theory | .2f }} |
| downstream B_y / B₀ | {{ r.piston.field_compression | .2f }} | {{ r.piston_null.field_compression | .2f }} | {{ r.piston.compression_theory | .2f }} |
| downstream (T_e + T_i) / jump conditions | {{ r.piston.temperature_ratio | .2f }} | {{ r.piston_null.temperature_ratio | .2f }} | 1 |
| downstream T_e, T_i [m_e c²] | {{ r.piston.Te_down | .2e }}, {{ r.piston.Ti_down | .2e }} | {{ r.piston_null.Te_down | .2e }}, {{ r.piston_null.Ti_down | .2e }} | sum {{ r.piston.temperature_theory | .2e }} |
| B_x mode ahead / behind (units of the wave) | {{ r.piston.Bx_ahead_mean | .2f }} / {{ r.piston.Bx_behind_mean | .2f }} | {{ r.piston_null.Bx_ahead_mean | .3f }} / {{ r.piston_null.Bx_behind_mean | .2f }} | frozen-in: behind = ahead |
| B_z mode ahead / behind (units of the wave) | {{ r.piston.Bz_ahead_mean | .2f }} / {{ r.piston.Bz_behind_mean | .2f }} | {{ r.piston_null.Bz_ahead_mean | .3f }} / {{ r.piston_null.Bz_behind_mean | .2f }} | frozen-in: behind = r × ahead |
| Gauss's law, max\|∇·B\| | {{ r.piston.gauss_max | .0e }}, {{ r.piston.divB_max | .0e }} | {{ r.piston_null.gauss_max | .0e }}, {{ r.piston_null.divB_max | .0e }} | ≤ 1e-10, ≤ 1e-12 |

The mode levels are means over the outputs at which both windows hold rows: ahead from x_sh + 4 to
x_sh + 20 d_e (the foot of the reflected ions lies there), behind from 0.2 x_sh to x_sh − 4 d_e.

## Notes

- **Behind the shock the k-mode is larger than a frozen-in field gives**: B_x grows by
  {{ r.piston.Bx_transmission | .1f }} and B_z by {{ r.piston.Bz_transmission | .1f }} across the shock
  (frozen-in: 1 and {{ r.piston.compression_theory | .2f }}). The null run shows the shock's own
  fluctuations at that k: {{ r.piston_null.Bx_behind_mean | .2f }} (B_x) and
  {{ r.piston_null.Bz_behind_mean | .2f }} (B_z) in units of the wave on average and
  {{ r.piston_null.Bx_behind_end | .2f }} and {{ r.piston_null.Bz_behind_end | .2f }} at the end, against
  {{ r.piston.Bx_behind_end | .2f }} and {{ r.piston.Bz_behind_end | .2f }} with the wave. One run each; this example does not
  separate the shock's own ripples from the processed wave.
- **The shock speed and the compression differ between the two runs** (with the wave
  {{ r.piston.shock_speed_ratio | .3f }}, without {{ r.piston_null.shock_speed_ratio | .3f }} of the jump
  conditions). Both are single runs measured at one time.
- **Why the envelope vanishes at the face too.** The engine refuses a vector potential tangential to a
  non-periodic face that is not zero there (a reflecting wall would zero the normal B). The face then
  starts without the wave over the last 9 d_e; that plasma is flushed in at V₀ and reaches the shock near
  t ≈ 850, after which the shock meets the wave of the face.
- **Why the null run keeps a table.** The engine allows the motional field E₀ tangential to the wall only
  when the open face of that axis carries a fluctuation table; the null run carries one of zero
  amplitude.
- The motional field of the wave along y (−V₀ δB_z), whose charge v₀ δv / c² the particles do not carry,
  is not loaded; the plasma sets it up.

## Reproduce

```bash
python -m pypic.launch -n 2 examples/alfven_inflow/run.py \
    --config examples/alfven_inflow/configs/gate.yaml --out output/wiki/alfven_inflow/gate
python -m pypic.launch -n 2 examples/alfven_inflow/run.py --config examples/alfven_inflow/configs/piston.yaml
python -m pypic.launch -n 2 examples/alfven_inflow/run.py --config examples/alfven_inflow/configs/piston.yaml \
    --set wave=false --out output/alfven_inflow/piston_null
python examples/alfven_inflow/analyze.py output/alfven_inflow/piston --movie
```

## Data provenance

{{ table:provenance }}

## References

- T. H. Stix, *Waves in Plasmas* (AIP, New York, 1992), chapter 2: the cold-plasma R and L waves.
- E. R. Priest, *Solar Magnetohydrodynamics* (Reidel, Dordrecht, 1982): jump conditions of a
  perpendicular MHD shock.
