# Bow shock of a magnetised planet in the solar wind

A magnetised electron–ion wind flows along +x through an open 3-D box and meets an absorbing sphere with a magnetic dipole at its centre. The dipole field stops the wind at the magnetopause; the wind is super-magnetosonic, so a bow shock stands upstream of it, with the sheath in between. The interplanetary field (IMF) is along +z and the dipole moment along −z, so the two fields are parallel on the day side; `south.yaml` and `oblique.yaml` turn the IMF.

- **Theory** (`theory.py`):
  - Magnetopause: the ram pressure of the wind equals the magnetic pressure of the compressed dipole field, (2 B_d)² / 2 = n m_i v², B_d = M / R³ (Chapman–Ferraro; Spreiter et al. 1966). The planner takes the distance R_mp and sets the moment M.
  - Bow shock on the subsolar line: Rankine–Hugoniot compression X of a perpendicular MHD shock (γ = 5/3), and the stand-off of gas-dynamic models, (R_bs − R_mp) / R_mp = 1.1 / X (Spreiter et al. 1966).
  - Oblique IMF: the compression of the fast MHD shock at the angle θ_Bn between the field and the shock normal (the form of Priest 1982, sec. 5.4); θ_Bn = 90° gives the perpendicular shock.
- **Known results** (not analytic):
  - a shock needs a magnetosphere larger than the ion scales, Dp = R_mp / d_i ≳ 1, and is Earth-like from Dp ≈ 20 (Omidi et al. 2004);
  - the PIC model of Mercury of Lavorenti et al. (2022) runs at Dp ≈ 8. The gate configurations here use Dp = 5, `dp10.yaml` Dp = 10.

| config | wind | grid | use |
|---|---|---|---|
| `gate.yaml` (gate) | m_i/m_e = 16, M_A = 4, v_sw = 0.25 c | 240 × 320 × 320, R_mp = 40 cells, 5760 steps | magnetopause, shock, conservation |
| `south.yaml` (gate) | the gate wind, IMF along −z (antiparallel to the dipole field on the day side) | as the gate | field reversal at the magnetopause |
| `oblique.yaml` (gate) | the gate wind, IMF along (1, 0, 1) / √2 (θ_Bn = 45° on the subsolar line) | as the gate | quasi-parallel and quasi-perpendicular flank |
| `dp10.yaml` | the gate wind, R_mp = 10 d_i | 480 × 640 × 640, R_mp = 80 cells, 11520 steps | larger magnetosphere |
| `full.yaml` | m_i/m_e = 100, M_A = 5, v_sw = 0.125 c | 600³, R_mp = 100 cells, 28800 steps | production |

All: ω_pe/Ω_ce = 4, v_te = 0.1 c (Debye length 0.2 Δx), T_i = T_e (β = 0.64), 8 particles per cell and species, three crossing times L_x / v_sw.

```bash
python -m pypic.launch -n 8 examples/bow_shock/run.py --config examples/bow_shock/configs/gate.yaml
python examples/bow_shock/analyze.py output/bow_shock/gate --movie
```

- **Set-up** (explicit scheme; engine features PIC-D12 to PIC-D16 in `docs/design/engine.md`):
  - The IMF, the motional field E = −v × B and the dipole are external fields: the grid holds the difference to them, so the absorbing layers of the open faces relax the fields to the wind values. Outputs hold the total fields.
  - All six faces are open. The wind enters through the upstream face (electrons and ions in pairs at the same place) and through the four side faces; the downstream face only removes particles.
  - The sphere (radius R_mp / 1.5) removes the particles that enter it and keeps their charge: an insulating body at floating potential. The box starts filled with the wind outside the sphere; the dipole is there from t = 0.
  - A compensated binomial current filter (4 passes) and the Marder zone at the open faces.
  - The planner refuses a body so small that Ω_ce Δt > 2 at its surface, and checks the particle ids against 2⁵³.
- **Gate** (definitions in `docs/plans/2026-11-10-bow-shock.md`, fixed before the first run): Gauss's law ≤ 1e-10; upstream density and u_x of both species within 5 % of the wind during the last crossing; ion density compression ≥ 2 in the sheath on the subsolar line; magnetopause (where the density falls back to the wind density) within 25 % of R_mp; shock upstream of the magnetopause. `south`: B_z reverses once between the shock and the body, at a distance x_rev from the centre with |x_rev / R_mp − 1| ≤ 0.25. `oblique`: the upstream state is taken at z − z_c ≥ 1.5 R_mp, the side the shocked field lines do not reach; its wind passed the null test without a planet first (`docs/plans/2026-11-11-more-examples.md`).
- **Reported**: the measured compression and shock distance against the MHD shock and Spreiter et al. (1966); the magnetopause by magnetic pressure; the sheath thickness; the charge of the body and the particles it absorbs per step; the largest fill ratio of the particle arrays.

## Result of the gate configuration (8 GPUs, four runs with the same seed)

| quantity | measured | reference |
|---|---|---|
| magnetopause distance | 0.92 to 0.93 R_mp | 1 (pressure balance) |
| shock distance | 1.38 R_mp | 1.37 (Spreiter et al. 1966, with the MHD compression) |
| largest density on the subsolar line | 3.8 n₀ | 2.97 (perpendicular MHD shock) |
| sheath thickness / magnetopause distance | 0.49 to 0.50 | 0.37 (Spreiter et al. 1966) |
| upstream density and u_x, last crossing | within 0.11 % of the wind | gate 5 % |
| Gauss's law | 4e-14 | gate 1e-10 |

- **Runs with the same seed give different third digits.** The deposit adds the particles' contributions in an order that differs from run to run, so two runs with the same seed differ in the last bit after the first step, and the flow amplifies the difference: 1e-6 of the field energy after 1440 steps in a smaller run, and the ranges of the table after the 5760 steps of this one (ANA-D3; roadmap R10). The gate criteria are far from these differences. The four runs are two before the commit of the example, the gate record (`bench/results/`) and the run of the wiki page; the plan lists each.
- **The density maximum is not the mean sheath density.** It lies above the MHD compression. Supercritical perpendicular shocks are known to overshoot behind the ramp (Leroy et al. 1982), and the sheath here is 9 d_e ≈ 2 d_i thick, too thin for a plateau behind an overshoot; this run does not separate the two.
- **The sheath is thicker than the gas-dynamic value** because the magnetopause lies 7 to 8 % inside R_mp while the shock is at the predicted distance.
- **The body charges positively** (the absorbed charge at the end is +730 to +749 in the four runs, in units of e n₀ d_e³; it is the difference of 2.4e6 absorbed ions and as many electrons) and the rate of absorption falls from 1480 to 420 particles per step between the first and the third crossing; it has not become constant by the end of the run. The particle arrays reach 57 % of their capacity (capacity factor 2).
- **Energy** is monitored only (explicit scheme, energy of the difference fields without the work of the external field): the total changes by +2.5 % over the run.
- **Cost**: 46 ms per step (0.85 ns per particle and step, times 8 GPUs), 264 s; 10.1 GB of device memory per GPU at the most. The planner's figure (11.0 GB) is for full particle arrays; for `full.yaml` it is 86 GB per GPU on 8 GPUs.

## Southward and oblique IMF (8 GPUs, one run each)

| quantity | `south` | `oblique` | northward gate |
|---|---|---|---|
| magnetopause distance (density) | 0.91 R_mp | 0.96 R_mp | 0.92 to 0.93 R_mp |
| shock distance | 1.18 R_mp | 1.25 R_mp | 1.38 R_mp |
| largest density on the subsolar line | 2.48 n₀ | 3.26 n₀ | 3.8 n₀ |
| MHD compression | 2.97 | 3.19 (θ_Bn = 45°) | 2.97 |
| sheath thickness / magnetopause distance | 0.30 | 0.30 | 0.49 to 0.50 |
| shock on the flank (line x = x_c, plane y = y_c), +z and −z | 2.24 and 2.30 R_mp | 2.38 and 1.57 R_mp | — |
| upstream density and u_x, last crossing | within 0.13 % | within 0.12 % | within 0.11 % |
| Gauss's law | 6e-14 | 3e-14 | 4e-14 |

- **Southward IMF: one reversal of B_z, at 1.005 R_mp**, 0.09 R_mp (1.9 d_e) upstream of the magnetopause by density: the density falls back to the wind value inside the field reversal. The shock stands 0.20 R_mp closer to the planet than with the northward field and the sheath is 40 % thinner. The run does not show what removes the flux (the reconnection rate at the magnetopause was not measured).
- **Oblique IMF: the flank is not symmetric.** On the side z > z_c (quasi-perpendicular) the shock is at 2.38 R_mp, on the side z < z_c (quasi-parallel) at 1.57 R_mp. The moments ahead of the quasi-parallel side deviate from the wind by 0.09 %: no foreshock is seen in density and mean velocity at this size and length (22.5 / Ω_ci). The compression is 2.4 % above the oblique MHD jump.
- **Cost**: as the gate (46 to 47 ms per step, 264 to 269 s).

## Result at R_mp = 10 d_i (8 GPUs, one run)

`dp10.yaml`: the gate wind with the magnetopause at 10 d_i = 80 cells; 480 × 640 × 640 cells, 3.1e9 particles at the start, 11520 steps (45 / Ω_ci). One recorded run (`bench/results/2026-11-11-4d55109b5-example-bow_shock-dp10-run-4b073ac2.json`); the criteria of the gate are evaluated and all hold.

| quantity | measured | reference | gate configuration (5 d_i) |
|---|---|---|---|
| magnetopause distance | 0.94 R_mp | 1 (pressure balance) | 0.92 to 0.93 |
| shock distance | 1.44 R_mp | 1.37 (Spreiter et al. 1966, with the MHD compression) | 1.38 |
| largest density on the subsolar line | 3.8 n₀ | 2.97 (perpendicular MHD shock) | 3.8 |
| sheath thickness / magnetopause distance | 0.53 | 0.37 (Spreiter et al. 1966) | 0.49 to 0.50 |
| shock on the flank, +z and −z | 2.35 and 2.34 R_mp | — | — |
| upstream density and u_x, last crossing | within 0.054 % of the wind | gate 5 % | 0.11 % |
| Gauss's law | 5e-14 | gate 1e-10 | 4e-14 |

- **The sheath is 20 d_e = 5 d_i thick** (gate: 9 d_e ≈ 2 d_i), and the largest density is the same 3.8 n₀: the maximum does not come down to the MHD compression when the sheath is 2.2 times thicker in d_i. Whether a plateau follows the maximum was not analysed; one run.
- **The shock stands 5 % farther out than the gas-dynamic value** (the gate is within 1 %), with the magnetopause at 0.94 R_mp.
- **Not stationary at the end**: the particle number grows by 10.5 % over the run and still by 0.4 % over the last 2400 steps; the body holds the charge +1670 and absorbs 5750, 2290 and 930 particles per step in the three crossings.
- **Fill**: the particle arrays reach 57 % of their capacity (capacity factor 2), with no reallocation.
- **Energy** (monitored): the total changes by +2.3 %; the budget with the outflow and the absorbed particles by +3.3 %.
- **Cost**: 1.0 hour; 313 ms per step at the end (293 at the start), 0.72 ns per particle and step, times 8 GPUs. Device memory: 78.7 GB on the fullest GPU at the most, 0.1 % above the planner's figure for full particle arrays (78.6 GB) with the arrays 57 % full — the planner's figure is not a bound (the value includes the blocks held by the memory pool).

## Result of the production configuration (8 GPUs, one run)

`full.yaml`: 600³ cells, 3.4e9 particles at the start, 28800 steps (18 / Ω_ci). One recorded run (`bench/results/2026-11-11-9ba2a6518-example-bow_shock-full-run-b3d1c288.json`); the numbers are reported, not gated.

| quantity | measured | reference |
|---|---|---|
| magnetopause distance | 0.94 R_mp | 1 (pressure balance) |
| shock distance | 1.32 R_mp | 1.34 (Spreiter et al. 1966, with the MHD compression) |
| largest density on the subsolar line | 4.3 n₀ | 3.25 (perpendicular MHD shock) |
| sheath thickness / magnetopause distance | 0.40 | 0.34 (Spreiter et al. 1966) |
| upstream density and u_x, last crossing | within 0.17 % of the wind | — |
| Gauss's law | 8e-14 | 1e-10 |

- **Not stationary at the end**: the particle number grows by 9.2 % over the run and still by 0.3 % over the last 2400 steps; the body holds the charge +7120 and absorbs 4240, 1460 and 840 particles per step in the three crossings. The sheath is 19 d_e = 1.9 d_i thick, about as thin in d_i as in the gate configuration; the note on the density maximum above holds here too.
- **The shock distance depends on its definition** (where the density reaches 2 n₀): the density rises over 0.4 R_mp = 2 d_i, exceeding 1.1 n₀ from 1.68 R_mp and 1.5 n₀ from 1.50 R_mp. The agreement with Spreiter et al. (1966) is for this definition.
- **Fill**: the particle arrays reach 62 % of their capacity (capacity factor 2), with no reallocation.
- **Energy** (monitored): the total changes by −3.2 %; the budget with the outflow and the absorbed particles by +0.4 %.
- **Cost**: 2.85 hours; 357 ms per step over the run (318 at the start, 375 at the end), 0.76 ns per particle and step, times 8 GPUs; 79.8 GB of device memory per GPU at the most (planner: 86 GB).
