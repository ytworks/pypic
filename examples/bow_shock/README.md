# Bow shock of a magnetised planet in the solar wind

A magnetised electron–ion wind flows along +x through an open 3-D box and meets an absorbing sphere with a magnetic dipole at its centre. The dipole field stops the wind at the magnetopause; the wind is super-magnetosonic, so a bow shock stands upstream of it, with the sheath in between. The interplanetary field (IMF) is along +z and the dipole moment along −z, so the two fields are parallel on the day side.

- **Theory** (`theory.py`):
  - Magnetopause: the ram pressure of the wind equals the magnetic pressure of the compressed dipole field, (2 B_d)² / 2 = n m_i v², B_d = M / R³ (Chapman–Ferraro; Spreiter et al. 1966). The planner takes the distance R_mp and sets the moment M.
  - Bow shock on the subsolar line: Rankine–Hugoniot compression X of a perpendicular MHD shock (γ = 5/3), and the stand-off of gas-dynamic models, (R_bs − R_mp) / R_mp = 1.1 / X (Spreiter et al. 1966).
- **Known results** (not analytic):
  - a shock needs a magnetosphere larger than the ion scales, Dp = R_mp / d_i ≳ 1, and is Earth-like from Dp ≈ 20 (Omidi et al. 2004);
  - the PIC model of Mercury of Lavorenti et al. (2022) runs at Dp ≈ 8. Both configurations here use Dp = 5.

| config | wind | grid | use |
|---|---|---|---|
| `gate.yaml` (gate) | m_i/m_e = 16, M_A = 4, v_sw = 0.25 c | 240 × 320 × 320, R_mp = 40 cells, 5760 steps | magnetopause, shock, conservation |
| `full.yaml` | m_i/m_e = 100, M_A = 5, v_sw = 0.125 c | 600³, R_mp = 100 cells, 28800 steps | production |

Both: ω_pe/Ω_ce = 4, v_te = 0.1 c (Debye length 0.2 Δx), T_i = T_e (β = 0.64), 8 particles per cell and species, three crossing times L_x / v_sw.

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
- **Gate** (definitions in `docs/plans/2026-11-10-bow-shock.md`, fixed before the first run): Gauss's law ≤ 1e-10; upstream density and u_x of both species within 5 % of the wind during the last crossing; ion density compression ≥ 2 in the sheath on the subsolar line; magnetopause (where the density falls back to the wind density) within 25 % of R_mp; shock upstream of the magnetopause.
- **Reported**: the measured compression and shock distance against the MHD shock and Spreiter et al. (1966); the magnetopause by magnetic pressure; the sheath thickness; the charge of the body and the particles it absorbs per step; the largest fill ratio of the particle arrays.

## Result of the gate configuration (8 GPUs, two runs with the same seed)

| quantity | measured | reference |
|---|---|---|
| magnetopause distance | 0.92 to 0.93 R_mp | 1 (pressure balance) |
| shock distance | 1.38 R_mp | 1.37 (Spreiter et al. 1966, with the MHD compression) |
| largest density on the subsolar line | 3.8 n₀ | 2.97 (perpendicular MHD shock) |
| sheath thickness / magnetopause distance | 0.49 to 0.50 | 0.37 (Spreiter et al. 1966) |
| upstream density and u_x, last crossing | within 0.11 % of the wind | gate 5 % |
| Gauss's law | 4e-14 | gate 1e-10 |

- **Two runs give different third digits.** The deposit adds the particles' contributions in an order that differs from run to run, so two runs with the same seed differ in the last bit after the first step, and the flow amplifies the difference: 1e-6 of the field energy after 1440 steps in a smaller run, and the ranges of the table after the 5760 steps of this one (ANA-D3; roadmap R10). The gate criteria are far from these differences.
- **The density maximum is not the mean sheath density.** It lies above the MHD compression. Supercritical perpendicular shocks are known to overshoot behind the ramp (Leroy et al. 1982), and the sheath here is 9 d_e ≈ 2 d_i thick, too thin for a plateau behind an overshoot; this run does not separate the two.
- **The sheath is thicker than the gas-dynamic value** because the magnetopause lies 7 to 8 % inside R_mp while the shock is at the predicted distance.
- **The body charges positively** (the absorbed charge is +742 at the end, in units of e n₀ d_e³) and the rate of absorption falls from 1480 to 420 particles per step between the first and the third crossing; it has not become constant by the end of the run. The particle arrays reach 57 % of their capacity (capacity factor 2).
- **Energy** is monitored only (explicit scheme, energy of the difference fields without the work of the external field): the total changes by +2.5 % over the run.
- **Cost**: 46 ms per step (0.85 ns per particle and step, times 8 GPUs), 264 s; 10.1 GB of device memory per GPU at the most. The planner's figure (11.0 GB) is for full particle arrays; for `full.yaml` it is 86 GB per GPU on 8 GPUs.
