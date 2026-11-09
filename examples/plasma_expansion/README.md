# Plasma expansion into vacuum: ion acceleration (Mora 2003)

A slab of hot electrons and cold ions expands into vacuum. The charge-separation sheath at the ion front accelerates the ions; this is the mechanism of target-normal sheath acceleration (TNSA).

**Isothermal self-similar solution** (`theory.py`), with τ = ω_pi t / √(2e):

- maximum energy: E_max = 2ZT_e [ln(τ + √(τ² + 1))]²
- front position: x_f = 2√(2e) λ_D0 [τ ln(τ + √(τ² + 1)) − √(τ² + 1) + 1]
- front field: E_f = 2E₀ / √(2e + ω_pi² t²)
- ion spectrum: dN/dE = n_i0 c_s t / √(2EZT_e) exp(−√(2E/ZT_e))

| config | set-up |
|---|---|
| `isothermal_hydrogen.yaml` (gate) | T_e = 5.1 keV, m_i/m_e = 1836, slab of 200 λ_D, to τ = 10 |
| `isothermal.yaml` | the same with m_i/m_e = 100 (fast; the model needs v_ion ≪ v_te, and here v_front ≈ 0.6 v_te) |
| `thin_foil.yaml` | slab of 10 λ_D: the electrons cool (adiabatic expansion) and E_max falls below the isothermal law |

```bash
python -m pypic.launch -n 8 examples/plasma_expansion/run.py --config examples/plasma_expansion/configs/isothermal_hydrogen.yaml
python examples/plasma_expansion/analyze.py output/plasma_expansion/isothermal_hydrogen --movie
```

- **Set-up**: a step density profile shared by electrons and ions (PIC-D9), so the plasma is exactly neutral. Both walls are reflecting.
- **Gate**:
  - the exact self-similar ion velocity v = c_s + x/t in the quasi-neutral rarefaction (−0.8 < x/c_s t < 3) within 10 % (median);
  - E_max(t) within 20 % (median over 2 ≤ τ ≤ τ_end).
- **Measured** (2 GPUs, proton mass):
  - velocity profile −8 %, density profile 7 %;
  - E_max −15 % at τ = 10 (13.5 ZT_e, against 18.0 for the isothermal model);
  - front position −8 %;
  - the ion spectrum follows the self-similar dN/dE from 0.01 to a few ZT_e.
- **Why the PIC result is lower**:
  - Collisionless electrons are not the model's ideal isothermal (Boltzmann, infinite heat bath) electrons: the expanding region is cooler than the slab.
  - All the deviations are consistent with a ~15 % lower effective T_e there.
  - The offset is unchanged by Δx = λ_D/2, 2× particles or a 2.5× longer vacuum region, so it is not numerical. It grows at m_i/m_e = 100 (E_max −18 %), as expected when v_front approaches v_te.
