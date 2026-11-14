# N7: the diffusion measurement on synthetic fields with a known D

Kind: null test N7 of the README (method). It tests the measurement code `diffusion.py`, not a plasma; no PIC
run is involved.

| | |
|---|---|
| Date | 2026-11-15 (04:30 JST) |
| Commit | 53cd3a902844617a6b136f77f6dedbf64dcfda40 (clean tree) |
| Command | `python n7_synthetic.py --json <file>` (defaults: `--n 20000 --seed 7`) |
| Seed | 7 (each case draws from `default_rng([7, theta, case])`) |
| GPUs | none (NumPy, one CPU process, 451 s) |
| Output | `2026-11-15-n7-synthetic.json` (the whole output, 8 kB) |

## Fixed before the run (README, N7)

The measured D is within 10 % of the known one for D between 0.1 and 10 D_B, and it is taken from an interval
of a factor 3 in time whose log-slope is within 1 ± 0.1.

## Design (chosen for this test; not in the README)

- Field with a known D: a uniform B0 at the angle θ to x plus a magnetic field without memory (drawn anew for
  every particle, component and step). The direction of a particle then diffuses on the sphere at the rate ν, and
  D / D_B = (2/3) [cos²θ Ω/ν + sin²θ (ν/Ω) / (1 + (ν/Ω)²)] in the conventions of `diffusion.py`
  (D = ⟨Δx²⟩/t, D_B = v r_L).
- p = m_i v_pt of V1 (u = 10), B0 = 0.0667; ν/Ω = 1/15, 2/9, 2/3, 20/9, 20/3; θ = 0° and 30°;
  dt = min(0.05, 0.01 Ω/ν)/Ω; the run lasts 100/ν; 20,000 particles (statistical error of ⟨Δx²⟩: 1 %).
- The value taken is ⟨Δx²⟩/t at the end of the last accepted interval. The slope of a straight line through
  ⟨Δx²⟩(t) over the same interval is given beside it.

## Numbers

| θ | ν/Ω | exact D/D_B | ⟨Δx²⟩/t | deviation | line slope | log-slope of the interval |
|---|---|---|---|---|---|---|
| 0° | 1/15 | 10.000 | 10.090 | +0.9 % | 10.217 | 1.025 |
| 0° | 2/9 | 3.000 | 2.973 | −0.9 % | 3.004 | 1.021 |
| 0° | 2/3 | 1.000 | 1.001 | +0.1 % | 1.024 | 1.038 |
| 0° | 20/9 | 0.300 | 0.295 | −1.6 % | 0.295 | 1.004 |
| 0° | 20/3 | 0.100 | 0.101 | +1.1 % | 0.103 | 1.037 |
| 30° | 1/15 | 7.511 | 7.430 | −1.1 % | 7.520 | 1.023 |
| 30° | 2/9 | 2.285 | 2.311 | +1.1 % | 2.360 | 1.036 |
| 30° | 2/3 | 0.827 | 0.827 | +0.1 % | 0.838 | 1.026 |
| 30° | 20/9 | 0.287 | 0.290 | +1.1 % | 0.294 | 1.021 |
| 30° | 20/3 | 0.0994 | 0.0978 | −1.7 % | 0.0992 | 1.017 |

- N7 passes: the nine rows inside 0.1–10 D_B deviate by at most 1.6 % (the criterion is 10 %); every value comes
  from an accepted interval. The tenth row (exact 0.0994) is just outside the range and deviates by 1.7 %.
- The line slope deviates by −1.6 % to +3.3 % from the exact value.
- |u| changes by at most 3e-14 of u (the random field is magnetic). ⟨Δx⟩ / rms(Δx) is at most 0.009.
- Ballistic control (the same uniform field at 30°, no random part): log-slope 2.000 over the last factor 3;
  the rule returns no interval, as it must.
- Fields that vary along x (no threshold fixed):
  - interpolation of a field linear in x, cut and mirrored: largest error 1.8e-15 for a component at the nodes
    and for one at the half cells;
  - constants of the motion in a static (0, B_y(x), B_z(x)) over five gyrations, largest change over u:

    | dt | 3.2 | 1.6 | 0.8 | 0.4 | 0.2 | 0.1 | 0.4, components read half a cell off |
    |---|---|---|---|---|---|---|---|
    | u_y − a_z | 4.1e-4 | 8.6e-5 | 2.0e-5 | 4.4e-6 | 5.5e-7 | 3.4e-7 | 6.7e-4 |
    | u_z + a_y | 9.9e-5 | 2.8e-5 | 5.7e-6 | 1.0e-6 | 1.8e-7 | 1.2e-7 | 3.3e-4 |

    The change falls by a factor 4 to 8 per halving of dt down to dt = 0.2, and by 1.6 for the last halving.
    With the stagger wrong at dt = 0.4 it is 150 and 330 times the value with the right stagger.

## Limits

- N7 covers the pusher, the bookkeeping of the displacement, the conventions (no factor 1/2 in D, no 1/3 in D_B)
  and the rule of the diffusive interval. The scattering here has no memory; in a field snapshot it comes from
  waves with a correlation length, where an interval that passes the rule may still be short of the limit of
  Eq. 9. The two estimators side by side show that case (they differ by about the correlation time over t).
- It does not test the approximation of a frozen snapshot (N9), the mirroring and the length of the domain
  (N10), or the electric field of the waves (the synthetic fields have none).
- One momentum (u = 10), one field strength, one seed; statistical error 1 %.
