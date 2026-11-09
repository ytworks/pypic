**English** · [日本語](Plasma-Waves-ja) · source: [`examples/plasma_waves`](https://github.com/ytworks/pypic/tree/main/examples/plasma_waves)

# Thermal plasma waves: dispersion relations from noise

A uniform Maxwellian plasma at rest is never quiet. Its thermal fluctuations excite every normal mode
of the plasma, and a space–time Fourier transform of the fields, |F(k, ω)|², shows the dispersion
branches ω(k) directly. No perturbation is imposed: the simulation is left alone, and the measured
spectrum is laid over linear theory. The gated branches agree to between
{{ r.perpendicular.branches.Ez: O mode.median_rel_error | .2% }} (O mode) and
{{ r.unmagnetized.branches.Ex: Langmuir (kinetic).median_rel_error | .1% }} (Langmuir waves), median over k.

## Physics

| config | magnetic field | branches compared (theory in `pypic.theory.waves`) |
|---|---|---|
| `unmagnetized` | none | E_x: Langmuir waves, the kinetic root of the Maxwellian dielectric function (Bohm–Gross ω² = ω_pe² + 3k²v_t² for kλ_D ≪ 1). E_y, E_z: light waves ω² = ω_pe² + c²k², with the dispersion of the Yee grid and leapfrog time step |
| `perpendicular` | B₀ ∥ z, k ⊥ B₀, ω_ce = {{ c.perpendicular.params.wce }} ω_pe | E_z: ordinary (O) mode. E_y: extraordinary (X) mode, upper and lower branches. E_x: electron Bernstein waves between the cyclotron harmonics |
| `parallel` | B₀ ∥ x, k ∥ B₀, ω_ce = {{ c.parallel.params.wce }} ω_pe | E_y, E_z: whistler, R and L waves. E_x: Langmuir waves |
| `ion_acoustic` | none, m_i/m_e = {{ c.ion_acoustic.params.mass_ratio }}, T_e/T_i = 50 | E_x: ion-acoustic waves, the kinetic root with electron and ion Landau damping |

The cold magnetised branches follow from the Stix dielectric tensor; for example the X mode solves

```math
n^2 = \frac{c^2k^2}{\omega^2} = \frac{(\omega^2-\omega_{pe}^2)^2 - \omega^2\omega_{ce}^2}
{\omega^2(\omega^2 - \omega_{uh}^2)}, \qquad \omega_{uh}^2 = \omega_{pe}^2 + \omega_{ce}^2 ,
```

and the electron Bernstein waves are the roots of the hot perpendicular (electrostatic) dielectric
function, with ρ_e = v_te/ω_ce the thermal Larmor radius:

```math
1 = \frac{2\omega_{pe}^2}{\omega_{ce}^2}\sum_{n=1}^{\infty}
\frac{n^2\, e^{-\lambda} I_n(\lambda)/\lambda}{\omega^2/\omega_{ce}^2 - n^2},
\qquad \lambda = k^2\rho_e^2 .
```

As λ → 0 the first band ends at the upper-hybrid frequency ω_uh, the others at the harmonics nω_ce.

## Experimental conditions

Normalised units: c = ε₀ = 1, lengths in c/ω_pe, times in 1/ω_pe, θ = T/(m c²). The box is quasi-1-D
along k ∥ x with Δx = λ_De; the lineouts are averaged over the 8 × 8 transverse cells. The frequency
resolution is Δω = 2π/t_end = {{ c.unmagnetized.planner.frequency_resolution | .4f }} ω_pe (for t_end = 300).

{{ table:conditions }}

<details><summary>Planner parameters</summary>

{{ table:params }}

</details>

## Animation

![plasma waves from thermal noise](images/plasma_waves/unmagnetized.gif)

`unmagnetized`. **Top row**: E_x(x, t) and E_y(x, t), drawn as time advances. The longitudinal field
E_x is a jumble of slow Langmuir wave packets. The transverse field E_y shows light waves crossing the
box in both directions at almost c (faint streaks at ±45° in the x–t plane). **Bottom left**: the
instantaneous lineouts. **Bottom right**: the energies. The kinetic energy stays constant. The fields
start from zero, rise within a few plasma periods to the level of the thermal fluctuations and then creep
up slowly; they stay more than two orders of magnitude below the kinetic energy.

## Comparison with linear theory

### ω–k spectra with the theory curves

The colour is log₁₀ of the spectral power; the dashed and dotted lines are linear theory. Every ridge
of power sits on a theory curve.

![unmagnetized](images/plasma_waves/unmagnetized-omega_k.png)

`unmagnetized`: Langmuir waves in E_x, light waves in E_y and E_z.

![perpendicular](images/plasma_waves/perpendicular-omega_k.png)

`perpendicular`: O mode (E_z), upper and lower X branches (E_y), and the Bernstein bands between the
harmonics nω_ce (E_x, out to kρ_e ≈ 3).

![parallel](images/plasma_waves/parallel-omega_k.png)

`parallel`: whistler, R and L waves (E_y, E_z), Langmuir waves (E_x).

![ion acoustic](images/plasma_waves/ion_acoustic-omega_k.png)

`ion_acoustic`: the ion-acoustic branch at low frequency.

### Error of every branch

![median error of every branch](images/plasma_waves/branch_errors.png)

At every k the peak of the spectrum is searched within ±15 % of the theoretical frequency. A branch
passes when it is found at half of the k points or more, and at least 90 % of the points found lie
within max(3 %, 2Δω/ω) of theory. Filled markers are gated branches, hollow ones are reported only.

| config | branch | found / k points | median error | max error | within tolerance | gated |
|---|---|---|---|---|---|---|
| `unmagnetized` | E_x: Langmuir (kinetic) | {{ r.unmagnetized.branches.Ex: Langmuir (kinetic).found }} / {{ r.unmagnetized.branches.Ex: Langmuir (kinetic).points }} | {{ r.unmagnetized.branches.Ex: Langmuir (kinetic).median_rel_error | .2% }} | {{ r.unmagnetized.branches.Ex: Langmuir (kinetic).max_rel_error | .2% }} | {{ r.unmagnetized.branches.Ex: Langmuir (kinetic).fraction_within_tolerance | .0% }} | yes |
| `unmagnetized` | E_y: EM (Yee) | {{ r.unmagnetized.branches.Ey: EM (Yee).found }} / {{ r.unmagnetized.branches.Ey: EM (Yee).points }} | {{ r.unmagnetized.branches.Ey: EM (Yee).median_rel_error | .2% }} | {{ r.unmagnetized.branches.Ey: EM (Yee).max_rel_error | .2% }} | {{ r.unmagnetized.branches.Ey: EM (Yee).fraction_within_tolerance | .0% }} | yes |
| `unmagnetized` | E_z: EM (Yee) | {{ r.unmagnetized.branches.Ez: EM (Yee).found }} / {{ r.unmagnetized.branches.Ez: EM (Yee).points }} | {{ r.unmagnetized.branches.Ez: EM (Yee).median_rel_error | .2% }} | {{ r.unmagnetized.branches.Ez: EM (Yee).max_rel_error | .2% }} | {{ r.unmagnetized.branches.Ez: EM (Yee).fraction_within_tolerance | .0% }} | yes |
| `perpendicular` | E_z: O mode | {{ r.perpendicular.branches.Ez: O mode.found }} / {{ r.perpendicular.branches.Ez: O mode.points }} | {{ r.perpendicular.branches.Ez: O mode.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ez: O mode.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ez: O mode.fraction_within_tolerance | .0% }} | yes |
| `perpendicular` | E_y: X upper | {{ r.perpendicular.branches.Ey: X upper.found }} / {{ r.perpendicular.branches.Ey: X upper.points }} | {{ r.perpendicular.branches.Ey: X upper.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ey: X upper.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ey: X upper.fraction_within_tolerance | .0% }} | yes |
| `perpendicular` | E_y: X lower | {{ r.perpendicular.branches.Ey: X lower.found }} / {{ r.perpendicular.branches.Ey: X lower.points }} | {{ r.perpendicular.branches.Ey: X lower.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ey: X lower.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ey: X lower.fraction_within_tolerance | .0% }} | reported |
| `perpendicular` | E_x: Bernstein 1 | {{ r.perpendicular.branches.Ex: Bernstein 1.found }} / {{ r.perpendicular.branches.Ex: Bernstein 1.points }} | {{ r.perpendicular.branches.Ex: Bernstein 1.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 1.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 1.fraction_within_tolerance | .0% }} | yes |
| `perpendicular` | E_x: Bernstein 2 | {{ r.perpendicular.branches.Ex: Bernstein 2.found }} / {{ r.perpendicular.branches.Ex: Bernstein 2.points }} | {{ r.perpendicular.branches.Ex: Bernstein 2.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 2.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 2.fraction_within_tolerance | .0% }} | yes |
| `perpendicular` | E_x: Bernstein 3 | {{ r.perpendicular.branches.Ex: Bernstein 3.found }} / {{ r.perpendicular.branches.Ex: Bernstein 3.points }} | {{ r.perpendicular.branches.Ex: Bernstein 3.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 3.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 3.fraction_within_tolerance | .0% }} | yes |
| `perpendicular` | E_x: Bernstein 4 | {{ r.perpendicular.branches.Ex: Bernstein 4.found }} / {{ r.perpendicular.branches.Ex: Bernstein 4.points }} | {{ r.perpendicular.branches.Ex: Bernstein 4.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 4.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 4.fraction_within_tolerance | .0% }} | reported |
| `parallel` | E_y: whistler | {{ r.parallel.branches.Ey: whistler.found }} / {{ r.parallel.branches.Ey: whistler.points }} | {{ r.parallel.branches.Ey: whistler.median_rel_error | .2% }} | {{ r.parallel.branches.Ey: whistler.max_rel_error | .2% }} | {{ r.parallel.branches.Ey: whistler.fraction_within_tolerance | .0% }} | — |
| `parallel` | E_y: R | {{ r.parallel.branches.Ey: R.found }} / {{ r.parallel.branches.Ey: R.points }} | {{ r.parallel.branches.Ey: R.median_rel_error | .2% }} | {{ r.parallel.branches.Ey: R.max_rel_error | .2% }} | {{ r.parallel.branches.Ey: R.fraction_within_tolerance | .0% }} | — |
| `parallel` | E_y: L | {{ r.parallel.branches.Ey: L.found }} / {{ r.parallel.branches.Ey: L.points }} | {{ r.parallel.branches.Ey: L.median_rel_error | .2% }} | {{ r.parallel.branches.Ey: L.max_rel_error | .2% }} | {{ r.parallel.branches.Ey: L.fraction_within_tolerance | .0% }} | — |
| `parallel` | E_x: Langmuir (kinetic) | {{ r.parallel.branches.Ex: Langmuir (kinetic).found }} / {{ r.parallel.branches.Ex: Langmuir (kinetic).points }} | {{ r.parallel.branches.Ex: Langmuir (kinetic).median_rel_error | .2% }} | {{ r.parallel.branches.Ex: Langmuir (kinetic).max_rel_error | .2% }} | {{ r.parallel.branches.Ex: Langmuir (kinetic).fraction_within_tolerance | .0% }} | — |
| `ion_acoustic` | E_x: ion acoustic (kinetic) | {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).found }} / {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).points }} | {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).median_rel_error | .2% }} | {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).max_rel_error | .2% }} | {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).fraction_within_tolerance | .0% }} | — |

"Gated" refers to the GPU gate runs: only `unmagnetized` and `perpendicular` are gate configurations
(`scripts/gpu_gates.py`); `parallel` and `ion_acoustic` are reported. The E_z rows of `parallel` repeat
E_y and are in the chart above.

Gauss's law holds to {{ r.unmagnetized.gauss_max | .1e }} (`unmagnetized`) and
{{ r.perpendicular.gauss_max | .1e }} (`perpendicular`); ∇·B to {{ r.perpendicular.divB_max | .1e }}.

## Notes

- **Lower X branch** (reported): it approaches the upper-hybrid resonance, where the thermal correction
  to the cold theory is no longer small, so it is shown but not gated.
- **Bernstein band 4** (reported): thermal noise excites it only weakly at these k, so the peak is found
  at fewer than half of the points.
- **Ion-acoustic waves** (reported): the branch is strongly Landau damped, so its spectral peak is broad
  compared with the frequency resolution (Δω = {{ c.ion_acoustic.planner.frequency_resolution | .4f }} ω_pe).
  The median error is {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).median_rel_error | .1% }}, and only
  {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).fraction_within_tolerance | .0% }} of the points are within
  the 3 % tolerance (90 % needed), so this configuration does not meet the analysis criterion.
- **Why θ_e = 10⁻³ in the magnetised cases**: at θ_e = 10⁻² the relativistic mass increase
  (⟨γ⟩ − 1 ≈ 1.5 θ_e) shifts the cyclotron harmonics, and the Bernstein bands come out 1.5–4 % below the
  non-relativistic theory.
- **Output cadence**: lineouts are written at the Nyquist rate of the highest frequency in the analysis
  window and time-averaged over each interval, which is an anti-aliasing filter for the high-k
  electromagnetic noise.

## Reproduce

```bash
python -m pypic.launch -n 4 examples/plasma_waves/run.py --config examples/plasma_waves/configs/unmagnetized.yaml \
    --out output/wiki/plasma_waves/unmagnetized
python examples/plasma_waves/analyze.py output/wiki/plasma_waves/unmagnetized --movie
```

## Data provenance

{{ table:provenance }}

## References

- B. D. Fried and S. D. Conte, *The Plasma Dispersion Function* (Academic Press, 1961).
- T. H. Stix, *Waves in Plasmas* (AIP, 1992).
- I. B. Bernstein, "Waves in a plasma in a magnetic field", Phys. Rev. 109, 10 (1958).
