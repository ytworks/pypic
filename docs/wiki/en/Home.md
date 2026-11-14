**English** · [日本語](Home-ja)

# pypic examples: simulations checked against exact and published results

[pypic](https://github.com/ytworks/pypic) is a GPU-accelerated, multi-GPU, 3-D special-relativistic
electromagnetic particle-in-cell (PIC) code. It has an explicit scheme (Buneman–Boris push, Esirkepov
charge-conserving deposition, Yee fields) and an implicit scheme that conserves the total energy exactly.
Gauss's law and ∇·B hold to round-off in every run.

Each example below is a plasma problem with an analytic solution or a well-known published result. Every
page shows the experimental conditions, an animation of the physics, and the key charts that compare the
simulation with the solution. All numbers, charts and animations are generated from the runs' own output
by [`scripts/build_wiki.py`](https://github.com/ytworks/pypic/blob/main/scripts/build_wiki.py); nothing is
copied by hand.

| example | physics | compared with | result |
|---|---|---|---|
| [Two-stream and Buneman instabilities](Two-Stream-Instability) | electrostatic beam instabilities (e-i, e-p, e-i-p, implicit scheme, Buneman) | warm-fluid linear dispersion relation | growth rate {{ two_stream:r.buneman.relative_error | +.1% }} … {{ two_stream:r.ei_implicit.relative_error | +.1% }}; implicit energy error {{ two_stream:r.ei_implicit.energy_error_max | .0e }} |
| [Thermal plasma waves](Plasma-Waves) | ω–k spectra of thermal noise: Langmuir, light, O/X, Bernstein, whistler, R/L, ion-acoustic waves | kinetic and cold magnetised dispersion relations | median error {{ plasma_waves:r.perpendicular.branches.Ez: O mode.median_rel_error | .2% }} (O mode) … {{ plasma_waves:r.unmagnetized.branches.Ex: Langmuir (kinetic).median_rel_error | .1% }} (Langmuir) |
| [Langmuir oscillations](Landau-Damping) | Landau damping, nonlinear trapping, Dawson oscillation, wave breaking, relativistic oscillation | kinetic root; exact cold-sheet solutions | damping rate {{ landau_damping:r.landau.gamma_rel_error | +.1% }}; wave-breaking time {{ landau_damping:r.wavebreaking.rel_error | +.1% }} |
| [Current filamentation](Filamentation-Instability) | Weibel instability of counter-streaming beams | cold relativistic dispersion relation | growth rates {{ filamentation:r.pair.median_rel_error_m2_10 | .1% }} (median over modes) |
| [Test-particle orbits](Particle-Orbits) | gyration, E×B, hyperbolic, E > B, plane-wave figure 8 | exact relativistic orbits | error ≤ {{ particle_orbits:r.gyration.max_rel_position_error | .2% }}; second order (ratio {{ particle_orbits:r.plane_wave.convergence_ratio | .2f }}–{{ particle_orbits:r.gyration.convergence_ratio | .2f }}) |
| [Plasma expansion into vacuum](Plasma-Expansion) | ion acceleration by the expansion sheath (TNSA) | Mora (2003) self-similar solution | velocity profile {{ plasma_expansion:r.isothermal_hydrogen.velocity_profile_median_rel_error | .1% }}; E_max {{ plasma_expansion:r.isothermal_hydrogen.E_max_median_rel_error | .0% }} below (kinetic electrons) |
| [Relativistic pair shock](Relativistic-Shock) | γ₀ = 15 collisionless shock, Fermi acceleration | relativistic jump conditions (Synge gas, 2-D and 3-D) | β_sh {{ relativistic_shock:r.pair_3d.beta_sh_rel_error | +.1% }}, n₂/n₁ {{ relativistic_shock:r.pair_3d.compression_rel_error | +.1% }} (3-D) |
| [Relativistic magnetic reconnection](Magnetic-Reconnection) | σ = 10 pair reconnection, plasmoids, power-law particles | Harris equilibrium; Werner et al. (2016) | p = {{ reconnection:r.pair_sigma10_large.power_law_index | .2f }} vs {{ reconnection:r.pair_sigma10_large.werner2016_index | .2f }} |
| [Bow shock of a magnetised planet](Bow-Shock) | magnetopause and bow shock of a dipole in a magnetised wind (open 3-D box, absorbing sphere) | pressure balance; perpendicular MHD shock with the stand-off of Spreiter et al. (1966) | magnetopause {{ bow_shock:r.gate.magnetopause | .2f }} R_mp; shock {{ bow_shock:r.gate.shock | .2f }} vs {{ bow_shock:r.gate.shock_spreiter | .2f }} R_mp |
| [Blast wave around a magnetised body](Blast-Wave) | debris shell expanding into a magnetised plasma around an absorbing sphere with a dipole (open 3-D box) | fast magnetosonic speed; debris extent scaling of Le et al. (2021) | front at {{ blast_wave:r.gate.front_mach_fast | .2f }} fast speeds; debris extent {{ blast_wave:r.gate.debris_extent_ratio | .2f }} × the scaling |
| [Plasma wakefield](Plasma-Wakefield) | wake of a laser pulse and of an electron bunch in a cold plasma (laser antenna, moving window); cavity of a strong pulse | linear theory of the wake (amplitude, wavenumber, radial structure); scalings of Lu et al. (2007) | amplitude {{ wakefield:r.laser_linear.amplitude_ratio | .3f }} (pulse) and {{ wakefield:r.beam_linear.amplitude_ratio | .3f }} (bunch) of the theory |
| [Alfvénic wave carried by a flow into a shock](Alfven-Inflow) | cold two-fluid R wave carried into an open box through a turbulent inflow face; the same wave meets a perpendicular piston shock (M_A ≈ 4.4) | cold two-fluid dispersion relation; perpendicular MHD jump conditions | k-mode within {{ alfven_inflow:r.gate.Bz_mean | .1% }} of the incident wave (mean); shock speed {{ alfven_inflow:r.piston.shock_speed_ratio | .2f }} × the jump conditions |

## Gallery

| | | |
|---|---|---|
| <img src="images/two_stream/ei.gif" width="280"><br>two-stream: phase-space vortices | <img src="images/filamentation/pair.gif" width="280"><br>filamentation: current filaments | <img src="images/reconnection/pair_sigma10_large_run2.gif" width="280"><br>reconnection: plasmoids |
| <img src="images/landau_damping/wavebreaking.gif" width="280"><br>wave breaking: sheets cross | <img src="images/relativistic_shock/pair_2d.gif" width="280"><br>relativistic shock | <img src="images/particle_orbits/orbits.gif" width="280"><br>test-particle orbits |
| <img src="images/bow_shock/gate.gif" width="280"><br>bow shock: ion density | <img src="images/blast_wave/gate.gif" width="280"><br>blast wave: field, ambient and debris ions | <img src="images/wakefield/laser_bubble.gif" width="280"><br>wakefield: cavity of a strong laser pulse |
| <img src="images/alfven_inflow/piston.gif" width="280"><br>Alfvénic wave: a shock runs into it | | |

## How the pages are made

- The gate-sized configurations are run fresh on 8 × RTX PRO 6000 Blackwell GPUs from a clean commit;
  the long production runs (shock, reconnection, expansion) are the recorded runs of `bench/results/`.
  Each page lists the commit, run id and GPU of every run it uses.
- `python scripts/build_wiki.py run` runs the missing configurations, `collect` analyses every run
  (charts, GIFs, numbers into `docs/wiki/data.json`), and `render` fills the page templates of
  `docs/wiki/`.
- What each step needs: `render` only needs the checked-in `docs/wiki/data.json` and templates. `run`
  makes only the fresh gate-size runs. `collect` also reads the recorded production runs, whose output
  directories (listed in each page's provenance table) are not in git: copy them into `output/` from the
  machine that made them. `collect` checks every recorded directory against its evidence record (run id,
  or commit and wall time) and stops if one is missing or is another run.
- Normalised units throughout: c = ε₀ = 1, lengths in c/ω_pe, times in 1/ω_pe, momenta u = γv/c,
  temperatures θ = T/(mc²).
