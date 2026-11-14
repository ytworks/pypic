# Set-up of the runs: the upstream state and how the shock meets it

Written 2026-11-14, before any recorded run. Numbers come from the planner (`plan.py`) through
`design/setup.py` → `design/setup-tables.md` ("Table n" below), from `design/estimates.md`, and from the
scratch measurements named in Sect. 7. Lengths in d_i, times in Ω_ci⁻¹, speeds in v_A of B0 unless stated.

**Purpose (owner, 2026-11-14).** The topic is not a replication of prior work; it does the next problem that
the prior kinetic work left. Read here as: the problem that `fulat2023kinetic` and `fulat2026pic` name as
future work, a fully kinetic shock in Alfvénic fluctuations (this reading was stated to the owner and is not
yet confirmed). `guo2010effect` is the source of the mechanism under test (field-line meandering, multiple
connection). Its set-up is not a target: Sect. 2 records the differences as facts.

## 1. What the prior kinetic work left open, and what this set-up has

Statements of the authors (verbatim, from the arXiv texts 2310.14388 and 2605.21986):

- `fulat2026pic`: "The level of pre-existing compressive turbulence that can be achieved on kinetic scales is
  limited by particle heating. To maintain sufficiently strong sonic Mach numbers, Ms ≳ 30, the maximum level
  of the density fluctuations should be on the order of δn/n ∼ 15%, with a few times lower amplitudes of the
  magnetic field fluctuations (Fulat et al. 2023)." and "Future studies will explore the effect of Alfvénic
  turbulence, in which the magnetic-field fluctuations dominate and may reach higher amplitudes without
  excessive heating."
- `fulat2023kinetic`: "Further studies should explore the effect of the Alfvénic turbulence, where the magnetic
  field fluctuations dominate."; "This may be different for turbulence on larger scales, as hybrid simulations
  suggest … Analogous studies in kinetic regime require using a wider simulation box and thus more
  computational resources."; "At lower Alfvénic Mach numbers the Weibel modes should be weaker, and
  correspondingly the impact of pre-existing turbulence larger."

| Statement | `fulat2023kinetic` / `fulat2026pic` | This set-up |
|---|---|---|
| magnetic fluctuations dominate | δn/n 3.5–15 %, δB/B0 ≈ 10 % at most; 2023 imposes velocity modes only | loaded δn/n = 0, rms δB_x/B0 = 0.6 (rms of \|δB\| 0.85), δv = 0.80–0.84 v_A (Table 4) |
| higher amplitude without excessive heating | amplitude capped by heating | not known: measured in the pilot (Sect. 8, heating relative to the laminar box) |
| larger scales, wider box | L_y ≈ 6 λ_si (2023), 12 λ_si (2026) | L_y = 50, 100, 200, 400 d_i |
| lower Alfvénic Mach number | M_A ≈ 32, M_s ≈ 36 (2023) | M_A = 4.5, M_s = 3.5, fast Mach 2.76 (Table 2) |
| geometry | 2D3V; θ_Bn = 90°, B0 in the plane (2023); 60°, B0 out of the plane (2026) | 2D3V; mean θ_Bn = 90°, B0 in the plane |
| mass ratio, time | m_i/m_e = 100; 10 and about 20 Ω_i⁻¹ | m_i/m_e = 100 (400 in one pair); 40 Ω_ci⁻¹ |
| nature of the fluctuations | decaying compressive turbulence, pre-evolved slabs appended upstream | waves along B0 with k_x = 0: one wave, or a few modes (Sect. 3) |

Facts that follow from the table:
- The regime differs: Fulat et al. keep M_s ≳ 30 (supernova remnants); this set-up has M_s = 3.5 (heliospheric
  shocks). The question fixed by the owner names M_A ≈ 4.5 and β ≈ 1.
- The upstream state here is not turbulence: no cascade, no oblique modes, one wavelength across the box. It
  holds what "Alfvénic" names in the quoted sentences (δB dominates, δn ≈ 0, δv ∥ δB) and nothing more. This
  is a limit of every claim of the topic (README, objections).

## 2. The source of the mechanism: the prescription of `guo2010effect`

Verbatim (arXiv:1003.5946, TeX source, Sect. 2):

> a random magnetic field is superposed on a mean field at the beginning of the simulation and is also injected
> continuously at the x=0 boundary during the simulation. The simplified one-dimensional fluctuations have the
> form B(z, t) = δB(z, t) + B_1, where B_1 is the averaged upstream magnetic field. The fluctuating component
> contains an equal mixture of right- and left-hand circularly polarized, forward and backward
> parallel-propagating plane Alfven waves. The amplitude of the fluctuations is determined from a Kolmogorov
> power spectrum P(k) ∝ 1/(1 + (k L_c)^{5/3}) … we take L_c = L_z … Note that in addition to magnetic
> fluctuations, there are also velocity perturbations with δv = v_{A1} δB/B_1 (Alfven waves). For most of the
> parts in the paper, we consider a turbulence variance σ = δB²/B_1² = δv²/v_{A1}² = 0.3

> the plasma is injected continuously from one end (x=0, in our case) of the simulation box, and reflected
> elastically at the other end (x=L_x). This boundary is also assumed to be a perfectly conducting barrier. …
> a shock propagating in the −x direction.

> in the model with at least one ignorable coordinate, the center of gyration of particles is confined to within
> one gyroradius of the original magnetic field line. The test-electrons can still move normal to the mean field
> in our model because of the field-line random walk.

| | `guo2010effect` | This set-up |
|---|---|---|
| fields | 2-D hybrid: kinetic ions, electron fluid | 2D3V, both species kinetic |
| electrons | test particles, a 100 eV shell released at 70 Ω_ci⁻¹, followed to 120 | the thermal pool, from t = 0; test species as in the README (Idea) |
| shock | piston, conducting wall; inflow Mach 4.0, about 5.6 in the shock frame | piston, conducting wall; V0 = 2.88, M_A = 4.5 in the shock frame |
| β_e, β_i | 0.5, 0.5 | 1, 1 |
| fluctuations depend on | the coordinate along the mean field, and t | the same (y, t; k_x = 0) |
| wave mix | equal R and L, forward and backward | primary state: one forward R wave (Sect. 3) |
| \|B\| upstream | not uniform (four kinds of waves) | uniform when loaded (one wave) |
| variance δB²/B² | 0.3, all components | 0.72 (0.36 in each of x and z) |
| spectrum | Kolmogorov shape, L_c = L_z = 1024 d_i | one wavelength L_y = 50–400 d_i |
| δv | v_A δB/B_1 | from the two-fluid wave: 0.94–0.99 of v_A δB/B0 for the ions (Table 4) |
| entry | loaded at t = 0, then injected at the inflow face | the same: loaded under an envelope, then carried by the face (Sect. 4) |
| Ξ = ξ_pp / L_sh | 38 (`design/estimates.md`, Sect. 2) | 4.4, 8.8, 17.7, 35.3 (Table 1) |
| box across the front | 1024 d_i | 50–400 d_i |

What is common is what the mechanism needs: a magnetic field that depends on the coordinate along the mean
field only, frozen into a flow that a piston stops, so that a field line crosses the front more than once
(Sect. 5). What differs is listed; none of it is matched on purpose.

## 3. The upstream state

The planner offers three states (`state`, `zeros` in `plan.py`). Each mode is a circularly polarized wave of
the cold two-fluid plasma that travels along +B0 (`examples/alfven_inflow/theory.py`); a circular polarization
does not rotate, a linear one does (`design/estimates.md`, Sect. 1).

1. **One wave** (`state: constant_b`, `zeros: 0`): δB_x + i s δB_z = √2 · db · B0 · e^{i(k y − ω t)}, s = +1
   for helicity R. An exact solution of the cold two-fluid equations at any amplitude, with uniform |B| and
   density.
2. **Several modes with uniform |B| at t = 0** (`zeros` ≥ 1): the same expression with a finite Blaschke
   product f(e^{i k_1 y}) in place of the exponential; |f| = 1, positive wavenumbers only. The modes disperse,
   so |B| loses its uniformity in time.
3. **Random phases** (`state: random`): weights from the slab spectrum P ∝ 1 / (1 + (k L_c)^{5/3}), random
   phases. |B| is not uniform from the start.

Facts (Table 1; std / mean of |B|² over y from the linear phases):

| State | L_y | spread at t = 0 | at t = 40 |
|---|---|---|---|
| one wave | any | 0 | 0 |
| several modes (1 zero, radius 0.3, k d_i ≤ 1) | 50 | 0.005 | 0.142 |
| several modes (1 zero, radius 0.3) | 100 | 0.005 | 0.036 |
| several modes (1 zero, radius 0.3) | 200 | 0.005 | 0.010 |
| several modes (2 zeros, radius 0.6) | 400 | 0.005 | 0.012 |
| random phases, 3 realizations | 50 | 0.21–0.36 | 0.25–0.38 |
| random phases, 3 realizations | 200 | 0.33–0.37 | 0.32–0.37 |

**Choice.** The primary state is one forward R wave at rms δB_x/B0 = 0.6. The facts behind it:
- the question says "incompressible Alfvénic fluctuations"; a |B|² spread of 0.2–0.4 (random phases) is a
  magnetic pressure fluctuation of that size at β = 2;
- the several-mode state keeps the spread below 0.05 over 40 Ω_ci⁻¹ only for L_y ≥ 200 d_i;
- at L_y = 50 d_i only one wave has a spread below 0.05 at t = 40.

The several-mode state is a test of the waveform at L_y ≥ 200 d_i. The random-phase state is a comparison in
the periodic pilot only.

**Helicity.** R is the whistler branch, L the ion-cyclotron branch. With the kinetic dispersion relation at
β_i = 1 (`design/estimates.md`, Sect. 1) an L mode has Im ω = −1.59e-3 Ω_ci at k d_i = 0.3 and −6.59e-2 at
0.5; an R mode has |Im ω| ≤ 4e-4 up to k d_i = 0.63. The planner notes every L mode above k d_i = 0.3
(`L_DAMPED_KDI`). The one wave of the runs has k d_i ≤ 0.126, where both branches are undamped in linear
theory. The pilot runs both helicities.

**Limits of the one-wave state.**
- It is not turbulence (Sect. 1).
- The distribution of the local θ_Bn is that of a circle: acute angles from 49.7° to 90° (Table 3).
- A "realization" is a shift of the wave along y, which is no new state in a periodic box. The spread between
  runs comes from the particle seed only (README, N4).
- The stability of the finite-amplitude wave in a warm plasma at β = 2 (parametric and modulational
  instabilities) is not established here: no source on it has been read. The pilot measures it (Sect. 8, C2
  and C3).

## 4. How the state is loaded and how it enters

- **At t = 0 (PIC-D24).** δB comes from a vector potential (div B = 0 bitwise). The motional E of the wave
  comes from a displacement of one species. Displacements are refused next to a non-periodic face, so the wave
  is loaded under an envelope h(x): 0 below 0.4 d_i from the wall, sin² rise to 1 at 1.2 d_i; sin² fall to 0
  over 0.8 d_i before the last d_e in front of the inflow face. Both species carry the wave's bulk velocity
  δu(y) h(x).
- **During the run (PIC-D26).** The face x+ injects the flow (flux face) with δu(y, t) on the injected
  particles and drives the wave's δB and E on a total-field / scattered-field surface. The modes have k_x = 0.
- **Known error of the face** (`examples/alfven_inflow`): a field error of 5–15 % next to the face; its cause is
  not established. Its size at this set-up's parameters is measured in P3 (Sect. 8).
- **The wall.** The wave amplitude is zero within 0.4 d_i of the wall. The plasma that the wall stops in the
  first Ω_ci⁻¹ carries the ramp of the envelope; this start-up layer stays at the wall, downstream.
- **Third direction.** nz = 2 is the engine's minimum; the runs are 2D3V with 2 · nz · ppc = 20 particles per
  column (x, y) at ppc = 5 per species.

## 5. How the wave meets the front

With k_x = 0 the upstream state depends on (y, t) only: there is no "age" along x. Every part of the front
meets the same wave; only the phase φ = k y − ω t differs along y.

**Field lines.** dx/dy = B_x / B_y gives x(y) = x_0 + (√2 db / k) sin φ: one arc per wavelength, with the
peak-to-peak excursion ξ_pp = √2 · db · L_y / π = 0.27 L_y at db = 0.6 (Table 4: 13.5, 27, 54, 108 d_i). The
field also has δB_z, out of the plane; it changes the local angles (Table 3) and not the excursion in x.

**Multiple connection.** A plane front that moves at U1 = 4.5 relative to the upstream plasma touches an arc
first at its leading point and is then connected to the same line at two points, until the trailing point has
passed: for ξ_pp / U1 = 3.0, 6.0, 12, 24 Ω_ci⁻¹ at L_y = 50, 100, 200, 400 d_i (Table 4, "connected"). The
thermal transit between the two points, (L_y / 2) / v_te, is 3.5, 7.0, 14, 28 Ω_ci⁻¹: the ratio is 0.85 at
m_i/m_e = 100 for every L_y, and it scales as √(m_i/m_e) (1.7 at 400). The README (Idea) gives 0.66, from the
random-phase estimate of ξ_pp.

**Local field.** B / B0 = (b cos φ, 1, s b sin φ), b = 0.849, |B| = 1.311 B0 at every phase (Table 3). The
angle to the normal goes from 49.7° (φ = 0) through 90° (φ = 90°) to 130.3° (φ = 180°).

**Local inflow.** The wave's bulk velocity is δv (cos φ, 0, s sin φ) with δv = −0.80 v_A for the ions at
L_y = 50 d_i (antiparallel to δB: the wave travels along +B0). Towards a front that moves at the shock speed
of the uniform state, the normal inflow speed is U1 − δv_x = 4.5 + 0.80 cos φ: 5.30 at φ = 0, 3.70 at
φ = 180°. Over the Alfvén speed of |B| that is 4.04 and 2.82 (Table 3).

**Consequence for the null (Sect. 6):** the four phases φ, −φ, 180° − φ and 180° + φ have the same acute
θ_Bn. φ and 180° ± φ differ in the inflow speed (by 1.6 |cos φ| v_A); φ and −φ differ by the sign of B_z,
which is no symmetry of an electron–ion plasma. The local state is a point on the full circle of φ, not a
value of the acute θ_Bn.

**Motion of the pattern.** The wave moves along the front at v_ph = ω / k = 1.058 (L_y = 50) to 1.003
(L_y = 400). The phase at a fixed point of the front goes once around the circle in 2π / ω = 47, 97, 198,
399 Ω_ci⁻¹ (Table 4). So within 40 Ω_ci⁻¹ a point of the front sees 0.85 of the circle at L_y = 50 d_i and
0.10 at L_y = 400 d_i. In the frame that moves along y at v_ph the upstream state is steady. In
`guo2010effect` forward and backward waves are mixed: no such frame exists there.

**The shock with the wave.** The planner's shock (V0 = 2.881, r = 2.780, V_sh = 1.619, M_A = 4.50; Table 2)
is that of the uniform state with B0 alone. The wave raises |B| to 1.311 B0 (M_A = 3.43 on |B|) and adds its
pressure and momentum. The shock speed, the compression and the thickness L_sh with the wave are not
predicted here; P2 and P3 measure them (Sect. 8). L_sh = 0.68 M_A d_i = 3.06 d_i in Ξ is an estimate
(`design/estimates.md`) until P2.

**The face and the length of the box** (Table 2, speeds of the uniform state):

| L_x | plasma that entered through the face reaches the shock | a fast-mode signal from the face does |
|---|---|---|
| 100 | 22.2 | 16.3 |
| 200 | 44.4 | 32.6 |

With L_x = 200 d_i the shock meets only preloaded plasma for 40 Ω_ci⁻¹, but a fast-mode signal from the face
(the face error above, the fall of the envelope) can reach it from 32.6 Ω_ci⁻¹. Either L_x ≥ 246 d_i
(200 × 40 / 32.6), or the evaluation ends at 32.6 Ω_ci⁻¹. P3 decides, from the measured size of the signal.

## 6. The laminar null

The null of the README (Idea a) asks what a laminar shock gives for the local states that the front of the
turbulent run meets. From Sect. 5:

- A local state is a phase φ on the full circle: B = B0 (b cos φ, 1, s b sin φ) and the bulk velocity
  −V0 x̂ + δv (cos φ, 0, s sin φ), uniform. The first plan (README: θ_Bn = 90, 80, 70, 60, 45°, B0 tilted in
  the plane, no bulk velocity) does not cover it: it has no B_z, no δv, and one sign of B_x.
- The null runs sample the circle of φ. Their number is fixed with the engine question below.
- The null is a stationary proxy: in the run the local state at a point of the front changes in time (Sect. 5,
  motion of the pattern) and along y within one arc. This is a limit of the null, to be stated with it.

**Engine.** A uniform B with a component normal to the reflecting wall cannot be given as `B0`: the wall
zeroes the normal B of the Yee fields (`docs/design/engine.md`, PIC-D12, "Unchanged"). Static external fields
(`ExternalFields`, PIC-D12, explicit solver) may be normal to any face. Whether a piston shock with an
external normal B runs correctly (the wall conditions on the difference fields, the flux face, the Gauss and
energy contracts) is not known. A smoke run on the present engine decides whether an engine change (E4, a
plan on main) is needed.

## 7. Resolution

**Grid heating, measured** (scratch runs on 458bfd9, not recorded runs: laminar periodic boxes, 4 × 10 d_i,
1.0 Ω_ci⁻¹, 5 cells per d_e, ppc 5, 4 filter passes; change of the kinetic energy of a species):

| ω_pe/Ω_ce | λ_De / dx | ΔT_e/T_e per Ω_ci⁻¹ | first half | × 40, linear | ΔT_i/T_i | budget drift |
|---|---|---|---|---|---|---|
| 10 | 0.35 | +0.0485 | +0.0214 | +1.94 | −0.0107 | +7.5e-3 |
| 4 | 0.88 | −0.0041 | −0.0034 | −0.16 | −0.0018 | +3.4e-4 |

- At ω_pe/Ω_ce = 10 the electrons would gain about twice their temperature in 40 Ω_ci⁻¹ if the rate stayed.
  The rate grows within the run (first half +0.021, second half +0.027).
- At ω_pe/Ω_ce = 4 the electrons lose 0.4 % in 1 Ω_ci⁻¹. A linear extrapolation gives −16 % in 40 Ω_ci⁻¹; one
  Ω_ci⁻¹ does not show whether the rate stays. The full-length pilot decides (Sect. 8, C1).

**What ω_pe/Ω_ce = 4 changes** (candidate):
- T_e = 0.03125 m_e c² (0.005 at ω_pe/Ω_ce = 10): 10 T_e = 0.31 m_e c², so the tail above 10 T_e is mildly
  relativistic. v_A = 0.025 c (0.00995 c).
- v_te / v_A, the access ratio and Ξ depend on β_e and the mass ratio only: unchanged.
- Steps for 40 Ω_ci⁻¹: 160000 instead of 400000 (Table 2). Cost at 1.1e10 particle-steps per second: about 2 h
  for 200 × 50 d_i (5 h), 16 h for 200 × 400 d_i (40 h).
- `fulat2023kinetic` uses 20 cells per λ_se and 20 particles per cell per species; this set-up 5 cells per d_e
  and 10 per column and species (Sect. 4).

## 8. The pilot

Order: P1a (below, fixed here), then P1b, P2, P3 (their configurations are written after P1a).

**P1a: periodic boxes without a shock** (`configs/p1a-*.yaml`; Table 5). 10 × 50 d_i, grid 500 × 2500 × 2,
ppc 5, one GPU each, 40 Ω_ci⁻¹ unless stated. The plasma drifts at −V0 as in the wall runs.

| Config | ω_pe/Ω_ce | State | rms δB_x/B0 | Ξ | Asks |
|---|---|---|---|---|---|
| `p1a-lam-w4` | 4 | laminar | 0 | 0 | grid heating over the full length (C1); noise floors |
| `p1a-lam-w4-ppc10` | 4 | laminar, ppc 10 | 0 | 0 | dependence of the heating on ppc |
| `p1a-lam-w10` | 10 | laminar, 10 Ω_ci⁻¹ | 0 | 0 | the rate at λ_De/dx = 0.35 in this box |
| `p1a-wave-r-w4` | 4 | one wave, R | 0.6 | 4.41 | stability and heating of the primary state (C2, C3) |
| `p1a-wave-l-w4` | 4 | one wave, L | 0.6 | 4.41 | the same for the other helicity |
| `p1a-weak-r-w4` | 4 | one wave, R | 0.1 | 0.74 | amplitude dependence; the state of the Ξ < 1 set |
| `p1a-multi-w4` | 4 | 5 modes (1 zero, radius 0.3, k d_i ≤ 1) | 0.6 | 2.69 | measured \|B\|² spread against the linear table |
| `p1a-random-w4` | 4 | random phases, 3 modes | 0.6 | 3.97 | comparison: heating and δn/n with a \|B\|² spread of 0.36 |

**Measured** (`analyze.py pilot`; plane z = 0 at 80 frames, averages over blocks of 1 d_i = 50 × 50 cells,
12500 particles per species and block): thermal energy per particle of each species in the local bulk frame;
δn/n; std / mean of |B|²; rms δB_x/B0 of the k_x = 0 part and its share in the loaded modes; the share of
δB_x² not in k_x = 0; ξ_pp and Ξ from the field lines; the budget; the share of electrons above 10 and 30 T_e0.

**Criteria, fixed before the runs:**
- **C1 (resolution).** `p1a-lam-w4`: |ΔT_e/T_e| ≤ 0.05 over 40 Ω_ci⁻¹. Reason for the number: a shift of T_e by
  5 % changes the Maxwellian share above 10 T_e by a factor 1.58 (1.70e-4 → 2.67e-4), against the factor 7.7
  that the model predicts for the excess (`design/model.md`). If C1 fails, ω_pe/Ω_ce = 4 at 5 cells per d_e is
  not used; the options of the README (more cells, the implicit scheme) are then tested the same way.
- **C2 (the wave stays incompressible).** `p1a-wave-r-w4`: δn_e/n ≤ 0.05 and std / mean of |B|² ≤ 0.05 at
  every frame, on 1 d_i blocks. The noise floor of both comes from `p1a-lam-w4` and is reported with them.
- **C3 (the wave stays).** `p1a-wave-r-w4`: rms δB_x/B0 of the k_x = 0 part ≥ 0.54 at t = 40 (0.9 of the
  loaded value).
- **Several modes.** `p1a-multi-w4`: measured |B|² spread within a factor 2 of the linear one (0.039, 0.076,
  0.142 at t = 10, 20, 40). If so, Table 1 is used for the several-mode state at L_y ≥ 200 d_i; if not, that
  state is measured directly before use.
- **Reported without a criterion:** the heating of T_e and T_i of every wave run relative to `p1a-lam-w4` (it
  fixes the tolerance of N6 and answers the "without excessive heating" of Sect. 1); the R–L differences; the
  ppc dependence; the share above 10 T_e0 relative to the laminar box (the no-shock part of N2).

If C2 or C3 fails for R and holds for L, the primary helicity changes to L, and this file records it. If both
fail, the primary state is not usable at this amplitude; the amplitude dependence (`p1a-weak-r-w4`) then shows
where it is.

**What P1a cannot show.** The box is 10 d_i wide: modes with 0 < k_x d_i < 0.63 do not fit. An instability
that needs them is not seen. P1b covers it.

**P1b.** The primary state in 50 × 50 d_i on 8 GPUs (k_x d_i ≥ 0.126), 40 Ω_ci⁻¹: the share of δB_x² not in
k_x = 0 and C2, C3 again.

**P2.** Laminar wall run (δB = 0), 200 × 50 d_i: shock speed, compression, L_sh (the foot and ramp widths that
define Ξ), ripples, the electron spectrum of the 90° shock; grid heating upstream with the wall.

**P3.** Wall run with the primary state, 200 × 50 d_i: the upstream state at arrival against Sect. 5 (N6); the
field error at the face and the time its signal reaches the shock (L_x, Sect. 5); shock speed and compression
with the wave; the output cadence that the field-line labelling needs.

After P2, the amplitude of the Ξ < 1 set is fixed: with one wave, rms δB_x/B0 = 0.15 at L_y = 50 d_i gives
ξ_pp = 3.38 d_i and Ξ = 1.10 with the estimated L_sh, not Ξ ≈ 0.85 (that number came from the random-phase
estimate); 0.1 gives Ξ = 0.74.

## 9. Open points

- Stability of the finite-amplitude wave at β = 2: no source read; P1a and P1b measure it.
- Grid heating over 40 Ω_ci⁻¹ at ω_pe/Ω_ce = 4 (C1).
- The face error (5–15 % in `examples/alfven_inflow`), cause not established; L_x against 32.6 Ω_ci⁻¹ (P3).
- Shock speed, compression and L_sh with the wave (P2, P3).
- The null on the circle of φ: engine route (external B or E4) and the number of runs (Sect. 6).
- The start-up layer at the wall (Sect. 4) and the fall of the envelope at the face: seen in P3.
- "A slab of waves is not turbulence" and the periodic return of every field line: limits that stay.
