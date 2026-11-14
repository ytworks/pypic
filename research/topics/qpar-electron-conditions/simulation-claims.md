# What a simulation can state — qpar-electron-conditions

The research of this topic is done by simulation. This file records what published simulations of collisionless
shocks have stated, on what basis, with which limits, and how they were tied to observations; then what the runs
planned here can and cannot state; then the boundary areas (neighbouring methods, laboratory shocks, border
regimes, domain boundaries). It complements `inputs.md` (the observed side) and `literature.md` (the map).

Rules of this file: a statement of a work is taken from its note in `research/literature/` (marks (S) and (A) as in
`literature.md`); a quotation was checked in the text read; "own" marks our estimate or reading, which is not a
result. Nothing here is a decision: criteria, definitions and the choice among set-ups are the owner's.

Read ranges: bodies of `gupta2024electron` (S), `gupta2025speed` (S), `kato2015particle` (S) (re-read for the
passages in §3 and §7), and of the nine works named "body" in the notebook entry of 2026-11-14 (second part);
abstracts only for twenty further works (each note says so in its `Read:` line). The search was 26 keyword queries
on OpenAlex plus the citations already collected; it is not exhaustive.

## 1. What published simulations state about electrons at quasi-parallel shocks

| Work | Method, dimension, frame | m_i/m_e | v_sh/c | M_A | θ_Bn | Run time | Statement on electrons | Tie to observations |
|---|---|---|---|---|---|---|---|---|
| `gupta2024electron` (S) | PIC, quasi-1-D (5 cells), upstream frame, moving piston | 16-1836 | 0.067-0.267 | 5-40 | 30° | 300 ω_ci⁻¹ (1836 shorter) | η_e ≲ 0.5 %, ε_e ~2 % (~1 % at 1836); slope ~4 only for M_A, M_s ≳ 10; stall near p ~ m_i v_pt at M_A ≲ 10 | consistency with SN 1006 (speed, M_A ~50, X-rays in quasi-parallel regions); no event |
| `gupta2025speed` (S) | same runs; test electrons in field snapshots; Monte Carlo | 100 (1836 in an appendix) | 0.133 | 5, 20 | 30° | ~300 ω_ci⁻¹ | injection threshold u > v_sh; D_1,e ~ D_Bohm(B_0), D_2,e ≈ 0.1 D_1,e at M_A 20; a formula for p_max(t) | none |
| `park2015simultaneous` (S) | PIC, 1-D, wall (downstream frame) | 100 (one 400) | u 0.1 (one 0.05) | 20 | 30° | 310 Ω_cp⁻¹ | both species p⁻⁴; K_ep 3.8e-3; 5.5e-3 (400); 1.2e-3 (0.05 c) | K_ep of Tycho |
| `kato2015particle` (S) | PIC, 1D3V, two walls (downstream frame) | 30 | 0.46 | ~28 | 30° | 590 ω_cp⁻¹ | both species accelerated; measured gain rates | — (M_A "still below remnant values") |
| `shalaby2022mechanism` (S) | PIC, 1D3V, wall | 1836, 100 | 0.133 | 5.3, 21.3 | 0° | — | efficient pre-acceleration only if the Mach number at the front is below √(m_i/m_e)/2; needs the real mass ratio | SN 1006 (post-shock M_A) |
| `zekovic2025slams` (S) | PIC, 1-D and narrow 2-D | 32 | 0.033-0.33 | 20-300 | 0° | 217 ω_ci⁻¹ (1-D, M_A 80), 61 (2-D) | quasi-periodic gain at SLAMS with B_⊥/B_0 ≳ 10 for M_A ≳ 80 | — |
| `fang2019early` (S), `fang2022kinetic`, `yu2024nonthermal` | PIC, 1-D | 30 | 0.27 | 16.5 (`fang2019early`) | 15-45° (`fang2019early`) | — | early pre-acceleration | — |
| `bessho2023electron`, `bessho2025electron` (abstracts) | PIC, 2-D, shock transition | — | — | Earth-like | quasi-parallel | ~20 ω_ci⁻¹ | reconnection and thinning current sheets heat and pre-accelerate; reach "the injection energy for diffusive shock acceleration" | parameters of Earth's bow shock; MMS current sheets |
| `tsubouchi2004full` (abstract) | full particle, 1-D | — | — | — | quasi-parallel | — | SLAMS with a whistler precursor that "behaves as a new shock front" | scales agree with hybrid runs and AMPTE/Cluster |
| `zhou2025planar` | semi-implicit PIC, 2-D (1-D for mass ratio), shock rest frame | 25 (100, 400 in 1-D) | — | — | quasi-parallel | 20 ω_ci⁻¹ | none on non-thermal electrons; whistler precursors at 100 and 400, not at 25 | — |
| `zhang2024kinetic` (abstract) | PIC, 2-D | — | — | — | quasi-parallel | — | formation of laboratory quasi-parallel shocks | planned NIF experiments |

"—" means: not in the text we read. No row measures D̃_e or an injection rate at the parameters of one of our four
cases, and no row confronts an event with numbers. That is the statement of the survey; the search is not
exhaustive.

## 2. Limits that recur, with their sources

1. **Dimension.**
   - `gupta2024electron` (S): "1D simulations tend to exaggerate this effect because they do not allow for a change
     in B_x and for shock rippling"; "Future investigation is necessary to validate whether the trend of precursor
     shortening with increasing shock speeds persists in multidimensions."
   - `zhou2025planar`: "at least two spatial dimensions are critical for … capturing the complex dynamics of
     quasi-parallel shocks, including surface rippling, shocklets, SLAMS, magnetic reconnection and jets".
   - `pohl2020pic` (S): "only 3D simulations can provide a complete picture"; the ranking of processes can change
     with dimension and mass ratio. `vanmarle2022diffusive`: "injection is, fundamentally, a 3D process".
   - `marcowith2020multi`: the Bell instability "can be triggered in 1D" at non-relativistic quasi-parallel
     shocks; the step from 2-D to 3-D matters for the scattering. `xu2020electron` (S): 1-D admits only waves
     along x. `morris2022preacceleration` (S): at θ_Bn 63° in 2-D the reflected share is ≈ 0.03 % with the
     field in the plane and ≈ 4 % out of it.
   - `orusa2025criteria`: a worked test of "when is a lower dimension enough" (perpendicular shocks, ions), with
     the particle number varied to exclude statistics.
   - `caprioli2015simulations`: reformation is "unlikely coherent along the shock surface" at real shocks; a
     quasi-1-D box makes it coherent by construction (own reading).
2. **Mass ratio.** `shalaby2022mechanism` (S) "necessarily requires realistic mass ratios"; `shalaby2024energy`: a
   reduced mass ratio changes electron acceleration in opposite directions at high and low M_A. `zhou2025planar`:
   precursor whistlers appear at 100 and 400, not at 25. `marcowith2020multi`: the whistler criterion
   M_A/√(m_i/m_e) < 1 and oblique whistlers depend on it. `gupta2024electron` (S): injection efficiency "does not
   show much dependence on shock speed and m_R" (their Figure 7); their K_ei at 1836 is an extrapolation.
3. **Speed v_sh/c.** All fully kinetic runs above have v_sh ≥ 0.033 c; the events have 0.0013-0.008 (`inputs.md`).
   - Published dependences: `gupta2024electron` (S), v_sh 0.067-0.267 c: η_nt,e without much dependence;
     ε_nt,e/ε_B "decreases with the shock speed"; "there are fewer energetic electrons for high shock speed",
     attributed to superluminal phases that 1-D exaggerates. `park2015simultaneous` (S): K_ep 3.8e-3 at 0.1 c and
     1.2e-3 at 0.05 c (one pair). `vanmarle2022diffusive` (ions, oblique): ε_inj ≈ 12 ε_sh at 0.133 c and ≈ 23 ε_sh at
     0.264 c, "not strongly dependent on the shock speed at non-relativistic speeds". `marcowith2020multi`: K_ep
     ~1e-3 "roughly proportional to u_sh/c".
   - The usual licence: `haggerty2019dhybridr`, "as long as there is a clear separation of scales … the underlying
     physics can be studied fruitfully" — stated for ions, not tested.
   - `gupta2025speed` (S): "the electron maximum momentum is smaller than protons' unless t ≫ (β_pt⁻²) ω_ci⁻¹".
   - Own reading: an electron of a given p/(m_i v_sh) moves at another u/v_sh in the run than in nature, so
     statements in p̃ do not carry over without a measured dependence (this is idea 2 and null test N4).
4. **Run time.** Full PIC covers "a couple of ion gyro-periods" to a few 100 ω_ci⁻¹; hybrid 10³
   (`marcowith2020multi`). `gupta2024electron` (S) limits its stall to the times explored; `pohl2020pic` (S): a few
   Ω_ci⁻¹. The PSP shock is ~10⁵ ω_ci⁻¹ old (`inputs.md`).
5. **Box, boundaries, launching, frame.**
   - `pohl2020pic` (S): a reflecting wall may cause unphysical initial reflection. `kumar2021nonthermal` (S): the
     result depends on the launching method. `vanmarle2022diffusive`: "an insufficient length of the simulation box
     in the direction of flow may change the outcome". `swisdak2023comparison` (abstract): the maximum energy is
     "limited by the size of the computational domain". `crumley2019kinetic` (S): the box must be wider than the
     gyroradius of the current-carrying ions.
   - Frame: `gupta2024electron` (S) runs in the upstream rest frame: "when the Debye length … is marginally resolved
     or under-resolved, the upstream rest frame is an excellent choice to avoid numerical heating for the parameter
     space used in this paper". Against `park2015simultaneous` (S) (same parameters, downstream frame): "our η_nt,e
     is slightly larger by a factor of ≲ 2 … it is likely that the choice of the upstream frame is less prone to
     the noise produced by drifting plasmas and hence more accurate for capturing the growth of the self-generated
     waves that scatter reflected particles".
   - `zhou2025planar` relaxes MHD jump states in the shock rest frame; `simon2026faux` replaces the shock by a
     prescribed injection and return probability and shrinks the upstream length in time.
6. **Noise.**
   - `kato2015particle` (S): "high-energy particles in the simulations generally suffer from the energy loss
     process due to the stopping power of the plasma"; "The energy loss rate is inversely proportional to the
     number of super-particles in the electron skin depth volume, N_e"; in 1-D, dE_kin/dt ~ −(1/2N_e) m_e c² ω_pe
     for relativistic particles. With N_e = 50 and 133 "the acceleration process becomes significantly inefficient
     and is almost completely nonfunctional"; "Even in the case of N_e = 667 (N_PPC = 40), the acceleration
     efficiency is still affected by the energy loss". Their run uses N_e = 2667.
   - `vanmarle2022diffusive`: "if a simulation is too long, artificial heating can substantially change the upstream
     conditions". `marcowith2020multi`: PIC noise is "less adapted" to tails; Vlasov codes are "virtually
     noise-free". `crumley2019kinetic` (S) uses 4 particles per cell; `gupta2024electron` (S) 200 per cell and
     species with 32 filter passes, "results … insensitive to them"; `gupta2025speed` (S): "We have validated our
     findings against different values of m_R, space resolution, and number of particles per cell".
7. **Extrapolation by formula.** K_ei at the real mass ratio (`gupta2024electron` (S)); K_ep ~1e-3 at the real mass
   ratio and p_max(t) at late times (`gupta2025speed` (S), whose Monte Carlo starts from an imposed state);
   PIC results as "robust scalings … as sub-grid prescriptions" (`marcowith2020multi`); injection in MHD-PIC "ad
   hoc", η = 2e-3 "encapsulates our ignorance" (`bai2015magnetohydrodynamic`).
8. **Ties to observations are mostly indirect** (§3).

## 3. How simulations have been confronted with observations

| Way | Works | What was compared | What it could not give |
|---|---|---|---|
| A global ratio at remnants | `park2015simultaneous` (S), `gupta2024electron` (S) | K_ep, shock speed and M_A of Tycho, SN 1006 | no event; no time or size |
| Virtual spacecraft / synthetic probe through the run | `bohdan2021magnetic` (S) (PIC, 422 Cassini crossings), `kropotina2023weibel` (abstract; hybrid), `haggerty2019dhybridr` (against Johlander 2018), `johlander2022quasi` (abstract; MMS and Vlasiator), `gingell2017mms` (abstract; MMS and hybrid) | field peaks, ripples, SLAMS reformation along a path | no electron spectra; `bohdan2021magnetic` (S) has no quantitative measure and rescales the peaks |
| Hybrid run at event parameters | `giacalone2025source` (abstract; the PSP event), `swisdak2023comparison` (abstract; Voyager 2 termination shock, PIC against hybrid) | proton or pickup-ion spectra | electrons |
| Test electrons in prescribed or hybrid fields | `guo2015acceleration` (S) (Saturn), `burgess2006simulations` (S), `hanusch2020electron`, `liseykina2019electron` (abstract), `trotta2019electron` (S), `battarbee2021vlasov` (abstract) | electron energies reached | fields at electron scales; the origin of the turbulence is an input |
| A model in prescribed fields fitted to an event | `liu2017fermi` (abstract; THEMIS), `lindberg2025fermi` (abstract; MMS), `shi2025compound` | electron energy gain at Earth | self-consistency; `shi2025compound` shows sufficiency, not necessity |
| A mechanism checked against in-situ data | `amano2020observational` (S), `katou2019theory` (S), `lindberg2024electron` (A), `amano2022nonthermal` (S) | stochastic shock drift acceleration against MMS | quasi-parallel shocks are not covered |
| Transport from data | `jebaraj2026wave`, `kouloumvakos2026wave` | proton κ_eff and wave amplitudes at the PSP event | electrons; a kinetic origin of κ |
| Laboratory | `fiuza2020electron` (abstract), `grassi2023electron`, `schaeffer2019direct` (abstract), `shankar2026particle` (A) | electron spectra and velocity distributions at laser or LAPD shocks | the heliospheric regime; quasi-parallel shocks "more difficult to form in the laboratory and to simulate" (`zhang2024kinetic`, abstract) |
| Monte Carlo with inputs | `nyberg2026electron` (abstract), `john2024effects` (abstract) | electron energies at coronal parameters | the mean free path profile and injection are inputs |

## 4. The planned statements: what a run can and cannot state

Each block lists the statement, what a run of the planned kind gives directly, what needs a transfer, and the
readings. No ranking.

### 4.1 Spectra and p_max(t) in a run
- Directly: whether a tail forms within the run time at the given M_A, θ, β, mass ratio and v_sh/c, in the run's
  dimension; its slope and p_max(t) in units of m_i v_sh.
- Not from the run alone: "the shock does not accelerate" from "no tail in the run". `gupta2024electron` (S)
  limits its stall to the times explored; the events are older (PSP) or of similar age (Earth, Saturn) in ω_ci t.
- Readings of a missing or weak tail: the physics at those parameters; too short a run; 1-D (superluminal
  phases, no oblique whistlers); numerical drag (§5); the wall or the frame (§2.5).

### 4.2 Injection rate
- Directly: a fraction under a named definition, in a named region, at a named time.
- Facts bearing on the definition: `gupta2025speed` (S) states the threshold u > v_sh. The reviewer found that at
  the events' electron sonic Mach numbers (PSP 0.9-1.3, Saturn ~1, Earth 0.3-0.4) thermal electrons already
  satisfy it, so that "the share passing p_inj" does not measure injection there; `gupta2024electron` (S) defines
  the non-thermal fraction by a momentum threshold downstream (p ≥ 5 p_th,e, as read by the reviewer). `caprioli2015simulations` shows the form of a
  finished injection result for ions (a threshold, a number of cycles, a loss probability, a fraction of 3.6 %).
- Readings: the definition is open and is the owner's to fix; two definitions side by side are possible.
- Transfer: to 3-D ("injection is, fundamentally, a 3D process"); to the real v_sh/c (published: weak dependence
  over 0.067-0.267 c for electrons in 1-D; one pair in `park2015simultaneous` (S) with a factor ~3).

### 4.3 D̃_e(x, p̃)
- Directly: the diffusion coefficient of test electrons in the saved fields of the run, where ⟨Δx²⟩ is linear in
  time. The test electrons do not suffer numerical drag; the fields they move in do carry the noise of the run.
- Not from the run alone: D̃_e at momenta whose scattering waves the run does not hold (above the stall if the
  tail stalls; oblique whistlers in 1-D); D̃_e far upstream of the box; D̃_e at a later age.
- Facts: published values exist only at M_A 5 and 20, θ 30° (`gupta2025speed` (S)); no seed scatter is given there;
  D_1,e is 0.6 D_B0 at one time and ~1 as a time average (reviewer, read in the source). `reville2013universal`
  (abstract): transport near parallel shocks is non-uniform and energy-dependent, below Bohm in places.
  `marcowith2020multi`: Bohm is favoured by habit.
- Readings on frozen snapshots: valid if the waves move slowly against the electron (the reviewer proposed a
  test against tracked PIC electrons at the same x and p̃).

### 4.4 Dependence on v_sh/c and the transfer to the events
- Directly: the measured quantities at the speeds run (planned 0.033, 0.067, 0.133 c).
- Transfer: a power law over 1.2-1.8 decades beyond the lowest point is an assumption. The reviewer's
  propagation: three points give an extrapolated 1σ of 2.2 σ_pt at PSP and 3.9 σ_pt at Earth and Saturn.
- Facts against a simple law: the superluminal condition θ* = acos(v_sh/c) and the escape condition
  v_∥ > (B_⊥/B_0) v_sh (`zekovic2025slams` (S)) change with speed as thresholds, not as powers; `gupta2024electron` (S)
  sees fewer energetic electrons at higher speed for that reason.
- Own: the slow runs are also the ones most exposed to numerical drag (§5), so a measured trend in v_sh/c mixes
  physics and noise unless the particle number is varied at each speed.

### 4.5 Confronting the events through a transport model
- Directly: a predicted spectrum for stated inputs (D̃, injection, compression, age, size).
- What such a chain had to concede elsewhere: `vanmarle2022diffusive` chose the injection rate "in between" two
  2-D values; `bai2015magnetohydrodynamic` and `simon2026faux` take injection as an input; `galishnikova2026self`
  (preprint, abstract) finds the long-term state "largely insensitive to the details of CR injection".
- Not from the chain: necessity. A model that reproduces an event shows that its ingredients suffice
  (`shi2025compound`). An intensity needs a density and a background (`inputs.md`: Saturn's density and T_e are not
  measured; the reviewer: three of the four "no" classes are detection-limited by their own authors).
- Age: the reviewer's estimate is that a half-run comparison cannot license carrying D̃ from 200 to 10⁵ ω_ci⁻¹.

### 4.6 The candidates H1-H8
- H1a, H1b (time, size): through 4.5 only. H2, H3 (M_A bounds): orderings among runs that differ in M_A alone, at
  the run's mass ratio and dimension; the two predecessors disagree (`gupta2024electron` (S) against
  `shalaby2022mechanism` (S)) at different θ (30° against 0°) and temperature.
- H4 (seeds): the planned seeds are a passive species. `caprioli2018diffusive`: the efficiency is "not an
  intrinsic property of the shock, but rather it is regulated by the amount of seeds available", through the
  current of reflected seeds (ions; seed density "much larger than … in the solar wind"); "the same
  re-acceleration mechanisms should apply also to … seed electrons". A passive species shows how seeds gain energy
  in the shock's own turbulence, not what the seeds drive.
- H5, H7 (upstream structure, history): outside a uniform-upstream run by construction; only "the shock alone
  suffices or not". `caprioli2018diffusive`: extrinsic turbulence "would need to be assessed on a case-by-case
  basis". `kang2019electron` (abstract): weak shocks need fossil electrons.
- H6 (non-DSA gains): the energy budget of tracked electrons separates field components in the run; reconnection
  in the transition (`bessho2023electron`, `bessho2025electron`) and oblique whistlers are absent in 1-D.
- H8 (detection): needs the predicted intensity, so 4.5.

## 5. Numerical drag against the acceleration rate (own estimate, unverified)

- Basis: the 1-D loss rate of `kato2015particle` (S), dE/dt ~ −(1/2N_e) m_e c² ω_pe, N_e = electrons per electron
  skin depth. Gain at the Bohm rate with compression 4 and D_2 = D_1/10: dE/dt ≈ (v_sh/c)³ √(m_i/m_e) /
  (5.6 M_A D̃) in m_e c² ω_pe (conventions leave a factor 2-4).
- Check on a published run: Gupta's B1 (M_A 20, mass ratio 100, v_sh 0.133 c) gives 2.1e-4 from this formula; the
  published p_max ≈ 9 m_i v_pt at 300 ω_ci⁻¹ is a mean rate of 2.0e-4.
- Particle numbers: `gupta2024electron` (S) has 200 per cell and species × 5 cells × 10 cells per d_e = 10⁴ electrons
  per d_e (loss 5e-5). The plan has ~1000 per cell column over both species = 5000 electrons per d_e (loss 1e-4).
  `kato2015particle` (S) has 2667 and a measured gain of ~1e-2 (loss/gain ~2 %); their "still affected" case (N_e 667) is ~7 %.

| Run | v_sh/c, M_A, m_i/m_e | Bohm gain (D̃ = 1) | loss/gain at N_e = 5000 |
|---|---|---|---|
| W (speed not fixed in the plan; 0.133 c taken) | 0.133, 3, 1836 | 6.0e-3 | 0.02 |
| V2 | 0.133, 5, 1836 | 3.6e-3 | 0.03 |
| P0 | 0.133, 7.5, 1836 | 2.4e-3 | 0.04 |
| P2 | 0.133, 20, 1836 | 9.0e-4 | 0.11 |
| P5a | 0.067, 7.5, 1836 | 3.0e-4 | 0.33 |
| V1 | 0.133, 20, 100 | 2.1e-4 | 0.48 |
| E (plan: 0.05-0.07 c) | 0.06, 8, 1836 | 2.1e-4 | 0.48 |
| S1 | 0.133, 100, 100 | 4.2e-5 | 2.4 |
| P5b | 0.033, 7.5, 1836 | 3.7e-5 | 2.7 |
| P6 | 0.067, 20, 100 | 2.6e-5 | 3.8 |
| S2 | 0.033, 100, 100 | 6.6e-7 | 150 |

- Readings:
  - Taken at face value, the same formula gives loss/gain ≈ 0.24 for Gupta's own B1, while they report p_max in
    agreement with their formula and insensitivity to the particle number. So either the 1-D loss formula
    overstates the drag in a filtered quasi-1-D run, or its effect is inside their agreement. Not checked.
  - The formula is for relativistic particles in 1-D. A binomial filter of 32 passes at Δx = d_e/10 lowers the
    amplitude at k = 1/d_e by ~8 % only (own calculation), so it does not remove the wake scale.
  - D̃ is in units of Bohm in B_0. Where the field is amplified (M_A 100: B_⊥/B_0 ≳ 10 in `zekovic2025slams` (S)) the
    gain is faster than the table says; where D̃ > 1 (low M_A, far upstream) it is slower.
  - The gain scales as (v_sh/c)³ at fixed N_e, so the low-speed points of the v_sh/c scan and the Saturn runs are
    the most exposed. The cost of raising N_e is linear in N_e.
  - The measurement of D̃_e with test electrons in saved fields is not subject to this drag; the PIC spectrum,
    p_max(t) and the population that reaches high energy are.
- A second difference of method, a fact: the plan runs in the downstream frame against a wall; the direct
  predecessor chose the upstream frame against numerical heating and noise (§2.5), and differs from the
  downstream-frame run of `park2015simultaneous` (S) by a factor ≲ 2 in η_nt,e.

## 6. Boundary areas

### 6.1 Neighbouring methods and what each leaves to the others
| Method | Works | Holds | Takes as input |
|---|---|---|---|
| Hybrid (kinetic ions, fluid electrons) | `caprioli2014simulationsb`, `caprioli2015simulations`, `caprioli2018diffusive`, `haggerty2019dhybridr`, `orusa2025criteria`, `kylivnyk2026hybrid` (abstract), `giacalone2025source` (abstract) | ion-scale turbulence to 10³ ω_ci⁻¹, 2-D and 3-D | electron physics (an equation of state) |
| Hybrid + test electrons | `burgess2006simulations` (S), `hanusch2020electron`, `liseykina2019electron` (abstract), `trotta2019electron` (S) | electron orbits in ion-scale fields | electron-scale waves; injection energy |
| Hybrid-Vlasov with Vlasov electrons in snapshots | `battarbee2021vlasov` (abstract), `johlander2022quasi` (abstract) | noise-free distributions, global geometry | electron feedback on the fields |
| MHD-PIC / PIC-MHD | `bai2015magnetohydrodynamic`, `vanmarle2022diffusive`, `galishnikova2026self` (abstract) | long times and large boxes with kinetic energetic particles | the injection fraction; "equally applicable to sufficiently energetic CR electrons" |
| MHD + Vlasov-Fokker-Planck | `reville2013universal` (abstract) | transport with self-generated fields | injection |
| Faux shock | `simon2026faux` | late-time acceleration without the shock structure | injection, energy gain and escape probability per cycle, from full runs |
| Monte Carlo, focused transport | `nyberg2026electron` (abstract), `john2024effects` (abstract), the Monte Carlo of `gupta2025speed` (S) | time and size of real events | mean free path profile, injection |
| Data-constrained transport | `shi2025compound`, `jebaraj2026wave`, `kouloumvakos2026wave` | the event's own fields and κ | a kinetic origin |
| Semi-implicit PIC | `zhou2025planar` | larger cells (0.5-4 d_e) and 2-D quasi-parallel structure | electron tails not studied; noise ~5 % |

- Each method of this table either has no kinetic electrons or takes the electron injection (or the fields at
  electron scales) as an input. `vanmarle2022diffusive`: "the exact
  fraction of injected electrons still should be defined with further PIC numerical experiments";
  `caprioli2015simulations`: a characterization "without free parameters … can be achieved only with kinetic plasma
  simulations". This is the place of the planned measurement in the chain of methods.
- The engine has an implicit scheme (exactly energy conserving). Whether it serves the long, slow runs is not
  examined here.

### 6.2 Laboratory shocks
- `schaeffer2019direct` (abstract; OMEGA, velocity distributions in a shock precursor; ~750 km/s, M_A 15, ~4 proton
  gyrations per `orusa2025criteria`); `fiuza2020electron` (abstract; NIF, electrons to ~500 keV, more than 100 times
  thermal); `grassi2023electron` (D from tracked electrons, agreement with NIF); `vanthieghem2024electron`
  (abstract; heating model in Weibel-mediated shocks); `lezhnin2021kinetic` (abstract; PIC with collisions);
  `shankar2026particle` (A) (PIC tied to LAPD); `zhang2024kinetic` (abstract; quasi-parallel experiments planned at
  NIF, shocks "more difficult to form in the laboratory and to simulate").
- A fact: the electron results are for perpendicular or weakly magnetized shocks; for quasi-parallel shocks we
  read a plan, no measurement.

### 6.3 Border regimes in parameter space
- Quasi-parallel against oblique: ions are injected up to θ ≈ 45-50° and the self-generated field makes the
  effective inclination ~45° (`caprioli2015simulations`); electron effects at 1 au are mostly at oblique and
  quasi-perpendicular shocks (`literature.md`); ripples at θ_Bn ~45° (`gingell2017mms`, abstract).
- M_A ~10 (`gupta2024electron` (S): growth or stall) and M_A ≲ √(m_i/m_e)/4-√(m_i/m_e)/2 (`shalaby2022mechanism` (S));
  the published PSP values (~4 to 9.1, `inputs.md`) lie across the first.
- Weak shocks: the overshoot vanishes at M ≲ 4 (`caprioli2015simulations`); a critical M_s ≈ 2.3 and no injection
  into DSA from the electron firehose instability alone at weak quasi-perpendicular high-β shocks
  (`kang2019electron`, abstract); 3-D hybrid at high β (`kylivnyk2026hybrid`, abstract). The 1-au case lies here.
- Subluminal against superluminal: θ* = acos(v_sh/c); escape needs v_∥ > (B_⊥/B_0) v_sh (`zekovic2025slams` (S)). A
  run at 0.133 c with B_⊥/B_0 ≳ 10 is on the other side of this border than Saturn's real shock.
- Thermal injection against re-acceleration of seeds (`caprioli2018diffusive`; `raptis2025revealing` (S)).
- Shock-scale against transient-scale: foreshock transients at Earth (`literature.md`), SLAMS as local shocks
  (`tsubouchi2004full`, `johlander2022quasi`, abstracts).

### 6.4 Domain boundaries of the run
- Wall or piston (downstream or upstream frame), injection face, shock-rest-frame relaxation
  (`zhou2025planar`), faux shock (`simon2026faux`), receding upstream boundary (`gupta2024electron` (S)).
- A free-escape boundary against focusing in a diverging field (`john2024effects`, abstract); the maximum energy
  set by the domain (`swisdak2023comparison`, abstract).

## 7. Passages re-read for this file
- `kato2015particle` (S): the section on the energy loss of high-energy particles (their Eq. 4 and the runs with
  N_e 50 to 2667).
- `gupta2024electron` (S): set-up (frame, filter, particle number), the dependence on shock speed, the comparison
  with Park et al., the closing discussion.
- `gupta2025speed` (S): the p_max formula against the runs, the validation statement, the time condition.

## 8. Open, for the owner
- The definition of the injection rate and its region.
- Frame and particle number of the runs, given §5 and §2.5.
- Whether the v_sh/c scan varies the particle number at each speed.
- How H4 is stated with passive seeds.
- Which border regimes of §6.3 belong to the question and which are named as limits.
