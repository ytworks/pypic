# F. Guo, Y.-H. Liu, S. Zenitani, M. Hoshino (2023): Magnetic Reconnection and Associated Particle Acceleration in High-energy Astrophysics

- Ref: Space Science Reviews, "in press" per the arXiv v2 comment (Apr 2024; ISSI 2022 reconnection workshop collection); arXiv:2309.13382. Journal volume and year not checked.
- Read: 2026-11-12, skimmed; topics: (none yet)

## Claim
A review of reconnection in the magnetically dominated (σ ≫ 1) regime. Relativistic collisionless reconnection is fast
(R ~ 0.1–0.3, bounded near 0.3 by outflow geometry) and gives power-law spectra that approach p ~ 1 at large σ. The
mechanisms of power-law formation, injection and energy partition are still debated. 3D effects, onset, multiscale
coupling and radiation are the frontier problems.

## Method
- A literature synthesis. Sec. 1: astrophysical sites and key issues (1.1, 1.2). Sec. 2: relativistic resistive MHD and
  two-fluid models (Petschek, Sweet–Parker, resistivity models, shocks). Sec. 3: collisionless PIC (relativistic Ohm's
  law, rate theory, bursty plasmoids, 3D turbulent reconnection). Sec. 4: heating and acceleration (mechanisms, spectra,
  injection, roll-over, 3D, scale separation). Sec. 5: outlook.
- Hoshino-related parts: accretion disks (Sec. 1.1, citing Hoshino 2015); the relativistic-drift stabilization of tearing
  (Sec. 3.3, Hoshino 2020); the Zenitani–Hoshino 2001 power-law model (Sec. 4.3.1); thermal/nonthermal energy partition
  (Sec. 4.3.2, Fig. 16, Hoshino 2022, 2023); drift-kink and guide field (Sec. 4.4).
- Read fully: Secs. 1, 3.2–5. Skimmed: Secs. 2 and 3.1.

## Numbers to compare against
- σ = B²/(4πw) and vA/c = √(σ/(σ+1)) (Eq. 1, Sec. 1; Eq. 26).
- Maximum plausible rate ≈ 0.3 relativistic and ≈ 0.2 nonrelativistic (Sec. 3.2.1, Fig. 10). 3D turbulent PIC runs stay
  below 0.3 (Sec. 3.4, Fig. 13).
- X-line pressure bound ≈ 2(2/σx0)^{1/2} of the asymptotic magnetic pressure (Eq. 30, Fig. 12a).
- Relativistic Sweet–Parker gives R ∝ S^{-1/2}. Above S ≳ 1e4 the sheet is plasmoid-dominated with R ~ 0.01 (Sec. 2.4, Fig. 6).
- The outflow is supersonic for σ > 1/2 (Sec. 2.6).
- Zenitani–Hoshino power law N ∝ E^{-2B̄z/(πErec)} with index O(1) (Eqs. 32–34). Fokker–Planck index
  p = 1 + 1/(αacc τesc) (Eqs. 35–36).
- Indices: p ~ 1 for pair plasma at high σ, p ≳ 2 for electron–proton plasma at σi, p ≳ 3–4 nonrelativistic (Sec. 4.2).
- Energy partition (Hoshino 2023, Fig. 16, Sec. 4.3.2): with no guide field, nonthermal particles hold more than 90 % of the
  particle energy in relativistic runs, while nonrelativistic runs mostly heat. The nonthermal fraction falls with guide field,
  except in the weak-guide-field nonrelativistic case.
- E > B regions supply only ~ 10 % of γinj (Sec. 4.3.3, from Guo et al. 2023). Roll-over γc ~ [2δσ(2 − p)]^{1/(2−p)} (Sec. 4.3.4).
- Scale separation L/de ~ 1e13–1e17 for PWNe and jets (Sec. 4.5).

## Assessment
- As a review, its numbers come from the cited works and were not checked here. The authors themselves mark several
  disputes: X-line vs Fermi origin of the power law (4.3.1), which mechanism injects (4.3.3), and the claim that escape is
  needed for a power law, called "incorrect or at least misleading" (4.3.1).
- Open questions named in Sec. 1.2 and Sec. 5:
  - the rate problem;
  - partition, injection, power-law formation and roll-over;
  - whether reconnection in 3D runs through kinetic-scale layers or broad resistive ones, and whether 3D acceleration
    differs fundamentally from 2D;
  - onset and current-sheet formation, including whether a Harris sheet is a valid initial condition;
  - multiscale coupling (no relativistic MHD-PIC coupling exists yet);
  - radiation and QED effects;
  - observational links (minijets, Crab flares, polarization signatures surviving in 3D).
- Stated limits: the energy-partition studies are 2D (4.3.2); PIC conclusions need extrapolation across a huge scale gap
  (4.5); GR effects on local reconnection physics are unclear (2.7); there are few relativistic resistive-MHD studies (2.7).
- Citekey caveat: arXiv year used. If the published version is dated 2024, rename to guo2024magnetic.
