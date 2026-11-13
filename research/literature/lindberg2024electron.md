# Lindberg, Vaivads, Amano, Raptis & Joshi (2024): Electron Acceleration at Earth's Bow Shock Due to Stochastic Shock Drift Acceleration

- Ref: Geophys. Res. Lett. 51, e2023GL106612 (2024); DOI 10.1029/2023GL106612; peer-reviewed, open access.
  Note: arXiv:2404.07404 (given in the task brief) is a different paper (Amano et al. 2024, whistler statistics,
  Phys. Plasmas); no arXiv id was found for this GRL letter.
- Read: 2026-11-13, full (publisher PDF text); topics: shock-alfvenic-turbulence

## Claim
The SSDA cutoff-energy scaling of Amano et al. (2020), tested before on one MMS crossing, holds for seven more
quasi-perpendicular bow-shock crossings with different upstream parameters (χ² = 0.55, p ≈ 0.8). The pitch-angle
scattering rate at the cutoff decreases with increasing M_A and with θBn closer to 90°.

## Method
- MMS burst data (FPI, FEEPS, FGM, SCM), electrons averaged over MMS1–3; shocks from the Lalti et al. (2022)
  database, chosen to span θBn and M_A; all show energetic electrons near ramp/foot up to ≳ 60 keV.
- u0, θBn: three mixed methods (Paschmann & Daly 2000), averages over several 45–60 s upstream/downstream window
  pairs outside the burst interval.
- f(E) ∝ E^-p exp(−E/Ec) fitted above 1 keV; Dμμ/Ωce from e-folding times of energetic-electron PSD profiles
  (Eq. 2), threshold D*μμ (Eq. 3) with E_sh = me u0²/(2cos²θBn).
- Test statistic X = Ec,obs/Ec,th (Eq. 4) with Ec,th from Eq. 1; uncertainties from Monte Carlo over Gaussian
  errors of Ec, u0, θBn, Dμμ.

## Numbers to compare against
- Scaling (Eq. 1): Ec = 20 keV (u0/600 km/s)² (cos85°/cosθBn)² (Dμμ/0.03Ωce).
- Table 1 (7 shocks): M_A 5.2–11.1; θBn 76–88°; u0 332–714 km/s; Ec 9–24 keV; Dμμ/Ωce 0.003–0.45; X (MC mean)
  0.3–2.5 with σ_X 0.3–6.1 (a second X column, 0.2–1.8, is unlabelled in the extracted text).
- Spectral index p = 3.4 ± 0.2 on average (Sect. 3), consistent with Oka et al. (2018) and Gosling et al. (1989).
- Two groups: shocks 1–2 (u0 ≈ 700 km/s, M_A ≈ 11) show energetic flux over much of the STL; shocks 3–7 only
  near the ramp, close to noise (Sect. 4).
- Lowest M_A/cosθBn (my arithmetic): shock 5, 5.2/cos76° ≈ 21; highest: shock 2, 11.1/cos87° ≈ 210.
- Whistler critical Mach number used: M_wh = 0.5 (mi/me)^½ |cosθBn|; all shocks have M_A ≳ M_wh (Sect. 4).

## Assessment
- Shown: the scaling is consistent with 7 crossings within large errors; σ_X is often larger than X, so the
  test has weak power (several X are consistent with both 1 and 0).
- The M_A/θBn trend of Dμμ rests on 8 shocks; the authors note observational bias cannot be ruled out.
- FPI data compression at 10–27 keV removes points where the cutoff sits, inflating Dμμ errors.
- Upstream conditions are characterised only by u0, θBn, M_A; no upstream turbulence level, β or δB/B is
  reported, so the role of pre-existing solar-wind fluctuations is not addressed.
- Relevance: gives the bow-shock baseline (M_A ~ 5–11, θBn ≳ 75°, Ec ~ 10–20 keV, p ≈ 3.4) against which an
  injection change due to upstream Alfvénic turbulence would have to be measured.
