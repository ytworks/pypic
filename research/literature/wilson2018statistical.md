# Wilson, Stevens, Kasper, Klein, Maruca, Bale, Bowen, Pulupa & Salem (2018): The Statistical Properties of Solar Wind Temperature Parameters Near 1 au

- Ref: Astrophys. J. Suppl. Ser. 236, 41 (2018); DOI 10.3847/1538-4365/aab71c; arXiv:1802.08585; peer-reviewed.
- Read: 2026-11-13, skimmed (ar5iv text to ~100k characters, via extraction; tables 1–3); topics: shock-alfvenic-turbulence

## Claim
Ten years of Wind data (1995–2004) give long-term distributions of electron, proton and alpha temperatures and
betas near 1 au. Only magnetic obstacles (ICME flux ropes) change them significantly; T_e is the least affected.

## Method
- Wind SWE Faraday cups (p, α), 3DP EESA (electrons), TNR (n_e), MFI (B); 2-min grid, > 10⁶ samples
  (~8.2×10⁵ electron, ~1.1×10⁶ proton).
- Categories: all good data; excluding 239 IP shocks; inside / excluding magnetic obstacles (170 ICMEs).
- No fast/slow wind split; electron core/halo not separated.

## Numbers to compare against
- All data, median (IQR): T_e,tot = 11.9 eV (10.0–14.0); T_p,tot = 8.6 eV (4.7–15.9); T_α,tot = 10.8 eV.
- β_e,tot = 1.09 (0.64–1.99), mean 2.31; β_p,tot = 1.05 (0.54–1.77), mean 1.79; β_α = 0.05.
- T_e/T_p = 1.27 (0.78–2.14); T_α/T_p = 1.94.
- Inside magnetic obstacles: β_e = 0.40 (0.23–0.62); β_p = 0.13 (0.05–0.31); T_e/T_p = 2.06.
- Not tabulated: n_e, |B|, f_pe/f_ce, V_A, d_i.

## Assessment
- The cleanest available 1 au β statistics: β_e ≈ β_p ≈ 1 typical, β ≈ 0.1–0.5 in ICME flux ropes (the
  drivers of most fast-forward IP shocks, per Kilpua et al. 2015 — but shocks propagate in the ambient wind
  ahead of the driver, so the "all data" values apply upstream).
- PIC shock studies at β ~ 0.1–1 per species cover the 1 au range; β_e ≈ 1 with T_e/T_p ≈ 1.3 is the median.
- Long-tail distributions (mean ≫ median for β_e); a single "typical" value hides a factor ~3 spread.
- No upstream-of-shock subset is given; the IP-shock category excludes the shocks rather than characterising
  their upstream.
