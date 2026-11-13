# Liu, Song, Cao, Liu, Li, Russell, Pollock, Lindqvist, Burch, Lu, et al. (2026): Efficient acceleration of energetic electrons upstream of Earth's bow shock

- Ref: Nat. Commun. 17, 5546 (2026); DOI 10.1038/s41467-026-72197-y; published online 2026-04-22; peer-reviewed,
  open access. Author list as in the page metadata (first ten read; first author Y. Y. Liu).
- Read: 2026-11-13, full (main text and Methods; Supplementary not read); topics: shock-alfvenic-turbulence

## Claim
When a solar-wind rotational discontinuity (RD) intersects the bow shock, field lines take a V shape with both
ends on the shock. Shock-reflected electrons bounce between the two converging field-line–shock intersections
and gain parallel energy by repeated fast-Fermi reflections (conservation of J = m v∥ L). MMS saw a bidirectional
1–30 keV electron burst inside an RD upstream of the shock, raising ~17 eV solar-wind electrons to > 200 k_B T_e.

## Method
- MMS-1, 2019-01-05; shock crossing 14:23 UT at (13.9, 0.3, 5.5) R_E GSE; RD with burst at 15:28:10 UT.
- RD identified from Rankine–Hugoniot jump conditions (Eqs. 1–5); normal (0.64, 0.61, −0.47).
- Field-line tracing to a Shue-model shock gives the magnetic connection state (forward/backward/both).
- Test particles (Boris) in a 2.5-D tanh model of RD + planar shock (Eqs. 11–12), E = −V×B.
- Liouville mapping of the measured upstream distribution with E∥2/E∥1 = (L1/L2)² (Eq. 10).

## Numbers to compare against
- Upstream: |B| = 4.3 nT, n = 4 cm⁻³, flow ≈ 520 km/s; downstream ≈ 190 km/s; B_msh = 25 nT (loss cone ≈ 25°);
  model compression 5.8. Test-particle run uses V_sw = 332 km/s (Methods) — differs from the 520 km/s quoted.
- RD rotation 96°; field 71.5° from RD normal; half-thicknesses L_RD = 801 km, L_shock = 143 km (Methods).
- Test electron 500 eV, PA 50°: 18 shock reflections in ~23 s, 13,262 → 24,095 km/s.
- Liouville mapping: f∥ = 3.66, i.e. field lines shortened to 52 % (Fig. 4c); fit good for PA < 140°.
- RD occurrence near Earth ≈ 17 per hour (their ref. 44).
- Derived (my arithmetic): d_i = 228 km/√4 ≈ 114 km; L_RD ≈ 7 d_i; f_pe/f_ce ≈ 18 kHz / 120 Hz ≈ 150.
- θBn and M_A of the local shock are not given in the text read.

## Assessment
- Shown: one event with consistent field-line topology, pitch-angle signatures and a Liouville reconstruction.
  A second, weaker candidate at 14:33:45 UT is mentioned, not analysed. No statistics.
- The mechanism is an upstream, large-scale (≫ d_i) field-line geometry effect of a pre-existing solar-wind
  structure; it needs neither scattering nor shock-layer waves. Maximum energy limited by the loss cone; the
  authors argue scattering in the ramp may extend it (not shown).
- Directly relevant: large-amplitude Alfvénic turbulence makes RD-like rotations ubiquitous; whether the same
  converging-intersection trapping operates for a turbulent upstream field at an IP shock is open. A 2-D/3-D PIC
  box needs field-line excursions of order the box to contain both intersections.
