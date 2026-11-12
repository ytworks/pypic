# Wang & Amano (2026): Generation of Whistler Waves by Reflected Electrons and Their Self-Confinement at Quasi-Perpendicular Shocks

- Ref: Physics of Plasmas, accepted (arXiv comment); volume/article not in the metadata, journal reference not
  checked; DOI 10.1063/5.0317681; arXiv:2603.19018 (v1, read here); accepted after peer review per the comment
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
Electrons mirror-reflected at a quasi-perpendicular shock drive two whistler instabilities at the upstream edge of
the transition layer: a downstream-directed mode through the normal cyclotron resonance at the loss cone (DCW) and
an oblique upstream-directed mode through the anomalous cyclotron resonance of the beam (UAW). Both need a high
HTF Alfvén Mach number MA/cosθBn and are helped by a high βe; for bow-shock parameters the threshold is about 50.
Such shocks exceed the whistler critical Mach number, so the waves stay in the layer, and scattering by them opens
further whistler instabilities. The authors suggest this self-confinement supplies the scattering SSDA needs.

## Method
- Linear theory only. Electron distribution by Liouville mapping in the HTF through a steady 1-D shock: tanh
  profile for Bz (Eq. 6) and for the potential (Eq. 7), compression r = 4, normalized HTF potential
  2eΔΦ/(mi V1²) = 0.05 (Eq. 8); adiabatic invariants; causality rule for the doubly accessible region (Eq. 13).
- Upstream drifting Maxwellian; downstream Maxwellian from thermalized transmitted electrons (Eqs. 14–16).
- Ad hoc pitch-angle diffusion in the local plasma frame with Dμμ t as a free parameter (Eqs. 17–19).
- Growth rate from the Kennel–Wong weak-growth formula (Eqs. 1–5) with cold-plasma whistler dispersion;
  resonances n = −1, 0, +1 of electrons only; evaluated at x = −3 (upstream edge).
- Parameters: mi/me = 1836, Ωce/ωpe = 0.01 (ωpe/Ωce = 100), non-relativistic. Baseline MA = 10 (normal incidence
  frame), βe = 2, θBn = 85° (Sect. III D). Scans: Dμμ t = 0.01, 0.05, 0.2, 0.5 (Fig. 5); θBn = 75°, 80°, 87° at
  Dμμ t = 0.1 (Fig. 6); βe = 1, 3, 4 at θBn = 85°, Dμμ t = 0.1 (Fig. 7). No simulation, no new data.

## Numbers to compare against
- HTF flow speed for the baseline: U1/VA,e ≈ 2.68 (Sect. III D). Resonance speeds: ≈ −6 VA,e for the DCW and
  ≈ −4 VA,e for the UAW (Sect. IV A).
- Dμμ t = 0.05 removes the DCW; the UAW survives to Dμμ t = 0.2 with a reduced rate (Sect. IV B, Fig. 5).
- Secondary modes: Landau-driven (ULW) near k∥c/ωpe ≃ −0.2, k⊥c/ωpe ≃ 0.2 at Dμμ t = 0.2; cyclotron-driven
  (UCW) near k∥c/ωpe ≃ −0.3, k⊥ ≃ 0 at Dμμ t = 0.5 (Sect. IV B, App. A, Table I).
- Beam speed v_b/VA ~ 2 MA/cosθBn (Eq. 21); UAW condition v_b ≳ VA,e/2, i.e. MA/cosθBn ≳ (α/4)(mi/me)^(1/2)
  (Eq. 22). Requiring γ/Ωci ≳ 10 (γ/Ωce ≳ 0.01) makes 75° and 80° marginal and 85° unstable at MA = 10, hence
  α ~ 5 and MA/cosθBn ≳ 50 (Sect. IV C).
- Not monotonic in obliquity: at 87° the growth rate drops because fewer electrons are reflected (Sect. IV C).
- Reflected density: n_b ∝ exp[−(MA/cosθBn)² (me/mi) / ((B2/B1 − 1) βe)] (Eq. 24), from
  v⊥,min = [(v∥² + 2eΔΦ/me)/(B2/B1 − 1)]^(1/2) (Eq. 23); growth increases from βe = 1 to 4 (Fig. 7).
- Whistler critical Mach number (HTF) (mi/me)^(1/2)/2 ≈ 21 (Sect. V A; printed "2158", i.e. 21 and citation 58).
- MMS statistics quoted: SSDA wave-power threshold met for MA/cosθBn ≳ 30–60 (Sect. V B, refs. 21, 11, 61).

## Assessment
- Shown (linear theory for a model distribution): two primary and two secondary unstable whistler branches with
  their driving resonances; trends with θBn, βe and the amount of isotropization.
- The threshold 50 rests on a chosen growth criterion and one MA; α is fitted, not derived. Dμμ is never
  estimated, so confinement and the sequence of instabilities are suggestions (Sect. V C).
- Stated/conjectured: wave accumulation and self-confinement; positive feedback between acceleration and waves;
  one-hertz foreshock whistlers as the sub-threshold case; SNR shocks pass the threshold at lower obliquity;
  reflection may fail at very high MA/cosθBn but Buneman/Weibel heating should help; high-βe cluster shocks
  favourable at low Mach number (Sect. V B).
- Stated limits: steady monotonic profile (no overshoot, reformation or rippling); potential a free parameter,
  so electron bulk flow and assumed frame differ deep in the layer; cold-plasma dispersion; scattering does not
  feed back on transport; no saturation level. PIC tests need a mass ratio and ωpe/Ωce that keep the beam above
  VA,e yet non-relativistic; existing PIC work on reflected-electron whistlers used relativistic beams (Sect. V B).
- Relation to other work: the DCW is the loss-cone mechanism of Amano & Hoshino (2010) with a fuller
  distribution, and is the more fragile mode here; the UAW is tied to Levinson (1992, 1996). Here a higher βe
  raises the growth rate (more reflected electrons), whereas the Amano & Hoshino (2010) threshold rises as
  (βe)^(1/2) (thermal cyclotron damping). The scans stay above that threshold and the paper does not compare
  the two trends. Minor: the text announces four βe values and lists three.
