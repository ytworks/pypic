# Caprioli, Zhang & Spitkovsky (2018): Diffusive shock re-acceleration

- Ref: J. Plasma Phys. 84 (3), 715840301; DOI 10.1017/S0022377818000478; arXiv:1801.01510. Peer reviewed.
  Read in the arXiv version.
- Read: 2026-11-14, arXiv PDF text: abstract, set-up, the results on reflection, the current of reflected
  seeds and the efficiencies, the sections on heliospheric shocks and on electrons, conclusions; topics:
  qpar-electron-conditions

## Claim
Pre-existing energetic particles (seeds) are reflected and re-accelerated at a shock of any inclination.
Their current drives the Bell instability, which "can eventually trigger the injection of thermal protons
even at oblique shocks that ordinarily cannot inject thermal particles". "the acceleration efficiency is
not an intrinsic property of the shock, but rather it is regulated by the amount of seeds available; this
means that the level of self-generated magnetic turbulence is not universal, either, which may in turn
limit the maximum energy achievable".

## Method
- 2-D hybrid, ions only, with an isotropic seed population upstream. Benchmark: θ 60°, M 30, seed speed
  50 v_A, n_CR = 0.01; runs at θ 0-80°.
- The seed density is "much larger than those expected in the interstellar medium or in the solar wind".

## Numbers to compare against
- Current of reflected seeds J_CR ≃ e n_CR v_sh for v_CR ≫ v_sh; "no threshold for the effect to be
  relevant", growth limited only by damping.
- With seeds, thermal-proton efficiency ~3 % against ≲ 0.5 % without (M 30, 60°).
- Without seeds, box width 1000 c/ω_p: "no sign of non-thermal activity or magnetic field amplification up
  to t ≈ 1000 ω_c⁻¹".
- The maximum energy of the seeds grows linearly in time. The standard DSA slope fails at θ ≳ 70°.

## Assessment
- Heliospheric shocks: "the level of pre-existing seeds in heliospheric shocks varies greatly from event to
  event, so it is impossible to provide a universal estimate for their re-acceleration efficiency"; the
  injection criterion and the seed phenomenology "should suffice for interpreting a given heliospheric
  event once the far upstream conditions (shock strength and inclination, seed abundance and chemical
  composition) are measured". Extrinsic turbulence "would need to be assessed on a case-by-case basis".
- Electrons: "the same re-acceleration mechanisms should apply also to … seed electrons". SN 1006: the
  bilateral shape "favors a scenario where ion DSA is efficient only in quasi-parallel region, while
  electron acceleration (or re-acceleration) can occur regardless of the shock inclination".
- "the main results in this work do not depend on the prescription for the electron equation of state".
- For this topic: the simulation basis of candidate H4 (seeds), for ions. Two consequences. (1) If seeds
  set the efficiency, the answer to our question is not a property of the shock parameters alone: the seed
  abundance is an input that the events must supply (Saturn's ~100 keV electrons, the keV electrons at
  Earth and 1 au). (2) Our seeds are a passive test species (PIC-D9), which cannot drive waves. That is
  right as long as the seed electrons carry a negligible current; the run shows how seeds are accelerated
  in the shock's own turbulence, not what turbulence the seeds make. The limit must be named with H4.
