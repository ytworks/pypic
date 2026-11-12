# Kumar & Reville (2021): Non-thermal particle acceleration at highly oblique non-relativistic shocks

- Ref: Astrophys. J. Lett. 921, L14 (2021); DOI 10.3847/2041-8213/ac30e0; arXiv:2110.09939v1; peer-reviewed
  Letter (journal reference from the arXiv metadata; comment: accepted by ApJL)
- Read: 2026-11-12, full; topics: spherical-shock-injection

## Claim
A long 1-D PIC run of a 75 deg, M_A = 63 subluminal shock at m_i/m_e = 256 shows electrons pre-heated in the shock
foot to gyroradii above the shock thickness and then both electrons and ions starting diffusive acceleration;
the ion tail appears only after the electrons have reached ion-scale gyroradii.

## Method
- PIC, SMILEI, 1D3V; shock-normal frame; relaxation set-up: a Rankine-Hugoniot discontinuity of width 2 rho_i with
  a tanh profile at t = 0, downstream T_i/T_e = 4 imposed (other ratios gave comparable early results).
- m_i/m_e = 256; test-electron comparison runs at 100 and 25 (Fig. 4).
- Upstream drift = shock speed 0.15 c. M_A = 63; M_s = 70 with M_s = v/sqrt((T_e+T_i)/m_i); T_i = T_e upstream.
- beta_e = beta_i ~ 1.8 as stated ("intermediate regime beta ~ 1"); see the doubt under Assessment.
- theta_Bn = 75 deg, field in the x-z plane. omega_ce/omega_pe = v sqrt(m_i/m_e)/M_A (stated relation; = 0.038).
- Box ~570 rho_i, rho_i = (v/B_d)(m_i/m_e); dx = 0.1 c/omega_pe; dt = 0.095/omega_pe; 512 (upstream) and 1024
  (downstream) particles per cell per species; energy conserved to better than 1 %.
- Times are in 1/omega_ci defined with the DOWNSTREAM field. End time not stated in the text (events at
  omega_ci t ~ 30 and ~50 are mentioned). 1000 test electrons per group follow the self-consistent fields.

## Numbers to compare against
- Shock reformation at about the downstream ion gyrofrequency; for omega_ci t > 30 sharp electrostatic and
  electromagnetic structures run upstream from the foot (Sect. 3, Fig. 1).
- Electrons reach p_e ~ 10-40 m_e c, where their gyroradius exceeds the shock thickness; a distinct population
  above that appears at omega_ci t ~ 50, described as possibly transient (Sect. 3, Fig. 2).
- Downstream spectra summed over 2-16 rho_i behind the shock, shown as p^2 times the distribution; "flat" tails;
  no fitted index is quoted (Fig. 2).
- Efficiencies at the end: eta_e ~ 45 %, eta_i ~ 5 %, defined as kinetic energy of non-thermal over thermal
  particles downstream, non-thermal meaning p_e > 4 m_e c and p_i > 90 m_e c (Sect. 3, Fig. 2 insets). Relative
  to the incoming kinetic energy the electron value "would remain small"; eta_i is not saturated.
- Electron energisation slows above p_e ~ u_sh m_i/m_e; beyond it the electron and ion cut-offs advance together
  (Sect. 4, Fig. 3).
- Wave condition quoted: M_A >= (cos theta / 2) sqrt((m_i/m_e) beta_e) (Amano & Hoshino 2010) for Buneman and
  oblique whistler modes; also M_A,crit = sqrt(m_i/m_e) (White et al. 2020). The run exceeds both (Sect. 4).
- Mass ratio: test electrons gain more energy at 256 than at 100 or 25, and pitch-angle scattering is clearest
  at 256; electrons arriving later gain more than those arriving early (Fig. 4).
- Geometry remarks: mean obliquity 60 deg for random orientation; superluminal if |cos theta| < v_sh/c (Sect. 1).

## Assessment
- Shown by simulation (1-D): non-thermal electrons and, later, ions at 75 deg; the most energetic ions make long
  upstream excursions with diffusive guiding-centre motion (Fig. 3); pre-heating grows with mass ratio.
- Injection path as described: electrons meeting the reforming shock at a favourable phase are reflected, then
  heated by high-frequency waves in the foot that the reflected ions drive, and by the cross-shock potential.
  The authors say the first energisation cannot be assigned to a single instability.
- Not established (their words): whether electron acceleration starts the ion acceleration or the reverse.
  Speculation: the precursor keeps growing until transport near the shock no longer depends on the field angle.
- Limits stated: 1-D, forced by the particle number needed against numerical drag; early evolution depends on
  the launching method; a 1-D box in another frame samples another wave vector; no parameter survey; ion
  efficiency not saturated; run times too short for higher energies.
- Disagreement: earlier work at theta > 60 deg (Riquelme & Spitkovsky 2011; Caprioli & Spitkovsky 2014; Xu et al.
  2020) found no significant ion acceleration; here it sets in late. Intermediate-time results are said to
  resemble Xu et al. (2020), but that paper's electron-driven upstream waves are not analysed here.
- The literature disagrees on the whistler Mach-number threshold (Amano & Hoshino 2010 against Riquelme &
  Spitkovsky 2011, who find whistlers at low M_A); the paper notes this.
- Doubt: with the stated definitions and M_A = 63, M_s = 70, the relation beta_e = 2 (M_A/M_s)^2 gives 1.62, not
  1.8, and the definition printed for beta_e (4 pi n k T / B^2) lacks the usual factor 2.
