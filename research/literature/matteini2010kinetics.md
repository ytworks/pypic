# Matteini, Landi, Velli & Hellinger (2010): Kinetics of parametric instabilities of Alfven waves: Evolution of ion distribution functions

- Ref: J. Geophys. Res. Space Phys. 115, A09106 (2010) (Crossref); DOI 10.1029/2009JA014987; no arXiv
  found; peer-reviewed.
- Read: 2026-11-13, abstract; topics: shock-alfvenic-turbulence. The preprint PDF (asu.cas.cz) was
  fetched but could not be converted to text here.

## Claim
In hybrid simulations, large-amplitude Alfven waves undergo modulational and decay instabilities whose
evolution and saturation are changed by ion kinetics: resonant protons are trapped, accelerated along B,
and in some cases form a parallel beam faster than the core. The beam is driven by the parallel electric
field of the density fluctuations of the ion-acoustic daughter modes.

## Method
- Hybrid code (kinetic protons, fluid electrons); per the search snippet of the full text, 1-D runs
  with the Matthews (1994) hybrid code.
- Monochromatic and incoherent (spectrum) pumps; several polarizations and amplitudes; beta < 1 to
  beta > 1 (abstract). Box, resolution, amplitude values: not read.

## Numbers to compare against
- None read (growth rates, saturation times, beam speed are in the full text only).
- Companion: Matteini et al. 2010 GRL (DOI 10.1029/2010GL044806, linearly polarized oblique waves,
  1-D/2-D hybrid): decay of the parallel component into an ion-acoustic wave plus a backward Alfven
  wave, then broadband oblique compressive spectra (search snippet only, not a separate note).

## Assessment
- Shows (per abstract) that an Alfvenic pump of large amplitude is not a stationary background in a
  kinetic plasma at any beta in their range: ion-acoustic daughters with E_par develop. For us this is
  the main reason an imposed |B| = const upstream will evolve on its own, and a reason why a "turbulence
  without shock" null run is required over the same time window.
- Hybrid: says nothing on electron heating or electron acceleration; electrons are a fluid.
- Open for our design: the growth time at our beta and delta B/B0 (0.3-1) must be taken from linear
  theory (Derby 1978; Goldstein 1978) or from gonzalez2023particle, not from this abstract.
