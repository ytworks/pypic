# Notebook — spherical-shock-injection

Dated lab notes, newest first. Write down what was tried, what was seen, and what was decided, failed attempts included.
Never rewrite a past entry. To correct one, add a new entry that refers to it.

## 2026-11-12 — Citation check of the survey: 16 corrections
- An independent check compared about 290 located statements of Sects. 1–8 of the survey with the paper texts.
  About 280 were confirmed. Three were wrong and 13 imprecise or mis-located; all 16 are corrected in the
  English and Japanese surveys, and one wording in the summaries.
- The three wrong statements:
  - "only one paper gives Mach numbers for SNR shocks": four others quote values (M_A ≈ 20–1000);
  - "K_ep is degenerate with the magnetic field": the Tycho paper says the X-ray flux fixes K_ep independently
    of the field (the note `morlino2012strong` is corrected too);
  - "the 2-D runs last 6–52 Ω_ci⁻¹": two 2-D studies of parallel shocks run longer.
- Not checked: statements that rest on `hoshino2002nonthermal`, `amano2010critical` and `amano2022theory`
  (no text at hand), and three values that are garbled or only in figures; the survey marks them.

## 2026-11-12 — One-page summary of the survey, in English and Japanese
- Asked by the owner. `survey-snr-electron-injection-summary.md` and its `.ja.md` give the whole picture on one
  page: the problem, what observations fix, the mechanisms by θ_Bn, what decides the answer in a simulation,
  what is supported, disputed and unknown, and what it means for this topic.
- Every statement is taken from the survey, with its citekeys.

## 2026-11-12 — Japanese version of the survey
- Asked by the owner. `survey-snr-electron-injection.ja.md` has the same sections as the English survey, which
  stays the original.
- Its reference lists (Sect. 10) are taken from the English file by script; only headings and remarks are
  translated.
- Checked by script, section by section: the citekeys are identical, and no number of the English version is
  missing.

## 2026-11-12 — Survey: electron injection at supernova remnant shocks
- Asked by the owner: survey the prior work on the electron injection problem at SNR shocks, in markdown, with
  the citations kept.
- Done:
  - arXiv searches by keyword (12 queries, 238 new results), then by author and title, then reference lists;
  - 54 works selected, full texts from arXiv, read by nine readers in groups (observations; surfing line;
    perpendicular SNR series; foreshock and reviews; quasi-parallel; thresholds 2021–2026; shock drift theory;
    geometry);
  - one note per work in `research/literature/`;
  - the headline numbers of 20 works checked against the paper texts;
  - the survey written: `survey-snr-electron-injection.md`.
- Main findings:
  - mechanisms are sorted by obliquity, and each is known from planar runs at one or a few angles;
  - no study covers the range of θ_Bn with one parameter set;
  - global models of remnants assume the efficiency law ε(θ_Bn);
  - SN 1006 points to quasi-parallel acceleration and the radio population to quasi-perpendicular;
  - no kinetic simulation of a curved front at SNR-like parameters was found.
- Not done:
  - the bow-shock works on curved shocks (not on arXiv) are not read, so the last finding is not certain;
  - relativistic shocks and wakefield acceleration were not covered;
  - the survey has not been reviewed by the research reviewer.
- Corrected along the way: the reference of `amano2022theory` (ApJ 927, 132).
- No decision on the research question; it stays with the owner.

## 2026-11-12 — First review (purpose stage)
- Review: `reviews/2026-11-12-purpose.md`. Verdict:
  - "every mechanism over all θ_Bn" needs rethinking (a programme, not one question);
  - "curvature-only channels" needs the controlling number R/L_s and a target system first;
  - "the switch-on angle θ_c and its scaling with m_i/m_e and β_e" is answerable.
- Fixed in the README: the Background said that known mechanisms are framed only for planar
  shocks. amano2022theory Sect. 4.5 (Eq. 22) already has a curvature limit.
- Added to `literature.md`: the curved-shock works the review found by search (not read).
- Points to keep whatever the question becomes:
  - the mass ratio moves every threshold, with different powers;
  - θ_Bn is to be taken from the measured front normal;
  - the superluminal band;
  - Buneman waves against resolution and the current filter;
  - tail statistics per θ_Bn band;
  - tracked sets fixed at the start.
- Open: the owner's decision on the research question.

## 2026-11-12 — Topic opened
- The owner set the purpose: the electron injection problem, studied with PIC of a spherical
  shock, considering every acceleration mechanism, unknown ones included.
  - Mechanisms named: surfing, drift, wakefield, first- and second-order Fermi, combinations,
    geometric effects of the curved surface, and the shape changes as θ_Bn goes from
    perpendicular to parallel.
- Branch cut from research/trunk at 6f00c11.
- README: purpose, background and a draft research question. Prior work lists the four notes
  read so far; the survey list is in `literature.md`.
- Open: the survey, a narrowed question, the idea, and the driver and its parameters.
  `examples/blast_wave` (M_A = 3) is likely below the Mach numbers that surfing needs.
