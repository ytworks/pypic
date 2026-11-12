# Notebook — spherical-shock-injection

Dated lab notes, newest first. Write down what was tried, what was seen, and what was decided, failed attempts included.
Never rewrite a past entry. To correct one, add a new entry that refers to it.

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
