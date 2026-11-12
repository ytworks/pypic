# Research

This is the cross-session memory of the research branches. Rules: `.claude/rules/research/`.
Layout: `research/README.md`. Engine state: `PROGRESS.md` (read-only here).

- Update this file on `research/trunk` only. A topic keeps its own state in its `README.md` Status and in `notebook.md`.
- Budget: at most 150 lines and 10 Change Log entries (`research/tools/check_branch.py`).
- Older entries move verbatim to `research/archive/change-log.md`.

## Topics

| Topic | Branch | Question | Status | Latest record |
|---|---|---|---|---|
| (none yet) | | | | |

## Engine baseline

- main 15e5c2c (2026-11-12, merge of PR #42), taken at the creation of `research/trunk`.
- Its `pypic/` tree is that of 300c1a4. Gates all passed there: `bench/results/2026-11-12-300c1a498-gpu-gates.json`, `-mixed.json` and `-compact.json`.

## Change Log

### 2026-11-12 — [design] Research stance in the harness
- What: `.claude/rules/research/stance.md` is now loaded on every research branch. It summarises the owner's research stance (A. Morihata, "研究とは何か") in our own words, with a link to the source.
  - The topic README sections follow its elements: Background, Research question, Prior work, Idea, Related work and novelty, Approach, Claims and evidence, Publication, Status.
  - `check_branch.py` checks those headings.
  - `.claude/agents/research-*.md` is added to the research paths, for the research reviewer.
- Why: the owner asked for a stance that holds in every session without being restated.
- Refs: `.claude/rules/research/stance.md`, `research/templates/topic/README.md`.

### 2026-11-12 — [design] Research branches and their harness
- What: `research/trunk` was cut from main 15e5c2c. It is never merged into main and takes main only.
  - Topics are `research/<topic>`.
  - The rules are `.claude/rules/research/`, loaded automatically on research branches only.
  - Enforcement: the hooks in `.claude/settings.json` (`research/tools/hooks/`), `research/tools/check_branch.py`, and `.github/workflows/research.yml`.
  - Skills: `research-new-topic`, `research-sync-main`.
- Why: the owner wants research runs kept reproducible and prior work surveyed, apart from engine development.
  - Every merge of main must stay conflict-free. So a research branch changes only paths that main does not have, plus the hooks in `.claude/settings.json`.
  - Run records name commits. So research history is never rewritten.
- Refs: `.claude/rules/research/branching.md`, `research/README.md`.
