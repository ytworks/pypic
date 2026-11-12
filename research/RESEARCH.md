# Research

This is the cross-session memory of the research branches. Rules: `.claude/rules/research/`.
Layout: `research/README.md`. Engine state: `PROGRESS.md` (read-only here).

- Update this file on `research/trunk` only. A topic keeps its own state in its `README.md` Status and in `notebook.md`.
- Budget: at most 150 lines and 10 Change Log entries (`research/tools/check_branch.py`).
- Older entries move verbatim to `research/archive/change-log.md`.

## Topics

| Topic | Branch | Question | Status | Latest record |
|---|---|---|---|---|
| spherical-shock-injection | `research/spherical-shock-injection` | Which mechanisms, known or unknown, inject electrons at a spherical shock, where on the front (θ_Bn, curvature), and which exist only because the front is curved? (draft) | Purpose stated 2026-11-12; survey open | none |

## Engine baseline

- main 15e5c2c (2026-11-12, merge of PR #42), taken at the creation of `research/trunk`.
- Its `pypic/` tree is that of 300c1a4. Gates all passed there: `bench/results/2026-11-12-300c1a498-gpu-gates.json`, `-mixed.json` and `-compact.json`.

## Change Log

### 2026-11-12 — [progress] Topic spherical-shock-injection opened
- What: `research/spherical-shock-injection` was cut from trunk at 6f00c11.
  - Purpose (owner): the electron injection problem, studied with PIC of a spherical shock. It covers every acceleration mechanism, known and unknown: surfing, drift, wakefield, first- and second-order Fermi, combinations, geometric effects of the curved front, and the front's shape from perpendicular to parallel θ_Bn.
  - The README has the background and a draft question. The survey list is in the topic's `literature.md`.
- Refs: `research/topics/spherical-shock-injection/README.md` (on its branch).

### 2026-11-12 — [amendment] Research branches use no pull requests
- What: research branches use no pull requests at all, not even between research branches. Merges are local `git merge --no-ff` followed by a push.
  - The PreToolUse guard blocks every PR creation, through the GitHub MCP tool or the gh CLI. Its gh check now parses commands, so the words in a commit message no longer trigger it.
  - The CI job `no-pull-requests` fails any PR that carries the research files.
- Why: owner's decision. Nothing on a research branch is merged by a PR, so a PR serves no purpose. PR #44 (research/trunk → main) was closed without a merge.
- Refs: `.claude/rules/research/branching.md`, `research/tools/hooks/guard.py`, `.github/workflows/research.yml`; amends the entry "Research branches and their harness".

### 2026-11-12 — [design] Research reviewer agent; literature notes on Hoshino's work
- What: `.claude/agents/research-reviewer.md` is a reviewer for research questions, plans, results and drafts.
  - It is read-only, and is for research only: never for code, the engine or engineering.
  - It judges work by the standards visible in 14 papers by M. Hoshino and his group, and by the stance. It does not claim to be the author.
  - Notes on those papers are in `research/literature/`.
  - The stance requires a review before every milestone merge into trunk.
- Why: the owner asked for a reviewer modelled on Hoshino, kept apart from code review.
- Refs: `.claude/agents/research-reviewer.md`, `research/literature/hoshino*.md`, `amano*.md`, `guo2023magnetic.md`.

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
