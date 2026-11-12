# Research branches

These rules exist only on research branches: `research/trunk`, `research/<topic>`, and session branches
cut from them. Research mode means `research/RESEARCH.md` exists. Layout and procedures: `research/README.md`.

## Branch model
- `research/trunk` holds the shared research base: tools, rules, `research/RESEARCH.md`, the
  literature. It is the only branch that takes main.
- `research/<topic>` (slug `[a-z0-9-]+`) holds one research theme. It is cut from trunk and takes trunk.
  Start one with the skill `research-new-topic`.
- A session branch (e.g. `claude/*`) cut from a research branch belongs to that research branch.
  Land its commits there, and nowhere else.
- Changes flow main → trunk → topic, always with `git merge --no-ff`. A topic merges into trunk
  (`--no-ff`) at a milestone: a recorded result or a finished survey.
- Nothing on `research/*` ever reaches main: no push, merge or cherry-pick.
- Research branches use no pull requests, not even between research branches. Integrate with a
  local `git merge --no-ff` and a push.
- Engine changes go the other way: a branch from main, its PR and gates, the merge to main. Then
  take main into trunk with the skill `research-sync-main`, only at a gated point.
- Never rebase, squash, amend a pushed commit, or force-push a research branch. Never delete a
  research branch that has commits not merged into trunk. Run records name commits, and those
  commits must stay reachable.

## What a research branch may change (keeps every merge of main conflict-free)
- Only these paths: `research/**`, `.claude/rules/research/**`, `.claude/skills/research-*/**`,
  `.claude/agents/research-*.md`, `.github/workflows/research.yml`, and the research hooks in
  `.claude/settings.json`. Every other file stays as main has it. `research/tools/check_branch.py`
  checks this, and the PreToolUse hook enforces it.
- A topic branch changes only `research/topics/<topic>/**` and `research/literature/*.md`. Shared
  research files change on trunk; the topic then takes trunk.
- If a merge of main ever conflicts, take main's side. The exception is `.claude/settings.json`: keep the hooks of both sides.

## Main's rules on a research branch
- The engine contracts (`numerics-invariants.md`, the test policy) hold unchanged.
- These paths are read-only here: `pypic/`, `examples/`, `scripts/`, `tests/`, `docs/`, `bench/`,
  `PROGRESS.md`, `CLAUDE.md`, `AGENTS.md`.
- Instead of `PROGRESS.md`, update `research/RESEARCH.md`, and only on trunk. A topic keeps its state
  in the Status section of its `README.md` and in `notebook.md`.
- Instead of `bench/results/`, record runs in `research/topics/<topic>/runs/`.
- At session start, read `PROGRESS.md` (the engine), then `research/RESEARCH.md`, then the topic's `README.md`.

## Runs and literature
- A recorded run starts from a clean, committed tree, `research/` included. Its record holds the
  commit, config, seed, GPUs and the analysis result.
- Big outputs stay out of git: the record holds their path and size. Every file under `research/` is at most 1 MB.
- Literature: one file per work, `research/literature/<citekey>.md` (format in that folder's
  README). A topic's `literature.md` places the topic against those works.
