# Research branches

The research branches collect research done with the pypic engine:

- runs recorded so that they can be reproduced from history;
- surveys of prior work.

They are kept apart from engine development, which happens on main. The rules that Claude Code
loads on these branches are in `.claude/rules/research/`.

## Branches

```
main ──(merge, at a gated point)──▶ research/trunk ──(merge)──▶ research/<topic>
                                         ▲                            │
                                         └──── merge --no-ff ─────────┘  (at a milestone)
research/* ──✕──▶ main                   (no PR, push, merge or cherry-pick)
```

- `research/trunk` holds what all topics share: tools, rules, `RESEARCH.md`, the literature.
  It is the only branch that takes main.
- `research/<topic>` holds one research theme. It changes only `research/topics/<topic>/` and `research/literature/`.
- Main is taken only when the engine has changed, and only at a commit whose `pypic/` tree has an
  all-passed gate record in `bench/results/index.json`. The skill `research-sync-main` does this.
- An engine change needed by research is made on a branch from main and goes through its PR and gates.
  Research takes it after the merge to main.

## Why merges of main never conflict

A research branch changes only paths that main does not have:

- `research/`
- `.claude/rules/research/`
- `.claude/skills/research-*/`
- `.claude/agents/research-*.md`
- `.github/workflows/research.yml`

The one exception is the hooks it adds to `.claude/settings.json`. Every other file stays as main has it, so git
takes main's changes to those files as they are. `research/tools/check_branch.py` checks this
against `merge-base(HEAD, origin/main)`.

## Layout

```
research/RESEARCH.md            cross-session memory: topics, engine baseline, change log (trunk only)
research/archive/change-log.md  older change log entries, verbatim
research/literature/            one file per work (<citekey>.md; format in its README)
research/templates/topic/       README.md, notebook.md, literature.md of a new topic
research/topics/<topic>/        README.md, notebook.md, literature.md, configs/, analyze.py, runs/, data/
research/tools/check_branch.py  the branch rules (CI, SessionStart hook, by hand)
research/tools/hooks/           SessionStart summary and PreToolUse guard (.claude/settings.json)
research/tests/                 tests of the tools (pytest research/tests)
```

## Enforcement

- **Rules:** `.claude/rules/research/*.md` load at every session start on a research branch, with no
  change to `CLAUDE.md`. Main has no such files, so they never apply there.
- **PreToolUse hook (`research/tools/hooks/guard.py`):** blocks these:
  - edits of main-owned files;
  - pushes to main, force pushes and branch deletes;
  - rebases;
  - commits of main-owned files outside a merge;
  - pull requests from a research branch into a non-research branch.

  It matches strings, so it is a tripwire, not a proof.
- **SessionStart hook (`research/tools/hooks/session_start.py`):** reports the following. It reads local refs only.
  - the branch and topic;
  - the main commit that was taken and its gate record;
  - how far `origin/main` is ahead;
  - any rule violations.
- **CI (`.github/workflows/research.yml`):** runs `check_branch.py`, `ruff` and the tool tests on
  every push to `research/**`. A pull request that carries the research files into a branch that is not a research branch fails.
  GitHub reads a pull request's workflows from its merge commit, so this check works without any change to main.

## Limits and settings outside the repository

- **Codex:** reads `AGENTS.md` → `CLAUDE.md`, which do not mention these rules; they are Claude Code only.
- **Hooks on main:** they run only where their scripts exist, so they do nothing in a main checkout. A merge of research into
  main done there is caught at the next commit or push (the research files are then in the tree),
  or by the CI of its pull request.
- **Recommended:** a GitHub ruleset on `research/**` that forbids force pushes and deletion.
- **Cloud sessions:** start a session for research with `research/trunk` or `research/<topic>` as the source
  branch. A session started from main does not load these rules.
