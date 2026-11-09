# Implementation plans

One file per non-trivial change: `YYYY-MM-DD-<topic>.md`, written before implementation and
kept as the record of what was planned. Suggested layout:

```markdown
# <topic>
Status: Planned | In progress | Done (PR #N) | Dropped
Refs: roadmap ID, D-IDs

## Context      — the problem and the owner's request
## Approach     — the chosen design (alternatives only if rejected for a recorded reason)
## Invariants   — which contracts in .claude/rules/numerics-invariants.md are touched
## Verification — CPU tests, GPU tests, gates to run, expected evidence
```

Design decisions that outlive the plan go to `docs/design/` (with a D-ID) and `PROGRESS.md`.
