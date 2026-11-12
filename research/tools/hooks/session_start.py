#!/usr/bin/env python3
"""SessionStart hook of the research branches (.claude/settings.json): a few lines of state for a
session that starts in research mode (research/RESEARCH.md in the tree); silent otherwise. Reads
local refs only (no fetch) and never fails the session.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import check_branch as cb  # noqa: E402


def report(root: Path) -> list[str]:
    branch = cb.current_branch(root)
    topic = cb.topic_of(branch)
    lines = [f"Research mode on {branch} (.claude/rules/research/). Read PROGRESS.md (engine, read-only here), "
             "then research/RESEARCH.md" + (f", then research/topics/{topic}/README.md." if topic else ".")]
    if branch == "main":
        lines.append("WARNING: research files are on main. Research is never merged into main: undo the merge.")
    main = cb.resolve(root, "origin/main", "main")
    if main:
        base = cb.git(root, "merge-base", "HEAD", main).strip()
        gated = cb.gate_record(root, base)
        lines.append(f"Engine: main {base[:9]} taken; pypic/ " + (f"gated ({gated})." if gated else "NOT gated."))
        ahead = cb.git(root, "rev-list", "--count", f"{base}..{main}").strip()
        engine = cb.git(root, "rev-list", "--count", f"{base}..{main}", "--", "pypic").strip()
        if ahead != "0":
            lines.append(f"{main} has {ahead} commits not taken ({engine} touch pypic/; refs as of the last fetch): "
                         "skill research-sync-main.")
    problems = cb.run_checks(root, None, None, branch, worktree=True)
    lines.append("Branch rules: OK." if not problems else
                 f"Branch rules: {len(problems)} violation(s), e.g. {problems[0]} "
                 "(python research/tools/check_branch.py --worktree).")
    return lines


def main() -> int:
    root = Path(os.environ.get("CLAUDE_PROJECT_DIR") or Path(__file__).resolve().parents[3])
    if not cb.research_mode(root):
        return 0
    try:
        print("\n".join(report(root)))
    except Exception as exc:  # noqa: BLE001 - a hook must not stop the session
        print(f"research session_start: {exc}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
