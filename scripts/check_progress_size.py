#!/usr/bin/env python3
"""Check PROGRESS.md against the size budget in .claude/rules/progress-log.md.

Run by the SessionStart hook in .claude/settings.json, so every agent session starts with a
warning when PROGRESS.md needs archiving, and by CI. Exits non-zero when the budget is
exceeded or the required section layout is broken.
"""
from __future__ import annotations

import sys
from pathlib import Path

MAX_TOTAL_LINES = 200
MAX_CHANGELOG_ENTRIES = 10
SECTIONS = ["## Product", "## Current Status", "## Change Log"]


def main() -> int:
    progress = Path(__file__).resolve().parent.parent / "PROGRESS.md"
    if not progress.is_file():
        print(f"check_progress_size: {progress} not found", file=sys.stderr)
        return 1
    lines = progress.read_text(encoding="utf-8").splitlines()

    sections = [ln.rstrip() for ln in lines if ln.startswith("## ")]
    in_changelog = False
    entries = 0
    for line in lines:
        if line.startswith("## "):
            in_changelog = line.rstrip() == "## Change Log"
        elif in_changelog and line.startswith("### "):
            entries += 1

    problems = []
    if sections != SECTIONS:
        problems.append(f"sections {sections} (required: {SECTIONS})")
    if len(lines) > MAX_TOTAL_LINES:
        problems.append(f"{len(lines)} lines (budget: {MAX_TOTAL_LINES})")
    if entries > MAX_CHANGELOG_ENTRIES:
        problems.append(f"{entries} Change Log entries (budget: {MAX_CHANGELOG_ENTRIES})")
    if problems:
        print("PROGRESS.md violates its rules: " + "; ".join(problems)
              + ". See .claude/rules/progress-log.md; to archive, follow docs/progress/archiving.md"
              " before other work updates PROGRESS.md.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
