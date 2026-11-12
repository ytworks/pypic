#!/usr/bin/env python3
"""Check a research branch against .claude/rules/research/branching.md.

    python research/tools/check_branch.py [--main origin/main] [--trunk origin/research/trunk]
                                          [--branch research/<topic>] [--worktree]

Run by CI (.github/workflows/research.yml), in summary form by the SessionStart hook
(research/tools/hooks/session_start.py), and by hand. Exits 1 on a violation.

1. Main-owned files: every path that differs between HEAD and merge-base(HEAD, main) is a
   research path. This is what keeps every `git merge origin/main` free of conflicts.
2. Topic scope: on research/<topic>, every path that differs from merge-base(HEAD, trunk) is in
   research/topics/<topic>/ or research/literature/.
3. research/RESEARCH.md: its three sections, at most MAX_LINES lines and MAX_ENTRIES entries.
4. Every topic README has the headings of research/templates/topic/README.md.
5. Run records (research/topics/*/runs/): strict JSON, index and files one to one, every recorded
   commit an ancestor of HEAD (so it stays reachable); files under research/ at most MAX_BYTES.

Pure Python + git; no CuPy, no package install.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

RESEARCH_PREFIXES = ("research/", ".claude/rules/research/", ".claude/skills/research-")
# The only file research branches share with main: the research hooks are added to it.
RESEARCH_FILES = (".github/workflows/research.yml", ".claude/settings.json")
TRUNK = "research/trunk"
TOPIC_SLUG = re.compile(r"[a-z0-9][a-z0-9-]*")
SECTIONS = ["## Topics", "## Engine baseline", "## Change Log"]
MAX_LINES = 150
MAX_ENTRIES = 10
MAX_BYTES = 1_000_000


def git(root: Path, *args: str, check: bool = True) -> str:
    r = subprocess.run(["git", *args], cwd=root, capture_output=True, text=True, timeout=60)
    if check and r.returncode:
        raise RuntimeError(f"git {' '.join(args)}: {r.stderr.strip()}")
    return r.stdout


def is_research_path(path: str) -> bool:
    return path.startswith(RESEARCH_PREFIXES) or path in RESEARCH_FILES


def research_mode(root: Path) -> bool:
    return (root / "research" / "RESEARCH.md").is_file()


def current_branch(root: Path) -> str:
    return git(root, "rev-parse", "--abbrev-ref", "HEAD", check=False).strip()


def topic_of(branch: str) -> str | None:
    """The topic of research/<topic>; None for trunk and for any other branch."""
    if not branch.startswith("research/") or branch == TRUNK:
        return None
    slug = branch.removeprefix("research/")
    return slug if TOPIC_SLUG.fullmatch(slug) else None


def resolve(root: Path, *candidates: str) -> str | None:
    for ref in candidates:
        if ref and git(root, "rev-parse", "--verify", "--quiet", f"{ref}^{{commit}}", check=False).strip():
            return ref
    return None


def changed_since(root: Path, base: str, worktree: bool) -> list[str]:
    """Paths that differ from ``base``: committed (HEAD), or with ``worktree`` also the uncommitted
    and the untracked, not ignored ones."""
    if worktree:
        paths = git(root, "diff", "--name-only", base).splitlines()
        paths += git(root, "ls-files", "--others", "--exclude-standard").splitlines()
    else:
        paths = git(root, "diff", "--name-only", base, "HEAD").splitlines()
    return sorted(set(p for p in paths if p))


def check_main_owned(root: Path, main: str, worktree: bool) -> list[str]:
    base = git(root, "merge-base", "HEAD", main).strip()
    bad = [p for p in changed_since(root, base, worktree) if not is_research_path(p)]
    return [f"{p}: main-owned file changed on a research branch (it must stay as main has it; "
            f"change it on a branch from main and take main with the skill research-sync-main)" for p in bad]


def check_topic_scope(root: Path, topic: str, trunk: str, worktree: bool) -> list[str]:
    base = git(root, "merge-base", "HEAD", trunk).strip()
    own = (f"research/topics/{topic}/", "research/literature/")
    bad = [p for p in changed_since(root, base, worktree) if not p.startswith(own)]
    return [f"{p}: outside research/topics/{topic}/ and research/literature/ (change it on {TRUNK}, "
            f"then merge {TRUNK} into research/{topic})" for p in bad]


def check_research_md(root: Path) -> list[str]:
    lines = (root / "research" / "RESEARCH.md").read_text(encoding="utf-8").splitlines()
    sections = [ln.rstrip() for ln in lines if ln.startswith("## ")]
    entries, in_log = 0, False
    for ln in lines:
        if ln.startswith("## "):
            in_log = ln.rstrip() == "## Change Log"
        elif in_log and ln.startswith("### "):
            entries += 1
    problems = []
    if sections != SECTIONS:
        problems.append(f"research/RESEARCH.md: sections {sections} (required: {SECTIONS})")
    if len(lines) > MAX_LINES:
        problems.append(f"research/RESEARCH.md: {len(lines)} lines (budget {MAX_LINES}); move the oldest "
                        "Change Log entries verbatim to research/archive/change-log.md")
    if entries > MAX_ENTRIES:
        problems.append(f"research/RESEARCH.md: {entries} Change Log entries (budget {MAX_ENTRIES})")
    return problems


def required_headings(root: Path) -> list[str]:
    template = root / "research" / "templates" / "topic" / "README.md"
    return [ln.rstrip() for ln in template.read_text(encoding="utf-8").splitlines() if ln.startswith("## ")]


def topics(root: Path) -> list[Path]:
    base = root / "research" / "topics"
    return sorted(p for p in base.iterdir() if p.is_dir() and TOPIC_SLUG.fullmatch(p.name)) if base.is_dir() else []


def check_topic_readmes(root: Path) -> list[str]:
    need = required_headings(root)
    problems = []
    for t in topics(root):
        readme = t / "README.md"
        if not readme.is_file():
            problems.append(f"research/topics/{t.name}/README.md: missing")
            continue
        have = {ln.rstrip() for ln in readme.read_text(encoding="utf-8").splitlines()}
        missing = [h for h in need if h not in have]
        if missing:
            problems.append(f"research/topics/{t.name}/README.md: missing headings {missing}")
    return problems


def loads_strict(text: str):
    def reject(token):
        raise ValueError(f"non-standard JSON constant {token}")
    return json.loads(text, parse_constant=reject)


def check_runs(root: Path) -> list[str]:
    problems = []
    for t in topics(root):
        runs = t / "runs"
        if not runs.is_dir():
            continue
        rel = f"research/topics/{t.name}/runs"
        files = {p.name for p in runs.glob("*.json") if p.name != "index.json"}
        index = runs / "index.json"
        if not index.is_file():
            if files:
                problems.append(f"{rel}: records without index.json")
            continue
        try:
            entries = loads_strict(index.read_text(encoding="utf-8")).get("entries", [])
        except ValueError as exc:
            problems.append(f"{rel}/index.json: not strict JSON ({exc})")
            continue
        listed = set()
        for e in entries:
            name = e.get("file")
            if not name or not e.get("commit"):
                problems.append(f"{rel}/index.json: entry without file or commit: {e}")
                continue
            listed.add(name)
            if name not in files:
                problems.append(f"{rel}/{name}: listed in index.json, not found")
                continue
            try:
                loads_strict((runs / name).read_text(encoding="utf-8"))
            except ValueError as exc:
                problems.append(f"{rel}/{name}: not strict JSON ({exc})")
            for key in ("commit", "analysis_commit"):
                c = str(e.get(key) or "").removesuffix("-dirty")
                if c and subprocess.run(["git", "merge-base", "--is-ancestor", c, "HEAD"], cwd=root,
                                        capture_output=True).returncode:
                    problems.append(f"{rel}/{name}: {key} {c} is not an ancestor of HEAD (rebased or "
                                    "unknown commit: the run can no longer be reproduced from history)")
        for name in sorted(files - listed):
            problems.append(f"{rel}/{name}: not listed in index.json")
    return problems


def check_sizes(root: Path) -> list[str]:
    problems = []
    for ln in git(root, "ls-tree", "-r", "-l", "HEAD", "--", "research").splitlines():
        meta, path = ln.split("\t", 1)
        size = meta.split()[3]
        if size.isdigit() and int(size) > MAX_BYTES:
            problems.append(f"{path}: {int(size)} bytes (limit {MAX_BYTES}); big outputs stay out of git")
    return problems


def gate_record(root: Path, commit: str) -> str | None:
    """The newest all-passed GPU gate record (bench/results/index.json) whose commit has the
    pypic/ tree of ``commit``; None when there is none."""
    index = root / "bench" / "results" / "index.json"
    if not index.is_file():
        return None
    gates = [e for e in json.loads(index.read_text())["entries"]
             if e.get("type", "gates") == "gates" and e.get("all_passed") and not e.get("dirty")]
    if not gates:
        return None
    out = git(root, "rev-parse", f"{commit}:pypic", *[f"{e['commit']}:pypic" for e in gates], check=False).split()
    if len(out) != len(gates) + 1:
        return None
    target = out[0]
    return next((e["file"] for e, tree in zip(gates, out[1:]) if tree == target), None)


def run_checks(root: Path, main: str | None, trunk: str | None, branch: str, worktree: bool) -> list[str]:
    problems = []
    main_ref = resolve(root, main, "origin/main", "main")
    if main_ref is None:
        problems.append("main not found (fetch it: git fetch origin main)")
    else:
        problems += check_main_owned(root, main_ref, worktree)
    topic = topic_of(branch)
    if topic:
        trunk_ref = resolve(root, trunk, f"origin/{TRUNK}", TRUNK)
        if trunk_ref is None:
            problems.append(f"{TRUNK} not found (fetch it: git fetch origin {TRUNK})")
        else:
            problems += check_topic_scope(root, topic, trunk_ref, worktree)
    problems += check_research_md(root)
    problems += check_topic_readmes(root)
    problems += check_runs(root)
    problems += check_sizes(root)
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--main", help="main ref (default: origin/main, else main)")
    ap.add_argument("--trunk", help=f"trunk ref (default: origin/{TRUNK}, else {TRUNK})")
    ap.add_argument("--branch", help="branch name (default: the checked-out branch)")
    ap.add_argument("--worktree", action="store_true", help="include uncommitted and untracked changes")
    a = ap.parse_args()
    if not research_mode(ROOT):
        print("not a research branch (no research/RESEARCH.md)", file=sys.stderr)
        return 1
    branch = a.branch or current_branch(ROOT)
    problems = run_checks(ROOT, a.main, a.trunk, branch, a.worktree)
    for p in problems:
        print(p)
    if problems:
        print(f"{len(problems)} violation(s) of .claude/rules/research/branching.md on {branch}", file=sys.stderr)
        return 1
    print(f"{branch}: research branch rules OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
