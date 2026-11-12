#!/usr/bin/env python3
"""PreToolUse hook of the research branches (.claude/settings.json): a tripwire for the rules of
.claude/rules/research/branching.md. Exit 2 blocks the tool call and shows the reason to Claude;
exit 0 lets it run. CI (research/tools/check_branch.py) is the backstop: the command checks here
are string matches and see one command at a time.

Active only in research mode (research/RESEARCH.md in the working tree). Blocks:
- Edit / Write / NotebookEdit of a main-owned file (anything but the research paths);
- git push to main (by refspec, or by an upstream on main), force pushes, --all / --mirror, deletes;
- git rebase; git commit of staged main-owned files outside a merge; any commit or push on main;
- git merge / cherry-pick in a command that also switches to main;
- any pull request: research branches use none (merges are local `git merge --no-ff` + push).
"""
from __future__ import annotations

import json
import os
import re
import shlex
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from check_branch import current_branch, git, is_research_path, research_mode  # noqa: E402

FILE_TOOLS = {"Edit": "file_path", "Write": "file_path", "MultiEdit": "file_path", "NotebookEdit": "notebook_path"}
PR_TOOLS = {"mcp__github__create_pull_request", "mcp__github__update_pull_request"}
SEPARATORS = re.compile(r"&&|\|\||[;|\n&]")
FORCE = {"-f", "--force", "--mirror", "--all", "-d", "--delete"}
RULES = ".claude/rules/research/branching.md"
NO_PR = ("pull request: research branches use none, since nothing on them is merged by a pull request; "
         "integrate with a local `git merge --no-ff` and push")


def block(reason: str) -> int:
    print(f"Blocked by the research-branch guard: {reason} (rules: {RULES})", file=sys.stderr)
    return 2


def segments(command: str) -> list[list[str]]:
    """The simple commands of a shell command line, as word lists (quotes respected)."""
    try:
        lex = shlex.shlex(command.replace("\n", " ; "), posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        out, cur = [], []
        for tok in lex:
            if tok and set(tok) <= set("&|;()"):
                out.append(cur)
                cur = []
            else:
                cur.append(tok)
        return [*out, cur]
    except ValueError:  # unbalanced quotes: split on the operators, then on blanks
        return [s.split() for s in SEPARATORS.split(command)]


def git_calls(command: str) -> list[list[str]]:
    """The argument lists of the git invocations in a shell command, global options dropped."""
    calls = []
    for words in segments(command):
        while words and "=" in words[0] and not words[0].startswith("-"):  # VAR=value prefixes
            words = words[1:]
        if not words or Path(words[0]).name != "git":
            continue
        args, i = words[1:], 0
        while i < len(args) and args[i].startswith("-"):
            i += 2 if args[i] in ("-C", "-c", "--git-dir", "--work-tree") else 1
        if i < len(args):
            calls.append(args[i:])
    return calls


def short(ref: str) -> str:
    return ref.removeprefix("refs/heads/")


def check_push(root: Path, args: list[str], branch: str) -> str | None:
    flags = [a for a in args[1:] if a.startswith("-")]
    if any(f in FORCE or f.startswith("--force") for f in flags):
        return "force push, --all, --mirror or a branch delete: research history must stay as it is"
    positional = [a for a in args[1:] if not a.startswith("-")]
    refspecs = positional[1:]
    if any(r.startswith(("+", ":")) for r in refspecs):
        return "forced (+) or deleting (:) refspec: research history must stay as it is"
    if branch == "main":
        return "push from main while research files are in the tree: research is never merged into main"
    for r in refspecs:
        dst = short(r.split(":", 1)[1] if ":" in r else r)
        if dst == "main":
            return f"push to main ({r}): research branches are never pushed to main"
    if not refspecs:
        upstream = git(root, "rev-parse", "--abbrev-ref", "@{upstream}", check=False).strip()
        if upstream.split("/", 1)[-1] == "main":
            return (f"this branch tracks {upstream}: a bare `git push` would push it to main; "
                    f"push with `git push -u origin {branch}`")
    return None


def check_commit(root: Path, args: list[str], branch: str) -> str | None:
    if branch == "main":
        return "commit on main with research files in the tree: research is never merged into main"
    gitdir = Path(git(root, "rev-parse", "--absolute-git-dir", check=False).strip() or root / ".git")
    if any((gitdir / m).exists() for m in ("MERGE_HEAD", "CHERRY_PICK_HEAD", "REVERT_HEAD")):
        return None  # a merge of main brings main's files; check_branch.py checks the result
    staged = git(root, "diff", "--cached", "--name-only").splitlines()
    if "-a" in args or "--all" in args or any(re.fullmatch(r"-[a-zA-Z]*a[a-zA-Z]*", a) for a in args[1:]):
        staged += git(root, "diff", "--name-only").splitlines()
    bad = sorted({p for p in staged if p and not is_research_path(p)})
    if bad:
        return (f"commit of main-owned files {bad[:5]}: they stay as main has them; change them on a "
                "branch from main, then take main with the skill research-sync-main")
    return None


def check_bash(root: Path, command: str) -> str | None:
    if any(w[:3] == ["gh", "pr", "create"] for w in segments(command)):
        return NO_PR
    calls = git_calls(command)
    if not calls:
        return None
    branch = current_branch(root)
    to_main = any(c[0] in ("checkout", "switch") and "main" in c[1:] for c in calls)
    for args in calls:
        sub = args[0]
        if sub == "push":
            reason = check_push(root, args, branch)
        elif sub == "rebase" and not {"--abort", "--quit"} & set(args):
            reason = "rebase rewrites research history (run records name its commits); merge with --no-ff instead"
        elif sub in ("merge", "cherry-pick") and (to_main or branch == "main"):
            reason = f"git {sub} on main: research is never merged into main"
        elif sub == "commit":
            reason = check_commit(root, args, branch)
        else:
            reason = None
        if reason:
            return reason
    return None


def check_file(root: Path, path: str) -> str | None:
    p = Path(path)
    p = (p if p.is_absolute() else root / p).resolve()
    try:
        rel = p.relative_to(root.resolve()).as_posix()
    except ValueError:
        return None  # outside the repository (scratchpad, plans)
    if rel.startswith(".git/") or is_research_path(rel):
        return None
    return (f"{rel} is main-owned: it stays as main has it on research branches. Engine and harness "
            "changes go on a branch from main (PR, gates), then into research with the skill research-sync-main")


def main() -> int:
    try:
        event = json.load(sys.stdin)
    except ValueError:
        return 0
    root = Path(os.environ.get("CLAUDE_PROJECT_DIR") or Path(__file__).resolve().parents[3])
    if not research_mode(root):
        return 0
    tool, tool_input = event.get("tool_name", ""), event.get("tool_input") or {}
    if tool in FILE_TOOLS:
        reason = check_file(root, str(tool_input.get(FILE_TOOLS[tool], "")))
    elif tool == "Bash":
        reason = check_bash(root, str(tool_input.get("command", "")))
    elif tool in PR_TOOLS:
        reason = NO_PR
    else:
        reason = None
    return block(reason) if reason else 0


if __name__ == "__main__":
    sys.exit(main())
