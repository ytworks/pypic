"""research/tools/hooks/guard.py (PreToolUse): exit 2 blocks, 0 lets the call run. One case per rule."""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest
from conftest import GIT_ENV, ROOT, git, write

GUARD = ROOT / "research" / "tools" / "hooks" / "guard.py"


def guard(root, tool, tool_input, script=GUARD):
    env = {**GIT_ENV, "CLAUDE_PROJECT_DIR": str(root)}
    event = json.dumps({"tool_name": tool, "tool_input": tool_input})
    return subprocess.run([sys.executable, str(script)], input=event, env=env, capture_output=True,
                          text=True).returncode


def bash(command):
    return ("Bash", {"command": command})


@pytest.mark.parametrize("tool, tool_input, code", [
    ("Edit", {"file_path": "pypic/a.py"}, 2),
    ("Write", {"file_path": "PROGRESS.md"}, 2),
    ("Edit", {"file_path": "research/topics/foo/README.md"}, 0),
    ("Write", {"file_path": ".claude/rules/research/stance.md"}, 0),
    ("Write", {"file_path": os.devnull}, 0),  # outside the repository
    (*bash("git push origin HEAD:main"), 2),
    (*bash("git push -u origin research/trunk"), 0),
    (*bash("git push --force origin research/trunk"), 2),
    (*bash("git push origin :research/trunk"), 2),
    (*bash("git rebase origin/main"), 2),
    (*bash('git checkout main && git merge --no-ff -m "x; y" research/trunk'), 2),
    (*bash("git merge --no-ff origin/main"), 0),
    ("mcp__github__create_pull_request", {"head": "research/trunk", "base": "main"}, 2),
    ("mcp__github__create_pull_request", {"head": "research/foo", "base": "research/trunk"}, 0),
])
def test_guard(repo, tool, tool_input, code):
    assert guard(repo, tool, tool_input) == code


def test_bare_push_with_upstream_on_main(repo):
    git(repo, "branch", "-u", "origin/main")
    assert guard(repo, *bash("git push")) == 2


@pytest.mark.parametrize("path, code", [("pypic/a.py", 2), ("research/topics/foo/notebook.md", 0)])
def test_commit_of_staged_files(repo, path, code):
    write(repo, path, "staged\n")
    git(repo, "add", path)
    assert guard(repo, *bash("git commit -m 'c'")) == code


def test_commit_during_merge_of_main_passes(repo):
    git(repo, "switch", "-q", "main")
    write(repo, "pypic/a.py", "engine change\n")
    git(repo, "commit", "-q", "-am", "engine")
    git(repo, "switch", "-q", "research/trunk")
    git(repo, "merge", "-q", "--no-ff", "--no-commit", "main")
    assert guard(repo, *bash("git commit -m 'take main'")) == 0


def test_inactive_outside_research_mode(repo):
    git(repo, "switch", "-q", "main")  # no research/RESEARCH.md in this tree
    assert guard(repo, *bash("git push origin main")) == 0


def test_settings_hook_without_the_script_lets_calls_run(tmp_path):
    """On a branch without research/ (main) the configured command must exit 0, not python's 2."""
    settings = json.loads((ROOT / ".claude" / "settings.json").read_text())
    command = settings["hooks"]["PreToolUse"][0]["hooks"][0]["command"]
    r = subprocess.run(["bash", "-c", command], input="{}", env={**os.environ, "CLAUDE_PROJECT_DIR": str(tmp_path)},
                       capture_output=True, text=True)
    assert r.returncode == 0
