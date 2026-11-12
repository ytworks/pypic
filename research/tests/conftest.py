"""A throw-away repository shaped like pypic with a research trunk: main with main-owned files,
research/trunk with the research memory and the topic template of this checkout."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research" / "tools"))

GIT_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1",
           "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.invalid",
           "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@example.invalid"}


def git(cwd: Path, *args: str) -> str:
    return subprocess.run(["git", *args], cwd=cwd, env=GIT_ENV, check=True, capture_output=True, text=True).stdout


def write(root: Path, rel: str, text: str = "x\n") -> None:
    p = root / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)


def commit(root: Path, message: str = "c") -> str:
    git(root, "add", "-A")
    git(root, "commit", "-q", "-m", message)
    return git(root, "rev-parse", "HEAD").strip()


def add_topic(root: Path, slug: str) -> None:
    shutil.copytree(ROOT / "research" / "templates" / "topic", root / "research" / "topics" / slug)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """research/trunk checked out; origin (bare) has main, and trunk tracks nothing."""
    remote, root = tmp_path / "remote.git", tmp_path / "work"
    git(tmp_path, "init", "-q", "--bare", "-b", "main", str(remote))
    git(tmp_path, "init", "-q", "-b", "main", str(root))
    for rel in ("pypic/a.py", "PROGRESS.md", "CLAUDE.md", ".claude/settings.json"):
        write(root, rel)
    commit(root, "main")
    git(root, "remote", "add", "origin", str(remote))
    git(root, "push", "-q", "origin", "main")
    git(root, "switch", "-q", "-c", "research/trunk")
    write(root, "research/RESEARCH.md", (ROOT / "research" / "RESEARCH.md").read_text())
    shutil.copytree(ROOT / "research" / "templates", root / "research" / "templates")
    commit(root, "research trunk")
    return root
