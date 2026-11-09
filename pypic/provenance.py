"""Provenance of the code that produces a result (run.json, GPU gate and example records).

A result is reproducible from its commit only if the code in the working tree is that commit:
a commit is reported "-dirty" when a tracked file under the code paths differs from it, or when
an untracked, not ignored file lies there (an untracked helper module imported by an example,
say). Documentation and recorded evidence are not code: an uncommitted bench/results/index.json
does not make the next run dirty. :func:`code_state` also gives a digest of the differences,
so a caller can check that the code did not change while it ran. Pure Python (git), no CuPy.
"""
from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

CODE_PATHS = ("pypic", "examples", "scripts", "tests", "pyproject.toml")


def _git(root, *args: str, text: bool = True):
    r = subprocess.run(["git", *args], cwd=root, capture_output=True, text=text, timeout=30)
    if r.returncode:
        raise RuntimeError(f"git {' '.join(args)} failed")
    return r.stdout


def repo_root(start) -> Path | None:
    """Top of the git working tree that contains ``start``, or None (e.g. an installed package)."""
    try:
        return Path(_git(start, "rev-parse", "--show-toplevel").strip())
    except Exception:  # noqa: BLE001 - provenance is best effort
        return None


def code_state(root, paths=CODE_PATHS) -> tuple[str, str]:
    """(HEAD, digest): the digest is "" when the code is the commit, else a hash of the tracked
    differences and of the untracked, not ignored files (names and contents) under ``paths``."""
    root = Path(root)
    head = _git(root, "rev-parse", "HEAD").strip()
    diff = _git(root, "diff", "HEAD", "--binary", "--", *paths, text=False)
    untracked = sorted(p for p in _git(root, "ls-files", "--others", "--exclude-standard", "-z", "--", *paths)
                       .split("\0") if p)
    if not diff and not untracked:
        return head, ""
    h = hashlib.sha256(diff)
    for p in untracked:
        h.update(b"\0untracked\0" + p.encode() + b"\0")
        f = root / p
        h.update(f.read_bytes() if f.is_file() else b"")
    return head, h.hexdigest()


def commit_label(root, paths=CODE_PATHS) -> str:
    """HEAD, with "-dirty" appended when the code differs from it (:func:`code_state`)."""
    head, changes = code_state(root, paths)
    return head + ("-dirty" if changes else "")
