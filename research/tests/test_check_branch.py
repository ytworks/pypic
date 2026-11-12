"""research/tools/check_branch.py: one case per rule it enforces."""
from __future__ import annotations

import json

import check_branch as cb
from conftest import add_topic, commit, git, write


def problems(root, branch="research/trunk"):
    return cb.run_checks(root, None, None, branch, worktree=False)


def test_clean_trunk_passes(repo):
    assert problems(repo) == []


def test_main_owned_file_changed(repo):
    write(repo, "pypic/a.py", "changed\n")
    commit(repo)
    assert [p for p in problems(repo) if p.startswith("pypic/a.py: main-owned")]


def test_merge_of_main_is_not_a_change(repo):
    git(repo, "switch", "-q", "main")
    write(repo, "pypic/a.py", "engine change on main\n")
    commit(repo)
    git(repo, "push", "-q", "origin", "main")
    git(repo, "switch", "-q", "research/trunk")
    git(repo, "merge", "-q", "--no-ff", "-m", "take main", "main")
    assert problems(repo) == []


def test_topic_changes_another_topic(repo):
    git(repo, "switch", "-q", "-c", "research/foo")
    add_topic(repo, "foo")
    add_topic(repo, "bar")
    commit(repo)
    found = problems(repo, "research/foo")
    assert [p for p in found if p.startswith("research/topics/bar/") and "outside" in p]
    assert not [p for p in found if p.startswith("research/topics/foo/")]


def test_topic_readme_without_required_heading(repo):
    add_topic(repo, "foo")
    readme = repo / "research/topics/foo/README.md"
    readme.write_text(readme.read_text().replace("## Claims and evidence", "## Results"))
    commit(repo)
    assert [p for p in problems(repo) if "missing headings ['## Claims and evidence']" in p]


def test_recorded_commit_must_stay_in_history(repo):
    add_topic(repo, "foo")
    commit(repo)
    git(repo, "switch", "-q", "-c", "scratch")
    write(repo, "research/topics/foo/notebook.md", "a commit that is rebased away\n")
    lost = commit(repo)
    git(repo, "switch", "-q", "research/trunk")
    runs = "research/topics/foo/runs"
    write(repo, f"{runs}/r1.json", json.dumps({"commit": lost}))
    write(repo, f"{runs}/index.json", json.dumps({"entries": [{"file": "r1.json", "commit": lost}]}))
    commit(repo)
    assert [p for p in problems(repo) if f"commit {lost} is not an ancestor of HEAD" in p]
