"""The wiki articles of the examples (scripts/build_wiki.py render: docs/wiki templates + data.json).

Every example page, in English and Japanese, must render with no placeholder left, show the experimental
conditions, at least one animation and one comparison chart, and reference only images the builder made.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
IMAGE = re.compile(r"images/[\w\-./]+\.(?:png|gif)")


def test_every_example_page_renders_in_both_languages(tmp_path):
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "build_wiki.py"), "render", "--wiki", str(tmp_path)],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    data = json.loads((ROOT / "docs" / "wiki" / "data.json").read_text())
    made = {f"images/{ex}/{name}" for ex, e in data["examples"].items() for name in e["images"]}
    sidebar = (tmp_path / "_Sidebar.md").read_text()
    for ex, e in data["examples"].items():
        pages = {lang: (tmp_path / f"{e['page']}{suffix}.md").read_text()
                 for lang, suffix in (("en", ""), ("ja", "-ja"))}
        refs = {lang: set(IMAGE.findall(text)) for lang, text in pages.items()}
        assert refs["en"] == refs["ja"], ex
        assert refs["en"] <= made, (ex, refs["en"] - made)
        assert any(x.endswith(".gif") for x in refs["en"]), f"{ex}: no animation"
        assert any(x.endswith(".png") for x in refs["en"]), f"{ex}: no comparison chart"
        assert "**grid (cells)**" in pages["en"] and "**格子数**" in pages["ja"], f"{ex}: no conditions table"
        assert f"({e['page']})" in sidebar and f"({e['page']}-ja)" in sidebar


def _builder():
    import importlib.util

    spec = importlib.util.spec_from_file_location("build_wiki", ROOT / "scripts" / "build_wiki.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_a_recorded_case_must_be_the_run_of_its_evidence(tmp_path, monkeypatch):
    B = _builder()
    monkeypatch.setattr(B, "ROOT", tmp_path)
    (tmp_path / "bench").mkdir()
    for name, ev in (("new.json", {"example": "reconnection", "config": "large", "run_id": "e879cfd6",
                                   "commit": "d9637cd9f", "run": {"wall_time_s": 686.2}}),
                     ("old.json", {"example": "reconnection", "config": "large", "commit": "200af1306",
                                   "run": {"wall_time_s": 635.9}})):
        (tmp_path / "bench" / name).write_text(json.dumps(ev))

    def case(evidence, **run):
        d = tmp_path / "out" / evidence
        d.mkdir(parents=True, exist_ok=True)
        meta = {"config_file": "examples/reconnection/configs/large.yaml", "run": {"wall_time_s": run.pop("wall")}}
        (d / "run.json").write_text(json.dumps({**run, "meta": meta}))
        return {"path": f"out/{evidence}", "evidence": f"bench/{evidence}"}

    assert B.check_recorded("reconnection", "x", case("new.json", run_id="e879cfd6", wall=686.2)) is None
    assert B.check_recorded("reconnection", "x", case("new.json", run_id="different_realization", wall=686.2))
    assert B.check_recorded("reconnection", "x", case("old.json", git_commit="200af1306", wall=635.9)) is None
    assert B.check_recorded("reconnection", "x", case("old.json", run_id="abc", git_commit="200af1306", wall=635.9))
    assert B.check_recorded("reconnection", "x", {"path": "out/missing", "evidence": "bench/new.json"})


def test_each_page_names_the_analysis_commit_of_its_own_data(tmp_path):
    B = _builder()
    data = json.loads((ROOT / "docs" / "wiki" / "data.json").read_text())
    data["examples"]["reconnection"]["analysis_commit"] = "b" * 40
    B.render(data, tmp_path)
    for ex, e in data["examples"].items():
        for suffix in ("", "-ja"):
            assert f"`{e['analysis_commit'][:9]}`" in (tmp_path / f"{e['page']}{suffix}.md").read_text(), ex
    assert "`bbbbbbbbb`" not in (tmp_path / "Two-Stream-Instability.md").read_text()


def test_a_partial_collect_keeps_the_history_of_the_other_examples(tmp_path):
    B = _builder()
    data = json.loads((ROOT / "docs" / "wiki" / "data.json").read_text())
    legacy = {"analysis_commit": "a" * 40, "generated": "2026-01-01",
              "examples": {ex: {k: v for k, v in e.items() if k not in ("analysis_commit", "generated")}
                           for ex, e in data["examples"].items()}}  # the format with one commit for all
    new = B.merge_collected(legacy, {"particle_orbits": legacy["examples"]["particle_orbits"]},
                            ("b" * 40, ""), ("b" * 40, ""), "2026-02-02")
    assert new["examples"]["particle_orbits"]["analysis_commit"] == "b" * 40
    assert new["examples"]["reconnection"]["analysis_commit"] == "a" * 40
    assert new["examples"]["reconnection"]["generated"] == "2026-01-01"
    B.render(new, tmp_path)


def test_nothing_is_merged_when_the_code_changed_during_the_analysis():
    B = _builder()
    with pytest.raises(RuntimeError):
        B.merge_collected({"examples": {}}, {}, ("a" * 40, ""), ("b" * 40, ""), "2026-02-02")


def test_the_wiki_images_are_kept_when_the_collect_is_refused(tmp_path, monkeypatch):
    import argparse
    import subprocess as sp

    import pypic.provenance

    B = _builder()
    wiki = tmp_path / "wiki"
    (wiki / "images" / "reconnection").mkdir(parents=True)
    (wiki / "images" / "reconnection" / "power_law.png").write_bytes(b"old")
    monkeypatch.setattr(B, "ROOT", tmp_path)
    monkeypatch.setattr(B, "DATA", tmp_path / "data.json")
    monkeypatch.setattr(B, "check_recorded", lambda *a: None)
    monkeypatch.setattr(B, "_complete", lambda d: True)

    def worker(cmd, **kw):  # a collect worker: writes its charts and its JSON
        out = Path(cmd[cmd.index("--wiki") + 1]) / "images" / cmd[3]
        out.mkdir(parents=True, exist_ok=True)
        (out / "power_law.png").write_bytes(b"new")
        Path(cmd[cmd.index("--json") + 1]).write_text(json.dumps({"page": "Magnetic-Reconnection", "cases": {}}))
        return sp.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(B.subprocess, "run", worker)
    states = iter([("a" * 40, ""), ("b" * 40, "")])  # the code changes during the analysis
    monkeypatch.setattr(pypic.provenance, "code_state", lambda root: next(states))
    args = argparse.Namespace(wiki=str(wiki), example=["reconnection"], jobs=1)
    assert B.cmd_collect(args) == 1
    assert (wiki / "images" / "reconnection" / "power_law.png").read_bytes() == b"old"
    assert not (tmp_path / "data.json").exists()

    states = iter([("a" * 40, ""), ("a" * 40, "")])  # unchanged code: the new charts are published
    assert B.cmd_collect(args) == 0
    assert (wiki / "images" / "reconnection" / "power_law.png").read_bytes() == b"new"


def _fake_collect(B, tmp_path, monkeypatch, after_worker=None):
    """cmd_collect with fake workers that write a chart reading the --wiki argument's run label."""
    import subprocess as sp

    import pypic.provenance

    monkeypatch.setattr(B, "ROOT", tmp_path)
    monkeypatch.setattr(B, "DATA", tmp_path / "data.json")
    monkeypatch.setattr(B, "check_recorded", lambda *a: None)
    monkeypatch.setattr(B, "_complete", lambda d: True)
    monkeypatch.setattr(pypic.provenance, "code_state", lambda root: ("a" * 40, ""))
    label = {}

    def worker(cmd, **kw):
        out = Path(cmd[cmd.index("--wiki") + 1]) / "images" / cmd[3]
        out.mkdir(parents=True, exist_ok=True)
        (out / "power_law.png").write_bytes(label["current"])
        Path(cmd[cmd.index("--json") + 1]).write_text(json.dumps({"page": "Magnetic-Reconnection", "cases": {}}))
        if after_worker:
            after_worker()
        return sp.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(B.subprocess, "run", worker)
    return label


def test_concurrent_collects_keep_their_own_images(tmp_path, monkeypatch):
    import argparse

    B = _builder()
    wiki_a, wiki_b = tmp_path / "wiki_a", tmp_path / "wiki_b"
    args = {w: argparse.Namespace(wiki=str(w), example=["reconnection"], jobs=1) for w in (wiki_a, wiki_b)}
    nested = []

    def run_b_meanwhile():  # collect B runs to completion while A has its results but has not published
        if not nested:
            nested.append(1)
            label["current"] = b"B"
            assert B.cmd_collect(args[wiki_b]) == 0

    label = _fake_collect(B, tmp_path, monkeypatch, run_b_meanwhile)
    label["current"] = b"A"
    assert B.cmd_collect(args[wiki_a]) == 0
    assert (wiki_a / "images" / "reconnection" / "power_law.png").read_bytes() == b"A"
    assert (wiki_b / "images" / "reconnection" / "power_law.png").read_bytes() == b"B"


@pytest.mark.skipif(os.name != "posix" or os.geteuid() == 0, reason="POSIX permissions; root ignores them")
def test_a_failed_replacement_keeps_the_old_images(tmp_path, monkeypatch):
    import argparse

    B = _builder()
    wiki = tmp_path / "wiki"
    old = wiki / "images" / "reconnection"
    old.mkdir(parents=True)
    (old / "power_law.png").write_bytes(b"old")
    label = _fake_collect(B, tmp_path, monkeypatch)
    label["current"] = b"new"
    (wiki / "images").chmod(0o555)  # the example's directory stays writable, its parent does not
    try:
        assert B.cmd_collect(argparse.Namespace(wiki=str(wiki), example=["reconnection"], jobs=1)) == 1
    finally:
        (wiki / "images").chmod(0o755)
    assert (old / "power_law.png").read_bytes() == b"old"
    assert not (tmp_path / "data.json").exists()


def test_conditions_of_a_run_with_a_laser_a_bunch_and_a_moving_window(tmp_path):
    B = _builder()
    shape = [{"shape": "gauss", "axis": ax, "width": 1.4, "center": 12.0} for ax in "xyz"]
    beam = {"density": 0.05, "drift": [9.95, 0.0, 0.0], "theta": 0.0, "profile": shape}
    cfg = {"ppc": 8, "boundaries": [["open", "reflecting"], ["periodic", "periodic"], ["periodic", "periodic"]],
           "species": [{"name": "bunch", "charge": -1.0, "mass": 1e6, "beams": [beam]}],
           "lasers": [{"plane_x": 2.0, "a0": 0.1, "omega": 10.0, "waist": 4.0, "focus_x": 31.4}],
           "window": {"start": 20.7, "margin": 1.0}}
    run = {"solver": "explicit", "n_cells": [8, 4, 4], "L": [8.0, 4.0, 4.0], "dx": [1.0, 1.0, 1.0], "dt": 0.5,
           "n_steps": 10, "config": cfg, "meta": {}, "window": {"start": 20.7, "window_shift_end": 3}}
    (tmp_path / "run.json").write_text(json.dumps(run))
    c = B.conditions(tmp_path)
    assert c["species"] == ["bunch (q=-1, m=1e+06): n=0.05, u=(9.95, 0, 0), θ=0, "
                            "gauss along x × gauss along y × gauss along z"]
    assert c["fields"] == ["laser antenna at x = 2: a0 = 0.1, ω = 10; waist 4; focus at x = 31.4"]
    assert c["boundaries"][-1] == "moving window: speed c from t = 20.7" and c["window"]["window_shift_end"] == 3


def test_conditions_name_shaped_perturbations_potentials_and_an_inflow_fluctuation(tmp_path):
    """PIC-D24 / PIC-D26 runs: the table names the shapes instead of printing their dicts (a table of
    an envelope h(x) holds a hundred numbers)."""
    B = _builder()
    h = {"shape": "table", "x": [0.0, 1.0], "values": [0.0, 1.0]}
    cos = {"shape": "cos", "mode": 1, "length": 32.0}
    beam = {"density": 1.0, "drift": [0.0, 0.0, 0.0], "theta": 0.0,
            "perturbation": [{"shape": {"y": cos, "x": h}, "velocity": 0.005, "velocity_component": "x"},
                             {"shape": {"y": cos}, "displacement": 0.0005, "displacement_component": "z"}]}
    cfg = {"ppc": 8, "boundaries": [["reflecting", "open"], ["periodic", "periodic"], ["periodic", "periodic"]],
           "species": [{"name": "e", "charge": -1.0, "mass": 1.0, "beams": [beam]}],
           "initial_fields": [{"kind": "potential", "component": "Az", "shape": {"y": cos, "x": h}, "scale": 0.1}],
           "injection": {"faces": {"x+": "flux"}, "fluctuation": {"face": "x+", "modes": [{}]}}}
    run = {"solver": "explicit", "n_cells": [8, 4, 4], "L": [8.0, 4.0, 4.0], "dx": [1.0, 1.0, 1.0], "dt": 0.5,
           "n_steps": 10, "config": cfg, "meta": {}}
    (tmp_path / "run.json").write_text(json.dumps(run))
    c = B.conditions(tmp_path)
    assert c["species"] == ["e (q=-1, m=1): n=1, θ=0, perturbation: u_x 0.005 cos(y) × table(x), ξ_z 0.0005 cos(y)"]
    assert c["fields"] == ["Az potential cos(y) × table(x)", "inflow fluctuation at x+: 1 mode(s)"]


def test_analysis_figures_are_written_into_the_collects_own_directory(tmp_path, monkeypatch):
    import sys
    import types

    B = _builder()
    run = tmp_path / "run"
    run.mkdir()
    seen = {}

    def analyze(run_dir, show_plot=True, fig_dir=None):  # a stand-in for examples/<ex>/analyze.py
        seen["fig_dir"] = fig_dir
        path = Path(fig_dir or run_dir) / "summary.png"
        path.write_bytes(b"png")
        return {"plot": str(path)}

    monkeypatch.setitem(sys.modules, "analyze", types.SimpleNamespace(analyze=analyze))
    monkeypatch.setitem(B.EXAMPLES, "toy", {"page": "Toy", "title": {}, "cases": {"a": {"path": "run"}}})
    monkeypatch.setattr(B, "ROOT", tmp_path)
    monkeypatch.setattr(B, "conditions", lambda d: {})
    monkeypatch.setattr(B, "provenance", lambda d: {})
    work = tmp_path / "work"
    out = B.collect_example("toy", work)
    assert Path(seen["fig_dir"]).is_relative_to(work)
    assert not (run / "summary.png").exists()  # the shared run directory is left alone
    assert (work / "images" / "toy" / "a-summary.png").read_bytes() == b"png" and "a-summary.png" in out["images"]


@pytest.mark.skipif(os.name != "posix", reason="POSIX file locks")
def test_publishing_holds_the_wiki_lock(tmp_path, monkeypatch):
    import argparse
    import fcntl

    B = _builder()
    wiki = tmp_path / "wiki"
    (wiki / ".git").mkdir(parents=True)
    label = _fake_collect(B, tmp_path, monkeypatch)
    label["current"] = b"new"
    held = []
    real = B.publish

    def publish(*a):  # another checkout's collect cannot take the same wiki's lock meanwhile
        with open(wiki / ".git" / "build_wiki.lock", "w") as other:
            try:
                fcntl.flock(other, fcntl.LOCK_EX | fcntl.LOCK_NB)
                held.append(False)
            except BlockingIOError:
                held.append(True)
        return real(*a)

    monkeypatch.setattr(B, "publish", publish)
    assert B.cmd_collect(argparse.Namespace(wiki=str(wiki), example=["reconnection"], jobs=1)) == 0
    assert held == [True]


def test_an_interrupted_publish_restores_the_replaced_images(tmp_path, monkeypatch):
    B = _builder()
    monkeypatch.setattr(B, "DATA", tmp_path / "data.json")
    src, images = tmp_path / "src", tmp_path / "wiki" / "images"
    for ex in ("a", "b"):
        (src / ex).mkdir(parents=True)
        (src / ex / "x.png").write_bytes(b"new")
        (images / ex).mkdir(parents=True)
        (images / ex / "x.png").write_bytes(b"old")
    real, calls = B.shutil.copytree, []

    def copytree(s_, d_, *a, **k):  # Ctrl+C while copying the second example
        calls.append(s_)
        if len(calls) == 2:
            raise KeyboardInterrupt
        return real(s_, d_, *a, **k)

    monkeypatch.setattr(B.shutil, "copytree", copytree)
    with pytest.raises(KeyboardInterrupt):
        B.publish(src, images, ["a", "b"], "{}")
    assert (images / "a" / "x.png").read_bytes() == b"old" and (images / "b" / "x.png").read_bytes() == b"old"
    assert sorted(p.name for p in images.iterdir()) == ["a", "b"] and not (tmp_path / "data.json").exists()


def test_collect_refuses_before_any_analysis_without_file_locks(tmp_path, monkeypatch):
    import argparse
    import sys

    B = _builder()
    label = _fake_collect(B, tmp_path, monkeypatch)
    label["current"] = b"new"
    monkeypatch.setitem(sys.modules, "fcntl", None)  # as on Windows

    def no_worker(*a, **k):
        raise AssertionError("an analysis was started")

    monkeypatch.setattr(B.subprocess, "run", no_worker)
    assert B.cmd_collect(argparse.Namespace(wiki=str(tmp_path / "wiki"), example=["reconnection"], jobs=1)) == 1
