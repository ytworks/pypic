"""The wiki articles of the examples (scripts/build_wiki.py render: docs/wiki templates + data.json).

Every example page, in English and Japanese, must render with no placeholder left, show the experimental
conditions, at least one animation and one comparison chart, and reference only images the builder made.
"""
from __future__ import annotations

import json
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
