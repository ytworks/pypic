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
