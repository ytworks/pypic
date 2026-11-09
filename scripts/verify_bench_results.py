#!/usr/bin/env python3
"""Validate bench/results/index.json against the result files (run on CPU CI).

Every index entry must point to an existing result file whose commit, verdict and gate
records agree with the entry; every result file must be indexed. This keeps the committed
GPU evidence (scripts/gpu_gates.py) consistent and reviewable.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

RES = Path(__file__).resolve().parent.parent / "bench" / "results"
ENTRY_KEYS = {"file", "date", "commit", "dirty", "gpus", "gpu_name", "all_passed"}
GATE_KEYS = {"name", "passed", "metrics", "criteria"}


def main() -> int:
    index_path = RES / "index.json"
    if not index_path.exists():
        print("no bench/results/index.json (no GPU evidence recorded yet)")
        return 0
    errors = []
    index = json.loads(index_path.read_text())
    seen = set()
    for e in index.get("entries", []):
        missing = ENTRY_KEYS - set(e)
        if missing:
            errors.append(f"index entry {e.get('file')}: missing keys {sorted(missing)}")
            continue
        seen.add(e["file"])
        path = RES / e["file"]
        if not path.exists():
            errors.append(f"index entry {e['file']}: file not found")
            continue
        rec = json.loads(path.read_text())
        for k in ("commit", "dirty", "gpus", "all_passed"):
            if rec.get(k) != e[k]:
                errors.append(f"{e['file']}: {k} differs between index ({e[k]}) and file ({rec.get(k)})")
        gates = rec.get("gates", [])
        if not gates:
            errors.append(f"{e['file']}: no gates")
        for g in gates:
            if GATE_KEYS - set(g):
                errors.append(f"{e['file']}: gate {g.get('name')} missing keys {sorted(GATE_KEYS - set(g))}")
        if rec.get("all_passed") != all(g.get("passed") for g in gates):
            errors.append(f"{e['file']}: all_passed does not match the gate verdicts")
    for f in RES.glob("*-gpu-gates.json"):
        if f.name not in seen:
            errors.append(f"{f.name}: not listed in index.json")
    for msg in errors:
        print(msg, file=sys.stderr)
    if not errors:
        print(f"bench/results: {len(seen)} indexed result file(s) consistent")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
