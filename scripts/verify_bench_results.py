#!/usr/bin/env python3
"""Validate bench/results/index.json against the result files (run on CPU CI).

Two record types: "gates" (scripts/gpu_gates.py; the default for entries without a type; a
`precision: mixed` gate record, PIC-D11, says so in file, index entry and file name) and
"example" (scripts/record_example.py: reported results of example runs, e.g. production
configurations). Every index entry must point to an existing result file whose commit and
verdicts agree with the entry; every result file must be indexed; every file is strict JSON
(NaN / Infinity are rejected: unmeasured values are null). An example record with a run id
(one record per run) carries it in its file name, and no two records share a run id. This keeps
the committed GPU evidence consistent and reviewable.
"""
from __future__ import annotations

import sys
from pathlib import Path

from bench_records import loads_strict  # scripts/bench_records.py (this directory is on sys.path)

RES = Path(__file__).resolve().parent.parent / "bench" / "results"
ENTRY_KEYS = {"file", "date", "commit", "dirty", "gpus", "gpu_name", "all_passed"}
EXAMPLE_KEYS = {"file", "type", "date", "commit", "dirty", "gpus", "gpu_name", "example", "config"}
GATE_KEYS = {"name", "passed", "metrics", "criteria"}


def load(path: Path, errors: list):
    try:
        return loads_strict(path.read_text())
    except ValueError as exc:
        errors.append(f"{path.name}: not strict JSON ({exc})")
        return None


def check_example(e: dict, errors: list, seen: set, run_ids: set):
    missing = EXAMPLE_KEYS - set(e)
    if missing:
        errors.append(f"index entry {e.get('file')}: missing keys {sorted(missing)}")
        return
    seen.add(e["file"])
    path = RES / e["file"]
    if not path.exists():
        errors.append(f"index entry {e['file']}: file not found")
        return
    rec = load(path, errors)
    if rec is None:
        return
    for k in ("type", "commit", "dirty", "gpus", "example", "config", "run_id"):
        if rec.get(k) != e.get(k):
            errors.append(f"{e['file']}: {k} differs between index ({e.get(k)}) and file ({rec.get(k)})")
    rid = e.get("run_id")
    if rid is not None:  # records made before run ids have none
        if not e["file"].endswith(f"-run-{rid}.json"):
            errors.append(f"{e['file']}: file name does not carry run id {rid}")
        if rid in run_ids:
            errors.append(f"{e['file']}: run id {rid} is recorded more than once")
        run_ids.add(rid)
    if not isinstance(rec.get("metrics"), dict) or not rec["metrics"]:
        errors.append(f"{e['file']}: no metrics")


def main() -> int:
    index_path = RES / "index.json"
    if not index_path.exists():
        print("no bench/results/index.json (no GPU evidence recorded yet)")
        return 0
    errors = []
    index = load(index_path, errors) or {}
    seen, run_ids = set(), set()
    for e in index.get("entries", []):
        if e.get("type", "gates") == "example":
            check_example(e, errors, seen, run_ids)
            continue
        missing = ENTRY_KEYS - set(e)
        if missing:
            errors.append(f"index entry {e.get('file')}: missing keys {sorted(missing)}")
            continue
        seen.add(e["file"])
        path = RES / e["file"]
        if not path.exists():
            errors.append(f"index entry {e['file']}: file not found")
            continue
        rec = load(path, errors)
        if rec is None:
            continue
        for k in ("commit", "dirty", "gpus", "all_passed"):
            if rec.get(k) != e[k]:
                errors.append(f"{e['file']}: {k} differs between index ({e[k]}) and file ({rec.get(k)})")
        # a gate record runs in fp64 unless it says precision: mixed (PIC-D11), in file, index and name
        prec = rec.get("precision", "fp64")
        if prec != e.get("precision", "fp64"):
            errors.append(f"{e['file']}: precision differs between index "
                          f"({e.get('precision', 'fp64')}) and file ({prec})")
        if e["file"].endswith("-gpu-gates-mixed.json") != (prec == "mixed"):
            errors.append(f"{e['file']}: file name does not match precision {prec}")
        gates = rec.get("gates", [])
        if not gates:
            errors.append(f"{e['file']}: no gates")
        for g in gates:
            if GATE_KEYS - set(g):
                errors.append(f"{e['file']}: gate {g.get('name')} missing keys {sorted(GATE_KEYS - set(g))}")
        if rec.get("all_passed") != all(g.get("passed") for g in gates):
            errors.append(f"{e['file']}: all_passed does not match the gate verdicts")
    for pattern in ("*-gpu-gates.json", "*-gpu-gates-mixed.json", "*-example-*.json"):
        for f in RES.glob(pattern):
            if f.name not in seen:
                errors.append(f"{f.name}: not listed in index.json")
    for msg in errors:
        print(msg, file=sys.stderr)
    if not errors:
        print(f"bench/results: {len(seen)} indexed result file(s) consistent")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
