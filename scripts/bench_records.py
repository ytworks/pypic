"""Shared helpers for the GPU evidence in bench/results/ (gpu_gates.py, record_example.py,
verify_bench_results.py).

Records are strict JSON: a quantity that could not be measured (NaN, an infinity) is stored as
null, so any JSON parser reads the files, not only Python's.
"""
from __future__ import annotations

import json
from pathlib import Path


def finite(x):
    """Copy of x with numpy values made plain Python and every non-finite number replaced by None."""
    if isinstance(x, dict):
        return {str(k): finite(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [finite(v) for v in x]
    if hasattr(x, "tolist") and not isinstance(x, (str, bytes)):  # numpy scalar or array
        return finite(x.tolist())
    if isinstance(x, float) and (x != x or x in (float("inf"), float("-inf"))):
        return None
    return x


def dumps(record) -> str:
    return json.dumps(finite(record), indent=2, allow_nan=False) + "\n"


def _reject(token):
    raise ValueError(f"non-standard JSON constant {token}")


def loads_strict(text: str):
    """json.loads that rejects NaN / Infinity (Python accepts them, JSON does not)."""
    return json.loads(text, parse_constant=_reject)


def write(res_dir: Path, fname: str, record: dict, entry: dict, replaces=lambda e: False) -> None:
    """Write a record and list it first in index.json. Index entries for the same file, and those
    for which ``replaces(entry)`` is true (their files are deleted), are dropped."""
    res_dir.mkdir(parents=True, exist_ok=True)
    index_path = res_dir / "index.json"
    index = loads_strict(index_path.read_text()) if index_path.exists() else {"entries": []}
    keep = []
    for e in index["entries"]:
        if e["file"] == fname:
            continue
        if replaces(e):
            (res_dir / e["file"]).unlink(missing_ok=True)
            continue
        keep.append(e)
    (res_dir / fname).write_text(dumps(record))
    index["entries"] = [{"file": fname, **entry}] + keep
    index_path.write_text(dumps(index))
