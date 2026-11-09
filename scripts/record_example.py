#!/usr/bin/env python3
"""Record the analysis of a finished example run (typically a production configuration, whose
results are reported rather than gated) as evidence in bench/results/.

    python scripts/record_example.py output/relativistic_shock/pair_2d_long

The example's own analyze.py runs on the output directory (in a separate interpreter); its result,
the run's provenance from run.json (code commit at run time, dirty flag, GPUs, wall time), the
commit of the analysis code and the configuration are written to
bench/results/<date>-<commit>-example-<example>-<config>.json and listed in
bench/results/index.json with type "example" (validated by verify_bench_results.py).
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "bench" / "results"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--note", default="", help="free text stored with the record")
    a = ap.parse_args()
    run_dir = Path(a.run_dir).resolve()
    info = json.loads((run_dir / "run.json").read_text())
    meta = info.get("meta", {})
    example = meta.get("example")
    if not example:
        print("run.json has no meta.example: was the run started by an example's run.py?", file=sys.stderr)
        return 2
    config = Path(meta.get("config_file", "unknown")).stem
    code = ("import json, sys; sys.path.insert(0, sys.argv[1]); from analyze import analyze; "
            "print(json.dumps(analyze(sys.argv[2], show_plot=False), default=float))")
    r = subprocess.run([sys.executable, "-c", code, str(ROOT / "examples" / example), str(run_dir)], cwd=ROOT,
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-3000:], file=sys.stderr)
        return 1
    result = json.loads(r.stdout.strip().splitlines()[-1])
    commit_full = info.get("git_commit") or "unknown"
    dirty = commit_full.endswith("-dirty")
    commit = commit_full.removesuffix("-dirty")
    date = dt.date.today().isoformat()
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    an_dirty = bool(subprocess.run(["git", "status", "--porcelain", "--untracked-files=no", "examples", "pypic"],
                                   cwd=ROOT, capture_output=True, text=True).stdout.strip())
    record = {"type": "example", "date": date, "commit": commit, "dirty": dirty, "gpus": info.get("world_size"),
              "analysis_commit": head + ("-dirty" if an_dirty else ""),
              "gpu_name": info.get("gpu_name"), "example": example, "config": config, "params": meta.get("params"),
              "run": meta.get("run", {}), "passed": result.get("passed"), "metrics": result, "note": a.note}
    fname = f"{date}-{commit[:9]}{'-dirty' if dirty else ''}-example-{example}-{config}.json"
    RES.mkdir(parents=True, exist_ok=True)
    (RES / fname).write_text(json.dumps(record, indent=2, default=float) + "\n")
    index_path = RES / "index.json"
    index = json.loads(index_path.read_text()) if index_path.exists() else {"entries": []}
    index["entries"] = [e for e in index["entries"] if e["file"] != fname]
    index["entries"].insert(0, {"file": fname, "type": "example", "date": date, "commit": commit, "dirty": dirty,
                                "gpus": record["gpus"], "gpu_name": record["gpu_name"], "example": example,
                                "config": config})
    index_path.write_text(json.dumps(index, indent=2) + "\n")
    print(f"recorded bench/results/{fname}")
    if dirty:
        print("NOTE: the run was made from a dirty working tree.")
    return 0


if __name__ == "__main__":
    os.environ.setdefault("PYTHONUNBUFFERED", "1")
    sys.exit(main())
