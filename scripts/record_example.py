#!/usr/bin/env python3
"""Record the analysis of a finished example run (typically a production configuration, whose
results are reported rather than gated) as evidence in bench/results/.

    python scripts/record_example.py output/relativistic_shock/pair_2d_long

The example's own analyze.py runs on the output directory (in a separate interpreter); its result,
the run's provenance from run.json (run id, code commit at run time, dirty flag, GPUs, wall time),
the commit of the analysis code and the configuration are written to
bench/results/<date>-<commit>-example-<example>-<config>-run-<run_id>.json and listed in
bench/results/index.json with type "example" (validated by verify_bench_results.py).

- One record per run: separate runs of a configuration (realizations) are kept side by side;
  recording the same run again replaces its record.
- The analysis commit is taken before the analysis code is loaded, "-dirty" when the code differs
  from it, untracked code files included (pypic.provenance); if the code changes while the
  analysis runs, nothing is recorded.
- The record is strict JSON: a metric that could not be measured (NaN) is stored as null.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
from pathlib import Path

import bench_records  # scripts/bench_records.py (this directory is on sys.path)

from pypic.provenance import code_state

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
    run_id = info.get("run_id")
    if not isinstance(run_id, str) or not run_id.isalnum():
        print("run.json has no run_id (a run made before run ids were written): records need it to keep "
              "separate runs apart", file=sys.stderr)
        return 2
    # run.json is data: its example name only selects one of this repository's examples, whose
    # analyze.py is loaded by its path (nothing from the run directory is put on sys.path)
    known = {p.name: p / "analyze.py" for p in (ROOT / "examples").iterdir() if (p / "analyze.py").is_file()}
    if example not in known:
        print(f"meta.example {example!r} is not an example of this repository ({sorted(known)})", file=sys.stderr)
        return 2
    config = Path(meta.get("config_file", "unknown")).stem
    code = ("import importlib.util, json, sys\n"
            "spec = importlib.util.spec_from_file_location('analyze', sys.argv[1])\n"
            "mod = importlib.util.module_from_spec(spec)\n"
            "spec.loader.exec_module(mod)\n"
            "print(json.dumps(mod.analyze(sys.argv[2], show_plot=False), default=float))")
    # the analysis code is identified before it is loaded and must be unchanged when it ends
    before = code_state(ROOT)
    r = subprocess.run([sys.executable, "-c", code, str(known[example]), str(run_dir)], cwd=ROOT,
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-3000:], file=sys.stderr)
        return 1
    if code_state(ROOT) != before:
        print("the analysis code (HEAD, or tracked or untracked files under the code paths) changed while "
              "the analysis ran: nothing recorded", file=sys.stderr)
        return 1
    result = json.loads(r.stdout.strip().splitlines()[-1])
    commit_full = info.get("git_commit") or "unknown"
    dirty = commit_full.endswith("-dirty")
    commit = commit_full.removesuffix("-dirty")
    date = dt.date.today().isoformat()
    head, changes = before
    record = {"type": "example", "date": date, "commit": commit, "dirty": dirty, "run_id": run_id,
              "run_started": info.get("created"), "gpus": info.get("world_size"),
              "analysis_commit": head + ("-dirty" if changes else ""),
              "gpu_name": info.get("gpu_name"), "example": example, "config": config, "params": meta.get("params"),
              "run": meta.get("run", {}), "passed": result.get("passed"), "metrics": result, "note": a.note}
    fname = f"{date}-{commit[:9]}{'-dirty' if dirty else ''}-example-{example}-{config}-run-{run_id}.json"
    entry = {"type": "example", "date": date, "commit": commit, "dirty": dirty, "run_id": run_id,
             "gpus": record["gpus"], "gpu_name": record["gpu_name"], "example": example, "config": config}
    # one record per run: re-recording a run replaces its earlier record; other runs are kept
    bench_records.write(RES, fname, record, entry, replaces=lambda e: e.get("run_id") == run_id)
    print(f"recorded bench/results/{fname}")
    if dirty:
        print("NOTE: the run was made from a dirty working tree.")
    return 0


if __name__ == "__main__":
    os.environ.setdefault("PYTHONUNBUFFERED", "1")
    sys.exit(main())
