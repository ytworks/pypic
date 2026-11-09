#!/usr/bin/env python3
"""GPU verification gates. Run on a GPU machine before merging any change to the engine.

Gates (criteria are the numerical contracts in .claude/rules/numerics-invariants.md):
  gpu-tests              pytest tests/gpu (includes the 1-vs-2-GPU equivalence test)
  two-stream:<config>    the 3-D two-stream examples on N GPUs: growth rate of the k_max mode vs
                         linear theory, Gauss law, div B and (implicit) total energy

Wall time and throughput are recorded as diagnostics only; they never decide a verdict.
Results go to bench/results/<date>-<commit>-gpu-gates.json and are listed in
bench/results/index.json (validated on CPU CI by scripts/verify_bench_results.py).

    python scripts/gpu_gates.py --gpus 8            # full gates
    python scripts/gpu_gates.py --gpus 2 --quick    # shorter runs while iterating (not recorded)
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "examples" / "two_stream"))

CONFIGS = ["ei", "ep", "eip", "ei_implicit"]
CRITERIA = {
    "explicit": {"growth_rel_error_max": 0.10, "gauss_max": 1e-10, "divB_max": 1e-12},
    "implicit": {"growth_rel_error_max": 0.10, "gauss_max": 1e-10, "divB_max": 1e-12, "energy_error_max": 1e-8},
}


def sh(cmd, **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, **kw)


def git_state() -> tuple[str, bool]:
    commit = sh(["git", "rev-parse", "HEAD"]).stdout.strip()
    dirty = bool(sh(["git", "status", "--porcelain", "--untracked-files=no"]).stdout.strip())
    return commit, dirty


def gate_gpu_tests() -> dict:
    r = sh([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "tests/gpu"], timeout=7200)
    m = re.search(r"(\d+) passed", r.stdout)
    f = re.search(r"(\d+) failed", r.stdout)
    s = re.search(r"(\d+) skipped", r.stdout)
    metrics = {"passed": int(m.group(1)) if m else 0, "failed": int(f.group(1)) if f else 0,
               "skipped": int(s.group(1)) if s else 0}
    ok = r.returncode == 0 and metrics["failed"] == 0 and metrics["skipped"] == 0
    tail = "\n".join(r.stdout.strip().splitlines()[-15:])
    return {"name": "gpu-tests", "passed": ok, "metrics": metrics, "criteria": {"failed": 0, "skipped": 0},
            "log_tail": None if ok else tail}


def gate_two_stream(cfg: str, gpus: int, quick: bool, workdir: Path) -> dict:
    from analyze import analyze

    out = workdir / cfg
    cmd = [sys.executable, "-m", "pypic.launch", "-n", str(gpus), "examples/two_stream/run.py",
           "--config", f"examples/two_stream/configs/{cfg}.yaml", "--out", str(out), "--progress", "1000000"]
    if quick:
        cmd += ["--set", "growth_times=15"]
    r = sh(cmd, timeout=7200)
    name = f"two-stream:{cfg}"
    if r.returncode != 0:
        return {"name": name, "passed": False, "metrics": {}, "criteria": {},
                "log_tail": "\n".join((r.stdout + r.stderr).strip().splitlines()[-20:])}
    res = analyze(str(out), show_plot=False)
    meta = json.loads((out / "meta.json").read_text())
    solver = "implicit" if "implicit" in cfg else "explicit"
    crit = CRITERIA[solver]
    metrics = {
        "growth_fit": res["gamma_fit_kmax"], "growth_theory": res["gamma_theory"],
        "growth_rel_error": res["relative_error"], "gauss_max": res["gauss_max"], "divB_max": res["divB_max"],
        "energy_error_max": res["energy_error_max"],
        # diagnostics (never part of the verdict)
        "wall_time_s": meta.get("wall_time_s"), "ns_per_particle_step_x_gpus": meta.get("ns_per_particle_step"),
    }
    checks = {
        "growth_rel_error_max": abs(metrics["growth_rel_error"]) <= crit["growth_rel_error_max"],
        "gauss_max": metrics["gauss_max"] <= crit["gauss_max"],
        "divB_max": metrics["divB_max"] <= crit["divB_max"],
    }
    if "energy_error_max" in crit:
        checks["energy_error_max"] = metrics["energy_error_max"] <= crit["energy_error_max"]
    return {"name": name, "passed": all(checks.values()), "metrics": metrics, "criteria": crit,
            "failed_checks": [k for k, v in checks.items() if not v]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gpus", type=int, default=8)
    ap.add_argument("--quick", action="store_true", help="shorter runs; results are printed, not recorded")
    ap.add_argument("--skip-tests", action="store_true", help="only run the two-stream gates")
    ap.add_argument("--configs", nargs="*", default=CONFIGS)
    a = ap.parse_args()

    import cupy as cp

    ndev = cp.cuda.runtime.getDeviceCount()
    if ndev < a.gpus:
        print(f"need {a.gpus} GPUs, found {ndev}", file=sys.stderr)
        return 2
    gpu_name = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
    commit, dirty = git_state()
    gates = []
    if not a.skip_tests:
        print("gate gpu-tests ...", flush=True)
        gates.append(gate_gpu_tests())
        print(f"  -> {'PASS' if gates[-1]['passed'] else 'FAIL'} {gates[-1]['metrics']}", flush=True)
    with tempfile.TemporaryDirectory(prefix="pypic-gates-") as td:
        for cfg in a.configs:
            print(f"gate two-stream:{cfg} on {a.gpus} GPU(s) ...", flush=True)
            g = gate_two_stream(cfg, a.gpus, a.quick, Path(td))
            gates.append(g)
            print(f"  -> {'PASS' if g['passed'] else 'FAIL'} "
                  + json.dumps({k: g["metrics"].get(k) for k in ("growth_rel_error", "gauss_max", "divB_max",
                                                                  "energy_error_max")}), flush=True)
    all_passed = all(g["passed"] for g in gates)
    record = {"date": dt.date.today().isoformat(), "commit": commit, "dirty": dirty, "gpus": a.gpus,
              "gpu_name": gpu_name, "quick": a.quick, "all_passed": all_passed, "gates": gates}
    print(f"ALL {'PASSED' if all_passed else 'FAILED'}")
    if a.quick or a.skip_tests or a.configs != CONFIGS:
        print("(partial/quick run: not recorded in bench/results)")
        return 0 if all_passed else 1
    res_dir = ROOT / "bench" / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    fname = f"{record['date']}-{commit[:9]}{'-dirty' if dirty else ''}-gpu-gates.json"
    (res_dir / fname).write_text(json.dumps(record, indent=2, default=float) + "\n")
    index_path = res_dir / "index.json"
    index = json.loads(index_path.read_text()) if index_path.exists() else {"entries": []}
    index["entries"] = [e for e in index["entries"] if e["file"] != fname]
    index["entries"].insert(0, {"file": fname, "date": record["date"], "commit": commit, "dirty": dirty,
                                "gpus": a.gpus, "gpu_name": gpu_name, "all_passed": all_passed})
    index_path.write_text(json.dumps(index, indent=2) + "\n")
    print(f"recorded bench/results/{fname}")
    if dirty:
        print("NOTE: the working tree was dirty; re-run on a clean commit for evidence that counts.")
    return 0 if all_passed else 1


if __name__ == "__main__":
    os.environ.setdefault("PYTHONUNBUFFERED", "1")
    sys.exit(main())
