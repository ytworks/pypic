#!/usr/bin/env python3
"""GPU verification gates. Run on a GPU machine before merging any change to the engine.

Gates (criteria are the numerical contracts in .claude/rules/numerics-invariants.md):
  gpu-tests              pytest tests/gpu (includes the 1-vs-2-GPU equivalence test)
  <example>:<config>     analytic-solution examples (plasma waves, Landau / Dawson, filamentation, ...):
                         each examples/<name>/analyze.py verdict vs theory, and Gauss's law
  two-stream:<config>    the 3-D two-stream / Buneman examples on N GPUs: growth rate of the k_max mode vs
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
import signal
import subprocess
import sys
import tempfile
from pathlib import Path

import bench_records  # scripts/bench_records.py (this directory is on sys.path)

from pypic.provenance import code_state

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "examples" / "two_stream"))

CONFIGS = ["ei", "ep", "eip", "ei_implicit", "buneman", "ei_b0"]
# precision: mixed (PIC-D11) applies to the explicit solver only; its record is the explicit gates
MIXED_CONFIGS = [c for c in CONFIGS if "implicit" not in c]
PRECISION = "fp64"


def precision_args() -> list:
    return ["--set", "overrides.precision=mixed"] if PRECISION == "mixed" else []
# analytic-solution examples: (example directory, config). Each examples/<name>/analyze.py returns a
# verdict ("passed") against its theory with the criteria documented in its README and in
# docs/design/output-analysis.md; Gauss's law is checked here for every one of them.
EXAMPLES = [
    ("plasma_waves", "unmagnetized"), ("plasma_waves", "perpendicular"),
    ("landau_damping", "landau"), ("landau_damping", "dawson"), ("landau_damping", "wavebreaking"),
    ("landau_damping", "relativistic"),
    ("filamentation", "pair"), ("filamentation", "pair_b0"), ("filamentation", "pair_b0_stable"),
    ("particle_orbits", "gyration"), ("particle_orbits", "exb"), ("particle_orbits", "hyperbolic"),
    ("particle_orbits", "e_gt_b"), ("particle_orbits", "plane_wave"),
    ("plasma_expansion", "isothermal_hydrogen"), ("relativistic_shock", "pair_2d"),
    ("relativistic_shock", "pair_2d_sigma"), ("reconnection", "pair_sigma10"),
    ("bow_shock", "gate"), ("bow_shock", "south"), ("bow_shock", "oblique"), ("blast_wave", "gate"),
]
EXAMPLE_NAMES = [f"{e}:{c}" for e, c in EXAMPLES]
# examples verified by a convergence check: run again with dt_scale=0.5, on one GPU (few particles)
REFINED = {"particle_orbits"}
CRITERIA = {
    "explicit": {"growth_rel_error_max": 0.10, "gauss_max": 1e-10, "divB_max": 1e-12},
    "implicit": {"growth_rel_error_max": 0.10, "gauss_max": 1e-10, "divB_max": 1e-12, "energy_error_max": 1e-8},
}


def sh(cmd, timeout=None) -> subprocess.CompletedProcess:
    """Run cmd in its own process group; a timeout kills the whole group (a hung launcher's ranks
    would otherwise keep spinning on the GPUs)."""
    p = subprocess.Popen(cmd, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         start_new_session=True)
    try:
        out, err = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        out, err = p.communicate()
        return subprocess.CompletedProcess(cmd, -9, out, err + f"\nTIMEOUT after {timeout} s")
    return subprocess.CompletedProcess(cmd, p.returncode, out, err)


def git_state() -> tuple[str, bool]:
    """Commit and whether the code differs from it, untracked code files included (pypic.provenance)."""
    commit, changes = code_state(ROOT)
    return commit, bool(changes)


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
           "--config", f"examples/two_stream/configs/{cfg}.yaml", "--out", str(out), "--progress", "1000000",
           *precision_args()]
    if quick:
        cmd += ["--set", "growth_times=15"]
    r = sh(cmd, timeout=7200)
    name = f"two-stream:{cfg}"
    if r.returncode != 0:
        return {"name": name, "passed": False, "metrics": {}, "criteria": {},
                "log_tail": "\n".join((r.stdout + r.stderr).strip().splitlines()[-20:])}
    res = analyze(str(out), show_plot=False)
    meta = json.loads((out / "run.json").read_text())["meta"].get("run", {})
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


def gate_example(example: str, cfg: str, gpus: int, workdir: Path) -> dict:
    name = f"{example}:{cfg}"
    out = workdir / example / cfg
    refine = example in REFINED
    runs = [(out, [])] + ([(workdir / example / f"{cfg}_half", ["--set", "dt_scale=0.5"])] if refine else [])
    for o, extra in runs:
        cmd = [sys.executable, "-m", "pypic.launch", "-n", str(1 if refine else gpus), f"examples/{example}/run.py",
               "--config", f"examples/{example}/configs/{cfg}.yaml", "--out", str(o), "--progress", "1000000", *extra,
               *precision_args()]
        r = sh(cmd, timeout=7200)
        if r.returncode != 0:
            return {"name": name, "passed": False, "metrics": {}, "criteria": {},
                    "log_tail": "\n".join((r.stdout + r.stderr).strip().splitlines()[-20:])}
    # each analysis in its own interpreter (examples have their own theory / plan modules)
    kw = f", refined={str(runs[1][0])!r}" if refine else ""
    a = sh([sys.executable, "-c", "import json, sys; sys.path.insert(0, sys.argv[1]); from analyze import analyze; "
            f"print(json.dumps(analyze(sys.argv[2], show_plot=False{kw}), default=float))",
            f"examples/{example}", str(out)], timeout=1800)
    if a.returncode != 0:
        return {"name": name, "passed": False, "metrics": {}, "criteria": {},
                "log_tail": "\n".join((a.stdout + a.stderr).strip().splitlines()[-20:])}
    res = json.loads(a.stdout.strip().splitlines()[-1])
    meta = json.loads((out / "run.json").read_text())["meta"].get("run", {})
    res.update(wall_time_s=meta.get("wall_time_s"), ns_per_particle_step_x_gpus=meta.get("ns_per_particle_step"))
    checks = {"analysis": bool(res.get("passed")), "gauss_max": res.get("gauss_max", 0.0) <= 1e-10}
    return {"name": name, "passed": all(checks.values()), "metrics": res,
            "criteria": {"analysis": "examples/" + example + "/README.md", "gauss_max": 1e-10},
            "failed_checks": [k for k, v in checks.items() if not v]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gpus", type=int, default=8)
    ap.add_argument("--quick", action="store_true", help="shorter runs; results are printed, not recorded")
    ap.add_argument("--skip-tests", action="store_true", help="only run the two-stream gates")
    ap.add_argument("--configs", nargs="*", default=CONFIGS)
    ap.add_argument("--examples", nargs="*", default=EXAMPLE_NAMES, help="example:config gates to run")
    ap.add_argument("--precision", choices=["fp64", "mixed"], default="fp64",
                    help="mixed: the explicit gates with precision: mixed (PIC-D11), without the GPU tests "
                         "(precision-independent, part of the fp64 record); recorded as *-gpu-gates-mixed.json")
    a = ap.parse_args()
    global PRECISION
    PRECISION = a.precision
    full_configs = MIXED_CONFIGS if a.precision == "mixed" else CONFIGS
    if a.precision == "mixed":
        a.skip_tests = True
        if a.configs == CONFIGS:
            a.configs = MIXED_CONFIGS

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
        for ex in a.examples:
            example, cfg = ex.split(":")
            print(f"gate {ex} on {a.gpus} GPU(s) ...", flush=True)
            g = gate_example(example, cfg, a.gpus, Path(td))
            gates.append(g)
            print(f"  -> {'PASS' if g['passed'] else 'FAIL'} {g.get('failed_checks') or ''}", flush=True)
    all_passed = all(g["passed"] for g in gates)
    mixed = a.precision == "mixed"
    record = {"date": dt.date.today().isoformat(), "commit": commit, "dirty": dirty, "gpus": a.gpus,
              "gpu_name": gpu_name, "quick": a.quick, "all_passed": all_passed, "gates": gates}
    if mixed:
        record["precision"] = "mixed"
    print(f"ALL {'PASSED' if all_passed else 'FAILED'}")
    if a.quick or (a.skip_tests and not mixed) or a.configs != full_configs or a.examples != EXAMPLE_NAMES:
        print("(partial/quick run: not recorded in bench/results)")
        return 0 if all_passed else 1
    fname = f"{record['date']}-{commit[:9]}{'-dirty' if dirty else ''}-gpu-gates{'-mixed' if mixed else ''}.json"
    entry = {"date": record["date"], "commit": commit, "dirty": dirty, "gpus": a.gpus, "gpu_name": gpu_name,
             "all_passed": all_passed}
    if mixed:
        entry["precision"] = "mixed"
    # strict JSON: a metric that could not be measured is null (bench_records)
    bench_records.write(ROOT / "bench" / "results", fname, record, entry)
    print(f"recorded bench/results/{fname}")
    if dirty:
        print("NOTE: the working tree was dirty; re-run on a clean commit for evidence that counts.")
    return 0 if all_passed else 1


if __name__ == "__main__":
    os.environ.setdefault("PYTHONUNBUFFERED", "1")
    sys.exit(main())
