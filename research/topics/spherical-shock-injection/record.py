"""Record of a run: runs/<date>-<commit>-<name>.json and runs/index.json (research/README.md, branching rules).

A record holds the commit, config, seed and GPUs of the run, its wall time, the place and size of its output
(which stays out of git), the conservation checks, and the summaries of the analysis files. analysis.commit is
the commit at which the record was written; an analysis file can be older, so its hash and time are kept too.
A record is written from a clean tree, and the commit of the run must be an ancestor of HEAD.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import census
import h5py
import numpy as np

TOPIC = Path(__file__).resolve().parent
LIMITS = {"gauss": 1e-10, "divB": 1e-12, "fill": 0.95}  # numerics-invariants.md; fill: the stop rule of the plan
# summaries taken from the analysis files (small values only; the files stay with the output)
METRICS = {
    "front": ("n_outputs", "bins", "detected_bins", "valid_bins", "usable_bins", "fit_rms", "fit_rms_B",
              "flags_in_valid", "flow_within_30pct"),
    "null": ("T0", "N1_max_per_dump>=1.5"),
    "cadence": ("n_tracked", "n_stay", "coarsest_passing"),
    "census": ("n_dumps", "b_max", "b_cell_max", "tau_dwell", "N3b"),
    "tracks": ("n_tracked", "stride", "all_by_category", "n_gamma_max>=1.25", "n_gamma_max>=1.5",
               "gamma_max>=1.5_by_category", "N4", "reflections, all energies", "reflections",
               "K2 several encounters, gamma_out >= 1.5", "K3 superluminal foot, mu not conserved, gamma_out >= 1.5",
               "returned from downstream, gamma_out >= 1.5", "census residuals on tracks (gamma >= 1.5)",
               "K1 on tracks"),
    "n3c": ("null", "commit_null", "largest ratio", "note on N3 and N5"),
    "phi": ("null", "bins on the detected front", "measured (error below 50 %)", "bins of the valid domain",
            "measured, valid domain", "significantly negative", "mean over the outputs with a detected front"),
    "baseline": ("phi", "phi_cells", "N3", "share reflected", "rule gamma >= 1.5", "rule gamma 1.25-1.5 (extra)"),
}


def _git(*args) -> str:
    return subprocess.run(["git", *args], cwd=TOPIC, capture_output=True, text=True, check=True).stdout.strip()


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _size(path: Path) -> int:
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())


def scalar_checks(scalars: dict) -> dict:
    """Largest value of each checked scalar against its limit, and the change of the energy budget (monitored
    only in the explicit scheme)."""
    out = {}
    for key, limit in LIMITS.items():
        if key in scalars and len(scalars[key]):
            worst = float(np.nanmax(np.abs(scalars[key])))
            out[key] = {"max": worst, "limit": limit, "ok": bool(worst < limit if key == "fill" else worst <= limit)}
    if "budget" in scalars and len(scalars["budget"]) > 1 and scalars["budget"][0] != 0:
        b = np.asarray(scalars["budget"], float)
        out["energy budget, relative change (monitored)"] = float((b[-1] - b[0]) / b[0])
    return out


def census_summary(cen: dict) -> dict:
    """Largest numbers of the census over the dumps: upstream electrons by result, and K1."""
    out = {}
    for row in cen.get("upstream", []):
        d = out.setdefault("upstream, largest dump, " + row["domain"], {"t": None, "n": -1, "results": {}})
        n = int(sum(row["results"].values()))
        if n > d["n"]:
            d.update({"t": row["t"], "n": n, "results": row["results"]})
    k1 = [(r["whole front"].get("k1", 0), r["whole front"].get("expected", 0.0), r["t"]) for r in cen.get("K1", [])]
    if k1:
        n, expected, t = max(k1)
        out["K1, largest dump, whole front"] = {"t": t, "k1": n, "expected": expected}
    return out


def build(run_dir: Path, name: str, root: Path) -> dict:
    run = census._open(run_dir)
    info, cfg = run.info, run.config
    meta = info.get("meta", {})
    commit = info.get("git_commit")
    head = _git("rev-parse", "HEAD")
    ancestor = subprocess.run(["git", "merge-base", "--is-ancestor", commit, head], cwd=TOPIC).returncode == 0
    config = meta.get("config_file")
    with h5py.File(run_dir / "scalars.h5", "r") as f:
        scalars = {k: f[k][...] for k in ("gauss", "divB", "fill", "budget", "n_particles") if k in f}
    params = meta.get("params", {})
    rec = {
        "name": name,
        "created": info.get("created"),
        "run": {
            "commit": commit, "commit is an ancestor of the commit of this record": ancestor,
            "example": meta.get("example"), "config": config,
            "config_sha256": _sha(root / config) if config and (root / config).exists() else None,
            "seed": params.get("seed", cfg.get("seed")), "ppc": params.get("ppc"),
            "mass_ratio": params.get("mass_ratio"), "wpe_wce": params.get("wpe_wce"),
            "mach_alfven": params.get("mach_alfven"),
            "species": [sp.get("name") if isinstance(sp, dict) else sp for sp in info.get("species", [])],
            "n_cells": list(info.get("n_cells", [])), "steps": info.get("n_steps"), "dt": info.get("dt"),
            "gpus": meta.get("run", {}).get("gpus", info.get("world_size")), "gpu_name": info.get("gpu_name"),
            "wall_time_s": meta.get("run", {}).get("wall_time_s"),
            "ns_per_particle_step": meta.get("run", {}).get("ns_per_particle_step"),
            "particles_at_the_end": int(scalars["n_particles"][-1]) if "n_particles" in scalars else None,
            "run_id": info.get("run_id"), "pypic_version": info.get("pypic_version"),
        },
        "output": {"path": str(run_dir), "bytes": _size(run_dir), "in git": False},
        "checks": scalar_checks(scalars),
        "analysis": {"commit": head, "files": {}, "metrics": {}},
    }
    for path in sorted((run_dir / "analysis").glob("*.json")) if (run_dir / "analysis").exists() else []:
        data = json.loads(path.read_text())
        rec["analysis"]["files"][path.name] = {
            "sha256": _sha(path), "written": datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds")}
        keep = {k: data[k] for k in METRICS.get(path.stem, ()) if k in data}
        if path.stem == "census":
            keep.update(census_summary(data))
        if keep:
            rec["analysis"]["metrics"][path.stem] = keep
    return rec


def cmd_record(args) -> int:
    root = Path(_git("rev-parse", "--show-toplevel"))
    if _git("status", "--porcelain"):
        print("the tree is not clean: commit the analysis code before a record is written")
        return 1
    run_dir = Path(args.run)
    rec = build(run_dir, args.name, root)
    if not rec["run"]["commit is an ancestor of the commit of this record"]:
        print(f"the commit of the run {rec['run']['commit']} is not an ancestor of HEAD")
        return 1
    runs = TOPIC / "runs"
    runs.mkdir(exist_ok=True)
    date = str(rec["created"])[:10]
    path = runs / f"{date}-{rec['run']['commit'][:9]}-{args.name}.json"
    text = json.dumps(rec, indent=1)
    if len(text) > 1_000_000:
        print(f"the record is {len(text)} bytes: above the limit of 1 MB of research/")
        return 1
    path.write_text(text + "\n")
    index_path = runs / "index.json"
    index = json.loads(index_path.read_text()) if index_path.exists() else []
    entry = {"file": path.name, "name": args.name, "created": rec["created"], "commit": rec["run"]["commit"],
             "config": rec["run"]["config"], "wall_time_s": rec["run"]["wall_time_s"],
             "output_bytes": rec["output"]["bytes"],
             "checks_ok": all(c["ok"] for c in rec["checks"].values() if isinstance(c, dict)),
             "analysis_commit": rec["analysis"]["commit"]}
    index = [e for e in index if e["file"] != path.name] + [entry]
    index.sort(key=lambda e: (e["created"] or "", e["name"]))
    index_path.write_text(json.dumps(index, indent=1) + "\n")
    print(f"wrote {path.relative_to(root)} ({len(text)} bytes); checks "
          + ", ".join(f"{k} {c['max']:.2e} ({'ok' if c['ok'] else 'FAIL'})" for k, c in rec["checks"].items()
                      if isinstance(c, dict)))
    return 0


# ------------------------------------------------------------------ check on synthetic data
def selfcheck(check) -> None:
    good = scalar_checks({"gauss": np.array([1e-15, 5e-15]), "divB": np.array([2e-14]), "fill": np.array([0.4, 0.67]),
                          "budget": np.array([100.0, 110.0])})
    bad = scalar_checks({"gauss": np.array([1e-15, 3e-10]), "divB": np.array([np.nan, 2e-12]),
                         "fill": np.array([0.95])})
    cen = census_summary({"upstream": [{"t": 1.0, "domain": "valid domain", "results": {"S0": 2, "residual": 0}},
                                       {"t": 2.0, "domain": "valid domain", "results": {"S0": 5, "residual": 1}}],
                          "K1": [{"t": 1.0, "whole front": {"k1": 3, "expected": 1.5}},
                                 {"t": 2.0, "whole front": {"k1": 1, "expected": 0.2}}]})
    ok = all(c["ok"] for c in good.values() if isinstance(c, dict)) \
        and abs(good["energy budget, relative change (monitored)"] - 0.1) < 1e-12 \
        and not any(c["ok"] for c in bad.values()) \
        and cen["upstream, largest dump, valid domain"]["n"] == 6 and cen["K1, largest dump, whole front"]["t"] == 1.0
    check("record, checks of the scalars and the summary of the census", ok,
          "limits of the Gauss law, div B and the fill; the largest dump of the census")
