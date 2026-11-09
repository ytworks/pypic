#!/usr/bin/env python3
"""Build the GitHub wiki articles of the examples (English and Japanese) from the runs themselves.

Every number, chart and animation on a page comes from a run's own output: nothing is copied by hand.

    python scripts/build_wiki.py run [--devices 0,1,2,3]       # fresh runs missing under output/wiki/
                                                               #   (from a clean checkout: run.json
                                                               #   records the commit, "-dirty" otherwise)
    python scripts/build_wiki.py collect --wiki ../pypic.wiki  # analysis, charts, GIFs -> <wiki>/images,
                                                               #   numbers -> docs/wiki/data.json
    python scripts/build_wiki.py render --wiki ../pypic.wiki   # pages from docs/wiki/{en,ja} + data.json

What each step needs:
    render   only the checked-in docs/wiki/data.json and templates (the images already in the wiki stay valid).
    collect  the output of every run: the fresh runs that `run` makes, and the recorded production runs
             (EXAMPLES[...]["cases"] with a "path"). Their output directories are not in git: copy them into
             output/ from the machine that made them. collect checks each one against its evidence record in
             bench/results/ (run id, or commit and wall time for records older than run ids) and stops if
             one is missing or is a different run.

Templates (docs/wiki/en/<Page>.md, docs/wiki/ja/<Page>.md; the Japanese page is published as <Page>-ja)
use {{ ... }} placeholders:
    {{ r.<case>.<key>[.<index>] | <filter> | <format spec> }}   analysis result of a case
    {{ c.<case>.<key> }}, {{ p.<case>.<key> }}                  conditions, provenance of a case
    {{ <example>:r.<case>.<key> }}                              the same, from another example (Home)
    {{ table:conditions }}, {{ table:params }}, {{ table:provenance }}
filters: abs, neg, half; the format spec is a Python one (".3f", "+.1%", ".1e"); a missing key is an error.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEMPLATES = ROOT / "docs" / "wiki"
DATA = TEMPLATES / "data.json"
FRESH = "output/wiki"
LANGS = ("en", "ja")

# GIF settings: small enough for a wiki page (target <= 8 MB each)
GIF = {"fps": 10, "n_frames": 90, "dpi": 70}


def fresh(config: str, gpus: int = 4, refined: bool = False) -> dict:
    return {"config": config, "gpus": gpus, "refined": refined}


def recorded(path: str, evidence: str) -> dict:
    return {"path": path, "evidence": evidence}


# example -> page, title, cases (fresh runs of a config, or recorded production runs), animations, summary chart
EXAMPLES: dict[str, dict] = {
    "two_stream": {
        "page": "Two-Stream-Instability",
        "title": {"en": "Two-stream and Buneman instabilities", "ja": "二流体不安定性とブネマン不安定性"},
        "cases": {c: fresh(c) for c in ("ei", "ep", "eip", "ei_implicit", "buneman")},
        "movies": ["ei", "buneman"],
        "summary": {"kind": "runs"},
    },
    "plasma_waves": {
        "page": "Plasma-Waves",
        "title": {"en": "Thermal plasma waves (dispersion relations)", "ja": "熱的プラズマ波（分散関係）"},
        "cases": {c: fresh(c) for c in ("unmagnetized", "perpendicular", "parallel", "ion_acoustic")},
        "movies": ["unmagnetized"],
        "summary": {"kind": "results", "file": "branch_errors.png"},
    },
    "landau_damping": {
        "page": "Landau-Damping",
        "title": {"en": "Langmuir oscillations: Landau damping, Dawson, wave breaking",
                  "ja": "ラングミュア振動：ランダウ減衰・Dawson 振動・波の砕波"},
        "cases": {c: fresh(c) for c in ("landau", "nonlinear", "dawson", "wavebreaking", "relativistic")},
        "movies": ["nonlinear", ("wavebreaking", "sheets_movie")],
        "summary": {"kind": "results", "file": "errors.png"},
    },
    "filamentation": {
        "page": "Filamentation-Instability",
        "title": {"en": "Current filamentation (Weibel) instability",
                  "ja": "電流フィラメンテーション（ワイベル）不安定性"},
        "cases": {c: fresh(c) for c in ("pair", "electron_ion")},
        "movies": ["pair"],
    },
    "particle_orbits": {
        "page": "Particle-Orbits",
        "title": {"en": "Test-particle orbits in prescribed fields", "ja": "与えられた場の中のテスト粒子軌道"},
        "cases": {c: fresh(c, gpus=1, refined=True) for c in ("gyration", "exb", "hyperbolic", "e_gt_b", "plane_wave")},
        "movies": "all",
        "summary": {"kind": "results", "file": "convergence.png"},
    },
    "plasma_expansion": {
        "page": "Plasma-Expansion",
        "title": {"en": "Plasma expansion into vacuum (Mora 2003)", "ja": "真空へのプラズマ膨張（Mora 2003）"},
        "cases": {
            "isothermal_hydrogen": fresh("isothermal_hydrogen"),
            "isothermal": fresh("isothermal"),
            "thin_foil": fresh("thin_foil"),
            "isothermal_hydrogen_long": recorded("output/plasma_expansion/isothermal_hydrogen_long",
                                                 "bench/results/2026-11-09-200af1306-example-plasma_expansion-"
                                                 "isothermal_hydrogen_long.json"),
        },
        "movies": ["isothermal_hydrogen_long"],
    },
    "relativistic_shock": {
        "page": "Relativistic-Shock",
        "title": {"en": "Relativistic pair shock", "ja": "相対論的電子陽電子衝撃波"},
        "cases": {
            "pair_2d": fresh("pair_2d"),
            "pair_2d_long": recorded("output/relativistic_shock/pair_2d_long",
                                     "bench/results/2026-11-09-200af1306-example-relativistic_shock-pair_2d_long.json"),
            "pair_2d_long_run2": recorded("output/relativistic_shock/pair_2d_long_run2",
                                          "bench/results/2026-11-10-20b2e74c1-example-relativistic_shock-"
                                          "pair_2d_long-run-40c0ca85.json"),
            "pair_3d": recorded("output/relativistic_shock/pair_3d",
                                "bench/results/2026-11-09-3bec889ef-example-relativistic_shock-pair_3d.json"),
        },
        "movies": ["pair_2d"],
        "summary": {"kind": "results", "file": "jump_conditions.png"},
    },
    "reconnection": {
        "page": "Magnetic-Reconnection",
        "title": {"en": "Relativistic magnetic reconnection", "ja": "相対論的磁気リコネクション"},
        "cases": {
            "pair_sigma10": fresh("pair_sigma10"),
            "pair_sigma10_large": recorded("output/reconnection/pair_sigma10_large",
                                           "bench/results/2026-11-09-200af1306-example-reconnection-"
                                           "pair_sigma10_large.json"),
            "pair_sigma10_large_run2": recorded("output/reconnection/pair_sigma10_large_run2",
                                                "bench/results/2026-11-10-d9637cd9f-example-reconnection-"
                                                "pair_sigma10_large-run-e879cfd6.json"),
        },
        "movies": ["pair_sigma10_large_run2"],
        "gif": {"dpi": 52},
        "figures": [("pair_sigma10_large_run2", "fieldlines", "fieldlines.png")],
        "summary": {"kind": "results", "file": "power_law.png"},
    },
}


def run_dir(ex: str, case: str, half: bool = False) -> str:
    c = EXAMPLES[ex]["cases"][case]
    if "path" in c:
        return c["path"]
    return f"{FRESH}/{ex}/{case}" + ("_half" if half else "")


# ----------------------------------------------------------------------------------------------- run

def _complete(d: Path) -> bool:
    try:
        return "wall_time_s" in json.loads((d / "run.json").read_text())["meta"].get("run", {})
    except (OSError, ValueError, KeyError):
        return False


def cmd_run(args) -> int:
    devices = args.devices.split(",")
    for ex, spec in EXAMPLES.items():
        if args.example and ex not in args.example:
            continue
        for case, c in spec["cases"].items():
            if "config" not in c:
                continue
            for half in ([False, True] if c["refined"] else [False]):
                out = ROOT / run_dir(ex, case, half)
                if _complete(out):
                    continue
                n = min(c["gpus"], len(devices))
                cmd = [sys.executable, "-m", "pypic.launch", "-n", str(n), "--devices", ",".join(devices[:n]),
                       f"examples/{ex}/run.py", "--config", f"examples/{ex}/configs/{c['config']}.yaml",
                       "--out", str(out), "--overwrite", "--progress", "1000000"]
                if half:
                    cmd += ["--set", "dt_scale=0.5"]
                print(" ".join(cmd), flush=True)
                r = subprocess.run(cmd, cwd=ROOT)
                if r.returncode != 0:
                    return r.returncode
    return 0


# ------------------------------------------------------------------------------------------- collect

def _clean(v):
    """Strict-JSON copy: numpy scalars to Python, NaN/inf to None."""
    if isinstance(v, dict):
        return {str(k): _clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_clean(x) for x in v]
    if hasattr(v, "tolist"):
        return _clean(v.tolist())
    if isinstance(v, float):
        return v if math.isfinite(v) else None
    return v


def _perturbation(p: dict) -> str:
    amp = " ".join(f"{k}={v:g}" if isinstance(v, (int, float)) else f"{k}={v}"
                   for k, v in p.items() if k not in ("axis", "mode"))
    return f"{amp} (mode {p.get('mode')} along {p.get('axis')})"


def _species(sp: dict) -> str:
    q, m = sp["charge"], sp["mass"]
    head = f"{sp['name']} (q={q:+g}, m={m:g})"
    if sp.get("test"):
        return f"{head}: {len(sp.get('particles', []))} test particles"
    beams = []
    for b in sp["beams"]:
        u = b["drift"]
        s = f"n={b['density']:g}"
        if any(abs(x) > 0 for x in u):
            s += ", u=(" + ", ".join(f"{x:.4g}" for x in u) + ")"
        s += f", θ={b['theta']:.3g}"
        if b.get("profile"):
            pr = b["profile"]
            s += f", {pr.get('shape')} along {pr.get('axis')}"
        if b.get("perturbation"):
            s += ", perturbation: " + ", ".join(_perturbation(p) for p in b["perturbation"])
        beams.append(s)
    return head + ": " + "; ".join(beams)


def conditions(d: Path) -> dict:
    """Experimental conditions of a run, from its run.json and the particle count at t = 0."""
    from pypic.analysis import open_run

    r = json.loads((d / "run.json").read_text())
    cfg = r["config"]
    try:
        npart = int(open_run(str(d)).scalars["n_particles"][0])
    except (OSError, KeyError, IndexError):  # test-particle runs write no scalars
        npart = sum(len(s.get("particles", [])) for s in cfg["species"]) or None
    bnd = []
    for ax, (lo, hi) in zip("xyz", cfg["boundaries"]):
        bnd.append(f"{ax}: {lo}" if lo == hi else f"{ax}: {lo}/{hi}")
    filt = cfg.get("current_filter") or {}
    fields = []
    if any(abs(x) > 0 for x in cfg.get("B0", [0, 0, 0])):
        fields.append("B0=(" + ", ".join(f"{x:.4g}" for x in cfg["B0"]) + ")")
    if any(abs(x) > 0 for x in cfg.get("E0", [0, 0, 0])):
        fields.append("E0=(" + ", ".join(f"{x:.4g}" for x in cfg["E0"]) + ")")
    for f in cfg.get("initial_fields") or []:
        fields.append(f"{f.get('component', '')} {f.get('kind', '')} {f.get('shape', '')}".strip())
    meta_run = r["meta"].get("run", {})
    return {
        "solver": r["solver"], "n_cells": r["n_cells"], "L": r["L"], "dx": r["dx"], "dt": r["dt"],
        "n_steps": r["n_steps"], "t_end": r["dt"] * r["n_steps"], "ppc": cfg["ppc"],
        "particles": npart,
        "species": [_species(s) for s in cfg["species"]], "boundaries": bnd,
        "fields": fields, "filter_passes": int(filt.get("passes", 0) or 0),
        "gpus": meta_run.get("gpus", r.get("world_size")), "wall_time_s": meta_run.get("wall_time_s"),
        "params": r["meta"].get("params", {}),
        "planner": {k: v for k, v in r["meta"].get("planner", {}).items() if k != "params"},
    }


def provenance(d: Path) -> dict:
    r = json.loads((d / "run.json").read_text())
    return {"git_commit": r.get("git_commit"), "run_id": r.get("run_id"), "created": r.get("created"),
            "gpu_name": r.get("gpu_name"), "world_size": r.get("world_size"), "run_dir": str(d.relative_to(ROOT))}


def collect_example(ex: str, wiki: Path) -> dict:
    """Analyse every case of one example (run in its own interpreter: examples have their own modules)."""
    sys.path.insert(0, str(ROOT / "examples" / ex))
    import analyze as A  # noqa: E402  (examples/<ex>/analyze.py)

    spec = EXAMPLES[ex]
    img = wiki / "images" / ex
    img.mkdir(parents=True, exist_ok=True)
    out = {"page": spec["page"], "title": spec["title"], "cases": {}, "images": {}}
    runs = {case: ROOT / run_dir(ex, case) for case in spec["cases"]}
    results = {}
    for case, d in runs.items():
        kw = {"refined": str(ROOT / run_dir(ex, case, half=True))} if spec["cases"][case].get("refined") else {}
        res = A.analyze(str(d), show_plot=True, **kw)
        plots = [res.pop("plot")] if "plot" in res else res.pop("plots", [])
        for p in plots:
            name = f"{case}-{Path(p).name}"
            shutil.copy(p, img / name)
            out["images"][name] = (img / name).stat().st_size
        results[case] = res
        c = spec["cases"][case]
        out["cases"][case] = {"source": "recorded" if "path" in c else "fresh", "evidence": c.get("evidence"),
                              "results": _clean(res), "conditions": _clean(conditions(d)),
                              "provenance": provenance(d)}
        if c.get("refined"):
            out["cases"][case]["provenance_refined"] = provenance(ROOT / run_dir(ex, case, half=True))
    summ = spec.get("summary")
    if summ:
        if summ["kind"] == "runs":
            paths = A.summary({k: str(v) for k, v in runs.items()}, results, str(img))
        else:
            paths = [A.summary(results, str(img / summ["file"]))]
        for p in paths:
            out["images"][Path(p).name] = Path(p).stat().st_size
    for case, fn, name in spec.get("figures", []):  # extra charts: (case, function of the example, file)
        p = getattr(A, fn)(str(runs[case]), str(img / name))
        out["images"][Path(p).name] = Path(p).stat().st_size
    movies = spec.get("movies", [])
    if movies == "all":
        p = A.movie({k: str(v) for k, v in runs.items()}, str(img / "orbits.gif"))
        out["images"][Path(p).name] = Path(p).stat().st_size
    else:
        for m in movies:  # a case, or (case, name of the movie function of the example)
            case, fn = (m, "movie") if isinstance(m, str) else m
            p = getattr(A, fn)(str(runs[case]), out=str(img / f"{case}.gif"), **{**GIF, **spec.get("gif", {})})
            out["images"][Path(p).name] = Path(p).stat().st_size
    return out


def check_recorded(ex: str, case: str, c: dict) -> str | None:
    """None when the output directory of a recorded case is the run of its evidence record, else why not."""
    d = ROOT / c["path"]
    try:
        run = json.loads((d / "run.json").read_text())
    except (OSError, ValueError):
        return (f"{ex}:{case}: no run output at {c['path']} (the run of {c['evidence']}); copy that run's "
                f"output directory there")
    ev = json.loads((ROOT / c["evidence"]).read_text())
    config = Path(run["meta"].get("config_file", "")).stem
    if ev.get("example") != ex or ev.get("config") != config:
        return f"{ex}:{case}: {c['path']} is config {config!r}, the evidence is {ev.get('example')}:{ev.get('config')}"
    if ev.get("run_id"):
        same = run.get("run_id") == ev["run_id"]
    else:  # records made before run ids: the commit and the wall time identify the run
        same = (run.get("run_id") is None and (run.get("git_commit") or "").removesuffix("-dirty") == ev["commit"]
                and run["meta"].get("run", {}).get("wall_time_s") == ev["run"]["wall_time_s"])
    if not same:
        return (f"{ex}:{case}: {c['path']} is not the run of {c['evidence']} (run id "
                f"{run.get('run_id')}, commit {run.get('git_commit')}); restore that run's output")
    return None


def merge_collected(data: dict, collected: dict, before: tuple, after: tuple, today: str) -> dict:
    """data.json with the newly collected examples, each stamped with the code state of its analysis.
    ``before`` and ``after`` are pypic.provenance.code_state() around the analysis: if the code changed
    meanwhile the results belong to no single commit, and nothing is merged (RuntimeError)."""
    if before != after:
        raise RuntimeError("the code changed while collecting: nothing written; collect again")
    data = {**data, "examples": dict(data.get("examples", {}))}
    # older data.json files have one analysis commit for every example: hand it to each before dropping it
    legacy = {k: data.pop(k) for k in ("analysis_commit", "generated") if k in data}
    for ex, e in data["examples"].items():
        data["examples"][ex] = {**legacy, **e}
    commit = before[0] + ("-dirty" if before[1] else "")
    for ex, e in collected.items():  # each example keeps the commit and date of its own analysis
        data["examples"][ex] = {**e, "analysis_commit": commit, "generated": today}
    data["examples"] = {ex: data["examples"][ex] for ex in EXAMPLES if ex in data["examples"]}
    return data


def cmd_collect(args) -> int:
    wiki = Path(args.wiki).resolve()
    names = [e for e in EXAMPLES if not args.example or e in args.example]
    problems = [m for ex in names for case, c in EXAMPLES[ex]["cases"].items() if "path" in c
                for m in [check_recorded(ex, case, c)] if m]
    problems += [f"{ex}:{case}: no fresh run at {run_dir(ex, case)}; run `build_wiki.py run` first"
                 for ex in names for case, c in EXAMPLES[ex]["cases"].items()
                 if "config" in c and not _complete(ROOT / run_dir(ex, case))]
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    from pypic.provenance import code_state

    state = code_state(ROOT)  # the code the analysis runs with, checked again before anything is written
    base = ROOT / "output" / "wiki" / "_collect"
    base.mkdir(parents=True, exist_ok=True)
    # a private working directory per collect: charts, GIFs and result JSONs wait here until everything
    # checked out, so concurrent collects never see each other's files
    work = Path(tempfile.mkdtemp(prefix="collect-", dir=base))
    try:
        return _collect(args, names, wiki, work, base, state, code_state)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _collect(args, names, wiki, work, base, state, code_state) -> int:
    def one(ex):
        cmd = [sys.executable, str(Path(__file__).resolve()), "_collect_one", ex, "--wiki", str(work),
               "--json", str(work / f"{ex}.json")]
        r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        return ex, r

    failed = False
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for ex, r in pool.map(one, names):
            print(f"{ex}: {'ok' if r.returncode == 0 else 'FAILED'}", flush=True)
            if r.returncode != 0:
                failed = True
                print(r.stderr[-4000:], file=sys.stderr)
    if failed:
        return 1
    collected = {ex: json.loads((work / f"{ex}.json").read_text()) for ex in names}
    import fcntl

    with open(base / "publish.lock", "w") as lock:  # one publish at a time: images and data.json together
        fcntl.flock(lock, fcntl.LOCK_EX)
        data = json.loads(DATA.read_text()) if DATA.exists() else {"examples": {}}
        try:
            data = merge_collected(data, collected, state, code_state(ROOT), dt.date.today().isoformat())
            publish(work / "images", wiki / "images", names,
                    json.dumps(data, indent=1, ensure_ascii=False, allow_nan=False) + "\n")
        except (RuntimeError, OSError) as e:
            print(f"{e}; the wiki images and {DATA.name} are unchanged", file=sys.stderr)
            return 1
    print(f"wrote {DATA.relative_to(ROOT)}")
    return 0


def publish(src: Path, images: Path, names: list[str], data_text: str) -> None:
    """Replace images/<example> of every collected example and data.json, all or nothing: each new set is
    first copied next to its destination (same file system), then swapped in by renames; data.json is
    written to a temporary file and renamed. On any failure the renames are undone and the old images
    stay."""
    swapped = []  # (destination, old copy or None)
    try:
        for ex in names:
            dst = images / ex
            new, old = images / f".{ex}.new", images / f".{ex}.old"
            for leftover in (new, old):
                shutil.rmtree(leftover, ignore_errors=True)
            shutil.copytree(src / ex, new)
            if dst.exists():
                dst.rename(old)
                swapped.append((dst, old))
            else:
                swapped.append((dst, None))
            new.rename(dst)
        tmp = DATA.with_name(DATA.name + ".tmp")
        tmp.write_text(data_text)
        os.replace(tmp, DATA)
    except OSError:
        for dst, old in reversed(swapped):
            if dst.exists():
                dst.rename(dst.with_name(f".{dst.name}.failed"))
            if old is not None:
                old.rename(dst)
        for ex in names:
            for junk in (images / f".{ex}.new", images / f".{ex}.failed"):
                shutil.rmtree(junk, ignore_errors=True)
        raise
    for _, old in swapped:
        if old is not None:
            shutil.rmtree(old, ignore_errors=True)


# -------------------------------------------------------------------------------------------- render

LABELS = {
    "en": {
        "case": "case", "solver": "solver", "grid": "grid (cells)", "box": "box [c/ω_pe]", "dx": "Δx [c/ω_pe]",
        "dt": "Δt [1/ω_pe]", "steps": "steps", "t_end": "t_end [1/ω_pe]", "ppc": "particles per cell",
        "particles": "macro-particles", "species": "species", "boundaries": "boundaries",
        "fields": "external / initial fields", "filter": "current filter", "gpus": "GPUs",
        "wall": "wall time", "source": "data", "fresh": "fresh run", "recorded": "recorded production run",
        "commit": "simulation commit", "run_id": "run id", "created": "started", "gpu": "GPU",
        "evidence": "recorded evidence", "dir": "output directory", "none": "none", "passes": "binomial, {n} passes",
        "param": "planner parameter", "refined": "refined run (Δt/2)",
    },
    "ja": {
        "case": "ケース", "solver": "解法", "grid": "格子数", "box": "計算領域 [c/ω_pe]", "dx": "Δx [c/ω_pe]",
        "dt": "Δt [1/ω_pe]", "steps": "ステップ数", "t_end": "終了時刻 [1/ω_pe]", "ppc": "セルあたり粒子数",
        "particles": "マクロ粒子数", "species": "粒子種", "boundaries": "境界条件",
        "fields": "外部場・初期場", "filter": "電流フィルタ", "gpus": "GPU 数",
        "wall": "計算時間", "source": "データ", "fresh": "新規実行", "recorded": "記録済み本番実行",
        "commit": "シミュレーションのコミット", "run_id": "run id", "created": "開始時刻", "gpu": "GPU",
        "evidence": "記録済みエビデンス", "dir": "出力ディレクトリ", "none": "なし", "passes": "二項フィルタ {n} パス",
        "param": "プランナーのパラメータ", "refined": "細分化実行（Δt/2）",
    },
}


def _num(v) -> str:
    if v is None:
        return "—"
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, int):
        return f"{v:,}"
    if isinstance(v, float):
        return f"{v:.4g}"
    if isinstance(v, (list, tuple)):
        return " × ".join(_num(x) for x in v)
    return str(v)


def _wall(s) -> str:
    if s is None:
        return "—"
    return f"{s:.0f} s" if s < 600 else f"{s / 60:.1f} min"


def _boundaries(bnd: list[str]) -> str:
    """["x: periodic", "y: periodic", ...] -> "x, y, z: periodic" when every axis has the same kind."""
    kinds = {b.split(": ", 1)[1] for b in bnd}
    if len(kinds) == 1:
        return ", ".join(b.split(":")[0] for b in bnd) + ": " + kinds.pop()
    return ", ".join(bnd)


def _table(header: list[str], rows: list[list[str]]) -> str:
    esc = [[str(c).replace("|", "\\|") for c in r] for r in rows]
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(r) + " |" for r in esc]
    return "\n".join(lines)


def table_conditions(ex: dict, lang: str) -> str:
    L = LABELS[lang]
    cases = ex["cases"]
    rows = []

    def row(label, fn):
        rows.append([f"**{label}**"] + [fn(c) for c in cases.values()])

    row(L["source"], lambda c: L[c["source"]])
    row(L["solver"], lambda c: c["conditions"]["solver"])
    row(L["grid"], lambda c: _num(c["conditions"]["n_cells"]))
    row(L["box"], lambda c: _num([round(x, 4) for x in c["conditions"]["L"]]))
    row(L["dx"], lambda c: _num(c["conditions"]["dx"][0]) if len(set(c["conditions"]["dx"])) == 1
        else _num(c["conditions"]["dx"]))
    row(L["dt"], lambda c: _num(c["conditions"]["dt"]))
    row(L["steps"], lambda c: _num(c["conditions"]["n_steps"]))
    row(L["t_end"], lambda c: _num(c["conditions"]["t_end"]))
    row(L["ppc"], lambda c: _num(c["conditions"]["ppc"]))
    row(L["particles"], lambda c: _num(c["conditions"]["particles"]))
    row(L["species"], lambda c: "<br>".join(c["conditions"]["species"]))
    row(L["boundaries"], lambda c: _boundaries(c["conditions"]["boundaries"]))
    row(L["fields"], lambda c: ", ".join(c["conditions"]["fields"]) or L["none"])
    row(L["filter"], lambda c: L["passes"].format(n=c["conditions"]["filter_passes"])
        if c["conditions"]["filter_passes"] else L["none"])
    row(L["gpus"], lambda c: _num(c["conditions"]["gpus"]))
    row(L["wall"], lambda c: _wall(c["conditions"]["wall_time_s"]))
    return _table([L["case"]] + [f"`{k}`" for k in cases], rows)


def table_params(ex: dict, lang: str) -> str:
    L = LABELS[lang]
    cases = ex["cases"]
    keys = []
    for c in cases.values():
        for k in c["conditions"]["params"]:
            if k not in keys:
                keys.append(k)
    rows = [[f"`{k}`"] + [_num(c["conditions"]["params"].get(k)) for c in cases.values()] for k in keys]
    return _table([L["param"]] + [f"`{k}`" for k in cases], rows)


def table_provenance(ex: dict, lang: str) -> str:
    L = LABELS[lang]
    rows = []
    for name, c in ex["cases"].items():
        p = c["provenance"]
        commit = (p["git_commit"] or "—")
        short = commit[:9] + ("-dirty" if commit.endswith("-dirty") else "")
        ev = f"`{c['evidence']}`" if c.get("evidence") else "—"
        rows.append([f"`{name}`", L[c["source"]], f"`{short}`", p.get("run_id") or "—", p.get("created") or "—",
                     f"{p.get('world_size')} × {p.get('gpu_name') or 'GPU'}", f"`{p.get('run_dir')}`", ev])
        if c.get("provenance_refined"):
            q = c["provenance_refined"]
            rows.append([f"`{name}` ({L['refined']})", L["fresh"], f"`{(q['git_commit'] or '—')[:9]}`",
                         q.get("run_id") or "—", q.get("created") or "—",
                         f"{q.get('world_size')} × {q.get('gpu_name') or 'GPU'}", f"`{q.get('run_dir')}`", "—"])
    return _table([L["case"], L["source"], L["commit"], L["run_id"], L["created"], L["gpu"], L["dir"], L["evidence"]],
                  rows)


TABLES = {"conditions": table_conditions, "params": table_params, "provenance": table_provenance}
PLACEHOLDER = re.compile(r"\{\{\s*(.*?)\s*\}\}")


def _lookup(data: dict, ex_name: str, path: str):
    if ":" in path.split(".")[0]:
        ex_name, path = path.split(":", 1)
    ex = data["examples"][ex_name]
    parts = path.split(".")
    kind, case, keys = parts[0], parts[1], parts[2:]
    section = {"r": "results", "c": "conditions", "p": "provenance"}[kind]
    v = ex["cases"][case][section]
    for k in keys:
        v = v[int(k)] if isinstance(v, list) else v[k]
    return v


def _format(v, filters: list[str]) -> str:
    for f in filters:
        if f == "abs":
            v = abs(v)
        elif f == "neg":
            v = -v
        elif f == "half":
            v = v / 2
        else:
            return "—" if v is None else format(v, f)
    if v is None:
        return "—"
    return _num(v) if not isinstance(v, str) else v


def render_text(text: str, data: dict, ex_name: str | None, lang: str) -> str:
    def sub(m):
        expr = m.group(1)
        if expr.startswith("table:"):
            return TABLES[expr[6:]](data["examples"][ex_name], lang)
        path, *filters = [s.strip() for s in expr.split("|")]
        return _format(_lookup(data, ex_name, path), filters)

    out = PLACEHOLDER.sub(sub, text)
    if "{{" in out:  # (a lone "}}" is LaTeX, e.g. v_{te}})
        raise ValueError(f"unresolved placeholder in {ex_name or 'Home'} ({lang})")
    return out


def page_name(page: str, lang: str) -> str:
    return page if lang == "en" else f"{page}-ja"


def sidebar(data: dict) -> str:
    lines = ["**English**", "", "- [Home](Home)"]
    lines += [f"- [{e['title']['en']}]({e['page']})" for e in data["examples"].values()]
    lines += ["", "**日本語**", "", "- [ホーム](Home-ja)"]
    lines += [f"- [{e['title']['ja']}]({e['page']}-ja)" for e in data["examples"].values()]
    return "\n".join(lines) + "\n"


def footer(data: dict) -> str:
    return ("Generated by [`scripts/build_wiki.py`](https://github.com/ytworks/pypic/blob/main/scripts/build_wiki.py) "
            "from the runs' own output. 自動生成ページ：数値・図・動画はすべて計算結果から生成しています。\n")


def page_stamp(ex: dict, lang: str) -> str:
    """The analysis commit and date of one example's data, appended to its pages."""
    c = ex["analysis_commit"]
    short = c[:9] + ("-dirty" if c.endswith("-dirty") else "")
    if lang == "ja":
        return (f"\n---\n*この記事の数値・図・動画はコミット `{short}` の解析コードで "
                f"{ex['generated']} に生成しました。*\n")
    return f"\n---\n*Numbers, charts and animations of this article: analysis commit `{short}`, {ex['generated']}.*\n"


def render(data: dict, out: Path) -> list[str]:
    """Write every page into out; returns the page file names."""
    out.mkdir(parents=True, exist_ok=True)
    written = []
    by_page = {e["page"]: name for name, e in data["examples"].items()}
    for lang in LANGS:
        for tpl in sorted((TEMPLATES / lang).glob("*.md")):
            page = tpl.stem
            ex_name = by_page.get(page)
            if page != "Home" and ex_name is None:
                raise ValueError(f"template {tpl} names no example page")
            text = render_text(tpl.read_text(), data, ex_name, lang)
            if ex_name:
                text += page_stamp(data["examples"][ex_name], lang)
            name = page_name(page, lang) + ".md"
            (out / name).write_text(text)
            written.append(name)
    (out / "_Sidebar.md").write_text(sidebar(data))
    (out / "_Footer.md").write_text(footer(data))
    return written + ["_Sidebar.md", "_Footer.md"]


def cmd_render(args) -> int:
    data = json.loads(Path(args.data).read_text())
    pages = render(data, Path(args.wiki))
    print(f"rendered {len(pages)} pages into {args.wiki}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run")
    p.add_argument("--devices", default="0,1,2,3")
    p.add_argument("--example", nargs="*")
    p = sub.add_parser("collect")
    p.add_argument("--wiki", required=True)
    p.add_argument("--example", nargs="*")
    p.add_argument("--jobs", type=int, default=8)
    p = sub.add_parser("render")
    p.add_argument("--wiki", required=True)
    p.add_argument("--data", default=str(DATA))
    p = sub.add_parser("_collect_one")
    p.add_argument("example")
    p.add_argument("--wiki", required=True)
    p.add_argument("--json", required=True)
    a = ap.parse_args()
    if a.cmd == "run":
        return cmd_run(a)
    if a.cmd == "collect":
        return cmd_collect(a)
    if a.cmd == "render":
        return cmd_render(a)
    res = collect_example(a.example, Path(a.wiki).resolve())
    Path(a.json).write_text(json.dumps(_clean(res), allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
