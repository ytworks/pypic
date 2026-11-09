"""Common command line of the examples: YAML parameters -> planner -> run -> output directory.

    python examples/<name>/run.py --config examples/<name>/configs/<cfg>.yaml
    python -m pypic.launch -n 8 examples/<name>/run.py --config ...
    options: --set key=value (parameters) / --set overrides.key=value (any SimConfig field)
             --set output.key=value (OutputConfig field), --out DIR, --overwrite, --plan-only, --progress N

A planner is ``plan(params, n_gpus) -> (SimConfig, report: str, info: dict)``; ``info`` (theory
predictions and chosen numbers) goes to run.json -> meta.planner for the analysis.
"""
from __future__ import annotations

import argparse
import dataclasses
import os
import sys

import yaml

from ..config import OutputConfig


def parse_value(v: str):
    try:
        return yaml.safe_load(v)
    except yaml.YAMLError:
        return v


def load_params(path: str, sets: list[str]) -> dict:
    with open(path) as f:
        d = yaml.safe_load(f) or {}
    for kv in sets:
        k, v = kv.split("=", 1)
        if k.startswith("overrides.") or k.startswith("output."):
            group, key = k.split(".", 1)
            d.setdefault(group, {})[key] = parse_value(v)
        else:
            d[k] = parse_value(v)
    return d


def run_example(example: str, plan, params_from_dict, argv=None, doc: str | None = None) -> int:
    ap = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--set", action="append", default=[],
                    help="key=value (parameters), overrides.key=value (SimConfig), output.key=value")
    ap.add_argument("--out", default=None, help="output directory (default output/<example>/<config name>)")
    ap.add_argument("--overwrite", action="store_true",
                    help="replace the run already in the output directory (its output files are removed)")
    ap.add_argument("--plan-only", action="store_true")
    ap.add_argument("--progress", type=int, default=100)
    a = ap.parse_args(argv)

    d = load_params(a.config, a.set)
    world = int(os.environ.get("PYPIC_WORLD_SIZE", "1"))
    rank = int(os.environ.get("PYPIC_RANK", "0"))
    out_over = d.pop("output", {})
    params = params_from_dict(d)
    cfg, report, info = plan(params, n_gpus=world)
    name = os.path.splitext(os.path.basename(a.config))[0]
    # output.* (YAML section or --set): merged into the planner's output and validated like any
    # OutputConfig (diagnostics given as dicts become Diagnostic objects)
    out = {f.name: getattr(cfg.output, f.name) for f in dataclasses.fields(cfg.output)}
    out.update(out_over)
    out["dir"] = a.out or out_over.get("dir") or os.path.join("output", example, name)
    if a.overwrite:
        out["overwrite"] = True
    cfg.output = OutputConfig.from_dict(out)
    from ..output.schedule import report as output_report

    full_report = report + "\n" + output_report(cfg)
    if rank == 0:
        print(full_report, flush=True)
    if a.plan_only:
        return 0

    from ..simulation import Simulation

    sim = Simulation(cfg)
    sim.output.update_meta({"example": example, "config_file": a.config, "params": d, "planner": info,
                            "report": full_report})
    sim.run(progress=a.progress)
    stats = {"wall_time_s": sim.wall_time, "ns_per_particle_step": sim.push_rate(), "gpus": world,
             "steps": sim.step_count}
    sim.output.update_meta({"run": stats})
    sim.close()
    if rank == 0:
        print(f"output in {cfg.output.dir}; wall {sim.wall_time:.1f} s, "
              f"{sim.push_rate():.2f} ns/particle/step (x GPUs)", flush=True)
    return 0


def main_for(example: str, plan, params_from_dict, doc=None):
    """Entry point used by examples/<name>/run.py."""
    sys.exit(run_example(example, plan, params_from_dict, doc=doc))
