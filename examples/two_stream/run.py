"""3-D two-stream instability experiment.

Single GPU:   python examples/two_stream/run.py --config examples/two_stream/configs/ei.yaml
Multi GPU:    python -m pypic.launch -n 8 examples/two_stream/run.py --config .../ep.yaml
Plan only:    python examples/two_stream/run.py --config ... --plan-only
Overrides:    --set u_drift=0.3 --set solver=implicit --set overrides.ppc=32
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import yaml


def parse_value(v: str):
    try:
        return yaml.safe_load(v)
    except yaml.YAMLError:
        return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--set", action="append", default=[], help="key=value (params) or overrides.key=value")
    ap.add_argument("--out", default=None, help="output directory (default output/<config name>)")
    ap.add_argument("--plan-only", action="store_true")
    ap.add_argument("--progress", type=int, default=100)
    a = ap.parse_args()

    with open(a.config) as f:
        d = yaml.safe_load(f) or {}
    for kv in a.set:
        k, v = kv.split("=", 1)
        if k.startswith("overrides."):
            d.setdefault("overrides", {})[k.split(".", 1)[1]] = parse_value(v)
        else:
            d[k] = parse_value(v)

    from pypic.setup.planner import params_from_dict, plan_two_stream

    world = int(os.environ.get("PYPIC_WORLD_SIZE", "1"))
    rank = int(os.environ.get("PYPIC_RANK", "0"))
    params = params_from_dict(d)
    cfg, report, info = plan_two_stream(params, n_gpus=world)
    name = os.path.splitext(os.path.basename(a.config))[0]
    cfg.output_dir = a.out or os.path.join("output", name)
    if rank == 0:
        print(report, flush=True)
    if a.plan_only:
        return

    from pypic.simulation import Simulation

    sim = Simulation(cfg)
    sim.run(progress=a.progress)
    if rank == 0:
        path = sim.save()
        meta = {"params": d, "linear": {k: v for k, v in info.items() if k != "dims"},
                "wall_time_s": sim.wall_time, "ns_per_particle_step": sim.push_rate(), "gpus": world}
        with open(os.path.join(cfg.output_dir, "meta.json"), "w") as f:
            json.dump(meta, f, indent=2, default=str)
        print(f"saved {path}; wall {sim.wall_time:.1f} s, {sim.push_rate():.2f} ns/particle/step (x GPUs)")


if __name__ == "__main__":
    sys.exit(main())
