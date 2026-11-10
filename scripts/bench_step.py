"""Step-time benchmark of an example configuration (diagnostic; never a test verdict).

    python scripts/bench_step.py --example two_stream --config examples/two_stream/configs/ei.yaml \
        --set cells_per_wavelength=64 --set t_end=4 --set overrides.diag_interval=10
    python -m pypic.launch -n 8 scripts/bench_step.py --example ... [--phases] [--with-output DIR]

The configuration comes from the example's planner, as in ``examples/<name>/run.py``; file output
is off unless ``--with-output`` is given. After ``--warmup`` steps (kernel compilation, first
sort) ``--steps`` steps are timed with ``Simulation.run``; rank 0 prints one line ``BENCH {json}``
with ms/step and GPU-ns per particle per step.

``--phases`` brackets the parts of a step with device synchronisations and reports the time of
each (mean and maximum over the ranks, ms per step). The synchronisations remove the overlap of
host and device work, so the step time of such a run is not the step time of the code: take the
total from a run without ``--phases``. A phase called inside another one is counted in the outer
phase; ``other`` is the rest of the step.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import runpy
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# phases bracketed by --phases: attribute of the object -> phase name
GRID_PHASES = {"zero_J": "zero_J", "reduce_ghosts": "reduce_ghosts", "fill_ghosts": "fill_ghosts",
               "apply_bc_J": "apply_bc_J", "push_B": "field_push", "push_E": "field_push",
               "field_energy": "field_energy", "poynting_out": "poynting_out"}
SIM_PHASES = {"_filter_J": "filter_J", "_migrate": "migrate", "_record": "record", "_cg": "cg",
              "_half_fields": "half_fields", "_gather_E": "gather_E", "_kinetic_now": "kinetic_now",
              "_remove_mean_current": "remove_mean_current"}


def load_planner(example: str):
    """(plan, params_from_dict) of examples/<example>/run.py (its imports, not its main)."""
    ns = runpy.run_path(str(REPO / "examples" / example / "run.py"), run_name="bench_step")
    plans = [v for k, v in ns.items() if k.startswith("plan_") and callable(v)]
    if len(plans) != 1 or "params_from_dict" not in ns:
        raise SystemExit(f"examples/{example}/run.py does not import one plan_* and params_from_dict")
    return plans[0], ns["params_from_dict"]


class Phases:
    """Wall time of the bracketed calls of one rank, per phase."""

    def __init__(self, sync):
        self.sync = sync
        self.depth = 0
        self.on = False
        self.time: dict[str, float] = {}
        self.calls: dict[str, int] = {}
        self.longest: dict[str, float] = {}
        self.series: dict[str, list] = {}
        self.step_time = 0.0

    def name(self, name: str):
        self.time.setdefault(name, 0.0)
        self.calls.setdefault(name, 0)
        self.longest.setdefault(name, 0.0)
        self.series.setdefault(name, [])

    def wrap(self, name: str, fn):
        self.name(name)

        def timed(*a, **k):
            if self.depth or not self.on:
                return fn(*a, **k)
            self.sync()
            self.depth += 1
            t0 = time.perf_counter()
            try:
                return fn(*a, **k)
            finally:
                self.sync()
                dt = time.perf_counter() - t0
                self.depth -= 1
                self.time[name] += dt
                self.calls[name] += 1
                self.longest[name] = max(self.longest[name], dt)
                self.series[name][-1] += dt

        return timed

    def wrap_step(self, step):
        def timed():
            if not self.on:
                return step()
            for s in self.series.values():
                s.append(0.0)
            self.sync()
            t0 = time.perf_counter()
            step()
            self.sync()
            self.step_time += time.perf_counter() - t0

        return timed

    def install(self, sim):
        import pypic.simulation as simulation

        for obj, table in ((sim.grid, GRID_PHASES), (sim, SIM_PHASES)):
            for attr, name in table.items():
                setattr(obj, attr, self.wrap(name, getattr(obj, attr)))
        sim.output.on_step = self.wrap("output", sim.output.on_step)
        for sp in sim.species:
            sp.sort = self.wrap("sort", sp.sort)
        # kernels launched by the step itself (the particle push), by kernel name
        kernels = {id(v): k[2:] for k, v in vars(sim).items() if k.startswith("k_")}
        timed = {name: self.wrap(f"kernel_{name}", simulation.launch) for name in kernels.values()}

        def launch(kernel, *a, **k):
            name = kernels.get(id(kernel))
            return timed[name](kernel, *a, **k) if name else simulation_launch(kernel, *a, **k)

        simulation_launch, simulation.launch = simulation.launch, launch
        sim.step = self.wrap_step(sim.step)

    def report(self, comm, steps: int, series: bool) -> dict:
        """ms per step of every phase, mean and maximum over the ranks (collective)."""
        import cupy as cp

        names = sorted(self.time)
        local = cp.asarray([self.time[n] for n in names] + [self.step_time])
        mean = cp.asnumpy(comm.allreduce(local.copy())) / comm.size
        peak = cp.asnumpy(comm.allreduce(local.copy(), "max"))
        ms = 1e3 / steps
        out = {"step_ms_synchronised": {"mean": mean[-1] * ms, "max": peak[-1] * ms}, "phases_ms": {}}
        for i, n in enumerate(names):
            if self.calls[n]:
                out["phases_ms"][n] = {"mean": mean[i] * ms, "max": peak[i] * ms, "calls": self.calls[n],
                                       "longest_call_rank0": self.longest[n] * 1e3}
        out["phases_ms"]["other"] = {"mean": (mean[-1] - mean[:-1].sum()) * ms}
        if series:
            out["series_ms_rank0"] = {n: [v * 1e3 for v in self.series[n]] for n in names if self.calls[n]}
        return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--example", required=True, help="directory name under examples/")
    ap.add_argument("--config", required=True)
    ap.add_argument("--set", action="append", default=[],
                    help="key=value (parameters), overrides.key=value (SimConfig), output.key=value")
    ap.add_argument("--steps", type=int, default=50, help="timed steps")
    ap.add_argument("--warmup", type=int, default=10, help="steps before the timed ones")
    ap.add_argument("--phases", action="store_true", help="time the parts of a step (synchronised)")
    ap.add_argument("--series", action="store_true", help="with --phases: per-step times of rank 0")
    ap.add_argument("--with-output", metavar="DIR", default=None,
                    help="write the example's output to DIR (replacing a run already there)")
    a = ap.parse_args(argv)

    from pypic.config import OutputConfig
    from pypic.setup.runner import load_params

    plan, params_from_dict = load_planner(a.example)
    d = load_params(a.config, a.set)
    out_over = d.pop("output", {})
    world = int(os.environ.get("PYPIC_WORLD_SIZE", "1"))
    cfg, _, _ = plan(params_from_dict(d), n_gpus=world)
    if a.with_output:
        out = {f.name: getattr(cfg.output, f.name) for f in dataclasses.fields(cfg.output)}
        out.update(out_over)
        out.update(dir=a.with_output, overwrite=True)
        cfg.output = OutputConfig.from_dict(out)
    else:
        cfg.output = OutputConfig()

    import cupy as cp

    from pypic.simulation import Simulation

    sim = Simulation(cfg, verbose=False)
    phases = Phases(cp.cuda.Device().synchronize) if a.phases else None
    if phases:
        phases.install(sim)
    sim.run(a.warmup)
    if phases:
        phases.on = True
    sim.run(a.steps)
    result = {
        "example": a.example, "config": a.config, "set": a.set, "gpus": sim.comm.size,
        "solver": cfg.solver, "precision": cfg.precision, "n_cells": list(cfg.n_cells),
        "particles": sim._n_total_end, "diag_interval": cfg.diag_interval,
        "sort_interval": cfg.sort_interval, "output": bool(a.with_output), "warmup": a.warmup,
        "steps": a.steps, "synchronised": bool(phases),
        "ms_per_step": sim.wall_time / a.steps * 1e3, "ns_per_particle_step": sim.push_rate(),
    }
    if phases:
        result.update(phases.report(sim.comm, a.steps, a.series))
    sim.close()
    if sim.rank == 0:
        print("BENCH " + json.dumps(result), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
