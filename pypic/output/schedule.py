"""Output cadence from physics (OUT-D4) and storage estimates. Pure Python (CPU-testable).

Every diagnostic states *what it must resolve* rather than a step count:

    every: 25                                  25 steps
    every: {t: 0.5}                            every 0.5 / omega_pe
    every: {count: 100}                        100 outputs over the run
    every: {per: growth_time, samples: 10}     10 samples per e-folding time
    every: {nyquist: shortest_period, oversample: 2}   Nyquist of that period, 2x oversampled
    every: auto                                kind-specific rule below

Named time scales (``OutputConfig.timescales``, in 1/omega_pe) are set by each experiment
planner: ``plasma_period`` (2 pi), ``shortest_period`` (2 pi / highest frequency the analysis
must resolve), ``growth_time`` (1 / growth rate), ``gyro_period``, ``t_end`` ...

``resolve`` turns every cadence into an integer step interval, estimates bytes per output,
checks aliasing / sampling, and enforces ``budget_GB`` by widening only ``auto`` cadences.
The same function runs in the planner (report) and in the engine (actual schedule).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

AUTO_RULES = {
    # kind: requirements; every rule whose timescale is set applies and the strictest (shortest
    # interval) wins; ("count", n) is the fallback when none applies
    "lineout": [("nyquist", "shortest_period", 2.0), ("per", "growth_time", 10), ("count", 200)],
    "modes": [("nyquist", "shortest_period", 2.0), ("per", "growth_time", 10), ("count", 200)],
    "slice": [("count", 100)],
    "kspectrum": [("count", 50)],
    "histogram": [("count", 100)],
    "fields": [("count", 20)],
    "particles": [("count", 5)],
    "tracks": [("per", "shortest_period", 16), ("count", 500)],
}
SPECTRAL_KINDS = ("lineout", "modes", "slice", "tracks")


@dataclass
class Resolved:
    name: str
    kind: str
    every: int  # steps
    auto: bool
    rule: str
    first: int  # first output step
    last: int  # last output step allowed
    bytes_per_output: int
    average: bool = False
    notes: list = field(default_factory=list)

    @property
    def n_outputs(self) -> int:
        if self.last < self.first:
            return 0
        lo = -(-self.first // self.every) * self.every
        return 0 if lo > self.last else (self.last - lo) // self.every + 1

    @property
    def total_bytes(self) -> int:
        return self.n_outputs * self.bytes_per_output

    def due(self, step: int) -> bool:
        return self.first <= step <= self.last and step % self.every == 0


def _component_count(names, species) -> int:
    n = 0
    for f in names:
        base = f.split(":", 1)[0]
        if base in ("E", "B", "J", "flux"):
            n += 3
        elif base in ("Ex", "Ey", "Ez", "Bx", "By", "Bz", "Jx", "Jy", "Jz", "rho", "n", "ekin"):
            n += 1
        else:
            raise ValueError(f"unknown field {f!r}")
        if ":" in f and f.split(":", 1)[1] not in species:
            raise ValueError(f"field {f!r}: unknown species")
    return n


def region_cells(cfg, region) -> list[tuple[int, int]]:
    """Global cell index range [i0, i1) per axis of a physical region (None = whole axis)."""
    out = []
    for a in range(3):
        n = cfg.n_cells[a]
        r = None if region is None else region[a]
        if r is None:
            out.append((0, n))
        else:
            i0 = max(0, int(math.floor(r[0] / cfg.dx[a] + 1e-9)))
            i1 = min(n, int(math.ceil(r[1] / cfg.dx[a] - 1e-9)))
            if i1 <= i0:
                raise ValueError(f"empty region along axis {'xyz'[a]}: {r}")
            out.append((i0, i1))
    return out


def particle_estimate(cfg, species=()) -> int:
    n = cfg.n_cells[0] * cfg.n_cells[1] * cfg.n_cells[2] * cfg.ppc
    tot = 0
    for s in cfg.species:
        if species and s.name not in species:
            continue
        tot += n * len(s.beams) + len(getattr(s, "particles", ()) or ())
    return tot


def bytes_per_output(d, cfg) -> int:
    """Uncompressed size of one output of diagnostic d (upper bound for the compressed file)."""
    item = 4 if d.dtype == "float32" else 8
    spn = [s.name for s in cfg.species]
    nc = list(cfg.n_cells)
    if d.kind in ("fields", "slice", "lineout", "modes", "kspectrum"):
        ncomp = _component_count(d.fields, spn)
    if d.kind == "fields":
        rc = region_cells(cfg, d.region)
        cells = math.prod(-(-(i1 - i0) // c) for (i0, i1), c in zip(rc, d.coarsen))
        return cells * ncomp * item
    if d.kind == "slice":
        a = "xyz".index(d.axis)
        cells = math.prod(-(-nc[b] // d.coarsen[b]) for b in range(3) if b != a)
        return cells * ncomp * item
    if d.kind == "lineout":
        return nc["xyz".index(d.axis)] * ncomp * item
    if d.kind == "modes":
        return len(d.modes) * ncomp * 2 * item
    if d.kind == "kspectrum":
        n = max(-(-nc[b] // d.coarsen[b]) for b in range(3))
        return (n // 2 + 1) * ncomp * item
    if d.kind == "histogram":
        nsp = len(d.species) if d.species else len(spn)
        return math.prod(int(b) for b in d.bins) * nsp * item
    if d.kind == "particles":
        n = particle_estimate(cfg, d.species) // d.stride
        frac = 1.0 if not d.select else 0.1  # a filter keeps a fraction; rough planning figure
        return int(n * frac) * (7 * item + 8)
    if d.kind == "tracks":
        n = particle_estimate(cfg, d.species) // d.stride
        return n * (6 + (6 if d.with_fields else 0)) * item
    raise ValueError(d.kind)


def _interval(spec, kind, ts: dict, dt: float, n_steps: int):
    """(steps, rule text, average default) for one cadence spec."""
    if isinstance(spec, bool):
        raise ValueError("every must be a number of steps, a dict or 'auto'")
    if isinstance(spec, int):
        if spec < 1:
            raise ValueError("every must be >= 1 step")
        return spec, f"{spec} steps", False
    if spec == "auto":
        cands = []
        for rule in AUTO_RULES[kind]:
            if rule[0] == "count":
                if not cands:
                    return max(1, n_steps // rule[1]), f"auto: {rule[1]} outputs", False
            elif rule[1] in ts:
                sub = {"nyquist": rule[1], "oversample": rule[2]} if rule[0] == "nyquist" else \
                    {"per": rule[1], "samples": rule[2]}
                st, txt, avg = _interval(sub, kind, ts, dt, n_steps)
                cands.append((st, txt, avg))
        st = min(c[0] for c in cands)
        txt = " & ".join(c[1] for c in cands)
        # time averaging (anti-aliasing) whenever a Nyquist requirement is part of the rule
        return st, "auto: " + txt, any(c[2] for c in cands)
    if not isinstance(spec, dict) or len(spec) not in (1, 2):
        raise ValueError(f"invalid cadence {spec!r}")
    if "t" in spec:
        return max(1, int(math.floor(float(spec["t"]) / dt + 1e-9))), f"every {spec['t']:g}/omega_pe", False
    if "count" in spec:
        return max(1, n_steps // int(spec["count"])), f"{spec['count']} outputs", False
    if "per" in spec or "nyquist" in spec:
        key = spec.get("per", spec.get("nyquist"))
        if key not in ts:
            raise ValueError(f"cadence {spec!r}: timescale {key!r} not set (have {sorted(ts)})")
        tau = float(ts[key])
        if "per" in spec:
            n = float(spec.get("samples", 10))
            return max(1, int(math.floor(tau / n / dt + 1e-9))), f"{n:g} per {key}", False
        s = float(spec.get("oversample", 2.0))
        return max(1, int(math.floor(tau / (2.0 * s) / dt + 1e-9))), f"Nyquist of {key} x{s:g}", True
    raise ValueError(f"invalid cadence {spec!r}")


def resolve(cfg) -> list[Resolved]:
    """Resolve every diagnostic of cfg.output into integer step intervals (and enforce the budget)."""
    out_cfg = cfg.output
    ts = dict(out_cfg.timescales)
    dt, n_steps = cfg.dt, cfg.n_steps
    res = []
    for d in out_cfg.diagnostics:
        every, rule, avg_default = _interval(d.every, d.kind, ts, dt, n_steps)
        first = int(math.ceil(d.t_start / dt - 1e-9))
        last = n_steps if d.t_end is None else min(n_steps, int(math.floor(d.t_end / dt + 1e-9)))
        average = bool(d.average) if d.average is not None else (avg_default and d.kind in ("lineout", "modes"))
        r = Resolved(d.name, d.kind, every, d.every == "auto", rule, first, last, bytes_per_output(d, cfg), average)
        tdt = every * dt
        if d.kind in SPECTRAL_KINDS and "shortest_period" in ts and not average and every > 1:
            if tdt > ts["shortest_period"] / 2:
                r.notes.append(f"WARNING aliasing: interval {tdt:.3g} > half the shortest period "
                               f"{ts['shortest_period'] / 2:.3g} (set average: true or a shorter cadence)")
        if d.kind in ("lineout", "modes") and "growth_time" in ts and tdt > ts["growth_time"] / 5:
            r.notes.append(f"WARNING sampling: {ts['growth_time'] / tdt:.1f} samples per growth time (< 5)")
        res.append(r)
    budget = out_cfg.budget_GB
    if budget is not None:
        lim = budget * 1e9
        total = lambda: sum(x.total_bytes for x in res) + scalars_bytes(cfg)  # noqa: E731
        while total() > lim:
            cand = [x for x in res if x.auto and x.n_outputs > 1]
            if not cand:
                raise ValueError(f"output exceeds budget_GB={budget}: {total() / 1e9:.3g} GB with no 'auto' "
                                 "cadence left to widen\n" + report(cfg, res))
            x = max(cand, key=lambda r: r.total_bytes)
            x.every *= 2
            if "widened" not in " ".join(x.notes):
                x.notes.append("widened to fit budget_GB")
    return res


def scalars_bytes(cfg) -> int:
    if not cfg.output.scalars or cfg.diag_interval <= 0:
        return 0
    return (cfg.n_steps // cfg.diag_interval + 1) * (20 + len(cfg.species)) * 8


def _fmt_bytes(b: float) -> str:
    for unit in ("B", "kB", "MB", "GB", "TB"):
        if b < 1000 or unit == "TB":
            return f"{b:.3g} {unit}"
        b /= 1000.0
    return ""


def report(cfg, res: list[Resolved] | None = None) -> str:
    res = resolve(cfg) if res is None else res
    ts = cfg.output.timescales
    lines = ["output plan (uncompressed sizes; zstd typically saves 10-50 % more):",
             f"  {'diagnostic':16s} {'kind':10s} {'every':>7s} {'dt_out':>9s} {'per period':>10s} "
             f"{'outputs':>7s} {'each':>9s} {'total':>9s}  rule"]
    tot = scalars_bytes(cfg)
    ref = ts.get("shortest_period", ts.get("plasma_period"))
    for r in res:
        tdt = r.every * cfg.dt
        per = f"{ref / tdt:.1f}" if ref else "-"
        lines.append(f"  {r.name:16s} {r.kind:10s} {r.every:7d} {tdt:9.4g} {per:>10s} {r.n_outputs:7d} "
                     f"{_fmt_bytes(r.bytes_per_output):>9s} {_fmt_bytes(r.total_bytes):>9s}  "
                     f"{r.rule}{' (time-averaged)' if r.average else ''}")
        tot += r.total_bytes
        for n in r.notes:
            lines.append(f"      {n}")
    lines.append(f"  scalars every {cfg.diag_interval} steps; total ~{_fmt_bytes(tot)}"
                 + (f" (budget {cfg.output.budget_GB:g} GB)" if cfg.output.budget_GB else ""))
    return "\n".join(lines)
