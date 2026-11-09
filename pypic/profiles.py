"""1-D shape functions for density profiles (Beam.profile) and initial field profiles
(SimConfig.initial_fields). Pure array code: works with NumPy and CuPy (``xp``).

A shape is a dict: {"shape": name, ...parameters} with positions in c/omega_pe:
    uniform                                    1
    step          x0, x1                       1 on [x0, x1), 0 elsewhere
    ramp          x0, x1                       linear 0 -> 1 on [x0, x1], 1 after
    tanh          center, width                tanh((x - center) / width)
    sech2         center, width                sech^2((x - center) / width)
    double_tanh   centers [c1, c2], width      tanh((x-c1)/w) - tanh((x-c2)/w) - 1  (periodic double Harris sheet)
    double_sech2  centers [c1, c2], width      sech^2((x-c1)/w) + sech^2((x-c2)/w)
    gauss         center, width                exp(-((x - center) / width)^2)
    cos           mode, length, phase          cos(2 pi mode x / length + phase)
    table         x [...], values [...]        linear interpolation (constant beyond the ends)
plus an optional "scale" (multiplier) and "offset" (added) for every shape.
"""
from __future__ import annotations

import math

import numpy as np

SHAPES = ("uniform", "step", "ramp", "tanh", "sech2", "double_tanh", "double_sech2", "gauss", "cos", "table")
_PARAMS = {"uniform": set(), "step": {"x0", "x1"}, "ramp": {"x0", "x1"}, "tanh": {"center", "width"},
           "sech2": {"center", "width"}, "double_tanh": {"centers", "width"}, "double_sech2": {"centers", "width"},
           "gauss": {"center", "width"}, "cos": {"mode", "length", "phase"}, "table": {"x", "values"}}


def validate(spec: dict, extra=()) -> dict:
    spec = dict(spec)
    name = spec.get("shape")
    if name not in SHAPES:
        raise ValueError(f"profile shape must be one of {SHAPES} (got {name!r})")
    allowed = _PARAMS[name] | {"shape", "scale", "offset"} | set(extra)
    unknown = set(spec) - allowed
    if unknown:
        raise KeyError(f"profile {name}: unknown key(s) {sorted(unknown)}")
    need = _PARAMS[name] - {"phase"}
    missing = need - set(spec)
    if missing:
        raise KeyError(f"profile {name}: missing {sorted(missing)}")
    return spec


def evaluate(spec: dict, x, xp=np):
    """Shape value at positions x (array)."""
    name = spec["shape"]
    x = xp.asarray(x, dtype=float)
    if name == "uniform":
        f = xp.ones_like(x)
    elif name == "step":
        f = ((x >= spec["x0"]) & (x < spec["x1"])).astype(float)
    elif name == "ramp":
        f = xp.clip((x - spec["x0"]) / (spec["x1"] - spec["x0"]), 0.0, 1.0)
    elif name == "tanh":
        f = xp.tanh((x - spec["center"]) / spec["width"])
    elif name == "sech2":
        f = 1.0 / xp.cosh((x - spec["center"]) / spec["width"]) ** 2
    elif name == "double_tanh":
        c1, c2 = spec["centers"]
        w = spec["width"]
        f = xp.tanh((x - c1) / w) - xp.tanh((x - c2) / w) - 1.0
    elif name == "double_sech2":
        c1, c2 = spec["centers"]
        w = spec["width"]
        f = 1.0 / xp.cosh((x - c1) / w) ** 2 + 1.0 / xp.cosh((x - c2) / w) ** 2
    elif name == "gauss":
        f = xp.exp(-(((x - spec["center"]) / spec["width"]) ** 2))
    elif name == "cos":
        f = xp.cos(2 * math.pi * spec["mode"] * x / spec["length"] + spec.get("phase", 0.0))
    elif name == "table":
        f = xp.interp(x, xp.asarray(spec["x"], dtype=float), xp.asarray(spec["values"], dtype=float))
    else:  # pragma: no cover - validate() rejects it
        raise ValueError(name)
    return spec.get("scale", 1.0) * f + spec.get("offset", 0.0)
