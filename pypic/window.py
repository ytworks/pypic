"""Moving window (PIC-D22): the box follows a driver along +x at the speed of light.

Positions stay those of the laboratory. After a step that ends at the time t the window has moved
by ``s = max(0, floor((t - start) / dx))`` cells, at most one more than after the step before
(dt < dx): the fields move one plane down, the particles below the new x- face are dropped and the
laboratory plane that enters behind the vacuum margin at x+ is loaded. Everything here is a
function of the configuration and of the number of steps, so every rank gets the same values.

Pure Python: the schedule and the checks of the configuration.
"""
from __future__ import annotations

import math

from . import antenna
from .profiles import support, support_interval

# the plasma ends this many cells (plus the reach of the current filter) before the x+ wall at
# least: the shape of a particle at the front and the current of one that streams ahead of it stay
# off the wall node and off the edges that the moving wall sets
MIN_MARGIN = 4
# a laser must have fallen to exp(-ENVELOPE_END^2) of its peak when its antenna leaves the window
ENVELOPE_END = 5.0


def margin_cells(cfg) -> int:
    """Cells of vacuum between the plasma and the x+ wall."""
    return round(cfg.window.margin / cfg.dx[0])


def _raw(cfg, steps: int) -> int:
    return math.floor((steps * cfg.dt - cfg.window.start) / cfg.dx[0])


def shift_after(cfg, steps: int) -> int:
    """Cells that the window has moved when ``steps`` steps are done."""
    return max(0, _raw(cfg, steps))


def steps_since_shift(cfg, steps: int) -> int:
    """For a number of steps at which the shift has just become s >= 1: the steps done since the
    shift became s - 1 (for s = 1: since the time ``start``). The plane that is loaded now was the
    front for this long."""
    s = _raw(cfg, steps)
    j = steps
    while _raw(cfg, j - 1) >= s - 1:
        j -= 1
    return steps - j


def front_cell(cfg) -> int:
    """The first cell along x that the first load leaves empty; the shift to s loads the
    laboratory plane ``front_cell + s - 1``."""
    return cfg.n_cells[0] - margin_cells(cfg)


def beam_cells(cfg, beam) -> tuple[int, int]:
    """Laboratory cells ``(lo, hi)`` along x that ``beam`` loads during the run: those that its
    profile reaches among the cells below ``front_cell + s_end``, with ``s_end`` the shift after
    ``n_steps`` steps (no profile along x: all of them). The first load takes the cells below
    ``front_cell``, every later plane is loaded when it enters."""
    dx = cfg.dx[0]
    end = front_cell(cfg) + shift_after(cfg, cfg.n_steps)
    length = [end * dx] + [cfg.n_cells[a] * cfg.dx[a] for a in (1, 2)]
    lo, hi = 0.0, length[0]
    p = beam.profile
    if isinstance(p, dict) and p.get("axis", "x") == "r":
        rs = support(p, math.sqrt(sum(max(abs(o), abs(length[a] - o)) ** 2 for a, o in enumerate(beam.origin))))
        lo, hi = beam.origin[0] - rs, beam.origin[0] + rs
    else:
        for q in p if isinstance(p, list) else [p] if p else []:
            if q.get("axis", "x") == "x":
                lo, hi = support_interval(q, length[0])
    lo = min(max(math.floor(lo / dx), 0), end)
    return lo, max(min(math.ceil(hi / dx), end), lo)


def antenna_off_shift(cfg, laser, retreat: int = 0) -> int:
    """The shift at which the antenna of ``laser`` is switched off: its cell has reached the nodes
    next to x- where Gauss's law is not checked."""
    return antenna.plane_cell(laser, cfg.dx) - antenna.first_cell(cfg, retreat) + 1


def check_config(cfg, retreat: int = 0, reach: int = 0) -> bool:
    """True if ``cfg`` has a moving window; raises on what is not supported. ``retreat``: nodes
    next to an open face that are left out of the Gauss check; ``reach``: the reach of the current
    filter."""
    w = cfg.window
    if w is None:
        return False
    if cfg.solver != "explicit":
        raise ValueError("the moving window is for the explicit solver only (PIC-D22)")
    if tuple(cfg.boundaries[0]) != ("open", "reflecting"):
        raise ValueError("the moving window needs boundaries (open, reflecting) along x (PIC-D22)")
    refused = {"injection": cfg.injection.active, "a body": cfg.body is not None,
               "external_fields": cfg.external_fields.active, "B0 / E0": any(cfg.B0) or any(cfg.E0),
               "initial_fields": bool(cfg.initial_fields),
               "the compact storage (precision: compact, PIC-D19)": cfg.precision == "compact",
               "a perturbation of a beam": any(b.perturbation for sp in cfg.species for b in sp.beams)}
    for name, on in refused.items():
        if on:
            raise ValueError(f"the moving window does not support {name} (PIC-D22)")
    n, dx = cfg.n_cells[0], cfg.dx
    m = margin_cells(cfg)
    if m < MIN_MARGIN + reach:
        raise ValueError(f"window.margin must be at least {MIN_MARGIN + reach} cells "
                         "(4 plus the reach of the current filter)")
    if n - m <= antenna.first_cell(cfg, retreat):
        raise ValueError("window.margin leaves no cell for the plasma")
    # a beam that is loaded plane by plane is loaded as a plasma at rest along x that is the same
    # on both sides of the plane; whatever else lies in the cells of the first load
    front = front_cell(cfg)
    for sp in cfg.species:
        for b in sp.beams:
            if beam_cells(cfg, b)[1] <= front:
                continue
            radial = b.drift_radial != 0.0 or (isinstance(b.profile, dict) and b.profile.get("axis", "x") == "r")
            if b.drift[0] != 0.0 or radial:
                what = "a radial profile or drift" if radial else "a drift along x"
                raise ValueError(f"species {sp.name}: a beam with {what} must end within the first load of a moving "
                                 f"window (x < {front * dx[0]:g}): the planes that enter later are loaded at rest "
                                 "along x (PIC-D22)")
    for laser in cfg.lasers:
        env = laser.envelope
        if env["shape"] != "gaussian":
            raise ValueError("a laser in a moving window needs a gaussian envelope: its antenna leaves the box")
        # the pulse has passed everywhere on the plane when the antenna is switched off; a front
        # that is curved away from the focus (advance < 0) is late by |c| r^2 at the radius r
        c = antenna.front_advance(laser, dx, cfg.dt)
        late = 0.0
        if c < 0.0:
            center = laser.center if laser.center is not None else tuple(0.5 * cfg.n_cells[a] * dx[a] for a in (1, 2))
            r2 = sum(max(center[a - 1], cfg.n_cells[a] * dx[a] - center[a - 1]) ** 2 for a in (1, 2))
            late = -c * r2
        t_off = w.start + antenna_off_shift(cfg, laser, retreat) * dx[0]
        if t_off - late < env["t_peak"] + ENVELOPE_END * env["duration"]:
            raise ValueError("the moving window takes the antenna out of the box before the laser pulse has ended: "
                             f"its cell leaves at t = {t_off:g}, the envelope needs t_peak + {ENVELOPE_END:g} duration"
                             f"{' and the delay of the curved front' if late else ''} = "
                             f"{env['t_peak'] + ENVELOPE_END * env['duration'] + late:g} (raise window.start)")
    return True
