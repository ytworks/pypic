"""Example analysis helpers that the wiki charts rely on (examples/<name>/analyze.py), on synthetic fields."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("scipy")
ROOT = Path(__file__).resolve().parents[2]


def _reconnection():
    spec = importlib.util.spec_from_file_location("reconnection_analyze", ROOT / "examples/reconnection/analyze.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _grid():
    x = np.linspace(0, 2 * np.pi, 256, endpoint=False)
    y = np.arange(-20.0, 20.0, 0.25)
    return x, y, np.meshgrid(x, y, indexing="ij")


def test_reconnected_flux_is_taken_between_the_x_and_o_points_of_a_bent_sheet():
    R = _reconnection()
    x, y, (X, Y) = _grid()
    A = -(Y - np.cos(X)) ** 2 + 0.1 * np.cos(X)  # null line y = cos x; A on it = 0.1 cos x
    xs, os_ = R.null_points(A, x, y, 0.0, 10.0, 0.01)
    assert len(xs) == 1 and len(os_) == 1
    assert R.xo_flux(A, x, y, 0.0, 10.0, 0.01) == pytest.approx(0.2, abs=0.01)  # grid: dy = 0.25


def test_no_null_points_on_the_edge_of_the_search_band():
    R = _reconnection()
    x, y, (X, Y) = _grid()
    A = 0.5 * (Y - 12.0) ** 2 + np.cos(X)  # the null line y = 12 lies outside |y| < 10
    xs, os_ = R.null_points(A, x, y, 0.0, 10.0, 0.01)
    assert xs == [] and os_ == []


def test_no_flux_is_reported_without_an_x_and_an_o_point():
    R = _reconnection()
    x, y, (X, Y) = _grid()
    A = np.cos(np.pi * (Y - 12 * np.cos(X)) / 40) + 0.1 * np.cos(X)  # the null line leaves |y| < 10 in places
    assert R.null_points(A, x, y, 0.0, 10.0, 0.01) == ([], [])
    assert np.isnan(R.xo_flux(A, x, y, 0.0, 10.0, 0.01))


def test_the_flux_panel_draws_frames_without_a_measurable_flux():
    from types import SimpleNamespace

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    R = _reconnection()
    x, y, (X, Y) = _grid()
    frames = [0.01 * Y ** 2 + np.cos(X), 0.01 * (Y - 12) ** 2 + np.cos(X)]  # the second has no X/O point
    run = SimpleNamespace(meta={"planner": {"y_sheets": [0.0, 40.0], "delta": 1.0, "B0": 1.0, "Lx": 2 * np.pi}})
    panel = R.FluxPanel((run, SimpleNamespace(t=np.array([0.0, 1.0])), x, y, frames))
    assert np.isfinite(panel.psi[0]) and np.isnan(panel.psi[1])
    fig, ax = plt.subplots()
    panel.init(ax)
    panel.update(1.0)
    plt.close(fig)


def _bow_shock():
    """examples/bow_shock/analyze.py; its ``import theory`` must not meet the theory module of another example."""
    import sys

    folder = ROOT / "examples/bow_shock"
    spec = importlib.util.spec_from_file_location("bow_shock_analyze", folder / "analyze.py")
    mod = importlib.util.module_from_spec(spec)
    other = sys.modules.pop("theory", None)
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.modules.pop("theory", None)
        if other is not None:
            sys.modules["theory"] = other
        if str(folder) in sys.path:
            sys.path.remove(str(folder))
    return mod


def test_bow_shock_upstream_deviation_is_nan_without_a_measurement():
    from types import SimpleNamespace

    B = _bow_shock()
    v, n_cells = 0.25, (48, 64, 64)  # R_mp = 20 cells: the upstream region is cells 4 to 16 along x
    info = {"r_mp": 10.0, "v_sw": v, "crossing_steps": 100}

    def state(block, density):
        """Upstream state of a uniform wind with the given density at each output step."""
        shape = tuple(-(-n // block) for n in n_cells)
        data = {(it, name): np.full(shape, d * (v if name.startswith("flux") else 1.0))
                for it, d in density.items() for name in ("n:e", "n:i", "fluxx:e", "fluxx:i")}
        f = SimpleNamespace(iterations=sorted(density), get=lambda it, name: SimpleNamespace(data=data[it, name]))
        run = SimpleNamespace(fields=lambda name: f, n_cells=n_cells, dx=(0.5, 0.5, 0.5),
                              config={"output": {"diagnostics": [{"name": "coarse", "coarsen": [block] * 3}]}})
        return B.upstream_state(run, info)

    up = state(4, {0: 1.0, 100: 1.0, 200: 0.9})
    assert up["blocks"] == 3 * 12 * 12 and up["outputs"] == 2
    assert up["deviation"] == pytest.approx(0.1) and up["e"]["u_x"] == pytest.approx(0.0, abs=1e-12)
    # an output without a finite density in the last crossing: not "no deviation"
    assert np.isnan(state(4, {0: 1.0, 100: float("nan"), 200: 1.0})["deviation"])
    # blocks of 12 cells: none lies inside the region, so nothing was measured
    up = state(12, {0: 1.0, 100: 1.0, 200: 1.0})
    assert up["blocks"] == 0 and np.isnan(up["deviation"])
