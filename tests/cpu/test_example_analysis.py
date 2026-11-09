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
    assert R.xo_flux(A, y, 0.0, 10.0) == pytest.approx(0.2, abs=0.01)  # grid: dy = 0.25


def test_no_null_points_on_the_edge_of_the_search_band():
    R = _reconnection()
    x, y, (X, Y) = _grid()
    A = 0.5 * (Y - 12.0) ** 2 + np.cos(X)  # the null line y = 12 lies outside |y| < 10
    xs, os_ = R.null_points(A, x, y, 0.0, 10.0, 0.01)
    assert xs == [] and os_ == []
