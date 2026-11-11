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


def _example(name):
    """examples/<name>/analyze.py; its ``import theory`` must not meet the theory module of another example."""
    import sys

    folder = ROOT / "examples" / name
    spec = importlib.util.spec_from_file_location(f"{name}_analyze", folder / "analyze.py")
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

    B = _example("bow_shock")
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


def test_bow_shock_upstream_region_follows_the_side_chosen_for_an_oblique_imf():
    from types import SimpleNamespace

    B = _example("bow_shock")
    v, n_cells, block = 0.25, (48, 64, 64), 4  # box 24 x 32 x 32, centre of z at 16; blocks of 2
    shape = tuple(n // block for n in n_cells)
    z = (np.arange(shape[2]) + 0.5) * 2.0
    wind = np.ones(shape)
    disturbed = wind * np.where(z < 16.0, 1.3, 1.0)  # 30 % more density on the side z < z_c only
    data = {}
    for it, d in ((0, wind), (100, disturbed), (200, disturbed)):
        for sp in "ei":
            data[it, f"n:{sp}"], data[it, f"fluxx:{sp}"] = d, d * v
    f = SimpleNamespace(iterations=[0, 100, 200], get=lambda it, name: SimpleNamespace(data=data[it, name]))
    run = SimpleNamespace(fields=lambda name: f, n_cells=n_cells, dx=(0.5, 0.5, 0.5),
                          config={"output": {"diagnostics": [{"name": "coarse", "coarsen": [block] * 3}]}})
    info = {"r_mp": 10.0, "v_sw": v, "crossing_steps": 100, "center": [12.0, 16.0, 16.0]}
    assert B.upstream_state(run, info)["deviation"] == pytest.approx(0.15)  # half of the whole region
    up = B.upstream_state(run, {**info, "upstream_z_range": (21.0, 32.0)})  # z - z_c >= 0.5 R_mp
    assert up["blocks"] == 3 * 12 * 3 and up["deviation"] == pytest.approx(0.0, abs=1e-12)
    assert up["foreshock"]["blocks"] == up["blocks"] and up["foreshock"]["deviation"] == pytest.approx(0.3)


def test_bow_shock_field_reversal_is_counted_between_the_shock_and_the_body():
    B = _example("bow_shock")
    x = np.linspace(0.0, 50.0, 501)
    bz = np.interp(x, [0.0, 29.0, 31.0, 50.0], [-1.0, -3.0, 5.0, 8.0])  # crosses zero at 29.75
    one = B.reversal(x, bz, 20.0, 45.0)
    assert one["count"] == 1 and one["x"] == pytest.approx(29.75, abs=1e-9)
    # a second and a third change inside the interval: not a single reversal, no position
    three = B.reversal(x, bz * np.where((x > 35.0) & (x < 40.0), -1.0, 1.0), 20.0, 45.0)
    assert three["count"] == 3 and np.isnan(three["x"])
    # outside the interval they do not count; without a shock (nan) nothing lies between
    assert B.reversal(x, bz, 20.0, 34.0)["count"] == 1 and B.reversal(x, bz, 31.0, 45.0)["count"] == 0
    assert B.reversal(x, bz, float("nan"), 45.0)["count"] == 0


def test_bow_shock_null_history_measures_density_velocity_and_temperature_per_region():
    from types import SimpleNamespace

    B = _example("bow_shock")
    v, t0, m, shape = 0.25, 0.01, 16.0, (8, 6, 6)

    def moments(density, ux, uz, temp):
        """Moments of a wind with the given density, velocity (ux, 0, uz) and temperature per block."""
        gu = 1.0 / np.sqrt(1.0 - ux * ux - uz * uz)
        return {"n": density, "fluxx": density * ux, "fluxy": 0.0 * density, "fluxz": density * uz,
                "ekin": density * (m * (gu - 1.0) + 1.5 * temp)}

    one = np.ones(shape)
    late = np.ones(shape)
    late[:2] = 1.1  # the slab at x-: 10 % more density, a velocity of 0.02 |v| along z, twice the temperature
    data = {0: moments(one, v * one, 0.0 * one, t0 * one),
            480: moments(late, v * one, np.where(late > 1.0, 0.02 * v, 0.0), np.where(late > 1.0, 2.0, 1.0) * t0)}
    f = SimpleNamespace(iterations=[0, 480],
                        get=lambda it, name: SimpleNamespace(data=data[it][name.split(":")[0]]))
    run = SimpleNamespace(fields=lambda name: f, species_info={"e": {"mass": m}, "i": {"mass": m}})
    h = B.null_history(run, {"v_sw": v})
    assert h["regions"] == ["interior", "x-", "x+", "y-", "y+", "z-", "z+"]
    k = h["regions"].index("x-")
    assert h["e"]["density"][k, 1] == pytest.approx(1.1) and h["e"]["density"][0, 1] == pytest.approx(1.0)
    assert h["e"]["velocity"][k, 1] == pytest.approx(0.02) and h["e"]["velocity"][0, 1] == pytest.approx(0.0, abs=1e-12)
    assert h["e"]["temperature"][k, 1] == pytest.approx(2.0)  # the drift does not count as heat
    assert h["i"]["temperature"][0, 1] == pytest.approx(1.0)
    # the slab at y- crosses the disturbed one: 2 of its 8 planes along x
    assert h["e"]["density"][h["regions"].index("y-"), 1] == pytest.approx(1.0 + 0.1 * 2 / 8)


def test_blast_wave_front_and_debris_extent_are_recovered_from_synthetic_planes():
    from scipy.special import erfinv

    A = _example("blast_wave")
    a = 0.25 + 0.5 * np.arange(256)  # the staggered positions of B_z; centre of the plane at 64
    c, dr = (64.0, 64.0), 0.5
    x, y = np.meshgrid(a - c[0], a - c[1], indexing="ij")
    # a front with four corners: B falls from 2 to 1 between the 4-norm radii 20 and 40, 1.25 at 35,
    # which is the radius 35 on the axes and 35 * 2^(1/4) on the diagonals
    bz = np.clip(2.0 - ((x ** 4 + y ** 4) ** 0.25 - 20.0) / 20.0, 1.0, 2.0)
    fronts = A.ray_fronts(bz, a, a, c, 4.0, dr, 8, A.FRONT_LEVEL)
    assert np.max(np.abs(fronts[::2] - 35.0)) < 0.05 and np.max(np.abs(fronts[1::2] - 35.0 * 2.0 ** 0.25)) < 0.05
    assert np.all(A.ray_fronts(np.ones_like(bz), a, a, c, 4.0, dr, 8, A.FRONT_LEVEL) == 0.0)  # no front
    # a history: the diagonals 10 % ahead of the axes, speed 0.1; R_B = 100, so the window is 30 to 70
    t = np.arange(0.0, 1000.0, 50.0)
    f = 0.1 * t[:, None] * np.where(np.arange(8) % 2, 1.1, 1.0)
    peak = np.where(t < 500.0, 2.0, 3.0)
    m = A.window_measures(t, f, peak, 100.0)
    assert list(t[m["window"]]) == [300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0]
    assert m["front_speed"] == pytest.approx(0.105) and m["asymmetry"] == pytest.approx(0.1 / 1.05)
    assert m["compression"] == 3.0
    one = A.window_measures(t[:7], f[:7], peak[:7], 100.0)  # a single sample in the window: no speed
    assert len(one["window"]) == 1 and np.isnan(one["front_speed"]) and np.isfinite(one["asymmetry"])
    none = A.window_measures(t[:5], f[:5], peak[:5], 100.0)
    assert all(np.isnan(none[k]) for k in ("front_speed", "compression", "asymmetry"))
    # Eq. 3 of Le et al. 2021 on a Gaussian cloud exp(-r^2 / L^2): D98 = L erfinv(0.98)
    n = 0.5 * np.arange(257)
    xn, yn = np.meshgrid(n - 64.0, n - 64.0, indexing="ij")
    r, ring = A.ring_mean(np.exp(-(xn * xn + yn * yn) / 12.0 ** 2), n, n, c, dr)
    assert r[0] == 0.25 and r[-1] == 63.75
    assert A.extent(r, ring) == pytest.approx(12.0 * erfinv(0.98), rel=5e-3)
    assert np.isnan(A.extent(r, np.zeros_like(r))) and np.isnan(A.extent(r, np.where(r < 1.0, np.nan, ring)))


def _blast_wave_slices(cs, n_cells=384, dx=0.5):
    """A run with the slice z = centre written as means over blocks of cs cells: B_z with a cavity inside
    the radius 10 and a ramp from 2 to 1 between 20 and 40, debris ions in a Gaussian shell (radius 20,
    width 5), E_r = 1 / r. The centre is at 96."""
    from types import SimpleNamespace

    stagger = {"Bz": (0.5, 0.5), "n:di": (0.0, 0.0), "Ex": (0.5, 0.0), "Ey": (0.0, 0.5)}

    def read(comp, j):
        x, y = np.meshgrid(*((np.arange(n_cells) + s) * dx - 96.0 for s in stagger[comp]), indexing="ij")
        r = np.maximum(np.hypot(x, y), 1e-12)
        plane = {"Bz": np.where(r < 10.0, 0.0, np.clip(2.0 - (r - 20.0) / 20.0, 1.0, 2.0)),
                 "n:di": np.exp(-((r - 20.0) / 5.0) ** 2), "Ex": x / r ** 2, "Ey": y / r ** 2}[comp]
        plane = plane.reshape(n_cells // cs, cs, n_cells // cs, cs).mean(axis=(1, 3))
        return plane if isinstance(j, int) else np.stack([plane] * 2)

    def plane_axes(comp):
        return tuple((np.arange(n_cells // cs) + ((cs - 1) / 2 + s) / cs) * dx * cs for s in stagger[comp])

    series = SimpleNamespace(layout={"axis": "z"}, step=np.array([0, 10, 20]), read=read, plane_axes=plane_axes)
    return SimpleNamespace(dx=(dx, dx, dx), dt=0.25, series=lambda name: series)


def test_blast_wave_analysis_follows_the_resolution_of_the_slices():
    """The rings are as wide as the samples of the slice are apart, so that every ring holds a sample.
    Rings of the cell width gave nan for slices written as means over blocks of 3 or 4 cells."""
    from scipy.special import erfinv

    A = _example("blast_wave")
    info = {"center": (96.0, 96.0, 96.0), "b0": 1.0, "body": True, "body_radius": 64.0 / 15.0,
            "dipole_radius": 6.4, "d_i": 4.0, "params": {"vte": 0.1}}
    for cs in (1, 3, 4):
        run = _blast_wave_slices(cs)
        h = A.plane_history(run, info)
        assert h["front"] == pytest.approx(35.0, abs=0.5 * cs), cs  # B_z = 1.25 on the ramp; one sample
        assert h["extent"] == pytest.approx(20.0 + 5.0 * erfinv(0.96), abs=0.5 * cs), cs  # Eq. 3 on the shell
        # E_r = 1 / r from the body radius to two dipole radii, three times as far out
        assert A.body_potential(run, info, "slice_z") * 0.1 ** 2 == pytest.approx(np.log(3.0), rel=0.05), cs


def test_blast_wave_gate_fails_without_a_measurement():
    from types import SimpleNamespace

    A = _example("blast_wave")
    good = {"gauss_max": 1e-12, "divB_max": 1e-14, "ambient_deviation": 0.01, "front_speed": 0.1, "fast_speed": 0.08,
            "compression": 2.0, "asymmetry": 0.02, "debris_extent": 30.0, "debris_extent_le": 32.6,
            "debris_extent_step": 880, "last_step": 1760}
    assert all(A.gate_checks(good).values())
    for key, check in (("gauss_max", "gauss"), ("divB_max", "divB"), ("ambient_deviation", "ambient"),
                       ("front_speed", "front_speed"), ("compression", "compression"), ("asymmetry", "symmetry"),
                       ("debris_extent", "debris_extent")):
        checks = A.gate_checks({**good, key: float("nan")})
        assert [k for k, ok in checks.items() if not ok] == [check]
    # the extent still grows at the end of the run: no maximum was seen
    assert not A.gate_checks({**good, "debris_extent_step": 1760})["debris_extent"]

    # ahead of the front: 16^3 blocks of 4 d_e around the centre (32, 32, 32); R_B = 20, so the region
    # begins 4 d_e beyond the front
    info = {"center": (32.0, 32.0, 32.0), "body": False, "r_b": 20.0}
    coords = tuple(2.0 + 4.0 * np.arange(16) for _ in range(3))

    def ahead(front, density):
        """A uniform plasma with the given density at each 3-D output and the front radius at each slice sample."""
        data = {(it, q): np.full((16, 16, 16), d if q.startswith("n") else 0.25 * (q == "Bz"))
                for it, d in density.items() for q in ("n:e", "n:i", "Bx", "By", "Bz")}
        f = SimpleNamespace(iterations=sorted(density),
                            get=lambda it, q: SimpleNamespace(data=data[it, q], coords=coords))
        step = np.array(sorted(front))
        return A.ahead_of_front(SimpleNamespace(fields=lambda name: f), info, step, np.array([front[k] for k in step]))

    out = ahead({0: 0.0, 50: 10.0, 100: 20.0}, {0: 1.0, 100: 1.1})
    assert out["blocks"][0] > out["blocks"][1] > 0 and out["deviation"] == pytest.approx(0.1)
    # the region does not come back when the front is lost later
    assert ahead({0: 0.0, 50: 20.0, 100: 0.0}, {0: 1.0, 100: 1.0})["blocks"] == out["blocks"]
    # a front so far out that no block lies ahead of it: nothing was measured
    out = ahead({0: 0.0, 50: 60.0, 100: 60.0}, {0: 1.0, 100: 1.0})
    assert out["blocks"][1] == 0 and np.isnan(out["deviation"])


def test_filamentation_stable_ratio_compares_the_modes_2_to_10_from_a_tenth_of_the_run():
    spec = importlib.util.spec_from_file_location("filamentation_analyze", ROOT / "examples/filamentation/analyze.py")
    F = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(F)
    t = np.linspace(0.0, 100.0, 101)
    gated = np.where(t < 10.0, 5.0, np.exp(0.5 * (t - 10.0) / 90.0))  # amplitude e^0.5 from t = 10 to the end
    amps = [np.exp(t)] + [gated] * 9 + [np.exp(t)] * 6  # modes 1 and 11..16 are not part of the criterion
    assert F.stable_ratio(t, amps, 100.0) == pytest.approx(np.e)
