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
    """examples/<name>/analyze.py; its ``import theory`` and ``import plan`` must not meet the modules of
    another example."""
    import sys

    folder = ROOT / "examples" / name
    spec = importlib.util.spec_from_file_location(f"{name}_analyze", folder / "analyze.py")
    mod = importlib.util.module_from_spec(spec)
    others = {m: sys.modules.pop(m, None) for m in ("theory", "plan")}
    try:
        spec.loader.exec_module(mod)
    finally:
        for m, other in others.items():
            sys.modules.pop(m, None)
            if other is not None:
                sys.modules[m] = other
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


def _wake_planes(W, p, t, ex, er, x_low=0.0):
    """A stand-in for the output of the plane through the axis: E_x(x, r) on its edges and E_z = +-E_r, in a
    box that begins at x_low."""
    from types import SimpleNamespace

    dx, dz, _ = W.PL.cells(p)
    xs = {"Ex": x_low + (np.arange(p.nx) + 0.5) * dx, "Ez": x_low + np.arange(p.nx) * dx}
    zs = {"Ex": np.arange(p.n_perp) * dz, "Ez": (np.arange(p.n_perp) + 0.5) * dz}
    zc = 0.5 * p.width
    planes = {"Ex": ex(xs["Ex"][:, None], np.abs(zs["Ex"][None, :] - zc)),
              "Ez": er(xs["Ez"][:, None], np.abs(zs["Ez"][None, :] - zc)) * np.sign(zs["Ez"][None, :] - zc)}
    return SimpleNamespace(t=np.array([t]), plane_axes=lambda c, i=None: (xs[c], zs[c]), read=lambda c, i: planes[c])


def test_wakefield_measurement_recovers_a_planted_wake_of_a_pulse():
    """The wake of the theory with 1.03 times its amplitude, 1.004 times its wavenumber and 0.95 times its
    radial field: the fits return these three numbers (the lines at waist / 2 lie between samples)."""
    W = _example("wakefield")
    p = W.PL.WakeParams()
    g = W.PL.laser_geometry(p)
    k, phase = 1.004 / g["v_g"], 0.7

    def potential(x, r):  # written out here, not taken from the planner
        f2 = 1.0 + ((x - g["focus"]) / g["rayleigh"]) ** 2
        a2 = p.a0 ** 2 / f2 / np.sqrt(1.0 - 1.0 / p.omega0 ** 2)
        return (a2 * 0.25 * p.duration * np.sqrt(0.5 * np.pi) * np.exp(-p.duration ** 2 / 8.0)
                * np.exp(-2.0 * r ** 2 / (p.waist ** 2 * f2))), p.waist ** 2 * f2

    def ex(x, r):
        return 1.03 * potential(x, r)[0] / g["v_g"] * np.cos(k * x + phase) + 2e-5

    def er(x, r):
        phi, w2 = potential(x, r)
        return -0.95 * 1.03 * 4.0 * r / w2 * phi * np.sin(k * x + phase) - 1e-5

    m = W._measure(p, _wake_planes(W, p, W.PL.arrival_time(p, p.x_end), ex, er), 0)
    assert m["window"] == pytest.approx([9.0, 48.0], abs=1e-6)
    assert m["amplitude_ratio"] == pytest.approx(1.03, rel=1e-3)
    assert m["k_times_speed"] == pytest.approx(1.004, rel=1e-5)
    assert m["phase"] == pytest.approx(phase, abs=1e-3) and m["offset"] == pytest.approx(2e-5, rel=1e-2)
    assert m["Er_over_Ex_ratio"] == pytest.approx(0.95, rel=2e-3)
    # the same wake seen by a moving window whose lower face is at x = 39.74: the fit keeps 1 c/omega_pe
    # from the face and has 1.2 wavelengths left, from which it returns the same numbers
    p.nx, p.focus = 1152, g["focus"]
    m = W._measure(p, _wake_planes(W, p, W.PL.arrival_time(p, p.x_end), ex, er, x_low=39.74), 0)
    assert m["window"] == pytest.approx([40.74, 48.0], abs=1e-6)
    assert m["amplitude_ratio"] == pytest.approx(1.03, rel=1e-3)
    assert m["k_times_speed"] == pytest.approx(1.004, rel=1e-4)
    assert m["Er_over_Ex_ratio"] == pytest.approx(0.95, rel=2e-3)


def test_wakefield_measurement_recovers_a_planted_wake_of_a_bunch():
    """A wake with 0.97 times the amplitude of the theory on the axis, the phase 0.1 behind the bunch
    and 0.7 of the axis value at r = sigma_r."""
    W = _example("wakefield")
    p = W.PL.WakeParams(driver="bunch", nx=900, n_perp=192, dx=0.125, ramp=(2.0, 5.0, 107.0, 110.0), x_end=100.0)
    beta = W.PL.speed(p)
    t = 88.5
    xb = p.start + beta * t
    amp = 0.97 * W.PL.wake_amplitude(p, 0.0)[0][0]

    def ex(x, r):
        return amp * 0.7 ** (r / p.sigma_r) ** 2 * np.cos(0.99 / beta * (x - xb) + 0.1)

    m = W._measure(p, _wake_planes(W, p, t, ex, lambda x, r: 0.0 * x * r), 0)
    assert m["window"] == pytest.approx([xb - 5.0 - 10.0 * np.pi * beta, xb - 5.0])
    assert m["amplitude_ratio"] == pytest.approx(0.97, rel=1e-4) and m["k_times_speed"] == pytest.approx(0.99, rel=1e-5)
    assert m["phase"] == pytest.approx(0.1, abs=1e-4) and m["R_ratio"] == pytest.approx(0.7, rel=1e-4)


@pytest.mark.parametrize("config, wavelengths", [("laser_linear", 1.0), ("beam_linear", 5.0)])
def test_wakefield_gate_keeps_the_driver_and_the_fit_inside_the_moving_window(config, wavelengths):
    """The gate configurations run in a moving window: the engine accepts them, the window starts when the
    driver is ``window_lead`` behind the x+ wall, and at the two outputs that are analysed the interval of the
    fit lies in the box behind the driver and holds at least this many wavelengths of the wake."""
    import yaml

    from pypic import antenna, window

    W = _example("wakefield")
    with open(ROOT / f"examples/wakefield/configs/{config}.yaml") as f:
        p = W.PL.params_from_dict(yaml.safe_load(f))
    cfg, report, info = W.PL.plan_wakefield(p, n_gpus=8)
    assert window.check_config(cfg) and antenna.check_config(cfg) == (p.driver == "laser")
    length = cfg.n_cells[0] * cfg.dx[0]
    assert W.PL.driver_position(p, cfg.window.start) == pytest.approx(length - p.window_lead, abs=1e-9)
    every = cfg.n_steps // p.n_out
    for steps in (cfg.n_steps - every, cfg.n_steps):
        low = window.shift_after(cfg, steps) * cfg.dx[0]
        x_driver = W.PL.driver_position(p, steps * cfg.dt)
        lo, hi = W._window(p, x_driver, low)
        assert low + 1.0 <= lo < hi < x_driver < low + length - p.window_margin
        assert hi - lo >= wavelengths * 2.0 * np.pi * W.PL.speed(p) - 1e-9
        assert p.ramp[1] < lo and hi < p.ramp[2]  # on the plateau
    assert info["x_low_end"] == pytest.approx(window.shift_after(cfg, cfg.n_steps) * cfg.dx[0])
