"""Frozen-in initial states (PIC-D24): a vector potential and shaped perturbations."""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys

import cupy as cp
import numpy as np
import pytest

from pypic.parallel.comm import SelfComm
from pypic.simulation import Simulation

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mgpu_frozen as M  # noqa: E402


def _mean_drift(sim):
    out = []
    for sp in sim.species:
        P = sp.rows()
        sel = (P[0] > 30.0) & (P[0] < M.LX - 30.0)  # outside the light cones of the walls up to t = 24
        out.append(float(cp.mean(P[3][sel])))
    return out


def test_potential_is_the_discrete_curl_and_divergence_free():
    """B += curl A_z with A_z = -(b/k) cos(k y) h(x): div B = 0 to round-off, B_x equal to the
    difference of A_z over a cell, and zero on the walls, where h (and so A_z) vanishes."""
    sim = Simulation(M.config(steps=0), comm=SelfComm(), verbose=False)
    g = sim.grid
    b = M.BFRAC * M.B
    assert float(cp.max(cp.abs(g.div_B()))) < 1e-12
    Bx = cp.asnumpy(g.Bx[g.physical_slice()])
    y = np.arange(M.NY) * M.DX
    ref = -(b / M.K) * (np.cos(M.K * (y + M.DX)) - np.cos(M.K * y)) / M.DX
    assert np.max(np.abs(Bx[M.NX // 2, :, 0] - ref)) < 1e-12 * b
    assert np.max(np.abs(Bx[0])) == 0.0
    sim.close()


def test_displacement_loads_the_motional_field_and_the_flow_keeps_its_drift():
    """The electrons displaced along z by xi = -delta E_z / (q n) leave delta E_z = -delta v_x B where h
    is flat; Gauss's law and div B hold over the run, and outside the light cones of the walls each
    species keeps its mean drift -u0 (E x B balance of the uniform part) within 1 %."""
    sim = Simulation(M.config(steps=240), comm=SelfComm(), verbose=False)
    g = sim.grid
    b = M.BFRAC * M.B
    dvx = -M.VA / M.B * b
    Ez = cp.asnumpy(g.Ez[g.physical_slice()])
    y = np.arange(M.NY) * M.DX
    ref = M.V0 * M.B - dvx * M.B * np.sin(M.K * y)
    assert np.max(np.abs(Ez[M.NX // 2, :, 0] - ref)) < 1e-3 * abs(dvx * M.B)
    u0 = M.V0 / math.sqrt(1 - M.V0 * M.V0)
    assert all(abs(d / -u0 - 1) < 1e-3 for d in _mean_drift(sim))
    sim.run()
    assert np.nanmax(sim.history["gauss"]) < 1e-10
    assert np.nanmax(sim.history["divB"]) < 1e-12
    assert all(abs(d / -u0 - 1) < 1e-2 for d in _mean_drift(sim))
    sim.close()


def test_potential_on_a_wall_is_refused():
    """A_z tangential to the reflecting x walls and not zero on them would leave a normal B that the
    wall zeroes: refused."""
    with pytest.raises(ValueError, match="potential"):
        Simulation(M.config(steps=0, a_on_wall=True), comm=SelfComm(), verbose=False)


def test_shaped_displacement_into_a_wall_is_refused():
    """A displacement along the wall normal that does not vanish next to the wall is refused."""
    with pytest.raises(ValueError, match="shaped displacement"):
        Simulation(M.config(steps=0, xi_x=True), comm=SelfComm(), verbose=False)


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 2, reason="needs >= 2 GPUs")
def test_frozen_state_is_decomposition_independent():
    """The shapes use global coordinates only: 1 GPU and 2 GPUs split along the wall axis give the
    same run, and Gauss's law and div B hold."""
    script = os.path.join(HERE, "mgpu_frozen.py")
    outs = {}
    for n in (1, 2):
        r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", str(n), script],
                           capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stderr[-2000:]
        outs[n] = json.loads(r.stdout.strip().splitlines()[-1])
    for k in ("kinetic", "electric", "magnetic"):
        a, b = np.asarray(outs[1][k]), np.asarray(outs[2][k])
        assert np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-30)) < 1e-10, k
    assert np.nanmax(np.asarray(outs[2]["gauss"], dtype=float)) < 1e-10
    assert np.nanmax(np.asarray(outs[2]["divB"], dtype=float)) < 1e-12
