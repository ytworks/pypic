"""Linear theory (pypic.theory) against literature values and exact limits."""
from __future__ import annotations

import math

import numpy as np

from pypic.theory import streaming as S
from pypic.theory import waves as W


def test_landau_damping_roots_match_literature():
    # Maxwellian Langmuir roots, omega / omega_pe at k lambda_De = 0.3 / 0.4 / 0.5
    # (e.g. Canosa 1973; Birdsall & Langdon): 1.1598-0.0126i, 1.2850-0.0661i, 1.4157-0.1534i
    vt = 0.05
    w = W.langmuir(np.array([0.3, 0.4, 0.5]) / vt, vt)
    ref = np.array([1.1598 - 0.0126j, 1.2850 - 0.0661j, 1.4157 - 0.1534j])
    assert np.max(np.abs(w - ref)) < 2e-4


def test_cold_magnetised_cutoffs_and_bernstein_limit():
    wce = 0.5
    p = W.cold_perpendicular(np.array([1e-6]), 1.0, wce)
    # X-mode cutoffs: omega_R,L = (+-wce + sqrt(wce^2 + 4 wpe^2)) / 2
    assert abs(p["X_upper"][0] - (wce + math.sqrt(wce ** 2 + 4)) / 2) < 1e-9
    assert abs(p["X_lower"][0] - (-wce + math.sqrt(wce ** 2 + 4)) / 2) < 1e-9
    # long-wavelength Bernstein branch tends to the upper hybrid frequency
    b = W.bernstein(np.array([1e-3]), 0.05, 1.0, wce, n_bands=2)
    assert abs(b[1, 0] - math.sqrt(1 + wce ** 2)) < 1e-6
    # whistler: omega -> k^2 c^2 wce / wpe^2 for k << wpe/c (corrections O(omega/wce))
    w = W.cold_parallel(np.array([0.05]), 1.0, 2.0)["whistler"][0]
    assert abs(w / (0.05 ** 2 * 2.0) - 1) < 0.02


def test_streaming_limits():
    assert abs(S.buneman_cold(1836.15267) - 0.0561) < 1e-4
    u = 3.0  # filamentation growth tends to beta wp / sqrt(gamma) for k >> wp/c
    g = math.sqrt(1 + u * u)
    assert abs(S.filamentation_cold(1e4, u) - (u / g) / math.sqrt(g)) < 1e-6
