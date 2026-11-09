"""The GPU suite needs CuPy and a CUDA device; without them it is not collected at all
(the reason is shown in the pytest header), so `pytest` on a CPU machine runs tests/cpu only."""
try:
    import cupy as cp

    _NDEV = cp.cuda.runtime.getDeviceCount()
    _REASON = None if _NDEV > 0 else "no CUDA device"
except ImportError:
    _REASON = "CuPy not installed (pip install -e '.[gpu]')"
except Exception as e:  # driver missing or not initialisable
    _REASON = f"CUDA not usable ({type(e).__name__})"

collect_ignore_glob = ["*.py"] if _REASON else []


def pytest_report_header(config):
    if _REASON:
        return f"tests/gpu not collected: {_REASON}"
