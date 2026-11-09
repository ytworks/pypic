"""GPU backend: CuPy import and compilation of the CUDA kernels with injected constants."""
from __future__ import annotations

import os
from pathlib import Path

import cupy as cp

KERNEL_DIR = Path(__file__).parent / "kernels"
_cache: dict = {}


def _fmt(v) -> str:
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(float(v))
    return str(v)


def compile_module(sources: list[str], defines: dict) -> cp.RawModule:
    """Compile the given kernel files (from pypic/kernels) with ``defines`` as macros."""
    header = (KERNEL_DIR / "common.cuh").read_text()
    body = "\n".join((KERNEL_DIR / s).read_text() for s in sources)
    macro = "\n".join(f"#define {k} {_fmt(v)}" for k, v in sorted(defines.items()))
    code = macro + "\n" + header.replace("#pragma once", "") + "\n" + body
    key = (cp.cuda.Device().id, code)
    if key not in _cache:
        opts = ("--std=c++17", "-use_fast_math") if os.environ.get("PYPIC_FAST_MATH") else ("--std=c++17",)
        _cache[key] = cp.RawModule(code=code, options=opts)
    return _cache[key]


def launch(kernel, n: int, args, block: int = 256):
    if n <= 0:
        return
    grid = (int((n + block - 1) // block),)
    kernel(grid, (block,), args)
