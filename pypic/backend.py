"""GPU backend: CuPy import and compilation of the CUDA kernels with injected constants."""
from __future__ import annotations

import os

import cupy as cp

from .kernelsrc import COMPILE_OPTIONS, assemble_source

_cache: dict = {}


def compile_module(sources, defines: dict) -> cp.RawModule:
    """Compile the given kernel files (from pypic/kernels) with ``defines`` as macros."""
    code = assemble_source(sources, defines)
    key = (cp.cuda.Device().id, code)
    if key not in _cache:
        opts = COMPILE_OPTIONS + (("-use_fast_math",) if os.environ.get("PYPIC_FAST_MATH") else ())
        _cache[key] = cp.RawModule(code=code, options=opts)
    return _cache[key]


def launch(kernel, n: int, args, block: int = 256):
    if n <= 0:
        return
    grid = (int((n + block - 1) // block),)
    kernel(grid, (block,), args)
