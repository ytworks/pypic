#!/usr/bin/env python3
"""Compile-check every CUDA kernel with NVRTC, without a GPU or CUDA driver.

The engine specialises its kernels at run time with ``#define`` constants that describe the
local grid (``pypic.kernelsrc``). This script builds the same sources for a set of
representative configurations (each boundary type, a decomposed rank, CPML on/off, the
implicit path-segment bound) and compiles them for the supported architectures, so a kernel
that would fail to build on a GPU machine fails CI instead.

NVRTC comes from the ``nvidia-cuda-nvrtc`` wheel (``pip install -e '.[kernelcheck]'``) or a
CUDA toolkit; set PYPIC_NVRTC to a libnvrtc path to override. Exits non-zero on any error.
"""
from __future__ import annotations

import ctypes
import ctypes.util
import importlib.util
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pypic.config import Body  # noqa: E402
from pypic.kernelsrc import (  # noqa: E402
    COMPILE_OPTIONS,
    FIELD_SOURCES,
    PARTICLE_SOURCES,
    assemble_source,
    grid_defines,
    particle_defines,
)
from pypic.parallel.decomp import Decomposition  # noqa: E402

ARCHS = ("compute_80", "compute_120")  # Ampere (portability) and Blackwell (target hardware)
P, O, R = "periodic", "open", "reflecting"
# name, n_cells, rank grid, rank, boundaries, solver, dx, dt, paired injection faces (bit mask), body,
# moving window, lattice of the particle places (PIC-D23; optional)
CASES = [
    ("explicit periodic, 1 rank", (16, 8, 8), (1, 1, 1), 0, ((P, P), (P, P), (P, P)),
     "explicit", 0.1, 0.05, 0, None, False),
    ("explicit open+reflecting, CPML, decomposed rank", (32, 8, 8), (2, 1, 2), 1, ((O, O), (P, P), (R, R)),
     "explicit", 0.1, 0.05, 0, None, False),
    ("implicit reflecting+open, c dt/dx = 2", (32, 16, 16), (2, 1, 1), 0, ((R, R), (P, P), (O, O)),
     "implicit", 0.1, 0.2, 0, None, False),
    ("explicit mixed-precision gather (PIC-D11), decomposed rank", (32, 8, 8), (2, 2, 1), 3,
     ((P, P), (P, P), (R, R)), "mixed", 0.1, 0.05, 0, None, False),
    ("explicit open, paired injection face x- (PIC-D13), decomposed rank", (32, 8, 8), (2, 1, 1), 0,
     ((O, O), (O, O), (O, O)), "explicit", 0.1, 0.05, 1, None, False),
    ("explicit open, absorbing sphere (PIC-D14), mixed-precision gather", (32, 16, 16), (2, 1, 1), 1,
     ((O, O), (O, O), (P, P)), "mixed", 0.1, 0.05, 0, Body(center=(1.6, 0.8, 0.8), radius=0.3), False),
    ("explicit moving window (PIC-D22), open x- and reflecting x+, split along x", (32, 8, 8), (2, 1, 1), 1,
     ((O, R), (R, R), (R, R)), "explicit", 0.1, 0.05, 0, None, True),
    ("explicit moving window (PIC-D22), mixed-precision gather, split along x and y", (32, 8, 8), (2, 2, 1), 0,
     ((O, R), (P, P), (R, R)), "mixed", 0.1, 0.05, 0, None, True),
    ("explicit periodic, particles loaded on a 2 x 2 x 1 lattice (PIC-D23)", (16, 8, 8), (1, 1, 1), 0,
     ((P, P), (P, P), (P, P)), "explicit", 0.1, 0.05, 0, None, False, (2, 2, 1)),
]


def find_nvrtc() -> str:
    env = os.environ.get("PYPIC_NVRTC")
    if env:
        return env
    spec = importlib.util.find_spec("nvidia")
    if spec and spec.submodule_search_locations:
        for base in spec.submodule_search_locations:
            for cand in sorted(Path(base).glob("cu*/lib/libnvrtc.so*")):
                return str(cand)
    lib = ctypes.util.find_library("nvrtc")
    if lib:
        return lib
    raise SystemExit("check_kernels: libnvrtc not found (pip install -e '.[kernelcheck]' or set PYPIC_NVRTC)")


class NVRTC:
    def __init__(self, path: str):
        self.lib = ctypes.CDLL(path)
        self.lib.nvrtcGetErrorString.restype = ctypes.c_char_p
        major, minor = ctypes.c_int(), ctypes.c_int()
        self._check(self.lib.nvrtcVersion(ctypes.byref(major), ctypes.byref(minor)))
        self.version = (major.value, minor.value)

    def _check(self, rc):
        if rc != 0:
            raise RuntimeError(self.lib.nvrtcGetErrorString(rc).decode())

    def compile(self, source: str, name: str, options) -> tuple[bool, str]:
        prog = ctypes.c_void_p()
        self._check(self.lib.nvrtcCreateProgram(ctypes.byref(prog), source.encode(), name.encode(), 0, None, None))
        try:
            opts = (ctypes.c_char_p * len(options))(*[o.encode() for o in options])
            rc = self.lib.nvrtcCompileProgram(prog, len(options), opts)
            size = ctypes.c_size_t()
            self._check(self.lib.nvrtcGetProgramLogSize(prog, ctypes.byref(size)))
            buf = ctypes.create_string_buffer(size.value)
            self._check(self.lib.nvrtcGetProgramLog(prog, buf))
            return rc == 0, buf.value.decode(errors="replace").strip()
        finally:
            self.lib.nvrtcDestroyProgram(ctypes.byref(prog))


def main() -> int:
    nvrtc = NVRTC(find_nvrtc())
    print(f"NVRTC {nvrtc.version[0]}.{nvrtc.version[1]}")
    failures = 0
    for name, n_cells, dims, rank, bc, solver, h, dt, paired, body, window, *lattice in CASES:
        implicit = solver == "implicit"
        ghost = 2 if not implicit else int(math.ceil(dt / h)) + 2
        decomp = Decomposition.build(n_cells, dims, rank, bc, npml=6)
        gdef = grid_defines(decomp, (h, h, h), dt, ghost, use_pml=not implicit and any(O in b for b in bc))
        pdef = particle_defines(gdef, mixed=solver == "mixed", paired_faces=paired, body=body, window=int(window),
                                x_low=decomp.x_low(ghost), lattice=lattice[0] if lattice else None)
        units = [("fields", FIELD_SOURCES, gdef), ("particles", PARTICLE_SOURCES, pdef)]
        for unit, sources, defs in units:
            src = assemble_source(sources, defs)
            for arch in ARCHS:
                ok, log = nvrtc.compile(src, f"{unit}.cu", [*COMPILE_OPTIONS, f"-arch={arch}"])
                status = "ok" if ok else "FAILED"
                print(f"[{status}] {name}: {unit} ({arch})")
                if log:
                    print("    " + log.replace("\n", "\n    "))
                failures += not ok
    if failures:
        print(f"{failures} kernel compilation(s) failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
