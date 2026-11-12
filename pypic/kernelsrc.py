"""CUDA kernel sources and their compile-time constants (pure Python, no GPU needed).

The engine compiles its kernels with NVRTC at run time, specialised by ``#define`` macros
that describe the local grid. Building the source here, without CuPy, lets CPU-only CI
compile-check exactly what a GPU run would compile (``scripts/check_kernels.py``).
"""
from __future__ import annotations

from pathlib import Path

from .body import defines as body_defines
from .parallel.decomp import BC_CODES, Decomposition

KERNEL_DIR = Path(__file__).parent / "kernels"
FIELD_SOURCES = ("fields.cu",)
PARTICLE_SOURCES = ("particles.cu", "implicit.cu", "migrate.cu")
COMPILE_OPTIONS = ("--std=c++17",)
CG_BLOCK = 256  # threads per block of cg_matvec_end, which sums over each block (a power of two)


def _fmt(v) -> str:
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(float(v))
    return str(v)


def assemble_source(sources, defines: dict) -> str:
    """Macro block + common header + the given kernel files, as passed to NVRTC."""
    header = (KERNEL_DIR / "common.cuh").read_text()
    body = "\n".join((KERNEL_DIR / s).read_text() for s in sources)
    macro = "\n".join(f"#define {k} {_fmt(v)}" for k, v in sorted(defines.items()))
    return macro + "\n" + header.replace("#pragma once", "") + "\n" + body


def local_shape(decomp: Decomposition, ghost: int) -> tuple:
    return tuple(a.n_int + 2 * ghost for a in decomp.axes)


def grid_defines(decomp: Decomposition, dx, dt: float, ghost: int, use_pml: bool) -> dict:
    """Constants describing this rank's local Yee grid (field kernels)."""
    ax = decomp.axes
    shape = local_shape(decomp, ghost)
    L = tuple(n * h for n, h in zip(decomp.n_cells, dx))
    D = {"NX": shape[0], "NY": shape[1], "NZ": shape[2], "GH": ghost,
         "DX": float(dx[0]), "DY": float(dx[1]), "DZ": float(dx[2]), "DT": float(dt),
         "LX": L[0], "LY": L[1], "LZ": L[2], "CG_BLOCK": CG_BLOCK}
    for a, c in enumerate("XYZ"):
        D[f"O{c}"] = (ax[a].cell_lo - ghost) * float(dx[a])
        D[f"S{c}0"] = ax[a].phys_lo * float(dx[a])
        D[f"S{c}1"] = (ax[a].phys_lo + ax[a].n_phys) * float(dx[a])
        D[f"BC{c}0"] = BC_CODES[decomp.bc[a][0]]
        D[f"BC{c}1"] = BC_CODES[decomp.bc[a][1]]
        D[f"DEC{c}"] = 1 if ax[a].nproc > 1 else 0
        D[f"I0{c}"] = ghost
        D[f"I1{c}"] = ghost + ax[a].n_int
    if use_pml:
        D["USE_PML"] = 1
    return D


def particle_defines(grid_defs: dict, mixed: bool = False, paired_faces: int = 0, body=None,
                     window: int = 0, x_low: int = 0, lattice=None) -> dict:
    """Particle kernels share the grid constants (no CPML terms) plus the path-segment bound.
    ``mixed`` (PIC-D11): the explicit push gathers E and B from float32 copies; OXI.. are the
    global cell indices of local index 0, so the gather weights come from global coordinates and
    are the same on every rank layout. ``paired_faces`` (PIC-D13): bit ``2 * axis + side`` per
    paired injection face; 0 compiles the push without the return of the neutraliser. ``body``
    (PIC-D14): the absorbing sphere of the configuration, or None. ``window`` (PIC-D22): 0 compiles
    the kernels of a box at rest; a moving window passes the number of its simulation, because the
    kernels then read the shift and the faces along x from the device array ``window_state`` of
    their module, and two simulations of one process must not share a module. ``x_low`` (with a
    window): the cells that the arrays of the first rank hold below the box along x, the same
    number on every rank; with OXI it gives the kernels an index coordinate along x that depends on
    the rank only by a whole number of cells (XIDX in ``common.cuh``). ``lattice``
    (PIC-D23): the places per axis of the particles of a cell, or None for random places."""
    D = dict(grid_defs)
    D.pop("USE_PML", None)
    D["MAXSEG"] = 3 * (D["GH"] + 2) + 4
    D["PAIRED_FACES"] = int(paired_faces)
    D.update(body_defines(body))
    if window:
        D["MOVING_WINDOW"] = int(window)
        D["OXI"] = int(round(D["OX"] / D["DX"]))  # the global cell of local index 0
        D["XLOW"] = int(x_low)
        if D["OXI"] + D["XLOW"] < 0:
            raise ValueError("x_low is below the cells that this rank holds below the box along x")
    if lattice is not None:
        for c, v in zip("XYZ", lattice):
            D[f"LATTICE_{c}"] = int(v)
    if mixed:
        D["MIXED_GATHER"] = 1
        for c in "XYZ":
            D[f"O{c}I"] = int(round(D[f"O{c}"] / D[f"D{c}"]))
    return D
