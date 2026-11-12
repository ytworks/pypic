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


def particle_defines(
    grid_defs: dict, mixed: bool = False, paired_faces: int = 0, body=None, compact: bool = False
) -> dict:
    """Particle kernels share the grid constants (no CPML terms) plus the path-segment bound.
    ``mixed`` (PIC-D11): the explicit push gathers E and B from float32 copies; OXI.. are the
    global cell indices of local index 0, so the gather weights come from global coordinates and
    are the same on every rank layout. ``paired_faces`` (PIC-D13): bit ``2 * axis + side`` per
    paired injection face; 0 compiles the push without the return of the neutraliser. ``body``
    (PIC-D14): the absorbing sphere of the configuration, or None. ``compact`` (PIC-D19): the
    particles are stored as cell, position within the cell and float32 momenta; the gather is the
    one of ``mixed``, GNX.. are the global cells per axis and S?0I / S?1I the physical range of
    the rank in global cells."""
    D = dict(grid_defs)
    D.pop("USE_PML", None)
    D["MAXSEG"] = 3 * (D["GH"] + 2) + 4
    D["PAIRED_FACES"] = int(paired_faces)
    D.update(body_defines(body))
    if mixed or compact:
        D["MIXED_GATHER"] = 1
        for c in "XYZ":
            D[f"O{c}I"] = int(round(D[f"O{c}"] / D[f"D{c}"]))
    if compact:
        D["COMPACT_STORE"] = 1
        if D["NX"] * D["NY"] * D["NZ"] >= 2 ** 31:
            raise ValueError("precision='compact': the cells of one rank (with ghosts) must be fewer than 2**31, "
                             "the range of the stored cell; use more ranks")
        for c in "XYZ":
            D[f"GN{c}"] = int(round(D[f"L{c}"] / D[f"D{c}"]))
            if D[f"GN{c}"] + 2 * D["GH"] > 2 ** 20:
                raise ValueError(f"precision='compact': more than 2**20 cells along {c.lower()} (with ghosts); the "
                                 "index coordinate would not be exact in 2**-32 of a cell")
            D[f"S{c}0I"] = int(round(D[f"S{c}0"] / D[f"D{c}"]))
            D[f"S{c}1I"] = int(round(D[f"S{c}1"] / D[f"D{c}"]))
    return D
