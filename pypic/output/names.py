"""Field names, openPMD records and Yee staggering (pure Python, shared by the GPU-side
diagnostics and the CPU-side reader)."""
from __future__ import annotations

# staggering offsets (x, y, z) of each Yee component, in cells
STAGGER = {
    "Ex": (0.5, 0, 0), "Ey": (0, 0.5, 0), "Ez": (0, 0, 0.5),
    "Bx": (0, 0.5, 0.5), "By": (0.5, 0, 0.5), "Bz": (0.5, 0.5, 0),
}
VECTORS = {"E": ("Ex", "Ey", "Ez"), "B": ("Bx", "By", "Bz"), "J": ("Jx", "Jy", "Jz")}
MOMENT_KIND = {"n": 0, "fluxx": 1, "fluxy": 2, "fluxz": 3, "ekin": 4}


def expand_fields(names) -> list[str]:
    """'E' -> Ex Ey Ez, 'flux:e' -> fluxx:e fluxy:e fluxz:e; components keep their order."""
    out = []
    for f in names:
        if f in VECTORS:
            out += list(VECTORS[f])
        elif f.startswith("flux:"):
            out += [f"flux{c}:{f[5:]}" for c in "xyz"]
        else:
            out.append(f)
    return out


def record_of(comp: str) -> tuple[str, str]:
    """openPMD (record, component) of an expanded field name ('' for scalar records)."""
    if len(comp) == 2 and comp[0] in "EBJ" and comp[1] in "xyz":
        return comp[0], comp[1]
    if comp == "rho":
        return "rho", ""
    base, sp = comp.split(":", 1)
    if base.startswith("flux"):
        return f"flux_{sp}", base[4]
    return f"{base}_{sp}", ""


def stagger_of(comp: str) -> tuple:
    """Staggering of a component in cells (J like E; densities and moments at nodes)."""
    if comp in STAGGER:
        return STAGGER[comp]
    if len(comp) == 2 and comp[0] == "J":
        return STAGGER["E" + comp[1]]
    return (0.0, 0.0, 0.0)
