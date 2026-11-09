"""pypic: GPU-accelerated, multi-GPU 3-D relativistic electromagnetic Particle-in-Cell engine."""
from .config import Beam, ImplicitConfig, SimConfig, SpeciesConfig

__all__ = ["SimConfig", "SpeciesConfig", "ImplicitConfig", "Beam", "Simulation"]
__version__ = "0.1.0"


def __getattr__(name):
    if name == "Simulation":
        from .simulation import Simulation

        return Simulation
    raise AttributeError(name)
