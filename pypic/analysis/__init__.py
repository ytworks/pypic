"""Analysis and visualisation of pypic output (CPU only: NumPy, SciPy, h5py, matplotlib).

    from pypic.analysis import open_run
    run = open_run("output/<run>")
    print(run.summary())

Modules: run (reader), spectral (omega-k, growth/damping fits), particles (spectra, power laws,
thermal fits, work along tracks), fields (flux function, fronts), plot (figures), animate
(movies), xdmf (ParaView). Command line: ``python -m pypic.analysis {info,plot,movie,xdmf} RUN``.
"""
from .run import Run, open_run


def __getattr__(name):
    if name == "Movie":
        from .animate import Movie

        return Movie
    raise AttributeError(name)


__all__ = ["Run", "open_run", "Movie"]
