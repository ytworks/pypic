"""Command line access to a pypic run directory.

    python -m pypic.analysis info  output/run
    python -m pypic.analysis plot  output/run [--out figs/]
    python -m pypic.analysis movie output/run --panel phase:phase:electron --panel xt:lineout:Ex \
                                              --panel history -o movie.mp4 [--fps 20 --frames 300]
    python -m pypic.analysis xdmf  output/run --diag fields

Panel specs: slice:<diag>:<comp>  cut:<diag>:<comp>[:<axis>]  lineout:<diag>:<comp>[,<comp>]
             xt:<diag>:<comp>  phase:<diag>:<species>  spectrum:<diag>:<species>  history
             tracks:<diag>:<species>[:<qx>:<qy>]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

from .run import open_run


def cmd_info(a):
    print(open_run(a.run).summary())


def cmd_plot(a):
    from . import plot as P
    from .particles import spectrum

    run = open_run(a.run)
    out = a.out or os.path.join(a.run, "figures")
    os.makedirs(out, exist_ok=True)
    P.style()
    made = []
    fig, ax = P.plt.subplots(1, 2, figsize=(11, 3.8))
    P.energy_history(ax[0], run)
    P.conservation(ax[1], run)
    made.append(P.save(fig, os.path.join(out, "history.png")))
    for name, d in run.diag.items():
        k = d["kind"]
        if k == "lineout":
            s = run.series(name)
            for c in s.components:
                fig, ax = P.plt.subplots(figsize=(6, 4))
                x = s.x(c)
                t = s.t
                P.field_map(ax, s[c].T, [x[0], x[-1], t[0], t[-1]], signed=c[:1] in "EBJ", label=c,
                            xlabel=f"{d['axis']} [c/ωpe]", ylabel=r"$\omega_{pe} t$")
                ax.set_title(f"{name}: {c}(x, t)")
                made.append(P.save(fig, os.path.join(out, f"{name}_{c.replace(':', '_')}_xt.png")))
        elif k == "histogram":
            s = run.series(name)
            for sp in s.keys():
                h = s.read(sp, -1)
                fig, ax = P.plt.subplots(figsize=(6, 4))
                if h.ndim == 2:
                    P.phase_space(ax, h, s.edges, xlabel=s.quantities[0], ylabel=s.quantities[1])
                else:
                    E, f = spectrum(s.edges[0], h)
                    P.spectrum_plot(ax, E, f, label=f"t={s.t[-1]:.1f}", xlabel=s.quantities[0])
                ax.set_title(f"{name}: {sp} at t={s.t[-1]:.1f}")
                made.append(P.save(fig, os.path.join(out, f"{name}_{sp}.png")))
        elif k == "kspectrum":
            s = run.series(name)
            fig, ax = P.plt.subplots(figsize=(6, 4))
            for i, c in enumerate(s.keys()):
                v = s.read(c, -1)
                ok = v > 0
                ax.loglog(s.k[ok], v[ok], color=P.color(i), label=c)
            ax.set_xlabel(r"$|k| c/\omega_{pe}$")
            ax.set_ylabel("power")
            ax.legend()
            made.append(P.save(fig, os.path.join(out, f"{name}.png")))
    print("\n".join(made))


def parse_panel(movie, spec: str):
    parts = spec.split(":")
    kind = parts[0]
    if kind == "history":
        return movie.history()
    if kind == "slice":
        return movie.slice(parts[1], parts[2])
    if kind == "cut":
        return movie.cut(parts[1], parts[2], parts[3] if len(parts) > 3 else "z")
    if kind == "lineout":
        return movie.lineout(parts[1], parts[2].split(","))
    if kind == "xt":
        return movie.xt(parts[1], parts[2])
    if kind == "phase":
        return movie.phase(parts[1], parts[2])
    if kind == "spectrum":
        return movie.spectrum(parts[1], parts[2])
    if kind == "tracks":
        q = parts[3:5] if len(parts) >= 5 else ("x", "y")
        return movie.tracks(parts[1], parts[2], *q)
    raise SystemExit(f"unknown panel kind {kind!r}")


def cmd_movie(a):
    from .animate import Movie

    run = open_run(a.run)
    m = Movie(run, title=a.title or os.path.basename(os.path.normpath(a.run)))
    for p in a.panel:
        parse_panel(m, p)
    print(m.save(a.o, fps=a.fps, n_frames=a.frames))


def cmd_xdmf(a):
    from .xdmf import write_xdmf

    run = open_run(a.run)
    names = [a.diag] if a.diag else run.names("fields")
    for n in names:
        print(write_xdmf(run, n))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m pypic.analysis", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("info")
    p.add_argument("run")
    p.set_defaults(fn=cmd_info)
    p = sub.add_parser("plot")
    p.add_argument("run")
    p.add_argument("--out")
    p.set_defaults(fn=cmd_plot)
    p = sub.add_parser("movie")
    p.add_argument("run")
    p.add_argument("--panel", action="append", required=True)
    p.add_argument("-o", default="movie.mp4")
    p.add_argument("--fps", type=int, default=15)
    p.add_argument("--frames", type=int, default=300)
    p.add_argument("--title")
    p.set_defaults(fn=cmd_movie)
    p = sub.add_parser("xdmf")
    p.add_argument("run")
    p.add_argument("--diag")
    p.set_defaults(fn=cmd_xdmf)
    a = ap.parse_args(argv)
    np.seterr(all="ignore")
    a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
