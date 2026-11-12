"""Animations of pypic runs: several synchronised panels, one frame per output time.

    from pypic.analysis import open_run, Movie
    run = open_run("output/two_stream_ei")
    m = Movie(run, title="e-i two-stream")
    m.phase("phase", "electron")        # (x, u_x) histogram
    m.xt("lineout", "Ex")               # space-time diagram growing with time
    m.history()                         # energy history with a time cursor
    m.save("two_stream.mp4", fps=20)    # MP4 (imageio-ffmpeg) or .gif (Pillow)

Diagnostics written at different cadences are shown at the latest output not after the frame
time. Colour limits are fixed over the whole movie (robust percentiles) so frames compare.
"""
from __future__ import annotations

import math
import os
import warnings

import numpy as np

from . import plot as P

plt = P.plt


class Panel:
    title = ""

    def times(self):  # output times this panel can show (None: time independent)
        return None

    def prepare(self, frame_times):
        pass

    def init(self, ax):
        raise NotImplementedError

    def update(self, t):
        raise NotImplementedError

    @staticmethod
    def _index(times, t):
        i = int(np.searchsorted(times, t + 1e-9, side="right")) - 1
        return max(0, min(i, len(times) - 1))


class SlicePanel(Panel):
    """A component of a slice diagnostic. `scale` multiplies the data. `blank=(centre, radius)` leaves
    out the part of the plane inside a sphere (an absorbing body: the dipole field in there would set
    the colour scale). `vlim=(low, high)` fixes the colour range; the default is the robust range of
    the frames. With a moving window every frame shows the plane at its laboratory place (OUT-D5)."""

    def __init__(self, run, diag, comp, signed=None, log=False, label=None, scale=1.0, blank=None, vlim=None):
        self.s = run.series(diag)
        self.moving = bool(np.any(self.s.window_shift))
        self.comp = comp
        self.signed = (comp[:1] in "EBJ" and len(comp) == 2) if signed is None else signed
        self.log = log
        self.label = label or comp
        self.t = self.s.t
        axis = self.s.layout["axis"]
        self.xl, self.yl = [c for c in "xyz" if c != axis]
        self.scale, self.vlim, self.mask = scale, vlim, None
        if blank is not None:
            centre, radius = blank
            a, b = self.s.plane_axes(comp)
            ca, cb = (c for c, ax in zip(centre, "xyz") if ax != axis)
            cn = centre["xyz".index(axis)]
            off = self.s.layout.get("plane_position", {}).get(comp, cn) - cn
            self.mask = (a[:, None] - ca) ** 2 + (b[None, :] - cb) ** 2 + off ** 2 < radius ** 2

    def times(self):
        return self.t

    def _read(self, i):
        d = self.s.read(self.comp, i)
        if self.scale != 1.0 or self.mask is not None:
            d = np.array(d, dtype=float) * self.scale
            if self.mask is not None:
                d[self.mask] = np.nan
        return d

    def prepare(self, frame_times):
        if self.vlim is None:
            idx = sorted({self._index(self.t, t) for t in frame_times[:: max(1, len(frame_times) // 12)]})
            self.vlim = P.robust_limits([self._read(i) for i in idx], self.signed, self.log)

    def init(self, ax):
        self.ax = ax
        d = self._read(0)
        self.im = P.field_map(ax, d, self.s.extent(self.comp), signed=self.signed, log=self.log, vlim=self.vlim,
                              label=self.label, xlabel=f"{self.xl} [c/ωpe]", ylabel=f"{self.yl} [c/ωpe]")

    def update(self, t):
        i = self._index(self.t, t)
        d = self._read(i)
        if self.log:
            d = np.where(d > 0, d, np.nan)
        self.im.set_data(np.asarray(d).T)
        if self.moving:
            ext = self.s.extent(self.comp, i)
            self.im.set_extent(ext)
            self.ax.set_xlim(ext[0], ext[1])
        self.ax.set_title(f"{self.label}  t={self.t[i]:.1f}")


class FieldCutPanel(Panel):
    """A plane of a 3-D fields dump (fields diagnostic), e.g. FieldCutPanel(run, 'fields', 'Bz', 'z')."""

    def __init__(self, run, diag, comp, axis="z", position=None, signed=None, log=False):
        self.f = run.fields(diag)
        self.comp, self.axis = comp, "xyz".index(axis)
        self.position = position
        self.signed = (comp[:1] in "EBJ" and len(comp) == 2) if signed is None else signed
        self.log = log
        self.its = self.f.iterations
        self.t = np.array([it * run.dt for it in self.its])
        self.cache = {}

    def times(self):
        return self.t

    def _cut(self, i):
        if i not in self.cache:
            fd = self.f.get(self.its[i], self.comp)
            c = fd.coords[self.axis]
            j = len(c) // 2 if self.position is None else int(np.argmin(np.abs(c - self.position)))
            self.cache[i] = (np.take(fd.data, j, axis=self.axis),
                             [fd.coords[b] for b in range(3) if b != self.axis])
        return self.cache[i]

    def prepare(self, frame_times):
        idx = sorted({self._index(self.t, t) for t in frame_times[:: max(1, len(frame_times) // 8)]})
        self.vlim = P.robust_limits([self._cut(i)[0] for i in idx], self.signed, self.log)

    def init(self, ax):
        self.ax = ax
        d, (a, b) = self._cut(0)
        ext = [a[0], a[-1], b[0], b[-1]]
        lab = [c for k, c in enumerate("xyz") if k != self.axis]
        self.im = P.field_map(ax, d, ext, signed=self.signed, log=self.log, vlim=self.vlim, label=self.comp,
                              xlabel=f"{lab[0]} [c/ωpe]", ylabel=f"{lab[1]} [c/ωpe]")

    def update(self, t):
        i = self._index(self.t, t)
        d, _ = self._cut(i)
        self.im.set_data(np.asarray(d).T)
        self.ax.set_title(f"{self.comp}  t={self.t[i]:.1f}")


class LineoutPanel(Panel):
    def __init__(self, run, diag, comps, ylim=None):
        self.s = run.series(diag)
        self.comps = [comps] if isinstance(comps, str) else list(comps)
        self.t = self.s.t
        self.data = {c: self.s[c] for c in self.comps}
        self.ylim = ylim

    def times(self):
        return self.t

    def init(self, ax):
        self.ax = ax
        self.lines = [ax.plot(self.s.x(c), self.data[c][0], color=P.color(k), label=c)[0]
                      for k, c in enumerate(self.comps)]
        lo = min(np.nanmin(v) for v in self.data.values())
        hi = max(np.nanmax(v) for v in self.data.values())
        pad = 0.05 * (hi - lo or 1.0)
        ax.set_ylim(self.ylim or (lo - pad, hi + pad))
        ax.set_xlabel(f"{self.s.layout['axis']} [c/ωpe]")
        if len(self.comps) > 1:
            ax.legend(loc="upper right")

    def update(self, t):
        i = self._index(self.t, t)
        for c, ln in zip(self.comps, self.lines):
            ln.set_ydata(self.data[c][i])
        self.ax.set_title(f"{', '.join(self.comps)}  t={self.t[i]:.1f}")


class XTPanel(Panel):
    """Space-time diagram of a lineout, revealed up to the frame time."""

    def __init__(self, run, diag, comp, signed=None):
        self.s = run.series(diag)
        self.comp = comp
        self.t = self.s.t
        self.data = self.s[comp]
        self.signed = (comp[:1] in "EBJ" and len(comp) == 2) if signed is None else signed

    def times(self):
        return self.t

    def init(self, ax):
        self.ax = ax
        x = self.s.x(self.comp)
        ext = [x[0], x[-1], self.t[0], self.t[-1]]
        self.full = np.asarray(self.data, dtype=float)
        vlim = P.robust_limits([self.full], self.signed)
        self.im = P.field_map(ax, np.full_like(self.full, np.nan).T, [ext[0], ext[1], ext[2], ext[3]],
                              signed=self.signed, vlim=vlim, label=self.comp,
                              xlabel=f"{self.s.layout['axis']} [c/ωpe]", ylabel=r"$\omega_{pe} t$")
        self.ax.set_title(f"{self.comp}(x, t)")

    def update(self, t):
        i = self._index(self.t, t)
        img = np.full_like(self.full, np.nan)
        img[: i + 1] = self.full[: i + 1]
        self.im.set_data(img)


class PhasePanel(Panel):
    def __init__(self, run, diag, species, log=True):
        self.s = run.series(diag)
        self.sp = species
        self.t = self.s.t
        self.edges = self.s.edges
        self.q = self.s.quantities
        self.log = log
        # x of a histogram is the place in the window (OUT-D5): a frame shows its laboratory place, on
        # whichever axis of the picture x is
        self.ix = list(self.q).index("x") if "x" in self.q and np.any(self.s.window_shift) else None
        self.shift = self.s.window_shift * run.dx[0]

    def times(self):
        return self.t

    def prepare(self, frame_times):
        idx = sorted({self._index(self.t, t) for t in frame_times[:: max(1, len(frame_times) // 12)]})
        self.vlim = P.robust_limits([self.s.read(self.sp, i) for i in idx], False, self.log)

    def init(self, ax):
        self.ax = ax
        units = {"x": "[c/ωpe]", "y": "[c/ωpe]", "z": "[c/ωpe]", "ux": "/ c", "uy": "/ c", "uz": "/ c"}
        self.im = P.phase_space(ax, self.s.read(self.sp, 0), self.edges, log=self.log, vlim=self.vlim,
                                xlabel=f"{self.q[0]} {units.get(self.q[0], '')}",
                                ylabel=f"{self.q[1]} {units.get(self.q[1], '')}")

    def update(self, t):
        i = self._index(self.t, t)
        d = self.s.read(self.sp, i)
        if self.log:
            d = np.where(d > 0, d, np.nan)
        self.im.set_data(np.asarray(d).T)
        if self.ix is not None:
            ext = [[e[0], e[-1]] for e in self.edges]
            ext[self.ix] = [v + self.shift[i] for v in ext[self.ix]]
            self.im.set_extent(ext[0] + ext[1])
            (self.ax.set_ylim if self.ix else self.ax.set_xlim)(*ext[self.ix])
        self.ax.set_title(f"{self.sp} ({self.q[0]}, {self.q[1]})  t={self.t[i]:.1f}")


class SpectrumPanel(Panel):
    """Energy spectrum (1-D histogram) at the frame time, with the initial spectrum for reference."""

    def __init__(self, run, diag, species, xlabel=None):
        from .particles import spectrum

        self.s = run.series(diag)
        self.sp = species
        self.t = self.s.t
        self.E, f0 = spectrum(self.s.edges[0], self.s.read(species, 0))
        self.f = [spectrum(self.s.edges[0], self.s.read(species, i))[1] for i in range(len(self.t))]
        self.f0 = f0
        q = self.s.quantities[0]
        self.xlabel = xlabel or {"ekin": r"$E_{kin}$ [$m_e c^2$]", "gamma": r"$\gamma$", "u": r"$|u|$"}.get(q, q)

    def times(self):
        return self.t

    def init(self, ax):
        self.ax = ax
        P.style()
        ok = self.f0 > 0
        ax.loglog(self.E[ok], self.f0[ok], color=P.INK["muted"], lw=1.0, label="t = 0")
        (self.line,) = ax.loglog(self.E, np.where(self.f[0] > 0, self.f[0], np.nan), color=P.color(0), label="t")
        allf = np.concatenate([v[v > 0] for v in self.f])
        ax.set_ylim(max(allf.min(), allf.max() * 1e-12), allf.max() * 3)
        # x range over every frame (a cold initial spectrum fills one bin and would set it alone)
        seen = np.any(np.asarray(self.f) > 0, axis=0) | ok
        if seen.any():
            ax.set_xlim(self.E[seen].min(), self.E[seen].max())
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel("dN/dE")
        ax.legend(loc="lower left")

    def update(self, t):
        i = self._index(self.t, t)
        self.line.set_ydata(np.where(self.f[i] > 0, self.f[i], np.nan))
        self.line.set_label(f"t = {self.t[i]:.1f}")
        self.ax.legend(loc="lower left")
        self.ax.set_title(f"{self.sp} spectrum")


class HistoryPanel(Panel):
    def __init__(self, run, log=True):
        self.run = run
        self.log = log

    def init(self, ax):
        self.ax = ax
        P.energy_history(ax, self.run, log=self.log)
        self.cursor = ax.axvline(0.0, color=P.INK["secondary"], lw=0.8)
        ax.set_title("energy")

    def update(self, t):
        self.cursor.set_xdata([t, t])


class TracksPanel(Panel):
    """Trajectories of tracked particles in a (qx, qy) plane up to the frame time."""

    def __init__(self, run, diag, species, qx="x", qy="y", n=40, color_by="gamma"):
        tr = run.tracks(diag, species)
        self.t = tr.t
        self.X, self.Y = tr[qx], tr[qy]
        self.qx, self.qy = qx, qy
        g = tr["gamma"]
        # the n particles with the largest final energy (most interesting for acceleration)
        last = np.nan_to_num(np.nanmax(g, axis=0), nan=0.0)
        self.cols = np.argsort(last)[::-1][:n]
        self.sp = species

    def times(self):
        return self.t

    def init(self, ax):
        P.style()
        self.ax = ax
        self.lines = [ax.plot([], [], lw=0.8, color=P.color(k % 8))[0] for k in range(len(self.cols))]
        x, y = self.X[:, self.cols], self.Y[:, self.cols]
        ax.set_xlim(np.nanmin(x), np.nanmax(x))
        ax.set_ylim(np.nanmin(y), np.nanmax(y))
        ax.set_xlabel(self.qx)
        ax.set_ylabel(self.qy)
        ax.set_title(f"{self.sp}: {len(self.cols)} most energetic tracked particles")

    def update(self, t):
        i = self._index(self.t, t)
        for ln, c in zip(self.lines, self.cols):
            ln.set_data(self.X[: i + 1, c], self.Y[: i + 1, c])


class Movie:
    def __init__(self, run, title: str | None = None, ncols: int | None = None, panel_size=(5.0, 3.6)):
        self.run = run
        self.title = title
        self.panels = []
        self.ncols = ncols
        self.panel_size = panel_size

    # builder API ---------------------------------------------------------
    def add(self, panel: Panel):
        self.panels.append(panel)
        return self

    def slice(self, diag, comp, **kw):
        return self.add(SlicePanel(self.run, diag, comp, **kw))

    def cut(self, diag, comp, axis="z", **kw):
        return self.add(FieldCutPanel(self.run, diag, comp, axis, **kw))

    def lineout(self, diag, comps, **kw):
        return self.add(LineoutPanel(self.run, diag, comps, **kw))

    def xt(self, diag, comp, **kw):
        return self.add(XTPanel(self.run, diag, comp, **kw))

    def phase(self, diag, species, **kw):
        return self.add(PhasePanel(self.run, diag, species, **kw))

    def spectrum(self, diag, species, **kw):
        return self.add(SpectrumPanel(self.run, diag, species, **kw))

    def history(self, **kw):
        return self.add(HistoryPanel(self.run, **kw))

    def tracks(self, diag, species, qx="x", qy="y", **kw):
        return self.add(TracksPanel(self.run, diag, species, qx, qy, **kw))

    # rendering -------------------------------------------------------------
    def frame_times(self, n_frames: int | None = None):
        ts = [p.times() for p in self.panels if p.times() is not None]
        if not ts:
            raise ValueError("no time-dependent panel")
        base = max(ts, key=len)  # the most frequent diagnostic sets the frame times
        lo = max(t[0] for t in ts)
        base = base[base >= lo - 1e-9]
        if n_frames and len(base) > n_frames:
            base = base[np.linspace(0, len(base) - 1, n_frames).round().astype(int)]
        return np.asarray(base)

    def save(self, path: str, fps: int = 15, n_frames: int | None = 300, dpi: int = 110) -> str:
        from matplotlib import animation

        P.style()
        times = self.frame_times(n_frames)
        for p in self.panels:
            p.prepare(times)
        n = len(self.panels)
        nc = self.ncols or (2 if n == 4 else min(n, 3))
        nr = math.ceil(n / nc)
        fig, axes = plt.subplots(nr, nc, figsize=(self.panel_size[0] * nc, self.panel_size[1] * nr), squeeze=False)
        for ax in axes.ravel()[n:]:
            ax.set_visible(False)
        for p, ax in zip(self.panels, axes.ravel()):
            p.init(ax)
        sup = fig.suptitle(self.title or "", color=P.INK["primary"])

        def frame(k):
            t = times[k]
            for p in self.panels:
                p.update(t)
            sup.set_text(f"{self.title + '   ' if self.title else ''}ω_pe t = {t:.2f}")
            return []

        frame(0)  # the titles of a frame take part in the layout
        fig.tight_layout()

        anim = animation.FuncAnimation(fig, frame, frames=len(times), blit=False)
        writer, path = _writer(path, fps)
        anim.save(path, writer=writer, dpi=dpi)
        plt.close(fig)
        return path


def _writer(path: str, fps: int):
    from matplotlib import animation

    if path.endswith(".gif"):
        return animation.PillowWriter(fps=fps), path
    try:
        import imageio_ffmpeg

        plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
        return animation.FFMpegWriter(fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p"]), path
    except Exception:  # noqa: BLE001 - no ffmpeg: fall back to an animated GIF
        if animation.FFMpegWriter.isAvailable():
            return animation.FFMpegWriter(fps=fps), path
        gif = os.path.splitext(path)[0] + ".gif"
        warnings.warn(f"ffmpeg not available (pip install imageio-ffmpeg); writing {gif} instead", stacklevel=3)
        return animation.PillowWriter(fps=fps), gif
