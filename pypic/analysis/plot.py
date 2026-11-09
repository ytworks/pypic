"""Matplotlib figure helpers with one consistent style.

Colour by job: categorical series in a fixed order (validated palette: adjacent pairs pass the
colour-vision-deficiency checks), magnitudes (densities, histograms) on a single-hue
sequential map, signed fields (E, B, J) on a blue-gray-red diverging map centred at zero.
Every plot with two or more series has a legend; theory curves are dashed black.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK = {"primary": "#0b0b0b", "secondary": "#52514e", "muted": "#898781", "grid": "#e1e0d9",
       "axis": "#c3c2b7", "surface": "#fcfcfb"}
THEORY = {"color": INK["primary"], "ls": "--", "lw": 1.2}
_SEQ = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SEQUENTIAL = LinearSegmentedColormap.from_list("pypic_sequential", _SEQ)
DIVERGING = LinearSegmentedColormap.from_list(
    "pypic_diverging", ["#0d366b", "#2a78d6", "#9ec5f4", "#f0efec", "#f2a7a6", "#e34948", "#7d1f1e"])


def style():
    """Apply the pypic figure style (idempotent)."""
    from cycler import cycler

    plt.rcParams.update({
        "axes.prop_cycle": cycler(color=SERIES), "lines.linewidth": 1.5, "lines.markersize": 4,
        "axes.edgecolor": INK["axis"], "axes.labelcolor": INK["secondary"], "axes.titlesize": 10,
        "axes.labelsize": 9, "axes.titlecolor": INK["primary"], "xtick.color": INK["muted"],
        "ytick.color": INK["muted"], "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.grid": True, "grid.color": INK["grid"], "grid.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": INK["surface"], "axes.facecolor": INK["surface"], "savefig.facecolor": INK["surface"],
        "legend.frameon": False, "legend.fontsize": 7.5, "font.family": "sans-serif",
        "image.cmap": "pypic_sequential",
    })
    for cm in (SEQUENTIAL, DIVERGING):
        if cm.name not in matplotlib.colormaps:
            matplotlib.colormaps.register(cm)


def color(i: int) -> str:
    """Categorical colour of series i (fixed order; a 9th series should be folded or faceted)."""
    return SERIES[i % len(SERIES)]


def robust_limits(arrays, signed: bool, log: bool = False, pct: float = 99.5):
    """Colour limits shared by several frames: symmetric for signed data, robust percentiles."""
    vals = np.concatenate([np.ravel(a)[np.isfinite(np.ravel(a))] for a in arrays]) if arrays else np.zeros(1)
    if signed:
        m = float(np.percentile(np.abs(vals), pct)) if vals.size else 1.0
        return -m or -1.0, m or 1.0
    if log:
        pos = vals[vals > 0]
        if not pos.size:
            return 1e-10, 1.0
        hi = float(np.percentile(pos, 100.0))
        return max(float(np.percentile(pos, 100 - pct)), hi * 1e-6), hi
    return float(np.percentile(vals, 100 - pct)), float(np.percentile(vals, pct))


def field_map(ax, data2d, extent, signed=True, log=False, vlim=None, label=None, xlabel="", ylabel=""):
    """2-D map of data2d[i, j] with i along the horizontal axis (as pypic arrays are ordered)."""
    style()
    if vlim is None:
        vlim = robust_limits([data2d], signed, log)
    if log:
        norm = LogNorm(vmin=vlim[0], vmax=vlim[1])
        img = np.where(np.asarray(data2d) > 0, data2d, np.nan)
    else:
        norm = Normalize(vmin=vlim[0], vmax=vlim[1])
        img = data2d
    im = ax.imshow(np.asarray(img).T, origin="lower", aspect="auto", extent=extent, norm=norm,
                   cmap=DIVERGING if signed else SEQUENTIAL, interpolation="nearest")
    ax.grid(False)
    cb = ax.figure.colorbar(im, ax=ax, pad=0.01, fraction=0.05)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, colors=INK["muted"])
    if label:
        cb.set_label(label, color=INK["secondary"], fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return im


def energy_history(ax, run, species=True, log=True):
    style()
    s = run.scalars
    t = s.t
    i = 0
    if species:
        for k in s.keys():
            if k.startswith("kinetic_"):
                ax.plot(t, s[k], color=color(i), label=f"K {k[8:]}")
                i += 1
    else:
        ax.plot(t, s["kinetic"], color=color(i), label="kinetic")
        i += 1
    ax.plot(t, s["electric"], color=color(i), label="electric")
    ax.plot(t, s["magnetic"], color=color(i + 1), label="magnetic")
    if log:
        ax.set_yscale("log")
    ax.set_xlabel(r"$\omega_{pe} t$")
    ax.set_ylabel(r"energy [$m_e c^2 n_0 (c/\omega_{pe})^3$]")
    ax.legend()


def conservation(ax, run):
    style()
    s = run.scalars
    t = s.t
    b = s["budget"]
    ax.semilogy(t, np.abs(b - b[0]) / abs(b[0]) + 1e-17, color=color(0), label=r"$|\Delta E_{tot}|/E_0$")
    g, d = s["gauss"], s["divB"]
    fin = np.isfinite(g)
    ax.semilogy(t[fin], g[fin] + 1e-20, "o", ms=3, color=color(1), label=r"max$|\nabla\cdot E-\rho|$")
    ax.semilogy(t[fin], d[fin] + 1e-20, "s", ms=3, color=color(2), label=r"max$|\nabla\cdot B|$")
    ax.set_xlabel(r"$\omega_{pe} t$")
    ax.legend()


def omega_k_map(ax, k, w, P, curves=(), k_max=None, w_max=None, floor=1e-8):
    """log10 power of an omega-k diagram with theory curves [(label, k, w), ...] (dashed)."""
    style()
    Pn = P / max(P.max(), 1e-300)
    ext = [k[0], k[-1], w[0], w[-1]]
    im = ax.imshow(np.log10(np.maximum(Pn, floor)), origin="lower", aspect="auto", extent=ext, cmap=SEQUENTIAL,
                   vmin=np.log10(floor), vmax=0, interpolation="nearest")
    ax.grid(False)
    cb = ax.figure.colorbar(im, ax=ax, pad=0.01, fraction=0.05)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, colors=INK["muted"])
    cb.set_label(r"$\log_{10}$ power (normalised)", color=INK["secondary"], fontsize=8)
    lss = ["--", ":", "-.", (0, (5, 2, 1, 2)), (0, (1, 1))]
    for i, (lab, kc, wc) in enumerate(curves):
        ax.plot(kc, wc, color=INK["primary"], ls=lss[i % len(lss)], lw=1.1, label=lab)
    if curves:
        ax.legend(loc="upper left")
    ax.set_xlim(-(k_max or abs(k).max()), k_max or abs(k).max())
    ax.set_ylim(0, w_max or w.max())
    ax.set_xlabel(r"$k c/\omega_{pe}$")
    ax.set_ylabel(r"$\omega/\omega_{pe}$")


def spectrum_plot(ax, E, f, label=None, i=0, fit=None, xlabel=r"$\gamma - 1$", ylabel=r"$dN/d\gamma$"):
    style()
    ok = f > 0
    ax.loglog(E[ok], f[ok], color=color(i), label=label)
    if fit is not None and np.isfinite(fit.get("p", np.nan)):
        sel = (E >= fit["lo"]) & (E <= fit["hi"])
        ax.loglog(E[sel], fit["norm"] * E[sel] ** (-fit["p"]), **THEORY,
                  label=rf"$\propto E^{{-{fit['p']:.2f}\pm{fit['p_err']:.2f}}}$")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if label or fit is not None:
        ax.legend()


def phase_space(ax, hist2d, edges, log=True, xlabel="", ylabel="", vlim=None):
    ext = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]]
    return field_map(ax, hist2d, ext, signed=False, log=log, vlim=vlim, xlabel=xlabel, ylabel=ylabel,
                     label="weight")


def save(fig, path, dpi=130):
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path
