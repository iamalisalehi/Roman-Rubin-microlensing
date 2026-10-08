#!/usr/bin/env python3
"""One figure style for every plot this project publishes.

Figures sit side by side at a fixed column width in a paper, and the caption, not the image,
carries the words. This module sets that up:

  Vector output. save_figure() writes a PDF for the paper and a PNG beside it for quick viewing.
  Embedded fonts. pdf.fonttype = 42 embeds TrueType so text stays selectable (journals ask for it).
  Real column widths. WIDTH holds the standard single (3.5 in) and double (7.2 in) column measures,
    so text sized here prints at that size.
  No title inside the figure; panels get a short (a)/(b) tag via panel_label().
  One colourblind-safe palette. The survey colours (deep blue, burnt orange, teal) separate under
    deuteranopia and protanopia and stay distinct in greyscale by lightness. Check both before
    adding a fourth.

Usage:

    import plotstyle as ps
    ps.use_paper_style()
    fig, ax = ps.figure(width="single")
    ax.plot(x, y, color=ps.SURVEY["joint"], label=ps.SURVEY_LABEL["joint"])
    ps.legend(ax)
    ps.stamp(fig, provenance_line)      # git commit, weighting, N_eff -- small, grey, bottom
    ps.save_figure(fig, "figures/p5_efficiency_vs_mass")   # writes .pdf and .png
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Palette ----
# Survey partitions.
SURVEY = {"joint": "#0d366b", "roman": "#c2410c", "rubin": "#0e7490"}
SURVEY_LABEL = {"joint": "joint fit", "roman": "Roman alone", "rubin": "Rubin alone"}

# Lens populations. Distinct in hue and lightness; the survey colours are not reused, since a
# figure may show both dimensions.
POPULATION = {"bulge": "#3f3f46", "bh": "#1d4ed8", "ns": "#b45309", "besancon": "#0f766e"}
# Labels are short to fit a 3.5-inch panel; mass ranges and distribution parameters go in the caption.
POPULATION_LABEL = {
    "bulge": "bulge: Kroupa + remnants",
    "bh":    "black holes: log-uniform",
    "ns":    "neutron stars",
    "besancon": "stars: Besancon list",
}

# Galactic components of the stellar catalogue (Besancon Pop codes: thin 1-7, bulge 10, thick 8+11,
# halo 9). All pairs checked on a white surface: CVD dE >= 9.2, normal-vision dE >= 16.3. The aqua is
# below 3:1 contrast on white, so figures using these must label the lines, not rely on colour alone.
COMPONENT = {"thin": "#2a78d6", "bulge": "#eb6834", "thick": "#1baf7a", "halo": "#4a3aa7"}
COMPONENT_LABEL = {"thin": "thin disc", "bulge": "bulge", "thick": "thick disc", "halo": "halo"}

INK = "#1a1a1a"      # axis labels, tick labels, data text
MUTED = "#6b6b6b"    # annotations, secondary text, the provenance stamp
GRID = "#e8e8e8"
SURFACE = "#ffffff"  # white, not off-white: journals composite onto white

# Column widths in inches: ApJ/AAS single column 3.5 in, double 7.2 in (MNRAS is within a few per cent).
WIDTH = {"single": 3.5, "double": 7.2, "wide": 7.2}


def use_paper_style():
    """Set the rcParams. Call once, before creating any figure."""
    plt.rcParams.update({
        # Vector output with embedded, selectable text.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        # DejaVu ships with matplotlib, so this renders identically on any machine.
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Nimbus Roman"],
        "mathtext.fontset": "dejavuserif",
        # Sized for a 3.5-inch column: 8 pt labels, 7 pt ticks.
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "figure.titlesize": 9,
        # Thin frame; ticks inside on both sides (astronomy convention).
        "axes.edgecolor": "#c8c8c8",
        "axes.labelcolor": INK,
        "axes.linewidth": 0.7,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "lines.linewidth": 1.4,
        "legend.frameon": False,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


def figure(width="single", height=None, nrows=1, ncols=1, **kw):
    """A figure at a real column width. `height` defaults to a 4:3 panel."""
    w = WIDTH.get(width, width if isinstance(width, (int, float)) else 3.5)
    if height is None:
        height = 0.75 * w * nrows / ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(w, height), **kw)
    return fig, axes


def panel_label(ax, text, loc="upper left"):
    """A short (a)/(b) tag inside the axes -- the paper equivalent of a title.

    Placed inside the axes so it adds no vertical space between stacked panels.
    """
    x, y, ha, va = (0.03, 0.97, "left", "top")
    if loc == "upper right":
        x, ha = 0.97, "right"
    elif loc == "lower left":
        y, va = 0.03, "bottom"
    elif loc == "lower right":
        x, y, ha, va = 0.97, 0.03, "right", "bottom"
    ax.text(x, y, text, transform=ax.transAxes, color=INK, fontsize=8,
            ha=ha, va=va, fontweight="bold")


def legend(ax, **kw):
    """Frameless legend with the project's text colour, or nothing if there is nothing to list.

    A panel can legitimately have no labelled series (every bin too thin to draw), and matplotlib
    would warn about missing artists.
    """
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return None
    kw.setdefault("frameon", False)
    kw.setdefault("handlelength", 1.6)
    kw.setdefault("borderpad", 0.2)
    kw.setdefault("labelcolor", INK)
    return ax.legend(**kw)


def stamp(fig, text):
    """The provenance line: commit, population, weighting, N_eff.

    Small and grey at the bottom of the figure, so weighting and population can be read off it.
    Strip it only for camera-ready submission, where the same facts belong in the caption.
    """
    # wrap=True: with bbox "tight" an unwrapped long stamp would widen the saved image.
    fig.text(0.0, -0.015, text, color=MUTED, fontsize=5.5, ha="left", va="top", wrap=True)


def plain_log_ticks(ax, lo, hi, axis="y"):
    """Plain numbers on a log axis that spans only a decade or two.

    Matplotlib labels log minor ticks on short ranges, which collide at column width, and
    suppressing them can leave a sub-decade axis with no numbers. Explicit ticks with a plain
    formatter avoid both. Above ~2.2 decades the default decade ticks are fine and this does nothing.

    The range is passed in from the data because the axes have not autoscaled yet when a figure
    function calls this, so get_ylim() would return provisional limits.
    """
    import numpy as np
    import matplotlib.ticker as mt
    a = ax.yaxis if axis == "y" else ax.xaxis
    if not (lo > 0 and hi > lo) or np.log10(hi / lo) > 2.2:
        return
    ticks = np.geomspace(lo, hi, 5)
    a.set_major_locator(mt.FixedLocator(ticks))
    a.set_minor_locator(mt.NullLocator())
    a.set_major_formatter(mt.FuncFormatter(
        lambda v, _: f"{v:.2f}".rstrip("0").rstrip(".") if v < 10 else f"{v:.0f}"))


def save_figure(fig, path_without_extension, dpi=300, formats=("pdf", "png")):
    """Write the figure once per format. Returns the paths written."""
    written = []
    for ext in formats:
        p = f"{path_without_extension}.{ext}"
        fig.savefig(p, dpi=dpi)
        written.append(p)
    plt.close(fig)
    return written
