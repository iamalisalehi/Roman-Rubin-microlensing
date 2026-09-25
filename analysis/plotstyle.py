#!/usr/bin/env python3
"""One figure style for every plot this project publishes (Step P5).

WHY THIS EXISTS. The figure scripts grew their own styles: three different ink/grid palettes,
two different survey colour sets, DPIs from 150 to 200, figure sizes chosen per script, PNG
only, and the explanatory text baked into the image. That is fine for a working plot and wrong
for a paper, where figures sit side by side at a fixed column width and the caption -- not the
image -- carries the words.

WHAT A PUBLICATION FIGURE NEEDS, AND WHAT THIS DOES ABOUT IT.

  Vector, not raster. A PDF stays sharp at any zoom and prints properly; a 200-dpi PNG does
  not. save_figure() writes BOTH: the PDF for the paper, a PNG beside it for quick viewing.

  Embedded, editable fonts. pdf.fonttype = 42 embeds TrueType rather than converting glyphs
  to outlines, so the text stays selectable and searchable. Journals ask for this.

  Real column widths. A figure drawn at 16 inches and shrunk to 3.5 has 4-point labels. The
  widths here are the standard single (3.5 in) and double (7.2 in) column measures, so text
  sized here is the size it prints at.

  No title inside the figure. set_title() is for working plots; in a paper the caption does
  that job and a baked-in title is duplicated text you cannot edit at proof stage. Panels get
  a short (a)/(b) tag instead, via panel_label().

  One palette, colourblind-safe. The survey colours are the ones F4 already used -- deep blue,
  burnt orange, teal -- which separate under deuteranopia and protanopia, the two common forms,
  and stay distinct in greyscale by lightness. Do not add a fourth without checking both.

USAGE

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

# ---------------------------------------------------------------------------------------
# Palette
# ---------------------------------------------------------------------------------------
# Survey partitions. Same three colours F4 has used since Phase F, kept so the new population
# figures read as part of the same set as the existing ones.
SURVEY = {"joint": "#0d366b", "roman": "#c2410c", "rubin": "#0e7490"}
SURVEY_LABEL = {"joint": "joint fit", "roman": "Roman alone", "rubin": "Rubin alone"}

# Lens populations. Distinct in hue AND in lightness, so a greyscale print still separates
# them; deliberately not reusing the survey colours, since a figure may show both dimensions.
POPULATION = {"bulge": "#3f3f46", "bh": "#1d4ed8", "ns": "#b45309"}
# Labels are kept short because they sit inside a 3.5-inch panel; the mass ranges and the
# distribution parameters belong in the caption, where there is room for them.
POPULATION_LABEL = {
    "bulge": "bulge: Kroupa + remnants",
    "bh":    "black holes: log-uniform",
    "ns":    "neutron stars",
}

INK = "#1a1a1a"      # axis labels, tick labels, data text
MUTED = "#6b6b6b"    # annotations, secondary text, the provenance stamp
GRID = "#e8e8e8"
SURFACE = "#ffffff"  # white, not off-white: journals composite onto white

# Column widths in inches. ApJ/AAS single column is 3.5 in (246 pt), double 7.2 in (513 pt);
# MNRAS is within a few per cent of both, so one set serves.
WIDTH = {"single": 3.5, "double": 7.2, "wide": 7.2}


def use_paper_style():
    """Set the rcParams. Call once, before creating any figure."""
    plt.rcParams.update({
        # Vector output with embedded, selectable text.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        # A serif face pairs with most journal body text; DejaVu ships with matplotlib, so
        # this renders identically on any machine rather than silently falling back.
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Nimbus Roman"],
        "mathtext.fontset": "dejavuserif",
        # Sized for a 3.5-inch column: 8 pt labels, 7 pt ticks. Text set here is text as
        # printed, which is the whole point of drawing at the real width.
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "figure.titlesize": 9,
        # Thin, unobtrusive frame; ticks inside, both sides, the convention in astronomy.
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

    Inside rather than above, so it survives the tight bounding box and does not add vertical
    space between stacked panels.
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

    The guard is not cosmetic. A panel legitimately ends up with no labelled series whenever
    every bin was too thin to draw -- a partial run, or a population whose events are rare in
    that panel's cut -- and matplotlib answers that with a UserWarning about missing artists.
    A warning that fires on normal sparse data trains the reader to ignore warnings.
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

    Small and grey at the bottom of the figure. It is NOT decoration -- a figure whose
    weighting or population cannot be read off it is a figure nobody can check. Strip it only
    for a camera-ready submission, where the same facts belong in the caption.
    """
    # wrap=True: a long stamp must wrap at the figure width, not widen the saved image --
    # with bbox "tight" an unwrapped line made p7_precision twice as wide as its panels.
    fig.text(0.0, -0.015, text, color=MUTED, fontsize=5.5, ha="left", va="top", wrap=True)


def plain_log_ticks(ax, lo, hi, axis="y"):
    """Plain numbers on a log axis that spans only a decade or two.

    Matplotlib labels log MINOR ticks on short ranges, which at column width collides into
    mush ("6x10^0 4x10^0 3x10^0 ..."); but simply suppressing the minor labels can leave a
    sub-decade axis with no numbers at all. Explicit ticks with a plain formatter is the only
    option that avoids both. Above ~2.2 decades the default decade ticks are fine and this
    does nothing.

    The range is passed IN, from the data, rather than read off the axes: at the point a
    figure function calls this, matplotlib has not autoscaled yet, so get_ylim() returns
    provisional limits and the span test silently takes the wrong branch. That mistake makes
    this function look correct while changing nothing.
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
