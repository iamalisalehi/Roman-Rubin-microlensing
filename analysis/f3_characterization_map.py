#!/usr/bin/env python3
"""Step F3 -- the (tE, piE) characterization map.

Deliberately in the format of Abrams et al. 2025, Figures 11-14: the log(tE) - log(piE)
plane, cells coloured by a RATIO of characterized fractions. The community already reads
that plot; matching the format is worth the effort.

Why this plane
--------------
tE (the Einstein-crossing time, days) and piE (the microlensing parallax, dimensionless)
are the two parameters a photometric light curve can actually deliver, and together they
pin the lens down: with the angular Einstein radius tetE from astrometry,

    Ml = tetE / (kappa * piE),     kappa = 8.144 mas / Msun

and even without astrometry the pair (tE, piE) separates the populations. Since
tE ~ sqrt(Ml) and piE ~ 1/sqrt(Ml), heavy lenses -- black holes -- sit toward LONG tE and
SMALL piE (lower right), and low-mass lenses toward SHORT tE and LARGE piE (upper left).
So the plane is a mass axis running along the diagonal, which is why the dashed iso-mass
lines are drawn: at the sample's median relative proper motion they are exactly lines of
constant lens mass.

What is coloured
----------------
Per cell, the fraction of events meeting the Abrams characterization criterion
(tE > 2 sigma_tE AND piE > 2 sigma_piE) with the joint fit, divided by the same fraction
with one survey alone. The ratio is >= 1 by construction: the joint Fisher matrix is the
sum of the per-survey ones, so sigma_joint <= sigma_single, so any event a single survey
characterizes the joint fit characterizes too. A cell below 1 is a bug, not a result.

Two panels, two different questions
-----------------------------------
(a) joint / Roman-alone, restricted to events Roman actually observed (ndw_R > 0). This is
    "what does adding Rubin buy Roman". Only ~39 of 1706 sightlines are inside Roman's
    footprint, so this is the small, expensive sample -- but it is the only one where the
    Roman-alone denominator means anything. Outside the footprint Roman characterizes
    nothing because it never pointed there, which is true and vacuous.

(b) joint / Rubin-alone, over every joint-detected event. This is "what does adding Roman
    buy Rubin", and it is the whole-survey statistic: Rubin sees the entire bulge region.

Cells where the single survey characterizes NOTHING but the joint fit characterizes
something are pure rescue -- the ratio is infinite, not large. They are hatched rather than
coloured, because painting infinity at the top of a colour ramp would understate them.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R

# Sequential single hue, light -> dark, for a magnitude (the ratio). Same blue family as
# the F2 tE ramp so the two figures read as one set. A ratio is an ordered quantity: it
# gets one hue, never a rainbow.
RAMP = LinearSegmentedColormap.from_list(
    "romanblue", ["#eef4fc", "#c3d9f5", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])

GREY_SPARSE = "#e4e4e4"   # cell with too few events to quote a fraction
INK = "#1a1a1a"
MUTED = "#6b6b6b"
SURFACE = "#ffffff"

KAPPA = 8.144             # mas / Msun


def cell_stats(df, xedges, yedges, alone, min_n):
    """Per-cell N, joint fraction, single-survey fraction, and their ratio."""
    x = np.log10(df["tE"].to_numpy())
    y = np.log10(df["piE"].to_numpy())
    cj = R.characterized(df, "joint").to_numpy()
    cs = R.characterized(df, alone).to_numpy()

    n, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])
    nj, _, _ = np.histogram2d(x, y, bins=[xedges, yedges], weights=cj.astype(float))
    ns, _, _ = np.histogram2d(x, y, bins=[xedges, yedges], weights=cs.astype(float))

    with np.errstate(invalid="ignore", divide="ignore"):
        fj = np.where(n >= min_n, nj / n, np.nan)
        fs = np.where(n >= min_n, ns / n, np.nan)
        ratio = np.where(fs > 0, fj / fs, np.nan)

    # Pure rescue: the single survey characterizes nothing here, the joint fit does.
    rescue = (n >= min_n) & (ns == 0) & (nj > 0)
    return n, fj, fs, ratio, rescue


def draw_panel(ax, xedges, yedges, n, ratio, rescue, vmax, title, subtitle,
               murel, masses):
    ax.set_facecolor(SURFACE)

    # Populated cells that get no colour (too few events, or no single-survey baseline)
    # are painted grey first, so the reader can tell "no data" from "no gain".
    sparse = (n > 0) & np.isnan(ratio) & ~rescue
    ax.pcolormesh(xedges, yedges, np.where(sparse, 1.0, np.nan).T,
                  cmap=LinearSegmentedColormap.from_list("g", [GREY_SPARSE, GREY_SPARSE]),
                  shading="flat", zorder=1)

    mesh = ax.pcolormesh(xedges, yedges, np.ma.masked_invalid(ratio).T, cmap=RAMP,
                         norm=Normalize(vmin=1.0, vmax=vmax), shading="flat", zorder=2)

    # Rescue cells: hatched, not coloured. Infinity is not a large number.
    ry, rx = np.where(rescue.T)
    for iy, ix in zip(ry, rx):
        ax.add_patch(plt.Rectangle((xedges[ix], yedges[iy]),
                                   xedges[ix + 1] - xedges[ix],
                                   yedges[iy + 1] - yedges[iy],
                                   facecolor="none", edgecolor="#0d366b",
                                   hatch="///", linewidth=0.9, zorder=3))

    # Iso-mass diagonals. theta_E = murel * tE / 365.25, Ml = theta_E / (kappa * piE),
    # so log10 piE = log10 tE + log10(murel / (365.25 * kappa * Ml)) -- slope +1.
    x0, x1 = xedges[0], xedges[-1]
    y0, y1 = yedges[0], yedges[-1]
    for Ml in masses:
        c = np.log10(murel / (365.25 * KAPPA * Ml))
        ax.plot([x0, x1], [x0 + c, x1 + c], ls="--", lw=0.9, color=MUTED, zorder=4)
        # Anchor the label to the line's exit point, pulled inside both limits so it
        # never lands in the margin.
        xl = min(x1, y1 - c) - 0.22
        yl = xl + c
        if x0 < xl < x1 and y0 < yl < y1:
            ax.text(xl, yl + 0.04, f"{Ml:g} " + r"M$_\odot$", color=MUTED, fontsize=7.5,
                    ha="center", va="bottom", rotation=45, rotation_mode="anchor",
                    zorder=5)

    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_xlabel(r"$\log_{10}\,t_{\rm E}$  [days]", color=INK, fontsize=10)
    ax.set_ylabel(r"$\log_{10}\,\pi_{\rm E}$", color=INK, fontsize=10)
    # Title and subtitle in one text object: a separate transAxes label collides with
    # the neighbouring panel's title at this figure width.
    ax.set_title(f"{title}\n{subtitle}", color=INK, fontsize=10.5, loc="left", pad=8,
                 linespacing=1.6)
    ax.tick_params(colors=MUTED, labelsize=8)
    for sp in ax.spines.values():
        sp.set_color("#d0d0d0")
    ax.grid(True, color="#f0f0f0", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    return mesh


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("events", help="per-event table (test5.dat)")
    ap.add_argument("-o", "--out", default="f3_characterization_map.png")
    ap.add_argument("--provenance", default=None,
                    help="run_provenance.txt; found automatically if omitted")
    ap.add_argument("--min-n", type=int, default=12,
                    help="cells with fewer events are greyed, not coloured (default 12)")
    ap.add_argument("--dex", type=float, default=0.5, help="cell size in dex (default 0.5)")
    ap.add_argument("--te-range", type=float, nargs=2, default=(-0.5, 3.5))
    ap.add_argument("--pie-range", type=float, nargs=2, default=(-2.5, 1.0))
    args = ap.parse_args()

    df = R.load_events(args.events)
    print(R.describe(args.events, args.provenance))

    bad = R.check_monotonicity(df)
    if bad:
        print("WARNING: sigma_joint > sigma_single somewhere -- see check_monotonicity:")
        for b in bad:
            print("   ", b)

    xedges = np.arange(args.te_range[0], args.te_range[1] + 1e-9, args.dex)
    yedges = np.arange(args.pie_range[0], args.pie_range[1] + 1e-9, args.dex)

    det = df[R.detected(df, "joint")]
    panels = [
        ("roman", det[det["ndw_R"] > 0],
         "(a)  What Rubin adds to Roman",
         "joint / Roman-alone characterized fraction, Roman footprint only"),
        ("rubin", det,
         "(b)  What Roman adds to Rubin",
         "joint / Rubin-alone characterized fraction, all joint-detected events"),
    ]

    murel = float(np.median(det["murel_yr"]))
    masses = [0.1, 1.0, 10.0]

    # Compute BOTH panels before drawing either: the two share one colourbar, so they
    # must share one colour scale. A per-panel vmax under a single bar would make the
    # same shade mean two different numbers.
    stats = []
    for alone, sub, title, subtitle in panels:
        stats.append(cell_stats(sub, xedges, yedges, alone, args.min_n))
    allr = np.concatenate([st[3][np.isfinite(st[3])] for st in stats])
    vmax = float(np.ceil(np.nanmax(allr) * 10.0) / 10.0) if allr.size else 1.5
    vmax = max(vmax, 1.1)
    print(f"  colour scale 1.0 -> {vmax:.2f} (shared by both panels)")

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.9), facecolor=SURFACE)
    rows = []
    meshes = []
    for ax, (alone, sub, title, subtitle), (n, fj, fs, ratio, rescue) in zip(
            axes, panels, stats):
        meshes.append(draw_panel(ax, xedges, yedges, n, ratio, rescue, vmax,
                                 title, subtitle, murel, masses))
        print(f"  panel {alone}: {len(sub)} events, {int(np.nansum(n))} in range, "
              f"{int(np.isfinite(ratio).sum())} cells coloured, "
              f"{int(rescue.sum())} rescue cells, "
              f"max ratio {np.nanmax(ratio) if np.isfinite(ratio).any() else float('nan'):.2f}")
        for i in range(len(xedges) - 1):
            for j in range(len(yedges) - 1):
                if n[i, j] == 0:
                    continue
                rows.append(dict(panel=alone,
                                 log_tE_lo=round(xedges[i], 3),
                                 log_piE_lo=round(yedges[j], 3),
                                 n=int(n[i, j]),
                                 frac_char_joint=fj[i, j],
                                 frac_char_alone=fs[i, j],
                                 ratio=ratio[i, j],
                                 rescue=bool(rescue[i, j])))

    fig.subplots_adjust(left=0.065, right=0.885, bottom=0.235, top=0.87, wspace=0.22)
    cax = fig.add_axes([0.905, 0.235, 0.018, 0.635])
    cbar = fig.colorbar(meshes[-1], cax=cax)
    cbar.set_label("characterized fraction, joint / single survey", color=INK, fontsize=9)
    cbar.ax.tick_params(colors=MUTED, labelsize=8)
    cbar.outline.set_color("#d0d0d0")

    fig.text(0.065, 0.125,
             f"grey = fewer than {args.min_n} events in the cell, or the single survey "
             f"characterizes none and neither does the joint fit\n"
             f"hatched = the single survey characterizes none but the joint fit does "
             f"(ratio infinite, not large)     "
             f"dashed = constant lens mass at the sample median "
             + r"$\mu_{\rm rel}$ = " + f"{murel:.2f} mas/yr\n"
             + R.describe(args.events, args.provenance),
             color=MUTED, fontsize=7.5, va="top", linespacing=1.7)

    fig.savefig(args.out, dpi=200, facecolor=SURFACE)
    print(f"wrote {args.out}")

    import pandas as pd
    csv = os.path.splitext(args.out)[0] + ".csv"
    pd.DataFrame(rows).to_csv(csv, index=False, float_format="%.4f")
    print(f"wrote {csv}")


if __name__ == "__main__":
    main()
