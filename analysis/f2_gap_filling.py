#!/usr/bin/env python3
"""Step F2 -- the gap-filling figure.

THE result this project exists to produce. Roman observes the bulge in ~70-day seasons
separated by ~110-day gaps forced by the Sun angle. An event peaking in a gap is one Roman
sees poorly or not at all, while Rubin -- which observes the bulge all year at low cadence --
keeps taking data straight through. The claim advocated in white papers for years, and never
simulated, is that Rubin's coverage rescues those events. This figure tests it.

What is plotted
---------------
x: dt_edge, days from t0 to the nearest Roman season boundary. NEGATIVE inside a season,
   positive in a gap, so the plot reads left-to-right as "deeper in season" -> "deeper in gap".

Left panel  -- PRECISION. Median sigma_joint / sigma_Roman per event, split by tE bin.
               Expected flat at ~1 mid-season and dropping into the gap, the drop deepening
               with tE. Only events Roman characterises ALONE can appear here, because the
               ratio needs a denominator.

Right panel -- YIELD. Fraction of joint-detected events that the joint fit characterises but
               Roman alone does not. This is where the gap-filling shows up once Roman fails
               outright: sigma_Roman is undefined, the ratio vanishes, and a precision-only
               plot would silently drop exactly the events that make the case.

Both panels are needed. The plan (Phase F) is explicit that short-tE gain is a yield statistic
and long-tE gain is a precision statistic, and that forcing one metric across both regimes
produces a misleadingly flat answer.

Scope -- two restrictions, both load-bearing
--------------------------------------------
1. Events peaking WITHIN Roman's mission (t0zone in-season or in-gap). Off-mission events --
   t0 before Roman's first epoch or after its last -- are Rubin-only by construction and are
   not a gap-filling result. That three-way distinction is why T0Zone exists rather than a
   simple in/out flag.

2. Events Roman ACTUALLY OBSERVED (ndw_R > 0). This one is easy to get wrong and it inverts
   the answer. Only ~39 of 1706 sightlines fall inside Roman's footprint, so the large
   majority of joint-detected events have no Roman epochs at all. Counting those as "Roman
   alone could not characterise it" is true but vacuous -- Roman missed them because it never
   pointed there, not because of a season gap. Including them buried the real signal under a
   flat ~20% floor that was really just footprint coverage. The gap-filling question is only
   meaningful where Roman has data and the timing is what limits it.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R

# Ordinal ramp, one hue light->dark, for the ORDERED tE bins. Validated with
# validate_palette.js --ordinal: monotone lightness, all gaps >= 0.06, light end 2.06:1
# on the surface. An ordered quantity must not get a categorical rainbow.
TE_COLORS = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
TE_EDGES = [10.0, 30.0, 100.0, 300.0, np.inf]
TE_LABELS = ["10-30 d", "30-100 d", "100-300 d", "> 300 d"]

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e6e5e1"
MIN_PER_BIN = 8   # below this a median and its quartiles are not worth drawing


def bin_edges(dt_max, width):
    lo = np.floor(-40.0 / width) * width
    return np.arange(lo, dt_max + width, width)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("events", help="per-event table (test2.dat)")
    ap.add_argument("-o", "--out", default="f2_gap_filling.png")
    ap.add_argument("--provenance", default="files/MONTLMC/files/run_provenance.txt")
    ap.add_argument("--dt-max", type=float, default=60.0,
                    help="largest dt_edge to plot [d]; dt_edge is distance to the NEAREST "
                         "edge, so a ~110 d gap tops out near 55 d")
    ap.add_argument("--width", type=float, default=15.0, help="dt_edge bin width [d]")
    a = ap.parse_args()

    df = R.load_events(a.events)

    # See the module docstring: both restrictions matter, and dropping the second one
    # (ndw_R > 0) silently turns this figure into a plot of Roman's footprint.
    d = df[R.detected(df, "joint")
           & df["t0zone"].isin([0, 1])
           & (df["ndw_R"] > 0)].copy()
    d = d[d["dt_edge"] <= a.dt_max]
    if d.empty:
        sys.exit("no in-mission joint-detected events in range -- nothing to plot")

    d["ratio_tE"] = R.ratio_joint_over(d, "tE", "roman")
    d["char_joint"] = R.characterized(d, "joint")
    d["char_roman"] = R.characterized(d, "roman")
    d["rescued"] = d["char_joint"] & ~d["char_roman"]
    d["teBin"] = pd.cut(d["tE"], TE_EDGES, labels=TE_LABELS, right=False)

    edges = bin_edges(a.dt_max, a.width)
    d["dtBin"] = pd.cut(d["dt_edge"], edges)
    centres = {iv: iv.mid for iv in d["dtBin"].cat.categories}

    fig, (axP, axY) = plt.subplots(1, 2, figsize=(12.5, 5.0), constrained_layout=True)
    rows = []

    for label, colour in zip(TE_LABELS, TE_COLORS):
        sub = d[d["teBin"] == label]
        if sub.empty:
            continue

        # ---- precision panel: median ratio, IQR band ----
        xs, med, q25, q75 = [], [], [], []
        for iv, g in sub.groupby("dtBin", observed=True):
            r = g["ratio_tE"].dropna()
            if len(r) < MIN_PER_BIN:
                continue
            xs.append(centres[iv]); med.append(r.median())
            q25.append(r.quantile(0.25)); q75.append(r.quantile(0.75))
            rows.append(dict(panel="precision", teBin=label, dt_centre=centres[iv],
                             n=len(r), median=r.median(),
                             q25=r.quantile(0.25), q75=r.quantile(0.75)))
        if xs:
            # Quartiles as thin error bars, not a filled band: with three series the bands
            # overlapped so heavily they hid the medians they were supposed to qualify.
            lo = np.array(med) - np.array(q25)
            hi = np.array(q75) - np.array(med)
            axP.errorbar(xs, med, yerr=[lo, hi], color=colour, linewidth=2.0,
                         elinewidth=1.0, capsize=3, marker="o", markersize=5,
                         markeredgecolor="#fcfcfb", markeredgewidth=1.0, label=label)

        # ---- yield panel: rescue fraction ----
        xs2, frac = [], []
        for iv, g in sub.groupby("dtBin", observed=True):
            if len(g) < MIN_PER_BIN:
                continue
            xs2.append(centres[iv]); frac.append(g["rescued"].mean())
            rows.append(dict(panel="yield", teBin=label, dt_centre=centres[iv],
                             n=len(g), rescue_fraction=g["rescued"].mean()))
        if xs2:
            axY.plot(xs2, frac, color=colour, linewidth=2.0, marker="o", markersize=5,
                     markeredgecolor="#fcfcfb", markeredgewidth=1.0, label=label)

    for ax in (axP, axY):
        ax.axvspan(edges[0], 0.0, color="#f2f1ed", zorder=0)   # inside a Roman season
        ax.axvline(0.0, color=INK2, linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=INK2, labelsize=9)
        ax.set_xlabel("days from t$_0$ to nearest Roman season edge\n"
                      "$\\leftarrow$ inside season      in gap $\\rightarrow$",
                      color=INK2, fontsize=10)

    axP.axhline(1.0, color=INK2, linewidth=1.0, alpha=0.5)
    axP.set_ylabel(r"median  $\sigma_{\rm joint}(t_E)\ /\ \sigma_{\rm Roman}(t_E)$",
                   color=INK, fontsize=10)
    axP.set_title("Precision: what the joint fit adds where Roman still measures",
                  color=INK, fontsize=11, loc="left")

    axY.set_ylabel("fraction characterised jointly but NOT by Roman alone",
                   color=INK, fontsize=10)
    axY.set_title("Yield: events Roman alone cannot characterise at all",
                  color=INK, fontsize=11, loc="left")
    axY.set_ylim(0, 1)

    leg = axP.legend(title="$t_E$", frameon=False, fontsize=9, title_fontsize=9,
                     loc="lower left", bbox_to_anchor=(0.0, 0.0))
    leg.get_title().set_color(INK2)
    for t in leg.get_texts():
        t.set_color(INK2)

    stamp = R.describe(a.events, a.provenance)
    fig.suptitle("Roman season gaps: what Rubin's year-round coverage recovers",
                 color=INK, fontsize=13, x=0.005, ha="left")
    fig.text(0.005, -0.02, stamp, color=INK2, fontsize=7.5, ha="left")

    fig.savefig(a.out, dpi=160, bbox_inches="tight", facecolor="#fcfcfb")
    print(f"wrote {a.out}")

    if rows:
        csv = os.path.splitext(a.out)[0] + ".csv"
        pd.DataFrame(rows).to_csv(csv, index=False)
        print(f"wrote {csv}")
        print(f"\nevents in scope: {len(d)}  "
              f"(ratio defined on {int(d['ratio_tE'].notna().sum())}, "
              f"rescued {int(d['rescued'].sum())})")
        print(f"characterised: joint {int(d['char_joint'].sum())}, "
              f"Roman alone {int(d['char_roman'].sum())}")


if __name__ == "__main__":
    sys.exit(main())
