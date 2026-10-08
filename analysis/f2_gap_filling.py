#!/usr/bin/env python3
"""The gap-filling figure: what Rubin's year-round coverage recovers in Roman's season gaps.

Roman observes the bulge in ~70-day seasons separated by ~110-day gaps forced by the Sun angle. An
event peaking in a gap is seen poorly or not at all by Roman, while Rubin keeps taking data at low
cadence through the gap.

x: dt_edge, days from t0 to the nearest Roman season boundary. Negative inside a season, positive in
   a gap.

Left panel  -- PRECISION. Median sigma_joint / sigma_Roman per event, split by tE bin, for one fitted
               parameter chosen with --param. Only events Roman characterises alone appear, because
               the ratio needs a denominator.
               --param tE (default): the drop deepens with short tE, because a 20-day event peaking
               in a 110-day gap is missed by Roman entirely while a 500-day one is still magnified
               when the next season opens.
               --param piE: parallax is measured from the distortion Earth's orbital motion
               imprints on the light curve, so it needs sampling across a substantial fraction of a
               year, which only long events provide and which Roman's gaps interrupt.

Right panel -- YIELD. Fraction of joint-detected events that the joint fit characterises but Roman
               alone does not. Independent of --param (characterisation is the two-parameter
               criterion tE > 2 sigma_tE AND piE > 2 sigma_piE). Where Roman fails outright the
               ratio is undefined and the precision panel drops those events; this panel keeps them.

Scope (both restrictions are load-bearing):
1. Events peaking within Roman's mission (t0zone in-season or in-gap). Off-mission events are
   Rubin-only by construction.
2. Events Roman actually observed (ndw_R > 0). Most sightlines lie outside Roman's footprint, and
   counting those as "Roman could not characterise it" turns the figure into a plot of footprint
   coverage rather than of season gaps.
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

# Ordinal ramp (one hue, light to dark) for the ordered tE bins.
TE_COLORS = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
TE_EDGES = [10.0, 30.0, 100.0, 300.0, np.inf]
TE_LABELS = ["10-30 d", "30-100 d", "100-300 d", "> 300 d"]

# Symbol used for each fittable parameter in the panel title (piE dimensionless, tE in days).
PARAM_TEX = {"tE": "t_E", "piE": r"\pi_E"}

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
    ap.add_argument("--param", choices=("tE", "piE"), default="tE",
                    help="which forecast sigma the precision panel plots the ratio of. "
                         "tE = Einstein crossing time [d]; piE = microlensing parallax "
                         "[dimensionless], the parameter the plan's long-tE argument is about")
    ap.add_argument("--provenance", default="files/MONTLMC/files/run_provenance.txt")
    ap.add_argument("--dt-max", type=float, default=60.0,
                    help="largest dt_edge to plot [d]; dt_edge is distance to the NEAREST "
                         "edge, so a ~110 d gap tops out near 55 d")
    ap.add_argument("--width", type=float, default=15.0, help="dt_edge bin width [d]")
    ap.add_argument("--map", default=None,
                    help="MapLMC5.dat -- needed for the event-rate weight")
    ap.add_argument("--log", action="append", default=[],
                    help="run log(s), for sightlines whose map rows a killed run lost")
    ap.add_argument("--chunksize", type=int, default=500_000,
                    help="rows per read chunk; holds peak memory to one chunk plus survivors")
    ap.add_argument("--unweighted", action="store_true",
                    help="deliberately report the raw sample, with no event-rate weight")
    a = ap.parse_args()

    # Stream the table, dropping rows of barren sightlines (no nsim, weight 0) as they are read.
    # Undetected events are kept: the denominator is every event in the tE bin.
    df = R.load_events(a.events, keep=R.keep_weightable(a.map, a.log), chunksize=a.chunksize)
    w, wlabel = R.attach_weight(df, a.map, a.log, a.unweighted)
    df["W"] = w
    print(f"weighting: {wlabel}")

    # Scope restrictions: see the module docstring.
    d = df[R.detected(df, "joint")
           & df["t0zone"].isin([0, 1])
           & (df["ndw_R"] > 0)].copy()
    d = d[d["dt_edge"] <= a.dt_max]
    if d.empty:
        sys.exit("no in-mission joint-detected events in range -- nothing to plot")

    d["ratio"] = R.ratio_joint_over(d, a.param, "roman")
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
        xs, med, q25, q75, med_raw = [], [], [], [], []
        for iv, g in sub.groupby("dtBin", observed=True):
            ok = g["ratio"].notna()
            r = g.loc[ok, "ratio"].to_numpy()
            rw = g.loc[ok, "W"].to_numpy()
            # Gate on the effective sample, not the raw count.
            neff = R.kish_neff(rw)
            if len(r) < MIN_PER_BIN or neff < MIN_PER_BIN:
                continue
            xs.append(centres[iv]); med.append(R.weighted_median(r, rw))
            q25.append(R.weighted_quantile(r, rw, 0.25))
            q75.append(R.weighted_quantile(r, rw, 0.75))
            med_raw.append(float(np.median(r)))
            rows.append(dict(panel="precision", param=a.param, teBin=label,
                             dt_centre=centres[iv],
                             n=len(r), n_eff=neff, median=R.weighted_median(r, rw),
                             median_unweighted=float(np.median(r)),
                             q25=R.weighted_quantile(r, rw, 0.25),
                             q75=R.weighted_quantile(r, rw, 0.75)))
        if xs:
            axP.plot(xs, med_raw, color=colour, lw=0.9, ls=(0, (3, 2)), alpha=0.55, zorder=2)
            # Quartiles as error bars rather than a band, so overlapping series stay readable.
            lo = np.array(med) - np.array(q25)
            hi = np.array(q75) - np.array(med)
            axP.errorbar(xs, med, yerr=[lo, hi], color=colour, linewidth=2.0,
                         elinewidth=1.0, capsize=3, marker="o", markersize=5,
                         markeredgecolor="#fcfcfb", markeredgewidth=1.0, label=label)

        # ---- yield panel: rescue fraction ----
        xs2, frac = [], []
        for iv, g in sub.groupby("dtBin", observed=True):
            gw = g["W"].to_numpy()
            neff = R.kish_neff(gw)
            if len(g) < MIN_PER_BIN or neff < MIN_PER_BIN:
                continue
            xs2.append(centres[iv])
            frac.append(R.weighted_fraction(g["rescued"], gw))
            rows.append(dict(panel="yield", param=a.param, teBin=label,
                             dt_centre=centres[iv],
                             n=len(g), n_eff=neff,
                             rescue_fraction=R.weighted_fraction(g["rescued"], gw),
                             rescue_fraction_unweighted=float(g["rescued"].mean())))
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
    axP.set_ylabel(f"median  $\\sigma_{{\\rm joint}}({PARAM_TEX[a.param]})"
                   f"\\ /\\ \\sigma_{{\\rm Roman}}({PARAM_TEX[a.param]})$",
                   color=INK, fontsize=10)
    axP.set_title(f"Precision in ${PARAM_TEX[a.param]}$: "
                  "what the joint fit adds where Roman still measures",
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

    stamp = (f"param={a.param}  solid: {wlabel}; dashed: raw sample  "
             + R.describe(a.events, a.provenance))
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
              f"(ratio defined on {int(d['ratio'].notna().sum())}, "
              f"rescued {int(d['rescued'].sum())})")
        print(f"characterised: joint {int(d['char_joint'].sum())}, "
              f"Roman alone {int(d['char_roman'].sum())}")


if __name__ == "__main__":
    sys.exit(main())
