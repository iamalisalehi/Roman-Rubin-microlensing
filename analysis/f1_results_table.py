#!/usr/bin/env python3
"""Step F1 -- the results table generator.

Produces, per (Roman field, tE bin), the quantities Phase F of the refactor plan asks for.
This is the numeric backbone the figures illustrate; F2 and F3 are pictures of rows in here.

The three currencies
--------------------
The plan is emphatic that one metric must not be forced across both regimes:

  short tE  -> a YIELD statistic. Roman literally cannot see an event that peaks and ends
               inside a season gap, so sigma_Roman does not exist and any ratio is undefined.
               What matters is how many events Rubin recovers at all.
  long tE   -> a PRECISION statistic. Roman sees the event; the question is how much Rubin's
               year-round baseline sharpens the fit, especially the annual-parallax signal
               that turns tE into a lens mass.

So the table reports both, side by side, and never averages one into the other.

Characterization criterion: tE > 2*sigma_tE AND piE > 2*sigma_piE, deliberately Abrams et al.
2025's, so the Rubin-alone column is directly comparable to their published numbers.

Fields
------
Events are assigned to the nearest of the six Roman GBTDS field centres read from
RomanBaseline.dat, within one Roman FoV radius. Everything else is "outside" -- Rubin-only
sky, which is most of the survey region and must not be pooled with the Roman fields.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R

FOV_ROMAN = 0.3003          # deg, radius -- Bulge.h
DEFAULT_EDGES = "10,30,100,300,inf"


def roman_fields(path):
    lb = pd.read_csv(path, sep=r"\s+", comment="#", header=None,
                     names=["ID", "RA", "Dec", "l", "b", "time", "sig5"])
    return lb[["l", "b"]].drop_duplicates().reset_index(drop=True)


def assign_field(df, fields):
    """Nearest field centre within one Roman FoV radius, else 'outside'."""
    name = pd.Series("outside", index=df.index, dtype=object)
    best = pd.Series(np.inf, index=df.index)
    for i, row in fields.iterrows():
        d = np.hypot(df["lon"] - row["l"], df["lat"] - row["b"])
        hit = (d < FOV_ROMAN) & (d < best)
        best = best.where(~hit, d)
        name = name.where(~hit, f"F{i}({row['l']:+.3f},{row['b']:+.3f})")
    return name


def summarise(g):
    """One output row from one (field, tE bin) group."""
    detL, detR = R.detected(g, "rubin"), R.detected(g, "roman")
    out = {
        "N_events": len(g),
        # ---- yield: who saw it ----
        "N_rubin_only": int((detL & ~detR).sum()),
        "N_roman_only": int((detR & ~detL).sum()),
        "N_both": int((detL & detR).sum()),
        "N_neither": int((~detL & ~detR).sum()),
    }

    # ---- yield: gap-peaking events recovered by Rubin ----
    # Only meaningful where Roman actually observes; elsewhere Roman missed the event
    # because it never pointed there, which is a footprint fact, not a cadence one.
    gap = g[(g["t0zone"] == 1) & (g["ndw_R"] > 0)]
    out["N_gap_peaking"] = len(gap)
    out["frac_gap_seen_by_rubin"] = float(R.detected(gap, "rubin").mean()) if len(gap) else np.nan

    # ---- precision: per-event ratio first, then the median. Never a ratio of means. ----
    for p in ("tE", "piE"):
        r = R.ratio_joint_over(g, p, "roman").dropna()
        out[f"N_ratio_{p}"] = len(r)
        out[f"med_ratio_{p}"] = float(r.median()) if len(r) else np.nan
        out[f"q25_ratio_{p}"] = float(r.quantile(0.25)) if len(r) else np.nan
        out[f"q75_ratio_{p}"] = float(r.quantile(0.75)) if len(r) else np.nan

    # ---- characterisation gain ----
    cj, cr = R.characterized(g, "joint"), R.characterized(g, "roman")
    out["N_char_joint"] = int(cj.sum())
    out["N_char_roman"] = int(cr.sum())
    out["dN_char"] = int(cj.sum() - cr.sum())

    # ---- how often Rubin alone cannot be inverted at all ----
    det = g[R.detected(g, "joint")]
    out["frac_rubin_singular"] = float((det["okA_L"] == 0).mean()) if len(det) else np.nan
    return pd.Series(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("events")
    ap.add_argument("-o", "--out", default="f1_results_table.csv")
    ap.add_argument("--baseline", default="Baseline/RomanBaseline.dat")
    ap.add_argument("--provenance", default="files/MONTLMC/files/run_provenance.txt")
    ap.add_argument("--te-edges", default=DEFAULT_EDGES,
                    help=f"comma-separated tE bin edges in days (default {DEFAULT_EDGES})")
    ap.add_argument("--fields-only", action="store_true",
                    help="drop the 'outside' row (Rubin-only sky)")
    a = ap.parse_args()

    df = R.load_events(a.events)
    edges = [float(x) for x in a.te_edges.split(",")]
    labels = [f"{edges[i]:g}-{edges[i+1]:g} d" for i in range(len(edges) - 1)]
    df["teBin"] = pd.cut(df["tE"], edges, labels=labels, right=False)
    df["field"] = assign_field(df, roman_fields(a.baseline))

    # A binning that dumps nearly everything in one bin describes the population badly and
    # makes every per-bin number a restatement of the whole sample. Say so rather than let
    # a reader assume the bins were chosen to fit the data.
    share = df["teBin"].value_counts(normalize=True, dropna=True)
    if len(share) and share.max() > 0.5:
        print(f"WARNING: {share.idxmax()} holds {100*share.max():.0f}% of events -- "
              f"the tE bins do not resolve this population. Override with --te-edges.",
              file=sys.stderr)

    if a.fields_only:
        df = df[df["field"] != "outside"]

    tab = (df.groupby(["field", "teBin"], observed=True)
             .apply(summarise, include_groups=False)
             .reset_index())
    tab = tab[tab["N_events"] > 0]
    tab.to_csv(a.out, index=False)

    print(R.describe(a.events, a.provenance))
    print(f"\n{len(tab)} (field, tE bin) rows -> {a.out}\n")
    show = ["field", "teBin", "N_events", "N_rubin_only", "N_roman_only", "N_both",
            "N_char_joint", "N_char_roman", "dN_char", "med_ratio_tE", "N_ratio_tE"]
    with pd.option_context("display.width", 200, "display.max_columns", 40,
                           "display.max_rows", 60):
        print(tab[show].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
