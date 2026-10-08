#!/usr/bin/env python3
"""Results table: per (Roman field, tE bin) yield and precision statistics.

Yield and precision are reported side by side and never averaged into one another:

  short tE  -> YIELD. Roman cannot see an event that peaks and ends inside a season gap, so
               sigma_Roman does not exist and any ratio is undefined; what matters is how many
               events Rubin recovers at all.
  long tE   -> PRECISION. Roman sees the event; the question is how much Rubin's year-round
               baseline sharpens the fit, especially the annual-parallax signal that turns tE
               into a lens mass.

Characterization criterion: tE > 2*sigma_tE AND piE > 2*sigma_piE (Abrams et al. 2025), so the
Rubin-alone column is directly comparable to their published numbers.

Fields: events are assigned to the GBTDS field whose detectors image their sightline, from the
visit list given by --baseline, using the simulator's own coverage test (gbtds_geometry): 'F<i>',
or 'F<i>/F<j>' where the spring and autumn rolls put a different field there. Everything else is
"outside" (Rubin-only sky, including chip gaps) and is not pooled with the Roman fields. A run made
with the old six-centre layout needs that run's visit list
(Baseline/legacy_layout40395/RomanBaseline.dat): it is recognised and the old rule applied
(nearest of six centres within 0.3003 deg).
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R
import gbtds_geometry as G

DEFAULT_EDGES = "10,30,100,300,inf"


def assign_field(df, visits):
    """Field label per event, from its sightline (gbtds_geometry.field_label)."""
    keys = df[["lon", "lat"]].drop_duplicates()
    lab = dict(zip(zip(keys["lon"], keys["lat"]),
                   G.field_label(keys["lon"].to_numpy(), keys["lat"].to_numpy(), visits)))
    return pd.Series([lab[k] for k in zip(df["lon"], df["lat"])], index=df.index, dtype=object)


def summarise(g):
    """One output row from one (field, tE bin) group.

    Counts are raw counts (they describe the simulated sample). Fractions and medians are
    event-rate weighted, because they are statements about the sky, and each carries `N_eff`
    (the Kish effective sample size).
    """
    detL, detR = R.detected(g, "rubin"), R.detected(g, "roman")
    w = g["W"].to_numpy()
    out = {
        "N_events": len(g),
        "N_eff": R.kish_neff(w),
        # ---- yield: who saw it ----
        "N_rubin_only": int((detL & ~detR).sum()),
        "N_roman_only": int((detR & ~detL).sum()),
        "N_both": int((detL & detR).sum()),
        "N_neither": int((~detL & ~detR).sum()),
    }

    # ---- yield: gap-peaking events recovered by Rubin ----
    # Only where Roman observes; elsewhere it missed the event because it never pointed there.
    gap = g[(g["t0zone"] == 1) & (g["ndw_R"] > 0)]
    out["N_gap_peaking"] = len(gap)
    out["frac_gap_seen_by_rubin"] = (R.weighted_fraction(R.detected(gap, "rubin"),
                                                         gap["W"].to_numpy())
                                     if len(gap) else np.nan)
    out["frac_gap_seen_by_rubin_unw"] = (float(R.detected(gap, "rubin").mean())
                                         if len(gap) else np.nan)

    # ---- precision: per-event ratio first, then the median. Never a ratio of means. ----
    # The ratio is per event and unweighted; its median over a set of events is weighted.
    for p in ("tE", "piE"):
        r = R.ratio_joint_over(g, p, "roman")
        ok = r.notna()
        rv, rw = r[ok].to_numpy(), w[ok.to_numpy()]
        out[f"N_ratio_{p}"] = int(ok.sum())
        out[f"Neff_ratio_{p}"] = R.kish_neff(rw)
        out[f"med_ratio_{p}"] = R.weighted_median(rv, rw)
        out[f"q25_ratio_{p}"] = R.weighted_quantile(rv, rw, 0.25)
        out[f"q75_ratio_{p}"] = R.weighted_quantile(rv, rw, 0.75)
        out[f"med_ratio_{p}_unw"] = float(np.median(rv)) if rv.size else np.nan

    # ---- characterisation gain ----
    cj, cr = R.characterized(g, "joint"), R.characterized(g, "roman")
    out["N_char_joint"] = int(cj.sum())
    out["N_char_roman"] = int(cr.sum())
    out["dN_char"] = int(cj.sum() - cr.sum())
    out["frac_char_joint"] = R.weighted_fraction(cj, w)
    out["frac_char_roman"] = R.weighted_fraction(cr, w)

    # ---- how often Rubin alone cannot be inverted at all ----
    det = g[R.detected(g, "joint")]
    out["frac_rubin_singular"] = (R.weighted_fraction(det["okA_L"] == 0, det["W"].to_numpy())
                                  if len(det) else np.nan)
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
    ap.add_argument("--map", default=None,
                    help="MapLMC5.dat -- needed for the event-rate weight")
    ap.add_argument("--log", action="append", default=[],
                    help="run log(s), for sightlines whose map rows a killed run lost")
    ap.add_argument("--chunksize", type=int, default=500_000,
                    help="rows per read chunk; holds peak memory to one chunk plus survivors")
    ap.add_argument("--unweighted", action="store_true",
                    help="deliberately report the raw sample, with no event-rate weight")
    a = ap.parse_args()

    # Stream the table and drop rows of barren sightlines (no nsim, weight 0, no detection) as
    # they are read; reading it whole can exhaust memory on a small machine.
    df = R.load_events(a.events, keep=R.keep_weightable(a.map, a.log), chunksize=a.chunksize)
    w, wlabel = R.attach_weight(df, a.map, a.log, a.unweighted)
    df["W"] = w
    print(f"weighting: {wlabel}")
    edges = [float(x) for x in a.te_edges.split(",")]
    labels = [f"{edges[i]:g}-{edges[i+1]:g} d" for i in range(len(edges) - 1)]
    df["teBin"] = pd.cut(df["tE"], edges, labels=labels, right=False)
    df["field"] = assign_field(df, G.read_roman_visits(a.baseline))

    # A binning that puts most events in one bin makes every per-bin number a restatement of
    # the whole sample.
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
