#!/usr/bin/env python3
"""Step W1 -- what the pooled event weight does to the published numbers.

Deviation 41 derives the weight; this script is its measurement, kept runnable so the table can
be regenerated rather than trusted. It reports each headline pooled statistic twice, unweighted
and weighted, with the Kish effective sample size beside the weighted one.

    .roman/bin/python analysis/w1_pooled_weight_check.py test5.dat \
        --map files/MONTLMC/files/MapLMC5.dat --log run.log --log run2.log

UNWEIGHTED IS THE CONTROL. The unweighted column must reproduce the published numbers exactly
(F4 20.4/14.3/6.3; F2 0.250/0.924/0.975/0.990) -- that is what shows the sample selection here
matches the figure scripts, so that any difference in the weighted column is the weight and not a
different event set.

--log is needed only for a damaged map file: the simulator never flushes it, so a killed run
loses its buffered tail (OPEN_ITEMS.md). The log prints `nsim` for every sightline.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import galaxy_model as G
import romanlib as R

TE_EDGES = [10.0, 30.0, 100.0, 300.0, np.inf]
TE_LABELS = ["10-30 d", "30-100 d", "100-300 d", ">300 d"]
TARGET = 0.1


def wfrac(mask, w):
    """Weighted fraction of a boolean mask."""
    m = np.asarray(mask.fillna(False) if hasattr(mask, "fillna") else mask, dtype=bool)
    w = np.asarray(w, dtype=float)
    return 100.0 * w[m].sum() / w.sum() if w.sum() else np.nan


def wmedian(v, w):
    v, w = np.asarray(v, float), np.asarray(w, float)
    o = np.argsort(v)
    c = np.cumsum(w[o]) / w.sum()
    return float(v[o][np.searchsorted(c, 0.5)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("events")
    ap.add_argument("--map", required=True, help="MapLMC5.dat, for per-sightline nsim")
    ap.add_argument("--log", action="append", default=[],
                    help="run log(s), to supply nsim for sightlines a damaged map file lost")
    ap.add_argument("--chunksize", type=int, default=500_000)
    a = ap.parse_args()

    df = R.load_events(a.events, keep=lambda c: c["detJ"] == 1, chunksize=a.chunksize)
    sl = R.load_sightlines(a.map)
    over = R.nsim_from_logs(a.log)

    worst, n_checked = G.check_against_map(sl)
    print(f"density port vs map file: worst {worst:.3f} dex over {n_checked} sightlines "
          f"(map writes 1 decimal, so 0.05 is the resolution)")

    w = R.event_weight(df, sl, nsim_override=over)
    df = df.assign(W=w)

    # ---- F4: fraction measured better than 10%, footprint sample ----
    fp = df[df["ndw_R"] > 0]
    print(f"\nF4 sample (joint-detected, Roman epochs): N = {len(fp)}, "
          f"Kish N_eff = {R.kish_neff(fp['W']):.0f}")
    print(f"{'':6s} {'unweighted joint/Roman/Rubin':>30s}   {'weighted':>22s}")
    for p in ("tE", "piE", "tetE"):
        cells = []
        for weighted in (False, True):
            ww = fp["W"].to_numpy() if weighted else np.ones(len(fp))
            fr = [wfrac(R.sigma(fp, p, s) / fp[p] < TARGET, ww)
                  for s in ("joint", "roman", "rubin")]
            cells.append("/".join(f"{f:.1f}" for f in fr))
        print(f"{p:6s} {cells[0]:>30s}   {cells[1]:>22s}")

    # ---- the tE distribution the sample is drawn from ----
    print("\nAll joint detections:")
    for name, ww in (("unweighted", np.ones(len(df))), ("weighted", df["W"].to_numpy())):
        print(f"  {name:10s} median tE {wmedian(df['tE'], ww):6.1f} d   "
              f"tE>200 d {wfrac(df['tE'] > 200, ww):5.2f}%   "
              f"N_eff {R.kish_neff(ww):7.0f}")

    # ---- F2: per-event ratio medians, in-gap and in-season ----
    d = df[df["t0zone"].isin([0, 1]) & (df["ndw_R"] > 0) & (df["dt_edge"] <= 60.0)].copy()
    d["teBin"] = pd.cut(d["tE"], TE_EDGES, right=False, labels=TE_LABELS)
    print("\nF2 median sigma_joint/sigma_Roman, unweighted | weighted [N, N_eff]:")
    for zone, zname in ((1, "gap"), (0, "season")):
        for p in ("tE", "piE"):
            cells = []
            for b in TE_LABELS:
                g = d[(d["t0zone"] == zone) & (d["teBin"] == b)]
                r = R.ratio_joint_over(g, p, "roman")
                ok = r.notna()
                if not ok.any():
                    cells.append("-")
                    continue
                v, ww = r[ok].to_numpy(), g.loc[ok[ok].index, "W"].to_numpy()
                cells.append(f"{np.median(v):.3f}|{wmedian(v, ww):.3f} "
                             f"[{len(v)},{R.kish_neff(ww):.0f}]")
            print(f"  {p:4s} peak in {zname:6s}: " + "   ".join(f"{c:>24s}" for c in cells))


if __name__ == "__main__":
    main()
