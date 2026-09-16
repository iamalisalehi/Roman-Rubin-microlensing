#!/usr/bin/env python3
"""Step W1 -- validate the event-rate weight against a published survey measurement.

The weighted median over DETECTED events answers nothing on its own: it carries this pipeline's
detection test, and no survey publishes that. The INTRINSIC distribution does compare. Applying
the weight to every draw, detected or not, estimates the same thing an efficiency-corrected
survey measurement does, and Mroz et al. 2019 (ApJS 244, 29; arXiv:1906.02210) publishes it for
8 yr of OGLE-IV over 121 bulge fields: each event reweighted by 1/eff(tE), giving a mean Einstein
timescale of 22 d in the central bins (within ~3 deg of the Galactic centre), rising to 32 d at
l ~ +8 deg and l ~ -6 deg.

    .roman/bin/python analysis/w1_intrinsic_te.py test5.dat

THE COMPARISON IS OF A POPULATION MEAN, NOT OF A MEASUREMENT. OGLE reweights its own detections
by its own efficiency and works at u0 < 1; this pipeline draws u0 < 3 and has its own mass
function and kinematics. What the check can catch -- and did -- is a sampler whose tE
distribution is wrong by a factor of 2.5 (Deviation 41).

Reads the whole table, so it uses `usecols` rather than the row filter every other script uses:
a statistic over all draws cannot throw rows away.
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import galaxy_model as G
import romanlib as R

COLS = ["tE", "Vt", "Ml", "Ds", "lon", "lat", "w_area", "detJ"]


def weighted_stats(tE, w):
    w = np.asarray(w, float)
    o = np.argsort(tE)
    c = np.cumsum(w[o]) / w.sum()
    return dict(mean=float(np.sum(w * tE) / w.sum()),
                median=float(tE[o][np.searchsorted(c, 0.5)]),
                f100=float(100 * w[tE > 100].sum() / w.sum()),
                f200=float(100 * w[tE > 200].sum() / w.sum()),
                n=int(tE.size),
                neff=R.kish_neff(w))


def show(label, tE, w):
    s = weighted_stats(tE, w)
    print(f"  {label:34s} mean {s['mean']:6.1f} d  median {s['median']:6.1f} d  "
          f">100 d {s['f100']:5.1f}%  >200 d {s['f200']:5.1f}%  "
          f"N {s['n']:>9,}  Neff {s['neff']:9.0f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("events")
    ap.add_argument("--chunksize", type=int, default=1_000_000)
    a = ap.parse_args()

    df = R.load_events(a.events, usecols=COLS, chunksize=a.chunksize)
    print(f"{len(df):,} draws, {int((df['detJ'] == 1).sum()):,} joint-detected")

    # nsim is the number of draws a sightline took, and every draw writes one row -- so count
    # them here. That also covers BARREN sightlines, which produce draws but never reach the
    # aggregation block, so they appear in neither the map file nor the log.
    key = pd.Series(list(zip(df["lon"].round(3), df["lat"].round(3))))
    nsim = key.value_counts().to_dict()

    W = np.empty(len(df))
    Ds, Ml, Vt, wa = (df[c].to_numpy() for c in ("Ds", "Ml", "Vt", "w_area"))
    worst_interp = 0.0
    rng = np.random.default_rng(0)
    for k, pos in pd.Series(np.arange(len(df))).groupby(key).groups.items():
        pos = np.asarray(pos)
        prof = G.density_profile(*k)
        ds = Ds[pos]
        # Z(Ds) is smooth in Ds; tabulating it on 200 points and interpolating avoids building
        # an (events x 9500) matrix per sightline. The error is checked, not assumed.
        grid = np.linspace(ds.min(), ds.max(), 200)
        Z = np.interp(ds, grid, G.lens_distance_norm(prof, grid))
        if worst_interp == 0.0:
            sub = rng.choice(len(ds), size=min(300, len(ds)), replace=False)
            worst_interp = float(np.max(np.abs(Z[sub] / G.lens_distance_norm(prof, ds[sub]) - 1)))
        W[pos] = wa[pos] * prof.Nstart / nsim[k] * np.sqrt(Ml[pos]) * Vt[pos] * Z
    print(f"Z interpolation: worst relative error {worst_interp:.1e}")

    tE = df["tE"].to_numpy()
    lat, lon = df["lat"].to_numpy(), df["lon"].to_numpy()
    print("\nINTRINSIC (all draws) -- the efficiency-corrected distribution a survey measures:")
    show("all sightlines, weighted", tE, W)
    show("all sightlines, unweighted", tE, np.ones(len(df)))
    print("   OGLE-IV (Mroz+2019): mean 22 d in the central bins, 32 d at l ~ 8 deg")

    print("\nby |b|, weighted:")
    for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 5), (5, 90)):
        m = (np.abs(lat) >= lo) & (np.abs(lat) < hi)
        if m.any():
            show(f"|b| in [{lo},{hi})", tE[m], W[m])

    print("\nby l, weighted (|b| in [1,5), to hold latitude roughly fixed):")
    band = (np.abs(lat) >= 1) & (np.abs(lat) < 5)
    for lo, hi in ((-5, -2), (-2, 0), (0, 2), (2, 5)):
        m = band & (lon >= lo) & (lon < hi)
        if m.any():
            show(f"l in [{lo:+d},{hi:+d})", tE[m], W[m])

    print("\nDETECTED (joint) -- carries the detection test, NOT comparable to the survey:")
    d = df["detJ"].to_numpy() == 1
    show("joint-detected, weighted", tE[d], W[d])
    show("joint-detected, unweighted", tE[d], np.ones(int(d.sum())))


if __name__ == "__main__":
    main()
