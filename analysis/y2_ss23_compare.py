#!/usr/bin/env python3
"""Step Y (follow-up): per-parameter Roman characterisation efficiency, matched to
Sajadian & Sahu 2023 (arXiv:2301.03812), AJ 165, 96.

WHY THIS EXISTS. The yield tables of y1_absolute_yield.py showed our Roman black-hole
DETECTION count sitting near SS23's while our CHARACTERISED count sat far below it. The
report's Table 4 suggested the deficit is not astrometric -- sigma(tetE) is within 20% of
theirs while sigma(piE) and sigma(tE) are ~4x worse -- but Table 4 is computed over the
3-1000 Msun population actually simulated, whereas SS23 use 2-50 Msun. Long events are
exactly what Roman's ~70-day seasons constrain badly, so the comparison had to be redone on
a MATCHED mass range before the cause could be named.

This script recomputes, for the Roman partition only, the event-rate-weighted fraction of
detections with sigma(X)/X below 1, 5 and 10% for X in {tE, piE, tetE, Ml}, over:
  - the full 3-1000 Msun population as run, and
  - the 3-50 Msun subset re-weighted to log-flat over that range (SS23's range is 2-50; we
    cannot go below the 3 Msun the run sampled, and say so rather than extrapolating).

The re-weighting is the same device as y1's mass_subset(): draws are log-flat in M, so
restricting to M < 50 and renormalising multiplies each surviving draw by
ln(1000/3)/ln(50/3). It changes the mass FUNCTION, not the run.

SS23 Table 1, dN/dM ~ M^-1 (their row matching our log-flat function), sparse observations
during the 2.3-yr gap, is hard-coded below for the comparison. Their "m" column requires all
three of mass, distance and proper motion to pass simultaneously, which is STRICTER than any
single-parameter column; ours is per parameter, so our numbers are the optimistic side of the
comparison and the deficit is if anything understated.

    .roman/bin/python analysis/y2_ss23_compare.py --run runs/prod_bh_20260917 -o figures/y2
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R          # noqa: E402

COLS = ["tE", "piE", "tetE", "Ml", "Vt", "Ds", "u0", "struc", "lon", "lat", "w_area",
        "detR", "ndw_R", "opt_1e6",
        "sigtE_R", "sigpiE_R", "sigtetE_R", "relMl_R", "okA_R", "okB_R"]

# SS23 Table 1, dN/dM ~ M^-1, "Sparse observations during the time gap" [%].
# Columns as printed: sigtE/tE, sigpiE/piE, sigtetE/tetE, sigMl/Ml, sigDl/Dl, sigmus/mus,
# sigmul/mul, m, Ne_BHs.  We compare the four we forecast.
SS23 = {
    "<=1%":  {"tE": 24.48, "piE": 7.55, "tetE": 67.89, "Ml": 5.30},
    "<=5%":  {"tE": 54.23, "piE": 21.38, "tetE": 96.75, "Ml": 20.81},
    "<=10%": {"tE": 66.95, "piE": 30.00, "tetE": 99.17, "Ml": 29.79},
}
SS23_M = {"<=1%": 3.92, "<=5%": 17.33, "<=10%": 25.61}   # all three params at once
SS23_N = {"<=1%": 3, "<=5%": 15, "<=10%": 22}            # events, their F1*F2 normalisation
SS23_NDET = 86.0   # = 22 / 0.2561, their detected total for this mass function

PARAMS = ["tE", "piE", "tetE", "Ml"]
THRESH = [("<=1%", 0.01), ("<=5%", 0.05), ("<=10%", 0.10)]


def rel_error(df, param):
    """sigma(param)/param for the Roman partition, NaN where not measured."""
    s = R.sigma(df, param, "roman")
    if param == "Ml":
        return s                      # relMl_* is already relative
    return s / df[param].astype(float)


def block(df, w, label, out):
    det = R.detected(df, "roman")
    wd = w.where(det, 0.0)
    tot = wd.sum()
    out.append(f"\n### {label}")
    out.append(f"\n- weighted Roman detections: {tot:.6g}   (raw draws {int(det.sum()):,})")
    out.append(f"- N_eff: {R.kish_neff(wd[det]):.6g}")
    out.append("\n| threshold | param | this work [%] | SS23 [%] | ratio |")
    out.append("|---|---|---|---|---|")
    for name, thr in THRESH:
        for p in PARAMS:
            r = rel_error(df, p)
            good = det & r.notna() & (r < thr)
            frac = 100.0 * w.where(good, 0.0).sum() / tot if tot > 0 else np.nan
            ss = SS23[name][p]
            out.append(f"| {name} | {p} | {frac:.2f} | {ss:.2f} | {frac/ss:.2f} |")
    return tot


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="run directory (the black-hole run)")
    ap.add_argument("-o", "--out", default="y2_ss23", help="output prefix")
    ap.add_argument("--chunksize", type=int, default=2_000_000)
    ap.add_argument("--mlo", type=float, default=3.0)
    ap.add_argument("--mhi", type=float, default=50.0)
    a = ap.parse_args()

    ev = R.find_events(a.run) if hasattr(R, "find_events") else None
    if ev is None:
        cands = [f for f in os.listdir(a.run) if f.endswith(".dat") and f.startswith("test")]
        if len(cands) != 1:
            sys.exit(f"cannot identify the event table in {a.run}: {cands}")
        ev = os.path.join(a.run, cands[0])
    mp = os.path.join(a.run, "files/MONTLMC/files")
    maps = [f for f in os.listdir(mp) if f.startswith("MapLMC")]
    if len(maps) != 1:
        sys.exit(f"cannot identify the map file in {mp}: {maps}")
    map_path = os.path.join(mp, maps[0])
    logs = [os.path.join(a.run, f) for f in os.listdir(a.run) if f.startswith("run") and f.endswith(".log")]

    print(f"events : {ev}\nmap    : {map_path}\nlogs   : {logs}", flush=True)

    df = R.load_events(ev, keep=R.keep_weightable(map_path, logs),
                       chunksize=a.chunksize, usecols=COLS)
    print(f"loaded {len(df):,} weightable draws", flush=True)

    cache = a.out + ".weights.npz"
    if os.path.exists(cache):
        z = np.load(cache)
        if len(z["w"]) != len(df):
            sys.exit(f"{cache} has {len(z['w'])} weights for {len(df)} rows; delete it and re-run")
        w = pd.Series(z["w"], index=df.index)
        wlabel = str(z["label"])
        print(f"weights from cache {cache}: {wlabel}", flush=True)
    else:
        w, wlabel = R.attach_weight(df, map_path=map_path, log_paths=logs)
        os.makedirs(os.path.dirname(cache) or ".", exist_ok=True)
        np.savez_compressed(cache, w=w.to_numpy(), label=wlabel)
        print(f"weights computed and cached -> {cache}: {wlabel}", flush=True)

    out = ["# Roman characterisation efficiency vs Sajadian & Sahu 2023",
           "",
           "Generated by `analysis/y2_ss23_compare.py`. SS23 = arXiv:2301.03812 Table 1, the",
           "`dN/dM ~ M^-1` block with sparse observations during the 2.3-yr gap -- the row whose",
           "mass function matches ours. Percentages are event-rate weighted fractions of Roman's",
           "own detections, so they are efficiencies and carry no assumption about abundance F.",
           "",
           "SS23's own `m` column (all of mass, distance and proper motion below the threshold",
           f"simultaneously) is {SS23_M['<=1%']}, {SS23_M['<=5%']}, {SS23_M['<=10%']}% at 1, 5, 10%,",
           f"giving {SS23_N['<=1%']}, {SS23_N['<=5%']}, {SS23_N['<=10%']} events out of ~{SS23_NDET:.0f} detections.",
           "Our columns are single-parameter and therefore the optimistic side of the comparison."]

    block(df, w, f"Full population as run, {df['Ml'].min():.3g}-{df['Ml'].max():.4g} Msun", out)

    sub = df[(df["Ml"] >= a.mlo) & (df["Ml"] <= a.mhi)].copy()
    boost = np.log(df["Ml"].max() / df["Ml"].min()) / np.log(a.mhi / a.mlo)
    ws = w.loc[sub.index] * boost
    out.append(f"\n(log-flat re-weighting factor for the subset: {boost:.4f})")
    block(sub, ws, f"Re-weighted to log-flat {a.mlo:g}-{a.mhi:g} Msun (SS23 use 2-50)", out)

    txt = "\n".join(out) + "\n"
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out + ".md", "w") as fh:
        fh.write(txt)
    print(txt)
    print(f"wrote {a.out}.md")


if __name__ == "__main__":
    main()
