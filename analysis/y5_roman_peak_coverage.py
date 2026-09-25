#!/usr/bin/env python3
"""Step Y5: how many of Roman's "detections" saw the peak at all.

WHY. OPEN_ITEMS.md ("Roman 'detects' events whose peak falls years outside its mission"): 50,000
Roman exposures can pass delta-chi2 >= 500 on the slow wing of an event whose peak Roman never
observed. The post-extinction-fix Roman yields doubled, so the share of such wing detections
decides how the Roman rows of the yield tables should be read.

For each run it splits the Roman-detected yield (per object, F = 1, the y1 weight exactly --
this reuses y1's Run) by
  t0zone   in-season / in-gap / off-mission (peak before Roman's first or after its last epoch)
  nepR_pk  Roman epochs within +-2 tE of t0; 0 = Roman never saw the magnified part
and writes the yield-weighted fractions with their Monte Carlo counts.

    .roman/bin/python analysis/y5_roman_peak_coverage.py \\
        --run bh=runs/prod_bh_20260924 --run ns=runs/prod_ns_20260924 \\
        --run bulge=runs/prod_bulge_20260924 -o figures/yield_20260925/y5_peak_coverage.md
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R                  # noqa: E402
import y1_absolute_yield as Y         # noqa: E402

Y.COLS = Y.COLS + [c for c in ("t0zone", "nepR_pk") if c not in Y.COLS]


def rows(run):
    d = run.df
    det = d["detR"].to_numpy() == 1
    y = d["y"].to_numpy()
    tot = y[det].sum()
    out = []
    for label, sel in [
            ("peak in a Roman season", d["t0zone"].to_numpy() == 0),
            ("peak in a gap between seasons", d["t0zone"].to_numpy() == 1),
            ("peak outside Roman's mission", d["t0zone"].to_numpy() == 2),
            ("no Roman epoch within +-2 tE of the peak", d["nepR_pk"].to_numpy() == 0)]:
        s = det & sel
        out.append((label, int(s.sum()), y[s].sum() / tot))
    return int(det.sum()), tot, out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", action="append", required=True, help="name=directory")
    ap.add_argument("-o", "--out", required=True, help="output .md")
    ap.add_argument("--chunksize", type=int, default=500_000)
    a = ap.parse_args()
    lines = ["# Roman detections by peak coverage (Step Y5)", "",
             "Fractions of the Roman-detected yield (per object, F = 1, y1's weight), whole scan. "
             "n_mc = Monte Carlo count in that class. The first three rows partition the total.", ""]
    for spec in a.run:
        name, directory = spec.split("=", 1)
        run = Y.Run(name, directory, a.chunksize)
        n, tot, out = rows(run)
        lines += [f"## {name}  (Roman detections: n_mc {n:,}, N_1 {tot:,.4g})", "",
                  "| class | n_mc | fraction of yield |", "|---|---:|---:|"]
        lines += [f"| {lab} | {k:,} | {f:.3f} |" for lab, k, f in out]
        lines.append("")
        del run
        print("\n".join(lines[-8:]), flush=True)
    open(a.out, "w").write("\n".join(lines))


if __name__ == "__main__":
    main()
