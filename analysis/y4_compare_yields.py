"""Step Y4: compare two y1_yields.csv files -- here the post-extinction-fix yields against the
pre-fix (September 22) ones.

Yields are linear in the abundance F (N(F) = F N_1), so every row is reduced to N_1 = N / F and
the F grids of the two files need not match. Per (population, scope, selection) it reports N_1
before and after, the ratio, and how many sigma apart they are (the two runs are independent
Monte Carlo samples, so the errors add in quadrature). The Monte Carlo counts n_mc are shown
too: they are sample sizes set by the run's stopping rule, not yields, and are NOT compared.

    .roman/bin/python analysis/y4_compare_yields.py \\
        figures/yield_prefix_20260922/y1_yields.csv figures/yield_20260925/y1_yields.csv \\
        -o figures/yield_20260925/compare_20260922.md
"""
import argparse

import numpy as np
import pandas as pd

KEYS = ["run", "scope", "selection"]


def per_unit_F(path):
    d = pd.read_csv(path)
    for c in ("N_per_object", "err_per_object", "N_all_stars", "err_all_stars"):
        d[c] = d[c] / d["F"]
    # One row per key: N_1 is the same at every F (checked, not assumed).
    g = d.groupby(KEYS, sort=False)
    spread = g["N_per_object"].agg(lambda x: np.ptp(x) / max(abs(x.mean()), 1e-300)).max()
    if spread > 1e-9:
        raise ValueError(f"{path}: N/F differs across F by {spread:.1e}; yields not linear in F")
    return g.first().reset_index()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    old, new = per_unit_F(a.old), per_unit_F(a.new)
    old["order"] = np.arange(len(old))
    m = old.merge(new, on=KEYS, suffixes=("_old", "_new"), how="outer", indicator=True)
    m = m.sort_values("order", kind="stable", na_position="last")   # the old file's row order

    lines = ["# Yields per unit F: post-extinction-fix vs 2026-09-22", "",
             f"old: `{a.old}` (commit {', '.join(old.commit.unique())})  ",
             f"new: `{a.new}` (commit {', '.join(new.commit.unique())})", "",
             "N_1 = N / F, events per resolved object ('per object'); ratio = new / old; "
             "z = (new - old) / sqrt(err_old^2 + err_new^2), Monte Carlo error only. "
             "n_mc = Monte Carlo sample sizes (set by the stopping rule), not yields.", ""]
    for c in ("N_all_stars",):
        m[f"ratio_{c}"] = m[f"{c}_new"] / m[f"{c}_old"]
    for (run, scope), g in m.groupby(["run", "scope"], sort=False):
        lines += [f"## {run} -- {scope}", "",
                  "| selection | n_mc old | n_mc new | N_1 old | N_1 new | ratio | z | ratio (all stars) |",
                  "|---|---:|---:|---:|---:|---:|---:|---:|"]
        for _, r in g.iterrows():
            if r["_merge"] != "both":
                lines.append(f"| {r.selection} | only in {r['_merge']} | | | | | | |")
                continue
            ratio = r.N_per_object_new / r.N_per_object_old
            z = (r.N_per_object_new - r.N_per_object_old) / np.hypot(r.err_per_object_old,
                                                                     r.err_per_object_new)
            lines.append(f"| {r.selection} | {int(r.n_mc_old):,} | {int(r.n_mc_new):,} | "
                         f"{r.N_per_object_old:,.4g} ± {r.err_per_object_old:.2g} | "
                         f"{r.N_per_object_new:,.4g} ± {r.err_per_object_new:.2g} | "
                         f"{ratio:.3f} | {z:+.1f} | {r.ratio_N_all_stars:.3f} |")
        lines.append("")
    open(a.out, "w").write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
