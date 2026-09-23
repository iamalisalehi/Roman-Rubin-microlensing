#!/usr/bin/env python3
"""Step Y3: the absolute yield as a function of the abundance F, and the per-unit-F coefficients.

WHY THIS EXISTS. y1_absolute_yield.py tabulates the yield at a handful of literature values of F.
But N is exactly linear in F (Deviation 54: F enters only through the prefactor F / <M>, and <M>
is fixed by the mass function, not by F), so the whole content of those tables is one
coefficient per selection -- the yield at F = 1, N_1 -- and the rule N(F) = F * N_1. This script
reads y1's CSV, recovers N_1 for every selection, checks that the linearity really holds across
the grid it was written at, writes the coefficients as a table, and draws N(F) with the
literature abundances marked, so a reader can take any F they prefer.

It also relates each population to the ordinary-lens (`bulge`) run, where F = 1 by definition:
eta = N_1(population) / N(bulge) is the fraction of ALL events a population would supply per
unit mass fraction. That ratio is the <sqrt(M)>/<M> suppression of the report, after detection
efficiency, in one number.

No new data: everything here is derived from y1_yields.csv.

    .roman/bin/python analysis/y3_yield_vs_F.py figures/yield_prefix_20260922/y1_yields.csv \\
        -o figures/yield_prefix_20260922/y3
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotstyle as ps        # noqa: E402

# Literature abundances (fraction of the Galaxy's stellar mass), as checked in Deviation 54.
# (label, F, populations it applies to)
LITERATURE_F = [
    # Labels are short because they sit inside the panel; the caption names the papers.
    ("SM26", 0.0045, ("bh",)),
    ("SS23", 0.019, ("bh",)),
    ("O20, L20, G00", 0.03, ("bh",)),
    ("S22", 0.01, ("ns",)),
    ("G00", 0.06, ("ns",)),
]

# Selections drawn in the figure: (CSV selection, survey colour key, linestyle, label).
CURVES = [
    ("Roman detects", "roman", "-", "Roman detects"),
    ("Rubin detects", "rubin", "-", "Rubin detects"),
    ("both detect", "joint", "-", "both detect"),
    ("Roman, sigma(M)/M < 10%", "roman", "--", r"Roman, $\sigma_M/M<10\%$"),
    ("joint, sigma(M)/M < 10%", "joint", "--", r"joint, $\sigma_M/M<10\%$"),
]

SCOPE = "Roman footprint only"


def coefficients(df):
    """N at F = 1 per (run, scope, selection), with a check that N / F is constant in F."""
    d = df.assign(k_obj=df.N_per_object / df.F, e_obj=df.err_per_object / df.F,
                  k_all=df.N_all_stars / df.F, e_all=df.err_all_stars / df.F)
    g = d.groupby(["run", "scope", "selection"], sort=False)
    spread = g.k_obj.agg(lambda s: s.max() / s.min() - 1.0)
    if (spread > 1e-9).any():
        raise SystemExit(f"N/F is not constant in F -- the CSV is not linear:\n{spread[spread > 1e-9]}")
    return g.first()[["n_mc", "k_obj", "e_obj", "k_all", "e_all"]].reset_index()


def figure(k, out):
    ps.use_paper_style()
    pops = [p for p in ("bh", "ns") if p in set(k.run)]
    fig, axes = ps.figure(width="double", height=3.0, ncols=len(pops), sharey=True)
    axes = np.atleast_1d(axes)
    F = np.geomspace(1e-3, 0.1, 50)
    for ax, pop, tag in zip(axes, pops, "ab"):
        sub = k[(k.run == pop) & (k.scope == SCOPE)].set_index("selection")
        for sel, col, ls, lab in CURVES:
            r = sub.loc[sel]
            ax.plot(F, F * r.k_obj, color=ps.SURVEY[col], ls=ls, label=lab)
            ax.fill_between(F, F * (r.k_obj - r.e_obj), F * (r.k_obj + r.e_obj),
                            color=ps.SURVEY[col], alpha=0.15, lw=0)
        for lab, f, applies in LITERATURE_F:
            if pop in applies:
                ax.axvline(f, color=ps.MUTED, lw=0.6, ls=":")
                ax.text(f * 1.06, 0.03, lab, rotation=90, transform=ax.get_xaxis_transform(),
                        ha="left", va="bottom", fontsize=5.5, color=ps.MUTED)
        ax.axhline(1.0, color=ps.MUTED, lw=0.6)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(F[0], F[-1])
        ax.set_ylim(0.03, 1e4)
        ax.set_xlabel(r"abundance $F$ (fraction of stellar mass)")
        ps.panel_label(ax, f"({tag}) " + ("black holes" if pop == "bh" else "neutron stars"),
                       loc="upper left")
    axes[0].set_ylabel("events in 10 yr, Roman footprint")
    # Above the panels, not inside: every corner of both panels is taken by a line or a label.
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, axes[0].get_position().y1), ncol=len(CURVES),
               frameon=False, fontsize=6.5)
    ps.stamp(fig, "y3_yield_vs_F.py from y1_yields.csv; per-object convention; band = MC 1-sigma; "
                  "pre-extinction-fix runs (Deviation 53)")
    return ps.save_figure(fig, os.path.join(out, "y3_yield_vs_F"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("csv", help="y1_yields.csv written by y1_absolute_yield.py")
    ap.add_argument("-o", "--out", required=True, help="output directory")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    k = coefficients(pd.read_csv(a.csv))
    lines = ["# Yield per unit abundance, N_1 = N(F) / F (Step Y3)", "",
             "Derived from `" + a.csv + "`. N(F) = F x N_1 exactly; N_1 is the yield at F = 1, "
             "over the 10-yr window. eta = N_1 / N(bulge) for the same scope and selection: the "
             "fraction of all events a population supplies per unit mass fraction.", ""]
    bulge = k[k.run == "bulge"].set_index(["scope", "selection"])
    for (run, scope), sub in k.groupby(["run", "scope"], sort=False):
        lines += [f"## {run}, {scope}", "",
                  "| selection | MC draws | N_1 per object | N_1 all stars | eta |",
                  "|---|---|---|---|---|"]
        for _, r in sub.iterrows():
            key = (scope, r.selection)
            eta = (r.k_obj / bulge.loc[key].k_obj) if (run != "bulge" and key in bulge.index) else np.nan
            lines.append(f"| {r.selection} | {r.n_mc:,} | {r.k_obj:.4g} ± {r.e_obj:.2g} "
                         f"| {r.k_all:.4g} ± {r.e_all:.2g} | {eta:.3g} |")
        lines.append("")
    with open(os.path.join(a.out, "y3_coefficients.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    k.to_csv(os.path.join(a.out, "y3_coefficients.csv"), index=False)
    print("\n".join(lines))
    print("wrote", figure(k, a.out))


if __name__ == "__main__":
    main()
