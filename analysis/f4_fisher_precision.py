#!/usr/bin/env python3
"""What the Fisher matrices forecast, survey by survey.

For each survey partition (joint, Roman alone, Rubin alone), how well each parameter is measured
over the whole detected sample: the direct read-out of the three Fisher matrices per event.

Top row, the parameters the light curve and centroid fit:
  (a) sigma_tE / tE      Einstein-crossing time; tE ~ sqrt(Ml), half of the mass measurement.
  (b) sigma_piE / piE    Microlensing parallax, from the Earth's orbital motion; piE ~ 1/sqrt(Ml).
                         The parameter Roman's season gaps are expected to degrade.
  (c) sigma_tetE / tetE  Angular Einstein radius, from the ASTROMETRIC matrix (source-centroid
                         wobble), so gated on okB, not okA.

Bottom row, the payoff and whether to believe it:
  (d) sigma_Ml / Ml      Lens mass, never fitted: Ml = tetE / (kappa * piE), kappa = 8.144 mas/Msun,
                         so its fractional error is the quadrature sum of those on tetE and piE.
  (e) joint vs single    Per-event sigma_Ml/Ml, joint against each survey alone. Every point must lie
                         on or below 1:1, since the joint information matrix is the sum of the
                         per-survey ones. Points above the line are round-off on ill-conditioned
                         matrices.
  (f) condition number   of the normalized photometric matrix. Above ~1e9 double precision has lost
                         most of its digits and the inverse, hence every sigma from it, is not
                         trustworthy.

The CDFs are normalized to the full sample, not the measured subset. An event whose matrix did not
invert, or whose parameter was never free, carries the -1.0 sentinel; dropping those rows and
renormalizing would flatter the survey that fails most often. The y axis is

    (events in the sample with sigma/theta < x) / (events in the sample)

so a curve that saturates below 1.0 shows that survey never constrained the remaining events. The
saturation level is printed in the legend.

Scope:
--scope footprint (default) keeps events Roman observed (ndw_R > 0), the only sample where the
Roman-alone curve is meaningful.
--scope all keeps every joint-detected event; ignore the Roman curve and compare joint against
Rubin-alone, the whole-survey statistic.
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

INK = "#1a1a1a"
MUTED = "#6b6b6b"
SURFACE = "#ffffff"
GRID = "#f0f0f0"

# One hue per survey, identical in every panel.
COLOR = {"joint": "#0d366b", "roman": "#c2410c", "rubin": "#0e7490"}
LABEL = {"joint": "joint fit", "roman": "Roman alone", "rubin": "Rubin alone"}

# Fractional-precision panels: parameter -> (denominator column, axis label, gate note)
PANELS = [
    ("tE",   "tE",   r"$t_{\rm E}$",
     "(a)  Einstein-crossing time", "photometric matrix"),
    ("piE",  "piE",  r"$\pi_{\rm E}$",
     "(b)  Microlensing parallax", "photometric matrix"),
    ("tetE", "tetE", r"$\theta_{\rm E}$",
     "(c)  Angular Einstein radius", "astrometric matrix"),
    ("Ml",   None,   r"$M_{\rm L}$",
     "(d)  Lens mass  (derived)", r"$M_{\rm L}=\theta_{\rm E}/(\kappa\,\pi_{\rm E})$"),
]

XLIM = (1e-3, 1e3)
TARGET = 0.1          # a 10% measurement: the conventional "this parameter is measured"


def frac_precision(df, param, denom, survey):
    """Fractional 1-sigma forecast, NaN where the parameter was not measured.

    relMl is already fractional; the others are absolute and divided by the true value drawn.
    R.sigma has already masked the -1.0 sentinels.
    """
    s = R.sigma(df, param, survey)
    return s if denom is None else s / df[denom]


def draw_cdf(ax, df, w, param, denom, symbol, title, note, surveys):
    """CDF of a fractional forecast. Solid curve weighted, dashed curve the raw sample.

    The weighted curve describes the sky (the draws are not distributed like events); the raw
    curve describes the simulated sample and stays visible as a control.
    """
    n_total = len(df)
    w = np.asarray(w, dtype=float)
    w_total = w.sum()
    ax.set_facecolor(SURFACE)
    ax.axvline(TARGET, color="#b0b0b0", lw=0.9, ls=":", zorder=1)
    ax.text(TARGET, 0.02, " 10%", color=MUTED, fontsize=7, ha="left", va="bottom", zorder=1)

    rows = []
    for s in surveys:
        vs = frac_precision(df, param, denom, s)
        ok = (vs.notna() & (vs > 0)).to_numpy()
        v, wv = vs.to_numpy()[ok], w[ok]
        o = np.argsort(v)
        v, wv = v[o], wv[o]
        # Normalized to the full sample, so unmeasured events cost the curve height.
        y = np.cumsum(wv) / w_total
        sat = wv.sum() / w_total
        below = float(wv[v < TARGET].sum()) / w_total
        below_raw = float(np.sum(v < TARGET)) / n_total
        ax.step(v, y, where="post", color=COLOR[s], lw=1.9, zorder=3,
                label=f"{LABEL[s]}   {below*100:4.1f}% < 10%   "
                      f"(raw {below_raw*100:.1f}%, max {sat*100:.0f}%)")
        ax.step(v, np.arange(1, len(v) + 1) / n_total, where="post", color=COLOR[s],
                lw=0.9, ls=(0, (3, 2)), alpha=0.55, zorder=2)
        rows.append(dict(param=param, survey=s, n_sample=n_total, n_measured=len(v),
                         n_eff=R.kish_neff(w), frac_measured=sat,
                         frac_below_10pct=below, frac_below_10pct_unweighted=below_raw,
                         p10=R.weighted_quantile(v, wv, 0.10),
                         median=R.weighted_median(v, wv),
                         p90=R.weighted_quantile(v, wv, 0.90),
                         median_unweighted=float(np.median(v)) if len(v) else np.nan))

    ax.set_xscale("log")
    ax.set_xlim(*XLIM)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel(r"fractional 1$\sigma$ forecast", color=INK, fontsize=9)
    ax.set_ylabel("cumulative fraction of sample", color=INK, fontsize=9)
    ax.set_title(f"{title}\n{note}", color=INK, fontsize=10, loc="left", pad=7,
                 linespacing=1.5)
    ax.legend(loc="upper left", fontsize=7.2, frameon=False, labelcolor=INK,
              handlelength=1.5, borderpad=0.2)
    style(ax)
    # Parameter symbol as a faint panel watermark.
    ax.text(0.985, 0.03, symbol, transform=ax.transAxes, color="#d8d8d8", fontsize=21,
            ha="right", va="bottom", zorder=0)
    return rows


def draw_invariant(ax, df, w, surveys):
    """sigma_Ml/Ml, joint against each single survey. The 1:1 line is a physics bound."""
    ax.set_facecolor(SURFACE)
    w = np.asarray(w, dtype=float)
    j = frac_precision(df, "Ml", None, "joint")
    lo, hi = 1e-3, 1e3
    ax.plot([lo, hi], [lo, hi], color="#b0b0b0", lw=1.0, ls="--", zorder=2)
    ax.text(hi, hi, "  1:1", color=MUTED, fontsize=7.5, ha="left", va="center")

    rows = []
    nviol = 0
    for s in surveys:
        if s == "joint":
            continue
        x = frac_precision(df, "Ml", None, s)
        m = j.notna() & x.notna()
        r = (j[m] / x[m])
        wm = w[m.to_numpy()]
        med = R.weighted_median(r.to_numpy(), wm)
        ax.scatter(x[m], j[m], s=3.5, alpha=0.28, linewidths=0, color=COLOR[s], zorder=3,
                   label=f"vs {LABEL[s]}   median "
                         + r"$\sigma_{\rm joint}/\sigma_{\rm single}$ = "
                         + f"{med:.2f}   (n={int(m.sum())})")
        nviol += int((r > 1.0 + 1e-9).sum())
        rows.append(dict(param="Ml", survey=s, n_sample=len(df), n_measured=int(m.sum()),
                         n_eff=R.kish_neff(wm), frac_measured=np.nan,
                         frac_below_10pct=np.nan, frac_below_10pct_unweighted=np.nan,
                         p10=np.nan, median=med, p90=np.nan,
                         median_unweighted=float(r.median())))
        print(f"    Ml joint/{s}: median ratio {med:.3f} weighted, {r.median():.3f} raw; "
              f"{int((r > 1.0 + 1e-9).sum())} above 1 (round-off, cond > 1e9)")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$\sigma_{M_{\rm L}}/M_{\rm L}$, single survey", color=INK, fontsize=9)
    ax.set_ylabel(r"$\sigma_{M_{\rm L}}/M_{\rm L}$, joint", color=INK, fontsize=9)
    ax.set_title("(e)  The Fisher invariant\n"
                 r"$\sigma_{\rm joint}\leq\sigma_{\rm single}$: everything on or below 1:1",
                 color=INK, fontsize=10, loc="left", pad=7, linespacing=1.5)
    ax.text(0.03, 0.72,
            f"{nviol} point{'' if nviol == 1 else 's'} above the line"
            + ("" if nviol == 0 else "  (round-off, cond > $10^{9}$)"),
            transform=ax.transAxes, color=MUTED, fontsize=7.5, va="top")
    leg = ax.legend(loc="upper left", fontsize=7.2, frameon=False, labelcolor=INK,
                    handlelength=1.0, borderpad=0.2)
    for h in leg.legend_handles:
        h.set_alpha(1.0); h.set_sizes([18])
    style(ax)
    return rows


def draw_condition(ax, df, w, surveys):
    """CDF of the normalized photometric matrix's condition number."""
    ax.set_facecolor(SURFACE)
    w = np.asarray(w, dtype=float)
    w_total = w.sum()
    ax.axvline(1e9, color="#b0b0b0", lw=0.9, ls=":", zorder=1)
    ax.text(1e9, 0.02, "  $10^{9}$", color=MUTED, fontsize=7, ha="left", va="bottom")

    rows = []
    for s in surveys:
        cs = df[f"condA_{R.SURVEYS[s]}"].astype(float)
        ok = (cs > 0).to_numpy()
        c, wc = cs.to_numpy()[ok], w[ok]
        o = np.argsort(c)
        c, wc = c[o], wc[o]
        y = np.cumsum(wc) / w_total
        good = float(wc[c < 1e9].sum()) / w_total
        good_raw = float(np.sum(c < 1e9)) / len(df)
        ax.step(c, y, where="post", color=COLOR[s], lw=1.9, zorder=3,
                label=f"{LABEL[s]}   {good*100:4.1f}% well-conditioned (raw {good_raw*100:.1f}%)")
        rows.append(dict(param="condA", survey=s, n_sample=len(df), n_measured=len(c),
                         n_eff=R.kish_neff(w), frac_measured=wc.sum() / w_total,
                         frac_below_10pct=good, frac_below_10pct_unweighted=good_raw,
                         p10=R.weighted_quantile(c, wc, 0.10),
                         median=R.weighted_median(c, wc),
                         p90=R.weighted_quantile(c, wc, 0.90),
                         median_unweighted=float(np.median(c))))

    ax.set_xscale("log")
    ax.set_xlim(1e2, 1e18)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("condition number of the normalized photometric matrix",
                  color=INK, fontsize=9)
    ax.set_ylabel("cumulative fraction of sample", color=INK, fontsize=9)
    ax.set_title("(f)  Are the inverses trustworthy?\n"
                 "beyond $10^{9}$, double precision has lost the answer",
                 color=INK, fontsize=10, loc="left", pad=7, linespacing=1.5)
    ax.legend(loc="upper left", fontsize=7.2, frameon=False, labelcolor=INK,
              handlelength=1.5, borderpad=0.2)
    style(ax)
    return rows


def style(ax):
    ax.tick_params(colors=MUTED, labelsize=8)
    for sp in ax.spines.values():
        sp.set_color("#d0d0d0")
    ax.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("events", help="per-event table (test5.dat), or a cached subset of it")
    ap.add_argument("-o", "--out", default="f4_fisher_precision.png")
    ap.add_argument("--provenance", default=None,
                    help="run_provenance.txt; found automatically if omitted")
    ap.add_argument("--chunksize", type=int, default=500_000,
                    help="rows per read chunk; lower it if memory is tight (default 5e5)")
    ap.add_argument("--scope", choices=("footprint", "all"), default="footprint",
                    help="footprint = events Roman observed (ndw_R>0); the only sample "
                         "where the Roman-alone curve means anything. Default footprint.")
    ap.add_argument("--map", default=None,
                    help="MapLMC5.dat -- needed for the event-rate weight")
    ap.add_argument("--log", action="append", default=[],
                    help="run log(s), for sightlines whose map rows a killed run lost")
    ap.add_argument("--unweighted", action="store_true",
                    help="deliberately report the raw sample, with no event-rate weight")
    args = ap.parse_args()

    # Filter inside the chunked reader: joint-detected events are ~1% of the table and loading
    # the rest first can exhaust memory. Nothing is lost, since FisherM runs only on events that
    # passed the detection test and detJ == 1 contains every event with a Fisher matrix.
    df = R.load_events(args.events, keep=lambda c: c["detJ"] == 1,
                       chunksize=args.chunksize)
    print(R.describe(args.events, args.provenance))

    bad = R.check_monotonicity(df)
    print(f"  monotonicity sigma_joint <= sigma_single: "
          f"{'OK' if not bad else 'VIOLATIONS ' + repr(bad)}")

    if args.scope == "footprint":
        df = df[df["ndw_R"] > 0]
        scope_note = ("events Roman actually observed (ndw_R > 0) -- the only sample in "
                      "which the Roman-alone curve is meaningful")
    else:
        scope_note = ("every joint-detected event; Roman never pointed at most of these, "
                      "so read the joint/Rubin pair and ignore the Roman curve")
    surveys = ("joint", "roman", "rubin")
    # Weight after the scope cut: every fraction is over the sample the scope selected.
    w, wlabel = R.attach_weight(df, args.map, args.log, args.unweighted)
    w = w.to_numpy()
    print(f"  scope={args.scope}: {len(df)} events; weighting: {wlabel}")

    fig, axes = plt.subplots(2, 3, figsize=(16.2, 9.4), facecolor=SURFACE)
    rows = []
    for ax, (param, denom, symbol, title, note) in zip(axes.flat, PANELS):
        rows += draw_cdf(ax, df, w, param, denom, symbol, title, note, surveys)
    rows += draw_invariant(axes.flat[4], df, w, surveys)
    rows += draw_condition(axes.flat[5], df, w, surveys)

    fig.suptitle("Fisher-matrix forecast precision, per survey partition",
                 color=INK, fontsize=13.5, x=0.035, ha="left", y=0.975)
    fig.text(0.035, 0.945,
             f"{len(df):,} joint-detected events -- {scope_note}\n"
             f"solid: {wlabel}    dashed: the raw sample, unweighted",
             color=MUTED, fontsize=9, va="top")

    fig.subplots_adjust(left=0.05, right=0.985, bottom=0.105, top=0.885,
                        wspace=0.26, hspace=0.42)
    fig.text(0.035, 0.052,
             "Curves are normalized to the FULL sample, not to the events each survey "
             "managed to measure: an unmeasurable event costs the curve height rather "
             "than leaving the denominator.\n"
             "A curve that saturates below 1.0 never constrained the rest of the sample "
             "at all. -1.0 sentinels are masked by analysis/romanlib.py and never enter a "
             "quantile.\n"
             + R.describe(args.events, args.provenance),
             color=MUTED, fontsize=7.5, va="top", linespacing=1.7)

    fig.savefig(args.out, dpi=200, facecolor=SURFACE)
    print(f"wrote {args.out}")

    csv = os.path.splitext(args.out)[0] + ".csv"
    pd.DataFrame(rows).to_csv(csv, index=False, float_format="%.5g")
    print(f"wrote {csv}")


if __name__ == "__main__":
    main()
