#!/usr/bin/env python3
"""Four panels for the astrometric-deflection result.

Companion to h5_astrometric_shift.py, which plots the shift itself in detail.

  (a) Where the signal sits relative to the noise: the centroid shift is far below a single
      exposure and clears the noise only after averaging tens of thousands of exposures.
  (b) Who measures theta_E: Rubin's ground-based astrometry is not a capability here; the
      Einstein radius is Roman's.
  (c) What it buys: the lens mass needs theta_E from astrometry and piE from photometry, so it
      requires both telescopes.
  (d) That the matrix responds to the deflection: precision must improve as the deflection
      grows; a flat panel would mean the other three measure something else.

Input: a whitespace-separated column extract of the per-event table with at least lon, lat, Ml,
Vt, Ds, w_area, u0, tetE, ndwR, magb_F146 and the okB_*/sigtetE_*/relMl_* columns.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R
from cparams import P

BG = "#fcfcfb"
C_J, C_L, C_R = "#1f4e79", "#c1121f", "#2a9d8f"   # joint, Rubin, Roman

PROV_PATH = None   # set in main(): the run's provenance says which astrometric error model it used


def roman_ast_error(mag):
    """errRomanA() via romanlib, by the model the run used. Per EXPOSURE, in mas."""
    return R.roman_ast_error_for_run(mag, PROV_PATH)


def cdf(ax, v, w=None, **kw):
    """Cumulative distribution, by weight when given."""
    v = np.asarray(v, dtype=float)
    if v.size == 0:
        return
    w = np.ones_like(v) if w is None else np.asarray(w, dtype=float)
    o = np.argsort(v)
    ax.plot(v[o], np.cumsum(w[o]) / w.sum(), **kw)


def panel_signal(ax, cov, w):
    sig1 = roman_ast_error(cov.magb_F146.to_numpy())
    med1 = R.weighted_median(sig1, w)
    n = R.weighted_median(cov.ndwR, w)
    stacked = med1 / np.sqrt(n)

    cdf(ax, cov.shift_max, w, color=C_J, lw=1.8, label="max centroid shift")
    ax.axvline(med1, color=C_L, lw=1.4, ls="--")
    ax.axvline(stacked, color=C_R, lw=1.4, ls="--")
    ax.text(med1, 0.55, f"  one exposure\n  {med1:.2f} mas", color=C_L, fontsize=7.5, va="center")
    ax.text(stacked, 0.22, f"{stacked:.4f} mas  \nall {int(n):,} exposures  ",
            color=C_R, fontsize=7.5, va="center", ha="right")
    ax.set_xscale("log")
    frac = R.weighted_fraction(cov.shift_max.to_numpy() > sig1, w)
    style(ax, f"(a) the signal lives between the two noise levels\n"
              f"only {frac:.2%} of events exceed a single exposure",
          "centroid shift [mas]", "cumulative fraction")
    ax.legend(fontsize=8, frameon=False, loc="upper left")


def panel_tetE(ax, cov, w):
    w_total = float(np.sum(w))
    for lab, sc, oc, c in (("joint", "sigtetE_J", "okB_J", C_J),
                           ("Rubin only", "sigtetE_L", "okB_L", C_L),
                           ("Roman only", "sigtetE_R", "okB_R", C_R)):
        m = (cov[oc] == 1) & (cov[sc] > 0)
        r = (cov.loc[m, sc] / cov.loc[m, "tetE"]).to_numpy()
        wm = np.asarray(w)[m.to_numpy()]
        # Normalised to the whole sample, so events this survey cannot measure lower the curve.
        f10 = float(wm[r < 0.10].sum()) / w_total
        cdf(ax, r, wm * (wm.sum() / w_total) if wm.sum() else wm,
            color=c, lw=1.8, label=f"{lab}  ({f10:.1%} < 10%)")
    ax.axvline(0.10, color="0.35", lw=1.0, ls=":")
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 1e2)
    style(ax, "(b) the Einstein radius is Roman's measurement",
          r"$\sigma(\theta_E)/\theta_E$", "cumulative fraction")
    ax.legend(fontsize=8, frameon=False, loc="upper left")


def panel_mass(ax, cov, w):
    thr = np.array([1.0, 0.3, 0.1])
    x = np.arange(len(thr))
    bw = 0.26
    got = {}
    for i, (lab, c, col) in enumerate((("joint", C_J, "relMl_J"),
                                       ("Rubin only", C_L, "relMl_L"),
                                       ("Roman only", C_R, "relMl_R"))):
        f = [R.weighted_fraction((cov[col] > 0) & (cov[col] < t), w) for t in thr]
        got[col] = f
        ax.bar(x + (i - 1) * bw, f, bw, color=c, label=lab)
    # events the joint fit measures that NEITHER survey measures alone
    for k, t in enumerate(thr):
        j = (cov.relMl_J > 0) & (cov.relMl_J < t)
        l = (cov.relMl_L > 0) & (cov.relMl_L < t)
        r = (cov.relMl_R > 0) & (cov.relMl_R < t)
        only = int((j & ~l & ~r).sum())
        ax.text(x[k], got["relMl_J"][k] + 0.012, f"+{only} joint-only",
                ha="center", fontsize=7.5, color=C_J)
    ax.set_xticks(x)
    ax.set_xticklabels([f"< {t:.0%}" for t in thr])
    style(ax, "(c) the lens mass needs BOTH matrices\n"
              r"$M_L=\theta_E/(\kappa\pi_E)$: $\theta_E$ from Roman, $\pi_E$ from the baseline",
          "fractional precision on $M_L$", "fraction of Roman-observed detections")
    ax.legend(fontsize=8, frameon=False)


def panel_monotone(ax, cov, w):
    m = (cov.okB_J == 1) & (cov.sigtetE_J > 0)
    c = cov[m]
    wm = np.asarray(w)[m.to_numpy()]
    rel = (c.sigtetE_J / c.tetE).to_numpy()
    ax.hexbin(c.shift_max, rel, xscale="log", yscale="log",
              gridsize=40, cmap="Blues", mincnt=1, linewidths=0)
    q = pd.qcut(c.shift_max, 5, labels=False, duplicates="drop")
    xs, ys = [], []
    for b in range(int(q.max()) + 1):
        s = (q == b).to_numpy()
        xs.append(R.weighted_median(c.shift_max.to_numpy()[s], wm[s]))
        ys.append(R.weighted_median(rel[s], wm[s]))
    ax.plot(xs, ys, "o-", color=C_L, lw=1.8, ms=5, label="quintile medians")
    r = float(np.corrcoef(np.log10(c.shift_max), np.log10(rel))[0, 1])
    style(ax, f"(d) precision tracks the deflection, as it must\n"
              f"corr(log $\\delta_{{max}}$, log $\\sigma/\\theta_E$) = {r:+.3f}",
          r"max centroid shift $\delta_{max}$ [mas]", r"$\sigma(\theta_E)/\theta_E$")
    ax.legend(fontsize=8, frameon=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("extract", help="column extract of test5.dat, with magb_F146")
    ap.add_argument("-o", "--out", default="figures/h5_astrometry_summary.png")
    ap.add_argument("--map", default=None,
                    help="MapLMC5.dat -- needed for the event-rate weight")
    ap.add_argument("--log", action="append", default=[],
                    help="run log(s), for sightlines whose map rows a killed run lost")
    ap.add_argument("--provenance", default=None,
                    help="the run's run_provenance.txt (default: searched for beside the extract)")
    ap.add_argument("--unweighted", action="store_true",
                    help="deliberately report the raw sample, with no event-rate weight")
    a = ap.parse_args()

    global PROV_PATH
    PROV_PATH = R.find_provenance(a.provenance, near=a.extract)
    print(f"Roman astrometric error model of this run: {R.roman_ast_model(PROV_PATH)}")
    d = pd.read_csv(a.extract, sep=r"\s+")
    if "magb_F146" not in d.columns:
        sys.exit("extract lacks magb_F146; re-extract it with column 70 of the event table")
    need = ("lon", "lat", "Ml", "Vt", "Ds", "w_area")
    missing = [c for c in need if c not in d.columns]
    if missing and not a.unweighted:
        sys.exit(f"extract lacks {missing}, which the event-rate weight needs. "
                 "Re-extract with these columns, or pass "
                 "--unweighted to report raw-sample numbers deliberately.")
    r2 = np.sqrt(2.0)
    u0, tetE = d.u0.to_numpy(), d.tetE.to_numpy()
    d["shift_max"] = np.where(u0 <= r2, tetE / np.sqrt(8.0), tetE * u0 / (u0 * u0 + 2.0))
    cov = d[d.ndwR > 0].copy()

    w, wlabel = R.attach_weight(cov, a.map, a.log, a.unweighted)
    w = w.to_numpy()
    print(f"weighting: {wlabel}")

    fig, axes = plt.subplots(2, 2, figsize=(12.6, 9.4))
    fig.patch.set_facecolor(BG)
    panel_signal(axes[0, 0], cov, w)
    panel_tetE(axes[0, 1], cov, w)
    panel_mass(axes[1, 0], cov, w)
    panel_monotone(axes[1, 1], cov, w)
    fig.suptitle(f"Astrometric microlensing in the Roman+Rubin joint fit — "
                 f"{len(cov):,} Roman-observed detections of {len(d):,}",
                 fontsize=12.5, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(a.out, dpi=150, facecolor=BG)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
