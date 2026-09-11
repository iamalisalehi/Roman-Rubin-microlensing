#!/usr/bin/env python3
"""Four panels for the astrometric-deflection result (PROGRESS.md 5f).

Companion to h5_astrometric_shift.py, which plots the shift itself in detail. This one plots
the four claims the result actually rests on, each with the check that makes it believable:

  (a) WHERE THE SIGNAL SITS relative to the noise. The centroid shift is far below a single
      exposure and only clears the noise after averaging ~50,000 of them. The two vertical
      lines are the whole story of the measurement, and the gap between them is the assumption
      recorded in OPEN_ITEMS.md.

  (b) WHO MEASURES theta_E. Rubin's ground-based astrometry is not a capability here; the
      Einstein radius is Roman's.

  (c) WHAT IT BUYS -- the lens mass, which needs theta_E from astrometry AND piE from
      photometry, so it is the one quantity that genuinely requires both telescopes.

  (d) THAT THE MATRIX IS RESPONDING TO THE DEFLECTION AT ALL: precision must improve as the
      deflection grows. If this panel were flat, every number in the other three would be
      measuring something else.

Input: the column extract described in h5_crosscheck.py (needs magb_F146).
Needs the project venv (.roman/bin/python) -- /usr/bin/python3 has no pandas.
"""
import argparse
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BG = "#fcfcfb"
C_J, C_L, C_R = "#1f4e79", "#c1121f", "#2a9d8f"   # joint, Rubin, Roman

ROMAN_AST_FLOOR, ROMAN_AST_MFLR, ROMAN_AST_MBKG = 1.1, 20.62, 23.5
ROMAN_AST_SLOPE_SRC, ROMAN_AST_SBKG, ROMAN_AST_SLOPE_BKG = 0.3329, 10.0, 0.4


def roman_ast_error(mag):
    mag = np.asarray(mag, dtype=float)
    out = np.full(mag.shape, ROMAN_AST_FLOOR)
    mid = (mag > ROMAN_AST_MFLR) & (mag <= ROMAN_AST_MBKG)
    out[mid] = ROMAN_AST_FLOOR * 10.0 ** (ROMAN_AST_SLOPE_SRC * (mag[mid] - ROMAN_AST_MFLR))
    hi = mag > ROMAN_AST_MBKG
    out[hi] = ROMAN_AST_SBKG * 10.0 ** (ROMAN_AST_SLOPE_BKG * (mag[hi] - ROMAN_AST_MBKG))
    return out


def style(ax, title, xlabel, ylabel):
    ax.set_facecolor(BG)
    ax.set_title(title, fontsize=10.5, loc="left", pad=8)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)
    ax.grid(alpha=0.25, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def cdf(ax, v, **kw):
    v = np.sort(np.asarray(v, dtype=float))
    if v.size == 0:
        return
    ax.plot(v, np.arange(1, v.size + 1) / v.size, **kw)


def panel_signal(ax, cov):
    sig1 = roman_ast_error(cov.magb_F146.to_numpy())
    med1 = float(np.median(sig1))
    n = float(np.median(cov.ndwR))
    stacked = med1 / np.sqrt(n)

    cdf(ax, cov.shift_max, color=C_J, lw=1.8, label="max centroid shift")
    ax.axvline(med1, color=C_L, lw=1.4, ls="--")
    ax.axvline(stacked, color=C_R, lw=1.4, ls="--")
    ax.text(med1, 0.55, f"  one exposure\n  {med1:.2f} mas", color=C_L, fontsize=7.5, va="center")
    ax.text(stacked, 0.22, f"{stacked:.4f} mas  \nall {int(n):,} exposures  ",
            color=C_R, fontsize=7.5, va="center", ha="right")
    ax.set_xscale("log")
    frac = float((cov.shift_max.to_numpy() > sig1).mean())
    style(ax, f"(a) the signal lives between the two noise levels\n"
              f"only {frac:.2%} of events exceed a single exposure",
          "centroid shift [mas]", "cumulative fraction")
    ax.legend(fontsize=8, frameon=False, loc="upper left")


def panel_tetE(ax, cov):
    for lab, sc, oc, c in (("joint", "sigtetE_J", "okB_J", C_J),
                           ("Rubin only", "sigtetE_L", "okB_L", C_L),
                           ("Roman only", "sigtetE_R", "okB_R", C_R)):
        m = (cov[oc] == 1) & (cov[sc] > 0)
        r = (cov.loc[m, sc] / cov.loc[m, "tetE"]).to_numpy()
        f10 = (r < 0.10).mean()
        cdf(ax, r, color=c, lw=1.8, label=f"{lab}  ({f10:.1%} < 10%)")
    ax.axvline(0.10, color="0.35", lw=1.0, ls=":")
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 1e2)
    style(ax, "(b) the Einstein radius is Roman's measurement",
          r"$\sigma(\theta_E)/\theta_E$", "cumulative fraction")
    ax.legend(fontsize=8, frameon=False, loc="upper left")


def panel_mass(ax, cov):
    thr = np.array([1.0, 0.3, 0.1])
    x = np.arange(len(thr))
    w = 0.26
    got = {}
    for i, (lab, c, col) in enumerate((("joint", C_J, "relMl_J"),
                                       ("Rubin only", C_L, "relMl_L"),
                                       ("Roman only", C_R, "relMl_R"))):
        f = [float(((cov[col] > 0) & (cov[col] < t)).mean()) for t in thr]
        got[col] = f
        ax.bar(x + (i - 1) * w, f, w, color=c, label=lab)
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


def panel_monotone(ax, cov):
    m = (cov.okB_J == 1) & (cov.sigtetE_J > 0)
    c = cov[m]
    rel = (c.sigtetE_J / c.tetE).to_numpy()
    ax.hexbin(c.shift_max, rel, xscale="log", yscale="log",
              gridsize=40, cmap="Blues", mincnt=1, linewidths=0)
    q = pd.qcut(c.shift_max, 5, labels=False, duplicates="drop")
    xs, ys = [], []
    for b in range(int(q.max()) + 1):
        s = q == b
        xs.append(float(np.median(c.shift_max[s])))
        ys.append(float(np.median(rel[s])))
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
    a = ap.parse_args()

    d = pd.read_csv(a.extract, sep=r"\s+")
    if "magb_F146" not in d.columns:
        sys.exit("extract lacks magb_F146; re-extract with $70 (see h5_crosscheck.py)")
    r2 = np.sqrt(2.0)
    u0, tetE = d.u0.to_numpy(), d.tetE.to_numpy()
    d["shift_max"] = np.where(u0 <= r2, tetE / np.sqrt(8.0), tetE * u0 / (u0 * u0 + 2.0))
    cov = d[d.ndwR > 0].copy()

    fig, axes = plt.subplots(2, 2, figsize=(12.6, 9.4))
    fig.patch.set_facecolor(BG)
    panel_signal(axes[0, 0], cov)
    panel_tetE(axes[0, 1], cov)
    panel_mass(axes[1, 0], cov)
    panel_monotone(axes[1, 1], cov)
    fig.suptitle(f"Astrometric microlensing in the Roman+Rubin joint fit — "
                 f"{len(cov):,} Roman-observed detections of {len(d):,}",
                 fontsize=12.5, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(a.out, dpi=150, facecolor=BG)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
