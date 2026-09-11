#!/usr/bin/env python3
"""The H3 satellite-parallax comparison, before and after the derivative-reference fix.

WHAT THIS SHOWS, AND WHY IT IS A FIGURE RATHER THAN A TABLE
-----------------------------------------------------------
DEVIATIONS.md 36: FisherM differenced every derivative against the model value cached when the
light curve was generated, which belongs to the observer the run actually used. Re-evaluating
with Roman moved to Earth made that reference inconsistent with the perturbed model, injecting a
constant dm_sat/Delta into every derivative. The artefact is proportional to dm_sat, so it grows
with the observer separation -- and it therefore leaves three fingerprints that this figure puts
side by side with the same quantities after the fix.

  (a) sigma_tE. Moving an observer re-weights information; it cannot destroy the event
      timescale, which is set by the shape of a light curve both observatories still sample at
      the same epochs. Before: a median ratio of 4.88 with 85.5% of events degraded. That is not
      an information statement, and it is the clearest single sign that the two matrices were
      not two descriptions of one event.

  (b) sigma_piE against the observer separation. The mechanism REQUIRES the gain to grow with
      separation: a wider baseline resolves the parallax better. Before, it shrank
      (corr = +0.227) -- the artefact running forwards, not the physics running backwards.

  (c) THE CONTROL, and the lesson in it. Events with no Roman epochs near the peak have
      dm_sat = 0 identically, so the artefact vanishes and the ratio is exactly 1. It was 1
      before the fix too. A control that passes tells you the plumbing runs; it does not tell
      you the measurement is sound, because here it passed precisely on the events the
      measurement was not about.

Needs the project venv (.roman/bin/python).
"""
import argparse
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BG = "#fcfcfb"
C_OLD, C_NEW = "#c1121f", "#1f4e79"

OLD_COLS = ("lon lat tE u0 piE du_sat okA_sat okA_nosat sigtE_sat sigtE_nosat "
            "sigpiE_sat sigpiE_nosat nepL_pk nepR_pk w_area").split()
NEW_COLS = ("lon lat tE u0 piE tetE du_sat okA_sat okA_nosat okB_sat okB_nosat "
            "sigtE_sat sigtE_nosat sigpiE_sat sigpiE_nosat sigpiER_sat sigpiER_nosat "
            "sigtetE_sat sigtetE_nosat sigpiEb_sat sigpiEb_nosat relMl_sat relMl_nosat "
            "condA_sat condA_nosat condB_sat condB_nosat nepL_pk nepR_pk w_area").split()


def load(path, cols):
    df = pd.read_csv(path, sep=r"\s+", comment="#", names=cols, engine="python")
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna()
    ok = ((df.okA_sat == 1) & (df.okA_nosat == 1) &
          (df.sigpiE_sat > 0) & (df.sigpiE_nosat > 0) &
          (df.sigtE_sat > 0) & (df.sigtE_nosat > 0))
    d = df[ok].copy()
    d["ratio"] = d.sigpiE_sat / d.sigpiE_nosat
    d["ratio_tE"] = d.sigtE_sat / d.sigtE_nosat
    return d


def style(ax, title, xlabel, ylabel):
    ax.set_facecolor(BG)
    ax.set_title(title, fontsize=10.5, loc="left", pad=8)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)
    ax.grid(alpha=0.25, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_tE(ax, o, n):
    bins = np.logspace(-1, 2, 46)
    for d, c, lab in ((o, C_OLD, "before"), (n, C_NEW, "after")):
        v = d.ratio_tE.to_numpy()
        if v.size == 0:
            continue
        worse = (v > 1.01).mean()
        ax.hist(np.clip(v, bins[0], bins[-1]), bins=bins, density=True, histtype="stepfilled",
                color=c, alpha=0.45, lw=1.4, edgecolor=c,
                label=f"{lab}: median {np.median(v):.3f}, {worse:.1%} worse  (n={v.size:,})")
    ax.axvline(1.0, color="0.3", lw=1.1, ls=":")
    ax.set_xscale("log")
    style(ax, "(a) moving an observer cannot destroy the timescale\n"
              r"$\sigma(t_E)$ ratio, satellite / no-satellite",
          r"$\sigma(t_E)_{\rm sat}\,/\,\sigma(t_E)_{\rm nosat}$", "density")
    ax.legend(fontsize=7.8, frameon=False, loc="upper right")


def panel_dusat(ax, o, n):
    for d, c, lab in ((o, C_OLD, "before"), (n, C_NEW, "after")):
        if len(d) < 5:
            continue
        ax.scatter(d.du_sat, d.ratio, s=6, color=c, alpha=0.18, linewidths=0)
        nb = min(6, max(2, len(d) // 25))
        q = pd.qcut(d.du_sat, nb, labels=False, duplicates="drop")
        xs, ys = [], []
        for b in range(int(q.max()) + 1):
            s = q == b
            xs.append(float(np.median(d.du_sat[s])))
            ys.append(float(np.median(d.ratio[s])))
        r = float(np.corrcoef(np.log10(d.du_sat), np.log10(d.ratio))[0, 1])
        ax.plot(xs, ys, "o-", color=c, lw=2.0, ms=5,
                label=f"{lab}: corr = {r:+.3f}  (n={len(d):,})")
    ax.axhline(1.0, color="0.3", lw=1.1, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    style(ax, "(b) a wider baseline must buy MORE, not less\n"
              "the gain has to fall as separation grows, so corr must be negative",
          r"observer separation $\Delta u_{\rm sat}$ [$\theta_E$]",
          r"$\sigma(\pi_E)_{\rm sat}\,/\,\sigma(\pi_E)_{\rm nosat}$")
    ax.legend(fontsize=7.8, frameon=False, loc="lower right")


def panel_control(ax, ob, nb):
    bins = np.linspace(0.97, 1.03, 61)
    for d, c, lab in ((ob, C_OLD, "before"), (nb, C_NEW, "after")):
        v = d.ratio.to_numpy()
        if v.size == 0:
            continue
        ax.hist(np.clip(v, bins[0], bins[-1]), bins=bins, histtype="stepfilled",
                color=c, alpha=0.45, lw=1.4, edgecolor=c, density=True,
                label=f"{lab}: median {np.median(v):.6f}  (n={v.size:,})")
    ax.axvline(1.0, color="0.3", lw=1.1, ls=":")
    style(ax, "(c) the control passed BEFORE the fix too\n"
              "no Roman epochs near peak, so the artefact is identically zero here",
          r"$\sigma(\pi_E)$ ratio", "density")
    ax.legend(fontsize=7.8, frameon=False, loc="upper left")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("old", help="h3_pair.dat from the pre-fix run (15 columns)")
    ap.add_argument("new", help="h3_pair.dat from the fixed run (30 columns)")
    ap.add_argument("-o", "--out", default="figures/h3_before_after.png")
    a = ap.parse_args()

    o = load(a.old, OLD_COLS)
    n = load(a.new, NEW_COLS)
    oc, nc = o[o.nepR_pk > 0], n[n.nepR_pk > 0]
    ob, nb = o[o.nepR_pk == 0], n[n.nepR_pk == 0]
    print(f"before: {len(oc):,} Roman-covered, {len(ob):,} control")
    print(f"after : {len(nc):,} Roman-covered, {len(nb):,} control")
    if len(nc) < 5:
        sys.exit("the fixed run has too few Roman-covered events yet")

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.4))
    fig.patch.set_facecolor(BG)
    panel_tE(axes[0], oc, nc)
    panel_dusat(axes[1], oc, nc)
    panel_control(axes[2], ob, nb)
    fig.suptitle("Step H3: the satellite-parallax comparison, before and after the "
                 "derivative-reference fix (DEVIATIONS 36)", fontsize=12.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(a.out, dpi=150, facecolor=BG)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
