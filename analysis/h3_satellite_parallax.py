#!/usr/bin/env python3
"""Step H3 -- what Roman's displacement to L2 buys for the parallax forecast.

THE EXPERIMENT
--------------
Every detected event is characterised TWICE inside one run: once with Roman at L2 and once
with `L2_OFFSET_AU` zeroed, the same draw, the same epochs, the same photometry, only the
observer moved (`--pair-satellite`). Each row of `h3_pair.dat` is therefore a genuine pair, and
the ratio below is a per-event ratio.

WHY NOT TWO RUNS, WHICH IS WHAT THE PLAN ASKS FOR
-------------------------------------------------
Because the physics forbids it (DEVIATIONS.md 35). The RNG is one mt19937_64 stream and the
per-event path draws conditionally on `acceptRubin or acceptRoman`; moving the observer changes
what is detectable, so the streams fork. Measured on the v3 pair: byte-identical for 1,740,091
events across 459 sightlines, then divergent from the first footprint sightline onward -- the
two runs stay matched only where the effect is identically zero.

The unpaired fallback was tried and its own control rejected it: the median `sigma_tE` ratio
came out 1.1503 [1.0805, 1.2170], i.e. the satellite apparently making the TIMESCALE forecast
15% worse. That is impossible for matched events -- a Fisher forecast at fixed parameters and
epochs does not degrade because the geometry changed -- so the 15% measured how much the two
DETECTED populations differ, which is larger than the ~6% effect being sought. Hence pairing.

GATING. A ratio needs both forecasts to exist: `okA_sat == 1 and okA_nosat == 1`, with -1.0
excluded as the not-measured sentinel. A sigma that never inverted is not a large sigma.

THE BUILT-IN CONTROL. Events with `nepR_pk == 0` have no Roman epochs near the peak, so the
satellite cannot act on them and their ratio must be 1. They are kept deliberately and reported
as a control; if they ever drift from 1, the measurement is wrong.
"""
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Temporal-baseline gain for H3c: median sigma_joint/sigma_Roman for piE on in-gap events,
# from Step F2 (DEVIATIONS.md 24.2). A DIFFERENT physical effect -- Rubin filling Roman's
# season gaps in time, not two observers separated in space. H3c exists to put the two on one
# axis at the right scale, not to declare a winner.
TEMPORAL_GAIN = {"10-30 d": 0.31, "30-100 d": 0.80, "100-300 d": 0.95, "> 300 d": 0.984}
TE_EDGES = [10.0, 30.0, 100.0, 300.0, np.inf]
TE_LABELS = list(TEMPORAL_GAIN.keys())
BG = "#fcfcfb"
COLS = ("lon lat tE u0 piE du_sat okA_sat okA_nosat sigtE_sat sigtE_nosat "
        "sigpiE_sat sigpiE_nosat nepL_pk nepR_pk w_area").split()


def load(path):
    df = pd.read_csv(path, sep=r"\s+", comment="#", names=COLS, engine="python")
    for c in COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pairs", help="h3_pair.dat from a --pair-satellite run")
    ap.add_argument("--out-prefix", default="figures/h3")
    a = ap.parse_args()

    df = load(a.pairs)
    print(f"== {len(df):,} detected events in {a.pairs}")

    ok = (df.okA_sat == 1) & (df.okA_nosat == 1) & \
         (df.sigpiE_sat > 0) & (df.sigpiE_nosat > 0) & \
         (df.sigtE_sat > 0) & (df.sigtE_nosat > 0)
    d = df[ok].copy()
    print(f"   both forecasts inverted (okA on both sides): {len(d):,}"
          f"   ({len(d)/max(len(df),1):.1%})")
    if len(d) == 0:
        sys.exit("no events with both forecasts; nothing to compare")

    d["ratio"] = d.sigpiE_sat / d.sigpiE_nosat
    d["ratio_tE"] = d.sigtE_sat / d.sigtE_nosat
    d["tb"] = pd.cut(d.tE, TE_EDGES, labels=TE_LABELS, right=False)

    covered = d[d.nepR_pk > 0]
    blind = d[d.nepR_pk == 0]

    # ---- the control comes first. If it is not ~1 the rest is meaningless.
    print("\n== CONTROL: events with no Roman epochs near the peak (nepR_pk == 0)")
    print("   the satellite cannot act on these, so the ratio must be 1")
    if len(blind):
        print(f"   n = {len(blind):,}   median sigma_piE ratio = {blind.ratio.median():.6f}")
        print(f"                       median sigma_tE  ratio = {blind.ratio_tE.median():.6f}")
        print(f"   fraction within 1% of unity: {(blind.ratio.between(0.99,1.01)).mean():.1%}")
    else:
        print("   none in this sample")

    print("\n== events where Roman covers the peak (nepR_pk > 0)")
    print(f"   n = {len(covered):,}")
    if len(covered) == 0:
        sys.exit("no Roman-covered events")
    q = covered.ratio.quantile([0.05, 0.25, 0.5, 0.75, 0.95])
    print(f"   sigma_piE ratio  median {covered.ratio.median():.4f}"
          f"   quartiles {q[0.25]:.4f} / {q[0.75]:.4f}"
          f"   5-95% {q[0.05]:.4f} / {q[0.95]:.4f}")
    print(f"   sigma_tE  ratio  median {covered.ratio_tE.median():.4f}   (not a control here:")
    print("                     tE and piE are correlated, so the geometry moves both)")
    better = (covered.ratio < 0.99).mean()
    big = (covered.ratio < 0.5).mean()
    print(f"   improved by >1%  : {better:.1%}")
    print(f"   improved by >2x  : {big:.1%}")

    figures(d, covered, blind, a.out_prefix)
    verdict(covered, blind)
    return 0


def figures(d, covered, blind, prefix):
    # ================================================================= H3a
    fig, ax = plt.subplots(figsize=(8.6, 5.4), facecolor=BG)
    ax.set_facecolor(BG)
    colors = plt.cm.viridis(np.linspace(0.12, 0.88, len(TE_LABELS)))
    for lab, col in zip(TE_LABELS, colors):
        g = covered[covered.tb == lab]
        if len(g) < 5:
            continue
        ax.scatter(g.du_sat, g.ratio, s=7, alpha=0.28, color=col, lw=0, label=f"tE {lab}")
    if len(blind):
        ax.scatter(blind.du_sat, blind.ratio, s=7, alpha=0.30, color="#b02020", lw=0,
                   marker="x", label="no Roman epochs at peak (control)")
    ax.axhline(1.0, color="#444", lw=1.0, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta u_{\rm sat} = \pi_E\,D_\perp/{\rm AU}$   [$\theta_E$]")
    ax.set_ylabel(r"$\sigma_{\pi_E}$(Roman at L2) / $\sigma_{\pi_E}$(Roman at Earth)")
    ax.set_title("H3a  Where the satellite baseline acts\n"
                 "per-event ratio, same event characterised both ways", fontsize=11)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.grid(alpha=0.25, lw=0.6)
    fig.text(0.01, -0.05,
             "Below 1 means Roman's displacement to L2 tightens the parallax forecast. Red "
             "crosses are the control: events with no Roman\nepochs near the peak, where the "
             "satellite cannot act and the ratio must sit at 1. Pairing is exact -- one run, "
             "two Fisher\nevaluations per event (DEVIATIONS.md 35).", fontsize=7.5, color="#444")
    fig.savefig(f"{prefix}a_where.png", dpi=160, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"\nwrote {prefix}a_where.png")

    # ================================================================= H3b
    both = covered[(covered.nepL_pk > 0) & (covered.nepR_pk > 0)]
    print(f"H3b contemporaneous coverage (both surveys at peak): {len(both):,}")
    te_e = np.geomspace(5, 400, 9)
    u0_e = np.linspace(0, 1.2, 9)
    M = np.full((len(u0_e) - 1, len(te_e) - 1), np.nan)
    for i in range(len(u0_e) - 1):
        for j in range(len(te_e) - 1):
            g = both[(both.u0 >= u0_e[i]) & (both.u0 < u0_e[i + 1]) &
                     (both.tE >= te_e[j]) & (both.tE < te_e[j + 1])]
            if len(g) >= 8:
                M[i, j] = g.ratio.median()

    fig, ax = plt.subplots(figsize=(8.2, 5.2), facecolor=BG)
    ax.set_facecolor(BG)
    if np.isfinite(M).any():
        span = max(np.nanmax(np.abs(np.log10(M))), 0.05)
        im = ax.pcolormesh(te_e, u0_e, np.log10(M), cmap="RdBu_r",
                           vmin=-span, vmax=span, shading="flat")
        cb = fig.colorbar(im, ax=ax)
        cb.set_label(r"$\log_{10}$ median ratio;  blue = satellite helps")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isfinite(M[i, j]):
                    ax.text(np.sqrt(te_e[j] * te_e[j + 1]), (u0_e[i] + u0_e[i + 1]) / 2,
                            f"{M[i, j]:.2f}", ha="center", va="center", fontsize=6.2)
    ax.set_xscale("log")
    ax.set_xlabel(r"$t_E$  [d]")
    ax.set_ylabel(r"$u_0$")
    ax.set_title("H3b  The corner, if there is one\n"
                 "median per-event ratio, contemporaneous coverage only", fontsize=11)
    fig.text(0.01, -0.05,
             "Cells need >=8 paired events; blank cells are unpopulated, not null results.",
             fontsize=7.5, color="#444")
    fig.savefig(f"{prefix}b_corner.png", dpi=160, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"wrote {prefix}b_corner.png")

    # ================================================================= H3c
    fig, ax = plt.subplots(figsize=(8.8, 5.0), facecolor=BG)
    ax.set_facecolor(BG)
    x = np.arange(len(TE_LABELS))
    w = 0.36
    sat = [covered[covered.tb == l].ratio.median() if (covered.tb == l).sum() >= 5 else np.nan
           for l in TE_LABELS]
    ax.bar(x - w / 2, sat, w, color="#3b6ea5",
           label="satellite baseline  (Roman at L2 vs at Earth), this work")
    ax.bar(x + w / 2, [TEMPORAL_GAIN[l] for l in TE_LABELS], w, color="#c0703a",
           label="temporal baseline  (Rubin filling Roman's gaps), Step F2")
    ax.axhline(1.0, color="#444", lw=1.0, ls="--")
    ax.set_xticks(x, TE_LABELS)
    ax.set_yscale("log")
    ax.set_xlabel(r"$t_E$ bin")
    ax.set_ylabel(r"median $\sigma_{\pi_E}$ ratio   (lower = bigger gain)")
    ax.set_title("H3c  Two different effects, at the same scale", fontsize=11)
    ax.legend(frameon=False, fontsize=8.5)
    ax.grid(alpha=0.25, axis="y", lw=0.6)
    for xi, v in zip(x - w / 2, sat):
        if np.isfinite(v):
            ax.text(xi, v * 1.06, f"{v:.3f}", ha="center", fontsize=8)
    for xi, l in zip(x + w / 2, TE_LABELS):
        ax.text(xi, TEMPORAL_GAIN[l] * 1.06, f"{TEMPORAL_GAIN[l]:.2f}", ha="center", fontsize=8)
    fig.text(0.01, -0.06,
             "DIFFERENT PHYSICAL EFFECTS. Pairing the bars compares their size, not their "
             "merit: one is two observers separated in\nspace, the other is one observer's "
             "gaps filled in time. The temporal numbers are in-gap medians from Step F2 "
             "(Deviations 24.2).", fontsize=7.5, color="#444")
    fig.savefig(f"{prefix}c_honest.png", dpi=160, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"wrote {prefix}c_honest.png")


def verdict(covered, blind):
    print("\n" + "=" * 74)
    print("THE RESULT, STATED PLAINLY")
    print("=" * 74)
    med = covered.ratio.median()
    better = (covered.ratio < 0.99).mean()
    print(f"  Roman-covered events: n = {len(covered):,}, median sigma_piE ratio {med:.4f}")
    print(f"  improved by more than 1%: {better:.1%} of them")
    for lab in TE_LABELS:
        g = covered[covered.tb == lab]
        if len(g) >= 5:
            print(f"    tE {lab:10s} n={len(g):6,}  median {g.ratio.median():.4f}"
                  f"  best 5% {g.ratio.quantile(0.05):.4f}")
    if len(blind):
        print(f"  control (no Roman epochs at peak): median {blind.ratio.median():.6f}"
              f" on n={len(blind):,}  <- must be 1")
    print()
    if better < 0.01:
        print("  THIS IS A NULL. Fewer than 1% of Roman-covered events gain anything")
        print("  measurable from Roman's displacement to L2.")
    else:
        strong = covered[covered.ratio < 0.5]
        print(f"  NOT a null. {better:.1%} of Roman-covered events improve by >1%, and")
        print(f"  {len(strong):,} ({len(strong)/len(covered):.1%}) improve by more than 2x.")
        if len(strong):
            print(f"  Those events sit at median tE {strong.tE.median():.1f} d, "
                  f"u0 {strong.u0.median():.3f}, piE {strong.piE.median():.3f}, "
                  f"du_sat {strong.du_sat.median():.4f}.")


if __name__ == "__main__":
    sys.exit(main())
