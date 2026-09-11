#!/usr/bin/env python3
"""Step H3 -- what Roman's displacement from Earth buys for the parallax forecast.

THE EXPERIMENT, AND WHY IT IS NOT THE ONE THE PLAN DESCRIBES
------------------------------------------------------------
The plan asks for two runs differing only in `L2_OFFSET_AU`, matched event by row, "every
event appears in both". That is impossible here, and the reason is physics rather than
bookkeeping (DEVIATIONS.md 35). The RNG is one mt19937_64 stream shared across the scan, and
the per-event path draws `RandR(0,100)` *conditionally* on `acceptRubin or acceptRoman`. Moving
Roman off the Earth-Sun line changes the trajectory, changes what is detectable, and so changes
whether that draw is taken. The two runs are byte-identical for 1,740,091 events across 459
sightlines -- and then fork at the first footprint sightline the scan reaches, which is exactly
the first place the satellite offset can do anything. They stay paired only where the effect is
identically zero.

So this is a POPULATION comparison, not a paired one. Both runs are gated identically, cut to
events where Roman actually contributes, and compared as distributions inside bins. The ratios
plotted are RATIOS OF BINNED MEDIANS, never per-event ratios; a reader who takes the scatter
for per-event dispersion will over-read it, which is why the captions say so.

The design that would restore true pairing -- call FisherM twice per detected event, once with
the offset and once with it zeroed, inside a single run -- is recorded in DEVIATIONS.md 35 and
OPEN_ITEMS.md. It costs about one extra Fisher call per detection.

GATING. Events must be jointly detected AND jointly characterised: `detJ == 1 and okA_J == 1`,
with -1.0 excluded as the not-measured sentinel. A sigma that was never measured is not a
large sigma. Applied in the extractor (h3_extract.py) and re-asserted here.

SELECTION. `nepR_pk > 0` -- Roman has epochs near the peak. This is the symmetric cut: it is a
property of the event and the schedule, identical in meaning in both runs, and it is the only
place the satellite offset can act. `du_sat` cannot be used to select, because it is
identically 0 in the no-satellite run by construction.

BINNING VARIABLE. `piE`, not `du_sat`, for the same reason: `du_sat` is 0 in one run. They are
near-equivalent -- du_sat = piE * D_perp/AU and D_perp/AU runs 0.87-0.99 of L2_OFFSET_AU, so
du_sat ~ 0.01 * piE to within ~13% -- and the top axis of H3a is labelled in du_sat accordingly.
"""
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L2_OFFSET_AU = 0.0100267

# Temporal-baseline gain, for H3c. Median sigma_joint/sigma_Roman for piE on in-gap events,
# measured in Step F2 and recorded in DEVIATIONS.md 24.2. This is a DIFFERENT physical effect
# from satellite parallax -- Rubin filling Roman's season gaps in time, not two observers
# separated in space -- and H3c exists to put the two on one axis at the right scale.
TEMPORAL_GAIN = {"10-30 d": 0.31, "30-100 d": 0.80, "100-300 d": 0.95, "> 300 d": 0.984}
TE_EDGES = [10.0, 30.0, 100.0, 300.0, np.inf]
TE_LABELS = list(TEMPORAL_GAIN.keys())

BG = "#fcfcfb"


def load(path):
    df = pd.read_csv(path, sep="\t")
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    n0 = len(df)
    df = df[(df["sigpiE_J"] > 0) & (df["sigtE_J"] > 0)].copy()
    if len(df) != n0:
        print(f"  dropped {n0 - len(df)} rows with a non-positive sigma (sentinel)")
    return df


def boot_median_ratio(a, b, n=400, rng=None):
    """Bootstrap CI for median(a)/median(b) with a and b INDEPENDENT samples.

    They are independent here -- the two runs stopped being paired at the first footprint
    sightline -- so the ratio's uncertainty is not the paired one and must not be computed as
    though the events matched.
    """
    rng = rng or np.random.default_rng(42)
    if len(a) < 20 or len(b) < 20:
        return np.nan, np.nan, np.nan
    ra = rng.choice(a, (n, len(a)), replace=True)
    rb = rng.choice(b, (n, len(b)), replace=True)
    r = np.median(ra, axis=1) / np.median(rb, axis=1)
    return np.median(a) / np.median(b), np.percentile(r, 16), np.percentile(r, 84)


def te_bin(df):
    return pd.cut(df["tE"], TE_EDGES, labels=TE_LABELS, right=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sat", required=True, help="cached subset of the satellite-parallax run")
    ap.add_argument("--nosat", required=True, help="cached subset of the twin (L2_OFFSET_AU=0)")
    ap.add_argument("--out-prefix", default="figures/h3")
    a = ap.parse_args()

    print("== loading")
    S, N = load(a.sat), load(a.nosat)
    print(f"  satellite run : {len(S):,} detected+characterised")
    print(f"  twin (no sat) : {len(N):,} detected+characterised")

    # ---- the satellite observable must be non-zero in one run and zero in the other.
    # If this ever fails the twin's flag did nothing and the whole comparison is empty.
    fs = float((S["du_sat"] > 0).mean())
    fn = float((N["du_sat"] != 0).mean())
    print(f"  du_sat > 0 in satellite run : {fs:6.1%}")
    print(f"  du_sat != 0 in twin         : {fn:6.1%}   (must be 0.0%)")
    if fn != 0.0:
        sys.exit("ABORT: the twin carries a non-zero satellite observable; --no-satellite-parallax did nothing")
    if fs == 0.0:
        sys.exit("ABORT: the satellite run carries no satellite observable")

    # ---- selection: Roman contributes near the peak
    S = S[S["nepR_pk"] > 0].copy()
    N = N[N["nepR_pk"] > 0].copy()
    print(f"\n== Roman contributes near peak (nepR_pk > 0)")
    print(f"  satellite run : {len(S):,}")
    print(f"  twin          : {len(N):,}")
    if len(S) == 0 or len(N) == 0:
        sys.exit("ABORT: no events with Roman coverage near the peak")

    # ---- are the two populations comparable? An unpaired comparison is only meaningful if
    # the underlying event populations match; these are the same draws from the same model,
    # so they should agree closely on quantities the satellite cannot influence.
    print("\n== population comparability (quantities the satellite should NOT shift)")
    for c in ("tE", "u0", "piE"):
        ms, mn = S[c].median(), N[c].median()
        print(f"  median {c:4s}  sat {ms:10.4f}   nosat {mn:10.4f}   "
              f"ratio {ms / mn:6.4f}")

    # ---- headline
    r, lo, hi = boot_median_ratio(S["sigpiE_J"].to_numpy(), N["sigpiE_J"].to_numpy())
    rt, tlo, thi = boot_median_ratio(S["sigtE_J"].to_numpy(), N["sigtE_J"].to_numpy())
    print("\n== headline, all Roman-covered detections")
    print(f"  median sigma_piE  sat/nosat = {r:.4f}  [{lo:.4f}, {hi:.4f}]  (68% boot)")
    print(f"  median sigma_tE   sat/nosat = {rt:.4f}  [{tlo:.4f}, {thi:.4f}]   <- control, expect ~1")

    figures(S, N, a.out_prefix, r, lo, hi)
    return 0


def figures(S, N, prefix, r_all, lo_all, hi_all):
    rng = np.random.default_rng(7)

    # =====================================================================  H3a
    # Where the effect lives: ratio of binned medians against du_sat (via piE).
    pi_edges = np.geomspace(max(1e-3, S["piE"].quantile(0.01)), S["piE"].quantile(0.999), 13)
    S["tb"], N["tb"] = te_bin(S), te_bin(N)

    fig, ax = plt.subplots(figsize=(8.4, 5.4), facecolor=BG)
    ax.set_facecolor(BG)
    colors = plt.cm.viridis(np.linspace(0.12, 0.88, len(TE_LABELS)))
    rows = []
    for lab, col in zip(TE_LABELS, colors):
        xs, ys, els, ehs = [], [], [], []
        for i in range(len(pi_edges) - 1):
            lo_e, hi_e = pi_edges[i], pi_edges[i + 1]
            sa = S[(S["tb"] == lab) & (S["piE"] >= lo_e) & (S["piE"] < hi_e)]["sigpiE_J"]
            nb = N[(N["tb"] == lab) & (N["piE"] >= lo_e) & (N["piE"] < hi_e)]["sigpiE_J"]
            if len(sa) < 30 or len(nb) < 30:
                continue
            rr, rlo, rhi = boot_median_ratio(sa.to_numpy(), nb.to_numpy(), 200, rng)
            xs.append(np.sqrt(lo_e * hi_e)); ys.append(rr)
            els.append(rr - rlo); ehs.append(rhi - rr)
            rows.append((lab, np.sqrt(lo_e * hi_e), rr, rlo, rhi, len(sa), len(nb)))
        if xs:
            ax.errorbar(xs, ys, yerr=[els, ehs], marker="o", ms=4.5, lw=1.4,
                        capsize=2.5, color=col, label=f"tE {lab}")
    ax.axhline(1.0, color="#444", lw=1.0, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\pi_E$   (satellite separation $\Delta u_{\rm sat}\simeq 0.010\,\pi_E$)")
    ax.set_ylabel(r"median $\sigma_{\pi_E}$  (with L2) / (without)")
    ax.set_title("H3a  Where the satellite baseline acts\n"
                 "ratio of BINNED MEDIANS -- not a per-event ratio", fontsize=11)
    ax2 = ax.secondary_xaxis("top", functions=(lambda x: x * L2_OFFSET_AU,
                                               lambda x: x / L2_OFFSET_AU))
    ax2.set_xlabel(r"$\Delta u_{\rm sat}$  [$\theta_E$]", fontsize=9)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25, lw=0.6)
    fig.text(0.01, -0.04,
             "Below 1 means the satellite baseline improves the parallax forecast. Bars are 68% "
             "bootstrap on the ratio of medians,\ncomputed for INDEPENDENT samples: the two runs "
             "stop being event-matched at the first footprint sightline (Deviations 35).",
             fontsize=7.5, color="#444")
    fig.savefig(f"{prefix}a_where.png", dpi=160, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"\nwrote {prefix}a_where.png  ({len(rows)} populated bins)")

    # =====================================================================  H3b
    # The corner: (tE, u0) with contemporaneous coverage from both observatories.
    Sc = S[(S["nepL_pk"] > 0) & (S["nepR_pk"] > 0)]
    Nc = N[(N["nepL_pk"] > 0) & (N["nepR_pk"] > 0)]
    print(f"H3b contemporaneous coverage: sat {len(Sc):,}, nosat {len(Nc):,}")
    te_e = np.geomspace(5, 400, 9)
    u0_e = np.linspace(0, 1.0, 9)
    M = np.full((len(u0_e) - 1, len(te_e) - 1), np.nan)
    Cnt = np.zeros_like(M)
    for i in range(len(u0_e) - 1):
        for j in range(len(te_e) - 1):
            sa = Sc[(Sc["u0"] >= u0_e[i]) & (Sc["u0"] < u0_e[i + 1]) &
                    (Sc["tE"] >= te_e[j]) & (Sc["tE"] < te_e[j + 1])]["sigpiE_J"]
            nb = Nc[(Nc["u0"] >= u0_e[i]) & (Nc["u0"] < u0_e[i + 1]) &
                    (Nc["tE"] >= te_e[j]) & (Nc["tE"] < te_e[j + 1])]["sigpiE_J"]
            if len(sa) >= 30 and len(nb) >= 30:
                M[i, j] = sa.median() / nb.median()
                Cnt[i, j] = min(len(sa), len(nb))

    fig, ax = plt.subplots(figsize=(8.0, 5.2), facecolor=BG)
    ax.set_facecolor(BG)
    span = np.nanmax(np.abs(np.log10(M))) if np.isfinite(M).any() else 0.1
    span = max(span, 0.02)
    im = ax.pcolormesh(te_e, u0_e, np.log10(M), cmap="RdBu_r",
                       vmin=-span, vmax=span, shading="flat")
    ax.set_xscale("log")
    ax.set_xlabel(r"$t_E$  [d]")
    ax.set_ylabel(r"$u_0$")
    ax.set_title("H3b  The corner, if there is one\n"
                 r"$\log_{10}$ of median $\sigma_{\pi_E}$ ratio, contemporaneous coverage only",
                 fontsize=11)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(r"$\log_{10}$ (with L2 / without);  blue = satellite helps")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(M[i, j]):
                ax.text(np.sqrt(te_e[j] * te_e[j + 1]), (u0_e[i] + u0_e[i + 1]) / 2,
                        f"{M[i, j]:.3f}", ha="center", va="center", fontsize=6.2,
                        color="#222")
    fig.text(0.01, -0.05,
             "Cells need >=30 characterised events on BOTH sides; blank cells are unpopulated, "
             "not null results.\nNumbers are ratios of medians over independent samples "
             "(Deviations 35).", fontsize=7.5, color="#444")
    fig.savefig(f"{prefix}b_corner.png", dpi=160, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"wrote {prefix}b_corner.png")

    # =====================================================================  H3c
    # Honest comparison: satellite baseline beside temporal baseline, same axis.
    sat_by_te, err_by_te = [], []
    for lab in TE_LABELS:
        sa = S[S["tb"] == lab]["sigpiE_J"].to_numpy()
        nb = N[N["tb"] == lab]["sigpiE_J"].to_numpy()
        rr, rlo, rhi = boot_median_ratio(sa, nb, 300, rng)
        sat_by_te.append(rr)
        err_by_te.append((rr - rlo, rhi - rr))

    fig, ax = plt.subplots(figsize=(8.6, 5.0), facecolor=BG)
    ax.set_facecolor(BG)
    x = np.arange(len(TE_LABELS))
    w = 0.36
    yerr = np.array(err_by_te).T
    ax.bar(x - w / 2, sat_by_te, w, yerr=yerr, capsize=3, color="#3b6ea5",
           label="satellite baseline  (Roman at L2 vs at Earth)")
    ax.bar(x + w / 2, [TEMPORAL_GAIN[l] for l in TE_LABELS], w, color="#c0703a",
           label="temporal baseline  (Rubin filling Roman's gaps, F2)")
    ax.axhline(1.0, color="#444", lw=1.0, ls="--")
    ax.set_xticks(x, TE_LABELS)
    ax.set_xlabel(r"$t_E$ bin")
    ax.set_ylabel(r"median $\sigma_{\pi_E}$ ratio  (lower = bigger gain)")
    ax.set_title("H3c  Two different effects, at the same scale", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25, axis="y", lw=0.6)
    for xi, v in zip(x - w / 2, sat_by_te):
        if np.isfinite(v):
            ax.text(xi, v + 0.012, f"{v:.3f}", ha="center", fontsize=8)
    for xi, l in zip(x + w / 2, TE_LABELS):
        ax.text(xi, TEMPORAL_GAIN[l] + 0.012, f"{TEMPORAL_GAIN[l]:.2f}", ha="center", fontsize=8)
    fig.text(0.01, -0.06,
             "These are DIFFERENT PHYSICAL EFFECTS and the pairing of bars is a scale "
             "comparison, not a competition. The temporal\nnumbers are in-gap medians from "
             "Step F2 (Deviations 24.2); the satellite numbers are all Roman-covered "
             "detections here.", fontsize=7.5, color="#444")
    fig.savefig(f"{prefix}c_honest.png", dpi=160, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"wrote {prefix}c_honest.png")

    # =====================================================================  the null
    print("\n" + "=" * 72)
    print("THE RESULT, STATED PLAINLY")
    print("=" * 72)
    print(f"  All Roman-covered joint detections: median sigma_piE ratio "
          f"{r_all:.4f} [{lo_all:.4f}, {hi_all:.4f}]")
    best = min((v for v in sat_by_te if np.isfinite(v)), default=np.nan)
    print("  by tE bin (satellite):")
    for lab, v, e in zip(TE_LABELS, sat_by_te, err_by_te):
        if np.isfinite(v):
            print(f"    {lab:10s} {v:7.4f}   -{e[0]:.4f}/+{e[1]:.4f}")
    if np.isfinite(M).any():
        k = np.unravel_index(np.nanargmin(M), M.shape)
        print(f"  most favourable (tE,u0) cell: ratio {M[k]:.4f} "
              f"at tE~{np.sqrt(te_e[k[1]]*te_e[k[1]+1]):.0f} d, u0~{(u0_e[k[0]]+u0_e[k[0]+1])/2:.2f}"
              f"  (n>={int(Cnt[k])})")
    print()
    if np.isfinite(best) and best > 0.99:
        print("  READ AS A NULL. No tE bin shows a median improvement better than 1%.")
    else:
        print(f"  Largest binned improvement: {(1 - best) * 100:.1f}% (median), in a single tE bin.")
    print("  Caveat that belongs in every caption: these are ratios of binned medians over")
    print("  INDEPENDENT samples, not per-event ratios. See DEVIATIONS.md 35.")


if __name__ == "__main__":
    sys.exit(main())
