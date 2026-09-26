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
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R

# Temporal-baseline gain for H3c: median sigma_joint/sigma_Roman for piE on in-gap events
# (t0zone == 1, Roman-covered), EVENT-RATE WEIGHTED, from the post-extinction-fix bulge run
# (runs/prod_bulge_20260924; figures/wp_20260926/w1.log, Deviation 55). The pre-fix v3 values were
# 0.433 / 0.936 / 0.978 / 0.991 (unweighted); the fix all but removed the pooled in-gap gain, which
# now survives only for peaks deep inside a gap (dt_edge > 45 d: tE ratio 0.077). A DIFFERENT
# physical effect -- Rubin filling Roman's season gaps in time, not two observers separated in
# space. H3c exists to put the two on one axis at the right scale, not to declare a winner.
# The in-season control for the same quantity is 0.999-1.000 in every bin.
TEMPORAL_GAIN = {"10-30 d": 0.984, "30-100 d": 0.999, "100-300 d": 1.000, "> 300 d": 1.000}
TE_EDGES = [10.0, 30.0, 100.0, 300.0, np.inf]
TE_LABELS = list(TEMPORAL_GAIN.keys())
BG = "#fcfcfb"
# The current h3_pair.dat layout. Order matters absolutely: the file's header line begins with
# a lone '#', so pandas is told comment="#" and the names are supplied here instead. Get this
# list out of step with Bulge_LSST.cpp and every column silently shifts by one -- which has
# already happened once in this project and cost a debugging session (PROGRESS.md traps).
COLS_LEGACY = ("lon lat tE u0 piE tetE du_sat "
               "okA_sat okA_nosat okB_sat okB_nosat "
               "sigtE_sat sigtE_nosat sigpiE_sat sigpiE_nosat "
               "sigpiER_sat sigpiER_nosat sigtetE_sat sigtetE_nosat "
               "sigpiEb_sat sigpiEb_nosat relMl_sat relMl_nosat "
               "condA_sat condA_nosat condB_sat condB_nosat "
               "nepL_pk nepR_pk w_area").split()

# Ml, Dl, Ds and Vt were appended on 2026-09-17 so a paired file can carry the event-rate
# weight (Deviation 41): W needs sqrt(Ml)*Vt*Z(Ds), and none of it can be reconstructed from
# the legacy columns -- theta_E and pi_E give Ml and theta_E/tE gives mu_rel, but
# pi_rel = 1/Dl - 1/Ds is one equation in two unknowns. A legacy file is therefore readable
# but NOT weightable, and this script says so rather than quietly reporting raw medians.
COLS = COLS_LEGACY + ["Ml", "Dl", "Ds", "Vt"]

# The superseded layout, recognised only so it can be refused by name.
COLS_OLD_15 = 15


def load(path):
    probe = pd.read_csv(path, sep=r"\s+", comment="#", header=None, nrows=1,
                        engine="python")
    n = probe.shape[1]
    if n == COLS_OLD_15:
        sys.exit(
            f"{path} has the superseded 15-column layout.\n"
            "That file was written before the derivative-reference fix (DEVIATIONS.md 36) and\n"
            "every forecast in it is corrupted on the no-satellite side. Refusing to read it\n"
            "rather than reporting numbers from it. Re-run with --pair-satellite on a binary\n"
            "built from 3a88180 or later.")
    if n == len(COLS):
        names, weightable = COLS, True
    elif n == len(COLS_LEGACY):
        names, weightable = COLS_LEGACY, False
    else:
        sys.exit(f"{path} has {n} columns; this script expects {len(COLS)} (or "
                 f"{len(COLS_LEGACY)} for a pre-2026-09-17 file). "
                 "If Bulge_LSST.cpp changed the row, update COLS to match it.")
    df = pd.read_csv(path, sep=r"\s+", comment="#", names=names, engine="python")
    for c in names:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(), weightable


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pairs", help="h3_pair.dat from a --pair-satellite run")
    ap.add_argument("--out-prefix", default="figures/h3")
    ap.add_argument("--map", default=None,
                    help="MapLMC5.dat FROM THE PAIRED RUN -- needed for the event-rate weight")
    ap.add_argument("--log", action="append", default=[],
                    help="the paired run's log(s), for sightlines its map file lost")
    ap.add_argument("--unweighted", action="store_true",
                    help="deliberately report raw-sample medians, with no event-rate weight")
    a = ap.parse_args()

    df, weightable = load(a.pairs)
    print(f"== {len(df):,} detected events in {a.pairs}")

    if not weightable and not a.unweighted:
        sys.exit(
            f"{a.pairs} is a pre-2026-09-17 paired file: it has no Ml/Dl/Ds/Vt columns, so the\n"
            "event-rate weight cannot be computed from it and every median below would be a\n"
            "statistic of the raw sample (Deviation 41: that sample over-represents long, slow,\n"
            "massive lenses roughly tenfold).\n"
            "Either re-run --pair-satellite with a binary built from 2026-09-17 or later, which\n"
            "writes the four columns, or pass --unweighted to say deliberately that you want the\n"
            "raw-sample numbers.")

    # Weights are attached to the FULL frame before any subsetting, so `covered`, `blind` and
    # every later selection carry their own weights with them.
    w, wlabel = R.attach_weight(df, a.map, a.log, a.unweighted or not weightable)
    df = df.assign(W=w)
    print(f"   weighting: {wlabel}")

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
        print(f"   n = {len(blind):,}   median sigma_piE ratio = "
              f"{R.weighted_median(blind.ratio, blind.W):.6f}")
        print(f"                       median sigma_tE  ratio = "
              f"{R.weighted_median(blind.ratio_tE, blind.W):.6f}")
        print(f"   fraction within 1% of unity: "
              f"{R.weighted_fraction(blind.ratio.between(0.99, 1.01), blind.W):.1%}")
    else:
        print("   none in this sample")

    print("\n== events where Roman covers the peak (nepR_pk > 0)")
    print(f"   n = {len(covered):,}")
    if len(covered) == 0:
        sys.exit("no Roman-covered events")
    q = {p: R.weighted_quantile(covered.ratio, covered.W, p)
         for p in (0.05, 0.25, 0.5, 0.75, 0.95)}
    print(f"   N_eff {R.kish_neff(covered.W):,.0f}")
    print(f"   sigma_piE ratio  median {q[0.5]:.4f}"
          f"   quartiles {q[0.25]:.4f} / {q[0.75]:.4f}"
          f"   5-95% {q[0.05]:.4f} / {q[0.95]:.4f}")
    print(f"   sigma_tE  ratio  median "
          f"{R.weighted_median(covered.ratio_tE, covered.W):.4f}   (not a control here:")
    print("                     tE and piE are correlated, so the geometry moves both)")
    better = R.weighted_fraction(covered.ratio < 0.99, covered.W)
    big = R.weighted_fraction(covered.ratio < 0.5, covered.W)
    print(f"   improved by >1%  : {better:.1%}")
    print(f"   improved by >2x  : {big:.1%}")

    astrometry(df, covered, blind)
    conditioning(covered)

    figures(d, covered, blind, a.out_prefix)
    passed = validate(df, covered, blind)
    verdict(covered, blind, passed)
    return 0


def astrometry(df, covered, blind):
    """What moving the observer does to the ASTROMETRIC forecast, and to the lens mass.

    Different physics from the photometric side, and worth separating. The magnification
    depends only on |u|, so the satellite baseline reaches piE through a change in separation.
    The centroid deflection theta_E u/(u^2+2) is a VECTOR, so moving the observer changes its
    direction as well as its magnitude, and the astrometric matrix carries both sigma(theta_E)
    and a route to piE independent of the photometric one.

    The lens mass is the quantity that matters: Ml = theta_E/(kappa piE) needs one observable
    from each matrix, so it is the only place where a change in either becomes a change in the
    science rather than a change in a nuisance parameter.
    """
    print("\n== ASTROMETRY: what the observer move does to theta_E and to the mass")
    okb = (df.okB_sat == 1) & (df.okB_nosat == 1) & \
          (df.sigtetE_sat > 0) & (df.sigtetE_nosat > 0)
    b = df[okb]
    print(f"   both astrometric matrices inverted: {len(b):,} of {len(df):,}")
    if len(b) == 0:
        print("   none; nothing to report")
        return
    for lab, sel in (("Roman covers the peak", b[b.nepR_pk > 0]),
                     ("CONTROL, no Roman epochs near peak", b[b.nepR_pk == 0])):
        if not len(sel):
            print(f"   {lab}: none")
            continue
        r = (sel.sigtetE_sat / sel.sigtetE_nosat)
        print(f"   {lab}: n = {len(sel):,}   median sigma(theta_E) ratio = "
              f"{R.weighted_median(r, sel.W):.6f}")

    m = (df.relMl_sat > 0) & (df.relMl_nosat > 0) & (df.nepR_pk > 0)
    if m.sum():
        rm = (df.loc[m, "relMl_sat"] / df.loc[m, "relMl_nosat"])
        wm = df.loc[m, "W"]
        print(f"   lens mass: n = {int(m.sum()):,}   median relMl ratio = "
              f"{R.weighted_median(rm, wm):.6f}")
        print(f"              improved by >1%: {R.weighted_fraction(rm < 0.99, wm):.1%}")

    mr = (covered.sigpiER_sat > 0) & (covered.sigpiER_nosat > 0)
    if mr.sum():
        rr = covered.loc[mr, "sigpiER_sat"] / covered.loc[mr, "sigpiER_nosat"]
        print(f"   Roman ALONE, sigma(piE): n = {int(mr.sum()):,}   median ratio = "
              f"{R.weighted_median(rr, covered.loc[mr, 'W']):.6f}")
        print("              the L2 offset is Roman's geometry, so an effect must show here first")


def conditioning(covered):
    """Answer the near-degeneracy question from the data instead of hypothesising it.

    OPEN_ITEMS recorded a guess that the satScale = 0 configuration might leave a
    near-degeneracy that clears the okA condition-number gate while leaving marginalised errors
    unstable. That guess was NOT the cause of the original H3 failure -- a stale derivative
    reference was, DEVIATIONS 36 -- but the question is real and the columns are now here.
    """
    print("\n== CONDITIONING of the two matrices (the OPEN_ITEMS question, answered)")
    for lab, a, b in (("photometric (condA)", "condA_sat", "condA_nosat"),
                      ("astrometric (condB)", "condB_sat", "condB_nosat")):
        m = (covered[a] > 0) & (covered[b] > 0)
        if not m.sum():
            print(f"   {lab}: no valid pairs")
            continue
        ca, cb, cw = covered.loc[m, a], covered.loc[m, b], covered.loc[m, "W"]
        print(f"   {lab}: n = {int(m.sum()):,}")
        print(f"      satellite    median {R.weighted_median(ca, cw):.4g}   "
              f"95th {R.weighted_quantile(ca, cw, 0.95):.4g}")
        print(f"      no-satellite median {R.weighted_median(cb, cw):.4g}   "
              f"95th {R.weighted_quantile(cb, cw, 0.95):.4g}")
        print(f"      median ratio no-sat/sat = {R.weighted_median(cb/ca, cw):.4f}"
              "   (far from 1 would mean the two geometries are not equally conditioned)")


def validate(df, c, b):
    """Checks the measurement must pass before any number from it can be quoted.

    These are statements about what the physics has to do, not statistical tests of a
    hypothesis, and they exist because the first full run failed them (DEVIATIONS 36).

    WHAT WAS REMOVED, AND WHY IT IS NOT A WEAKENING. The original gate required
    corr(log du_sat, log ratio) < 0, on the argument that a wider observer separation cannot
    buy less. That argument treats du_sat as a per-event measure of the separation. It is not:
    du_sat = piE * D_perp/AU, and D_perp is one observatory's orbit, identical for every event.
    Measured, du_sat/piE spans 0.00878 to 0.01002 -- a factor of 1.14 -- while piE spans a
    factor of 35.6 and corr(log piE, log du_sat) = +0.9985. du_sat is piE rescaled by a
    constant, so the check tested whether the gain correlates with piE, which is a different
    claim with competing effects on both sides. It could never have passed or failed for the
    right reason, and no sample size fixes it: at the measured correlation the n needed to put
    it two sigma from zero is ~2.5 million events. DEVIATIONS 37.

    nepR_pk was tried as a replacement monotone axis -- it IS independent of piE -- and shows no
    trend either, for a physical reason: the lowest tercile, median 43 Roman epochs near the
    peak, already shows the full gain. A simultaneous baseline is a geometric constraint, so it
    saturates as soon as a few epochs see the source from both positions at once. There is no
    slope to test because the effect has none.

    Check 2 below is the one that caught the bug, and is unchanged, so this gate would still
    refuse the corrupted run.
    """
    import numpy as np
    ok = True
    print("\n" + "=" * 74)
    print("VALIDATION -- must pass before any ratio here is quotable")
    print("=" * 74)

    # 1. The control must be EXACT, not merely close. These events have no Roman epochs near
    #    the peak, so moving Roman cannot change their light curves at all and the two Fisher
    #    matrices must be the same matrix.
    exact = float((b.ratio == 1.0).mean()) if len(b) else 0.0
    medb = float(b.ratio.median()) if len(b) else float("nan")
    print(f"  1. control (no Roman epochs at peak): n={len(b):,}  median {medb:.6f}  "
          f"{exact:.1%} bit-exactly 1")
    print("     moving Roman cannot touch these, so the ratio must be exactly 1")
    if not len(b) or abs(medb - 1.0) > 1e-6 or exact < 0.5:
        print("     *** FAILED: the control is not exact ***")
        ok = False
    else:
        print("     passed")

    # 2. THE CHECK THAT CAUGHT THE BUG. Unchanged.
    worse = float((c.ratio_tE > 1.01).mean())
    med = float(c.ratio_tE.median())
    print(f"  2. sigma_tE worse for {worse:.1%} of events, median ratio {med:.3f}")
    print("     moving an observer cannot destroy the timescale; expect ~1")
    if worse > 0.5 or med > 1.5:
        print("     *** FAILED: sigma_tE degrades systematically ***")
        ok = False
    else:
        print("     passed")

    # 3. theta_E comes from the deflection amplitude, not from a parallax baseline, so the
    #    astrometric Einstein radius must not move when the observer does.
    m = (c.okB_sat == 1) & (c.okB_nosat == 1) & (c.sigtetE_sat > 0) & (c.sigtetE_nosat > 0)
    if int(m.sum()):
        rt = float((c.loc[m, "sigtetE_sat"] / c.loc[m, "sigtetE_nosat"]).median())
        print(f"  3. sigma(theta_E) ratio median {rt:.6f}  (n={int(m.sum()):,})")
        print("     theta_E is set by the deflection, not the baseline; expect ~1")
        if abs(rt - 1.0) > 0.01:
            print("     *** FAILED: the observer move changed theta_E ***")
            ok = False
        else:
            print("     passed")

    # 4. If the two geometries were not equally conditioned, a difference in sigma could be an
    #    inversion artefact rather than an information statement.
    mc = (c.condA_sat > 0) & (c.condA_nosat > 0)
    if int(mc.sum()):
        rc = float((c.loc[mc, "condA_nosat"] / c.loc[mc, "condA_sat"]).median())
        print(f"  4. condition-number ratio no-sat/sat: median {rc:.4f}  (n={int(mc.sum()):,})")
        print("     the two geometries must be comparably conditioned")
        if not (0.5 < rc < 2.0):
            print("     *** FAILED: the two geometries are not equally conditioned ***")
            ok = False
        else:
            print("     passed")

    # 5. The gain must be a DECREASE, tested against the control rather than against a trend.
    #    The control pins the null at exactly 1, so a sign test is the right instrument.
    v = c.ratio.to_numpy()
    nz = v[v != 1.0]
    if nz.size:
        frac = float((nz < 1.0).mean())
        print(f"  5. {int((nz < 1.0).sum()):,} of {nz.size:,} non-tied events improve "
              f"({frac:.1%})")
        print("     a simultaneous baseline adds parallax information; it must not subtract it")
        if frac <= 0.5:
            print("     *** FAILED: the satellite does not improve sigma(piE) ***")
            ok = False
        else:
            print("     passed")

    if not ok:
        print("\n  VALIDATION FAILED. The ratios in this run are dominated by something other")
        print("  than satellite parallax and MUST NOT be quoted as an H3 result.")
    return ok


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
    # Weighted, like TEMPORAL_GAIN beside it: a raw bar next to a weighted one is not a comparison.
    sat = [R.weighted_median(covered.loc[covered.tb == l, "ratio"], covered.loc[covered.tb == l, "W"])
           if (covered.tb == l).sum() >= 5 else np.nan for l in TE_LABELS]
    # Plotted as the IMPROVEMENT, 100 (1 - ratio) %, on a linear axis from zero. A log axis of
    # the ratio with labels at 1.06x the bar only worked while the gains were large; once both
    # effects are percent-level (post-extinction-fix) the labels land far outside the axes.
    gs = [100.0 * (1.0 - v) if np.isfinite(v) else np.nan for v in sat]
    gt = [100.0 * (1.0 - TEMPORAL_GAIN[l]) for l in TE_LABELS]
    ax.bar(x - w / 2, gs, w, color="#3b6ea5",
           label="satellite baseline  (Roman at L2 vs at Earth), this work")
    ax.bar(x + w / 2, gt, w, color="#c0703a",
           label="temporal baseline  (Rubin filling Roman's gaps), Step F2")
    ax.axhline(0.0, color="#444", lw=1.0)
    ax.set_xticks(x, TE_LABELS)
    top = np.nanmax(gs + gt)
    ax.set_ylim(0, 1.25 * top if top > 0 else 1)
    ax.set_xlabel(r"$t_E$ bin")
    ax.set_ylabel(r"median improvement in $\sigma_{\pi_E}$  [%]")
    ax.set_title("H3c  Two different effects, at the same scale", fontsize=11)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    ax.grid(alpha=0.25, axis="y", lw=0.6)
    for xi, v in zip(x - w / 2, gs):
        if np.isfinite(v):
            ax.text(xi, v + 0.02 * top, f"{v:.1f}", ha="center", fontsize=8)
    for xi, v in zip(x + w / 2, gt):
        ax.text(xi, v + 0.02 * top, f"{v:.1f}", ha="center", fontsize=8)
    fig.text(0.01, -0.06,
             "DIFFERENT PHYSICAL EFFECTS. Pairing the bars compares their size, not their "
             "merit: one is two observers separated in\nspace, the other is one observer's "
             "gaps filled in time. The temporal numbers are in-gap medians from Step F2 "
             "(Deviations 24.2).", fontsize=7.5, color="#444")
    fig.savefig(f"{prefix}c_honest.png", dpi=160, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"wrote {prefix}c_honest.png")


def verdict(covered, blind, passed=True):
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
    if not passed:
        print("  NO RESULT. The validation above failed, so none of these numbers measure")
        print("  satellite parallax. What IS established, and does not depend on the Fisher")
        print("  comparison at all, is the size of the observable itself:")
        print(f"    du_sat median {covered.du_sat.median():.5f}, "
              f"95th pct {covered.du_sat.quantile(0.95):.5f}, "
              f"max {covered.du_sat.max():.5f}  [theta_E]")
        print("  The two observers are separated by a few thousandths of an Einstein radius,")
        print("  so the light-curve perturbation is tiny and a large precision gain would be")
        print("  surprising on physical grounds. See DEVIATIONS.md 35.")
        return
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
