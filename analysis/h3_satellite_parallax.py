#!/usr/bin/env python3
"""What Roman's displacement to L2 buys for the parallax forecast.

Every detected event is characterised twice inside one run (`--pair-satellite`): once with
Roman at L2 and once with `L2_OFFSET_AU` zeroed, with the same draw, epochs and photometry; only
the observer moves. Each row of `h3_pair.dat` is a genuine pair and the ratios are per-event.

Two separate runs cannot be compared instead: the RNG is a single mt19937_64 stream and the
per-event path draws conditionally on `acceptRubin or acceptRoman`, so moving the observer
changes what is detectable and the streams fork. An unpaired comparison gave a spurious
sigma_tE ratio of 1.15, which measured the difference between the two detected populations.

Gating: a ratio needs both forecasts, `okA_sat == 1 and okA_nosat == 1`, with -1.0 (not
measured) excluded.

Control: events with `nepR_pk == 0` have no Roman epochs near the peak, so the satellite cannot
act on them and their ratio must be 1. They are kept and reported; if they drift from 1 the
measurement is wrong.
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

# Temporal-baseline gain for panel c: median sigma_joint/sigma_Roman for piE on in-gap events
# (t0zone == 1, Roman-covered), event-rate weighted, from the bulge production run. This is a
# different effect (Rubin filling Roman's season gaps in time, not two observers separated in
# space); the panel puts the two on one axis at the right scale. The in-season control for the
# same quantity is 0.999-1.000 in every bin.
TEMPORAL_GAIN = {"10-30 d": 0.984, "30-100 d": 0.999, "100-300 d": 1.000, "> 300 d": 1.000}
TE_EDGES = [10.0, 30.0, 100.0, 300.0, np.inf]
TE_LABELS = list(TEMPORAL_GAIN.keys())
BG = "#fcfcfb"
# h3_pair.dat column layout. The header line begins with a lone '#', so pandas is told
# comment="#" and the names are supplied here. The order must match src/sim/record.cpp, or every
# column silently shifts.
COLS_LEGACY = ("lon lat tE u0 piE tetE du_sat "
               "okA_sat okA_nosat okB_sat okB_nosat "
               "sigtE_sat sigtE_nosat sigpiE_sat sigpiE_nosat "
               "sigpiER_sat sigpiER_nosat sigtetE_sat sigtetE_nosat "
               "sigpiEb_sat sigpiEb_nosat relMl_sat relMl_nosat "
               "condA_sat condA_nosat condB_sat condB_nosat "
               "nepL_pk nepR_pk w_area").split()

# Ml, Dl, Ds and Vt carry the event-rate weight, W ~ sqrt(Ml)*Vt*Z(Ds). They cannot be
# reconstructed from the other columns (pi_rel = 1/Dl - 1/Ds is one equation in two unknowns),
# so a file without them is readable but not weightable.
COLS = COLS_LEGACY + ["Ml", "Dl", "Ds", "Vt"]

# Superseded 15-column layout, recognised only so it can be refused.
COLS_OLD_15 = 15


def load(path):
    probe = pd.read_csv(path, sep=r"\s+", comment="#", header=None, nrows=1,
                        engine="python")
    n = probe.shape[1]
    if n == COLS_OLD_15:
        sys.exit(
            f"{path} has the superseded 15-column layout.\n"
            "That file was written before the derivative-reference fix and\n"
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
                 "If src/sim/record.cpp changed the row, update COLS to match it.")
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
            "statistic of the raw sample (that sample over-represents long, slow,\n"
            "massive lenses roughly tenfold).\n"
            "Either re-run --pair-satellite with a binary built from 2026-09-17 or later, which\n"
            "writes the four columns, or pass --unweighted to say deliberately that you want the\n"
            "raw-sample numbers.")

    # Weights are attached to the full frame before subsetting so every selection carries its own.
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

    # The control comes first; if it is not ~1 the rest is meaningless.
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
    """What moving the observer does to the astrometric forecast and to the lens mass.

    The magnification depends only on |u|, so the satellite baseline reaches piE through a change
    in separation. The centroid deflection theta_E u/(u^2+2) is a vector, so moving the observer
    changes its direction as well, giving the astrometric matrix a route to piE independent of
    the photometric one. The lens mass Ml = theta_E/(kappa piE) needs one observable from each
    matrix.
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
    """Compare the condition numbers of the with- and without-satellite matrices.

    Checks that the satScale = 0 configuration does not leave a near-degeneracy that passes the
    okA condition-number gate while leaving marginalised errors unstable.
    """
    print("\n== CONDITIONING of the two matrices")
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

    These are statements about what the physics has to do, not statistical tests.

    There is deliberately no check that the gain correlates with du_sat. du_sat = piE * D_perp/AU
    and D_perp is one observatory's orbit, identical for every event, so du_sat is just piE
    rescaled by a constant (corr(log piE, log du_sat) = +0.9985). The gain does not trend with
    nepR_pk either: a simultaneous baseline is geometric and saturates once a few epochs see the
    source from both positions.
    """
    import numpy as np
    ok = True
    print("\n" + "=" * 74)
    print("VALIDATION -- must pass before any ratio here is quotable")
    print("=" * 74)

    # 1. The control must be exact: with no Roman epochs near the peak, moving Roman cannot
    #    change the light curve and the two Fisher matrices must be identical.
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

    # 2. Moving an observer cannot degrade the timescale forecast.
    worse = float((c.ratio_tE > 1.01).mean())
    med = float(c.ratio_tE.median())
    print(f"  2. sigma_tE worse for {worse:.1%} of events, median ratio {med:.3f}")
    print("     moving an observer cannot destroy the timescale; expect ~1")
    if worse > 0.5 or med > 1.5:
        print("     *** FAILED: sigma_tE degrades systematically ***")
        ok = False
    else:
        print("     passed")

    # 3. theta_E comes from the deflection amplitude, not a parallax baseline, so it must not
    #    move with the observer.
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

    # 4. Unequal conditioning would make a sigma difference an inversion artefact.
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

    # 5. The gain must be a decrease; the control pins the null at exactly 1, so a sign test applies.
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
    # ---- panel a
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
             "two Fisher\nevaluations per event.", fontsize=7.5, color="#444")
    fig.savefig(f"{prefix}a_where.png", dpi=160, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"\nwrote {prefix}a_where.png")

    # ---- panel b
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

    # ---- panel c
    fig, ax = plt.subplots(figsize=(8.8, 5.0), facecolor=BG)
    ax.set_facecolor(BG)
    x = np.arange(len(TE_LABELS))
    w = 0.36
    # Weighted, like TEMPORAL_GAIN beside it.
    sat = [R.weighted_median(covered.loc[covered.tb == l, "ratio"], covered.loc[covered.tb == l, "W"])
           if (covered.tb == l).sum() >= 5 else np.nan for l in TE_LABELS]
    # Plotted as the improvement, 100 (1 - ratio) %, on a linear axis from zero.
    gs = [100.0 * (1.0 - v) if np.isfinite(v) else np.nan for v in sat]
    gt = [100.0 * (1.0 - TEMPORAL_GAIN[l]) for l in TE_LABELS]
    ax.bar(x - w / 2, gs, w, color="#3b6ea5",
           label="satellite baseline  (Roman at L2 vs at Earth), this work")
    ax.bar(x + w / 2, gt, w, color="#c0703a",
           label="temporal baseline  (Rubin filling Roman's gaps)")
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
             "gaps filled in time. The temporal numbers are in-gap medians of the gap-filling analysis."
             "", fontsize=7.5, color="#444")
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
        print("  surprising on physical grounds.")
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
