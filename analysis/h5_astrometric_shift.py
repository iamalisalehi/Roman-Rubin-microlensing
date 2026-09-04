#!/usr/bin/env python3
"""Step H5 -- the astrometric shift itself, not just the precision on theta_E.

WHY THIS EXISTS
---------------
The astrometric microlensing signal is computed by the simulator (`s.def1c`/`s.def2c` in
lightcurve(), stored per epoch in `l.soux`/`l.souy`) and it already feeds the astrometric
Fisher matrix. But until now nothing READ it. F4 plots sigma_tetE -- the forecast *precision*
on the angular Einstein radius -- which answers "can we fit theta_E?" and never answers "how
big is the wobble, and can anyone see it?". Those are different questions and only the second
one decides whether astrometric microlensing is a measurement or an extrapolation.

THE PHYSICS
-----------
A point lens shifts the centroid of the source's light by

    delta_theta(u) = theta_E * u / (u^2 + 2)          [mas]

where u is the lens-source separation in Einstein radii. This is NOT monotonic: it rises from
zero, peaks at u = sqrt(2), and falls again. So

    delta_theta_max = theta_E / sqrt(8) ~ 0.354 theta_E     if u0 <= sqrt(2)
                    = theta_E * u0 / (u0^2 + 2)             if u0 >  sqrt(2)

Two consequences the photometric intuition gets wrong:

1. **The astrometric peak is not the photometric peak.** Photometry peaks at u = u0 (closest
   approach, t = t0). Astrometry peaks at u = sqrt(2), which for u0 < sqrt(2) happens at

       |t - t_0| = tE * sqrt(2 - u0^2)

   i.e. TWICE per event, symmetrically, up to ~1.41 tE either side of t0. For a 300-day event
   that is well over a year from the photometric peak. Panel (c) asks what that does when the
   observatory has 72-day seasons: an event whose photometric peak Roman catches may have both
   its astrometric peaks in a gap, and vice versa.

2. **A high-magnification event is a poor astrometric event.** Small u0 gives a huge
   photometric peak, but the centroid shift at u0 -> 0 goes to ZERO. The astrometric signal
   comes from the wings.

BLENDING, AND WHY TWO CURVES IN PANEL (a)
-----------------------------------------
The simulator's astrometric positions (`s.pos1c`/`pos2c`) add the deflection undiluted:
the modelled centroid is the SOURCE's centroid. A real measurement sees the centroid of all
the light in the aperture, so an unlensed blend of fraction (1 - fb) drags the measured shift
down by roughly fb. Panel (a) therefore shows both: the shift as the simulator's Fisher matrix
sees it, and the same shift multiplied by the F146 source-flux fraction, which is closer to
what a real centroid measurement would deliver. The gap between them is the size of a
simplification that is currently in the code -- recorded in OPEN_ITEMS.md, not fixed here.

SCOPE
-----
Default is events Roman actually observed (ndw_R > 0). Roman is the only astrometric
instrument here worth the name, and outside its footprint the Roman-alone curves are vacuous.

CAVEAT ON ANY TABLE WRITTEN BEFORE STEP H4
------------------------------------------
Before H4, Roman's per-epoch astrometric error in the simulator was Rubin's model -- in fact a
STALE Rubin value from a different timestep (DEVIATIONS.md 29.2). Panels (d), (e) and (f) read
sigma_tetE and therefore inherit that. The script prints a warning when the provenance says
the run predates H4. Panels (a), (b) and (c) are unaffected: they are computed from theta_E,
u0 and tE, which are drawn quantities, not forecasts.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R

INK = "#1a1a1a"
MUTED = "#6b6b6b"
SURFACE = "#ffffff"
GRID = "#f0f0f0"

COLOR = {"joint": "#0d366b", "roman": "#c2410c", "rubin": "#0e7490"}
LABEL = {"joint": "joint fit", "roman": "Roman alone", "rubin": "Rubin alone"}
ACCENT = "#7c3aed"

# Roman WFI per-exposure astrometric precision, F146, in mas.
# KEEP IN SYNC WITH Bulge.h (ROMAN_AST_* constants, Step H4). Sources: Sanderson et al. 2019
# (arXiv:1712.05420) and arXiv:2608.24998. Duplicated here rather than parsed out of the
# header because parsing a C++ expression is more fragile than a cross-reference; if the C++
# constants move, move these.
ROMAN_AST_FLOOR = 1.1        # mas, 1% of the 110 mas pixel
ROMAN_AST_MFLR = 20.62
ROMAN_AST_MBKG = 23.5
ROMAN_AST_SBKG = 10.0
ROMAN_AST_SLOPE_SRC = 0.33285
ROMAN_AST_SLOPE_BKG = 0.4

U_AST_PEAK = np.sqrt(2.0)    # the separation at which the centroid shift is maximal


def roman_ast_error(mag):
    """errRomanA() from Bulge.h/helper.cpp, vectorized. Per EXPOSURE, in mas."""
    mag = np.asarray(mag, dtype=float)
    out = np.full(mag.shape, ROMAN_AST_FLOOR)
    mid = (mag > ROMAN_AST_MFLR) & (mag <= ROMAN_AST_MBKG)
    out[mid] = ROMAN_AST_FLOOR * 10.0 ** (ROMAN_AST_SLOPE_SRC * (mag[mid] - ROMAN_AST_MFLR))
    hi = mag > ROMAN_AST_MBKG
    out[hi] = ROMAN_AST_SBKG * 10.0 ** (ROMAN_AST_SLOPE_BKG * (mag[hi] - ROMAN_AST_MBKG))
    return out


def centroid_shift(theta_e, u):
    """delta_theta = theta_E * u / (u^2 + 2)  [mas]."""
    return theta_e * u / (u * u + 2.0)


def max_shift(df):
    """Maximum centroid shift reached by each event, in mas.

    theta_E/sqrt(8) when the trajectory passes through u = sqrt(2); otherwise the shift at
    closest approach, because delta_theta FALLS beyond u = sqrt(2).
    """
    te = df["tetE"].to_numpy(dtype=float)
    u0 = df["u0"].to_numpy(dtype=float)
    u_at_max = np.where(u0 <= U_AST_PEAK, U_AST_PEAK, u0)
    return centroid_shift(te, u_at_max)


def ast_peak_offset(df):
    """|t - t0| at which the centroid shift peaks, in days. Zero when u0 > sqrt(2)."""
    tE = df["tE"].to_numpy(dtype=float)
    u0 = df["u0"].to_numpy(dtype=float)
    inside = u0 <= U_AST_PEAK
    off = np.zeros_like(tE)
    off[inside] = tE[inside] * np.sqrt(2.0 - u0[inside] ** 2)
    return off


def style(ax):
    ax.tick_params(colors=MUTED, labelsize=8)
    for sp in ax.spines.values():
        sp.set_color("#d0d0d0")
    ax.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.set_facecolor(SURFACE)


# ---------------------------------------------------------------------------- panels
def panel_shift_cdf(ax, df, rows):
    """(a) How big is the wobble, and is it above Roman's per-epoch precision?"""
    dmax = max_shift(df)
    fb = df["blend_F146"].to_numpy(dtype=float)
    dmax_bl = dmax * fb

    prec = roman_ast_error(df["magb_F146"].to_numpy(dtype=float))
    n = len(df)

    fb_med = float(np.median(fb))
    for v, c, lab in ((dmax, ACCENT, "source centroid (what the Fisher matrix uses)"),
                      (dmax_bl, "#c026d3",
                       rf"$\times f_{{\rm b}}$ (blend-diluted; median $f_{{\rm b}}={fb_med:.2f}$)")):
        v = np.sort(v[np.isfinite(v) & (v > 0)])
        ax.step(v, np.arange(1, len(v) + 1) / n, where="post", color=c, lw=1.9, zorder=3,
                label=lab)

    ax.axvline(ROMAN_AST_FLOOR, color="#b0b0b0", lw=0.9, ls=":", zorder=1)
    ax.text(ROMAN_AST_FLOOR, 0.02, "  1.1 mas\n  1-exposure floor", color=MUTED, fontsize=7,
            ha="left", va="bottom", zorder=1)

    above = float(np.mean(dmax > prec))
    above_bl = float(np.mean(dmax_bl > prec))
    ax.set_xscale("log")
    ax.set_xlabel(r"maximum centroid shift  $\delta\theta_{\rm max}$  [mas]", color=INK, fontsize=9)
    ax.set_ylabel("cumulative fraction of sample", color=INK, fontsize=9)
    ax.set_title("(a)  The astrometric signal\n"
                 f"{above*100:.1f}% exceed Roman's own per-exposure precision "
                 f"({above_bl*100:.1f}% blend-diluted)",
                 color=INK, fontsize=10, loc="left", pad=7, linespacing=1.5)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right", labelcolor=INK)
    ax.set_ylim(0, 1)
    rows.append(dict(panel="a", metric="frac_shift_above_per_exposure_precision", value=above))
    rows.append(dict(panel="a", metric="frac_blenddiluted_above_precision", value=above_bl))
    rows.append(dict(panel="a", metric="median_max_shift_mas", value=float(np.median(dmax))))
    rows.append(dict(panel="a", metric="median_blend_F146", value=fb_med))
    rows.append(dict(panel="a", metric="median_roman_per_exposure_precision_mas",
                     value=float(np.median(prec))))


def panel_u0(ax, df, rows):
    """(b) How the wobble compares with what Roman can measure -- per exposure and stacked.

    The previous version of this panel plotted delta_theta_max/theta_E against u0, which is a
    deterministic function of u0 and therefore drew the analytic curve twice. It said nothing
    the formula did not. This one is the panel that reconciles (a) with (d): the signal sits
    far BELOW Roman's single-exposure precision and far ABOVE its stacked precision, which is
    why theta_E is forecastable from data in which no individual exposure sees the wobble.
    """
    dmax = max_shift(df)
    tE = df["tE"].to_numpy(dtype=float)
    teE = df["tetE"].to_numpy(dtype=float)
    nR = df["ndw_R"].to_numpy(dtype=float)
    prec = roman_ast_error(df["magb_F146"].to_numpy(dtype=float))

    sc = ax.scatter(tE, dmax, s=5, c=np.log10(teE), cmap="viridis", alpha=0.55, lw=0, zorder=3)
    cb = ax.figure.colorbar(sc, ax=ax, pad=0.015)
    cb.set_label(r"$\log_{10}(\theta_{\rm E}/{\rm mas})$", color=INK, fontsize=8)
    cb.ax.tick_params(colors=MUTED, labelsize=7)

    per_exp = float(np.median(prec))
    stacked = float(np.median(prec / np.sqrt(np.maximum(nR, 1.0))))
    ax.axhline(per_exp, color="#8a8a8a", lw=1.1, ls=":", zorder=2,
               label=f"single exposure ({per_exp:.1f} mas)")
    ax.axhline(stacked, color=COLOR["roman"], lw=1.1, ls="--", zorder=2,
               label=rf"stacked over $N_{{\rm epoch}}$ ({stacked*1e3:.0f} $\mu$as)")
    ax.legend(frameon=False, fontsize=7.5, loc="lower right", labelcolor=INK)

    reach = float(np.mean(df["u0"].to_numpy(dtype=float) <= U_AST_PEAK))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t_{\rm E}$  [d]", color=INK, fontsize=9)
    ax.set_ylabel(r"$\delta\theta_{\rm max}$  [mas]", color=INK, fontsize=9)
    ax.set_title("(b)  Below one exposure, far above the stack\n"
                 f"{reach*100:.1f}% of events pass through $u=\\sqrt{{2}}$ "
                 r"and reach $\theta_{\rm E}/\sqrt{8}$",
                 color=INK, fontsize=10, loc="left", pad=7, linespacing=1.5)
    rows.append(dict(panel="b", metric="frac_reaching_theta_over_sqrt8", value=reach))
    rows.append(dict(panel="b", metric="median_stacked_precision_mas", value=stacked))
    rows.append(dict(panel="b", metric="median_shift_over_stacked_precision",
                     value=float(np.median(dmax)) / stacked if stacked > 0 else np.nan))


def panel_season(ax, df, rows):
    """(c) The astrometric peak is displaced from t0 -- across a season edge, sometimes."""
    off = ast_peak_offset(df)
    dt_edge = df["dt_edge"].to_numpy(dtype=float)

    # dt_edge is signed days from t0 to the NEAREST Roman season boundary, negative when t0
    # fell inside a season. So |dt_edge| is the distance to that boundary, and the astrometric
    # peak crosses it when the offset exceeds it. Approximate -- it only knows about the
    # nearest edge, not the whole season structure -- and stated as such on the figure.
    crosses = off > np.abs(dt_edge)
    in_season = dt_edge < 0.0

    ax.scatter(np.abs(dt_edge[in_season]), off[in_season], s=5, alpha=0.30, lw=0,
               c=COLOR["roman"], zorder=3, label=r"$t_0$ in a Roman season")
    ax.scatter(np.abs(dt_edge[~in_season]), off[~in_season], s=5, alpha=0.30, lw=0,
               c=COLOR["rubin"], zorder=3, label=r"$t_0$ in a gap / off-mission")
    lim = np.array([1.0, 3.0e3])
    ax.plot(lim, lim, color=INK, lw=1.1, ls="--", zorder=4, label="astrometric peak on the edge")

    f_cross = float(np.mean(crosses))
    f_cross_in = float(np.mean(crosses[in_season])) if in_season.any() else np.nan
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lim)
    ax.set_ylim(1.0, 3.0e3)
    ax.set_xlabel(r"$|$days from $t_0$ to nearest Roman season edge$|$", color=INK, fontsize=9)
    ax.set_ylabel(r"$|t-t_0|$ of the astrometric peak  [d]", color=INK, fontsize=9)
    ax.set_title("(c)  Astrometry peaks at $u=\\sqrt{2}$, not at $t_0$\n"
                 f"{f_cross*100:.1f}% peak beyond a season edge (above the line)",
                 color=INK, fontsize=10, loc="left", pad=7, linespacing=1.5)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right", labelcolor=INK)
    rows.append(dict(panel="c", metric="frac_ast_peak_crosses_season_edge", value=f_cross))
    rows.append(dict(panel="c", metric="frac_crossing_given_t0_in_season", value=f_cross_in))


def panel_tetE_cdf(ax, df, surveys, rows):
    """(d) The forecast on theta_E, per survey -- F4's panel (c), kept for context."""
    n = len(df)
    ax.axvline(0.1, color="#b0b0b0", lw=0.9, ls=":", zorder=1)
    ax.text(0.1, 0.02, " 10%", color=MUTED, fontsize=7, ha="left", va="bottom")
    for s in surveys:
        v = (R.sigma(df, "tetE", s) / df["tetE"]).dropna().to_numpy()
        v = np.sort(v[v > 0])
        below = float(np.sum(v < 0.1)) / n
        ax.step(v, np.arange(1, len(v) + 1) / n, where="post", color=COLOR[s], lw=1.9,
                zorder=3, label=f"{LABEL[s]}   {below*100:4.1f}% < 10%")
        rows.append(dict(panel="d", metric=f"frac_sigtetE_below_10pct_{s}", value=below))
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 1e3)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\sigma_{\theta_{\rm E}}/\theta_{\rm E}$", color=INK, fontsize=9)
    ax.set_ylabel("cumulative fraction of sample", color=INK, fontsize=9)
    ax.set_title("(d)  Forecast precision on $\\theta_{\\rm E}$\n"
                 "astrometric matrix, gated on okB",
                 color=INK, fontsize=10, loc="left", pad=7, linespacing=1.5)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right", labelcolor=INK)


def panel_invariant(ax, df, rows):
    """(e) sigma_joint <= sigma_single on theta_E. A physics invariant, so a bug detector."""
    for s in ("roman", "rubin"):
        x = R.sigma(df, "tetE", s) / df["tetE"]
        y = R.sigma(df, "tetE", "joint") / df["tetE"]
        m = x.notna() & y.notna()
        ax.scatter(x[m], y[m], s=5, alpha=0.30, lw=0, c=COLOR[s], zorder=3,
                   label=f"joint vs {LABEL[s]}")
        r = (y[m] / x[m]).to_numpy()
        rows.append(dict(panel="e", metric=f"median_ratio_tetE_joint_over_{s}",
                         value=float(np.median(r)) if len(r) else np.nan))
        rows.append(dict(panel="e", metric=f"n_above_unity_{s}",
                         value=float(np.sum(r > 1.0 + 1e-9))))
    lim = np.array([1e-4, 1e3])
    ax.plot(lim, lim, color=INK, lw=1.1, ls="--", zorder=4)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_xlabel(r"$\sigma_{\theta_{\rm E}}/\theta_{\rm E}$, single survey", color=INK, fontsize=9)
    ax.set_ylabel(r"joint fit", color=INK, fontsize=9)
    ax.set_title("(e)  Adding data cannot hurt\n"
                 "every point must be on or below the 1:1 line",
                 color=INK, fontsize=10, loc="left", pad=7, linespacing=1.5)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left", labelcolor=INK)


def panel_mass_plane(ax, df, rows):
    """(f) The mass solution needs BOTH matrices: Ml = theta_E / (kappa piE)."""
    piE = df["piE"].to_numpy(dtype=float)
    teE = df["tetE"].to_numpy(dtype=float)
    ok = np.isfinite(piE) & np.isfinite(teE) & (piE > 0) & (teE > 0)
    ml = teE[ok] / (8.144 * piE[ok])

    sc = ax.scatter(piE[ok], teE[ok], s=5, c=np.log10(ml), cmap="viridis", alpha=0.55,
                    lw=0, zorder=3)
    cb = ax.figure.colorbar(sc, ax=ax, pad=0.015)
    cb.set_label(r"$\log_{10}(M_{\rm L}/M_\odot)$", color=INK, fontsize=8)
    cb.ax.tick_params(colors=MUTED, labelsize=7)
    for m in (0.1, 1.0, 10.0):
        g = np.logspace(-3, 1, 100)
        ax.plot(g, 8.144 * m * g, color="#b0b0b0", lw=0.8, ls=":", zorder=2)
        ax.text(g[-1], 8.144 * m * g[-1], f" {m:g}$M_\\odot$", color=MUTED, fontsize=6.5,
                ha="left", va="center")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\pi_{\rm E}$  (photometric matrix)", color=INK, fontsize=9)
    ax.set_ylabel(r"$\theta_{\rm E}$  [mas]  (astrometric matrix)", color=INK, fontsize=9)
    ax.set_title("(f)  Why both matrices are needed\n"
                 r"$M_{\rm L}=\theta_{\rm E}/(\kappa\pi_{\rm E})$, "
                 r"$\kappa=8.144$ mas/$M_\odot$",
                 color=INK, fontsize=10, loc="left", pad=7, linespacing=1.5)
    rows.append(dict(panel="f", metric="median_Ml_msun", value=float(np.median(ml))))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("events", help="per-event table (test5.dat), or a cached subset of it")
    ap.add_argument("-o", "--out", default="h5_astrometric_shift.png")
    ap.add_argument("--provenance", default=None)
    ap.add_argument("--chunksize", type=int, default=500_000)
    ap.add_argument("--scope", choices=("footprint", "all"), default="footprint")
    args = ap.parse_args()

    # Filtered chunked read: joint detections are 1.3% of the table and reading the rest
    # first gets the process OOM-killed, silently, on this machine (OPEN_ITEMS.md).
    df = R.load_events(args.events, keep=lambda c: c["detJ"] == 1, chunksize=args.chunksize)
    stamp = R.describe(args.events, args.provenance)
    print(stamp)

    prov_path = R.find_provenance(args.provenance, near=args.events)
    prov = R.load_provenance(prov_path) if prov_path else {}
    pre_h4 = "satellite_parallax" not in prov
    if pre_h4:
        print("  WARNING: this table predates Step H4. Roman's per-epoch astrometric error\n"
              "  was a stale Rubin value (DEVIATIONS.md 29.2), so panels (d), (e) and (f)\n"
              "  inherit it. Panels (a), (b), (c) are computed from drawn quantities and\n"
              "  are unaffected.")

    if args.scope == "footprint":
        df = df[df["ndw_R"] > 0]
        scope_note = "events Roman actually observed (ndw_R > 0)"
    else:
        scope_note = "every joint-detected event"
    df = df[(df["tetE"] > 0) & (df["u0"] > 0) & (df["tE"] > 0)]
    print(f"  scope={args.scope}: {len(df)} events")
    if len(df) == 0:
        print("  nothing in scope -- no figure written")
        return

    fig, axes = plt.subplots(2, 3, figsize=(16.2, 9.4), facecolor=SURFACE)
    rows = []
    panel_shift_cdf(axes.flat[0], df, rows)
    panel_u0(axes.flat[1], df, rows)
    panel_season(axes.flat[2], df, rows)
    panel_tetE_cdf(axes.flat[3], df, ("joint", "roman", "rubin"), rows)
    panel_invariant(axes.flat[4], df, rows)
    panel_mass_plane(axes.flat[5], df, rows)
    for ax in axes.flat:
        style(ax)

    fig.suptitle("The astrometric microlensing signal, and what it is worth",
                 color=INK, fontsize=13.5, x=0.035, ha="left", y=0.975)
    fig.text(0.035, 0.945, f"{len(df):,} joint-detected events -- {scope_note}",
             color=MUTED, fontsize=9, va="top")

    fig.subplots_adjust(left=0.05, right=0.985, bottom=0.115, top=0.885,
                        wspace=0.30, hspace=0.44)
    footer = (
        r"Centroid shift $\delta\theta=\theta_{\rm E}u/(u^2+2)$, maximal at $u=\sqrt{2}$ "
        r"where $\delta\theta=\theta_{\rm E}/\sqrt{8}$ -- so the astrometric peak is "
        r"displaced from $t_0$ by $t_{\rm E}\sqrt{2-u_0^2}$ and a high-magnification event "
        "is a poor astrometric one.\n"
        "Panel (c)'s season test uses the NEAREST season edge only, so it is a lower bound "
        "on how often the astrometric peak lands in different observing conditions than the "
        "photometric one.\n"
        + ("PANELS (d)-(f) REST ON A PRE-H4 TABLE: Roman's astrometric error was a stale "
           "Rubin value (DEVIATIONS.md 29.2). Re-run after the next production run.\n"
           if pre_h4 else "")
        + stamp)
    fig.text(0.035, 0.062, footer, color=MUTED, fontsize=7.5, va="top", linespacing=1.7)

    fig.savefig(args.out, dpi=190, facecolor=SURFACE)
    print(f"  wrote {args.out}")

    csv = os.path.splitext(args.out)[0] + ".csv"
    import pandas as pd
    pd.DataFrame(rows).to_csv(csv, index=False)
    print(f"  wrote {csv}")
    for r in rows:
        print(f"    {r['panel']}  {r['metric']:<46s} {r['value']:.6g}")


if __name__ == "__main__":
    main()
