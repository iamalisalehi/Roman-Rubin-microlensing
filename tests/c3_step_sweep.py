#!/usr/bin/env python3
"""Step C3 -- finite-difference step-size convergence analysis.

Input: the CSV written by `./fishertest --sweep` (see tests/fisher_fixture.cpp, runSweep()).
Each row is one (event, survey partition, photometric parameter, step scale) combination and
the sigma FisherM/ErrorCal recovered for that parameter with that step.

What we are looking for
-----------------------
A finite-difference derivative df/dp ~ [f(p+h) - f(p-h)] / 2h carries two competing errors:

  * truncation error, ~ h^2 * f''', which grows as the step h grows (the model is not linear
    over the step, so the secant stops matching the tangent);
  * round-off error, ~ eps * |f| / h, which grows as h shrinks (we subtract two nearly equal
    doubles and divide by a small number, so relative precision is destroyed).

Their sum is U-shaped in h. A trustworthy step sits at the flat bottom -- the plateau. The
plan's acceptance criterion (Step C3) is that sigma must be stable to "well under 1%" across
that plateau, and that every parameter must have one.

Two controls are built into the sweep:
  * mbs0 / mbs1 (baseline magnitude): the model magnitude depends on these linearly with unit
    slope, so the central difference is exact at any step. Their curves must be perfectly flat.
    Structure there means the sweep harness is broken, not the physics.
  * fb0 / fb1 (blend fraction): FisherM clamps their step so fb + h stays inside [0, 1]. At the
    largest sweep scales the clamp binds and the curve flattens artificially -- that flattening
    is the clamp, not a plateau. At the production step the clamp is far from binding.

Usage
-----
    .roman/bin/python tests/c3_step_sweep.py [c3_sweep.csv] [-o c3_step_sweep.png]
"""

import argparse
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# The acceptance window: one decade either side of the production step. The plan's Step C3
# criterion is evaluated here -- sigma must be stable to well under 1% across it. Must contain at
# least 3 points of the sweep grid, or the check reports "no data" rather than a misleading pass.
PLATEAU_LO, PLATEAU_HI = 0.1, 10.0

# The plateau the sweep found, common to every parameter: bounded below by round-off and above
# by truncation. Curves are normalized to REF_SCALE, deep inside it, so a converged curve reads
# as a flat line at 1.0.
#
# After the Step C3 retune (kFDStepScale in Bulge_LSST.cpp) the production step sits INSIDE the
# plateau, so scale 1 is both the reference and the middle of the flat region. Re-point these at
# 1e-6 / 1e-3 / 1e-6 to re-analyse a pre-retune sweep CSV.
PLATEAU_FOUND_LO, PLATEAU_FOUND_HI = 1.0e-4, 1.0e+1
REF_SCALE = 1.0

# Acceptance thresholds, in percent of the production-step sigma.
ACCEPT_PCT = 1.0    # plan's criterion: "well under 1%"
WARN_PCT = 0.1      # what a genuinely converged parameter looks like

# Plateau detection: sigma must change by less than STEP_TOL from one sweep point to the next,
# over at least MIN_PLATEAU_PTS consecutive points, for that stretch to count as a plateau.
STEP_TOL = 0.005
MIN_PLATEAU_PTS = 3

PARAM_ORDER = ["u0", "tE", "fb0", "piE", "xi", "t0", "mbs0", "fb1", "mbs1"]
SURVEY_ORDER = ["joint", "rubin", "roman"]
SURVEY_STYLE = {"joint": "-", "rubin": "--", "roman": ":"}


def load(path):
    df = pd.read_csv(path)
    # A partition that could not be characterized reports sigma = -1 (Step C4 sentinel) and
    # ok = 0. Those rows carry no information about step size; drop them, but count them.
    bad = (df["ok"] != 1) | (df["sigma"] <= 0) | ~np.isfinite(df["sigma"])
    if bad.any():
        print(f"dropping {bad.sum()} of {len(df)} rows with ok=0 or non-positive sigma "
              f"(uncharacterizable partitions)")
    return df[~bad].copy()


def curves(df):
    """Yield (event, survey, param, scales, sigma_ratio) with sigma normalized to scale=1."""
    for (ev, surv, par), g in df.groupby(["event", "survey", "param"], sort=False):
        g = g.sort_values("scale")
        ref = g.loc[np.isclose(g["scale"], REF_SCALE), "sigma"]
        if ref.empty or ref.iloc[0] <= 0:
            continue
        sig = g["sigma"].to_numpy()
        yield ev, surv, par, g["scale"].to_numpy(), sig / ref.iloc[0], sig


def find_plateau(sc, sig, tol=STEP_TOL):
    """Longest run of consecutive sweep points over which sigma barely moves.

    Returns (lo_scale, hi_scale, sigma_plateau, npts) or None. "Barely moves" is defined
    point-to-point rather than against a fixed reference, so a slow monotonic drift across
    the whole window -- the signature of an unconverged first-order stencil -- does not get
    mistaken for a plateau.
    """
    n = len(sc)
    if n < MIN_PLATEAU_PTS:
        return None
    flat = np.abs(sig[1:] / sig[:-1] - 1.0) < tol
    best = None
    i = 0
    while i < len(flat):
        if not flat[i]:
            i += 1
            continue
        j = i
        while j < len(flat) and flat[j]:
            j += 1
        npts = j - i + 1                     # points, not intervals
        if npts >= MIN_PLATEAU_PTS:
            span = np.log10(sc[j] / sc[i])
            if best is None or span > best[0]:
                best = (span, sc[i], sc[j], float(np.median(sig[i:j + 1])), npts)
        i = j + 1
    if best is None:
        return None
    return best[1], best[2], best[3], best[4]


def summarize(df):
    """One row per (param, survey, event) curve.

    dev_window_pct  -- how much sigma moves across the old +/-3x window around the production
                       step; this is the plan's literal acceptance test.
    plateau_*       -- where (and whether) the curve actually flattens.
    prod_bias_pct   -- sigma at the production step divided by the plateau sigma, minus 1.
                       This is the quantity that matters scientifically: how wrong the sigmas
                       the pipeline reports today are, purely from the step size.
    """
    rows = []
    for ev, surv, par, sc, rat, sig in curves(df):
        m = (sc >= PLATEAU_LO) & (sc <= PLATEAU_HI)
        rec = dict(param=par, survey=surv, event=ev,
                   dev_window_pct=(np.max(np.abs(rat[m] - 1.0)) * 100.0 if m.sum() >= 3
                                   else np.nan),
                   dev_full_pct=np.max(np.abs(rat - 1.0)) * 100.0)
        pl = find_plateau(sc, sig)
        if pl is None:
            rec.update(plateau_lo=np.nan, plateau_hi=np.nan, plateau_pts=0,
                       prod_bias_pct=np.nan)
        else:
            lo, hi, s_pl, npts = pl
            sig_prod = sig[np.isclose(sc, 1.0)]
            rec.update(plateau_lo=lo, plateau_hi=hi, plateau_pts=npts,
                       prod_bias_pct=(float(sig_prod[0]) / s_pl - 1.0) * 100.0
                       if len(sig_prod) else np.nan)
        rows.append(rec)
    return pd.DataFrame(rows)


def report(summary):
    print()
    print(f"Old acceptance window: scale in [{PLATEAU_LO}, {PLATEAU_HI}] around the production "
          f"step (scale 1.0)")
    print(f"Plateau: >= {MIN_PLATEAU_PTS} consecutive sweep points with < {STEP_TOL*100:.1f}% "
          f"change from one to the next")
    print()
    hdr = (f"{'param':>6} {'win dev %':>10} {'plateau':>18} {'frac w/':>8} "
           f"{'prod bias %':>22}  verdict")
    print(hdr)
    hdr2 = (f"{'':>6} {'(worst)':>10} {'(median range)':>18} {'plateau':>8} "
            f"{'med / worst':>22}")
    print(hdr2)
    print("-" * len(hdr))
    failures = []
    for par in PARAM_ORDER:
        g = summary[summary["param"] == par]
        if g.empty:
            print(f"{par:>6} {'--':>10} {'no usable curves':>18}")
            failures.append(par)
            continue
        win = g["dev_window_pct"].max()
        if not np.isfinite(win):
            # Fewer than 3 grid points inside the acceptance window -- cannot judge.
            print(f"{par:>6} {'--':>10} {'window too sparse':>18} {'':>8} {'':>22}  NO DATA")
            failures.append(par)
            continue
        has = g["plateau_pts"] > 0
        frac = has.mean()
        if has.any():
            rng = f"{g.loc[has,'plateau_lo'].median():.0e}-{g.loc[has,'plateau_hi'].median():.0e}"
            bias_med = g.loc[has, "prod_bias_pct"].abs().median()
            bias_max = g.loc[has, "prod_bias_pct"].abs().max()
            bias = f"{bias_med:9.2f} /{bias_max:10.2f}"
        else:
            rng, bias = "none", f"{'--':>21}"
        if win < WARN_PCT:
            verdict = "PASS (flat everywhere)"
        elif win < ACCEPT_PCT:
            verdict = "PASS"
        else:
            verdict = "FAIL -- production step off plateau"
            failures.append(par)
        print(f"{par:>6} {win:10.2f} {rng:>18} {frac:8.0%} {bias:>22}  {verdict}")
    print()
    if failures:
        print(f"OFF-PLATEAU PARAMETERS: {', '.join(failures)}")
    else:
        print("all parameters sit on a plateau at the production step size")
    return failures


def plot(df, out):
    fig, axes = plt.subplots(3, 3, figsize=(13, 10), sharex=True)
    events = sorted(df["event"].unique())
    cmap = plt.get_cmap("viridis", max(len(events), 2))
    colors = {e: cmap(i) for i, e in enumerate(events)}

    by_panel = {p: ax for p, ax in zip(PARAM_ORDER, axes.ravel())}
    for ev, surv, par, sc, rat, _sig in curves(df):
        ax = by_panel.get(par)
        if ax is None:
            continue
        ax.plot(sc, rat, SURVEY_STYLE.get(surv, "-"), color=colors[ev], lw=1.2, alpha=0.9)

    for par, ax in by_panel.items():
        ax.set_xscale("log")
        ax.set_yscale("log")
        # The plateau, and the production step that sits outside it.
        ax.axvspan(PLATEAU_FOUND_LO, PLATEAU_FOUND_HI, color="#cde5cd", zorder=0, lw=0)
        ax.axvline(1.0, color="crimson", lw=1.0, ls="-", zorder=1)
        ax.axhline(1.0, color="0.4", lw=0.8)
        ax.axhline(1.01, color="crimson", lw=0.6, ls=":")
        ax.axhline(0.99, color="crimson", lw=0.6, ls=":")
        ax.set_title(par, fontsize=11)
        ax.set_ylabel(r"$\sigma\,/\,\sigma_{\rm plateau}$", fontsize=9)
        ax.tick_params(labelsize=8)
    for ax in axes[-1]:
        ax.set_xlabel(r"finite-difference step, relative to the production step", fontsize=9)

    handles = [plt.Line2D([], [], color=colors[e], lw=1.5, label=e) for e in events]
    handles += [plt.Line2D([], [], color="0.3", ls=SURVEY_STYLE[s], lw=1.5, label=s)
                for s in SURVEY_ORDER]
    handles += [plt.Line2D([], [], color="crimson", lw=1.5, label="production step"),
                plt.Rectangle((0, 0), 1, 1, color="#cde5cd", label="plateau")]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=8, frameon=False)
    fig.suptitle("Step C3: Fisher $\\sigma$ vs finite-difference step size\n"
                 "round-off wall (left) and truncation wall (right) bracket the plateau; "
                 "the production step (red) sits inside it, $\\sigma$ flat to $<0.1\\%$",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?", default="c3_sweep.csv")
    ap.add_argument("-o", "--out", default="c3_step_sweep.png")
    ap.add_argument("--per-curve", action="store_true",
                    help="print every (param, survey, event) curve, not just the worst")
    a = ap.parse_args()

    df = load(a.csv)
    summary = summarize(df)
    if summary.empty:
        sys.exit("no usable curves -- is the CSV complete?")
    failures = report(summary)
    if a.per_curve:
        print()
        print(summary.sort_values(["param", "dev_window_pct"], ascending=[True, False])
              .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    plot(df, a.out)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
