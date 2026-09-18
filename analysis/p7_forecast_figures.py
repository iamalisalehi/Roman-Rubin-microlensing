#!/usr/bin/env python3
"""Step P7: the Fisher-forecast figures a microlensing forecast paper is expected to carry.

Four figures, each answering a question the P5/P6 set does not:

  p7_precision     -- HOW WELL each parameter is measured: the weighted cumulative
                      distribution of sigma(X)/X per survey, for tE, piE, theta_E and the lens
                      mass. This is the figure that says "N% of events are measured to better
                      than 10%", which is the currency of every forecast paper.
  p7_precision_te  -- the same precision as a function of tE, which is where the two surveys
                      differ: Rubin's decade baseline and Roman's dense cadence fail at
                      opposite ends of the timescale axis.
  p7_mass_distance -- WHERE the measurable lenses are, in the mass-distance plane. A forecast
                      that quotes a mass precision without saying which lenses it applies to is
                      quoting a number about a sample, not about the Galaxy.
  p7_sky           -- the sky dependence: event yield and characterised fraction against
                      Galactic coordinates, the map form used by comparable surveys.

WEIGHTING. Every pooled quantity is event-rate weighted (DEVIATIONS.md 41) and quoted with the
Kish N_eff. Per-event ratios carry no weight; only the aggregation does.

A NOTE ON WHAT IS **NOT** HERE. There is no detection-efficiency-versus-tE curve, though it is
the most standard figure of all. The honest denominator for an efficiency is every star drawn,
and this script drops the barren sightlines' rows (77% of a production table) to fit in memory.
Those rows are all non-detections, so including them changes an efficiency and excluding them
inflates it. The simulator computes the efficiency properly over every draw and writes it to
EfLMC<tag>.dat -- that file, not this script, is where an efficiency curve must come from.

USAGE
    .roman/bin/python analysis/p7_forecast_figures.py \
        --run bh=runs/prod_bh_20260917 --run ns=runs/prod_ns_20260917 -o figures/prod/p7
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R          # noqa: E402
import plotstyle as ps        # noqa: E402

COLS = ["tE", "piE", "tetE", "Ml", "Dl", "Ds", "u0", "Vt", "lon", "lat", "w_area",
        "detL", "detR", "detJ", "okA_J", "okA_L", "okA_R", "okB_J", "okB_L", "okB_R",
        "sigtE_J", "sigtE_L", "sigtE_R", "sigpiE_J", "sigpiE_L", "sigpiE_R",
        "sigtetE_J", "sigtetE_L", "sigtetE_R", "relMl_J", "relMl_L", "relMl_R",
        "ndw_L", "ndw_R"]

# (key, display, how to get the fractional error for survey S)
PARAMS = [("tE",   r"$t_{\rm E}$",        lambda d, s: d[f"sigtE_{s}"] / d["tE"]),
          ("piE",  r"$\pi_{\rm E}$",      lambda d, s: d[f"sigpiE_{s}"] / d["piE"]),
          ("tetE", r"$\theta_{\rm E}$",   lambda d, s: d[f"sigtetE_{s}"] / d["tetE"]),
          ("Ml",   r"$M_{\rm L}$",        lambda d, s: d[f"relMl_{s}"])]

SURVEYS = [("J", "joint"), ("L", "rubin"), ("R", "roman")]


class Run:
    def __init__(self, name, directory, chunksize):
        self.name = name
        self.dir = directory
        prov = os.path.join(directory, "files/MONTLMC/files/run_provenance.txt")
        self.prov = R.load_provenance(prov) if os.path.exists(prov) else {}
        self.population = R.assert_same_population([self.prov], what=f"run '{name}'")
        tag = self.prov.get("population_tag")
        if tag is None:
            cands = [f for f in os.listdir(directory)
                     if f.startswith("test") and f.endswith(".dat")]
            if len(cands) != 1:
                sys.exit(f"{directory}: cannot tell which table to read ({cands})")
            tag = cands[0][len("test"):-len(".dat")]
        self.tag = tag
        self.table = os.path.join(directory, f"test{tag}.dat")
        self.map = os.path.join(directory, f"files/MONTLMC/files/MapLMC{tag}.dat")
        self.logs = [os.path.join(directory, f) for f in ("run.log", "run2.log")
                     if os.path.exists(os.path.join(directory, f))]
        self.df = R.load_events(self.table, usecols=COLS, chunksize=chunksize,
                                keep=R.keep_weightable(self.map, self.logs))
        w, self.wlabel = R.attach_weight(self.df, self.map, self.logs)
        self.df["W"] = w
        self.neff = R.kish_neff(w)
        print(f"  {name}: {len(self.df):,} draws, "
              f"{int((self.df.detJ == 1).sum()):,} joint detections, N_eff={self.neff:,.0f}")

    @property
    def colour(self):
        return ps.POPULATION.get(self.name, ps.INK)

    @property
    def label(self):
        return ps.POPULATION_LABEL.get(self.name, self.name)

    def stamp(self):
        return f"{self.population}: {self.wlabel} · commit {self.prov.get('git_commit','?')}"


def stamp_for(runs):
    return "   |   ".join(r.stamp() for r in runs)


def ok_mask(d, surv, key):
    """Events where THIS survey measured THIS parameter.

    Two gates, not one. `okA`/`okB` says the matrix inverted; a positive sigma says this
    particular parameter was actually free (an inactive one carries the -1 sentinel, and a
    sentinel must never enter a distribution). theta_E and the mass come from the ASTROMETRIC
    matrix, so they are gated on okB rather than okA -- a row can have okA = 1 and okB = 0.
    """
    astrometric = key in ("tetE", "Ml")
    gate = d[f"okB_{surv}"] if astrometric else d[f"okA_{surv}"]
    return (gate == 1).to_numpy()


def wcdf(ax, x, w, **kw):
    """Weighted CDF, drawn as a step. Returns the fraction below 0.1, the usual headline."""
    good = np.isfinite(x) & (x > 0) & np.isfinite(w) & (w > 0)
    if good.sum() < 20:
        return np.nan
    x, w = x[good], w[good]
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w.sum()
    ax.plot(x[o], 100 * c, **kw)
    return 100 * float(w[x < 0.1].sum() / w.sum())


def fig_precision(runs, out):
    """sigma(X)/X for the four parameters, per survey. One panel per parameter."""
    fig, axes = ps.figure(width="double", height=5.0, nrows=2, ncols=2)
    axes = np.ravel(axes)
    headline = {}

    for ax, (key, disp, frac) in zip(axes, PARAMS):
        for r in runs:
            d = r.df
            for surv, skey in SURVEYS:
                m = ok_mask(d, surv, key)
                if m.sum() < 20:
                    continue
                x = frac(d, surv).to_numpy(float)[m]
                w = d["W"].to_numpy(float)[m]
                # One line style per population, one colour per survey: the figure has two
                # dimensions and giving each its own visual channel is the only way both stay
                # readable when four populations' worth of curves overlap.
                style = "-" if r.name == runs[0].name else "--"
                lab = f"{r.name} {ps.SURVEY_LABEL[skey].split()[0]}" if len(runs) > 1 \
                    else ps.SURVEY_LABEL[skey]
                f10 = wcdf(ax, x, w, color=ps.SURVEY[skey], ls=style, lw=1.3, label=lab)
                headline[(r.name, key, surv)] = f10
        ax.axvline(0.1, color=ps.MUTED, lw=0.7, ls=":")
        ax.set_xscale("log")
        ax.set_xlim(1e-3, 10)
        ax.set_xlabel(rf"$\sigma$({disp})$/${disp}")
        ax.set_ylabel("events below this [%]")
        ps.panel_label(ax, f"({'abcd'[list(PARAMS).index((key, disp, frac))]}) {disp}")

    # Figure-level legend above the grid. Six entries will not fit inside a panel without
    # sitting on either the curves or the panel label, and it belongs to all four panels
    # anyway.
    h, l = axes[0].get_legend_handles_labels()
    if h:
        fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 1.07), ncol=3,
                   frameon=False, fontsize=7, labelcolor=ps.INK, handlelength=1.8)

    # EACH CURVE HAS ITS OWN DENOMINATOR, and saying so is not pedantry: a reader comparing
    # "Roman 81%" against "joint 44%" for theta_E would conclude that adding Rubin makes the
    # forecast worse, which is impossible -- the joint information matrix is the sum of the
    # parts. What differs is the sample. Roman characterises only inside its footprint, where
    # every event is well covered; the joint curve also contains every Rubin-only event, whose
    # theta_E is poor. For the controlled, same-event comparison see the p6 synergy figure.
    counts = []
    for r in runs:
        for surv, skey in SURVEYS:
            n = int(ok_mask(r.df, surv, "tetE").sum())
            counts.append(f"{r.name}/{ps.SURVEY_LABEL[skey].split()[0]} {n:,}")
    ps.stamp(fig, stamp_for(runs) + "   |   dotted: the 10% criterion   |   each curve is over "
                  "ITS OWN characterised set, not a common one (theta_E counts: "
             + ", ".join(counts) + ")")
    return ps.save_figure(fig, f"{out}_precision"), headline


def fig_precision_te(runs, out):
    """Median sigma(X)/X against tE -- where each survey's strength lies."""
    fig, axes = ps.figure(width="double", height=5.0, nrows=2, ncols=2)
    axes = np.ravel(axes)

    for ax, (key, disp, frac) in zip(axes, PARAMS):
        for r in runs:
            d = r.df
            for surv, skey in SURVEYS:
                m = ok_mask(d, surv, key)
                if m.sum() < 50:
                    continue
                x = frac(d, surv).to_numpy(float)[m]
                te = d["tE"].to_numpy(float)[m]
                w = d["W"].to_numpy(float)[m]
                g = np.isfinite(x) & (x > 0) & (te > 0)
                if g.sum() < 50:
                    continue
                bins = np.geomspace(max(te[g].min(), 1.0), te[g].max(), 9)
                cen, med = [], []
                for i in range(len(bins) - 1):
                    b = g & (te >= bins[i]) & (te < bins[i + 1])
                    if b.sum() < 20 or w[b].sum() <= 0:
                        continue
                    cen.append(np.sqrt(bins[i] * bins[i + 1]))
                    med.append(R.weighted_median(x[b], w[b]))
                if len(cen) >= 3:
                    ax.plot(cen, med, "-" if r.name == runs[0].name else "--",
                            color=ps.SURVEY[skey], lw=1.3,
                            label=f"{r.name} {ps.SURVEY_LABEL[skey].split()[0]}"
                                  if len(runs) > 1 else ps.SURVEY_LABEL[skey])
        ax.axhline(0.1, color=ps.MUTED, lw=0.7, ls=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$t_{\rm E}$ [d]")
        ax.set_ylabel(rf"median $\sigma$({disp})$/${disp}")
        ps.panel_label(ax, f"({'abcd'[list(PARAMS).index((key, disp, frac))]}) {disp}")
    ps.legend(axes[0], loc="lower left")
    ps.stamp(fig, stamp_for(runs) + "   |   dotted line: the 10% criterion")
    return ps.save_figure(fig, f"{out}_precision_te")


def fig_mass_distance(runs, out):
    """Where the measurable lenses sit in the mass-distance plane.

    A mass precision quoted without this is a statement about a sample, not about the Galaxy:
    the lenses whose mass Roman can measure are not drawn uniformly from the population.
    """
    fig, axes = ps.figure(width="double", height=3.0, ncols=max(len(runs), 2))
    axes = np.ravel(axes)

    for ax, r in zip(axes, runs):
        d = r.df
        det = ((d.detL == 1) | (d.detR == 1) | (d.detJ == 1)).to_numpy()
        meas = det & (d["relMl_J"] > 0).to_numpy() & (d["relMl_J"] < 0.1).to_numpy()
        ml, dl = d["Ml"].to_numpy(float), d["Dl"].to_numpy(float)
        w = d["W"].to_numpy(float)

        # A MAP OF THE FRACTION, not two scatter clouds on top of each other. With ~10^5
        # detections the scatter version is a solid blob whose structure is entirely hidden by
        # overplotting, and the eye reads point DENSITY, which here is the Monte Carlo's
        # sampling rather than anything physical. The fraction measured per cell is the
        # quantity the panel is actually about, and it is weighted.
        xb = np.linspace(0.0, np.nanpercentile(dl[det], 99.5), 22)
        yb = np.geomspace(max(ml[det].min(), 1e-3), ml[det].max(), 22)
        num, _, _ = np.histogram2d(dl[meas], ml[meas], bins=[xb, yb], weights=w[meas])
        den, _, _ = np.histogram2d(dl[det], ml[det], bins=[xb, yb], weights=w[det])
        cnt, _, _ = np.histogram2d(dl[det], ml[det], bins=[xb, yb])
        with np.errstate(invalid="ignore", divide="ignore"):
            frac = np.where((den > 0) & (cnt >= 20), 100 * num / den, np.nan)

        im = ax.pcolormesh(xb, yb, frac.T, cmap="viridis", shading="auto", vmin=0)
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label(r"$\sigma(M_{\rm L})/M_{\rm L}<0.1$ [%]", fontsize=7)
        cb.ax.tick_params(labelsize=6)
        if meas.sum():
            print(f"  {r.name}: mass measured to 10% for "
                  f"{100 * w[meas].sum() / max(w[det].sum(), 1e-30):.3f}% of detections "
                  f"(weighted), median Dl {R.weighted_median(dl[meas], w[meas]):.2f} kpc")
        ax.set_yscale("log")
        ax.grid(False)
        ax.set_xlabel(r"lens distance $D_{\rm L}$ [kpc]")
        ax.set_ylabel(r"lens mass $M_{\rm L}$ [$M_\odot$]")
        ps.plain_log_ticks(ax, yb[0], yb[-1], "y")
        ps.panel_label(ax, f"({'ab'[list(runs).index(r)]}) {r.name}")
    for ax in axes[len(runs):]:
        ax.set_visible(False)
    ps.stamp(fig, stamp_for(runs)
             + "   |   blank cells: fewer than 20 detections, too few to quote a fraction")
    return ps.save_figure(fig, f"{out}_mass_distance")


def fig_sky(runs, out):
    """Event yield and characterised fraction across the sky.

    The yield is the SUM OF WEIGHTS per cell, which is the rate-weighted event count and the
    only version of "how many events here" that means anything; the raw count per cell is a
    statement about where the Monte Carlo spent its draws.
    """
    fig, axes = ps.figure(width="double", height=3.0, ncols=2)
    ax1, ax2 = np.ravel(axes)
    r = runs[0]
    d = r.df
    det = ((d.detL == 1) | (d.detR == 1) | (d.detJ == 1)).to_numpy()
    lon, lat = d["lon"].to_numpy(float), d["lat"].to_numpy(float)
    w = d["W"].to_numpy(float)
    # "Characterised" -- okA_J == 1 -- is true for very nearly every detection, so mapping it
    # gives a uniformly saturated panel that says nothing. The fraction measured TO 10% has
    # real dynamic range and is the quantity a forecast cares about.
    good = det & (d["okA_J"] == 1).to_numpy() & (d["sigtE_J"] > 0).to_numpy()
    good &= (d["sigtE_J"].to_numpy(float) / d["tE"].to_numpy(float) < 0.1)

    nb = 26
    lb = np.linspace(lon.min(), lon.max(), nb + 1)
    bb = np.linspace(lat.min(), lat.max(), nb + 1)
    yield_, _, _ = np.histogram2d(lon[det], lat[det], bins=[lb, bb], weights=w[det])
    nchar, _, _ = np.histogram2d(lon[good], lat[good], bins=[lb, bb], weights=w[good])
    cnt, _, _ = np.histogram2d(lon[det], lat[det], bins=[lb, bb])

    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where((yield_ > 0) & (cnt >= 20), nchar / yield_, np.nan)
    # Empty cells must be blank, not zero: a cell the scan never reached and a cell with no
    # events are different statements, and a zero would colour the second like the first.
    shown = np.where(yield_ > 0, yield_ / np.nanmax(yield_), np.nan)

    for ax, z, lab, tag in ((ax1, shown.T, "relative weighted event yield", "(a)"),
                            (ax2, 100 * frac.T,
                             r"joint $\sigma(t_{\rm E})/t_{\rm E}<0.1$ [%]", "(b)")):
        im = ax.pcolormesh(lb, bb, z, cmap="viridis", shading="auto")
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label(lab, fontsize=7)
        cb.ax.tick_params(labelsize=6)
        ax.set_xlabel(r"$\ell$ [deg]")
        ax.set_ylabel(r"$b$ [deg]")
        ax.grid(False)
        ps.panel_label(ax, tag)
    ps.stamp(fig, r.stamp() + "   |   blank cells: no detected events in the cell")
    return ps.save_figure(fig, f"{out}_sky")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, metavar="NAME=DIR")
    ap.add_argument("-o", "--out", default="figures/p7")
    ap.add_argument("--chunksize", type=int, default=500_000)
    a = ap.parse_args()

    ps.use_paper_style()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)

    runs = []
    print("loading:")
    for spec in a.run:
        name, directory = spec.split("=", 1)
        runs.append(Run(name, directory, a.chunksize))

    written = []
    w1, headline = fig_precision(runs, a.out)
    written += w1
    written += fig_precision_te(runs, a.out)
    written += fig_mass_distance(runs, a.out)
    written += fig_sky(runs, a.out)

    print("\n  fraction of measured events better than 10%:")
    for (pop, key, surv), v in sorted(headline.items()):
        if np.isfinite(v):
            print(f"    {pop:4s} {key:5s} {surv}: {v:6.2f}%")

    for p in written:
        print(f"wrote {p}")


if __name__ == "__main__":
    main()
