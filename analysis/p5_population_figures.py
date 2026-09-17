#!/usr/bin/env python3
"""Step P5 -- publication-ready figures comparing lens populations.

    .roman/bin/python analysis/p5_population_figures.py \
        --run bulge=/path/to/2026-xx_bulge_run \
        --run bh=/path/to/2026-xx_bh_run \
        --run ns=/path/to/2026-xx_ns_run \
        -o figures/p5

Each --run names a population and a RUN DIRECTORY holding that run's own outputs:
`test<tag>.dat`, `files/MONTLMC/files/MapLMC<tag>.dat`, `files/MONTLMC/files/run_provenance.txt`
and (optionally) `run.log`, which supplies `nsim` for any sightline whose map row a killed run
lost. Each population is weighted with its OWN map file: the weight is per sightline, and two
runs never share sightline statistics.

WHAT IS AND IS NOT POOLED. Populations are drawn one curve each and never summed into a single
statistic -- `romanlib.assert_same_population()` enforces that per table. The mass functions
differ, and the event-rate weight carries a sqrt(Ml) factor that is only valid for the mass
function actually sampled (Deviation 45), so a pooled BH+bulge number would be meaningless. The
BH:NS ratio, if a paper wants one, is applied to the finished per-population numbers.

THE FIGURES
  (1) mass_function  -- what was drawn against what was assumed. A validation plot, and the
      first thing to check on a new run: if the histogram does not lie on the analytic curve,
      nothing downstream means anything.
  (2) efficiency     -- detection efficiency against lens mass, per survey. The plot that turns
      a yield into a statement about which masses a survey can find.
  (3) mass_precision -- sigma(M_L)/M_L against lens mass, per survey.
  (4) timescale      -- the tE distribution each population produces, weighted.
  (5) gain           -- joint over single-survey sigma(tE) and sigma(piE) against lens mass:
      where the combination earns its keep.

Every panel is event-rate weighted and reports N_eff, because the weight costs precision and a
weighted number without it cannot be judged (Deviation 41).
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotstyle as ps
import romanlib as R

# Enough columns for every panel; reading all 92 over millions of draws is the documented
# OOM trap, and these figures need the UNDETECTED draws too (they are the efficiency
# denominator), so rows cannot be filtered away instead.
COLS = ["tE", "piE", "tetE", "Ml", "Vt", "Ds", "lon", "lat", "w_area",
        "detL", "detR", "detJ", "ndw_R",
        "okA_J", "okA_L", "okA_R", "sigtE_J", "sigtE_L", "sigtE_R",
        "sigpiE_J", "sigpiE_L", "sigpiE_R", "relMl_J", "relMl_L", "relMl_R"]

NBIN = 14


class Run:
    """One population's run: its table, its weights, its provenance."""

    def __init__(self, name, directory, chunksize):
        self.name = name
        self.dir = directory
        prov_path = os.path.join(directory, "files/MONTLMC/files/run_provenance.txt")
        self.prov = R.load_provenance(prov_path) if os.path.exists(prov_path) else {}
        self.population = R.assert_same_population([self.prov], what=f"run '{name}'")
        tag = self.prov.get("population_tag")
        if tag is None:
            # Pre-population runs wrote test5.dat; find whatever single table is there.
            cands = [f for f in os.listdir(directory)
                     if f.startswith("test") and f.endswith(".dat")]
            if len(cands) != 1:
                sys.exit(f"{directory}: cannot tell which table to read ({cands}); "
                         "expected exactly one test<tag>.dat")
            tag = cands[0][len("test"):-len(".dat")]
        self.tag = tag
        self.table = os.path.join(directory, f"test{tag}.dat")
        self.map = os.path.join(directory, f"files/MONTLMC/files/MapLMC{tag}.dat")
        logs = [os.path.join(directory, f) for f in ("run.log", "run2.log")]
        self.logs = [f for f in logs if os.path.exists(f)]

        self.df = R.load_events(self.table, usecols=COLS, chunksize=chunksize)
        w, self.wlabel = R.attach_weight(self.df, self.map, self.logs)
        self.df["W"] = w
        self.mass_range = self._mass_range()
        print(f"  {name}: {len(self.df):,} draws, {int((self.df.detJ == 1).sum()):,} joint "
              f"detections, population={self.population}, {self.wlabel}")

    def _mass_range(self):
        r = self.prov.get("lens_mass_range")
        if r:
            lo, hi = (float(x) for x in r.split()[:2])
            return lo, hi
        m = self.df["Ml"].to_numpy()
        return float(m.min()), float(m.max())

    @property
    def colour(self):
        return ps.POPULATION.get(self.name, ps.INK)

    @property
    def label(self):
        return ps.POPULATION_LABEL.get(self.name, self.name)

    def mass_bins(self):
        lo, hi = self.mass_range
        return np.geomspace(lo, hi, NBIN + 1)


def weighted_binned(x, w, bins, mask=None):
    """Sum of weights per bin, and per bin among `mask`. Returns (centres, total, selected)."""
    x = np.asarray(x, float)
    w = np.asarray(w, float)
    idx = np.digitize(x, bins) - 1
    ok = (idx >= 0) & (idx < len(bins) - 1)
    tot = np.zeros(len(bins) - 1)
    sel = np.zeros(len(bins) - 1)
    np.add.at(tot, idx[ok], w[ok])
    if mask is not None:
        m = np.asarray(mask, bool) & ok
        np.add.at(sel, idx[m], w[m])
    centres = np.sqrt(bins[:-1] * bins[1:])
    return centres, tot, sel


def analytic_mass_function(name, lo, hi, grid):
    """The mass function each population is defined by, normalised to unit area in log M.

    Only the two new populations have a closed form worth drawing: the bulge population's
    present-day mass function is the Kroupa IMF pushed through the initial-final mass
    relation, which is a piecewise map with a spike at the neutron-star mass and no tidy
    expression. For that one the histogram is the statement.
    """
    if name == "bh":
        # Flat in log M: constant density per decade.
        y = np.ones_like(grid)
    elif name == "ns":
        mu, sig = 1.35, 0.15
        y = np.exp(-0.5 * ((grid - mu) / sig) ** 2) * grid   # x grid: per log M
    else:
        return None
    inside = (grid >= lo) & (grid <= hi)
    y = np.where(inside, y, 0.0)
    area = np.trapezoid(y / grid, grid) if hasattr(np, "trapezoid") else np.trapz(y / grid, grid)
    return y / area if area > 0 else y


def fig_mass_function(runs, out):
    """Drawn mass function against the assumed one. The first check on any new run."""
    fig, ax = ps.figure(width="double", height=3.0)
    seen = []
    for r in runs:
        m = r.df["Ml"].to_numpy()
        w = r.df["W"].to_numpy()
        nb = 40 if len(m) > 2000 else max(8, int(np.sqrt(len(m))))
        bins = np.geomspace(max(m.min(), 1e-3), m.max(), nb + 1)
        widths = np.diff(np.log(bins))

        # THE SOLID CURVE IS UNWEIGHTED, AND THAT IS THE POINT OF THIS PANEL. It asks one
        # question -- did the sampler draw the mass function it was asked to? -- so it must
        # compare like with like: the raw draws against the assumed mass function.
        #
        # The event-rate weight would break that comparison rather than improve it. W carries
        # sqrt(Ml) (Deviation 41), so even a perfectly flat-in-log black-hole sample comes out
        # tilted upward by half a decade per decade once weighted, and a reader comparing it
        # against the flat dashed line would see a disagreement that is not there. The
        # weighted curve is still worth showing -- it is the population that actually produces
        # events -- so it is drawn thin, and labelled as a different thing.
        hu, _ = np.histogram(m, bins=bins)
        dens = np.where(hu > 0, hu / (hu.sum() * widths), np.nan)
        hw, _ = np.histogram(m, bins=bins, weights=w)
        dens_w = np.where(hw > 0, hw / (hw.sum() * widths), np.nan)
        #
        # Both densities are masked where their bin is empty: an empty bin has density zero,
        # and zero on a log axis is minus infinity, so plotted literally it draws a spike to
        # the bottom of the frame and drags the y-range down by ten orders of magnitude. An
        # empty bin is missing data, not a measured zero.
        centres = np.sqrt(bins[:-1] * bins[1:])
        ax.plot(centres, dens_w, color=r.colour, lw=0.7, alpha=0.55)
        ax.step(centres, dens, where="mid", color=r.colour, label=r.label)
        # Both curves set the y-range, or the thin weighted one can fall outside the frame.
        seen.append(dens[np.isfinite(dens)])
        seen.append(dens_w[np.isfinite(dens_w)])
        grid = np.geomspace(*r.mass_range, 200)
        a = analytic_mass_function(r.name, *r.mass_range, grid)
        if a is not None:
            ax.plot(grid, np.where(a > 0, a, np.nan),
                    color=r.colour, lw=0.9, ls=(0, (3, 2)), alpha=0.85)
    ax.set_xscale("log")
    ax.set_yscale("log")
    # Limits from the populated bins, with a little headroom for the legend.
    vals = np.concatenate([s for s in seen if s.size]) if any(s.size for s in seen) else None
    if vals is not None and vals.size:
        ax.set_ylim(0.5 * vals.min(), 30.0 * vals.max())
    ax.set_xlabel(r"lens mass  $M_{\rm L}$  [M$_\odot$]")
    ax.set_ylabel(r"drawn density  ${\rm d}N/{\rm d}\ln M_{\rm L}$")
    ps.panel_label(ax, "(a)")
    ps.legend(ax, loc="upper center", ncol=len(runs))
    ps.stamp(fig, "Solid: the raw draws, against the dashed assumed mass function -- the "
                  "sampler check.  Thin: the same population event-rate weighted, which "
                  "tilts by sqrt(M) and is not meant to match the dashed line.  "
                  "Empty bins are masked, not plotted as zero.")
    return ps.save_figure(fig, f"{out}_mass_function")


def fig_efficiency(runs, out):
    from matplotlib.ticker import LogLocator, NullFormatter

    fig, axes = ps.figure(width="double", height=2.7, ncols=len(runs), sharey=True)
    axes = np.atleast_1d(axes)
    finite = []
    for ax, r in zip(axes, runs):
        bins = r.mass_bins()
        for s, key in (("joint", "detJ"), ("rubin", "detL"), ("roman", "detR")):
            c, tot, sel = weighted_binned(r.df["Ml"], r.df["W"], bins, r.df[key] == 1)
            # A bin where this survey detected nothing is masked, not drawn as zero: zero on
            # a log axis is minus infinity and plots as a spike to the bottom of the frame.
            # "No detections in 30 draws" is an upper limit, not a measurement of zero, and a
            # figure should not assert more than the sample supports.
            eff = np.where((tot > 0) & (sel > 0), 100.0 * sel / np.maximum(tot, 1e-30), np.nan)
            finite.append(eff[np.isfinite(eff)])
            ax.plot(c, eff, color=ps.SURVEY[s], label=ps.SURVEY_LABEL[s], marker="o",
                    ms=2.2, lw=1.2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        # Tick labelling on a log axis needs two different answers here, because the
        # populations span wildly different ranges.
        #
        # Over decades (bulge, bh) the default labels every log MINOR tick, which at this
        # panel width overlaps into mush: decade majors only, minor labels off.
        #
        # Under a decade (neutron stars live in 1.1-2.2 Msun) that same rule leaves the axis
        # with NO tick labels at all, since the range contains no power of ten -- an axis with
        # no numbers on it. There, explicit ticks with a plain-number formatter.
        lo, hi = r.mass_range
        if np.log10(hi / lo) >= 1.0:
            ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=4))
            ax.xaxis.set_minor_formatter(NullFormatter())
        else:
            ticks = np.round(np.geomspace(lo, hi, 4), 2)
            ax.set_xticks(ticks)
            ax.set_xticklabels([f"{t:g}" for t in ticks])
            ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel(r"$M_{\rm L}$  [M$_\odot$]")
        ps.panel_label(ax, r.name, loc="lower left")
    vals = np.concatenate([f for f in finite if f.size]) if any(f.size for f in finite) else None
    if vals is not None and vals.size:
        axes[0].set_ylim(0.5 * vals.min(), 3.0 * vals.max())
    axes[0].set_ylabel("detection efficiency  [%]")
    ps.legend(axes[0], loc="upper left")
    ps.stamp(fig, "Efficiency = weighted detections / weighted draws per mass bin; bins with "
                  "no detections are masked.  "
                  + "  ".join(f"{r.name}: N_eff {R.kish_neff(r.df['W']):,.0f}" for r in runs))
    return ps.save_figure(fig, f"{out}_efficiency_vs_mass")


def fig_mass_precision(runs, out):
    fig, ax = ps.figure(width="single")
    for r in runs:
        bins = r.mass_bins()
        d = r.df[r.df["detJ"] == 1]
        v = d["relMl_J"].to_numpy()
        ok = np.isfinite(v) & (v > 0)
        idx = np.digitize(d["Ml"].to_numpy()[ok], bins) - 1
        w = d["W"].to_numpy()[ok]
        vv = v[ok]
        centres, med = [], []
        for b in range(len(bins) - 1):
            m = idx == b
            if m.sum() >= 5 and R.kish_neff(w[m]) >= 3:
                centres.append(np.sqrt(bins[b] * bins[b + 1]))
                med.append(R.weighted_median(vv[m], w[m]))
        if centres:
            ax.plot(centres, med, color=r.colour, label=r.label, marker="o", ms=2.5, lw=1.3)
    ax.axhline(0.1, color=ps.MUTED, lw=0.8, ls=":")
    ax.text(ax.get_xlim()[0], 0.1, " 10%", color=ps.MUTED, fontsize=6, va="bottom")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"lens mass  $M_{\rm L}$  [M$_\odot$]")
    ax.set_ylabel(r"median  $\sigma(M_{\rm L})/M_{\rm L}$,  joint fit")
    ps.panel_label(ax, "(b)")
    ps.legend(ax, loc="best")
    ps.stamp(fig, "Bins with fewer than 5 events or N_eff < 3 are not drawn.")
    return ps.save_figure(fig, f"{out}_mass_precision")


def fig_timescale(runs, out):
    fig, ax = ps.figure(width="single")
    for r in runs:
        d = r.df[r.df["detJ"] == 1]
        t = d["tE"].to_numpy()
        w = d["W"].to_numpy()
        o = np.argsort(t)
        ax.step(t[o], np.cumsum(w[o]) / w.sum(), where="post", color=r.colour,
                label=f"{r.label}\n  median {R.weighted_median(t, w):.0f} d")
    ax.set_xscale("log")
    ax.set_xlabel(r"$t_{\rm E}$  [days]")
    ax.set_ylabel("cumulative fraction of detections")
    ax.set_ylim(0, 1)
    ps.panel_label(ax, "(c)")
    ps.legend(ax, loc="lower right")
    ps.stamp(fig, "Joint detections, event-rate weighted.")
    return ps.save_figure(fig, f"{out}_timescale")


def fig_gain(runs, out):
    fig, axes = ps.figure(width="double", height=2.7, ncols=2, sharex=True)
    for ax, param in zip(axes, ("tE", "piE")):
        for r in runs:
            bins = r.mass_bins()
            d = r.df[r.df["detJ"] == 1]
            ratio = R.ratio_joint_over(d, param, "roman")
            ok = ratio.notna().to_numpy()
            idx = np.digitize(d["Ml"].to_numpy()[ok], bins) - 1
            w = d["W"].to_numpy()[ok]
            vv = ratio.to_numpy()[ok]
            centres, med = [], []
            for b in range(len(bins) - 1):
                m = idx == b
                if m.sum() >= 5 and R.kish_neff(w[m]) >= 3:
                    centres.append(np.sqrt(bins[b] * bins[b + 1]))
                    med.append(R.weighted_median(vv[m], w[m]))
            if centres:
                ax.plot(centres, med, color=r.colour, label=r.label, marker="o", ms=2.5, lw=1.3)
        ax.axhline(1.0, color=ps.MUTED, lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(r"$M_{\rm L}$  [M$_\odot$]")
        ax.set_ylabel(rf"median $\sigma_{{\rm joint}}/\sigma_{{\rm Roman}}$  (${param}$)"
                      .replace("piE", r"\pi_{\rm E}").replace("tE", r"t_{\rm E}"))
    ps.panel_label(axes[0], "(d)")
    ps.panel_label(axes[1], "(e)")
    ps.legend(axes[0], loc="lower left")
    ps.stamp(fig, "Below 1 means the joint fit is tighter than Roman alone. "
                  "Bins with fewer than 5 events or N_eff < 3 are not drawn.")
    return ps.save_figure(fig, f"{out}_gain_vs_mass")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, metavar="NAME=DIR",
                    help="population name and its run directory; repeatable")
    ap.add_argument("-o", "--out", default="figures/p5")
    ap.add_argument("--chunksize", type=int, default=500_000)
    a = ap.parse_args()

    ps.use_paper_style()
    runs = []
    print("loading:")
    for spec in a.run:
        if "=" not in spec:
            sys.exit(f"--run wants NAME=DIR, got '{spec}'")
        name, directory = spec.split("=", 1)
        runs.append(Run(name, directory, a.chunksize))

    written = []
    written += fig_mass_function(runs, a.out)
    written += fig_efficiency(runs, a.out)
    written += fig_mass_precision(runs, a.out)
    written += fig_timescale(runs, a.out)
    written += fig_gain(runs, a.out)
    for p in written:
        print(f"wrote {p}")


if __name__ == "__main__":
    main()
