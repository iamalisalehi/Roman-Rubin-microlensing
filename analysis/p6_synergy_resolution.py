#!/usr/bin/env python3
"""Step P6 figures: what the two surveys do for each other, and what the lens reveals.

Four figures, in the order the science reads:

  p6_synergy      -- how Rubin and Roman help EACH OTHER. Not one number: the help runs
                     both ways and for different reasons, so the figure shows both.
  p6_astrometry   -- the centroid shift, against the precision that has to measure it.
  p6_parallax     -- satellite parallax: what Roman-at-L2 buys over Roman-at-Earth, both
                     in precision and in breaking the degeneracy that limits the lens mass.
  p6_resolution   -- the probability of resolving the two lensing-induced images, after
                     Sajadian & Makler (arXiv:2608.16448).

EVERY POOLED NUMBER HERE IS EVENT-RATE WEIGHTED (DEVIATIONS.md 41), and the script will not
run without either --map (plus --log) or an explicit --unweighted. A raw fraction over this
table describes the sample, not the sky: the Monte Carlo draws lens mass from the number IMF
and velocities from plain Gaussians, so it is missing the rate's sqrt(Ml)*Vt factor. Every
weighted number is quoted with N_eff, the Kish effective sample size, because the weight costs
precision.

POPULATIONS ARE NEVER POOLED. The weight's sqrt(Ml) is valid only for the mass function that
was sampled, so mixing a black-hole run into a neutron-star run is a wrong number rather than
a style choice. Each population is a separate series in every panel.

USAGE
    .roman/bin/python analysis/p6_synergy_resolution.py \
        --run bh=runs/prod_bh_20260917 --run ns=runs/prod_ns_20260917 -o figures/p6
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R          # noqa: E402
import plotstyle as ps        # noqa: E402

# The paper's criterion: the images count as resolvable when at least three recorded data
# points have both images detectable AND separated by more than the bar. Three, not one,
# because a single qualifying epoch is as likely to be noise as signal.
RESOLVE_MIN_EPOCHS = 3

# Roman's per-exposure astrometric floor [mas] -- Bulge.h ROMAN_AST_FLOOR. Used only as a
# reference line; the per-event precision is a function of magnitude and is far worse than
# this for a typical bulge source (measured median 6.69 mas in the v3 run).
ROMAN_AST_FLOOR = 1.1

U_AST_PEAK = np.sqrt(2.0)   # the impact parameter at which the centroid shift is maximal

COLS = ["tE", "piE", "tetE", "Ml", "u0", "Vt", "Ds", "lon", "lat", "w_area",
        "detL", "detR", "detJ", "synClass",
        "okA_J", "okA_L", "okA_R", "okB_J", "okB_L", "okB_R",
        "sigtE_J", "sigtE_L", "sigtE_R",
        "sigpiE_J", "sigpiE_L", "sigpiE_R",
        "sigtetE_J", "sigtetE_L", "sigtetE_R",
        "relMl_J", "relMl_L", "relMl_R",
        "ndw_L", "ndw_R", "fb1",
        "nres5_L", "nres20_L", "nresPSF_L", "dsep_max_L",
        "nres5_R", "nres20_R", "nresPSF_R", "dsep_max_R"]

PAIR_COLS = ("lon lat tE u0 piE tetE du_sat okA_sat okA_nosat okB_sat okB_nosat "
             "sigtE_sat sigtE_nosat sigpiE_sat sigpiE_nosat sigpiER_sat sigpiER_nosat "
             "sigtetE_sat sigtetE_nosat sigpiEb_sat sigpiEb_nosat relMl_sat relMl_nosat "
             "condA_sat condA_nosat condB_sat condB_nosat nepL_pk nepR_pk w_area "
             "Ml Dl Ds Vt").split()


class Run:
    """One population's run: its table, its paired-satellite file, its weights."""

    def __init__(self, name, directory, chunksize, unweighted=False):
        self.name = name
        self.dir = directory
        prov_path = os.path.join(directory, "files/MONTLMC/files/run_provenance.txt")
        self.prov = R.load_provenance(prov_path) if os.path.exists(prov_path) else {}
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

        # Stream, and drop the barren sightlines' rows as they are read rather than after.
        # usecols alone is not enough on a 5M-row table: 43 columns x 5M x 8 bytes is ~1.7 GB
        # per population, and this script holds two at once. The dropped rows have no nsim,
        # hence weight 0, and carry no detection -- see romanlib.keep_weightable.
        self.df = R.load_events(self.table, usecols=COLS, chunksize=chunksize,
                                keep=None if unweighted else R.keep_weightable(self.map, self.logs))
        self.n_barren = self._drop_unweightable(unweighted)
        w, self.wlabel = R.attach_weight(self.df, None if unweighted else self.map,
                                         self.logs, unweighted=unweighted)
        self.df["W"] = w
        self.neff = R.kish_neff(w)

        # The paired satellite file is optional: only a --pair-satellite run has one.
        self.pair = None
        pair_path = os.path.join(directory, "h3_pair.dat")
        if os.path.exists(pair_path):
            self.pair = self._load_pair(pair_path, unweighted)

        print(f"  {name}: {len(self.df):,} draws, "
              f"{int((self.df.detJ == 1).sum()):,} joint detections, "
              f"N_eff={self.neff:,.0f}, population={self.population}, {self.wlabel}"
              + (f", pair={len(self.pair):,}" if self.pair is not None else ", no pair file")
              + (f", dropped {self.n_barren:,} barren rows" if self.n_barren else ""))

    def _drop_unweightable(self, unweighted):
        """Drop rows from BARREN sightlines, after proving they hold nothing we count.

        A sightline that draws stars but ends with no characterised event takes the barren
        branch in Bulge_LSST.cpp, which `continue`s past BOTH the map-row write and the
        `nsim:` print. Its rows are therefore in the table with no draw count to normalise
        them by, and `event_weight` refuses the whole table because of them -- correctly, since
        a missing nsim is indistinguishable from a map file truncated by a kill (Deviation 42),
        and guessing would silently mis-weight a real sightline.

        These rows are not a loss: barren means nothing survived to be aggregated, and at the
        scan's western edge such a sightline runs to the full --maxdraws cap, so a handful of
        them can be a million rows of nothing. But "should be empty" is not "is empty", so the
        rows are CHECKED before being dropped: if any carries a detection or a characterisation,
        this refuses instead, because then the drop would bias the very statistics the figures
        are made of.
        """
        self.df, n = self._drop_from(self.df, unweighted, what=self.table)
        return n

    def _known_sightlines(self):
        known = set()
        if os.path.exists(self.map):
            sl = R.load_sightlines(self.map)
            known |= {(round(float(a), 3), round(float(b), 3)) for a, b in
                      zip(sl["lon"], sl["lat"])}
        known |= set(R.nsim_from_logs(self.logs).keys())
        return known

    def _drop_from(self, df, unweighted, what):
        """The drop itself, applied to the event table and to the paired-satellite file alike.

        Both carry lon/lat and both are weighted the same way, so both hit the same wall on a
        sightline with no draw count. Doing it in one place keeps the two from drifting apart.
        """
        if unweighted:
            return df, 0
        known = self._known_sightlines()
        if not known:
            return df, 0

        lon = df["lon"].round(3)
        lat = df["lat"].round(3)
        keyed = list(zip(lon.to_numpy(), lat.to_numpy()))
        bad = np.array([k not in known for k in keyed], dtype=bool)
        if not bad.any():
            return df, 0

        # A run still in flight has one sightline part-written: its rows are in the table but
        # its nsim line and map row come only when it finishes. That one is legitimately
        # droppable -- it is incomplete, and an incomplete sightline must not be weighted as
        # if it stood for its full sky area. Anything ELSE carrying detections is not.
        inflight = self._last_entered()
        lost = df.loc[bad]
        if inflight is not None:
            not_inflight = bad & ~((lon == inflight[0]) & (lat == inflight[1])).to_numpy()
            lost = df.loc[not_inflight]

        # Only the event table has detection columns; the paired file's rows are all detected
        # by construction, so for it the presence of ANY non-in-flight row is the alarm.
        if {"detL", "detR", "detJ"} <= set(df.columns):
            carried = int(((lost.detL == 1) | (lost.detR == 1) | (lost.detJ == 1)
                           | (lost.okA_J == 1) | (lost.okB_J == 1)).sum())
        else:
            carried = len(lost)
        if carried:
            sys.exit(
                f"{what}: {len(lost):,} rows belong to sightlines with no nsim that are "
                f"not the one in flight, and {carried:,} of them carry a detection or a "
                "characterisation. These are NOT barren leftovers, so dropping them would bias "
                "every pooled number. Recover their draw counts from the run log first.")
        return df.loc[~bad].reset_index(drop=True), int(bad.sum())

    def _last_entered(self):
        """(lon, lat) of the last sightline the scan entered, or None.

        Read from the log rather than the table so it is the SCAN's last sightline, not merely
        the last one that happened to write a row.
        """
        import re
        last = None
        for p in self.logs:
            with open(p, errors="replace") as fh:
                for line in fh:
                    m = re.match(r"^longtitude:\s*(\S+)\s+latitude:\s*(\S+)", line)
                    if m:
                        last = (round(float(m.group(1)), 3), round(float(m.group(2)), 3))
        return last

    def _load_pair(self, path, unweighted):
        """The --pair-satellite side file, weighted the same way as the main table.

        Refuses the LEGACY 30-column layout unless explicitly unweighted: it has no Ml/Vt/Ds,
        so the event-rate weight cannot be formed and a weighted number off it would be a
        silently wrong one (DEVIATIONS.md 44).
        """
        import pandas as pd
        with open(path, errors="replace") as fh:
            header = fh.readline()
        names = header.lstrip("#").split()
        df = pd.read_csv(path, sep=r"\s+", comment="#", names=names, header=None,
                         engine="c", on_bad_lines="skip")
        if "Ml" not in df.columns:
            if not unweighted:
                sys.exit(f"{path}: legacy layout without Ml/Vt/Ds -- it cannot be weighted. "
                         "Re-run with a binary from 2026-09-17 or later, or pass --unweighted.")
            df["W"] = 1.0
            return df
        df, _ = self._drop_from(df, unweighted, what=path)
        w, _ = R.attach_weight(df, None if unweighted else self.map, self.logs,
                               unweighted=unweighted)
        df["W"] = w
        return df

    @property
    def colour(self):
        return ps.POPULATION.get(self.name, ps.INK)

    @property
    def label(self):
        return ps.POPULATION_LABEL.get(self.name, self.name)

    def stamp(self):
        # wlabel already carries N_eff; repeating it was printing the same number twice.
        commit = self.prov.get("git_commit", "?")
        return f"{self.population}: {self.wlabel} · commit {commit}"


def stamp_for(runs):
    return "   |   ".join(r.stamp() for r in runs)


def wfrac(mask, w):
    """Weighted fraction, returning nan rather than 0 when nothing qualifies to be counted."""
    w = np.asarray(w, float)
    tot = w.sum()
    if not np.isfinite(tot) or tot <= 0:
        return np.nan
    return float(w[np.asarray(mask, bool)].sum() / tot)


def mass_bin_count(lo, hi):
    """How many log bins a mass range deserves, from how many decades it spans.

    A fixed count is wrong for both ends of this project: nine bins across the black holes'
    2.5 decades is sensible, but the same nine across the neutron stars' 0.3 decades slices
    them so thinly that each bin is noise, and the curve dives off the left edge as a binning
    artifact that reads as a physical cutoff.
    """
    span = np.log10(hi / lo)
    return int(np.clip(round(4 * span) + 4, 5, 11))


def plain_log_ticks(ax, lo, hi, axis="y"):
    """Plain numbers on a log axis that spans only a decade or two.

    Matplotlib labels log MINOR ticks on short ranges, which at column width collides into
    mush ("6x10^0 4x10^0 3x10^0 ..."); but simply suppressing the minor labels can leave a
    sub-decade axis with no numbers at all. Explicit ticks with a plain formatter is the only
    option that avoids both.

    The range is passed in from the DATA rather than read off the axes: at the point this is
    called matplotlib has not autoscaled yet, so get_ylim() returns provisional limits and the
    span test silently takes the wrong branch. That is why the first attempt at this changed
    nothing.
    """
    import matplotlib.ticker as mt
    a = ax.yaxis if axis == "y" else ax.xaxis
    if not (lo > 0 and hi > lo) or np.log10(hi / lo) > 2.2:
        return
    ticks = np.geomspace(lo, hi, 5)
    a.set_major_locator(mt.FixedLocator(ticks))
    a.set_minor_locator(mt.NullLocator())
    a.set_major_formatter(mt.FuncFormatter(
        lambda v, _: f"{v:.2f}".rstrip("0").rstrip(".") if v < 10 else f"{v:.0f}"))


def max_centroid_shift(df):
    """Largest centroid shift each event reaches, in mas.

    delta(u) = theta_E * u / (u^2 + 2) peaks at u = sqrt(2), NOT at closest approach. So an
    event whose trajectory crosses u = sqrt(2) reaches theta_E/sqrt(8); one that does not
    peaks at its own u0. Taking theta_E/sqrt(8) for every event would overstate the shallow
    ones, and taking the value at u0 would understate the deep ones.
    """
    te = df["tetE"].to_numpy(float)
    u0 = df["u0"].to_numpy(float)
    u_at_max = np.where(u0 <= U_AST_PEAK, U_AST_PEAK, u0)
    return te * u_at_max / (u_at_max * u_at_max + 2.0)


# ---------------------------------------------------------------------------------------
# 1. Synergy: what each survey does for the other
# ---------------------------------------------------------------------------------------
def fig_synergy(runs, out):
    """Two panels, because the help runs in two directions for two different reasons.

    (a) Who detects what. The weighted share of detections that only Rubin saw, only Roman
        saw, or that needed the two together. A survey's value is not only what it finds
        alone -- an event Rubin finds and Roman characterises is a joint result.
    (b) What the joint fit buys, as a function of tE. sigma_joint/sigma_single is bounded
        above by 1 by construction (the joint information matrix is the sum of the parts),
        so the interesting quantity is HOW FAR below 1 it goes and WHERE in tE.
    """
    fig, (ax1, ax2) = ps.figure(width="double", height=3.0, ncols=2)

    # ---- (a) detection provenance ----
    cats = ["Rubin only", "Roman only", "both"]
    width = 0.8 / max(len(runs), 1)
    for k, r in enumerate(runs):
        d = r.df
        det = (d.detL == 1) | (d.detR == 1) | (d.detJ == 1)
        w = d["W"].to_numpy(float)[det.to_numpy()]
        dl = d.detL.to_numpy()[det.to_numpy()] == 1
        dr = d.detR.to_numpy()[det.to_numpy()] == 1
        shares = [wfrac(dl & ~dr, w), wfrac(dr & ~dl, w), wfrac(dl & dr, w)]
        x = np.arange(len(cats)) + (k - (len(runs) - 1) / 2) * width
        ax1.bar(x, [100 * s for s in shares], width=width * 0.9,
                color=r.colour, label=r.label, linewidth=0)
    ax1.set_xticks(np.arange(len(cats)))
    ax1.set_xticklabels(cats)
    ax1.set_ylabel("share of detected events [%]")
    ps.panel_label(ax1, "(a)")
    ps.legend(ax1, loc="upper right")   # upper LEFT is where the tallest bar and (a) both sit

    # ---- (b) what the joint fit buys, BOTH WAYS ----
    #
    # RESTRICTED TO EVENTS BOTH SURVEYS CHARACTERISED ON THEIR OWN, and that restriction is the
    # whole correctness of the panel. Roman has epochs for only a few per cent of draws -- its
    # footprint is 2% of the scanned area -- so over all characterised events the "joint" fit
    # IS the Rubin fit for the overwhelming majority, and the ratio is exactly 1 by
    # construction. Pooling those in buries the real gain under a tautology: measured on this
    # run, the unrestricted weighted median is 1.000, while the restricted one is 0.14.
    #
    # Two curves per population, because the help is not symmetric and the asymmetry is the
    # result: Roman's dense cadence sharpens what Rubin alone could do, and Rubin's decade-long
    # baseline sharpens what Roman alone could do, by very different factors.
    for r in runs:
        d = r.df
        both = ((d.okA_L == 1) & (d.okA_R == 1)).to_numpy()
        if both.sum() < 50:
            print(f"  fig_synergy(b): {r.name} has {both.sum()} events characterised by both "
                  "surveys -- too few to plot")
            continue
        te = d["tE"].to_numpy(float)[both]
        w = d["W"].to_numpy(float)[both]
        bins = np.geomspace(max(te.min(), 1e-2), te.max(), 9)
        for other, style, lab in (("L", "-", "vs Rubin alone"), ("R", "--", "vs Roman alone")):
            num = d["sigtE_J"].to_numpy(float)[both]
            den = d[f"sigtE_{other}"].to_numpy(float)[both]
            good = (num > 0) & (den > 0)          # -1 is the not-measured sentinel, never a value
            cen, med = [], []
            for i in range(len(bins) - 1):
                m = good & (te >= bins[i]) & (te < bins[i + 1])
                if m.sum() < 20 or w[m].sum() <= 0:
                    continue
                cen.append(np.sqrt(bins[i] * bins[i + 1]))
                med.append(R.weighted_median(num[m] / den[m], w[m]))
            if cen:
                ax2.plot(cen, med, style, color=r.colour, lw=1.4,
                         label=f"{r.name} {lab}")
    ax2.axhline(1.0, color=ps.MUTED, lw=0.7, ls=":")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$t_{\rm E}$ [d]")
    ax2.set_ylabel(r"median $\sigma_{t_{\rm E}}({\rm joint})/\sigma_{t_{\rm E}}({\rm single})$")
    ps.panel_label(ax2, "(b)", loc="lower left")
    ps.legend(ax2, loc="center left")   # upper right is where the "vs Roman" curves run

    ps.stamp(fig, stamp_for(runs))
    return ps.save_figure(fig, f"{out}_synergy")


# ---------------------------------------------------------------------------------------
# 2. The astrometric shift
# ---------------------------------------------------------------------------------------
def fig_astrometry(runs, out):
    """(a) the shift each event reaches; (b) how much of it survives blending.

    The shift is what makes theta_E -- and through it the lens MASS -- measurable without a
    degeneracy. But it is measured against a per-exposure precision of order mas, and it is
    diluted by blending: the centroid is of ALL the light in the aperture, so an unlensed
    blend of fraction (1 - fb) drags the measured shift down by roughly fb.
    """
    fig, (ax1, ax2) = ps.figure(width="double", height=3.0, ncols=2)
    yspan = [np.inf, -np.inf]     # data range of panel (b), for the tick formatter

    for r in runs:
        d = r.df
        det = ((d.detL == 1) | (d.detR == 1) | (d.detJ == 1)).to_numpy()
        shift = max_centroid_shift(d)[det]
        w = d["W"].to_numpy(float)[det]
        good = np.isfinite(shift) & (shift > 0) & (w > 0)
        if good.sum() < 10:
            continue
        s, ww = shift[good], w[good]
        order = np.argsort(s)
        cdf = np.cumsum(ww[order]) / ww.sum()
        ax1.plot(s[order], 100 * (1.0 - cdf), color=r.colour, label=r.label)

        # (b) how the signal scales with the lens mass. delta_theta_max is proportional to
        # theta_E, which goes as sqrt(Ml), so this panel is the mass reach of the astrometric
        # channel -- and it is the panel that says which lenses are worth chasing.
        #
        # Shown against the BLEND-DILUTED shift, because the centroid is of all the light in
        # the aperture: an unlensed blend of fraction (1 - fb) drags the measurement down by
        # roughly fb. On this run the F146 source fraction averages 0.95, so the dilution is
        # small in the median and matters only in the faint tail -- which is exactly why it is
        # drawn as a band rather than asserted to be negligible.
        ml = d["Ml"].to_numpy(float)[det][good]
        fb = np.clip(d["fb1"].to_numpy(float)[det][good], 0.0, 1.0)
        sd = s * fb
        lo, hi = ml.min(), ml.max()
        if hi > lo > 0:
            bins = np.geomspace(lo, hi, mass_bin_count(lo, hi))
            cen, med, q1, q3 = [], [], [], []
            for i in range(len(bins) - 1):
                m = (ml >= bins[i]) & (ml < bins[i + 1])
                if m.sum() < 50 or ww[m].sum() <= 0:
                    continue
                cen.append(np.sqrt(bins[i] * bins[i + 1]))
                med.append(R.weighted_median(sd[m], ww[m]))
                q1.append(R.weighted_quantile(sd[m], ww[m], 0.25))
                q3.append(R.weighted_quantile(sd[m], ww[m], 0.75))
            if cen:
                ax2.plot(cen, med, color=r.colour, label=r.label)
                ax2.fill_between(cen, q1, q3, color=r.colour, alpha=0.18, linewidth=0)
                yspan[0] = min(yspan[0], min(q1))
                yspan[1] = max(yspan[1], max(q3))

    ax1.axvline(ROMAN_AST_FLOOR, color=ps.MUTED, lw=0.7, ls=":")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"max centroid shift $\delta\theta_{\rm c}$ [mas]")
    ax1.set_ylabel("events above this shift [%]")
    ps.panel_label(ax1, "(a)")
    ps.legend(ax1, loc="lower left")

    ax2.axhline(ROMAN_AST_FLOOR, color=ps.MUTED, lw=0.7, ls=":")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"lens mass $M_{\rm L}$ [$M_\odot$]")
    ax2.set_ylabel(r"blend-diluted $\delta\theta_{\rm c}$ [mas]")
    plain_log_ticks(ax2, yspan[0], yspan[1], "y")
    ps.panel_label(ax2, "(b)", loc="lower right")
    ps.legend(ax2, loc="upper left")

    ps.stamp(fig, stamp_for(runs)
             + "   |   dotted: Roman per-exposure floor 1.1 mas; band: weighted 25-75%")
    return ps.save_figure(fig, f"{out}_astrometry")


# ---------------------------------------------------------------------------------------
# 3. Satellite parallax
# ---------------------------------------------------------------------------------------
def fig_parallax(runs, out):
    """What Roman at L2 buys over the same event seen from Earth.

    This is the one comparison in the project that is genuinely controlled: the SAME event is
    characterised twice, once with Roman at L2 and once with the offset zeroed, so the two
    forecasts differ only in the observer's position. Per-event ratios carry NO weight -- the
    weight enters only when the ratios are pooled.
    """
    have = [r for r in runs if r.pair is not None and "Ml" in r.pair.columns]
    if not have:
        print("  fig_parallax: no weightable h3_pair.dat in any run -- skipped")
        return []

    fig, (ax1, ax2) = ps.figure(width="double", height=3.0, ncols=2)

    for r in have:
        p = r.pair
        ok = ((p.okA_sat == 1) & (p.okA_nosat == 1)
              & (p.sigpiE_sat > 0) & (p.sigpiE_nosat > 0)).to_numpy()
        if ok.sum() < 20:
            continue
        gain = (p.sigpiE_nosat.to_numpy(float)[ok] / p.sigpiE_sat.to_numpy(float)[ok])
        w = p["W"].to_numpy(float)[ok]
        g = gain[np.isfinite(gain) & (gain > 0)]
        ww = w[np.isfinite(gain) & (gain > 0)]
        # SURVIVAL, not the CDF. Satellite parallax does nothing for most events -- the ratio
        # is 1 -- so a CDF is a vertical line at 1 that hides the entire result. What matters
        # is the tail: the share of events for which L2 helps by MORE than a given factor.
        order = np.argsort(g)
        surv = 100 * (1.0 - np.cumsum(ww[order]) / ww.sum())
        ax1.plot(g[order], surv, color=r.colour, label=r.label)
        print(f"  parallax {r.name}: share with sigma(piE) better by >1.1x = "
              f"{100 * wfrac(g > 1.1, ww):.2f}%, >2x = {100 * wfrac(g > 2.0, ww):.2f}%")

        # Degeneracy breaking, which is the physical point rather than the precision gain:
        # the lens mass follows from Ml = theta_E/(kappa piE), so a piE that is only bounded
        # leaves the mass unbounded however well tE is measured.
        okm = ok & (p.relMl_sat > 0).to_numpy() & (p.relMl_nosat > 0).to_numpy()
        if okm.sum() >= 20:
            rg = (p.relMl_nosat.to_numpy(float)[okm] / p.relMl_sat.to_numpy(float)[okm])
            wm = p["W"].to_numpy(float)[okm]
            fin = np.isfinite(rg) & (rg > 0)
            order = np.argsort(rg[fin])
            surv = 100 * (1.0 - np.cumsum(wm[fin][order]) / wm[fin].sum())
            ax2.plot(rg[fin][order], surv, color=r.colour, label=r.label)

    for ax, tag, xl in ((ax1, "(a)", r"$\sigma_{\pi_{\rm E}}({\rm Earth})/\sigma_{\pi_{\rm E}}({\rm L2})$"),
                        (ax2, "(b)", r"$\sigma_{M_{\rm L}}({\rm Earth})/\sigma_{M_{\rm L}}({\rm L2})$")):
        ax.axvline(1.0, color=ps.MUTED, lw=0.7, ls=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xl)
        ax.set_ylabel("share of events above this gain [%]")
        ps.panel_label(ax, tag, loc="lower left")
    ps.legend(ax1, loc="upper right")

    ps.stamp(fig, stamp_for(have) + "   |   right of the dotted line: L2 is better")
    return ps.save_figure(fig, f"{out}_parallax")


# ---------------------------------------------------------------------------------------
# 4. Resolving the two images
# ---------------------------------------------------------------------------------------
def fig_resolution(runs, out):
    """P(the two images can be told apart), after Sajadian & Makler (arXiv:2608.16448).

    (a) The probability per survey and per population, at all three bars, because the answer
        moves by more than a factor of four across the paper's own range of the D factor and
        a single bar would hide that.
    (b) The same probability against lens mass, at the most permissive bar. theta_E scales as
        sqrt(Ml), so this is the panel that shows the mass reach of the technique.
    """
    fig, (ax1, ax2) = ps.figure(width="double", height=3.0, ncols=2)

    bars = [("nres5", r"$D{=}5$"), ("nres20", r"$D{=}20$"), ("nresPSF", "PSF")]
    surveys = [("L", "rubin"), ("R", "roman")]

    # Short labels on purpose: "Rubin PSF FWHM" next to "Roman PSF FWHM" collided into mush at
    # column width. The bar meanings are spelled out in the stamp and the caption instead.
    labels, xs, k = [], [], 0
    for bar, blab in bars:
        for suf, skey in surveys:
            labels.append(f"{ps.SURVEY_LABEL[skey].split()[0]}\n{blab}")
            xs.append(k)
            k += 1

    width = 0.8 / max(len(runs), 1)
    for j, r in enumerate(runs):
        d = r.df
        vals = []
        for bar, _ in bars:
            for suf, _ in surveys:
                det = (d[f"det{suf}"] == 1).to_numpy()
                if det.sum() == 0:
                    vals.append(np.nan)
                    continue
                w = d["W"].to_numpy(float)[det]
                n = d[f"{bar}_{suf}"].to_numpy(float)[det]
                f = 100 * wfrac(n >= RESOLVE_MIN_EPOCHS, w)
                # A zero on a log axis is minus infinity, not a short bar. Drawing it would
                # give a spike to the bottom of the frame; nan simply leaves the slot empty,
                # which is the honest rendering of "none in this sample".
                vals.append(f if f > 0 else np.nan)
        x = np.asarray(xs, float) + (j - (len(runs) - 1) / 2) * width
        ax1.bar(x, vals, width=width * 0.9, color=r.colour, label=r.label, linewidth=0)
    ax1.set_xticks(xs)
    ax1.set_xticklabels(labels, fontsize=6)
    ax1.set_ylabel(r"$P(N_{\Delta\theta}\geq 3)$ [%]")
    ax1.set_yscale("log")
    ps.panel_label(ax1, "(a)", loc="lower left")
    ps.legend(ax1, loc="upper right")

    # ---- (b) vs lens mass, most permissive bar ----
    xspan = [np.inf, -np.inf]
    for r in runs:
        d = r.df
        for suf, skey in surveys:
            det = (d[f"det{suf}"] == 1).to_numpy()
            if det.sum() < 50:
                continue
            ml = d["Ml"].to_numpy(float)[det]
            w = d["W"].to_numpy(float)[det]
            n = d[f"nres5_{suf}"].to_numpy(float)[det]
            lo, hi = ml.min(), ml.max()
            if not (hi > lo > 0):
                continue
            bins = np.geomspace(lo, hi, mass_bin_count(lo, hi))
            cen, frac = [], []
            for i in range(len(bins) - 1):
                m = (ml >= bins[i]) & (ml < bins[i + 1])
                # 50, not 20: at 20 the neutron-star panel drew two-point curves that dived
                # vertically off the left edge -- a binning artifact reading as a physical
                # cutoff. A bin too thin to measure should be absent, not drawn.
                if m.sum() < 50 or w[m].sum() <= 0:
                    continue
                cen.append(np.sqrt(bins[i] * bins[i + 1]))
                frac.append(100 * wfrac(n[m] >= RESOLVE_MIN_EPOCHS, w[m]))
            cen, frac = np.asarray(cen), np.asarray(frac)
            vis = frac > 0     # a zero on a log axis is minus infinity, not a data point
            if vis.sum() >= 3:
                # Markers as well as a line: where a population resolves almost never, the
                # curve is a three-bin fragment, and a bare fragment reads as a rendering
                # artifact rather than as three honest measurements.
                ax2.plot(cen[vis], frac[vis],
                         "-o" if suf == "L" else "--s",
                         color=r.colour, lw=1.4, ms=2.5,
                         label=f"{r.name}: {ps.SURVEY_LABEL[skey].split()[0]}")
                xspan[0] = min(xspan[0], cen[vis].min())
                xspan[1] = max(xspan[1], cen[vis].max())
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"lens mass $M_{\rm L}$ [$M_\odot$]")
    ax2.set_ylabel(r"$P(N_{\Delta\theta}\geq 3)$ at $D=5$ [%]")
    # With a single narrow-mass population on the axis (neutron stars span 0.3 decades) the log
    # MINOR ticks get labelled and collide into mush. Same treatment as the astrometry panel.
    plain_log_ticks(ax2, xspan[0], xspan[1], "x")
    # Legend low-right, label top-left: the curves rise left-to-right, so that corner pair is
    # the one both populations leave empty. Upper left put the legend on top of the black-hole
    # curves.
    ps.panel_label(ax2, "(b)", loc="upper left")
    ps.legend(ax2, loc="lower right")

    ps.stamp(fig, stamp_for(runs)
             + "   |   resolvable = both images detectable AND separated by >= the bar, at >= 3"
               " epochs; bars: D=5 and D=20 are D*sigma_a, PSF is the filter's FWHM")
    return ps.save_figure(fig, f"{out}_resolution")


def summary(runs):
    """Print the numbers the figures are made of.

    A figure is for seeing a shape; a number is for quoting in a paper. Every value here is the
    weighted one, over the same cuts the corresponding panel uses, so the two cannot drift apart.
    """
    for r in runs:
        d = r.df
        det = ((d.detL == 1) | (d.detR == 1) | (d.detJ == 1)).to_numpy()
        w = d["W"].to_numpy(float)
        wd = w[det]
        print(f"\n=== {r.name}  ({r.population}) ===")
        print(f"  draws {len(d):,} | detected {int(det.sum()):,} | N_eff {r.neff:,.0f}"
              f" | dropped barren {r.n_barren:,}")

        # --- who detects what ---
        dl = (d.detL == 1).to_numpy()[det]
        dr = (d.detR == 1).to_numpy()[det]
        print(f"  detection share:  Rubin only {100*wfrac(dl & ~dr, wd):6.2f}%"
              f"   Roman only {100*wfrac(dr & ~dl, wd):6.2f}%"
              f"   both {100*wfrac(dl & dr, wd):6.2f}%")

        # --- joint gain, on events BOTH surveys characterised ---
        both = ((d.okA_L == 1) & (d.okA_R == 1)).to_numpy()
        print(f"  characterised by both surveys: {int(both.sum()):,}")
        if both.sum() >= 50:
            wb = w[both]
            for p in ("tE", "piE", "tetE"):
                num = d[f"sig{p}_J"].to_numpy(float)[both]
                out = []
                for other, lab in (("L", "Rubin"), ("R", "Roman")):
                    den = d[f"sig{p}_{other}"].to_numpy(float)[both]
                    g = (num > 0) & (den > 0)
                    out.append(f"vs {lab} {R.weighted_median(num[g]/den[g], wb[g]):.3f}"
                               if g.sum() >= 20 else f"vs {lab}   n/a")
                print(f"    median sigma({p:4s}) joint/single:  " + "   ".join(out))

        # --- astrometric shift ---
        shift = max_centroid_shift(d)[det]
        fb = np.clip(d["fb1"].to_numpy(float)[det], 0.0, 1.0)
        ok = np.isfinite(shift) & (shift > 0)
        if ok.sum():
            print(f"  max centroid shift [mas]: median {R.weighted_median(shift[ok], wd[ok]):.4f}"
                  f"   95th {R.weighted_quantile(shift[ok], wd[ok], 0.95):.4f}"
                  f"   above Roman floor {100*wfrac(shift[ok] > ROMAN_AST_FLOOR, wd[ok]):.2f}%"
                  f"   (blend-diluted {100*wfrac(shift[ok]*fb[ok] > ROMAN_AST_FLOOR, wd[ok]):.2f}%)")

        # --- resolving the two images ---
        for suf, name in (("L", "Rubin"), ("R", "Roman")):
            m = (d[f"det{suf}"] == 1).to_numpy()
            if m.sum() == 0:
                continue
            wm = w[m]
            vals = [100 * wfrac(d[f"{bar}_{suf}"].to_numpy(float)[m] >= RESOLVE_MIN_EPOCHS, wm)
                    for bar in ("nres5", "nres20", "nresPSF")]
            print(f"  P(resolvable) {name:5s}:  D=5 {vals[0]:7.4f}%   D=20 {vals[1]:7.4f}%"
                  f"   PSF {vals[2]:7.4f}%   (of {int(m.sum()):,} detections)")

        # --- satellite parallax ---
        if r.pair is not None and "Ml" in r.pair.columns:
            p = r.pair
            ok = ((p.okA_sat == 1) & (p.okA_nosat == 1)
                  & (p.sigpiE_sat > 0) & (p.sigpiE_nosat > 0)).to_numpy()
            if ok.sum() >= 20:
                g = (p.sigpiE_nosat.to_numpy(float)[ok] / p.sigpiE_sat.to_numpy(float)[ok])
                wp = p["W"].to_numpy(float)[ok]
                f2 = np.isfinite(g) & (g > 0)
                print(f"  satellite parallax: median gain {R.weighted_median(g[f2], wp[f2]):.4f}"
                      f"   >1.1x {100*wfrac(g[f2] > 1.1, wp[f2]):.3f}%"
                      f"   >2x {100*wfrac(g[f2] > 2.0, wp[f2]):.3f}%"
                      f"   (of {int(ok.sum()):,} paired)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, metavar="NAME=DIR",
                    help="population name and its run directory; repeatable")
    ap.add_argument("-o", "--out", default="figures/p6")
    ap.add_argument("--chunksize", type=int, default=500_000)
    ap.add_argument("--unweighted", action="store_true",
                    help="use the raw sample. Says you MEANT an unweighted number; there is "
                         "no silent fallback.")
    a = ap.parse_args()

    ps.use_paper_style()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)

    runs = []
    print("loading:")
    for spec in a.run:
        if "=" not in spec:
            sys.exit(f"--run wants NAME=DIR, got '{spec}'")
        name, directory = spec.split("=", 1)
        runs.append(Run(name, directory, a.chunksize, a.unweighted))

    summary(runs)

    written = []
    written += fig_synergy(runs, a.out)
    written += fig_astrometry(runs, a.out)
    written += fig_parallax(runs, a.out)
    written += fig_resolution(runs, a.out)
    for p in written:
        print(f"wrote {p}")


if __name__ == "__main__":
    main()
