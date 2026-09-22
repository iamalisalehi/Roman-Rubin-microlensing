#!/usr/bin/env python3
"""Step S2: draw one sample event -- light curves and astrometry, with and without parallax.

Reads the three files Step S1's `--dump-samples` writes per event (Deviation 50):

    <class>_<id>_epochs.dat   what Rubin and Roman actually recorded, one row per visit
    <class>_<id>_model.dat    a dense noise-free model, both observer frames, all 7 filters
    <class>_<id>_params.dat   the event's true parameters and its forecast sigmas

and draws one double-column figure per event, 3 x 2 panels under a parameter strip:

    (a) the full decade, dm = m - m_base per filter    (b) the peak, t0 +- 2 tE
    (c) the parallax signal: data - no-parallax model  (d) |centroid shift| vs time
    (e) the astrometric ellipse, shift_2 vs shift_1    (f) sky tracks, proper motion removed

WHAT "WITHOUT PARALLAX" MEANS HERE -- read before writing a caption. It removes the
MICROLENSING parallax only: the bending of the RELATIVE lens-source track by the observer's
orbit, measured by pi_E, which is what gives the lens distance. The source's own annual
parallax pi_s = 1/D_s is an ordinary astrometric wobble and stays in the sky track either way;
panel (f) shows it, labelled as such. Conflating the two in a caption is a physics error.

WHICH "WITHOUT PARALLAX" -- the gauge, and why the simulation's own version is not used.
lightcurve() measures the observer's displacement from Earth's position at t = 0, the START
of the simulation, and its no-parallax trajectory (u_noplx, def*a, mag0_*) is the straight
line in that gauge. For an event years later, that line is offset from the true trajectory by
pi_E times the Earth's displacement since t = 0 -- 1.75 AU projected for the first event this
script drew -- so the "no-parallax" curve peaks at a different time and height from the
event itself. Plotted as the parallax signal, it is dominated by the constant and linear
parts of the Earth's motion, which any real fit absorbs into u0, t0, tE and the direction of
motion; it would overstate parallax by orders of magnitude.

The default here (--noplx geocentric) is therefore the standard geocentric frame (Gould
2004): the no-parallax model is the straight line matching the TRUE trajectory's position and
velocity at its own peak, in each observer's frame. What remains between the two curves is the
observer's acceleration -- the only part of the parallax a fit can measure. The trajectory
vector is recovered exactly from the dump as u_vec = def_c (u^2 + 2) / thetaE. --noplx
simulation draws the simulation's own t = 0-gauge curves, for comparison only.

None of this touches a result: the simulation's no-parallax chi-square feeds only a dead
diagnostic file, and the marginalised sigma(pi_E) is invariant to the reparametrisation.

THREE DISPLAY CHOICES, AND WHY. Each was agreed before this was written.

  dm, not m. The bands sit at baselines several magnitudes apart, so on an absolute axis each
  is a flat line in its own strip. dm overlays them, and makes CHROMATIC BLENDING visible:
  a band's amplitude is set by the source's share of the light in it, blend[i], which differs
  per band. Different amplitudes per band are physics, not a plotting artefact.

  Binned astrometry. A per-visit astrometric error is ~5 mas for Roman and ~10 mas for Rubin
  at typical bulge magnitudes; a typical centroid shift is 0.1-4 mas. Drawn per epoch, the
  curve disappears inside its own error bars. Astrometric microlensing is detected by
  averaging many exposures, and the binned points show that averaging: inverse-variance mean
  per bin, error 1/sqrt(sum 1/sigma^2) -- i.e. sigma/sqrt(N) for equal errors. Two assumptions
  ride on that and belong in the caption: exposures are treated as independent (which
  OPEN_ITEMS.md questions for Roman's 1.1 mas floor), and the scalar `err_ast` is taken as the
  per-axis 1 sigma, as the simulation's own chi-square uses it.

  Astrometric points sit ON the model. The simulation never draws a 2D astrometric
  measurement -- its chi-square uses the scalar |position| and one scalar noise draw -- so there
  is no noisy datum to plot, and inventing one here would be a figure of noise we made up. The
  points are the model at the observed epochs, binned, with the instrument's error bars.
  Photometry is different: mag_obs IS the simulation's own draw, so panels (a)-(c) are real
  simulated data.

NOT WEIGHTED, DELIBERATELY. The event-rate weight exists so that a pooled fraction or median
describes the sky rather than the Monte Carlo sample. A figure of one event pools nothing, so
there is nothing to weight. The one romanlib rule that does apply is the sentinel: a sigma of
-1.0 is "not measured" and is printed as such, never as a number.

USAGE
    .roman/bin/python analysis/s2_sample_lightcurves.py samples/bh -o figures/samples/bh \\
        [--prov runs/<dir>/files/MONTLMC/files/run_provenance.txt] [--only both_003]
"""

import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R          # noqa: E402
import plotstyle as ps        # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402  (after plotstyle, which selects the backend)
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

YEAR = 365.2425
FILTERS = ["u", "g", "r", "i", "z", "y", "F146"]
ROMAN_FILT = 6

# Bin widths [days]. Roman observes a field every 12.1 min in season (~119 exposures a day),
# Rubin every few days in some band. The widths are chosen so a bin holds enough exposures to
# beat the per-visit error down, and short enough not to smear the peak of a ~10 d event.
BIN_ROMAN_PHOT = 0.5
# Astrometric bins are ADAPTIVE, per event: wide enough that a bin's error is ~1/4 of the peak
# centroid shift, capped at tE/3 (and 30 d) so the curve is not smeared out of existence. A
# faint source can hit the cap with error bars still larger than its signal -- which is the
# physics: that event's shift is not resolvable in time, only across a season. The width
# actually used is printed in the legend. Bins with fewer than this many exposures are
# dropped: a one-exposure bin at a season edge is a 5 mas error bar carrying no information.
ASTR_TARGET_FRACTION = 0.25
ASTR_BIN_MAX_DAYS = 30.0
ASTR_BIN_MIN_N = 3
# ...and no more than this many Roman bins across the window. A bright source with ~50,000
# exposures reaches the error target in a couple of days, and hundreds of individually
# informative crosses overplot into a solid block. Widening the bins MERGES information
# (each bin's error shrinks) rather than hiding any.
ASTR_MAX_BINS = 60
# A season observed at more than this many exposures a day is high-cadence (~119/d at the
# 12.1-min cycle); the ROTAC low-cadence seasons are one ~1.5 h unit every 5 days.
HIGH_CADENCE_PER_DAY = 10.0

# Consecutive Roman epochs further apart than this start a new season. The same threshold
# RomanSchedule uses in C++ (SEASON_GAP_MIN_DAYS in Bulge.h), so the shaded windows are the
# seasons the simulation itself believed in.
SEASON_GAP_MIN_DAYS = 20.0

# Largest possible point-lens centroid shift: u*thetaE/(u^2+2) peaks at u = sqrt(2).
SHIFT_MAX_OVER_TETE = 1.0 / (2.0 * np.sqrt(2.0))

EPOCH_COLS = ["t", "tele", "filt", "mag_obs", "mag_mod", "mag_mod0", "err_mag", "u", "u0",
              "A", "A0", "def1c", "def2c", "def1a", "def2a", "pos1b", "pos2b", "pos1c",
              "pos2c", "lens1", "lens2", "err_ast"]
MODEL_COLS = (["t", "frame", "u", "u0", "A", "A0"]
              + [f"mag_{f}" for f in FILTERS] + [f"mag0_{f}" for f in FILTERS]
              + ["def1c", "def2c", "def1a", "def2a", "pos1b", "pos2b", "pos1c", "pos2c",
                 "lens1", "lens2"])


def band_colours():
    """Rubin's six bands in wavelength order along viridis, F146 in the Roman survey colour.

    viridis is perceptually ordered and colourblind-safe, and has no orange in it, so F146
    cannot be mistaken for a Rubin band. Every band also gets its own marker, so the figure
    does not depend on colour alone.
    """
    rubin = plt.cm.viridis(np.linspace(0.0, 0.82, 6))
    cols = {f: rubin[i] for i, f in enumerate(FILTERS[:6])}
    cols["F146"] = ps.SURVEY["roman"]
    return cols


MARKERS = {"u": "v", "g": "^", "r": "o", "i": "s", "z": "D", "y": "P", "F146": "."}


# ---------------------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------------------
def read_params(path):
    """The key/value block. Three line shapes occur: `key value`, `key v1 v2 ...` (magb,
    blend: one value per FILTER), and `k1 v1 k2 v2 ...` (the per-survey flags and sigmas)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].split()
            if not line:
                continue
            if len(line) == 2:
                out[line[0]] = _num(line[1])
            elif len(line) > 2 and not _isnum(line[2]):
                for k, v in zip(line[0::2], line[1::2]):
                    out[k] = _num(v)
            else:
                out[line[0]] = np.array([float(v) for v in line[1:]])
    return out


def _isnum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def _num(s):
    return float(s) if _isnum(s) else s


def read_table(path, cols):
    a = np.loadtxt(path, comments="#", ndmin=2)
    if a.shape[1] != len(cols):
        sys.exit(f"{path}: {a.shape[1]} columns, expected {len(cols)}. Was it written by a "
                 f"different version of --dump-samples?")
    return {c: a[:, i] for i, c in enumerate(cols)}


def subset(tab, mask):
    return {k: v[mask] for k, v in tab.items()}


# ---------------------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------------------
def seasons_from_epochs(t_roman):
    """Roman's observing windows, as (start, end) pairs, from this sightline's own epochs."""
    if t_roman.size == 0:
        return []
    t = np.unique(t_roman)
    cut = np.where(np.diff(t) > SEASON_GAP_MIN_DAYS)[0]
    starts = np.r_[t[0], t[cut + 1]]
    ends = np.r_[t[cut], t[-1]]
    return list(zip(starts, ends))


def classify_seasons(t_roman):
    """Seasons as (start, end, is_high_cadence)."""
    out = []
    for a, b in seasons_from_epochs(t_roman):
        n = np.count_nonzero((t_roman >= a) & (t_roman <= b))
        out.append((a, b, n / max(b - a, 1.0) > HIGH_CADENCE_PER_DAY))
    return out


def trajectory_vector(tab, tetE):
    """(u1, u2): the lens-source separation VECTOR in Einstein radii, recovered from the
    centroid shift, which lightcurve() writes as u_vec * thetaE / (u^2 + 2)."""
    k = (tab["u"] ** 2 + 2.0) / tetE
    return tab["def1c"] * k, tab["def2c"] * k


class Geocentric:
    """The standard no-parallax model for one observer frame: the straight line tangent to the
    true (parallax-bent) trajectory at the moment of closest approach (Gould 2004). Built from
    the dense model of that frame; evaluable at any time."""

    def __init__(self, dense, tetE):
        u1, u2 = trajectory_vector(dense, tetE)
        t = dense["t"]
        i = int(np.argmin(dense["u"]))
        i = min(max(i, 1), t.size - 2)
        self.tref = t[i]
        self.w = np.array([u1[i], u2[i]])
        # central difference on the dense grid, which is sampled at tE/100 around the peak
        self.v = np.array([(u1[i + 1] - u1[i - 1]) / (t[i + 1] - t[i - 1]),
                           (u2[i + 1] - u2[i - 1]) / (t[i + 1] - t[i - 1])])
        self.tetE = tetE

    def u_vec(self, t):
        dt = np.asarray(t) - self.tref
        return self.w[0] + self.v[0] * dt, self.w[1] + self.v[1] * dt

    def u(self, t):
        return np.hypot(*self.u_vec(t))

    def shift(self, t):
        a, b = self.u_vec(t)
        k = self.tetE / (a * a + b * b + 2.0)
        return a * k, b * k

    def A(self, t):
        uu = self.u(t)
        return (uu * uu + 2.0) / np.sqrt(uu * uu * (uu * uu + 4.0))

    def mag(self, t, magb, blend):
        return magb - 2.5 * np.log10(self.A(t) * blend + 1.0 - blend)


def mag_to_A(m, sig_m, magb, blend):
    """A measured magnitude as the SOURCE's magnification, and its 1-sigma.

    m = m_base - 2.5 log10(A fb + 1 - fb), inverted. This is what lets every filter sit on one
    curve (Figure 2 of Sajadian & Makler plots magnification for exactly this reason): the
    band-to-band differences in dm are all blending, and dividing it out leaves A, which is
    achromatic for a point lens. The price is the error bar, which grows as 1/fb -- a band in
    which the source is 1% of the light has its magnification error inflated 100x.
    """
    f = 10.0 ** (-0.4 * (m - magb))                 # flux relative to the blended baseline
    A = (f - 1.0 + blend) / blend
    sA = 0.4 * np.log(10.0) * f * sig_m / blend
    return A, sA


def straight_track(t, pos, mu, t0):
    """A body's sky track WITHOUT its annual parallax: the straight line mu (t - t0) + c.

    The parallax term in lightcurve() is -pi * (P(X(t)) - P(X(0))): the observer's projected
    orbit, which for the circular orbit the code uses is a pure first harmonic in time, offset
    by the constant P(X(0)). So the residual pos - mu (t - t0) is exactly c + a cos wt + b sin wt,
    and a linear least-squares fit returns the centre c of the parallax ellipse exactly --
    unlike a plain mean, which is biased over a non-integer number of years.
    """
    w = 2.0 * np.pi / YEAR
    r = pos - mu * (t - t0)
    X = np.column_stack([np.ones_like(t), np.cos(w * t), np.sin(w * t)])
    c = np.linalg.lstsq(X, r, rcond=None)[0][0]
    return mu * (t - t0) + c


def emptiest_corner(ax, xs, ys):
    """The axes corner holding the fewest curve points, for an inset or a legend."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    fx = (np.asarray(xs) - x0) / (x1 - x0)
    fy = (np.asarray(ys) - y0) / (y1 - y0)
    boxes = {"upper left": (0, .5, .5, 1), "upper right": (.5, 1, .5, 1),
             "lower left": (0, .5, 0, .5), "lower right": (.5, 1, 0, .5)}
    count = {k: np.count_nonzero((fx >= a) & (fx < b) & (fy >= c) & (fy < d))
             for k, (a, b, c, d) in boxes.items()}
    return min(count, key=count.get)


def omitted_note(ax, n, what="points", loc="lower right"):
    """Say, on the panel, how many points were left out for carrying no information."""
    if n > 0:
        x, ha = (0.98, "right") if "right" in loc else (0.03, "left")
        y, va = (0.03, "bottom") if "lower" in loc else (0.97, "top")
        what = what if n != 1 else what.rstrip("s")
        ax.text(x, y, f"{n} {what} with $\\sigma$ > signal/3 not shown",
                transform=ax.transAxes, ha=ha, va=va, fontsize=5.5, color=ps.MUTED,
                zorder=10, bbox=dict(fc=ps.SURFACE, ec="none", alpha=0.85, pad=0.8))


def astr_bin_width(err_ast, rate_per_day, peak_shift, tE):
    """Adaptive astrometric bin width [days] -- see ASTR_TARGET_FRACTION."""
    if err_ast.size == 0 or peak_shift <= 0 or rate_per_day <= 0:
        return 1.0
    target = ASTR_TARGET_FRACTION * peak_shift
    n_needed = (np.median(err_ast) / target) ** 2
    cap = min(ASTR_BIN_MAX_DAYS, max(1.0, tE / 3.0))
    return float(np.clip(n_needed / rate_per_day, 1.0, cap))


def bin_ivar(t, y, err, width, min_n=1):
    """Inverse-variance mean of y in bins of `width` days. Returns (t, y, sigma, n) per bin."""
    if t.size == 0:
        return (np.array([]),) * 4
    key = np.floor(t / width).astype(np.int64)
    order = np.argsort(key, kind="stable")
    key, t, y, err = key[order], t[order], y[order], err[order]
    edges = np.r_[0, np.where(np.diff(key) != 0)[0] + 1, key.size]
    w = 1.0 / err ** 2
    tb, yb, sb, nb = [], [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        ws = w[a:b].sum()
        tb.append((t[a:b] * w[a:b]).sum() / ws)
        yb.append((y[a:b] * w[a:b]).sum() / ws)
        sb.append(1.0 / np.sqrt(ws))
        nb.append(b - a)
    tb, yb, sb, nb = (np.array(x) for x in (tb, yb, sb, nb))
    keep = nb >= min_n
    return tb[keep], yb[keep], sb[keep], nb[keep]


def fmt_sigma(p, key, true_value):
    """sigma/true as a percentage, or 'n/m' for the not-measured sentinel."""
    s = p.get(key, R.SENTINEL)
    if not np.isfinite(s) or s == R.SENTINEL or true_value <= 0:
        return "n/m"
    pct = 100.0 * s / true_value
    if pct >= 1000:
        return ">1000%"          # measured, in the sense that the matrix inverted; useless
    return f"{pct:.2g}%" if pct < 10 else f"{pct:.0f}%"


# ---------------------------------------------------------------------------------------
# Physics sanity checks -- asserted, not just drawn
# ---------------------------------------------------------------------------------------
def sanity_checks(p, ep, mo):
    """Three checks that the dump is internally consistent. Returns printable lines; raises
    on a hard violation. These are the checks a figure cannot show you are failing."""
    lines = []
    m0 = subset(mo, mo["frame"] == 0)
    m1 = subset(mo, mo["frame"] == 1)

    # 1. The epoch file's model magnitude must agree with the dense model at the same instant,
    #    in the same frame and filter. Two independent writers of one quantity: if they
    #    disagree, one of them read stale scratch state.
    worst = 0.0
    for tele, mf in ((0, m0), (1, m1)):
        e = subset(ep, ep["tele"] == tele)
        for fi in np.unique(e["filt"]).astype(int):
            sel = e["filt"] == fi
            interp = np.interp(e["t"][sel], mf["t"], mf[f"mag_{FILTERS[fi]}"])
            worst = max(worst, float(np.max(np.abs(interp - e["mag_mod"][sel]))))
    lines.append(f"epoch-vs-model magnitude, worst |diff| = {worst:.2e} mag "
                 f"(interpolation error only; limit 1e-2)")
    if worst > 1e-2:
        raise AssertionError(f"epoch and model files disagree by {worst:.3g} mag")

    # 2. No centroid shift may exceed the point-lens maximum, thetaE / (2 sqrt 2).
    cap = SHIFT_MAX_OVER_TETE * p["tetE"]
    smax = max(float(np.max(np.hypot(mo["def1c"], mo["def2c"]))),
               float(np.max(np.hypot(ep["def1c"], ep["def2c"]))) if ep["t"].size else 0.0)
    lines.append(f"max |shift| = {smax:.4g} mas vs point-lens cap {cap:.4g} mas")
    # Tolerance is the dump's own precision: positions are written to 8 significant digits
    # and thetaE to 10, so an event that passes through u = sqrt(2) -- every event with
    # u0 < sqrt(2) does -- sits ON the cap and can round a few parts in 1e8 above it.
    if smax > cap * (1.0 + 1e-6):
        raise AssertionError(f"centroid shift {smax} exceeds thetaE/(2 sqrt 2) = {cap}")

    # 3. The parallax gauge. lightcurve() measures the observer's displacement FROM EARTH'S
    #    POSITION AT t = 0, so in the geocentric frame u and u_noplx must coincide there. In
    #    the L2 frame they must NOT: the residual is the satellite-parallax offset, ~1e-3
    #    (the gauge-trap note in lightcurve()).
    i0 = int(np.argmin(np.abs(m0["t"])))
    i1 = int(np.argmin(np.abs(m1["t"])))
    g0 = abs(m0["u"][i0] - m0["u0"][i0])
    g1 = abs(m1["u"][i1] - m1["u0"][i1])
    lines.append(f"gauge at t=0: |u - u_noplx| = {g0:.2e} (Earth, must be ~0), "
                 f"{g1:.2e} (L2, the satellite offset)")
    if m0["t"][i0] == 0.0 and g0 > 1e-9:
        raise AssertionError(f"parallax gauge broken: |u - u0| = {g0} at t = 0 (Earth)")
    return lines


# ---------------------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------------------
# Colours for the sky-trajectory panel. Roles follow Figure 2 of Sajadian & Makler
# (arXiv:2608.16448) -- cyan undeflected source, blue deflected, magenta lens, red lens without
# parallax, dark red the relative track, black the deflection -- so the two figures can be read
# against each other. Line style carries the with/without-parallax distinction everywhere in
# THIS figure: solid = with, dashed = without. (Figure 2 uses the opposite convention on its
# light curves; consistency within one figure wins.)
TRACK = {"src_u": "#0891b2", "src_d": "#1e3a8a", "lens": "#c026d3", "lens0": "#dc2626",
         "rel": "#7f1d1d", "defl": "#1a1a1a"}
TOBS = 10.0 * YEAR     # Rubin's 10-year window; the tracks are dotted outside it

# ONE RULE FOR EVERY PANEL: a point or bin is drawn only if its 1-sigma is below this fraction
# of the signal it is plotted against -- the event's amplitude, the parallax signal, the peak
# centroid shift. At 1/3, every point shown distinguishes that signal at 3 sigma ON ITS OWN;
# a point that cannot is not information about this curve, and a bar that dwarfs the curve
# only hides it. Each panel states how many it left out.
INFO_FRACTION = 1.0 / 3.0


def draw_event(stem, out_dir, prov, noplx="geocentric"):
    p = read_params(stem + "_params.dat")
    ep = read_table(stem + "_epochs.dat", EPOCH_COLS)
    mo = read_table(stem + "_model.dat", MODEL_COLS)
    checks = sanity_checks(p, ep, mo)

    cols = band_colours()
    magb, blend, tetE = p["magb"], p["blend"], p["tetE"]
    t0, tE = p["t0"], p["tE"]
    rub = subset(ep, ep["tele"] == 0)
    rom = subset(ep, ep["tele"] == 1)
    frames = {0: subset(mo, mo["frame"] == 0), 1: subset(mo, mo["frame"] == 1)}
    geo = {f: Geocentric(frames[f], tetE) for f in (0, 1)}
    seasons = classify_seasons(rom["t"])

    # The no-parallax model, per observer frame. See the module docstring: geocentric is the
    # default and the only one to publish.
    def A_noplx(frame, t):
        if noplx == "geocentric":
            return geo[frame].A(t)
        return np.interp(t, frames[frame]["t"], frames[frame]["A0"])

    def m_noplx(frame, i, t, tab=None):
        if noplx == "geocentric":
            return geo[frame].mag(t, magb[i], blend[i])
        if tab is not None:
            return tab["mag_mod0"]
        return np.interp(t, frames[frame]["t"], frames[frame][f"mag0_{FILTERS[i]}"])

    def shift_noplx(frame, t, tab):
        if noplx == "geocentric":
            return geo[frame].shift(t)
        return tab["def1a"], tab["def2a"]

    # Reference time: the OBSERVED peak (Earth frame), not the parameter t0 (OPEN_ITEMS.md).
    tpk = geo[0].tref
    m0, m1 = frames[0], frames[1]
    t_lo, t_hi = float(mo["t"].min()), float(mo["t"].max())
    win = (max(tpk - 3.0 * tE, t_lo), min(tpk + 3.0 * tE, t_hi))
    zoom = (max(tpk - 2.0 * tE, t_lo), min(tpk + 2.0 * tE, t_hi))
    checks.append(f"observed peak (Earth) at t - t0 = {tpk - t0:+.2f} d, u_min = "
                  f"{m0['u'].min():.4f}; parameter u0 = {p['u0']:.4f}")
    # Rubin and Roman see different magnifications (the satellite parallax). The A panels
    # draw the Earth-frame curves; the L2 curve is added only if it would be distinguishable.
    dA = float(np.max(np.abs(m1["A"] - m0["A"]) / m0["A"]))
    draw_l2 = dA > 1e-2
    checks.append(f"satellite parallax: max |A_L2 - A_Earth| / A = {dA:.2e}"
                  + ("  -> L2 curve drawn separately" if draw_l2 else ""))

    # Panel letters, row by row. Rows 1-2 are Figure 2 of Sajadian & Makler, extended: the
    # magnification (left) beside the sky-plane trajectories (right), first over the decade
    # and then at the closest approach. Rows 3-5 are the diagnostics.
    #   (a) A(t), decade        (b) sky trajectories, decade
    #   (c) A(t), peak          (d) sky trajectories, closest approach
    #   (e) dm per band, peak   (f) the parallax signal
    #   (g) |shift| vs time     (h) the astrometric ellipse
    #   (i) proper motion off   (j) key for the trajectory panels
    fig = plt.figure(figsize=(ps.WIDTH["double"], 14.2))
    outer = GridSpec(2, 1, figure=fig, height_ratios=[0.40, 5.6], hspace=0.05,
                     top=0.995, bottom=0.03)
    head = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0], height_ratios=[3.0, 0.8],
                                   hspace=0.0)
    body = GridSpecFromSubplotSpec(5, 2, subplot_spec=outer[1],
                                   height_ratios=[1.25, 1.25, 1, 1.05, 1.05],
                                   hspace=0.33, wspace=0.26)
    strip = fig.add_subplot(head[0])
    key1 = fig.add_subplot(head[1])
    ax_a, ax_b = fig.add_subplot(body[0, 0]), fig.add_subplot(body[0, 1])
    ax_c, ax_d = fig.add_subplot(body[1, 0]), fig.add_subplot(body[1, 1])
    ax_e, ax_f = fig.add_subplot(body[2, 0]), fig.add_subplot(body[2, 1])
    ax_g, ax_h = fig.add_subplot(body[3, 0]), fig.add_subplot(body[3, 1])
    ax_i, key2 = fig.add_subplot(body[4, 0]), fig.add_subplot(body[4, 1])

    # ---- parameter strip: data, not a title --------------------------------------------
    strip.axis("off")
    yes = lambda v: "yes" if int(v) == 1 else "no"          # noqa: E731
    zone = {0: "in a Roman season", 1: "in a mid-mission gap", 2: "outside Roman's mission"}
    in_season = any(a <= tpk <= b for a, b, _ in seasons)
    pk_zone = "in a Roman season" if in_season else "outside Roman's seasons"
    l1 = (rf"$M_L = {p['Ml']:.3g}\,M_\odot$   $D_L = {p['Dl']:.2f}$ kpc   "
          rf"$D_S = {p['Ds']:.2f}$ kpc   $t_E = {tE:.3g}$ d   $\pi_E = {p['piE']:.3g}$   "
          rf"$\theta_E = {tetE:.3g}$ mas   $\mu_{{rel}} = {p['murel_yr']:.3g}$ mas/yr   "
          rf"$u_{{min}} = {m0['u'].min():.3g}$")
    l2 = (f"detected by Rubin: {yes(p['detL'])}   Roman: {yes(p['detR'])}   "
          f"joint: {yes(p['detJ'])}      observed peak {pk_zone} "
          f"(table's $t_0$: {zone.get(int(p['t0zone']), '?')})")
    l2b = (f"epochs Rubin {int(p['ndw_L'])}, Roman {int(p['ndw_R'])};  within "
           f"$\\pm2t_E$ of $t_0$: Rubin {int(p['nep_pk_L'])}, Roman {int(p['nep_pk_R'])}")
    l3 = ("forecast 1$\\sigma$ (joint / Rubin / Roman):   "
          f"$t_E$ {fmt_sigma(p, 'sigtE_J', tE)} / {fmt_sigma(p, 'sigtE_L', tE)} / "
          f"{fmt_sigma(p, 'sigtE_R', tE)}     "
          f"$\\pi_E$ {fmt_sigma(p, 'sigpiE_J', p['piE'])} / "
          f"{fmt_sigma(p, 'sigpiE_L', p['piE'])} / {fmt_sigma(p, 'sigpiE_R', p['piE'])}     "
          f"$\\theta_E$ {fmt_sigma(p, 'sigtetE_J', tetE)} / "
          f"{fmt_sigma(p, 'sigtetE_L', tetE)} / {fmt_sigma(p, 'sigtetE_R', tetE)}"
          "     (n/m = not measured)")
    for k, txt in enumerate((l1, l2, l2b, l3)):
        strip.text(0.0, 1.0 - 0.26 * k, txt, transform=strip.transAxes, fontsize=7,
                   color=ps.INK if k < 3 else ps.MUTED, va="top", ha="left")

    # ---- two keys: photometry, then the sky trajectories -------------------------------
    nolab = "without parallax" + (" (geocentric)" if noplx == "geocentric"
                                  else " (simulation t = 0 gauge)")
    for key in (key1, key2):
        key.axis("off")
    h1 = [Line2D([], [], color=cols[f], marker=MARKERS[f], ls="", ms=3.5, label=f)
          for f in FILTERS]
    h1 += [Line2D([], [], color=ps.INK, lw=0.9, label="model, with parallax"),
           Line2D([], [], color=ps.INK, lw=0.8, ls=(0, (3, 2)), label=nolab)]
    key1.legend(handles=h1, ncol=len(h1), loc="center", frameon=False, fontsize=6.5,
                handlelength=1.5, columnspacing=1.0, handletextpad=0.4, labelcolor=ps.INK)
    h2 = [Line2D([], [], color=TRACK["src_u"], lw=1.0, label="source, undeflected"),
          Line2D([], [], color=TRACK["src_d"], lw=1.0, label="source, deflected"),
          Line2D([], [], color=TRACK["lens"], lw=1.0, label="lens"),
          Line2D([], [], color=TRACK["lens0"], lw=0.9, ls=(0, (3, 2)),
                 label="lens, no parallax"),
          Line2D([], [], color=TRACK["rel"], lw=1.0, label="lens $-$ source"),
          Line2D([], [], color=TRACK["defl"], lw=1.0, label="deflection"),
          Line2D([], [], color=ps.MUTED, lw=0.9, ls=":", label="outside Rubin's 10 yr")]
    h2.insert(6, Line2D([], [], color=TRACK["defl"], lw=0.8, ls=(0, (3, 2)),
                        label="deflection, no microlensing parallax"))
    h2.insert(7, Line2D([], [], color=TRACK["src_d"], lw=0.8, ls=(0, (3, 2)),
                        label="source deflected, no microlensing parallax (d)"))
    h2.append(Line2D([], [], color=ps.SURVEY["roman"], marker="o", ls="", ms=3,
                     label="Roman positions, binned (d)"))
    h2.append(Line2D([], [], color=TRACK["lens0"], marker="x", ls="", ms=4, mew=1.0,
                     label="no parallax: loop centres (i)"))
    key2.legend(handles=h2, ncol=1, loc="center left", frameon=False, fontsize=6.5,
                handlelength=2.2, handletextpad=0.6, labelcolor=ps.INK,
                title="trajectory panels (b), (d), (i)", title_fontsize=6.5)

    def shade(ax, to_x=lambda t: t):
        for a, b, hi in seasons:
            ax.axvspan(to_x(a), to_x(b), color=ps.SURVEY["roman"], alpha=0.10 if hi else 0.04,
                       lw=0, zorder=0)

    # ---- (a) magnification over the decade, (c) at the peak ----------------------------
    # One curve pair for every filter: dividing out each band's blending leaves the source
    # magnification A, which a point lens makes achromatic. Error bars grow as 1/fb, so a
    # heavily blended band's points mostly fail INFO_FRACTION of the amplitude A_peak - 1.
    amp_A = float(m0["A"].max() - 1.0)

    def magnification(ax, xform, xlim):
        n_out = 0
        sel = (m0["t"] >= xlim[0]) & (m0["t"] <= xlim[1])
        tt = m0["t"][sel]
        ax.plot(xform(tt), m0["A"][sel], color=ps.INK, lw=1.0, zorder=4)
        ax.plot(xform(tt), A_noplx(0, tt), color=ps.INK, lw=0.9, ls=(0, (3, 2)), zorder=4)
        if draw_l2:
            s1 = (m1["t"] >= xlim[0]) & (m1["t"] <= xlim[1])
            ax.plot(xform(m1["t"][s1]), m1["A"][s1], color=cols["F146"], lw=0.8, zorder=4)
        for i, f in enumerate(FILTERS[:6]):
            s = (rub["filt"] == i) & (rub["t"] >= xlim[0]) & (rub["t"] <= xlim[1])
            if not s.any():
                continue
            A, sA = mag_to_A(rub["mag_obs"][s], rub["err_mag"][s], magb[i], blend[i])
            ok = sA < INFO_FRACTION * amp_A
            n_out += np.count_nonzero(~ok)
            if ok.any():
                # Drawn light: ~2,400 Rubin points, most of them informative only in bulk.
                ax.errorbar(xform(rub["t"][s][ok]), A[ok], yerr=sA[ok], fmt=MARKERS[f],
                            ms=1.8, mew=0, color=cols[f], ecolor=cols[f], elinewidth=0.3,
                            alpha=0.45, zorder=5)
        s = (rom["t"] >= xlim[0]) & (rom["t"] <= xlim[1])
        if s.any():
            A, sA = mag_to_A(rom["mag_obs"][s], rom["err_mag"][s], magb[ROMAN_FILT],
                             blend[ROMAN_FILT])
            ok = sA < INFO_FRACTION * amp_A
            n_out += np.count_nonzero(~ok)
            ax.scatter(xform(rom["t"][s][ok]), A[ok], s=0.3, color=cols["F146"], alpha=0.12,
                       lw=0, rasterized=True, zorder=2)
            tb, yb, sb, _ = bin_ivar(rom["t"][s][ok], A[ok], sA[ok], BIN_ROMAN_PHOT)
            ax.errorbar(xform(tb), yb, yerr=sb, fmt="o", ms=1.8, mew=0, color=cols["F146"],
                        ecolor=cols["F146"], elinewidth=0.5, zorder=6)
        top = float(m0["A"][sel].max()) if sel.any() else 1.0 + amp_A
        ax.set_ylim(1.0 - 0.06 * (top - 1.0), top + 0.10 * (top - 1.0))
        omitted_note(ax, n_out)

    magnification(ax_a, lambda t: t / YEAR, (t_lo, t_hi))
    shade(ax_a, lambda t: t / YEAR)
    ax_a.set_xlim(t_lo / YEAR, t_hi / YEAR)
    ax_a.set_xlabel("time since first Rubin bulge visit (yr)")
    ax_a.set_ylabel("magnification $A$")
    ps.panel_label(ax_a, "(a)", loc="upper left")

    magnification(ax_c, lambda t: t - tpk, zoom)
    shade(ax_c, lambda t: t - tpk)
    ax_c.set_xlim(zoom[0] - tpk, zoom[1] - tpk)
    ax_c.set_xlabel(r"$t - t_{\rm peak}$ (d)")
    ax_c.set_ylabel("magnification $A$")
    ps.panel_label(ax_c, "(c)", loc="upper left")

    # ---- (e) dm per band at the peak: chromatic blending made visible -------------------
    n_out = 0
    for i, f in enumerate(FILTERS):
        fr = 1 if i == ROMAN_FILT else 0
        mf = frames[fr]
        sel = (mf["t"] >= zoom[0]) & (mf["t"] <= zoom[1])
        tt = mf["t"][sel]
        ax_e.plot(tt - tpk, mf[f"mag_{f}"][sel] - magb[i], color=cols[f], lw=0.9, zorder=3)
        ax_e.plot(tt - tpk, m_noplx(fr, i, tt) - magb[i], color=cols[f], lw=0.8,
                  ls=(0, (3, 2)), zorder=3)
    for i, f in enumerate(FILTERS[:6]):
        amp = float(np.max(np.abs(m0[f"mag_{f}"] - magb[i])))
        s = (rub["filt"] == i) & (rub["t"] >= zoom[0]) & (rub["t"] <= zoom[1])
        ok = s & (rub["err_mag"] < INFO_FRACTION * amp)
        n_out += np.count_nonzero(s & ~ok)
        if ok.any():
            ax_e.errorbar(rub["t"][ok] - tpk, rub["mag_obs"][ok] - magb[i],
                          yerr=rub["err_mag"][ok], fmt=MARKERS[f], ms=2.2, mew=0,
                          color=cols[f], ecolor=cols[f], elinewidth=0.4, alpha=0.8, zorder=4)
    s = (rom["t"] >= zoom[0]) & (rom["t"] <= zoom[1])
    if s.any():
        amp = float(np.max(np.abs(m1["mag_F146"] - magb[ROMAN_FILT])))
        ok = s & (rom["err_mag"] < INFO_FRACTION * amp)
        n_out += np.count_nonzero(s & ~ok)
        ax_e.scatter(rom["t"][ok] - tpk, rom["mag_obs"][ok] - magb[ROMAN_FILT], s=0.3,
                     color=cols["F146"], alpha=0.12, lw=0, rasterized=True, zorder=2)
        tb, yb, sb, _ = bin_ivar(rom["t"][ok], rom["mag_obs"][ok], rom["err_mag"][ok],
                                 BIN_ROMAN_PHOT)
        ax_e.errorbar(tb - tpk, yb - magb[ROMAN_FILT], yerr=sb, fmt="o", ms=1.8, mew=0,
                      color=cols["F146"], ecolor=cols["F146"], elinewidth=0.5, zorder=5)
    ax_e.invert_yaxis()
    shade(ax_e, lambda t: t - tpk)
    ax_e.set_xlim(zoom[0] - tpk, zoom[1] - tpk)
    ax_e.set_xlabel(r"$t - t_{\rm peak}$ (d)")
    ax_e.set_ylabel(r"$\Delta m = m - m_{\rm base}$ (mag)")
    ps.panel_label(ax_e, "(e)", loc="upper left")
    omitted_note(ax_e, n_out)

    # ---- (f) the parallax signal itself ------------------------------------------------
    # Data minus the no-parallax model, with (model - model_noplx) drawn over it. In the
    # geocentric frame this is the observer's ACCELERATION: zero at the peak by construction.
    # A point is drawn only if its error is below INFO_FRACTION of the signal it is compared
    # with -- otherwise it cannot tell the two curves apart, and says so by its absence.
    sig = 0.0
    for i, f in enumerate(FILTERS):
        fr = 1 if i == ROMAN_FILT else 0
        mf = frames[fr]
        sel = (mf["t"] >= win[0]) & (mf["t"] <= win[1])
        d = mf[f"mag_{f}"][sel] - m_noplx(fr, i, mf["t"][sel])
        sig = max(sig, float(np.max(np.abs(d))) if d.size else 0.0)
        ax_f.plot(mf["t"][sel] - tpk, d, color=cols[f], lw=0.9, zorder=4)
    n_out = 0
    for i, f in enumerate(FILTERS[:6]):
        s = (rub["filt"] == i) & (rub["t"] >= win[0]) & (rub["t"] <= win[1])
        if s.any():
            e = subset(rub, s)
            ok = e["err_mag"] < INFO_FRACTION * sig
            n_out += np.count_nonzero(~ok)
            if ok.any():
                ax_f.errorbar(e["t"][ok] - tpk, (e["mag_obs"] - m_noplx(0, i, e["t"], e))[ok],
                              yerr=e["err_mag"][ok], fmt=MARKERS[f], ms=2.0, mew=0,
                              color=cols[f], ecolor=cols[f], elinewidth=0.35, alpha=0.7,
                              zorder=3)
    s = (rom["t"] >= win[0]) & (rom["t"] <= win[1])
    if s.any():
        e = subset(rom, s)
        tb, yb, sb, nb = bin_ivar(e["t"], e["mag_obs"] - m_noplx(1, ROMAN_FILT, e["t"], e),
                                  e["err_mag"], BIN_ROMAN_PHOT)
        ok = sb < INFO_FRACTION * sig
        n_out += int(nb[~ok].sum())
        ax_f.errorbar(tb[ok] - tpk, yb[ok], yerr=sb[ok], fmt="o", ms=1.8, mew=0,
                      color=cols["F146"], ecolor=cols["F146"], elinewidth=0.5, zorder=5)
    shade(ax_f, lambda t: t - tpk)
    ax_f.axhline(0.0, color=ps.MUTED, lw=0.6, zorder=1)
    lim = max(1.6 * sig, 1e-5)
    ax_f.set_ylim(lim, -lim)
    ax_f.set_xlim(win[0] - tpk, win[1] - tpk)
    ax_f.set_xlabel(r"$t - t_{\rm peak}$ (d)")
    ax_f.set_ylabel(r"$m - m_{\rm no\ parallax}$ (mag)")
    ps.panel_label(ax_f, "(f)", loc="upper left")
    omitted_note(ax_f, n_out)

    # ---- astrometric binning, chosen per event -----------------------------------------
    selm = (m1["t"] >= win[0]) & (m1["t"] <= win[1])
    peak = float(np.max(np.hypot(m1["def1c"], m1["def2c"])[selm])) if selm.any() else 0.0
    hi_days = sum(b - a for a, b, h in seasons if h)
    n_hi = sum(np.count_nonzero((rom["t"] >= a) & (rom["t"] <= b)) for a, b, h in seasons if h)
    w_rom = astr_bin_width(rom["err_ast"], n_hi / hi_days if hi_days > 0 else 0.0, peak, tE)
    covered = sum(max(0.0, min(b, win[1]) - max(a, win[0])) for a, b, _ in seasons)
    w_rom = max(w_rom, covered / ASTR_MAX_BINS)
    rate_rub = rub["t"].size / max(np.ptp(rub["t"]), 1.0) if rub["t"].size > 1 else 0.0
    w_rub = astr_bin_width(rub["err_ast"], rate_rub, peak, tE)

    def astro_bins(tab, width, key1="def1c", key2="def2c", tlim=win):
        s = (tab["t"] >= tlim[0]) & (tab["t"] <= tlim[1])
        if not s.any():
            return None
        e = subset(tab, s)
        tb, xb, sx, nb = bin_ivar(e["t"], e[key1], e["err_ast"], width, ASTR_BIN_MIN_N)
        _, yb, sy, _ = bin_ivar(e["t"], e[key2], e["err_ast"], width, ASTR_BIN_MIN_N)
        return tb, xb, yb, sx, sy, nb

    # ---- (g) |centroid shift| vs time --------------------------------------------------
    gx, gy = shift_noplx(1, m1["t"], m1)
    ax_g.plot(m1["t"][selm] - tpk, np.hypot(m1["def1c"], m1["def2c"])[selm],
              color=ps.SURVEY["roman"], lw=1.0, label="model, with parallax")
    ax_g.plot(m1["t"][selm] - tpk, np.hypot(gx, gy)[selm], color=ps.SURVEY["roman"], lw=0.9,
              ls=(0, (3, 2)), label="without")
    n_out = 0
    for tab, width, colour, mk, name in ((rom, w_rom, ps.SURVEY["roman"], "o", "Roman"),
                                         (rub, w_rub, ps.SURVEY["rubin"], "s", "Rubin")):
        b = astro_bins(tab, width)
        if b is None or not b[0].size:
            continue
        tb, xb, yb, sx, _, nb = b
        ok = sx < INFO_FRACTION * peak
        n_out += np.count_nonzero(~ok)
        if ok.any():
            ax_g.errorbar(tb[ok] - tpk, np.hypot(xb, yb)[ok], yerr=sx[ok], fmt=mk, ms=2.0,
                          mew=0, color=colour, ecolor=colour, elinewidth=0.5, alpha=0.85,
                          label=f"{name}, {width:.3g}-d bins", zorder=4)
    top = max(1.7 * peak, 0.05)
    if 1.1 < top:
        ax_g.axhline(1.1, color=ps.MUTED, lw=0.7, ls=":", zorder=1)
        ax_g.text(0.02, 1.1, " 1.1 mas: Roman per-exposure floor", ha="left", va="bottom",
                  fontsize=6, color=ps.MUTED, transform=ax_g.get_yaxis_transform())
    ax_g.set_ylim(0.0, top)
    shade(ax_g, lambda t: t - tpk)
    ax_g.set_xlim(win[0] - tpk, win[1] - tpk)
    ax_g.set_xlabel(r"$t - t_{\rm peak}$ (d)")
    ax_g.set_ylabel(r"$|\delta\theta_c|$ (mas)")
    ps.panel_label(ax_g, "(g)", loc="upper left")
    ps.legend(ax_g, loc="upper right", fontsize=6)
    omitted_note(ax_g, n_out, "bins")

    # ---- (h) the astrometric ellipse ---------------------------------------------------
    ax_h.plot(m1["def1c"][selm], m1["def2c"][selm], color=ps.SURVEY["roman"], lw=1.0,
              label="with parallax")
    ax_h.plot(gx[selm], gy[selm], color=ps.SURVEY["roman"], lw=0.9, ls=(0, (3, 2)),
              label="without")
    n_out = 0
    b = astro_bins(rom, w_rom)
    if b is not None and b[0].size:
        _, xb, yb, sx, sy, _ = b
        ok = np.maximum(sx, sy) < INFO_FRACTION * peak
        n_out = np.count_nonzero(~ok)
        if ok.any():
            ax_h.errorbar(xb[ok], yb[ok], xerr=sx[ok], yerr=sy[ok], fmt="o", ms=1.8, mew=0,
                          color=ps.SURVEY["roman"], ecolor=ps.SURVEY["roman"],
                          elinewidth=0.4, alpha=0.6, label=f"Roman, {w_rom:.3g}-d bins",
                          zorder=4)
    ax_h.plot([0], [0], marker="+", color=ps.INK, ms=6, mew=0.8, ls="")
    xs = np.r_[m1["def1c"][selm], gx[selm], 0.0]
    ys = np.r_[m1["def2c"][selm], gy[selm], 0.0]
    half = 0.62 * max(np.ptp(xs), np.ptp(ys), 1e-3)
    cx, cy = 0.5 * (xs.max() + xs.min()), 0.5 * (ys.max() + ys.min())
    ax_h.set_xlim(cx - half, cx + half)
    ax_h.set_ylim(cy - half, cy + half)
    ax_h.set_aspect("equal", adjustable="box")
    ax_h.set_xlabel(r"$\delta\theta_{c,1}$ (mas)")
    ax_h.set_ylabel(r"$\delta\theta_{c,2}$ (mas)")
    ps.panel_label(ax_h, "(h)", loc="upper left")
    ps.legend(ax_h, loc="upper right", fontsize=6)
    omitted_note(ax_h, n_out, "bins")

    # ---- (b) the sky-plane trajectories: Figure 2's right-hand panel -------------------
    # Absolute positions over the whole model span, Earth frame. The lens WITHOUT parallax
    # is its straight proper-motion line through the centre of its parallax ellipse. The
    # deflection is drawn about the origin, as in Figure 2. At this project's thetaE (0.1 to a
    # few mas, against 70-80 mas for Figure 2's LMC events) the deflected and undeflected
    # source coincide on the decade scale, so panel (d) zooms on the closest approach: the
    # source and lens are within a few thetaE of the origin at the peak, so one box holds
    # both of them and the deflection loop.
    t = m0["t"]
    inside = (t >= 0.0) & (t <= TOBS)
    su = (m0["pos1b"], m0["pos2b"])
    sd = (m0["pos1c"], m0["pos2c"])
    ln = (m0["lens1"], m0["lens2"])
    ln0 = (straight_track(t, m0["lens1"], p["mul1"], t0),
           straight_track(t, m0["lens2"], p["mul2"], t0))
    rel = (ln[0] - su[0], ln[1] - su[1])
    dfl = (m0["def1c"], m0["def2c"])
    dfl0 = shift_noplx(0, t, m0)

    def track(ax, xy, colour, ls="-", lw=1.0, z=3):
        x, y = xy
        ax.plot(np.where(inside, x, np.nan), np.where(inside, y, np.nan), color=colour,
                ls=ls, lw=lw, zorder=z)
        ax.plot(np.where(~inside, x, np.nan), np.where(~inside, y, np.nan), color=colour,
                ls=":", lw=lw, zorder=z)

    def all_tracks(ax, inset=False):
        track(ax, rel, TRACK["rel"], lw=0.9)
        track(ax, ln0, TRACK["lens0"], ls=(0, (3, 2)), lw=0.9)
        track(ax, ln, TRACK["lens"])
        track(ax, su, TRACK["src_u"])
        track(ax, sd, TRACK["src_d"])
        track(ax, dfl, TRACK["defl"], lw=1.0, z=5)
        track(ax, dfl0, TRACK["defl"], ls=(0, (3, 2)), lw=0.8, z=5)
        if inset:
            # the deflected source WITHOUT microlensing parallax: undeflected + no-plx shift
            track(ax, (su[0] + dfl0[0], su[1] + dfl0[1]), TRACK["src_d"], ls=(0, (3, 2)),
                  lw=0.8)

    all_tracks(ax_b)
    X = np.r_[su[0], sd[0], ln[0], ln0[0], rel[0], dfl[0]]
    Y = np.r_[su[1], sd[1], ln[1], ln0[1], rel[1], dfl[1]]
    span = 1.08 * max(np.ptp(X), np.ptp(Y))
    xc, yc = 0.5 * (X.max() + X.min()), 0.5 * (Y.max() + Y.min())
    ax_b.set_xlim(xc - span / 2, xc + span / 2)
    ax_b.set_ylim(yc - span / 2, yc + span / 2)
    ax_b.set_aspect("equal", adjustable="box")
    ax_b.set_xlabel("$x$ (mas)")
    ax_b.set_ylabel("$y$ (mas)")
    ps.panel_label(ax_b, "(b)", loc="upper left")

    # panel (d): the closest approach, the deflection loop, and Roman's binned positions
    # The zoom is set by the ANGULAR scale of the lensing, not by a time window: +-1.5 tE is
    # the whole decade for a slow event. Half-width max(1.5 thetaE, 3 x the peak shift),
    # centred between the deflected source and the lens at the observed peak and the origin
    # the deflection is drawn about.
    ipk = int(np.argmin(np.abs(t - tpk)))
    ih = max(1.5 * tetE, 3.0 * float(np.max(np.hypot(*dfl))), 1e-3)
    ixc = (sd[0][ipk] + ln[0][ipk] + 0.0) / 3.0
    iyc = (sd[1][ipk] + ln[1][ipk] + 0.0) / 3.0
    ins = ax_d
    all_tracks(ins, inset=True)
    n_out = 0
    b = astro_bins(rom, w_rom, "pos1c", "pos2c")
    if b is not None and b[0].size:
        _, xb, yb, sx, sy, _ = b
        ok = np.maximum(sx, sy) < INFO_FRACTION * peak
        n_out = np.count_nonzero(~ok)
        if ok.any():
            ins.errorbar(xb[ok], yb[ok], xerr=sx[ok], yerr=sy[ok], fmt="o", ms=1.6, mew=0,
                         color=ps.SURVEY["roman"], ecolor=ps.SURVEY["roman"], elinewidth=0.4,
                         alpha=0.75, zorder=6)
    ins.set_xlim(ixc - ih, ixc + ih)
    ins.set_ylim(iyc - ih, iyc + ih)
    ins.set_aspect("equal", adjustable="box")
    ins.set_xlabel("$x$ (mas)")
    ins.set_ylabel("$y$ (mas)")
    ps.panel_label(ins, "(d)", loc="upper left")
    ins.text(0.98, 0.97, "closest approach", transform=ins.transAxes,
             ha="right", va="top", fontsize=6, color=ps.MUTED)
    # Roman's binned positions outside the box are simply off-panel, not omitted
    omitted_note(ins, n_out, "Roman bins")
    # the zoomed region, marked on the decade panel
    ax_b.add_patch(plt.Rectangle((ixc - ih, iyc - ih), 2 * ih, 2 * ih, fill=False,
                                 ec=ps.MUTED, lw=0.6, zorder=7))

    # ---- (i) sky tracks with the linear proper motion removed --------------------------
    # pos - mu (t - t0) leaves each body's own annual-parallax loop (1/D in mas). These are the
    # SOURCE and LENS parallaxes, pi_s and pi_l -- not pi_E. Without parallax, each collapses
    # to the centre of its loop, marked x.
    dt = t - t0
    lx, ly = m0["lens1"] - p["mul1"] * dt, m0["lens2"] - p["mul2"] * dt
    ux, uy = m0["pos1b"] - p["mus1"] * dt, m0["pos2b"] - p["mus2"] * dt
    cxd, cyd = m0["pos1c"] - p["mus1"] * dt, m0["pos2c"] - p["mus2"] * dt
    ax_i.plot(lx, ly, color=TRACK["lens"], lw=0.9,
              label=rf"lens ($\pi_l = {1.0 / p['Dl']:.3g}$ mas)")
    ax_i.plot(ux, uy, color=TRACK["src_u"], lw=0.9, ls=(0, (3, 2)),
              label=rf"source, undeflected ($\pi_s = {1.0 / p['Ds']:.3g}$ mas)")
    ax_i.plot(cxd, cyd, color=TRACK["src_d"], lw=1.0, label="source, deflected")
    ax_i.plot(ln0[0][0] - p["mul1"] * dt[0], ln0[1][0] - p["mul2"] * dt[0], marker="x",
              color=TRACK["lens0"], ms=5, mew=1.0, ls="", label="no parallax (loop centres)")
    c_s = (straight_track(t, m0["pos1b"], p["mus1"], t0)[0] - p["mus1"] * dt[0],
           straight_track(t, m0["pos2b"], p["mus2"], t0)[0] - p["mus2"] * dt[0])
    ax_i.plot(*c_s, marker="x", color=TRACK["lens0"], ms=5, mew=1.0, ls="")
    allx = np.r_[lx, ux, cxd]
    ally = np.r_[ly, uy, cyd]
    x0, x1 = float(allx.min()), float(allx.max())
    y0, y1 = float(ally.min()), float(ally.max())
    y0 -= 0.55 * (y1 - y0)                     # the legend's space, below the loops
    span = 1.06 * max(x1 - x0, y1 - y0)
    xc, yc = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    ax_i.set_xlim(xc - span / 2, xc + span / 2)
    ax_i.set_ylim(yc - span / 2, yc + span / 2)
    ax_i.set_aspect("equal", adjustable="box")
    ax_i.set_xlabel(r"$x - \mu_x (t - t_0)$ (mas)")
    ax_i.set_ylabel(r"$y - \mu_y (t - t_0)$ (mas)")
    ps.panel_label(ax_i, "(i)", loc="upper left")
    ps.legend(ax_i, loc="lower right", fontsize=6)

    for ax in (ax_c, ax_e, ax_f, ax_g):
        ax.axvline(0.0, color=ps.GRID, lw=0.8, zorder=0)

    cls, sid = os.path.basename(stem).rsplit("_", 1)
    commit = prov.get("git_commit", "?") if prov else "?"
    ps.stamp(fig, (
        f"sample {cls} {sid} · population {p['population']} · commit {commit} · no-parallax "
        f"model: {noplx} · shading: Roman seasons at this sightline, dark = high cadence, "
        f"light = low\n"
        f"photometry is the simulated draw (Roman raw, and in {BIN_ROMAN_PHOT:g}-d bins) · "
        f"astrometric points: the model at observed epochs, inverse-variance binned, error bars "
        f"from the instrument model, exposures assumed independent\n"
        f"a point or bin whose 1$\\sigma$ exceeds the signal in its panel is not drawn; each "
        f"panel states how many were left out"))
    os.makedirs(out_dir, exist_ok=True)
    suffix = "" if noplx == "geocentric" else "_simgauge"
    paths = ps.save_figure(fig, os.path.join(out_dir, os.path.basename(stem) + suffix))
    return paths, checks


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("samples", help="directory holding <class>_<id>_*.dat from --dump-samples")
    ap.add_argument("-o", "--out", required=True, help="output directory for the figures")
    ap.add_argument("--prov", help="run_provenance.txt of the run that wrote the samples")
    ap.add_argument("--only", help="draw only this <class>_<id>")
    ap.add_argument("--noplx", choices=("geocentric", "simulation"), default="geocentric",
                    help="which no-parallax model to draw; see the module docstring. Only "
                         "'geocentric' is fit to publish")
    args = ap.parse_args()

    ps.use_paper_style()
    prov = R.load_provenance(args.prov) if args.prov else None
    stems = sorted(f[:-len("_params.dat")]
                   for f in glob.glob(os.path.join(args.samples, "*_params.dat")))
    if args.only:
        stems = [s for s in stems if os.path.basename(s) == args.only]
    if not stems:
        sys.exit(f"no *_params.dat under {args.samples}")

    for stem in stems:
        paths, checks = draw_event(stem, args.out, prov, args.noplx)
        print(f"{os.path.basename(stem)} -> {paths[0]}")
        for c in checks:
            print(f"    {c}")


if __name__ == "__main__":
    main()
