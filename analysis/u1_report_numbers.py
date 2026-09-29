#!/usr/bin/env python3
"""Step U1: every pooled number the overview report quotes, with a Monte Carlo error bar.

WHY THIS EXISTS. The analysis scripts print their headline numbers as bare values -- a weighted
median of 0.973, a share of 72.9% -- with the Kish N_eff beside them as the only statement of
precision. N_eff says how much the weight costs; it does not say how uncertain a MEDIAN of a
skewed ratio is, or a fraction of 0.2%. A report that wants to be read as science needs the
error on the number itself. This script recomputes each quoted number from the event tables and
attaches a bootstrap standard deviation, and it prints the value the original script logged so
the two can be compared: a recomputed central value that does not reproduce the log means the
selection here is not the selection there, and the error bar would belong to a different number.

It also computes the numbers the report did not yet have: the image-resolution probabilities of
all three populations (after Sajadian & Makler, arXiv:2608.16448) with their yields, and what
Roman does for the events RUBIN detects inside the footprint.

THE ERROR MODEL. A Poisson bootstrap over simulated events: each replicate multiplies every
event's weight by an independent Poisson(1) count and recomputes the statistic. That is the
standard large-sample stand-in for resampling with replacement, and it has one property that
matters here -- the counts can be drawn once per population and shared by every statistic over
that population, so the replicates are consistent across numbers that are later compared. The
error covers the Monte Carlo sampling only. It does not cover modelling assumptions, which the
report treats separately and which are larger.

Yields carry y1's own error, sqrt(sum y_i^2) over the selected draws (Poisson on the draws),
so a yield here must reproduce y1's value AND its error.

INPUTS. Detection-only tables, `test<tag>_detJ.dat` (every row with detJ == 1, header kept):
    awk 'NR==1 || /^#/ {print; next} $40==1' test<tag>.dat > test<tag>_detJ.dat
Exact for everything here, because a detected row's weight depends only on that row and its
sightline's nsim (from the map file), never on the other rows. The one quantity that needs the
undetected draws -- each lens component's mean mass <M>, which a yield divides by -- is read
from y1's own output (`--mean-mass`), where it was computed over every draw.
The intrinsic-timescale check needs every draw; it reads a column extract (`--draws`):
    awk '/^#/ {if (!h) {print "#tE Vt Ml Ds lon lat w_area detJ Ai_r ndw_R"; h=1}; next}
         {print $3,$7,$9,$12,$60,$61,$89,$40,$35,$37}' test5.dat > test5_w1cols.dat

    .roman/bin/python analysis/u1_report_numbers.py \\
        --run bulge=runs/prod_bulge_20260924 --run bh=runs/prod_bh_20260924 \\
        --run ns=runs/prod_ns_20260924 \\
        --mean-mass figures/yield_20260925 \\
        --draws bulge=runs/prod_bulge_20260924/test5_w1cols.dat \\
        -o figures/u1_20260929
"""

import argparse
import ast
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import galaxy_model as G      # noqa: E402
import romanlib as R          # noqa: E402

N_BOOT = 400
SEED = 20260929

# The paper's criterion (Sajadian & Makler, their sec. 3): resolvable when at least three
# recorded data points have both images detectable AND separated by more than the bar.
RESOLVE_MIN_EPOCHS = 3
ROMAN_AST_FLOOR = 1.1          # mas, Bulge.h ROMAN_AST_FLOOR
U_AST_PEAK = np.sqrt(2.0)
BARS = ("nres5", "nres20", "nresPSF")

COLS = ["tE", "piE", "tetE", "u0", "Ml", "Dl", "Ds", "Vt", "lon", "lat", "w_area", "struc",
        "detL", "detR", "detJ",
        "okA_J", "okA_L", "okA_R", "okB_J", "okB_L", "okB_R",
        "sigtE_J", "sigtE_L", "sigtE_R", "sigpiE_J", "sigpiE_L", "sigpiE_R",
        "sigtetE_J", "sigtetE_L", "sigtetE_R", "relMl_J", "relMl_L", "relMl_R",
        "ndw_L", "ndw_R", "t0zone", "dt_edge", "Ai_r", "magb_F146", "blend_F146", "nepR_pk",
        "nres5_L", "nres20_L", "nresPSF_L", "dsep_max_L",
        "nres5_R", "nres20_R", "nresPSF_R", "dsep_max_R"]

PAIR_COLS = ("lon lat tE u0 piE tetE okA_sat okA_nosat sigpiE_sat sigpiE_nosat "
             "nepR_pk w_area Ml Dl Ds Vt").split()

SURV = {"J": "joint", "L": "rubin", "R": "roman"}


# ---------------------------------------------------------------------------------------------
# The bootstrap
# ---------------------------------------------------------------------------------------------
class Boot:
    """Poisson(1) replicate counts for one sample of n events, shared by every statistic on it."""

    def __init__(self, n, B=N_BOOT, seed=SEED):
        rng = np.random.default_rng(seed)
        # uint8 is ample: P(Poisson(1) > 255) is zero to machine precision.
        self.C = rng.poisson(1.0, size=(n, B)).astype(np.uint8)
        self.B = B

    def _block(self, idx, lo, hi):
        return self.C[idx, lo:hi].astype(np.float64)

    def frac(self, num, den, w):
        """Weighted fraction sum(w[num]) / sum(w[den]); num must be a subset of den."""
        num, den = np.asarray(num, bool), np.asarray(den, bool)
        w = np.asarray(w, float)
        idx = np.flatnonzero(den)
        if idx.size == 0 or w[idx].sum() <= 0:
            return np.nan, np.nan, 0, 0.0
        v = w[num & den].sum() / w[idx].sum()
        wn = (w * num)[idx]
        wd = w[idx]
        reps = []
        for lo in range(0, self.B, 100):
            Cb = self._block(idx, lo, min(lo + 100, self.B))
            reps.append((wn @ Cb) / (wd @ Cb))
        reps = np.concatenate(reps)
        return float(v), float(np.nanstd(reps)), int(idx.size), R.kish_neff(wd)

    def quantile(self, values, sel, w, q=0.5):
        """Weighted quantile of `values` over `sel` (finite, positive-weight rows only)."""
        values = np.asarray(values, float)
        w = np.asarray(w, float)
        sel = np.asarray(sel, bool) & np.isfinite(values) & (w > 0)
        idx = np.flatnonzero(sel)
        if idx.size < 20:
            return np.nan, np.nan, int(idx.size), 0.0
        o = np.argsort(values[idx], kind="stable")
        idx = idx[o]
        v, ws = values[idx], w[idx]
        c = np.cumsum(ws)
        central = float(v[np.searchsorted(c, q * c[-1])])
        reps = []
        for lo in range(0, self.B, 50):
            cw = np.cumsum(ws[:, None] * self._block(idx, lo, min(lo + 50, self.B)), axis=0)
            k = [np.searchsorted(cw[:, j], q * cw[-1, j]) for j in range(cw.shape[1])]
            reps.append(v[np.minimum(k, len(v) - 1)])
        reps = np.concatenate(reps)
        return central, float(np.std(reps)), int(idx.size), R.kish_neff(ws)

    def ratio_sums(self, a, b, w):
        """sum(w[a]) / sum(w[b]) for two selections that need not be nested."""
        a, b = np.asarray(a, bool), np.asarray(b, bool)
        w = np.asarray(w, float)
        idx = np.flatnonzero(a | b)
        wa, wb = (w * a)[idx], (w * b)[idx]
        v = wa.sum() / wb.sum()
        reps = []
        for lo in range(0, self.B, 100):
            Cb = self._block(idx, lo, min(lo + 100, self.B))
            reps.append((wa @ Cb) / (wb @ Cb))
        return float(v), float(np.std(np.concatenate(reps)))


# ---------------------------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------------------------
def mean_mass_from_y1(path):
    """{struc: <M>} from y1's own markdown, where it was computed over every draw."""
    with open(path) as fh:
        for line in fh:
            m = re.search(r"mean lens mass <M> \[Msun\], per component: (\{.*\})", line)
            if m:
                return {int(k): float(v) for k, v in ast.literal_eval(m.group(1)).items()}
    raise ValueError(f"{path}: no 'mean lens mass' line")


class Run:
    def __init__(self, name, directory, mean_mass_dir):
        self.name, self.dir = name, directory
        prov = R.load_provenance(os.path.join(directory, "files/MONTLMC/files/run_provenance.txt"))
        self.prov = prov
        self.population = R.population(prov)
        tag = prov.get("population_tag", "5")
        full = os.path.join(directory, f"test{tag}.dat")
        det = os.path.join(directory, f"test{tag}_detJ.dat")
        self.table = det if os.path.exists(det) else full
        self.map = os.path.join(directory, f"files/MONTLMC/files/MapLMC{tag}.dat")
        self.logs = [os.path.join(directory, "run.log")]
        self.tobs = float(prov.get("Tobs_days", 3652.43))
        print(f"[{name}] reading {self.table}", flush=True)
        keep = (lambda c: c["detJ"] == 1) if self.table == full else None
        df = R.load_events(self.table, usecols=COLS, chunksize=500_000, keep=keep)
        if not (df["detJ"] == 1).all():
            raise ValueError(f"{self.table}: expected detections only")
        w, self.wlabel = R.attach_weight(df, self.map, self.logs)
        df["W"] = w.to_numpy()
        # Absolute yield per row, F = 1, per object: y1's formula with <M> per lens component.
        mm = mean_mass_from_y1(os.path.join(mean_mass_dir, name, "y1_yields.md"))
        Mbar = df["struc"].map(mm).to_numpy(float)
        if np.isnan(Mbar).any():
            raise ValueError(f"[{name}] a lens component has no <M> in y1's output")
        T = (self.tobs - 2.0 * R.T0_MARGIN_DAYS) * 86400.0
        df["y"] = T * R.RATE_UNIT * df["W"].to_numpy() / Mbar
        # Footprint = sightlines Roman observes, exactly as y1 defines it (any row with Roman
        # epochs marks the whole sightline).
        codes, _ = R.sightline_index(df)
        df["foot"] = np.isin(codes, np.unique(codes[df["ndw_R"].to_numpy() > 0]))
        self.df = df
        self.boot = Boot(len(df))
        print(f"[{name}] {len(df):,} detections, {int(df.foot.sum()):,} in the footprint, "
              f"{self.wlabel}", flush=True)


# ---------------------------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------------------------
ROWS = []


def rec(pop, section, quantity, value, err, n=None, neff=None, logged=None, scale=1.0,
        note=""):
    ROWS.append(dict(population=pop, section=section, quantity=quantity,
                     value=value * scale if value is not None else np.nan,
                     err=err * scale if err is not None else np.nan,
                     n=n, n_eff=neff, logged=logged, note=note))


def sig_over(df, p, s):
    """sigma(p)/p from survey s, NaN where not measured (romanlib gating)."""
    sv = R.sigma(df, p, SURV[s])
    if p == "Ml":
        return sv.to_numpy(float)                 # relMl is already relative
    return (sv / df[p]).to_numpy(float)


def max_centroid_shift(df):
    """Largest centroid shift [mas]: theta_E u/(u^2+2) at u = max(u0, sqrt 2) (p6's helper)."""
    te = df["tetE"].to_numpy(float)
    u0 = df["u0"].to_numpy(float)
    u = np.where(u0 <= U_AST_PEAK, U_AST_PEAK, u0)
    return te * u / (u * u + 2.0)


# ---------------------------------------------------------------------------------------------
# The numbers
# ---------------------------------------------------------------------------------------------
def yields_block(run):
    """Reproduce y1's footprint and whole-scan yields for the rows the report quotes."""
    d, pop = run.df, run.name
    y = d["y"].to_numpy()
    foot = d["foot"].to_numpy()
    L, Rr = (d.detL == 1).to_numpy(), (d.detR == 1).to_numpy()
    mR = (d.relMl_R > 0).to_numpy() & (d.relMl_R < 0.1).to_numpy()
    mJ = (d.relMl_J > 0).to_numpy() & (d.relMl_J < 0.1).to_numpy()
    sels = {
        "joint": np.ones(len(d), bool), "Rubin detects": L, "Roman detects": Rr,
        "both detect": L & Rr, "Rubin only": L & ~Rr, "Roman only": Rr & ~L,
        "Roman, sigma(M)/M<10%": Rr & mR, "joint, sigma(M)/M<10%": mJ,
        "Roman, sigma(M)/M<1%": Rr & (d.relMl_R > 0).to_numpy() & (d.relMl_R < .01).to_numpy(),
    }
    for scope, m in (("whole scan", np.ones(len(d), bool)), ("footprint", foot)):
        for name, s in sels.items():
            s = s & m
            rec(pop, f"yield N_1, {scope}", name, y[s].sum(), np.sqrt((y[s] ** 2).sum()),
                n=int(s.sum()))
    # What combining the surveys adds to the 10% masses, as a ratio with a correlated error.
    v, e = run.boot.ratio_sums(mJ & foot, Rr & mR & foot, y)
    rec(pop, "yield N_1, footprint", "joint/Roman, sigma(M)/M<10%", v, e)
    # Share of all detections that fall outside the footprint.
    f = run.boot.frac(~foot, np.ones(len(d), bool), y)
    rec(pop, "yield N_1, whole scan", "share of detections outside the footprint", f[0], f[1],
        n=f[2])
    # Share of Rubin's detections that fall inside the footprint.
    f = run.boot.frac(L & foot, L, y)
    rec(pop, "yield N_1, whole scan", "share of Rubin detections inside the footprint", f[0],
        f[1], n=f[2])
    # Of the events Rubin detects in the footprint, the share Roman detects too.
    f = run.boot.frac(L & Rr & foot, L & foot, y)
    rec(pop, "yield N_1, footprint", "Roman also detects, of Rubin's detections", f[0], f[1],
        n=f[2], neff=f[3])


def detection_shares(run, logged=None):
    d, b = run.df, run.boot
    W = d["W"].to_numpy()
    L, Rr = (d.detL == 1).to_numpy(), (d.detR == 1).to_numpy()
    allm = np.ones(len(d), bool)
    for lab, m in (("Rubin only", L & ~Rr), ("Roman only", Rr & ~L), ("both", L & Rr)):
        f = b.frac(m, allm, W)
        rec(run.name, "detection share (weighted, all detections)", lab, f[0], f[1], n=f[2],
            neff=f[3], scale=100, logged=(logged or {}).get(lab))


def joint_over_single(run, logged=None):
    """p6's per-event medians on events both surveys' photometric matrices constrain."""
    d, b = run.df, run.boot
    W = d["W"].to_numpy()
    both = ((d.okA_L == 1) & (d.okA_R == 1)).to_numpy()
    rec(run.name, "joint/single, events both surveys constrain", "N (okA_L & okA_R)",
        int(both.sum()), None, logged=(logged or {}).get("N"))
    for p in ("tE", "piE", "tetE"):
        num = d[f"sig{p}_J"].to_numpy(float)
        for other, lab in (("L", "Rubin"), ("R", "Roman")):
            den = d[f"sig{p}_{other}"].to_numpy(float)
            g = both & (num > 0) & (den > 0)
            with np.errstate(divide="ignore", invalid="ignore"):
                r = num / den
            q = b.quantile(r, g, W)
            rec(run.name, "joint/single, events both surveys constrain",
                f"median sigma_J/sigma_{lab} ({p})", q[0], q[1], n=q[2], neff=q[3],
                logged=(logged or {}).get((p, lab)))


def astrometric_shift(run, scope_mask, scope, logged=None):
    d, b = run.df, run.boot
    W = d["W"].to_numpy()
    s = max_centroid_shift(d)
    ok = scope_mask & np.isfinite(s) & (s > 0)
    q = b.quantile(s, ok, W)
    rec(run.name, f"centroid shift ({scope})", "median max shift [mas]", q[0], q[1], n=q[2],
        neff=q[3], logged=(logged or {}).get("median"))
    f = b.frac(ok & (s > ROMAN_AST_FLOOR), ok, W)
    rec(run.name, f"centroid shift ({scope})", "above Roman's 1.1 mas floor [%]", f[0], f[1],
        n=f[2], neff=f[3], scale=100, logged=(logged or {}).get("floor"))


def precision_10pct(run, logged=None):
    """p7's 'fraction of measured events better than 10%', per survey: each over its own
    measured set (ok flag for the matrix the parameter comes from, and a positive sigma)."""
    d, b = run.df, run.boot
    W = d["W"].to_numpy()
    for p in ("tE", "piE", "tetE", "Ml"):
        for s in ("J", "L", "R"):
            gate = d[f"okB_{s}"] if p in ("tetE", "Ml") else d[f"okA_{s}"]
            x = sig_over(d, p, s)
            m = (gate == 1).to_numpy() & np.isfinite(x) & (x > 0)
            f = b.frac(m & (x < 0.1), m, W)
            rec(run.name, "measured to 10% (own measured set, p7)", f"{p} {SURV[s]}", f[0],
                f[1], n=f[2], neff=f[3], scale=100, logged=(logged or {}).get((p, s)))


def resolution(run, logged=None):
    """P(images resolvable) per survey and bar over that survey's detections, and the yield."""
    d, b = run.df, run.boot
    W, y = d["W"].to_numpy(), d["y"].to_numpy()
    foot = d["foot"].to_numpy()
    for s, lab in (("L", "Rubin"), ("R", "Roman")):
        det = (d[f"det{s}"] == 1).to_numpy()
        for bar in BARS:
            res = d[f"{bar}_{s}"].to_numpy() >= RESOLVE_MIN_EPOCHS
            f = b.frac(det & res, det, W)
            rec(run.name, "image resolution: P over the survey's detections", f"{lab} {bar}",
                f[0], f[1], n=f[2], neff=f[3], scale=100,
                logged=(logged or {}).get((lab, bar)))
            for scope, m in (("whole scan", np.ones(len(d), bool)), ("footprint", foot)):
                sel = det & res & m
                rec(run.name, f"image resolution: yield N_1, {scope}", f"{lab} {bar}",
                    y[sel].sum(), np.sqrt((y[sel] ** 2).sum()), n=int(sel.sum()))
        # Rubin in the footprint only, for a like-for-like comparison with Roman.
        if s == "L":
            for bar in BARS:
                res = d[f"{bar}_L"].to_numpy() >= RESOLVE_MIN_EPOCHS
                f = b.frac(det & foot & res, det & foot, W)
                rec(run.name, "image resolution: P, Rubin detections in the footprint",
                    f"Rubin {bar}", f[0], f[1], n=f[2], neff=f[3], scale=100)
    # What kind of lens gets resolved: Roman at D = 5, the bar that is not empty for every
    # population.
    det = (d.detR == 1).to_numpy()
    res = det & (d["nres5_R"].to_numpy() >= RESOLVE_MIN_EPOCHS)
    for col, lab in (("tetE", "theta_E [mas]"), ("Ml", "M_L [Msun]"), ("Dl", "D_L [kpc]"),
                     ("tE", "t_E [d]"), ("dsep_max_R", "largest resolved separation [mas]")):
        v = d[col].to_numpy(float)
        for nm, m in (("resolvable (Roman, D=5)", res), ("all Roman detections", det)):
            ok = m & (v > 0)
            q = b.quantile(v, ok, W)
            rec(run.name, "image resolution: who is resolved", f"median {lab}, {nm}", q[0],
                q[1], n=q[2], neff=q[3])
    # theta_E threshold: the smallest theta_E at which half of Roman's detections resolve.
    te = d["tetE"].to_numpy(float)
    for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 5), (5, 10), (10, 30), (30, 1e9)):
        m = det & (te >= lo) & (te < hi)
        if m.sum() >= 20:
            f = b.frac(m & res, m, W)
            rec(run.name, "image resolution: P(Roman, D=5) by theta_E",
                f"theta_E in [{lo},{hi}) mas", f[0], f[1], n=f[2], neff=f[3], scale=100)


def peak_coverage(run):
    """Y5's 'no Roman epoch within +-2 tE of the peak', as a share of Roman's detected yield."""
    d, b = run.df, run.boot
    det = (d.detR == 1).to_numpy()
    f = b.frac(det & (d.nepR_pk == 0).to_numpy(), det, d["y"].to_numpy())
    rec(run.name, "Roman detections by peak coverage (yield-weighted)",
        "no Roman epoch within +-2 tE of the peak [%]", f[0], f[1], n=f[2], scale=100,
        logged={"bulge": 4.4, "bh": 2.9, "ns": 3.9}.get(run.name))


def roman_helps_rubin(run):
    """Inside the footprint, over the events RUBIN detects: what adding Roman does."""
    d, b = run.df, run.boot
    W, y = d["W"].to_numpy(), d["y"].to_numpy()
    foot = d["foot"].to_numpy()
    for base, lab in (("L", "Rubin"), ("R", "Roman")):
        other = "Roman" if base == "L" else "Rubin"
        den = foot & (d[f"det{base}"] == 1).to_numpy()
        sec = f"what {other} adds, over {lab}'s footprint detections"
        rec(run.name, sec, "N_1 of the denominator", y[den].sum(), np.sqrt((y[den] ** 2).sum()),
            n=int(den.sum()))
        charS = R.characterized(d, SURV[base]).to_numpy()
        charJ = R.characterized(d, "joint").to_numpy()
        crit = {"characterised (tE>2sig, piE>2sig)": (charS, charJ)}
        for p in ("tE", "piE", "tetE", "Ml"):
            xs, xj = sig_over(d, p, base), sig_over(d, p, "J")
            crit[f"{p} to 10%"] = ((xs > 0) & (xs < 0.1), (xj > 0) & (xj < 0.1))
        for name, (ms, mj) in crit.items():
            ms, mj = np.nan_to_num(ms).astype(bool), np.nan_to_num(mj).astype(bool)
            fs, fj = b.frac(den & ms, den, W), b.frac(den & mj, den, W)
            rec(run.name, sec, f"{name}: {lab} alone [%]", fs[0], fs[1], n=fs[2], neff=fs[3],
                scale=100)
            rec(run.name, sec, f"{name}: joint [%]", fj[0], fj[1], n=fj[2], neff=fj[3],
                scale=100)
            ys, yj = den & ms, den & mj
            rec(run.name, sec, f"{name}: {lab} alone, N_1", y[ys].sum(),
                np.sqrt((y[ys] ** 2).sum()), n=int(ys.sum()))
            rec(run.name, sec, f"{name}: joint, N_1", y[yj].sum(), np.sqrt((y[yj] ** 2).sum()),
                n=int(yj.sum()))
            if y[ys].sum() > 0:
                v, e = b.ratio_sums(yj, ys, y)
                rec(run.name, sec, f"{name}: joint / {lab} alone (yield ratio)", v, e)
        for p in ("tE", "piE", "tetE", "Ml"):
            xs, xj = sig_over(d, p, base), sig_over(d, p, "J")
            with np.errstate(divide="ignore", invalid="ignore"):
                r = xj / xs
            g = den & np.isfinite(r) & (xs > 0) & (xj > 0)
            q = b.quantile(r, g, W)
            rec(run.name, sec, f"median sigma_J/sigma_{lab} ({p})", q[0], q[1], n=q[2],
                neff=q[3])


def bulge_specific(run, logged):
    """The ordinary-lens numbers of the report's Sections 4.2-4.8."""
    d, b = run.df, run.boot
    W = d["W"].to_numpy()
    pop = run.name
    # --- the F4 sample: joint detections with Roman epochs (per event, as f4 selects) ---
    f4 = (d["ndw_R"] > 0).to_numpy()
    rec(pop, "F4 footprint sample", "N", int(f4.sum()), None, logged=12461)
    rec(pop, "F4 footprint sample", "N_eff", R.kish_neff(W[f4]), None, logged=3969)
    rec(pop, "F4 footprint sample", "N (per-sightline footprint, y1's definition)",
        int(d.foot.sum()), None, logged=12461)
    for p in ("tE", "piE", "tetE", "Ml"):
        for s in ("J", "R", "L"):
            x = sig_over(d, p, s)
            f = b.frac(f4 & (x > 0) & (x < 0.1), f4, W)
            rec(pop, "F4: measured to 10%, normalised to the whole footprint sample",
                f"{p} {SURV[s]}", f[0], f[1], n=f[2], neff=f[3], scale=100,
                logged=logged["f4"].get((p, s)))
    fu = float(np.mean((sig_over(d, "Ml", "J")[f4] > 0) & (sig_over(d, "Ml", "J")[f4] < 0.1)))
    rec(pop, "F4: measured to 10%, normalised to the whole footprint sample",
        "Ml joint, UNWEIGHTED", fu, None, scale=100, logged=9.9)
    # sigma(theta_E)/theta_E medians (h5 panel d / f4 csv).
    for s, lg in (("J", 0.098), ("R", 0.099), ("L", 2.30)):
        x = sig_over(d, "tetE", s)
        q = b.quantile(x, f4 & (x > 0), W)
        rec(pop, "F4: median sigma(theta_E)/theta_E", SURV[s], q[0], q[1], n=q[2], neff=q[3],
            logged=lg)
    # --- raw characterisation counts in the footprint (F3 text) ---
    cJ, cR, cL = (R.characterized(d, s).to_numpy() for s in ("joint", "roman", "rubin"))
    for lab, m, lg in (("joint", cJ, 3948), ("Roman alone", cR, 3412), ("Rubin alone", cL, 670),
                       ("joint and neither alone", cJ & ~cR & ~cL, 358)):
        rec(pop, "F3: characterised, raw counts, footprint", lab, int((f4 & m).sum()), None,
            logged=lg)
    for thr, lg in ((0.1, 143), (0.3, 319)):
        mj = (d.relMl_J > 0) & (d.relMl_J < thr)
        mr = (d.relMl_R > 0) & (d.relMl_R < thr)
        ml = (d.relMl_L > 0) & (d.relMl_L < thr)
        rec(pop, "H5: mass only from the combination, raw counts, footprint",
            f"sigma(M)/M<{thr:g} jointly and from neither alone",
            int((f4 & (mj & ~mr & ~ml).to_numpy()).sum()), None, logged=lg)
    # --- astrometric shift, footprint (h5) ---
    s = max_centroid_shift(d)
    q = b.quantile(s, f4 & (s > 0), W)
    rec(pop, "H5 (footprint)", "median max shift [mas]", q[0], q[1], n=q[2], neff=q[3],
        logged=0.1236)
    # --- dust and source brightness, footprint ---
    for lab, m in (("Roman detections", f4 & (d.detR == 1).to_numpy()), ("all footprint", f4)):
        q = b.quantile(d["Ai_r"].to_numpy(float), m, W)
        rec(pop, "dust (footprint)", f"median A_r [mag], {lab}", q[0], q[1], n=q[2], neff=q[3])
    src = (d["magb_F146"] - 2.5 * np.log10(d["blend_F146"].clip(lower=1e-12))).to_numpy(float)
    m = f4 & (d.detR == 1).to_numpy()
    q = b.quantile(src, m, W)
    rec(pop, "source brightness (Roman detections)", "median source F146 [mag]", q[0], q[1],
        n=q[2], neff=q[3], logged=21.8)
    q = b.quantile(d["magb_F146"].to_numpy(float), m, W)
    rec(pop, "source brightness (Roman detections)", "median baseline (blended) F146 [mag]",
        q[0], q[1], n=q[2], neff=q[3])
    # --- gap filling (Table 3): 10-30 d footprint events peaking in a gap or a season ---
    base = (d.t0zone.isin([0, 1]) & (d.ndw_R > 0) & (d.dt_edge <= 60.0)
            & (d.tE >= 10.0) & (d.tE < 30.0)).to_numpy()
    gap = base & (d.t0zone == 1).to_numpy()
    sea = base & (d.t0zone == 0).to_numpy()
    dt = d["dt_edge"].to_numpy(float)
    for p in ("tE", "piE"):
        r = R.ratio_joint_over(d, p, "roman").to_numpy(float)
        for lab, m in (("anywhere in a gap", gap), (">15 d past an edge", gap & (dt > 15)),
                       (">30 d past", gap & (dt > 30)), (">45 d past", gap & (dt > 45)),
                       ("in a season (control)", sea)):
            q = b.quantile(r, m & np.isfinite(r), W)
            rec(pop, "gap filling: median sigma_J/sigma_Roman, 10-30 d", f"{p}, {lab}", q[0],
                q[1], n=q[2], neff=q[3], logged=logged["gap"].get((p, lab)))
    # What Roman has on these events. "In a gap" does not mean Roman has no data: every event
    # here is on a sightline Roman observes all mission and has a Roman forecast (the ratio is
    # defined). What changes with depth is whether Roman saw the magnified part, and whether it
    # still detects the event from the wings. Same rows as above, over the tE-ratio sample.
    rt = R.ratio_joint_over(d, "tE", "roman").to_numpy(float)
    detR, detL = (d.detR == 1).to_numpy(), (d.detL == 1).to_numpy()
    nopk = (d.nepR_pk == 0).to_numpy()
    for lab, m in (("anywhere in a gap", gap), (">15 d past an edge", gap & (dt > 15)),
                   (">30 d past", gap & (dt > 30)), (">45 d past", gap & (dt > 45)),
                   ("in a season (control)", sea)):
        m = m & np.isfinite(rt)
        for what, sel in (("Roman detects", detR), ("Rubin detects", detL),
                          ("no Roman epoch within +-2 tE", nopk)):
            f = b.frac(m & sel, m, W)
            rec(pop, "gap filling: Roman's coverage of the same events, 10-30 d",
                f"{lab}: {what} [%]", f[0], f[1], n=f[2], neff=f[3], scale=100)
    # Improvement factor deep in the gap, 1/ratio, with its own bootstrap (not 1/err).
    r = R.ratio_joint_over(d, "tE", "roman").to_numpy(float)
    with np.errstate(divide="ignore"):
        inv = 1.0 / r
    q = b.quantile(inv, gap & (dt > 45) & np.isfinite(inv), W)
    rec(pop, "gap filling: median sigma_J/sigma_Roman, 10-30 d",
        "tE, >45 d past: median sigma_Roman/sigma_J", q[0], q[1], n=q[2], neff=q[3])
    # Rescue: characterised jointly and not by Roman, among in-gap footprint detections.
    ingap = (d.t0zone == 1).to_numpy() & (d.ndw_R > 0).to_numpy() & (dt <= 60)
    resc = cJ & ~cR
    f = b.frac(ingap & resc, ingap, W)
    rec(pop, "gap filling: rescue", "in-gap footprint detections, joint-only char. [%]", f[0],
        f[1], n=f[2], neff=f[3], scale=100, logged=2.0)
    for lo, hi in ((10, 30), (30, 100)):
        m = ingap & (dt > 45) & (d.tE >= lo).to_numpy() & (d.tE < hi).to_numpy()
        f = b.frac(m & resc, m, W)
        rec(pop, "gap filling: rescue", f"tE {lo}-{hi} d, >45 d past an edge [%]", f[0], f[1],
            n=f[2], neff=f[3], scale=100)
    # --- 'Rubin detects, Roman does not' in the footprint: where do they peak? ---
    m = d.foot.to_numpy() & (d.detL == 1).to_numpy() & (d.detR == 0).to_numpy()
    for z, lab in ((0, "in a Roman season"), (1, "in a gap"), (2, "outside Roman's mission")):
        rec(pop, "Rubin detects, Roman does not (footprint, raw)", lab,
            int((m & (d.t0zone == z).to_numpy()).sum()), None,
            logged={0: 1, 1: 116, 2: 800}[z])
    # --- astrometric maximum across a season edge (h5 panel c) is left as logged ---


def satellite(run, logged):
    """H3: sigma(piE) with Roman at L2 over Roman at Earth, events Roman covers at the peak."""
    path = os.path.join(run.dir, "h3_pair.dat")
    p = R.load_events(path, usecols=PAIR_COLS)
    w, _ = R.attach_weight(p, run.map, run.logs)
    p["W"] = w.to_numpy()
    ok = ((p.okA_sat == 1) & (p.okA_nosat == 1) & (p.sigpiE_sat > 0)
          & (p.sigpiE_nosat > 0)).to_numpy()
    r = (p.sigpiE_sat / p.sigpiE_nosat).to_numpy(float)
    cov = ok & (p.nepR_pk > 0).to_numpy()
    b = Boot(len(p), seed=SEED + 1)
    W = p["W"].to_numpy()
    q = b.quantile(r, cov, W)
    rec(run.name, "satellite parallax (Roman covers the peak)", "median sigma_L2/sigma_Earth",
        q[0], q[1], n=q[2], neff=q[3], logged=logged.get("median"))
    for lab, m, lg in (("improves at all [%]", r < 1.0, 86.8), ("improves >1% [%]", r < 0.99, 40.7),
                       ("improves >2x [%]", r < 0.5, 1.2)):
        f = b.frac(cov & m, cov, W)
        rec(run.name, "satellite parallax (Roman covers the peak)", lab, f[0], f[1], n=f[2],
            neff=f[3], scale=100, logged=lg)
    te = p["tE"].to_numpy(float)
    m = cov & (te >= 10) & (te < 30)
    q = b.quantile(r, m, W)
    rec(run.name, "satellite parallax (Roman covers the peak)", "median, tE 10-30 d", q[0], q[1],
        n=q[2], neff=q[3], logged=0.9914)
    ctrl = ok & (p.nepR_pk == 0).to_numpy()
    q = b.quantile(r, ctrl, W)
    rec(run.name, "satellite parallax (control, no Roman data at the peak)", "median", q[0],
        q[1], n=q[2], neff=q[3], logged=1.0)


def intrinsic_te(name, path):
    """W1's check on the post-fix draws: the event-rate-weighted mean tE over ALL draws.

    As in w1_intrinsic_te.py: nsim is the sightline's row count (every draw writes a row,
    barren sightlines included), and Z is interpolated on a 200-point Ds grid with its error
    measured. The comparison with OGLE-IV (Mroz et al. 2019, Fig. 13) is made both over the
    whole scan and over the region OGLE's central longitude bins cover, -6 <= b <= -1, |l| < 2.
    """
    print(f"[{name}] intrinsic tE: reading {path}", flush=True)
    df = R.load_events(path, chunksize=1_000_000, narrow=True)
    codes, keys = R.sightline_index(df)
    nsim = np.bincount(codes)
    W = np.empty(len(df))
    Ds = df["Ds"].to_numpy(np.float64)
    fac = (df["w_area"].to_numpy(np.float64) * np.sqrt(df["Ml"].to_numpy(np.float64))
           * df["Vt"].to_numpy(np.float64))
    rng = np.random.default_rng(0)
    worst = 0.0
    for c, pos in R.sightline_groups(codes):
        prof = G.density_profile(*keys[c])
        ds = Ds[pos]
        grid = np.linspace(ds.min(), ds.max(), 200)
        Z = np.interp(ds, grid, G.lens_distance_norm(prof, grid))
        if worst == 0.0:
            sub = rng.choice(len(ds), size=min(300, len(ds)), replace=False)
            worst = float(np.max(np.abs(Z[sub] / G.lens_distance_norm(prof, ds[sub]) - 1)))
        W[pos] = fac[pos] * prof.Nstart / nsim[c] * Z
    print(f"[{name}] Z interpolation worst relative error {worst:.1e}", flush=True)
    tE = df["tE"].to_numpy(np.float64)
    lon, lat = df["lon"].to_numpy(np.float64), df["lat"].to_numpy(np.float64)
    rng = np.random.default_rng(SEED + 2)
    regions = (("whole scan", np.ones(len(df), bool)),
               ("-6<=b<=-1, |l|<2 (OGLE central bins)",
                (lat >= -6) & (lat <= -1) & (np.abs(lon) < 2)),
               ("-6<=b<=-1, whole longitude range", (lat >= -6) & (lat <= -1)))
    for lab, m in regions:
        w, t = W[m], tE[m]
        mean_w = float((w * t).sum() / w.sum())
        mean_u = float(t.mean())
        reps_w, reps_u = [], []
        for _ in range(200):
            c = rng.poisson(1.0, size=t.size)
            reps_w.append((w * c * t).sum() / (w * c).sum())
            reps_u.append((c * t).sum() / c.sum())
        rec(name, "intrinsic tE, all draws", f"weighted mean tE [d], {lab}", mean_w,
            float(np.std(reps_w)), n=int(m.sum()), neff=R.kish_neff(w),
            logged=24.0 if lab == "whole scan" else None)
        rec(name, "intrinsic tE, all draws", f"unweighted mean tE [d], {lab}", mean_u,
            float(np.std(reps_u)), n=int(m.sum()), logged=56.2 if lab == "whole scan" else None)
    # Dust over every draw in the footprint (the report's 'footprint median A_V').
    foot = df["ndw_R"].to_numpy() > 0
    b = Boot(int(foot.sum()), B=200, seed=SEED + 3)
    q = b.quantile(df["Ai_r"].to_numpy(np.float64)[foot], np.ones(int(foot.sum()), bool),
                   W[foot])
    rec(name, "dust (footprint)", "median A_r [mag], all draws with Roman epochs", q[0], q[1],
        n=q[2], neff=q[3])


# ---------------------------------------------------------------------------------------------
# Values the source scripts logged, for the comparison column
# ---------------------------------------------------------------------------------------------
LOGGED = {
    "bh": dict(share={"Rubin only": 72.88, "Roman only": 21.45, "both": 5.63},
               jos={"N": 23961, ("tE", "Rubin"): 0.022, ("tE", "Roman"): 0.997,
                    ("piE", "Rubin"): 0.010, ("piE", "Roman"): 0.998,
                    ("tetE", "Rubin"): 0.041, ("tetE", "Roman"): 0.999},
               shift={"median": 3.6023, "floor": 90.43},
               res={("Rubin", "nres5"): 5.4540, ("Rubin", "nres20"): 0.0310,
                    ("Rubin", "nresPSF"): 0.0003, ("Roman", "nres5"): 65.5078,
                    ("Roman", "nres20"): 30.9002, ("Roman", "nresPSF"): 0.6788},
               p7={("Ml", "J"): 4.16, ("Ml", "L"): 0.27, ("Ml", "R"): 12.60,
                   ("piE", "J"): 5.32, ("piE", "L"): 1.46, ("piE", "R"): 12.85,
                   ("tE", "J"): 18.66, ("tE", "L"): 11.74, ("tE", "R"): 23.64,
                   ("tetE", "J"): 48.84, ("tetE", "L"): 37.53, ("tetE", "R"): 96.97}),
    "ns": dict(share={"Rubin only": 69.59, "Roman only": 26.27, "both": 3.92},
               jos={"N": 18895, ("tE", "Rubin"): 0.008, ("tE", "Roman"): 0.999,
                    ("piE", "Rubin"): 0.008, ("piE", "Roman"): 0.999,
                    ("tetE", "Rubin"): 0.035, ("tetE", "Roman"): 0.999},
               shift={"median": 0.2637, "floor": 0.58},
               res={("Rubin", "nres5"): 0.0002, ("Rubin", "nres20"): 0.0,
                    ("Rubin", "nresPSF"): 0.0, ("Roman", "nres5"): 14.8533,
                    ("Roman", "nres20"): 0.3516, ("Roman", "nresPSF"): 0.0011},
               p7={("Ml", "J"): 1.36, ("Ml", "L"): 0.00, ("Ml", "R"): 3.66,
                   ("piE", "J"): 2.40, ("piE", "L"): 0.78, ("piE", "R"): 4.55,
                   ("tE", "J"): 10.01, ("tE", "L"): 4.45, ("tE", "R"): 16.03,
                   ("tetE", "J"): 22.12, ("tetE", "L"): 0.43, ("tetE", "R"): 67.17}),
    "bulge": dict(share={"Rubin only": 67, "Roman only": 29, "both": 4},
                  f4={("tE", "J"): 13.2, ("tE", "R"): 11.6, ("tE", "L"): 0.8,
                      ("piE", "J"): 3.7, ("piE", "R"): 3.3, ("piE", "L"): 0.2,
                      ("tetE", "J"): 50.5, ("tetE", "R"): 50.3, ("tetE", "L"): 0.1,
                      ("Ml", "J"): 2.4, ("Ml", "R"): 2.2, ("Ml", "L"): 0.001},
                  gap={("tE", "anywhere in a gap"): 0.973, ("piE", "anywhere in a gap"): 0.984,
                       ("tE", ">15 d past an edge"): 0.826, ("piE", ">15 d past an edge"): 0.892,
                       ("tE", ">30 d past"): 0.436, ("piE", ">30 d past"): 0.703,
                       ("tE", ">45 d past"): 0.077, ("piE", ">45 d past"): 0.183,
                       ("tE", "in a season (control)"): 1.000,
                       ("piE", "in a season (control)"): 1.000},
                  sat={"median": 0.9911}),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", action="append", required=True, metavar="NAME=DIR")
    ap.add_argument("--mean-mass", required=True,
                    help="y1 output directory holding <name>/y1_yields.md")
    ap.add_argument("--draws", action="append", default=[], metavar="NAME=FILE",
                    help="all-draws column extract, for the intrinsic-timescale check")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    for spec in a.draws:
        name, path = spec.split("=", 1)
        intrinsic_te(name, path)

    for spec in a.run:
        name, directory = spec.split("=", 1)
        run = Run(name, directory, a.mean_mass)
        lg = LOGGED.get(name, {})
        yields_block(run)
        detection_shares(run, lg.get("share"))
        joint_over_single(run, lg.get("jos"))
        astrometric_shift(run, np.ones(len(run.df), bool), "all detections", lg.get("shift"))
        astrometric_shift(run, run.df["foot"].to_numpy(), "footprint detections")
        precision_10pct(run, lg.get("p7"))
        resolution(run, lg.get("res"))
        roman_helps_rubin(run)
        peak_coverage(run)
        if name == "bulge":
            bulge_specific(run, lg)
            satellite(run, lg["sat"])
        del run

    out = pd.DataFrame(ROWS)
    out.to_csv(os.path.join(a.out, "u1_numbers.csv"), index=False)
    lines = ["# Report numbers with Monte Carlo errors (Step U1)", "",
             "Generated by `analysis/u1_report_numbers.py`. `err` is the Poisson-bootstrap "
             f"standard deviation ({N_BOOT} replicates) for fractions and medians, and "
             "sqrt(sum y^2) for yields (y1's convention). `logged` is the value the original "
             "script printed, where there is one: the recomputed value must reproduce it.", ""]
    for (pop, sec), g in out.groupby(["population", "section"], sort=False):
        lines += [f"## {pop}: {sec}", "", "| quantity | value | err | n | N_eff | logged |",
                  "|---|---:|---:|---:|---:|---:|"]
        for _, r in g.iterrows():
            def f(x, fmt="{:.4g}"):
                return "" if x is None or (isinstance(x, float) and np.isnan(x)) else fmt.format(x)
            lines.append(f"| {r.quantity} | {f(r.value)} | {f(r.err, '{:.2g}')} | "
                         f"{f(r.n, '{:,}')} | {f(r.n_eff, '{:,.0f}')} | {f(r.logged)} |")
        lines.append("")
    with open(os.path.join(a.out, "u1_numbers.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
