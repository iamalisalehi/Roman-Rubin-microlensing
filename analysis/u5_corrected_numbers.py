#!/usr/bin/env python3
"""Step U5: the report's headline numbers with the dust corrected, and their uncertainties.

WHY THIS EXISTS. Step U4 (Deviation 61) found the model's extinction three to four times too low
within ~1 deg of the plane (optical Bayestar/DECaPS saturate; checked against the near-infrared
Marshall et al. 2006 map) and estimated the size of the effect. The user asked for the abstract
and summary to report the CORRECTED results, with uncertainties. This script produces them, for
every number the abstract and summary quote, each with a Monte Carlo error and a systematic error
of the correction itself.

THE DIAGNOSIS IT RESTS ON (Deviation 63, all checked, not assumed): maps.py follows the dustmaps
documentation's combination rule, Bayestar north of dec -30 and DECaPS south. 139 of Roman's 147
sightlines are north of -30, and Bayestar's OWN reliable-distance flag is false at >= 4 kpc on every
one of them: the tables hold Bayestar values beyond the range Bayestar vouches for, for every bulge
source. DECaPS (Zucker et al. 2025; 239 < l < 6, |b| < 10, which contains the whole scan) has data at
every sightline and table position, and its flag is true at 8 kpc on all of Roman's five-field block
and 71% of the Galactic-centre field. But DECaPS is sensitive to A_V ~ 12 only, and on the
Galactic-centre field it is saturated despite its flag: its column there is 1.07 times its
five-field column, where the VVV reddening map (Surot et al. 2020) measures 4.89 and Marshall 4.46
(Step U6, analysis/u6_vvv_check.py; the ratio is independent of the extinction law). 78 tables are
empty (Bayestar has no data there) and were simulated with zero dust.

THE CORRECTION. For every simulated draw i, dA_V,i = A_V,reference(l, b, Ds_i) - A_V,model(l, b, Ds_i)
at the source's OWN distance, the reference being the nominal "hybrid" dust of the Dust class
(DECaPS where it is reliable and not saturated; the near-infrared map, calibrated onto DECaPS's
scale where both are valid, where it is not). The extra dust only dims the
source: dF146 = 0.197 dA_V, dr = 0.854 dA_V (CCM89, R_V = 2.5). For an outcome S (detected by Roman,
mass to 10%, ...) that depends on the source brightness through a survey's magnitude m, the
efficiency eff_S(m) is measured from the same draws, and each draw with S_i = 1 is re-weighted by
eff_S(m_i + dm_i) / eff_S(m_i): its probability of still having outcome S under the thicker dust,
relative to now. Summing y_i over S-draws with that factor gives the corrected yield of S. Composite
outcomes use the product of the two surveys' factors (Rubin's in r, Roman's in F146), which assumes
the two respond independently to the dimming. Per-event fractions are ratios of corrected yields;
medians are taken over detections re-weighted by their detection factor. The two-image resolution
count is redone epoch by epoch with the dimmed F146 magnitudes (u2_resolution_depth.count_epochs, at
Roman's 5-sigma depth).

THE UNCERTAINTY of every corrected number is
    MC   : Poisson on the draws, sqrt(sum (y s)^2), for yields; a Poisson bootstrap for fractions
           and medians;
    syst : the spread of the correction under its own choices, added in quadrature --
           reference dust: the shift of the near-infrared-only variant (a DECaPS-only variant is
           computed and reported, but not counted: VVV rules it out on the Galactic-centre field),
           A_Ks/A_V: the larger shift of 0.0734 (VVV contrast) and 0.102 (CCM89, R_V 2.5) from the
           nominal 0.0805,
           efficiency measured per field block (nominal) vs pooled over the footprint.
What it does NOT cover: the law converting A_V to A_F146 and A_r is the simulator's CCM89 at
R_V = 2.5 throughout; its own uncertainty is not included.
Nor: that a surviving event's PRECISION also degrades when its source dims. For
yields of precision outcomes (masses to 10%) that is included, because they use their own
efficiency; for re-weighted medians of per-event ratios it is not.

INPUTS per population (run directory): test<tag>_foot.dat (all footprint draws; Step U4) and
test<tag>_rubincols.dat (all draws, lon lat w_area Ml Vt Ds struc magb_r detL detJ) for the
whole-scan Rubin numbers; test<tag>_detJ.dat for whole-scan per-event quantities.

    .roman/bin/python analysis/u5_corrected_numbers.py --run bulge=runs/prod_bulge_20260924 \\
        --run bh=runs/prod_bh_20260924 --run ns=runs/prod_ns_20260924 \\
        --mean-mass figures/yield_20260925 -o figures/u1_20260929
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import galaxy_model as G              # noqa: E402
import romanlib as R                  # noqa: E402
import u1_report_numbers as U         # noqa: E402

AF146_AV, AR_AV = 0.197, 0.854
# A_Ks/A_V puts the near-infrared (Marshall) map on the A_V scale; it matters only where the
# reference uses that map, i.e. where DECaPS cannot see (Step U6, analysis/u6_vvv_check.py):
#   nominal 0.0805: the median A_Ks(Marshall)/A_V(DECaPS) over Roman's five-field block at 8 kpc,
#                   where DECaPS is flagged reliable on all 119 sightlines -- the near-infrared map
#                   calibrated onto DECaPS's scale where both are valid (16-84%: 0.068-0.093);
#   aks_lo  0.0734: 0.0805 x 4.46/4.89, i.e. Marshall's column raised by the ~10% by which its
#                   Galactic-centre/five-field contrast (4.46) falls short of the VVV reddening
#                   map's (4.89, Surot et al. 2020) -- more dust; the GC-field column is then 4%
#                   above VVV's calibrated one, the nominal 5% below it;
#   aks_hi  0.102 : CCM89 at R_V 2.5, the law the simulator applies to bulge sources -- the
#                   uncalibrated alternative, less dust.
AKS_AV = {"nominal": 0.0805, "aks_lo": 0.0734, "aks_hi": 0.102}
# DECaPS's stated sensitivity limit (Zucker et al. 2025): beyond the distance at which the
# (calibrated) near-infrared column passes it, DECaPS no longer sees the stars behind the dust.
DECAPS_AV_MAX = 12.0
# (reference variant, A_Ks/A_V, pooled efficiency): nominal first; the rest set the systematic error,
# except "decaps", which is kept for the record only: VVV rules it out on the Galactic-centre field.
VARIANTS = {"nominal": ("hybrid", AKS_AV["nominal"], False),
            "marshall": ("marshall", AKS_AV["nominal"], False),
            "decaps": ("decaps", AKS_AV["nominal"], False),
            "aks_lo": ("hybrid", AKS_AV["aks_lo"], False),
            "aks_hi": ("hybrid", AKS_AV["aks_hi"], False),
            "pooled": ("hybrid", AKS_AV["nominal"], True)}
DEPTH_5SIG = 25.52
T_S = (3652.43 - 2 * R.T0_MARGIN_DAYS) * 86400.0
EFF_BINS = np.arange(10.0, 34.01, 0.25)
DGRID = np.arange(0.05, 20.0, 0.05)          # kpc, for per-sightline dust profiles
FOOT_COLS = ["lon", "lat", "w_area", "Ml", "Vt", "Ds", "struc", "detL", "detR", "detJ",
             "magb_F146", "magb_r", "blend_F146", "tE", "piE", "tetE", "u0", "t0",
             "okA_J", "okA_L", "okA_R", "okB_J", "okB_L", "okB_R",
             "sigtE_J", "sigtE_L", "sigtE_R", "sigpiE_J", "sigpiE_L", "sigpiE_R",
             "sigtetE_J", "sigtetE_L", "sigtetE_R", "relMl_J", "relMl_L", "relMl_R",
             "ndw_R", "Ai_r"]
ROWS = []
LOG_FOR_SUMMARY = None


# ---------------------------------------------------------------------------------------------
# Dust
# ---------------------------------------------------------------------------------------------
class Dust:
    """The simulator's extinction (its own tables, chosen exactly as the C++ chooses them) and the
    reference dust the correction moves each source to.

    MODEL. readBayestar() walks files/ext/ in directory order and takes each table's l, b from its
    first line; nearestSightline() keeps the first table with the strictly smallest flat distance.
    os.scandir() returns the same (readdir) order, and np.argmin returns the first minimum, so the
    choice is identical, ties included (verified: 99.92% of 300,000 draws' recorded A_r agree to the
    0.017-mag scatter the simulator adds; the rest were a tie this now reproduces). An empty table
    (all NaN; 78 of them) fails the C++ read and gives zero dust, and is modelled as zero.

    REFERENCE, three variants (nir = Marshall's A_Ks / A_Ks/A_V):
      "hybrid"   (nominal) -- what the fixed pipeline would use: DECaPS (A_V = 3.32 E(B-V), the
                 dustmaps convention) out to the largest distance its own flag calls reliable at
                 that sightline; beyond it, DECaPS's last reliable value plus nir's further increase;
                 and from the first distance at which nir reaches DECAPS_AV_MAX on, nir itself --
                 DECaPS is saturated there (on the Galactic-centre field it reads A_V ~ 7 at 8 kpc,
                 flat from ~3 kpc, while nir and the VVV reddening map read ~31-33; Step U6);
      "decaps"   DECaPS at every distance, ignoring its reliability flag and its saturation;
      "marshall" nir at every distance.
    self.dsat[(l, b, A_Ks/A_V)] records the distance at which the nominal switches to nir (inf if
    it never does), for the summary.
    """

    def __init__(self):
        from dustmaps.config import config
        config["data_dir"] = "dustmaps"
        from dustmaps.marshall import MarshallQuery
        from dustmaps.decaps import DECaPSQueryLite
        self.mq = MarshallQuery()
        self.dq = DECaPSQueryLite(mean_only=True)
        self.files = [e.path for e in os.scandir("files/ext")
                      if e.is_file() and e.name.endswith(".txt")]
        self.pos = np.array([np.loadtxt(f, max_rows=1, usecols=(0, 1)) for f in self.files])
        self._model, self._aks, self._dec = {}, {}, {}
        self.dsat = {}

    def _nearest(self, l, b):
        return int(np.argmin((self.pos[:, 0] - l) ** 2 + (self.pos[:, 1] - b) ** 2))

    def model_av(self, l, b, d):
        """A_V at distances d on the table the simulator used for sightline (l, b)."""
        i = self._nearest(l, b)
        if i not in self._model:
            a = np.loadtxt(self.files[i], usecols=(2, 3))
            self._model[i] = None if not np.isfinite(a[:, 1]).any() else a
        t = self._model[i]
        return np.zeros(np.size(d)) if t is None else np.interp(d, t[:, 0], t[:, 1])

    def _coords(self, l, b):
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        return SkyCoord(l=np.full(DGRID.size, l) * u.deg, b=np.full(DGRID.size, b) * u.deg,
                        distance=DGRID * u.kpc, frame="galactic")

    def aks_profile(self, l, b):
        """Marshall A_Ks on DGRID; NaN beyond its coverage is held at the last valid value."""
        k = (round(l, 3), round(b, 3))
        if k not in self._aks:
            a = np.asarray(self.mq(self._coords(l, b)), float)
            ok = np.isfinite(a)
            if ok.sum() >= 2:
                a = np.interp(DGRID, DGRID[ok], a[ok])
            elif ok.sum() == 1:
                a = np.full(DGRID.size, a[ok][0])
            else:
                a = np.full(DGRID.size, np.nan)
            self._aks[k] = a
        return self._aks[k]

    def decaps_profile(self, l, b):
        """(A_V on DGRID, reliable mask) from DECaPS, A_V = 3.32 E(B-V)."""
        k = (round(l, 3), round(b, 3))
        if k not in self._dec:
            v, fl = self.dq(self._coords(l, b), mode="mean", return_flags=True)
            self._dec[k] = (3.32 * np.asarray(v, float), np.asarray(fl["reliable_dist"], bool))
        return self._dec[k]

    def reference_profile(self, l, b, variant, aks_av):
        aks = self.aks_profile(l, b)
        nir = aks / aks_av
        if variant == "marshall":
            return nir
        dav, rel = self.decaps_profile(l, b)
        if not np.isfinite(dav).any():
            return nir                                # outside DECaPS (not the case in this scan)
        dav = np.interp(DGRID, DGRID[np.isfinite(dav)], dav[np.isfinite(dav)])
        if variant == "decaps":
            return dav
        if np.all(np.isnan(aks)):
            return dav
        if not rel.any():
            return nir
        imax = int(np.flatnonzero(rel).max())
        prof = dav.copy()
        beyond = np.arange(DGRID.size) > imax
        prof[beyond] = dav[imax] + np.maximum(aks[beyond] - aks[imax], 0.0) / aks_av
        hit = np.flatnonzero(nir >= DECAPS_AV_MAX)
        self.dsat[(round(l, 3), round(b, 3), aks_av)] = DGRID[hit[0]] if hit.size else np.inf
        if hit.size:
            prof[hit[0]:] = nir[hit[0]:]
        return prof

    def delta_av(self, lon, lat, ds, variant, aks_av, codes=None, keys=None):
        """Per-draw A_V,reference - A_V,model at each source's own distance."""
        lon, lat, ds = map(np.asarray, (lon, lat, ds))
        out = np.zeros(lon.size)
        if codes is None:
            sr = pd.Series(list(zip(np.round(lon, 3), np.round(lat, 3))))
            groups = [(k, np.asarray(p)) for k, p in sr.groupby(sr).groups.items()]
        else:
            groups = [(keys[c], pos) for c, pos in R.sightline_groups(codes)]
        for k, pos in groups:
            ref = np.interp(ds[pos], DGRID, self.reference_profile(*k, variant, aks_av))
            out[pos] = ref - self.model_av(*k, ds[pos])
        return out


def efficiency(mag, succ, y):
    """Yield-weighted P(success) in magnitude bins, interpolated; 0 fainter than the data."""
    idx = np.clip(np.digitize(mag, EFF_BINS) - 1, 0, len(EFF_BINS) - 2)
    num = np.bincount(idx, weights=y * succ, minlength=len(EFF_BINS) - 1)
    den = np.bincount(idx, weights=y, minlength=len(EFF_BINS) - 1)
    ok = den > 0
    cen = 0.5 * (EFF_BINS[1:] + EFF_BINS[:-1])[ok]
    eff = num[ok] / den[ok]
    return lambda m: np.interp(m, cen, eff, left=eff[0], right=0.0)


def survival(mag, dmag, succ, y, groups):
    """Per-draw factor eff(m + dm) / eff(m) for draws with succ, efficiency measured per group."""
    s = np.zeros(mag.size)
    for g in np.unique(groups):
        m = groups == g
        eff = efficiency(mag[m], succ[m], y[m])
        e0, e1 = eff(mag[m]), eff(mag[m] + dmag[m])
        s[m] = np.where(e0 > 0, np.minimum(e1 / np.where(e0 > 0, e0, 1), 10.0), 0.0)
    # A draw that did not have the outcome contributes nothing, even where its magnitude is
    # undefined (NaN x 0 would otherwise poison every sum). A NaN on a draw that DID have it is a
    # real problem, so it is counted and reported rather than silently zeroed.
    bad = (succ > 0) & ~np.isfinite(s)
    if bad.any():
        print(f"  WARNING: {int(bad.sum())} successful draws with an undefined correction "
              f"(magnitude or dust NaN); treated as unchanged", flush=True)
        s[bad] = 1.0
    return np.where(succ > 0, s, 0.0)


# ---------------------------------------------------------------------------------------------
# One variant of the correction on one population's footprint
# ---------------------------------------------------------------------------------------------
def footprint_variant(d, y, dav, groups):
    """All footprint quantities for one set of correction choices. Returns {name: (value, mc)}."""
    F, r = d["magb_F146"].to_numpy(float), d["magb_r"].to_numpy(float)
    dF, dr = AF146_AV * dav, AR_AV * dav
    detR = (d.detR == 1).to_numpy(); detL = (d.detL == 1).to_numpy()
    mR = detR & (d.relMl_R > 0).to_numpy() & (d.relMl_R < 0.1).to_numpy()
    mR1 = detR & (d.relMl_R > 0).to_numpy() & (d.relMl_R < 0.01).to_numpy()
    mJ = (d.relMl_J > 0).to_numpy() & (d.relMl_J < 0.1).to_numpy() & (d.detJ == 1).to_numpy()
    chL = R.characterized(d, "rubin").to_numpy() & detL
    chJ = R.characterized(d, "joint").to_numpy()
    mL = detL & (d.relMl_L > 0).to_numpy() & (d.relMl_L < 0.1).to_numpy()
    f = lambda succ, mag, dm: survival(mag, dm, succ.astype(float), y, groups)
    sR, sL = f(detR, F, dF), f(detL, r, dr)
    out = {}

    def yld(name, w):
        out[name] = (float(np.sum(y * w)), float(np.sqrt(np.sum((y * w) ** 2))))

    yld("Roman detects", sR)
    yld("Roman, sigma(M)/M<10%", f(mR, F, dF))
    yld("Roman, sigma(M)/M<1%", f(mR1, F, dF))
    yld("joint, sigma(M)/M<10%", f(mJ, F, dF))
    yld("Rubin detects (footprint)", sL)
    both = sR * sL                      # non-zero only where both detect (each is 0 elsewhere)
    yld("both detect", both)
    yld("Rubin only (footprint)", sL * ~detR)
    # Roman for Rubin: over Rubin's footprint detections
    # conditional Roman-side factors, measured among Rubin's detections
    def cond(succ):
        s = np.zeros(d.shape[0])
        m = detL
        s[m] = survival(F[m], dF[m], succ[m].astype(float), y[m], groups[m])
        return s
    wts = {}
    wts["RfR: characterised, Rubin alone"] = f(chL, r, dr)
    wts["RfR: characterised, joint"] = sL * cond(chJ & detL)
    wts["RfR: mass to 10%, Rubin alone"] = f(mL, r, dr)
    wts["RfR: mass to 10%, joint"] = sL * cond(mJ & detL)
    for k, v in wts.items():
        yld(k, v)

    def ratio(a, b, scale=1.0):
        # sum(y a) / sum(y b), MC error by the delta method: sum (y a - rho y b)^2 / (sum y b)^2
        A, B = y * a, y * b
        rho = A.sum() / B.sum()
        return (scale * rho, scale * float(np.sqrt(np.sum((A - rho * B) ** 2)) / B.sum()))

    for key in ("characterised", "mass to 10%"):
        for who in ("Rubin alone", "joint"):
            out[f"RfR: {key}, {who} [%]"] = ratio(wts[f"RfR: {key}, {who}"], sL, 100)
    out["joint/Roman, sigma(M)/M<10%"] = ratio(f(mJ, F, dF), f(mR, F, dF))
    out["Roman also detects, of Rubin's [%]"] = ratio(both, sL, 100)
    out["Roman sigma(M)/M<10% fraction of its detections [%]"] = ratio(f(mR, F, dF), sR, 100)
    # re-weighted medians
    W = d["W"].to_numpy()
    boot = U.Boot(d.shape[0], B=200, seed=U.SEED + 30)
    def wmed(vals, sel, w):
        q = boot.quantile(vals, sel, w)
        return (q[0], q[1])
    tE_J, tE_L = R.sigma(d, "tE", "joint").to_numpy(float), R.sigma(d, "tE", "rubin").to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rt = tE_J / tE_L
    out["median sigma_J/sigma_Rubin (tE), Rubin's detections"] = wmed(rt, detL & np.isfinite(rt), W * sL)
    tE_R = R.sigma(d, "tE", "roman").to_numpy(float)
    both_c = ((d.okA_L == 1) & (d.okA_R == 1)).to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        rr = tE_J / tE_R
    out["median sigma_J/sigma_Roman (tE), both constrain"] = wmed(rr, both_c & np.isfinite(rr), W * np.maximum(sR, sL))
    shift = U.max_centroid_shift(d)
    out["median max shift [mas], Roman's detections"] = wmed(shift, detR & (shift > 0), W * sR)
    # dust itself: median corrected A_V of Roman's footprint detections
    av_corr = d["Ai_r"].to_numpy(float) / AR_AV + dav
    out["median A_V, Roman's footprint detections (corrected)"] = wmed(av_corr, detR, W * sR)
    out["_sR"], out["_sL"] = sR, sL
    return out


def resolution_variant(d, y, dav, sR):
    """P(resolvable) for Roman at its 5-sigma depth with the dimmed magnitudes, over detections
    re-weighted by their detection factor."""
    import u2_resolution_depth as U2
    det = (d.detR == 1).to_numpy()
    ev = d[det].reset_index(drop=True).copy()
    ev["magb_F146"] = ev["magb_F146"].to_numpy(float) + AF146_AV * dav[det]
    times = U2.roman_times()
    cnt = U2.count_epochs(ev, times, (DEPTH_5SIG,))
    w = (y * sR)[det]
    out = {}
    for bar in ("5", "20", "PSF"):
        res = cnt[(DEPTH_5SIG, bar)] >= U.RESOLVE_MIN_EPOCHS
        # MC error of the weighted fraction by the delta method, as ratio() in footprint_variant
        rho = float(w[res].sum() / w.sum())
        mc = float(np.sqrt(np.sum((w * res - rho * w) ** 2)) / w.sum())
        out[f"P(resolvable) Roman D={bar} [%]"] = (100 * rho, 100 * mc)
        out[f"Roman resolves, D={bar}, N_1"] = (float(w[res].sum()), float(np.sqrt((w[res] ** 2).sum())))
    return out


# ---------------------------------------------------------------------------------------------
# Whole scan (Rubin)
# ---------------------------------------------------------------------------------------------
def whole_scan_context(name, directory, mean_mass_dir):
    """Load every draw once (lon lat w_area Ml Vt Ds struc magb_r detL detJ) and rebuild y."""
    import subprocess
    tag = R.load_provenance(os.path.join(directory, "files/MONTLMC/files/run_provenance.txt")
                            ).get("population_tag", "5")
    d = R.load_events(os.path.join(directory, f"test{tag}_rubincols.dat"), chunksize=1_000_000,
                      narrow=True)
    codes, keys = R.sightline_index(d)
    nsim = np.bincount(codes)
    Ds = d["Ds"].to_numpy(np.float64)
    fac = (d["w_area"].to_numpy(np.float64) * np.sqrt(d["Ml"].to_numpy(np.float64))
           * d["Vt"].to_numpy(np.float64))
    W = np.empty(len(d))
    for c, pos in R.sightline_groups(codes):
        prof = G.density_profile(*keys[c])
        ds = Ds[pos]
        grid = np.linspace(ds.min(), ds.max(), 200)
        W[pos] = fac[pos] * prof.Nstart / nsim[c] * np.interp(ds, grid, G.lens_distance_norm(prof, grid))
    mm = U.mean_mass_from_y1(os.path.join(mean_mass_dir, name, "y1_yields.md"))
    y = T_S * R.RATE_UNIT * W / d["struc"].map(mm).to_numpy(float)
    lines = subprocess.run(["grep", "-E", r"^longtitude:|^ndd \(", os.path.join(directory, "run.log")],
                           capture_output=True, text=True, errors="replace", check=True).stdout.splitlines()
    cov, key, nl = {}, None, None
    for line in lines:
        if line.startswith("longtitude:"):
            p = line.split(); key = (round(float(p[1]), 3), round(float(p[3]), 3))
        elif line.startswith("ndd (LSST)") and key:
            nl = int(line.split()[2])
        elif line.startswith("ndd (Roman)") and key:
            cov[key] = (nl, int(line.split()[2])); key = None
    nvis = np.array([cov.get(k, (0, 0))[0] for k in keys])
    foot = np.array([cov.get(k, (0, 0))[1] > 0 for k in keys])
    vbins = np.array([0, 15, 50, 100, 200, 400, 800, 1600, 5000])
    return dict(lon=d["lon"].to_numpy(float), lat=d["lat"].to_numpy(float), Ds=Ds, y=y,
                codes=codes, keys=keys, vis_group=np.digitize(nvis, vbins)[codes],
                ft=foot[codes], detL=(d["detL"] == 1).to_numpy(),
                r=d["magb_r"].to_numpy(float))


def whole_scan(ctx, dust, variant, aks_av, pooled):
    """Rubin's whole-scan yield with the dust corrected, for one variant of the correction."""
    y, detL, ft = ctx["y"], ctx["detL"], ctx["ft"]
    groups = np.zeros(y.size, int) if pooled else ctx["vis_group"]
    dav = dust.delta_av(ctx["lon"], ctx["lat"], ctx["Ds"], variant, aks_av,
                        codes=ctx["codes"], keys=ctx["keys"])
    s = survival(ctx["r"], AR_AV * dav, detL.astype(float), y, groups)
    tot0 = float(y[detL].sum())
    res = dict(total=(float((y * s).sum()), float(np.sqrt(((y * s) ** 2).sum()))),
               outside=(float((y * s)[~ft].sum()), float(np.sqrt(((y * s)[~ft] ** 2).sum()))),
               inside=float((y * s)[ft].sum()), uncorrected=tot0,
               outside_uncorrected=float(y[detL & ~ft].sum()))
    return res


def dust_summary(dust, out_dir):
    """Model vs DECaPS vs near-infrared vs nominal reference on Roman's sightlines, with each map's
    own reliability flag at 8 kpc -- the evidence the report's dust table quotes."""
    import subprocess
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from dustmaps.bayestar import BayestarQuery
    lines = subprocess.run(["grep", "-E", r"^longtitude:|^ndd \(Roman", LOG_FOR_SUMMARY],
                           capture_output=True, text=True).stdout.splitlines()
    rom, key = [], None
    for ln in lines:
        if ln.startswith("longtitude:"):
            p = ln.split(); key = (float(p[1]), float(p[3]))
        elif key is not None:
            if int(ln.split()[2]) > 0:
                rom.append(key)
            key = None
    rom = np.array(rom)
    c8 = SkyCoord(l=rom[:, 0] * u.deg, b=rom[:, 1] * u.deg, distance=np.full(len(rom), 8.0) * u.kpc,
                  frame="galactic")
    _, fb = BayestarQuery(max_samples=1)(c8, mode="mean", return_flags=True)
    rows = []
    for (l, b), brel in zip(rom, fb["reliable_dist"]):
        dav, drel = dust.decaps_profile(l, b)
        i8 = int(np.argmin(np.abs(DGRID - 8.0)))
        rows.append(dict(l=l, b=b, block="GC field" if b > -0.8 else "five-field",
                         model=float(dust.model_av(l, b, np.array([8.0]))[0]),
                         decaps=float(np.interp(8.0, DGRID, dav)),
                         aks=float(np.interp(8.0, DGRID, dust.aks_profile(l, b))),
                         nir=float(np.interp(8.0, DGRID, dust.aks_profile(l, b))) / AKS_AV["nominal"],
                         nominal=float(np.interp(8.0, DGRID, dust.reference_profile(l, b, "hybrid",
                                                                                      AKS_AV["nominal"]))),
                         nir_from_kpc=dust.dsat.get((round(l, 3), round(b, 3), AKS_AV["nominal"]),
                                                    np.nan),
                         bayestar_reliable_8kpc=bool(brel), decaps_reliable_8kpc=bool(drel[i8])))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(out_dir, "u5_dust_summary.csv"), index=False)
    g = t.groupby("block").agg(n=("l", "size"), model=("model", "median"), decaps=("decaps", "median"),
                               nir=("nir", "median"), nominal=("nominal", "median"),
                               saturated_by_8kpc=("nir_from_kpc", lambda x: float((x <= 8.0).mean())),
                               bayestar_reliable=("bayestar_reliable_8kpc", "mean"),
                               decaps_reliable=("decaps_reliable_8kpc", "mean"))
    print("Dust at 8 kpc on Roman's sightlines (median A_V; reliable = fraction flagged reliable):\n"
          + g.round(3).to_string(), flush=True)
    return g


# ---------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", action="append", required=True, metavar="NAME=DIR")
    ap.add_argument("--mean-mass", required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--no-resolution", action="store_true")
    a = ap.parse_args()
    dust = Dust()
    global LOG_FOR_SUMMARY
    LOG_FOR_SUMMARY = os.path.join(a.run[0].split("=", 1)[1], "run.log")
    dust_summary(dust, a.out)
    for spec in a.run:
        name, directory = spec.split("=", 1)
        tag = R.load_provenance(os.path.join(directory, "files/MONTLMC/files/run_provenance.txt")
                                ).get("population_tag", "5")
        d = R.load_events(os.path.join(directory, f"test{tag}_foot.dat"), usecols=FOOT_COLS,
                          chunksize=500_000)
        mapf = os.path.join(directory, f"files/MONTLMC/files/MapLMC{tag}.dat")
        w, _ = R.attach_weight(d, mapf, [os.path.join(directory, "run.log")])
        d["W"] = w.to_numpy()
        mm = U.mean_mass_from_y1(os.path.join(a.mean_mass, name, "y1_yields.md"))
        y = T_S * R.RATE_UNIT * d["W"].to_numpy() / d["struc"].map(mm).to_numpy(float)
        check = float(y[(d.detR == 1).to_numpy()].sum())
        print(f"[{name}] footprint draws {len(d):,}; uncorrected Roman N_1 {check:.4g} "
              f"(must equal y1's footprint 'Roman detects')", flush=True)
        blocks = (d["lat"].to_numpy() >= -0.8).astype(int)
        lon, lat, ds = (d[c].to_numpy(float) for c in ("lon", "lat", "Ds"))
        variants = {k: (dust.delta_av(lon, lat, ds, ref, aks), np.zeros(len(d), int) if pooled else blocks)
                    for k, (ref, aks, pooled) in VARIANTS.items()}
        res = {k: footprint_variant(d, y, dav, g) for k, (dav, g) in variants.items()}
        zero = footprint_variant(d, y, np.zeros(len(d)), blocks)      # must reproduce as simulated
        if not a.no_resolution:
            for k in VARIANTS:
                if k == "pooled":
                    continue
                res[k].update(resolution_variant(d, y, variants[k][0], res[k]["_sR"]))
            zero.update(resolution_variant(d, y, np.zeros(len(d)), zero["_sR"]))
            res["pooled"].update({q: res["nominal"][q] for q in res["nominal"]
                                  if q.startswith(("P(res", "Roman resolves"))})
        # whole scan
        ws_path = os.path.join(directory, f"test{tag}_rubincols.dat")
        if os.path.exists(ws_path):
            ctx = whole_scan_context(name, directory, a.mean_mass)
            wsv = {k: whole_scan(ctx, dust, ref, aks, pooled)
                   for k, (ref, aks, pooled) in VARIANTS.items()}
            print(f"[{name}] whole-scan Rubin: uncorrected {wsv['nominal']['uncorrected']:.4g} "
                  f"(must equal y1), corrected {wsv['nominal']['total'][0]:.4g}", flush=True)
            del ctx
            for k, v in wsv.items():
                res[k]["Rubin detects, whole scan"] = v["total"]
                res[k]["Rubin detects, outside footprint"] = v["outside"]
                tot_in = res[k]["Roman detects"][0] + res[k]["Rubin only (footprint)"][0]
                # MC: the outside and footprint draws are independent, and Roman's detections and
                # Rubin-only ones are disjoint, so the three Poisson variances simply add.
                o, so = v["outside"]
                si2 = res[k]["Roman detects"][1] ** 2 + res[k]["Rubin only (footprint)"][1] ** 2
                mc = float(np.sqrt((tot_in * so) ** 2 + o ** 2 * si2) / (o + tot_in) ** 2)
                res[k]["share of detections outside the footprint [%]"] = (
                    100 * o / (o + tot_in), 100 * mc)
            zero["Rubin detects, whole scan"] = (wsv["nominal"]["uncorrected"], np.nan)
        # collect
        for q in res["nominal"]:
            if q.startswith("_"):
                continue
            v0, mc = res["nominal"][q]
            get = lambda k: res[k].get(q, (np.nan,))[0]
            # reference map: the near-infrared-only variant's shift (DECaPS-only is excluded by
            # VVV, Step U6); A_Ks/A_V: the larger shift of its two alternatives (the nominal is not
            # central between them); efficiency modelling: the pooled variant's shift
            s_map = abs(get("marshall") - v0)
            s_aks = max(abs(get("aks_hi") - v0), abs(get("aks_lo") - v0))
            s_eff = abs(get("pooled") - v0)
            syst = float(np.sqrt(np.nansum([s_map ** 2, s_aks ** 2, s_eff ** 2])))
            ROWS.append(dict(population=name, quantity=q, corrected=v0, mc=mc, syst=syst,
                             total=np.sqrt(np.nan_to_num(mc) ** 2 + syst ** 2),
                             as_simulated=zero.get(q, (np.nan,))[0],
                             s_map=s_map, s_aks=s_aks, s_eff=s_eff,
                             marshall=get("marshall"), decaps=get("decaps"), aks_lo=get("aks_lo"),
                             aks_hi=get("aks_hi"), pooled=get("pooled")))
        print(pd.DataFrame([r for r in ROWS if r["population"] == name]).to_string(index=False),
              flush=True)
        del d
    ds = {k[:2]: v for k, v in dust.dsat.items() if k[2] == AKS_AV["nominal"]}
    v = np.array(list(ds.values()))
    print(f"Nominal reference over all {v.size} sightlines met: switches to the near-infrared map "
          f"(DECaPS saturated) by 8 kpc on {(v <= 8).sum()}, by 20 kpc on {np.isfinite(v).sum()}",
          flush=True)
    t = pd.DataFrame(ROWS)
    t.to_csv(os.path.join(a.out, "u5_corrected_numbers.csv"), index=False)
    lines = ["# Headline numbers with the dust corrected (Step U5)", "",
             "Generated by `analysis/u5_corrected_numbers.py`. corrected = nominal correction; "
             "mc = Monte Carlo; syst = quadrature sum of: reference dust (shift of the "
             "near-infrared-only variant from the nominal DECaPS + near-infrared where DECaPS is "
             "saturated), A_Ks/A_V (larger shift of 0.0734 and 0.102 from 0.0805), efficiency per "
             "block vs pooled; total = quadrature sum; "
             "as_simulated = the same code with zero extra dust (must reproduce the report).", ""]
    for pop, g in t.groupby("population", sort=False):
        lines += [f"## {pop}", "", "| quantity | corrected | mc | syst | total | as simulated |",
                  "|---|---:|---:|---:|---:|---:|"]
        for _, r in g.iterrows():
            f = lambda x: "" if not np.isfinite(x) else f"{x:.4g}"
            lines.append(f"| {r.quantity} | {f(r.corrected)} | {f(r.mc)} | {f(r.syst)} | "
                         f"{f(r.total)} | {f(r.as_simulated)} |")
        lines.append("")
    with open(os.path.join(a.out, "u5_corrected_numbers.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
