#!/usr/bin/env python3
"""Step U5: the report's headline numbers with the dust corrected, and their uncertainties.

WHY THIS EXISTS. Step U4 (Deviation 61) found the model's extinction three to four times too low
within ~1 deg of the plane (optical Bayestar/DECaPS saturate; checked against the near-infrared
Marshall et al. 2006 map) and estimated the size of the effect. The user asked for the abstract
and summary to report the CORRECTED results, with uncertainties. This script produces them, for
every number the abstract and summary quote, each with a Monte Carlo error and a systematic error
of the correction itself.

THE CORRECTION. For every simulated draw i, dA_V,i = A_V,Marshall(l, b, Ds_i) - A_V,model(l, b, Ds_i)
at the source's OWN distance (Step U4 used 8 kpc for every source). The extra dust only dims the
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
           A_Ks/A_V = 0.10 vs 0.114 (nominal 0.11; half the difference),
           efficiency measured per field block (nominal) vs pooled over the footprint,
           per-draw distance (nominal) vs every source at 8 kpc.
What it does NOT cover: that a surviving event's PRECISION also degrades when its source dims. For
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
AKS_AV = {"nominal": 0.11, "aks_lo": 0.10, "aks_hi": 0.114}
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


# ---------------------------------------------------------------------------------------------
# Dust
# ---------------------------------------------------------------------------------------------
class Dust:
    """Model (the simulator's extinction tables) and near-infrared (Marshall) A_V(d) profiles."""

    def __init__(self):
        from dustmaps.config import config
        config["data_dir"] = "dustmaps"
        from dustmaps.marshall import MarshallQuery
        self.mq = MarshallQuery()
        files = sorted(glob.glob("files/ext/bayestar_*.txt"))
        self.files = files
        self.pos = np.array([[float(x) for x in re.findall(r"bayestar_(-?[\d.]+)_(-?[\d.]+)\.txt", f)[0]]
                             for f in files])
        self._model, self._aks = {}, {}

    def model_profile(self, l, b):
        """A_V on DGRID from the table nearestSightline() would pick."""
        i = int(np.argmin((self.pos[:, 0] - l) ** 2 + (self.pos[:, 1] - b) ** 2))
        if i not in self._model:
            a = np.loadtxt(self.files[i])
            if not np.isfinite(a[:, 3]).any():
                # A "total dropout" table (78 of 2,518, all just north of dec -30 where Bayestar
                # has no data): every A_V is NaN. The simulator reads it with >>, the stream fails
                # on the first "nan", and the sightline is simulated with ~zero dust (measured:
                # median A_r 0.0002 mag on its detections). So the model's dust here IS zero.
                self._model[i] = np.zeros(DGRID.size)
            else:
                self._model[i] = np.interp(DGRID, a[:, 2], a[:, 3])
        return self._model[i]

    def aks_profile(self, l, b):
        """Marshall A_Ks on DGRID; NaN beyond its coverage is held at the last valid value."""
        k = (round(l, 3), round(b, 3))
        if k not in self._aks:
            import astropy.units as u
            from astropy.coordinates import SkyCoord
            c = SkyCoord(l=np.full(DGRID.size, l) * u.deg, b=np.full(DGRID.size, b) * u.deg,
                         distance=DGRID * u.kpc, frame="galactic")
            a = np.asarray(self.mq(c), float)
            ok = np.isfinite(a)
            if ok.sum() >= 2:
                a = np.interp(DGRID, DGRID[ok], a[ok])
            elif ok.sum() == 1:
                a = np.full(DGRID.size, a[ok][0])
            else:
                a = np.full(DGRID.size, np.nan)
            self._aks[k] = a
        return self._aks[k]

    def delta_av(self, lon, lat, ds, aks_av, fixed_d=None, codes=None, keys=None):
        """Per-draw A_V,Marshall - A_V,model; 0 where Marshall has no value. Pass the
        sightline codes/keys (romanlib.sightline_index) for large tables."""
        lon, lat, ds = map(np.asarray, (lon, lat, ds))
        out = np.zeros(lon.size)
        if codes is None:
            s = pd.Series(list(zip(np.round(lon, 3), np.round(lat, 3))))
            groups = [(k, np.asarray(p)) for k, p in s.groupby(s).groups.items()]
        else:
            groups = [(keys[c], pos) for c, pos in R.sightline_groups(codes)]
        for k, pos in groups:
            d = np.full(pos.size, fixed_d) if fixed_d else ds[pos]
            mod = np.interp(d, DGRID, self.model_profile(*k))
            aks = self.aks_profile(*k)
            if np.all(np.isnan(aks)):
                continue
            out[pos] = np.interp(d, DGRID, aks) / aks_av - mod
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
        out[f"P(resolvable) Roman D={bar} [%]"] = (100 * float(w[res].sum() / w.sum()), np.nan)
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


def whole_scan(ctx, dust, aks_av, pooled, fixed_d):
    """Rubin's whole-scan yield with the dust corrected, for one variant of the correction."""
    y, detL, ft = ctx["y"], ctx["detL"], ctx["ft"]
    groups = np.zeros(y.size, int) if pooled else ctx["vis_group"]
    dav = dust.delta_av(ctx["lon"], ctx["lat"], ctx["Ds"], aks_av, fixed_d,
                        codes=ctx["codes"], keys=ctx["keys"])
    s = survival(ctx["r"], AR_AV * dav, detL.astype(float), y, groups)
    tot0 = float(y[detL].sum())
    res = dict(total=(float((y * s).sum()), float(np.sqrt(((y * s) ** 2).sum()))),
               outside=(float((y * s)[~ft].sum()), float(np.sqrt(((y * s)[~ft] ** 2).sum()))),
               inside=float((y * s)[ft].sum()), uncorrected=tot0,
               outside_uncorrected=float(y[detL & ~ft].sum()))
    return res


# ---------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", action="append", required=True, metavar="NAME=DIR")
    ap.add_argument("--mean-mass", required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--no-resolution", action="store_true")
    a = ap.parse_args()
    dust = Dust()
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
        variants = {
            "nominal": (dust.delta_av(lon, lat, ds, AKS_AV["nominal"]), blocks),
            "aks_lo": (dust.delta_av(lon, lat, ds, AKS_AV["aks_lo"]), blocks),
            "aks_hi": (dust.delta_av(lon, lat, ds, AKS_AV["aks_hi"]), blocks),
            "pooled": (dust.delta_av(lon, lat, ds, AKS_AV["nominal"]), np.zeros(len(d), int)),
            "d8kpc": (dust.delta_av(lon, lat, ds, AKS_AV["nominal"], fixed_d=8.0), blocks),
        }
        res = {k: footprint_variant(d, y, dav, g) for k, (dav, g) in variants.items()}
        zero = footprint_variant(d, y, np.zeros(len(d)), blocks)      # must reproduce as simulated
        if not a.no_resolution:
            for k in ("nominal", "aks_lo", "aks_hi", "d8kpc"):
                res[k].update(resolution_variant(d, y, variants[k][0], res[k]["_sR"]))
            zero.update(resolution_variant(d, y, np.zeros(len(d)), zero["_sR"]))
            for k in ("pooled",):
                res[k].update({q: res["nominal"][q] for q in res["nominal"] if q.startswith(("P(res", "Roman resolves"))})
        # whole scan
        ws_path = os.path.join(directory, f"test{tag}_rubincols.dat")
        if os.path.exists(ws_path):
            ctx = whole_scan_context(name, directory, a.mean_mass)
            wsv = {"nominal": whole_scan(ctx, dust, 0.11, False, None),
                   "aks_lo": whole_scan(ctx, dust, 0.10, False, None),
                   "aks_hi": whole_scan(ctx, dust, 0.114, False, None),
                   "pooled": whole_scan(ctx, dust, 0.11, True, None),
                   "d8kpc": whole_scan(ctx, dust, 0.11, False, 8.0)}
            print(f"[{name}] whole-scan Rubin: uncorrected {wsv['nominal']['uncorrected']:.4g} "
                  f"(must equal y1), corrected {wsv['nominal']['total'][0]:.4g}", flush=True)
            del ctx
            for k, v in wsv.items():
                res[k]["Rubin detects, whole scan"] = v["total"]
                res[k]["Rubin detects, outside footprint"] = v["outside"]
                tot_in = res[k]["Roman detects"][0] + res[k]["Rubin only (footprint)"][0]
                res[k]["share of detections outside the footprint [%]"] = (
                    100 * v["outside"][0] / (v["outside"][0] + tot_in), np.nan)
            zero["Rubin detects, whole scan"] = (wsv["nominal"]["uncorrected"], np.nan)
        # collect
        for q in res["nominal"]:
            if q.startswith("_"):
                continue
            v0, mc = res["nominal"][q]
            get = lambda k: res[k].get(q, (np.nan,))[0]
            syst = np.sqrt(((get("aks_hi") - get("aks_lo")) / 2) ** 2 + (get("pooled") - v0) ** 2
                           + (get("d8kpc") - v0) ** 2)
            ROWS.append(dict(population=name, quantity=q, corrected=v0, mc=mc, syst=syst,
                             total=np.sqrt(np.nan_to_num(mc) ** 2 + syst ** 2),
                             as_simulated=zero.get(q, (np.nan,))[0],
                             aks_lo=get("aks_lo"), aks_hi=get("aks_hi"), pooled=get("pooled"),
                             d8kpc=get("d8kpc")))
        print(pd.DataFrame([r for r in ROWS if r["population"] == name]).to_string(index=False),
              flush=True)
        del d
    t = pd.DataFrame(ROWS)
    t.to_csv(os.path.join(a.out, "u5_corrected_numbers.csv"), index=False)
    lines = ["# Headline numbers with the dust corrected (Step U5)", "",
             "Generated by `analysis/u5_corrected_numbers.py`. corrected = nominal correction; "
             "mc = Monte Carlo; syst = spread of the correction (A_Ks/A_V 0.10-0.114, efficiency "
             "per block vs pooled, per-draw distance vs 8 kpc); total = quadrature sum; "
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
