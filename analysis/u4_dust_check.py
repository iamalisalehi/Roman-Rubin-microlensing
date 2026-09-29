#!/usr/bin/env python3
"""Step U4: the model's dust against an independent near-infrared map, and what it does to yields.

WHY THIS EXISTS. Checking the footprint-mismatch estimate (Step U3) showed the model's extinction
toward the inner footprint behaving backwards: the median A_r of footprint sources DROPS toward
the plane (4.4 mag at b = -1.24, 2.3 at b = -0.94) and the Galactic-centre field has A_V ~ 3 to
8 kpc. The dust comes from files/ext/ (maps.py): Bayestar19 north of dec = -30 deg, DECaPS south
of it. Bayestar is built from Pan-STARRS OPTICAL photometry; near the plane toward the Galactic
centre its stars cannot see through the dust lanes, its profiles saturate after a few kpc, and the
far values are lower limits. The Marshall et al. (2006) 3D map, built from 2MASS near-infrared
star counts for the inner Galaxy, is the independent check (dustmaps' MarshallQuery; A_Ks, turned
into A_V with A_Ks/A_V = 0.11, the CCM89 value for R_V 2.5-3.1, 0.10-0.114).

PART 1 -- the whole scan. Per sightline, model A_V to 8 kpc (the nearest extinction file, exactly as
nearestSightline() picks it) against Marshall's, and the share of each survey's yield that comes
from sightlines where the model has less than half, or more than 1.25x, Marshall's dust.

PART 2 -- the footprint, corrected. Every footprint draw (detected or not) is re-used: the dust
changes only how bright the source is, so the expected number of detections on a sightline with
Delta A extra magnitudes of extinction is sum_i y_i eff(m_i + Delta A), where eff(m) is the
yield-weighted detection efficiency against the blended baseline magnitude, measured from the same
draws (F146 for Roman, r for Rubin). The correction factor per sightline is that sum over the same
sum with Delta A = 0 (so the smoothing in eff cancels). Delta A = (A_V,Marshall - A_V,model) times
A_F146/A_V = 0.197 or A_r/A_V = 0.854 (R_V = 2.5, the bulge value). Then (a) the corrected yield
of the SIMULATED footprint, and (b) Step U3's real/simulated footprint ratio recomputed on the
corrected row densities.

ASSUMPTIONS, stated. eff(m) is taken as the same everywhere in the footprint (the timescale, impact
parameter and blending distributions do not depend on the dust); sources are taken at 8 kpc, so the
correction ignores that foreground disc sources see less of the extra dust (it overcorrects them);
A_Ks/A_V is uncertain by ~10%. The result is an estimate of size and direction, not a replacement
for a run with better dust.

INPUTS. Footprint draws extracted from the full tables (ndw_R > 0, column 37):
    awk 'NR==1 || /^#/ {print; next} $37>0' test<tag>.dat > test<tag>_foot.dat
and the detection-only tables of Step U1 for Part 1.

    .roman/bin/python analysis/u4_dust_check.py --run bulge=runs/prod_bulge_20260924 \\
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
import romanlib as R                  # noqa: E402
import u1_report_numbers as U         # noqa: E402
import u3_footprint_offset as U3      # noqa: E402

AK_AV = 0.11
AF146_AV, AR_AV = 0.197, 0.854       # CCM89, R_V = 2.5 (extinctiontest)
DIST = 8.0
FOOT_COLS = ["lon", "lat", "w_area", "Ml", "Vt", "Ds", "struc", "detL", "detR", "detJ",
             "magb_F146", "magb_r", "ndw_R"]


def ext_files():
    files = sorted(glob.glob("files/ext/bayestar_*.txt"))
    pos = np.array([[float(x) for x in re.findall(r"bayestar_(-?[\d.]+)_(-?[\d.]+)\.txt", f)[0]]
                    for f in files])
    return files, pos


class Dust:
    def __init__(self):
        from dustmaps.config import config
        config["data_dir"] = "dustmaps"
        from dustmaps.marshall import MarshallQuery
        self.mq = MarshallQuery()
        self.files, self.pos = ext_files()
        self.cache = {}

    def model_av(self, l, b):
        i = int(np.argmin((self.pos[:, 0] - l) ** 2 + (self.pos[:, 1] - b) ** 2))
        if i not in self.cache:
            a = np.loadtxt(self.files[i])
            self.cache[i] = float(np.interp(DIST, a[:, 2], a[:, 3]))
        return self.cache[i]

    def marshall_av(self, l, b):
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        c = SkyCoord(l=np.atleast_1d(l) * u.deg, b=np.atleast_1d(b) * u.deg,
                     distance=np.full(np.size(l), DIST) * u.kpc, frame="galactic")
        return np.asarray(self.mq(c), float) / AK_AV


def part1(dust, runs, out):
    """Model vs Marshall per sightline, weighted by each survey's yield."""
    keys = set()
    for run in runs:
        keys |= set(zip(run.df["lon"].round(3), run.df["lat"].round(3)))
    keys = sorted(keys)
    L = np.array([k[0] for k in keys]); B = np.array([k[1] for k in keys])
    mod = np.array([dust.model_av(l, b) for l, b in keys])
    mar = dust.marshall_av(L, B)
    ratio = pd.Series(mod / np.where(mar > 0, mar, np.nan), index=pd.MultiIndex.from_tuples(keys))
    print(f"Part 1: {len(keys)} sightlines with detections; model/Marshall A_V(8 kpc): median "
          f"{np.nanmedian(ratio):.2f}", flush=True)
    for lo, hi in ((0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.5), (2.5, 6)):
        m = (np.abs(B) >= lo) & (np.abs(B) < hi)
        if m.any():
            out.append(dict(part="1", population="scan", quantity=f"model/Marshall A_V, |b| {lo}-{hi}",
                            value=float(np.nanmedian(mod[m] / mar[m])), n=int(m.sum())))
    for run in runs:
        d = run.df
        key = pd.Series(list(zip(d["lon"].round(3), d["lat"].round(3))))
        r_ev = key.map(ratio.to_dict()).to_numpy(float)
        for surv, col in (("Rubin", "detL"), ("Roman", "detR")):
            sel = (d[col] == 1).to_numpy()
            y = d["y"].to_numpy()
            tot = y[sel].sum()
            for lab, m in (("model < 0.5x Marshall", r_ev < 0.5),
                           ("model 0.5-0.8x", (r_ev >= 0.5) & (r_ev < 0.8)),
                           ("model 0.8-1.25x", (r_ev >= 0.8) & (r_ev <= 1.25)),
                           ("model > 1.25x", r_ev > 1.25)):
                out.append(dict(part="1", population=run.name,
                                quantity=f"share of {surv} detected yield, {lab}",
                                value=float(y[sel & m].sum() / tot), n=int((sel & m).sum())))


def efficiency(mag, det, y, bins):
    """Yield-weighted detection efficiency in magnitude bins, monotone-cleaned at the faint end."""
    idx = np.clip(np.digitize(mag, bins) - 1, 0, len(bins) - 2)
    num = np.bincount(idx, weights=y * det, minlength=len(bins) - 1)
    den = np.bincount(idx, weights=y, minlength=len(bins) - 1)
    eff = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
    centres = 0.5 * (bins[1:] + bins[:-1])
    ok = np.isfinite(eff)
    return lambda m: np.interp(m, centres[ok], eff[ok], left=eff[ok][0], right=0.0)


def part2(dust, name, directory, mean_mass_dir, realpix, out):
    tag = R.load_provenance(os.path.join(directory, "files/MONTLMC/files/run_provenance.txt")
                            ).get("population_tag", "5")
    path = os.path.join(directory, f"test{tag}_foot.dat")
    d = R.load_events(path, usecols=FOOT_COLS, chunksize=500_000)
    mapf = os.path.join(directory, f"files/MONTLMC/files/MapLMC{tag}.dat")
    w, _ = R.attach_weight(d, mapf, [os.path.join(directory, "run.log")])
    mm = U.mean_mass_from_y1(os.path.join(mean_mass_dir, name, "y1_yields.md"))
    T = (3652.43 - 2 * R.T0_MARGIN_DAYS) * 86400.0
    y = T * R.RATE_UNIT * w.to_numpy() / d["struc"].map(mm).to_numpy(float)
    key = pd.Series(list(zip(d["lon"].round(3), d["lat"].round(3))))
    sl = sorted(set(key))
    dav = {k: dust.marshall_av(k[0], k[1])[0] - dust.model_av(*k) for k in sl}
    dA = key.map(dav).to_numpy(float)
    for surv, mcol, dcol, coef in (("Roman", "magb_F146", "detR", AF146_AV),
                                   ("Rubin", "magb_r", "detL", AR_AV)):
        mag = d[mcol].to_numpy(float)
        det = (d[dcol] == 1).to_numpy().astype(float)
        eff = efficiency(mag, det, y, np.arange(12.0, 32.01, 0.25))
        e0, e1 = eff(mag), eff(mag + coef * dA)
        g = pd.DataFrame({"key": key, "lat": d["lat"].round(2), "y_det": y * det,
                          "y2": (y * det) ** 2, "e0": y * e0, "e1": y * e1})
        per = g.groupby("key").agg(lat=("lat", "first"), y_det=("y_det", "sum"), y2=("y2", "sum"),
                                   e0=("e0", "sum"), e1=("e1", "sum"))
        per["c"] = np.where(per.e0 > 0, per.e1 / per.e0, 1.0)
        per["y_corr"] = per.y_det * per.c
        tot0, tot1 = per.y_det.sum(), per.y_corr.sum()
        out.append(dict(part="2", population=name,
                        quantity=f"{surv} detections, simulated footprint: dust-corrected / model",
                        value=tot1 / tot0, n=len(per)))
        for bname, bm in (("five-field", per.lat < U3.BLOCK_SPLIT), ("GC field", per.lat >= U3.BLOCK_SPLIT)):
            out.append(dict(part="2", population=name,
                            quantity=f"{surv} detections, {bname}: dust-corrected / model",
                            value=per.y_corr[bm].sum() / per.y_det[bm].sum(), n=int(bm.sum())))
        per["v1"] = per.y2 * per.c ** 2
        rows = per.groupby("lat").agg(n=("y_det", "size"), y0=("y_det", "sum"), y1=("y_corr", "sum"),
                                      v1=("v1", "sum"), c=("c", "median"))
        rows["area"] = rows.n * U3.FINE ** 2
        print(f"[{name}] {surv}: per-row yield density, model -> dust-corrected (median factor)\n"
              + (rows.assign(rho0=rows.y0 / rows.area, rho1=rows.y1 / rows.area)
                 [["n", "rho0", "rho1", "c"]].round(3).to_string()), flush=True)
        # U3's real/simulated ratio on the corrected densities
        bpix, lpix, wpix, dApix = realpix[:4]
        cells = np.array([k for k in per.index])
        num = den = 0.0
        for bname, bm_row, bm_pix, bm_cell in (
                ("five-field", rows.index < U3.BLOCK_SPLIT, bpix < U3.BLOCK_SPLIT,
                 cells[:, 1] < U3.BLOCK_SPLIT),
                ("GC field", rows.index >= U3.BLOCK_SPLIT, bpix >= U3.BLOCK_SPLIT,
                 cells[:, 1] >= U3.BLOCK_SPLIT)):
            r = rows[bm_row]
            t = pd.DataFrame({"b": r.index.to_numpy(float), "area": r.area.to_numpy(),
                              "rho": (r.y1 / r.area).to_numpy(),
                              "err": (np.sqrt(r.v1) / r.area).to_numpy()})
            blk = r.y1.sum()
            if (t.rho > 0).sum() < 3 or blk <= 0:
                # A block the corrected dust empties (Rubin in the Galactic-centre field): it adds
                # nothing on either layout, so it drops out of the ratio.
                out.append(dict(part="2", population=name,
                                quantity=f"{surv}: dust-corrected slope d ln(rho)/db, {bname}",
                                value=np.nan, n=len(r)))
                continue
            a, s, cov, chi2, dof = U3.fit_block(t)
            rho = lambda bb: np.exp(a + s * bb)
            n_sim_fit = (rho(cells[bm_cell, 1]) * U3.FINE ** 2).sum()
            n_real_fit = (wpix[bm_pix] * rho(bpix[bm_pix])).sum() * dApix
            num += blk * n_real_fit / n_sim_fit
            den += blk
            out.append(dict(part="2", population=name,
                            quantity=f"{surv}: dust-corrected slope d ln(rho)/db, {bname}",
                            value=s, n=len(r), err=float(np.sqrt(cov[1, 1]))))
        out.append(dict(part="2", population=name,
                        quantity=f"{surv} detections: real/simulated footprint, dust-corrected",
                        value=num / den, n=len(per)))


def log_coverage(path):
    """{(lon, lat): (Rubin visits, Roman epochs)} for every sightline the run entered."""
    import subprocess
    lines = subprocess.run(["grep", "-E", r"^longtitude:|^ndd \(", path], capture_output=True,
                           text=True, errors="replace", check=True).stdout.splitlines()
    out, key, nl = {}, None, None
    for line in lines:
        if line.startswith("longtitude:"):
            p = line.split()
            key = (round(float(p[1]), 3), round(float(p[3]), 3))
        elif line.startswith("ndd (LSST)") and key:
            nl = int(line.split()[2])
        elif line.startswith("ndd (Roman)") and key:
            out[key] = (nl, int(line.split()[2]))
            key = None
    return out


def part3(dust, name, directory, path, mean_mass_dir, roman_factor, out):
    """Whole-scan Rubin yield with the dust corrected, efficiency stratified by Rubin visits.

    Needs every draw (not only detections), for the efficiency curves and the weights: a column
    extract of the full table (lon lat w_area Ml Vt Ds struc magb_r detL detJ). The weight is
    rebuilt as in u1's intrinsic check (nsim = rows per sightline, Z interpolated), and the
    uncorrected Rubin yield it gives is checked against y1's before anything is corrected.
    """
    import galaxy_model as G
    d = R.load_events(path, chunksize=1_000_000, narrow=True)
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
    T = (3652.43 - 2 * R.T0_MARGIN_DAYS) * 86400.0
    y = T * R.RATE_UNIT * W / d["struc"].map(mm).to_numpy(float)
    detL = (d["detL"] == 1).to_numpy()
    total0 = float(y[detL].sum())
    print(f"[{name}] whole-scan Rubin yield rebuilt from all draws: {total0:.4g} (y1: see "
          f"y1_yields.csv)", flush=True)
    cov = log_coverage(os.path.join(directory, "run.log"))
    nvis = np.array([cov.get(k, (0, 0))[0] for k in keys])
    foot = np.array([cov.get(k, (0, 0))[1] > 0 for k in keys])
    L = np.array([k[0] for k in keys]); B = np.array([k[1] for k in keys])
    mar = dust.marshall_av(L, B)
    mod = np.array([dust.model_av(*k) for k in keys])
    dav = mar - mod
    nbad = int(np.sum(~np.isfinite(dav)))
    dav = np.where(np.isfinite(dav), dav, 0.0)
    dAr = AR_AV * dav[codes]
    mag = d["magb_r"].to_numpy(float)
    vbins = np.array([0, 15, 50, 100, 200, 400, 800, 1600, 5000])
    cls = np.digitize(nvis, vbins)[codes]
    e0 = np.zeros(len(d)); e1 = np.zeros(len(d))
    for k in np.unique(cls):
        m = cls == k
        eff = efficiency(mag[m], detL[m].astype(float), y[m], np.arange(12.0, 32.01, 0.25))
        e0[m], e1[m] = eff(mag[m]), eff(mag[m] + dAr[m])
    g = pd.DataFrame({"c": codes, "yd": y * detL, "e0": y * e0, "e1": y * e1})
    per = g.groupby("c").sum()
    per["corr"] = np.where(per.e0 > 0, per.e1 / per.e0, 1.0)
    per["foot"] = foot[per.index]
    Y0 = per.yd.sum(); Y1 = (per.yd * per["corr"]).sum()
    out0 = per.yd[~per.foot].sum(); out1 = (per.yd * per["corr"])[~per.foot].sum()
    in0 = per.yd[per.foot].sum(); in1 = (per.yd * per["corr"])[per.foot].sum()
    print(f"[{name}] Rubin whole scan: {Y0:.4g} -> {Y1:.4g} ({Y1/Y0:.3f}); outside the footprint "
          f"{out0:.4g} -> {out1:.4g} ({out1/out0:.3f}); inside {in0:.4g} -> {in1:.4g} ({in1/in0:.3f}); "
          f"{nbad} sightlines without a Marshall value (left uncorrected)", flush=True)
    for q, v in (("Rubin whole-scan yield, rebuilt (uncorrected)", Y0),
                 ("Rubin whole-scan yield: dust-corrected / model", Y1 / Y0),
                 ("Rubin outside the footprint: dust-corrected / model", out1 / out0),
                 ("Rubin inside the footprint (whole-scan eff.): dust-corrected / model", in1 / in0)):
        out.append(dict(part="3", population=name, quantity=q, value=v, n=len(per)))
    # share of all detections outside the footprint: Rubin-only outside + everything inside
    # (inside, approximated by Roman's corrected yield plus the Rubin-only part, which is small).
    out.append(dict(part="3", population=name,
                    quantity="note: outside-footprint detections are Rubin's by construction",
                    value=np.nan))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", action="append", required=True, metavar="NAME=DIR")
    ap.add_argument("--mean-mass", required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--whole-scan", action="append", default=[], metavar="NAME=FILE",
                    help="all-draws column extract (lon lat w_area Ml Vt Ds struc magb_r detL detJ)")
    a = ap.parse_args()
    dust = Dust()
    out, runs = [], []
    specs = [s.split("=", 1) for s in a.run]
    for name, directory in specs:
        runs.append(U.Run(name, directory, a.mean_mass))
    part1(dust, runs, out)
    del runs
    realpix = U3.real_pixels()
    for name, directory in specs:
        part2(dust, name, directory, a.mean_mass, realpix, out)
    dirs = dict(specs)
    for spec in a.whole_scan:
        name, path = spec.split("=", 1)
        part3(dust, name, dirs[name], path, a.mean_mass, None, out)
    t = pd.DataFrame(out)
    t.to_csv(os.path.join(a.out, "u4_dust_check.csv"), index=False)
    with pd.option_context("display.width", 200, "display.max_colwidth", 90, "display.max_rows", 200):
        print(t.to_string(index=False))


if __name__ == "__main__":
    main()
