#!/usr/bin/env python3
"""Step U3: what the simulated-vs-real Roman footprint mismatch, and the scan's corner cut, cost.

LEGACY ANALYSIS (Deviation 69). This script measures the runs made on the notional layout
(layout_40395 circles, box scan with a corner cut) against the real tiles; its FIELDS and scan
constants are that layout's on purpose. From Deviation 69 the simulator uses the adopted layout
and a distance-rule scan, so the mismatch it estimates no longer exists in new runs; the shared
geometry for those is analysis/gbtds_geometry.py.

WHY THIS EXISTS. Deviation 59 measured the simulated GBTDS fields against the real tile layout
(Whitepaper/roman_967_{spring,autumn,both}_aladinX.png): the modelled five-field block sits at
b = -1.2 where the real one is centred at -1.40, the Galactic-centre field ~0.1 deg closer to the
plane, and the real fields cover more sky. No run was made on the real layout, so this script
estimates what it would change from the runs that exist, and says how far the estimate can be
trusted. Two parts.

PART A -- Roman's footprint. Inside the simulated footprint the yield per deg^2 is measured row by
row in latitude (0.1 deg rows: six in the five-field block, seven in the Galactic-centre field) and
fitted as ln(rho) = a + s*b per block. The real layout is taken from the spring and autumn tile
images separately, each placed in (l, b) from its own Galactic grid (make_footprints.py's
calibration); each real pixel is weighted by the fraction of seasons it is observed (seasons
alternate spring/autumn, 5 + 5, with 3 + 3 high-cadence), which is what "the real footprint" means
for a ten-season yield. Then

    N_real / N_sim = sum_real_pixels w rho(b) dA / sum_sim_cells rho(b) dA

split into an AREA factor (rho constant) and a LATITUDE factor (the rest). Assumptions, stated:
rho depends on b only within a block (l is similar: real -0.60..1.58, simulated -0.72..1.52); it is
extrapolated below b = -1.5, where no simulated sightline is; the seasonal weight treats a patch seen
in half the seasons as half the yield. The fit's slope error is propagated by sampling.
Per-event quantities (fractions characterised, P(resolvable)) are also regressed on b: if they are
flat, the mismatch moves counts but not the conclusions drawn from fractions and ratios.

PART B -- the corner cut. The scan drops l < -0.9447 and b > 0.31 (Bulge.h lx, bx). For each 0.1-deg
cell of that corner and of the opposite upper corner, count the Rubin visits that image it (pointing
centre within 1.75 deg, as matchVisibleEpochs does) and, of those, the ones whose field also
overlaps a Roman field (pointing centre within 1.75 + 0.3003 deg of a Roman field centre). The
visit counter is first checked against the run log's own count on every simulated sightline. The
Rubin yield lost is then estimated from the simulated northern sightlines' yield density as a
function of Rubin visit count.

    .roman/bin/python analysis/u3_footprint_offset.py --run bulge=runs/prod_bulge_20260924 \\
        --run bh=runs/prod_bh_20260924 --run ns=runs/prod_ns_20260924 \\
        --mean-mass figures/yield_20260925 -o figures/u1_20260929
"""

import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R                 # noqa: E402
import u1_report_numbers as U        # noqa: E402

FIELDS = [(-0.417948, -1.2), (-0.008974, -1.2), (0.4, -1.2), (0.808974, -1.2),
          (1.217948, -1.2), (0.0, -0.125)]
FOV_ROMAN, FOV_RUBIN = 0.3003, 1.75
WID = 1.75
LON_MIN, LON_MAX = -0.219 - 0.2 - WID - WID, 1.4134 + 0.2 + WID + WID
LAT_MIN, LAT_MAX = -1.64 - 0.2 - WID - WID, -0.85 + 0.2 + WID + WID
LX, BX = 1.0053 - 0.2 - WID, -1.64 + 0.2 + WID
FINE = 0.1
BLOCK_SPLIT = -0.8          # b above: the Galactic-centre field; below: the five-field block
IMG = "Whitepaper/roman_967_{}_aladinX.png"
LOG = "runs/prod_bulge_20260924/run.log"
N_DRAW = 4000
OGLE_SLOPE = (0.0, 0.0, 0)


# ---------------------------------------------------------------------------------------------
# The real tiles, per season, in (l, b)
# ---------------------------------------------------------------------------------------------
def season_tiles(season):
    img = np.asarray(Image.open(IMG.format(season)).convert("RGB")).astype(int)
    R_, G, B = img[..., 0], img[..., 1], img[..., 2]
    green = (G > R_ + 20) & (G > B + 20)
    H, W = green.shape
    xs = [x for x in range(W) if green[:, x].sum() > 0.9 * H]
    ys = [y for y in range(H) if green[y, :].sum() > 0.9 * W]
    ppd = (np.diff(xs).mean() + np.diff(ys).mean()) / 2 / 0.5
    x0 = min(xs, key=lambda x: abs(x - 615))       # l = 0 ('000' label)
    y0 = min(ys, key=lambda y: abs(y - 124))       # b = 0 ('+00' label)
    edge = (G > 170) & (R_ < 190) & (B < 150) & (G > R_ + 40)
    edge[:45, :] = edge[790:, :] = False
    edge[:, :45] = edge[:, 935:] = False
    # Dilate to close the outlines, fill, then erode the dilation back off: without the erosion
    # every detector comes out ~9% too large (a 1-px rim on a ~45-px square).
    tiles = ndimage.binary_erosion(
        ndimage.binary_fill_holes(ndimage.binary_dilation(edge, iterations=1)), iterations=1)
    yy, xx = np.mgrid[0:H, 0:W]
    return tiles, -(xx - x0) / ppd, -(yy - y0) / ppd, ppd


# ---------------------------------------------------------------------------------------------
# Part A
# ---------------------------------------------------------------------------------------------
def row_table(d, sel, w_col):
    """Per 0.1-deg latitude row of the footprint: area, yield density and its error."""
    lat = d["lat"].round(2).to_numpy()
    key = pd.Series(list(zip(d["lon"].round(3), d["lat"].round(3))))
    foot = d["foot"].to_numpy()
    rows = []
    for b in np.unique(lat[foot]):
        m = foot & (lat == b)
        nsl = key[m].nunique()
        y = d[w_col].to_numpy()[m & sel]
        area = nsl * FINE ** 2
        rows.append(dict(b=b, n_sl=nsl, area=area, rho=y.sum() / area,
                         err=np.sqrt((y ** 2).sum()) / area))
    return pd.DataFrame(rows)


def fit_block(t):
    """Weighted fit ln(rho) = a + s b; returns (a, s, cov)."""
    t = t[t.rho > 0]
    x, y = t.b.to_numpy(), np.log(t.rho.to_numpy())
    sig = (t.err / t.rho).to_numpy()
    A = np.vstack([np.ones_like(x), x]).T / sig[:, None]
    coef, *_ = np.linalg.lstsq(A, y / sig, rcond=None)
    cov = np.linalg.inv(A.T @ A)
    chi2 = float(np.sum(((y - coef[0] - coef[1] * x) / sig) ** 2))
    dof = max(len(x) - 2, 1)
    cov = cov * max(chi2 / dof, 1.0)       # inflate if the rows scatter more than Poisson
    return coef[0], coef[1], cov, chi2, dof


def real_pixels():
    """(b, l, weight, pixel area) of the real tiles, weight = fraction of seasons observed."""
    sp, L, B, ppd = season_tiles("spring")
    au, L2, B2, ppd2 = season_tiles("autumn")
    assert np.allclose(L, L2) and np.allclose(B, B2)
    w = 0.5 * sp + 0.5 * au
    m = w > 0
    return B[m], L[m], w[m], 1.0 / ppd ** 2, sp.sum() / ppd ** 2, au.sum() / ppd ** 2


def ogle_slope_correction(path="figures/yield_20260925/bulge/y1_yields.md", lo=-1.95, hi=-0.75):
    """d ln(model/OGLE rate per star) / db over the footprint's latitudes, from y1's OGLE check.

    The model's rate is known to be too flat in latitude (OPEN_ITEMS). If the real rate per star
    rises toward the plane faster than the model's by this slope, the real yield-density slope is
    the model's MINUS it (the model/OGLE ratio falls toward the plane, so the slope is negative and
    the correction steepens the gradient).
    """
    b, r = [], []
    in_tab = False
    for line in open(path):
        if line.startswith("| b | draws |"):
            in_tab = True
            continue
        if in_tab:
            if not line.startswith("|"):
                break
            c = [x.strip() for x in line.strip("|\n").split("|")]
            if c[0].startswith("---"):
                continue
            bb = float(c[0])
            if lo <= bb <= hi:
                b.append(bb); r.append(float(c[4]))
    b, y = np.array(b), np.log(np.array(r))
    A = np.vstack([np.ones_like(b), b]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    sig2 = float(np.sum((y - A @ coef) ** 2)) / max(len(b) - 2, 1)
    err = np.sqrt(sig2 * np.linalg.inv(A.T @ A)[1, 1])
    return float(coef[1]), float(err), len(b)


def part_a(run, realpix, out):
    d = run.df
    bpix, lpix, wpix, dA, a_spring, a_autumn = realpix
    foot = d["foot"].to_numpy()
    cells = d.loc[foot, ["lon", "lat"]].round(3).drop_duplicates().to_numpy()
    sels = {"Roman detects": (d.detR == 1).to_numpy(),
            "Roman, sigma(M)/M<10%": ((d.detR == 1) & (d.relMl_R > 0) & (d.relMl_R < .1)).to_numpy(),
            "Rubin detects (footprint)": (d.detL == 1).to_numpy()}
    rng = np.random.default_rng(U.SEED + 20)
    for name, sel in sels.items():
        t = row_table(d, sel, "y")
        for bname, bm_t, bm_pix, bm_cell in (
                ("five-field", t.b < BLOCK_SPLIT, bpix < BLOCK_SPLIT, cells[:, 1] < BLOCK_SPLIT),
                ("GC field", t.b >= BLOCK_SPLIT, bpix >= BLOCK_SPLIT, cells[:, 1] >= BLOCK_SPLIT)):
            tb = t[bm_t]
            a, s, cov, chi2, dof = fit_block(tb)
            rho = lambda b, a=a, s=s: np.exp(a + s * b)
            area_sim = bm_cell.sum() * FINE ** 2
            area_real = (wpix[bm_pix]).sum() * dA
            n_sim = (rho(cells[bm_cell, 1]) * FINE ** 2).sum()
            n_real = (wpix[bm_pix] * rho(bpix[bm_pix])).sum() * dA
            # propagate the fit
            draws = rng.multivariate_normal([a, s], cov, size=N_DRAW)
            rr = [(wpix[bm_pix] * np.exp(aa + ss * bpix[bm_pix])).sum() * dA
                  / (np.exp(aa + ss * cells[bm_cell, 1]) * FINE ** 2).sum() for aa, ss in draws]
            lo, hi = np.percentile(rr, [16, 84])
            # the same, with the latitude slope steepened by OGLE-IV's (model too flat)
            s2 = s - OGLE_SLOPE[0]
            a2 = np.log((tb.rho * tb.area).sum() / (np.exp(s2 * tb.b) * tb.area).sum())
            ratio_ogle = ((wpix[bm_pix] * np.exp(a2 + s2 * bpix[bm_pix])).sum() * dA
                          / (np.exp(a2 + s2 * cells[bm_cell, 1]) * FINE ** 2).sum())
            out.append(dict(population=run.name, quantity=name, block=bname,
                            slope_per_deg=s, slope_err=np.sqrt(cov[1, 1]), chi2=chi2, dof=dof,
                            b_rows=f"{tb.b.min():.2f}..{tb.b.max():.2f}",
                            area_sim=area_sim, area_real_eff=area_real,
                            area_factor=area_real / area_sim,
                            latitude_factor=(n_real / n_sim) / (area_real / area_sim),
                            ratio=n_real / n_sim, ratio_lo=lo, ratio_hi=hi,
                            ratio_ogle_slope=ratio_ogle,
                            N1_sim_measured=float(tb.rho.mul(tb.area).sum())))
        print(f"[{run.name}] {name}: rows\n" + t.to_string(index=False), flush=True)
    # whole footprint (both blocks), per quantity
    df = pd.DataFrame([r for r in out if r["population"] == run.name])
    for name in sels:
        g = df[df.quantity == name]
        tot_sim = g.N1_sim_measured.sum()
        # scale each block's measured yield by its ratio; errors combined by sampling-free
        # linear propagation of the 16-84% half-widths (blocks independent)
        tot_real = (g.N1_sim_measured * g.ratio).sum()
        tot_ogle = (g.N1_sim_measured * g.ratio_ogle_slope).sum()
        half = np.sqrt(((g.N1_sim_measured * (g.ratio_hi - g.ratio_lo) / 2) ** 2).sum())
        out.append(dict(population=run.name, quantity=name, block="whole footprint",
                        area_sim=g.area_sim.sum(), area_real_eff=g.area_real_eff.sum(),
                        area_factor=g.area_real_eff.sum() / g.area_sim.sum(),
                        ratio=tot_real / tot_sim, ratio_lo=(tot_real - half) / tot_sim,
                        ratio_hi=(tot_real + half) / tot_sim, ratio_ogle_slope=tot_ogle / tot_sim,
                        latitude_factor=(tot_real / tot_sim) / (g.area_real_eff.sum() / g.area_sim.sum()),
                        N1_sim_measured=tot_sim))
    # per-event quantities against b, over Roman's footprint detections
    W = d["W"].to_numpy()
    det = foot & (d.detR == 1).to_numpy()
    charR = R.characterized(d, "roman").to_numpy()
    res5 = d["nres5_R"].to_numpy() >= U.RESOLVE_MIN_EPOCHS
    m10 = (d.relMl_R > 0).to_numpy() & (d.relMl_R < 0.1).to_numpy()
    lat = d["lat"].round(2).to_numpy()
    five = det & (lat < BLOCK_SPLIT)
    for qn, q in (("Roman-characterised fraction", charR), ("Roman sigma(M)/M<10% fraction", m10),
                  ("P(resolvable), Roman D=5 (as simulated)", res5)):
        bs, fs, es = [], [], []
        for b in np.unique(lat[five]):
            m = five & (lat == b)
            f = run.boot.frac(m & q, m, W)
            bs.append(b); fs.append(f[0]); es.append(f[1])
        bs, fs, es = map(np.asarray, (bs, fs, es))
        ok = es > 0
        wts = 1 / es[ok] ** 2
        A = np.vstack([np.ones(ok.sum()), bs[ok]]).T
        cov = np.linalg.inv(A.T @ (A * wts[:, None]))
        coef = cov @ (A.T @ (wts * fs[ok]))
        out.append(dict(population=run.name, quantity=qn, block="five-field, per-event trend",
                        slope_per_deg=coef[1], slope_err=np.sqrt(cov[1, 1]),
                        b_rows=" ".join(f"{b:.2f}:{100*f:.1f}+-{100*e:.1f}"
                                        for b, f, e in zip(bs, fs, es)),
                        ratio=(coef[0] + coef[1] * -1.40) / (coef[0] + coef[1] * -1.2)))


# ---------------------------------------------------------------------------------------------
# Part B
# ---------------------------------------------------------------------------------------------
def rubin_visits():
    v = np.loadtxt("Baseline/BulgeBaseline.dat", comments="#", usecols=(3, 4, 5))
    return v


def visits_at(l, b, v, full=False):
    """matchVisibleEpochs: flat distance <= FoV, strictly increasing times."""
    m = np.hypot(l - v[:, 0], b - v[:, 1]) <= FOV_RUBIN
    t = v[m, 2]
    if t.size == 0:
        return (0, 0, 0, 0) if full else (0, 0)
    # replicate "skip if cade <= 0 against the previous KEPT epoch"
    n, last = 0, -np.inf
    over = 0
    fr = np.array(FIELDS)
    dist = np.hypot(v[m, 0][:, None] - fr[:, 0], v[m, 1][:, None] - fr[:, 1])
    ov = dist.min(axis=1) <= FOV_RUBIN + FOV_ROMAN              # overlaps any Roman field
    ov_gc = dist[:, 5] <= FOV_RUBIN + FOV_ROMAN                 # overlaps the Galactic-centre field
    ov_5 = dist[:, :5].min(axis=1) <= FOV_RUBIN + FOV_ROMAN     # overlaps the five-field block
    gc = five = 0
    for ti, oi, og, o5 in zip(t, ov, ov_gc, ov_5):
        if ti > last:
            n += 1
            over += int(oi)
            gc += int(og)
            five += int(o5)
            last = ti
    return (n, over, gc, five) if full else (n, over)


def logged():
    lines = subprocess.run(["grep", "-E", r"^longtitude:|^ndd \(LSST", LOG], capture_output=True,
                           text=True, errors="replace", check=True).stdout.splitlines()
    out, key = {}, None
    for line in lines:
        if line.startswith("longtitude:"):
            p = line.split()
            key = (round(float(p[1]), 3), round(float(p[3]), 3))
        elif key:
            out[key] = int(line.split()[2])
            key = None
    return out


def in_scan(l, b):
    """The scanned region: the box minus the corner cut (Bulge.h)."""
    return ((l >= LON_MIN) & (l < LON_MAX + FINE) & (b >= LAT_MIN) & (b < LAT_MAX + FINE)
            & ~((l < LX) & (b > BX)))


def completeness(v, out):
    """Of the sky that Rubin images in exposures which also contain a Roman field, how much did
    the scan cover? This is the scan's design intent (Rubin's full field around Roman's)."""
    fr = np.array(FIELDS)
    over = (np.hypot(v[:, 0][:, None] - fr[:, 0], v[:, 1][:, None] - fr[:, 1]).min(axis=1)
            <= FOV_RUBIN + FOV_ROMAN)
    ptg = np.unique(np.round(v[over, :2], 4), axis=0)          # distinct overlapping pointings
    ls = np.arange(-6.5, 8.0, FINE) + FINE / 2
    bs = np.arange(-7.5, 5.5, FINE) + FINE / 2
    L, B = np.meshgrid(ls, bs)
    L, B = L.ravel(), B.ravel()
    hit = np.zeros(L.size, bool)
    for i in range(0, ptg.shape[0], 200):
        p = ptg[i:i + 200]
        hit |= (np.hypot(L[:, None] - p[:, 0], B[:, None] - p[:, 1]) <= FOV_RUBIN).any(axis=1)
    sc = in_scan(L, B)
    a_hit, a_in = hit.sum() * FINE ** 2, (hit & sc).sum() * FINE ** 2
    cut = hit & (L < LX) & (B > BX) & (L >= LON_MIN) & (B < LAT_MAX + FINE)
    outside_box = hit & ~((L >= LON_MIN) & (L < LON_MAX + FINE) & (B >= LAT_MIN) & (B < LAT_MAX + FINE))
    print(f"sky imaged by Rubin exposures that overlap a Roman field: {a_hit:.2f} deg^2 "
          f"({len(ptg)} distinct pointings); inside the scan {a_in:.2f} deg^2 = {100*a_in/a_hit:.0f}%; "
          f"lost to the corner cut {cut.sum()*FINE**2:.2f}; beyond the box {outside_box.sum()*FINE**2:.2f}",
          flush=True)
    out.append(dict(population="scan", quantity="Rubin sky co-imaged with Roman fields",
                    block="completeness", area_sim=a_in, area_real_eff=a_hit, ratio=a_in / a_hit,
                    b_rows=f"corner cut {cut.sum()*FINE**2:.2f} deg^2; beyond box "
                           f"{outside_box.sum()*FINE**2:.2f} deg^2; {len(ptg)} pointings"))


def part_b(runs, out):
    v = rubin_visits()
    log = logged()
    completeness(v, out)
    # 1. the visit counter against the simulator's own, on every simulated sightline
    keys = sorted(log)
    mism = [k for k in keys if visits_at(k[0], k[1], v)[0] != log[k]]
    print(f"visit counter vs run log: {len(keys) - len(mism)}/{len(keys)} sightlines agree"
          + (f"; e.g. {mism[:3]} -> {[(visits_at(*k, v)[0], log[k]) for k in mism[:3]]}"
             if mism else ""), flush=True)
    out.append(dict(population="scan", quantity="visit counter reproduces run log",
                    block="all sightlines", ratio=(len(keys) - len(mism)) / len(keys)))
    # 2. the two upper corners, cell by cell
    width = LX - LON_MIN
    corners = {"cut (l < -0.94, b > 0.31)": (LON_MIN, LX),
               "opposite (high l, b > 0.31)": (LON_MAX - width, LON_MAX)}
    cell_rows = []
    for cname, (l_lo, l_hi) in corners.items():
        for l in np.arange(l_lo + FINE / 2, l_hi, FINE):
            for b in np.arange(BX + FINE / 2, LAT_MAX, FINE):
                n, over, gc, five = visits_at(l, b, v, full=True)
                cell_rows.append(dict(corner=cname, l=l, b=b, visits=n, visits_overlapping_roman=over,
                                      visits_overlapping_gc=gc, visits_overlapping_five=five))
    cells = pd.DataFrame(cell_rows)
    for cname, g in cells.groupby("corner"):
        out.append(dict(population="scan", quantity=cname, block="corner",
                        area_sim=len(g) * FINE ** 2,
                        area_real_eff=(g.visits > 0).sum() * FINE ** 2,
                        N1_sim_measured=float((g.visits_overlapping_roman > 0).sum() * FINE ** 2),
                        b_rows=f"median visits {g.visits.median():.0f}, "
                               f"median overlapping-Roman visits {g.visits_overlapping_roman.median():.0f}"))
        print(f"{cname}: {len(g)} cells = {len(g)*FINE**2:.2f} deg^2; with Rubin visits "
              f"{(g.visits>0).sum()*FINE**2:.2f} deg^2; imaged by Rubin exposures that also contain "
              f"Roman's fields {(g.visits_overlapping_roman>0).sum()*FINE**2:.2f} deg^2; median visits "
              f"{g.visits.median():.0f}; imaged with the GC field {(g.visits_overlapping_gc>0).sum()*FINE**2:.2f} "
              f"deg^2 ({g.visits_overlapping_gc.sum():,} cell-visits), with the five-field block "
              f"{(g.visits_overlapping_five>0).sum()*FINE**2:.2f} deg^2 ({g.visits_overlapping_five.sum():,} "
              f"cell-visits)", flush=True)
        out.append(dict(population="scan", quantity=cname, block="corner, by Roman block",
                        area_sim=(g.visits_overlapping_gc > 0).sum() * FINE ** 2,
                        area_real_eff=(g.visits_overlapping_five > 0).sum() * FINE ** 2,
                        b_rows="area_sim = imaged with GC field; area_real_eff = with five-field block"))
    # 3. Rubin yield lost from the cut corner: density vs visit count on simulated northern
    #    sightlines (b > 0.31, outside the footprint), applied to the corner's cells
    for run in runs:
        d = run.df
        key = pd.Series(list(zip(d["lon"].round(3), d["lat"].round(3))))
        det_by_sl = d.loc[(d.detL == 1).to_numpy(), "y"].groupby(key[(d.detL == 1).to_numpy()]).sum()
        area_by_sl = d["w_area"].groupby(key).first()
        foot_sl = set(key[d["foot"].to_numpy()])
        # Every northern sightline OUTSIDE the footprint, detections or not (a sightline with no
        # detection is absent from the detection-only table but is a real zero). Area: the
        # table's w_area where known, else a full coarse block.
        rows = []
        for k, nv in log.items():
            if k[1] > BX and k not in foot_sl:
                a = float(area_by_sl.get(k, 0.04))
                rows.append(dict(nvis=nv, rho=float(det_by_sl.get(k, 0.0)) / a))
        per = pd.DataFrame(rows)
        bins = np.array([0, 15, 50, 100, 200, 400, 800, 1600, 5000])
        per["bin"] = np.digitize(per.nvis, bins)
        med = per.groupby("bin")["rho"].median()
        cut = cells[cells.corner.str.startswith("cut")].copy()
        cut["bin"] = np.digitize(cut.visits, bins)
        lost = float((cut["bin"].map(med).fillna(0) * FINE ** 2).sum())
        total = float(d.loc[d.detL == 1, "y"].sum())
        out.append(dict(population=run.name, quantity="Rubin detections lost to the corner cut",
                        block="estimate", N1_sim_measured=lost, ratio=lost / total,
                        b_rows=f"of whole-scan Rubin N_1 {total:.4g}; {len(per)} northern sightlines"))
        print(f"[{run.name}] Rubin yield in the cut corner ~ {lost:.4g} (N_1), "
              f"{100*lost/total:.1f}% of the whole-scan Rubin yield "
              f"(density from {len(per)} northern sightlines by visit bin: "
              f"{med.round(0).to_dict()})", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", action="append", required=True, metavar="NAME=DIR")
    ap.add_argument("--mean-mass", required=True)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    global OGLE_SLOPE
    OGLE_SLOPE = ogle_slope_correction()
    print(f"d ln(model/OGLE)/db over -1.95..-0.75: {OGLE_SLOPE[0]:+.3f} +- {OGLE_SLOPE[1]:.3f} per deg "
          f"({OGLE_SLOPE[2]} latitude bins)", flush=True)
    out_ogle = dict(population="scan", quantity="d ln(model/OGLE rate per star)/db",
                    block="footprint latitudes", slope_per_deg=OGLE_SLOPE[0],
                    slope_err=OGLE_SLOPE[1])
    realpix = real_pixels()
    print(f"real tiles: spring {realpix[4]:.3f} deg^2, autumn {realpix[5]:.3f} deg^2, "
          f"season-weighted {realpix[2].sum()*realpix[3]:.3f} deg^2, union "
          f"{(realpix[2] > 0).sum()*realpix[3]:.3f} deg^2", flush=True)
    out, runs = [out_ogle], []
    for spec in a.run:
        name, directory = spec.split("=", 1)
        run = U.Run(name, directory, a.mean_mass)
        part_a(run, realpix, out)
        runs.append(run)
    part_b(runs, out)
    t = pd.DataFrame(out)
    t.to_csv(os.path.join(a.out, "u3_footprint_offset.csv"), index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30,
                           "display.max_colwidth", 90):
        print(t.to_string(index=False))


if __name__ == "__main__":
    main()
