#!/usr/bin/env python3
"""Step U6: an independent check of the two dust maps on Roman's sightlines, against VVV.

WHY THIS EXISTS. The dust correction (Step U5) needs a reference extinction. Two maps are
available through dustmaps: DECaPS (Zucker et al. 2025; optical + near-infrared stars, sensitive to
A_V ~ 12 only) and Marshall et al. (2006; 2MASS giants, A_Ks). On Roman's five-field block they
agree; on the Galactic-centre field (|b| < 0.35) DECaPS gives A_V ~ 7 at 8 kpc and Marshall ~25.
They cannot both be right, and choosing between them by assumption would be a guess. This script
settles it with a third, independent measurement: the VVV reddening map of Surot et al. (2020,
A&A 644, A140; VizieR J/A+A/644/A140), E(J-Ks) of the bulge red clump and red giants on a grid of
10 arcsec to 2 arcmin -- i.e. the total column to the bulge, measured in the near-infrared where
A_V ~ 30 is still transparent.

THE TEST IS LAW-FREE. Every map is converted to A_V with some extinction law, and the law toward the
bulge is itself uncertain. So the test compares CONTRASTS, not values: the ratio of each map's
column in the Galactic-centre field to its column in the five-field block. The law cancels in the
ratio. A map that sees the whole column must reproduce VVV's contrast; a saturated one falls short.

It also measures the calibration the correction uses where it has to fall back on the
near-infrared map: the A_Ks/A_V at which Marshall's A_Ks matches DECaPS's A_V where DECaPS is valid
(the five-field block, reliable flag true on all 119 sightlines at 8 kpc).

Queries VizieR once (147 box searches, a few minutes) and caches the medians.

    .roman/bin/python analysis/u6_vvv_check.py -o figures/u1_20260929
"""

import argparse
import io
import os
import sys
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

VIZIER = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
HALF = 0.05                            # deg: the box is the 0.1-deg scan cell around the sightline


def surot_median(l, b, half=HALF):
    """Median VVV E(J-Ks) (and node count) in the box |l'-l|, |b'-b| <= half.

    A box on the catalogue's own GLON/GLAT columns, not a cone search: the table stores l < 0 as
    negative GLON, and VizieR's cone search around l < 0 returns nothing (verified 2026-09-30 at
    l = 359.9, b = -0.34, where the box query returns tile b333's nodes)."""
    q = {"-source": "J/A+A/644/A140/ejkmap", "-out": "GLON,GLAT,E(J-Ks)",
         "GLON": f"{l - half:.4f}..{l + half:.4f}", "GLAT": f"{b - half:.4f}..{b + half:.4f}",
         "-out.max": "1000000"}
    url = VIZIER + "?" + urllib.parse.urlencode(q)
    for attempt in range(4):
        try:
            txt = urllib.request.urlopen(url, timeout=120).read().decode()
            break
        except OSError:
            time.sleep(5 * (attempt + 1))
    else:
        raise RuntimeError(f"VizieR query failed at ({l}, {b})")
    rows = [r.split("\t") for r in txt.splitlines() if r and not r.startswith("#")]
    vals = [float(r[2]) for r in rows[3:] if len(r) >= 3 and r[2].strip()]
    return (float(np.median(vals)) if vals else np.nan), len(vals)


BBINS = [0.0, 0.5, 1.0, 1.5, 2.5, 6.0]      # |b| bins of the report's dust table


def scan_ejk(log, out_dir, workers=4):
    """VVV E(J-Ks) at every scan sightline (from a production run.log), cached."""
    import subprocess
    from concurrent.futures import ThreadPoolExecutor
    cache = os.path.join(out_dir, "u6_vvv_scan_ejk.csv")
    if os.path.exists(cache):
        return pd.read_csv(cache)
    lines = subprocess.run(["grep", "-E", "^longtitude:", log], capture_output=True,
                           text=True).stdout.splitlines()
    pos = sorted({(float(p.split()[1]), float(p.split()[3])) for p in lines})
    with ThreadPoolExecutor(workers) as ex:
        res = list(ex.map(lambda lb: surot_median(*lb), pos))
    t = pd.DataFrame([dict(l=l, b=b, ejk=e, n_nodes=n) for (l, b), (e, n) in zip(pos, res)])
    t.to_csv(cache, index=False)
    return t


def scan_table(ejk, out_dir, cal_aks, cal_ejk):
    """By |b| over the whole scan: each dust column at 8 kpc against VVV's, all on DECaPS's A_V
    scale (Marshall and VVV calibrated on the five-field block)."""
    import u5_corrected_numbers as U5
    dust = U5.Dust()
    i8 = int(np.argmin(np.abs(U5.DGRID - 8.0)))
    rows = []
    for l, b in zip(ejk.l, ejk.b):
        dav, rel = dust.decaps_profile(l, b)
        nom = dust.reference_profile(l, b, "hybrid", U5.AKS_AV["nominal"])
        rows.append(dict(l=l, b=b, model=float(dust.model_av(l, b, np.array([8.0]))[0]),
                         decaps=float(np.interp(8.0, U5.DGRID, dav)),
                         aks=float(np.interp(8.0, U5.DGRID, dust.aks_profile(l, b))),
                         nominal=float(np.interp(8.0, U5.DGRID, nom)),
                         decaps_reliable_8kpc=bool(rel[i8])))
    t = pd.DataFrame(rows).merge(ejk, on=["l", "b"], validate="one_to_one")
    t["vvv_av"] = t.ejk / cal_ejk
    t["nir_av"] = t.aks / cal_aks
    t["absb"] = t.b.abs()
    t["bbin"] = pd.cut(t.absb, BBINS, include_lowest=True)
    t.to_csv(os.path.join(out_dir, "u6_vvv_scan.csv"), index=False)
    out = t.groupby("bbin", observed=True).apply(lambda g: pd.Series(dict(
        n=len(g), vvv_av=g.vvv_av.median(),
        model=(g.model / g.vvv_av).median(), decaps=(g.decaps / g.vvv_av).median(),
        marshall=(g.nir_av / g.vvv_av).median(), nominal=(g.nominal / g.vvv_av).median(),
        decaps_reliable=g.decaps_reliable_8kpc.mean())))
    print("\nWhole scan, per |b| bin: median over sightlines of each map's A_V(8 kpc) / VVV's "
          "(all on DECaPS's scale):")
    print(out.round(3).to_string())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-o", "--out-dir", required=True)
    ap.add_argument("--summary", default=None,
                    help="u5_dust_summary.csv (default: <out-dir>/u5_dust_summary.csv)")
    ap.add_argument("--scan-log", default=None,
                    help="a production run.log: also compare the maps over the whole scan")
    ap.add_argument("--queries-only", action="store_true",
                    help="with --scan-log: only fetch and cache VVV (no dust maps loaded)")
    a = ap.parse_args()
    if a.scan_log and a.queries_only:
        t = scan_ejk(a.scan_log, a.out_dir)
        print(f"cached VVV at {len(t)} scan sightlines; {int((t.n_nodes == 0).sum())} without nodes")
        return
    summ = pd.read_csv(a.summary or os.path.join(a.out_dir, "u5_dust_summary.csv"))
    cache = os.path.join(a.out_dir, "u6_vvv_ejk.csv")
    if os.path.exists(cache):
        ejk = pd.read_csv(cache)
    else:
        rec = []
        for l, b in zip(summ.l, summ.b):
            e, n = surot_median(l, b)
            rec.append(dict(l=l, b=b, ejk=e, n_nodes=n))
        ejk = pd.DataFrame(rec)
        ejk.to_csv(cache, index=False)
    d = summ.merge(ejk, on=["l", "b"], validate="one_to_one")
    g = d.groupby("block")
    med = g[["model", "decaps", "aks", "ejk"]].median()
    print(f"Median column to 8 kpc on Roman's sightlines (VVV: median E(J-Ks) in the "
          f"{2 * HALF:.1f}-deg cell; {int(d.n_nodes.min())}-{int(d.n_nodes.max())} nodes each):")
    print(med.round(3).to_string())
    c = med.loc["GC field"] / med.loc["five-field"]
    print("\nContrast GC field / five-field block (law-free):")
    for k, name in [("ejk", "VVV E(J-Ks) (Surot+2020)"), ("aks", "Marshall A_Ks"),
                    ("decaps", "DECaPS A_V"), ("model", "simulator's tables A_V")]:
        print(f"  {name:28s} {c[k]:5.2f}")
    ff = d[d.block == "five-field"]
    cal = ff["aks"] / ff["decaps"]
    print(f"\nCalibration on the five-field block (DECaPS reliable at 8 kpc on "
          f"{ff.decaps_reliable_8kpc.mean():.0%}): A_Ks(Marshall)/A_V(DECaPS) median "
          f"{cal.median():.4f}, 16-84% {cal.quantile(.16):.4f}-{cal.quantile(.84):.4f}; "
          f"ratio of medians {ff.aks.median() / ff.decaps.median():.4f}")
    r = ff["ejk"] / ff["decaps"]
    print(f"  E(J-Ks)(VVV)/A_V(DECaPS) there: median {r.median():.4f}")
    rg = d[d.block == "GC field"]
    print(f"  same ratio on the GC field: {(rg.ejk / rg.decaps).median():.4f}  "
          f"(a saturated DECaPS makes this larger)")
    print(f"\nGC field column to 8 kpc on DECaPS's A_V scale (median): DECaPS itself "
          f"{rg.decaps.median():.1f}; Marshall calibrated {rg.aks.median() / cal.median():.1f}; "
          f"VVV calibrated {rg.ejk.median() / r.median():.1f}; the simulator's tables "
          f"{rg.model.median():.1f}")
    print(f"A_Ks/A_V that makes Marshall's contrast match VVV's: "
          f"{cal.median() * c['aks'] / c['ejk']:.4f}")
    d.to_csv(os.path.join(a.out_dir, "u6_vvv_check.csv"), index=False)
    if a.scan_log:
        scan_table(scan_ejk(a.scan_log, a.out_dir), a.out_dir, cal.median(), r.median())


if __name__ == "__main__":
    main()
