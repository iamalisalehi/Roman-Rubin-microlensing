#!/usr/bin/env python3
"""Step V1: validate the rebuilt extinction tables (files/ext/ext_tables.dat) against VVV.

WHY THIS EXISTS. Deviation 70 rebuilt the simulator's extinction tables from DECaPS + Marshall
(maps.py, analysis/dustref.py) because the old ones were 3-6x too thin within 1 deg of the plane
(Deviations 61-63). The OPEN_ITEMS fix plan sets the bar before any run: against the independent VVV
reddening map (Surot et al. 2020; E(J-Ks) of the bulge red clump and giants, i.e. the total column
to the bulge, measured in the near-infrared),
  (a) A_V(8 kpc) / VVV's A_V should be ~0.9-1.0 in every |b| bin, with no fall toward the plane
      (the old tables: 0.17 at |b| < 0.5, 0.31 at 0.5-1);
  (b) the Galactic-centre field / five-field block contrast should be ~4.5-4.9 (VVV 4.89, Marshall
      4.46; the old tables 0.57). The contrast is LAW-FREE: any A_V conversion cancels in it.

HOW. Every sightline of the scan the simulator now builds (analysis/gbtds_geometry, the production
grid --stride 10 --stride-roman 5), and for each the table the simulator would use (nearest
position, first minimum -- helper.cpp nearestSightline), A_V interpolated at 8 kpc. VVV: the median
E(J-Ks) in the 0.1-deg box around the sightline (u6_vvv_check.surot_median; cached). VVV is put on
the tables' A_V scale with E(J-Ks)/A_V measured on the five-field block (where DECaPS, which the
tables follow there, is reliable). The GC field / five-field blocks are the footprint sightlines on
a detector of the GC field / of fields 0-4, either roll.

    .roman/bin/python analysis/v1_ext_vvv.py -o figures/ext_20261001
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gbtds_geometry as G           # noqa: E402

TABLES = os.path.join(G.ROOT, "files", "ext", "ext_tables.dat")
BBINS = [0.0, 0.5, 1.0, 1.5, 2.5, 6.0]


def read_tables(path=TABLES):
    dist = None
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("# dist "):
                dist = np.array(line.split()[2:], float)
    t = np.loadtxt(path, comments="#", dtype=np.float32)
    return t[:, 0].astype(float), t[:, 1].astype(float), t[:, 2:], dist


def vvv_ejk(pos, out_dir, workers=6):
    """VVV E(J-Ks) at each (l, b), cached in out_dir/v1_vvv_ejk.csv (resumable)."""
    from concurrent.futures import ThreadPoolExecutor
    import u6_vvv_check as U6
    cache = os.path.join(out_dir, "v1_vvv_ejk.csv")
    have = pd.read_csv(cache) if os.path.exists(cache) else pd.DataFrame(columns=["l", "b", "ejk", "n_nodes"])
    done = {(round(l, 4), round(b, 4)) for l, b in zip(have.l, have.b)}
    todo = [p for p in pos if (round(p[0], 4), round(p[1], 4)) not in done]
    print(f"VVV: {len(pos) - len(todo)} cached, {len(todo)} to query", flush=True)
    for s in range(0, len(todo), 60):
        chunk = todo[s:s + 60]
        with ThreadPoolExecutor(workers) as ex:
            res = list(ex.map(lambda lb: U6.surot_median(*lb), chunk))
        add = pd.DataFrame([dict(l=l, b=b, ejk=e, n_nodes=n) for (l, b), (e, n) in zip(chunk, res)])
        have = pd.concat([have, add], ignore_index=True)
        have.to_csv(cache, index=False)
        print(f"  VVV {s + len(chunk)}/{len(todo)}", flush=True)
    key = {(round(l, 4), round(b, 4)): (e, n) for l, b, e, n in zip(have.l, have.b, have.ejk, have.n_nodes)}
    return np.array([key[(round(l, 4), round(b, 4))] for l, b in pos], float)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-o", "--out-dir", required=True)
    ap.add_argument("--vvv-only", action="store_true", help="only fetch and cache VVV")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    g = G.scan_sightlines(10, 5)
    lon, lat = g["lon"], g["lat"]
    ejk = vvv_ejk(list(zip(lon, lat)), a.out_dir)
    if a.vvv_only:
        return

    tl, tb, av, dist = read_tables()
    i8 = int(np.argmin(np.abs(dist - 8.0)))
    # nearestSightline: first strict minimum of the flat distance, in file order.
    near = np.array([int(np.argmin((tl - l) ** 2 + (tb - b) ** 2)) for l, b in zip(lon, lat)])
    t = pd.DataFrame(dict(l=lon, b=lat, fine=g["fine"], area=g["area"],
                          table_l=tl[near], table_b=tb[near],
                          av8=av[near, i8].astype(float),
                          ejk=ejk[:, 0], n_nodes=ejk[:, 1]))
    t["sep_deg"] = np.hypot(t.l - t.table_l, t.b - t.table_b)
    # Blocks: footprint sightlines on a detector of the GC field (field 5) / of fields 0-4.
    visits = G.read_roman_visits()
    pl = visits[["l", "b", "field", "layout"]].drop_duplicates().to_numpy()
    gc = np.zeros(len(t), bool)
    ff = np.zeros(len(t), bool)
    for x, y, f, k in pl:
        on = G.in_detector(t.l - x, t.b - y, int(k))
        (gc if int(f) == 5 else ff)[on] = True
    t["block"] = np.where(gc, "GC field", np.where(ff, "five-field", ""))
    ok = np.isfinite(t.ejk) & (t.n_nodes > 0)
    cal = (t.ejk / t.av8)[ok & ff].median()              # E(J-Ks) per A_V, five-field block
    t["vvv_av"] = t.ejk / cal
    t["ratio"] = t.av8 / t.vvv_av
    t["absb"] = t.b.abs()
    t["bbin"] = pd.cut(t.absb, BBINS, include_lowest=True)
    t.to_csv(os.path.join(a.out_dir, "v1_ext_vvv.csv"), index=False)

    print(f"{len(t)} scan sightlines; nearest table within {t.sep_deg.max():.3f} deg "
          f"(median {t.sep_deg.median():.3f}); VVV nodes found at {ok.sum()}")
    print(f"VVV -> A_V on the five-field block: E(J-Ks)/A_V = {cal:.4f} "
          f"({(ok & ff).sum()} sightlines)")
    out = t[ok].groupby("bbin", observed=True).apply(lambda q: pd.Series(dict(
        n=len(q), vvv_av=q.vvv_av.median(), tables_av=q.av8.median(), ratio=q.ratio.median(),
        ratio_p16=q.ratio.quantile(.16), ratio_p84=q.ratio.quantile(.84))))
    print("\n(a) per |b| bin, A_V(8 kpc) tables / VVV (both on the tables' scale):")
    print(out.round(3).to_string())
    m = t[ok].groupby("block")[["av8", "ejk"]].median()
    c_tab = m.loc["GC field", "av8"] / m.loc["five-field", "av8"]
    c_vvv = m.loc["GC field", "ejk"] / m.loc["five-field", "ejk"]
    print(f"\n(b) contrast GC field / five-field block (law-free): tables {c_tab:.2f}, VVV {c_vvv:.2f}"
          f"   [target ~4.5-4.9; old tables 0.57]")
    with open(os.path.join(a.out_dir, "v1_ext_vvv.md"), "w") as f:
        f.write("# V1: rebuilt extinction tables vs VVV (analysis/v1_ext_vvv.py)\n\n")
        f.write(f"Tables: `files/ext/ext_tables.dat`; {len(t)} scan sightlines (production grid); "
                f"E(J-Ks)/A_V = {cal:.4f} on the five-field block.\n\n")
        f.write("(a) A_V(8 kpc) tables / VVV per |b| bin (target ~0.9-1.0; old tables 0.17 at |b|<0.5)\n\n")
        f.write("```\n" + out.round(3).to_string() + "\n```\n\n")
        f.write(f"(b) GC field / five-field contrast: tables {c_tab:.2f}, VVV {c_vvv:.2f} "
                f"(target ~4.5-4.9; old tables 0.57)\n")


if __name__ == "__main__":
    main()
