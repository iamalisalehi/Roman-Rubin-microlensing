"""Build the simulator's extinction tables: files/ext/ext_tables.dat (Deviation 70).

WHAT CHANGED AND WHY. This script used to write one Bayestar19/DECaPS table per Rubin pointing
centre, choosing the map by declination (Bayestar north of -30, the dustmaps documentation's rule).
That put 139 of Roman's 147 sightlines on Bayestar beyond its own reliable distance, left 78 tables
all-NaN (read by the simulator as zero dust), and made the dust 3-6x too thin within 1 deg of the
plane against the VVV reddening map (Deviations 61-63, OPEN_ITEMS CRITICAL dust entry). Now:

  WHERE  a regular grid over the simulator's scan region (analysis/gbtds_geometry): 0.1 deg
         everywhere, 0.05 deg within |b| < 1.5 and over Roman's footprint, where the dust varies
         fastest -- no longer tied to Rubin's pointing pattern.
  WHAT   the reference profile of analysis/dustref.py: DECaPS where it can see, Marshall's
         near-infrared map (calibrated onto DECaPS's scale by k = A_Ks/A_V, re-measured here on
         Roman's five-field block) where it cannot; non-decreasing in distance. Both through the
         dustmaps library. Bayestar is not used.
  FORMAT one text file, one line per sky position: `l b A_V(d_1) ... A_V(d_n)`, the distance grid
         and the provenance in `#` header lines; src/galaxy/extinction.cpp readExtinction() reads it line by line
         into ~25 MB (float). Provenance also in ext_provenance.json (not .txt: nothing in files/ext
         is globbed any more, but the old reader read every .txt there).

The DECaPS queries cost ~0.13 s per sky position (~35 min for the grid) and are cached in
files/ext_raw/raw_<grid hash>.npz, so a rebuild with another k, or after a rule change, is
seconds:

    .roman/bin/python maps.py                 # query (or reuse the cache) and write the tables
    .roman/bin/python maps.py --k 0.0805      # force k instead of measuring it

Exit status is non-zero if any table value is non-finite or any profile decreases.
"""

import argparse
import datetime as dt
import hashlib
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
import dustref as D                     # noqa: E402
import gbtds_geometry as G              # noqa: E402

OUT_DIR = os.path.join(ROOT, "files", "ext")
OUT = os.path.join(OUT_DIR, "ext_tables.dat")
PROV = os.path.join(OUT_DIR, "ext_provenance.json")
CACHE_DIR = os.path.join(ROOT, "files", "ext_raw")
STEP = 0.05            # fine lattice [deg]; the coarse one is every second point
FINE_ABS_B = 1.5       # 0.05-deg tables within this |b| ...
MARGIN = 0.15          # ... and the grid reaches this far past the scan region [deg]


def table_positions():
    """Lattice points (multiples of STEP) inside the scan region + MARGIN: all at 0.1 deg, plus
    the 0.05-deg ones within |b| < FINE_ABS_B or over the footprint's detector outlines."""
    pl = G.placements()
    reach = G.scan_reach() + MARGIN
    lmin, lmax = min(p[0] for p in pl) - reach, max(p[0] for p in pl) + reach
    bmin, bmax = min(p[1] for p in pl) - reach, max(p[1] for p in pl) + reach
    i = np.arange(int(np.floor(lmin / STEP)), int(np.ceil(lmax / STEP)) + 1)
    j = np.arange(int(np.floor(bmin / STEP)), int(np.ceil(bmax / STEP)) + 1)
    I, J = np.meshgrid(i, j, indexing="ij")
    I, J = I.ravel(), J.ravel()
    l, b = I * STEP, J * STEP
    near = np.min([np.hypot(l - x, b - y) for x, y, _ in pl], axis=0) <= reach
    foot = np.zeros(l.size, bool)
    for x, y, k in pl:
        dl0, dl1, db0, db1 = G.outline_bbox(k)
        foot |= ((l - x >= dl0 - 0.1) & (l - x <= dl1 + 0.1)
                 & (b - y >= db0 - 0.1) & (b - y <= db1 + 0.1))
    coarse = (I % 2 == 0) & (J % 2 == 0)
    keep = near & (coarse | (np.abs(b) < FINE_ABS_B) | foot)
    # Spatial order (by b, then l) keeps DECaPS's disk reads local.
    order = np.lexsort((l[keep], b[keep]))
    return np.round(l[keep][order], 4), np.round(b[keep][order], 4)


def query_raw(l, b, chunk):
    key = hashlib.sha1(np.concatenate([l, b, D.DGRID]).tobytes()).hexdigest()[:12]
    path = os.path.join(CACHE_DIR, f"raw_{key}.npz")
    if os.path.exists(path):
        z = np.load(path)
        print(f"raw maps from cache {path}")
        return z["dav"], z["rel"], z["aks"], path
    os.makedirs(CACHE_DIR, exist_ok=True)
    dust = D.ReferenceDust()
    n, m = l.size, D.DGRID.size
    dav = np.empty((n, m), np.float32)
    rel = np.empty((n, m), bool)
    aks = np.empty((n, m), np.float32)
    t0 = time.time()
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        dav[s:e], rel[s:e], aks[s:e] = dust.batch(l[s:e], b[s:e])
        el = time.time() - t0
        print(f"  queried {e}/{n} positions, {el/60:.1f} min, ~{el/e*(n-e)/60:.0f} min left",
              flush=True)
    np.savez(path, l=l, b=b, dav=dav, rel=rel, aks=aks)
    print(f"raw maps cached in {path}")
    return dav, rel, aks, path


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--k", type=float, default=None,
                    help="A_Ks/A_V for the near-infrared map (default: measured on the five-field "
                         "block at 8 kpc)")
    ap.add_argument("--chunk", type=int, default=200, help="positions per DECaPS query")
    a = ap.parse_args()

    l, b = table_positions()
    print(f"{l.size} table positions (l {l.min():.2f}..{l.max():.2f}, b {b.min():.2f}..{b.max():.2f}),"
          f" {D.DGRID.size} distances")
    dav, rel, aks, cache = query_raw(l, b, a.chunk)
    dav, aks = dav.astype(float), aks.astype(float)

    # k on Roman's five-field block (fields 0-4, either roll), at 8 kpc, where DECaPS is reliable.
    i8 = int(np.argmin(np.abs(D.DGRID - 8.0)))
    block = np.zeros(l.size, bool)
    for x, y, k in G.placements():
        if y < -0.8:
            block |= G.in_detector(l - x, b - y, k)
    ok = block & rel[:, i8] & np.isfinite(dav[:, i8]) & np.isfinite(aks[:, i8]) & (dav[:, i8] > 0)
    ratio = aks[ok, i8] / dav[ok, i8]
    k_meas = (float(np.median(ratio)), float(np.percentile(ratio, 16)),
              float(np.percentile(ratio, 84)), int(ok.sum()), int(block.sum()))
    k = a.k if a.k is not None else k_meas[0]
    print(f"k = A_Ks/A_V: measured {k_meas[0]:.4f} (16-84%: {k_meas[1]:.4f}-{k_meas[2]:.4f}) on "
          f"{k_meas[3]} of {k_meas[4]} five-field positions; using {k:.4f}"
          f"{' (forced)' if a.k is not None else ''}; notional-layout value was {D.K_NOMINAL_LEGACY}")

    prof = np.empty((l.size, D.DGRID.size))
    dsat = np.empty(l.size)
    lifted = 0
    for n in range(l.size):
        raw, _ = D.combine(dav[n], rel[n], aks[n], "hybrid", k, monotone=False)
        prof[n], dsat[n] = D.combine(dav[n], rel[n], aks[n], "hybrid", k, monotone=True)
        lifted += not np.array_equal(raw, prof[n])

    bad = ~np.isfinite(prof)
    dec = np.diff(prof, axis=1) < 0
    if bad.any() or dec.any():
        print(f"ERROR: {bad.any(axis=1).sum()} positions with non-finite values, "
              f"{dec.any(axis=1).sum()} decreasing; no tables written.", file=sys.stderr)
        sys.exit(1)

    used_nir = np.isfinite(dsat)
    a8 = prof[:, i8]
    print(f"switched to Marshall outright at {used_nir.sum()} positions (median from "
          f"{np.median(dsat[used_nir]) if used_nir.any() else float('nan'):.2f} kpc); monotone fix "
          f"changed {lifted}; A_V(8 kpc) median {np.median(a8):.2f}, |b|<0.5: "
          f"{np.median(a8[np.abs(b) < 0.5]):.2f}, max {a8.max():.2f}")

    import dustmaps
    from importlib.metadata import version
    stamp = dt.datetime.now().isoformat(timespec="seconds")
    os.makedirs(OUT_DIR, exist_ok=True)
    tmp = OUT + ".tmp"
    with open(tmp, "w") as f:
        f.write(f"# ext_tables v1 -- built by maps.py {stamp} (Deviation 70)\n")
        f.write(f"# A_V(d): DECaPS (DECaPSQueryLite mean, R_V {D.RV_DECAPS}) where reliable and not "
                f"saturated, Marshall A_Ks/k beyond; non-decreasing; dustmaps {version('dustmaps')}\n")
        f.write(f"# k {k:.6f}\n")
        f.write(f"# n_tables {l.size}\n")
        f.write(f"# n_dist {D.DGRID.size}\n")
        f.write("# dist " + " ".join(f"{d:.2f}" for d in D.DGRID) + "\n")
        for n in range(l.size):
            f.write(f"{l[n]:.4f} {b[n]:.4f} " + " ".join(f"{v:.4f}" for v in prof[n]) + "\n")
    os.replace(tmp, OUT)
    prov = dict(built=stamp, script="maps.py", deviation=70, out=os.path.relpath(OUT, ROOT),
                dustmaps_version=version("dustmaps"), dustmaps_data_dir="dustmaps",
                decaps="DECaPSQueryLite(mean_only=True), A_V = 3.32 E(B-V)",
                marshall="MarshallQuery, A_V = A_Ks / k", decaps_av_max=D.DECAPS_AV_MAX,
                k_used=k, k_forced=a.k is not None,
                k_measured=dict(median=k_meas[0], p16=k_meas[1], p84=k_meas[2],
                                n_used=k_meas[3], n_block=k_meas[4], at_kpc=8.0),
                n_tables=int(l.size), n_dist=int(D.DGRID.size),
                grid=dict(step_deg=STEP, coarse_deg=2 * STEP, fine_abs_b=FINE_ABS_B,
                          margin_deg=MARGIN, scan_reach_deg=G.scan_reach()),
                n_switched_to_marshall=int(used_nir.sum()), n_monotone_lifted=int(lifted),
                raw_cache=os.path.relpath(cache, ROOT))
    with open(PROV, "w") as f:
        json.dump(prov, f, indent=1)
    print(f"wrote {OUT} ({os.path.getsize(OUT)/1e6:.1f} MB) and {PROV}")


if __name__ == "__main__":
    main()
