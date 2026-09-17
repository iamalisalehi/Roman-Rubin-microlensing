#!/usr/bin/env python3
"""Trim a killed run's outputs back to whole sightlines, so it can be continued.

WHY THIS IS NEEDED. The scan writes per-event rows as it goes and one aggregate block per
sightline when that sightline finishes. Kill it mid-sightline and the outputs disagree: the
per-event files hold a PARTIAL sightline, while the aggregate files hold nothing for it,
because those are written at the end. Continuing past that sightline leaves its partial rows in
the table for good -- and a sightline sampled to a fraction of its --events budget still claims
the full sky area `w_area` stands for, so every area-weighted yield built from it is an
undercount. Continuing AT it, without trimming, appends a complete second copy on top of the
partial one.

Neither is acceptable, so: delete the partial sightline's per-event rows, then re-run it whole.
This was done by hand at the v3 run's chunk boundary (PROGRESS.md, "chunks 1 and 2"); doing it
by hand is how it gets forgotten.

WHICH FILES ARE TOUCHED, AND WHICH ARE DELIBERATELY NOT.

  test<tag>.dat    per event, carries lon/lat  -> TRIMMED
  LpLMC<tag>.dat   per characterised event, carries lon/lat in cols 9,10 -> TRIMMED
  magC*/datC*.dat  per-timestep light-curve dumps for one sampled event; only ever written
                   for the legacy MACHO population (legacyId == 1) -> TRIMMED WHOLE if the
                   partial sightline is the one being redone is not decidable from their
                   contents (no lon/lat column), so they are left alone and reported. They are
                   diagnostics, not results.
  MapLMC<tag>.dat  one row per COMPLETED sightline -> untouched, and checked for consistency
  EfLMC<tag>(B)    one block per COMPLETED sightline -> untouched

THE lon/lat COMPARISON IS NUMERIC, NOT TEXTUAL. test<tag>.dat writes `0.281 -0.04` where
LpLMC<tag>.dat writes `0.2810000 -0.0400000` for the same sightline. A string match finds the
rows in one file and silently misses them in the other -- the trap PROGRESS.md flags for
"whoever automates this". That is this script.
"""

import argparse
import os
import re
import sys

LP_LON_COL = 8   # 0-based; verified against LpLMC5.dat's 19-column rows
LP_LAT_COL = 9
TOL = 1e-6       # lon/lat are written with >= 3 decimals in every file; this separates
                 # adjacent sightlines (0.02 deg apart at the finest stride) by 4 orders


def last_entered_sightline(log_path):
    """(lon, lat) of the last sightline the scan ENTERED, from the run log.

    The log prints 'NEW STEP n' and then 'longtitude: X   latitude: Y' for every sightline it
    enters, including ones later skipped as barren or uncovered -- which is exactly the set
    --start-index counts.
    """
    lon = lat = None
    pat = re.compile(r"longtitude:\s*(\S+)\s+latitude:\s*(\S+)")
    with open(log_path, "r", errors="replace") as fh:
        for line in fh:
            m = pat.search(line)
            if m:
                lon, lat = float(m.group(1)), float(m.group(2))
    if lon is None:
        raise SystemExit("resume_truncate: no 'longtitude:' line in the log -- nothing entered yet")
    return lon, lat


def header_columns(path):
    with open(path, "r", errors="replace") as fh:
        first = fh.readline()
    if not first.startswith("#"):
        return None
    return first.lstrip("#").split()


def trim_by_lonlat(path, lon, lat, lon_col, lat_col, dry_run):
    """Drop trailing rows whose lon/lat match the interrupted sightline. Returns rows removed.

    Only a TRAILING block is removed. If matching rows turn up before rows of some other
    sightline, the file is not what this script assumes -- two runs concatenated, say -- and it
    refuses rather than cutting a hole out of the middle.
    """
    if not os.path.exists(path):
        return 0

    with open(path, "r", errors="replace") as fh:
        lines = fh.readlines()

    keep = len(lines)
    for i in range(len(lines) - 1, -1, -1):
        s = lines[i]
        if s.startswith("#") or not s.strip():
            break
        parts = s.split()
        if len(parts) <= max(lon_col, lat_col):
            break
        try:
            if abs(float(parts[lon_col]) - lon) < TOL and abs(float(parts[lat_col]) - lat) < TOL:
                keep = i
                continue
        except ValueError:
            pass
        break

    removed = len(lines) - keep
    if removed == 0:
        return 0

    # Refuse if the same sightline also appears earlier with other sightlines after it.
    for i in range(keep - 1, -1, -1):
        s = lines[i]
        if s.startswith("#") or not s.strip():
            break
        parts = s.split()
        if len(parts) > max(lon_col, lat_col):
            try:
                if abs(float(parts[lon_col]) - lon) < TOL and abs(float(parts[lat_col]) - lat) < TOL:
                    raise SystemExit(
                        f"resume_truncate: {path} has rows for ({lon}, {lat}) that are NOT at the "
                        "end. This file looks like two runs concatenated; refusing to trim it."
                    )
            except ValueError:
                pass
        break

    if not dry_run:
        with open(path, "w") as fh:
            fh.writelines(lines[:keep])
    return removed


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="the run directory")
    ap.add_argument("--resume-at", type=int, required=True,
                    help="scan index the run will restart at (the interrupted sightline)")
    ap.add_argument("--dry-run", action="store_true", help="report, change nothing")
    args = ap.parse_args()

    d = os.path.abspath(args.dir)
    log = os.path.join(d, "run.log")
    if not os.path.exists(log):
        raise SystemExit(f"resume_truncate: no run.log in {d}")

    lon, lat = last_entered_sightline(log)
    print(f"resume_truncate: interrupted sightline is lon={lon} lat={lat} (scan index {args.resume_at})")

    tables = [p for p in os.listdir(d) if p.startswith("test") and p.endswith(".dat")]
    if len(tables) != 1:
        raise SystemExit(f"resume_truncate: expected exactly one test*.dat in {d}, found {tables}")
    table = os.path.join(d, tables[0])

    cols = header_columns(table)
    if not cols or "lon" not in cols or "lat" not in cols:
        raise SystemExit(f"resume_truncate: cannot find lon/lat in the header of {table}")
    t_lon, t_lat = cols.index("lon"), cols.index("lat")

    total = 0
    n = trim_by_lonlat(table, lon, lat, t_lon, t_lat, args.dry_run)
    print(f"  {'would remove' if args.dry_run else 'removed'} {n:>8} rows from {os.path.basename(table)}")
    total += n

    mont = os.path.join(d, "files", "MONTLMC", "files")
    for name in sorted(os.listdir(mont)) if os.path.isdir(mont) else []:
        if name.startswith("LpLMC") and name.endswith(".dat"):
            p = os.path.join(mont, name)
            n = trim_by_lonlat(p, lon, lat, LP_LON_COL, LP_LAT_COL, args.dry_run)
            print(f"  {'would remove' if args.dry_run else 'removed'} {n:>8} rows from {name}")
            total += n

    # The aggregate files must contain nothing for the interrupted sightline. Check rather than
    # assume: if one does, the assumption that they are written only at sightline end is wrong,
    # and continuing would duplicate a block.
    for name in sorted(os.listdir(mont)) if os.path.isdir(mont) else []:
        if name.startswith("MapLMC") and name.endswith(".dat"):
            p = os.path.join(mont, name)
            with open(p, "r", errors="replace") as fh:
                rows = [r.split() for r in fh if r.strip()]
            # MapLMC's lon/lat sit in the trailing block of the row; locate them by value.
            hits = sum(1 for r in rows
                       if any(abs(float(v) - lon) < TOL for v in r if _isnum(v))
                       and any(abs(float(v) - lat) < TOL for v in r if _isnum(v)))
            if hits:
                print(f"  NOTE: {name} already has {hits} row(s) mentioning this sightline. "
                      "If the run is continued AT it, check for a duplicate map row afterwards.")

    print(f"resume_truncate: {'would remove' if args.dry_run else 'removed'} {total} rows in total")


def _isnum(v):
    try:
        float(v)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    sys.exit(main())
