"""
generateRomanBaseline.py

Generates RomanBaseline.dat: a synthetic per-visit observation log for Roman's Galactic Bulge
Time Domain Survey (GBTDS), F146 only. Columns, as read by src/run/inputs.cpp:

    ID  RA  Dec  l  b  time  sig5  field  layout

`field` is the GBTDS field index (0-4 the contiguous block west to east in the layout file's
order, 5 the Galactic Centre field) and `layout` the roll angle of the visit (0 = spring,
1 = autumn), which selects the detector layout src/surveys/footprints.cpp places at (l, b).
`time` is in days on the simulation clock (day 0 = first Rubin bulge visit), RA/Dec and l/b in
degrees, as in BulgeBaseline.dat.

Roman has no public per-visit database, so the visits are generated from the ROTAC 2025 season
and cadence design (the Roman-side analogue of readbaselineBulge.py).

Sources:
- ROTAC 2025 Final Report (arXiv:2505.10574): 6 fields; 6 high-cadence seasons (the first three
  and last three of 10 bulge seasons over the 5-yr mission), ~72 d each, F146 every 12.1 min;
  4 intervening low-cadence seasons.
- Roman GBTDS user documentation (roman-docs.stsci.edu, galactic-bulge-time-domain-survey).
- Field centres: the adopted layout, mtpenny/gbtds_optimizer field_layouts/
  gbtds_{spring,autumn}_2026.4.3.centers, vendored with its detector layout in
  Baseline/gbtds_layout/. Five contiguous fields at b = -1.400 plus the Galactic Centre field at
  b = -0.221; the l centres differ between the spring and autumn rolls.

Assumptions to check against the survey revision being cited:
1. LOW_CADENCE_DAYS = 5: the GBTDS design page says "repeated every five days"; the ROTAC 2025
   overguide quotes 3-day F146/F213. These may be different design iterations.
2. Season windows (SEASON_PATTERN) follow STScI's published schedule for the first two years,
   not an exact half-year split; the sun-angle constraint sets them.
3. MISSION_START_DAY places the mission inside the simulation's Tobs window; a modelling choice.
4. SIG5_PLACEHOLDER is one fixed depth per visit. errRomanM() in src/surveys/noise.cpp does not
   use this column; it is kept for symmetry with BulgeBaseline.dat.
5. No dithering is modelled. Detector gaps and the roll angle are (see `layout` above).
6. Only F146 (filter index 6 in config/parameters.h) is generated, not F087/F213.

Usage:
    python3 generateRomanBaseline.py [--mission-start DAYS]
Writes ./Baseline/RomanBaseline.dat and prints the row count; run
`python3 tools/sync_data_products.py` from the repo root afterwards to update `NlRoman` in
config/data_products.h.
"""

import os
import argparse
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis"))
from cparams import P  # noqa: E402

# Configuration

OUTPUT_PATH = "./Baseline/RomanBaseline.dat"

YEAR_DAYS          = P.year                              # config: physical_constants.h
MISSION_YEARS      = 5
SEASONS_PER_YEAR   = 2                                   # bulge visible ~twice/year
N_SEASONS          = MISSION_YEARS * SEASONS_PER_YEAR    # 10 total bulge seasons
# Bulge visibility windows from STScI's published schedule for the first two years of science
# operations (roman-docs.stsci.edu, roman-observations-in-the-first-two-years-of-science-operations):
#     high  2027-02-11 -> 2027-04-20   69 d
#     high  2027-08-15 -> 2027-10-25   72 d   (gap before: 117 d)
#     high  2028-02-11 -> 2028-04-21   71 d   (gap before: 109 d)
#     low   2028-08-16 -> 2028-10-24   70 d   (gap before: 117 d)
# The gaps alternate (~116 d spring->fall, ~108 d fall->spring): the sun-angle constraint sets the
# windows, not a half-year split. The (69, 72) pattern averages the design page's ~70.5 d.
#
# Annual pattern as (offset from that year's spring-window start [d], window length [d]).
SEASON_PATTERN = [(0.0, 69.0), (185.0, 72.0)]   # spring, fall

# High-cadence = first three and last three of the 10 seasons (ROTAC 2025)
HIGH_CADENCE_SEASONS = {0, 1, 2, N_SEASONS - 3, N_SEASONS - 2, N_SEASONS - 1}

HIGH_CADENCE_DAYS = 12.1 / (24.0 * 60.0)  # 12.1 minutes -> days (ROTAC 2025)
# GBTDS design page: "~1.5 hour observing units that are repeated every five days".
# (The ROTAC 2025 overguide quotes 3 days; see assumption 1 in the module docstring.)
LOW_CADENCE_DAYS  = 5.0

# Where the GBTDS starts on the simulation clock (day 0 = first Rubin bulge visit,
# MJD 61141.312 = 2026-04-11, set by readbaselineBulge.py). 306 d = 2027-02-11, the start of
# the first published high-cadence window, so SEASON_PATTERN reproduces the real spring/autumn
# dates that the two rolls depend on. Leaves ~0.8 yr of Rubin-only baseline before Roman.
# Override with --mission-start DAYS.
MISSION_START_DAY = 306.0
TOBS_DAYS = P.Tobs  # config/parameters.h Tobs; the C++ read CHECKs against it

# Adopted GBTDS fields, one list per roll (layout 0 = spring, 1 = autumn), from the vendored
# layout files (Baseline/gbtds_layout/README.md).
LAYOUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gbtds_layout")
CENTERS_FILES = ("gbtds_spring_2026.4.3.centers", "gbtds_autumn_2026.4.3.centers")


def read_centers(name):
    """[(l, b), ...] in the file's order: fields 1-5 (the contiguous block), then GC."""
    out = []
    with open(os.path.join(LAYOUT_DIR, name)) as f:
        next(f)                                  # header: field l b fixed
        for line in f:
            p = line.split()
            if p:
                out.append((float(p[1]), float(p[2])))
    if len(out) != 6:
        raise SystemExit(f"{name}: expected 6 GBTDS fields, found {len(out)}")
    return out


FIELDS_BY_LAYOUT = [read_centers(n) for n in CENTERS_FILES]

SIG5_PLACEHOLDER = 24.0  # placeholder; replace with a per-visit depth model if needed


def galactic_to_radec(l_deg, b_deg):
    """Galactic (l, b) [deg] to ICRS (RA, Dec) [deg]."""
    c = SkyCoord(l=l_deg * u.deg, b=b_deg * u.deg, frame="galactic")
    icrs = c.icrs
    return icrs.ra.deg, icrs.dec.deg


def build_season_windows(mission_start=MISSION_START_DAY):
    """[(season_index, start_day, end_day, cadence_days), ...] for all N_SEASONS seasons.

    `mission_start` is where season 0 begins on the simulation clock. The whole mission must
    land inside [0, Tobs], or the C++ read's CHECK(ro->tim[i] <= Tobs) aborts."""
    last_year, last_w = divmod(N_SEASONS - 1, len(SEASON_PATTERN))
    last_off, last_len = SEASON_PATTERN[last_w]
    last_end = mission_start + last_year * YEAR_DAYS + last_off + last_len
    if mission_start < 0.0 or last_end > TOBS_DAYS:
        raise SystemExit(
            f"MISSION_START_DAY={mission_start:.1f} puts the Roman mission at "
            f"[{mission_start:.1f}, {last_end:.1f}] d, outside the simulation window "
            f"[0, {TOBS_DAYS:.1f}] d. src/run/inputs.cpp's CHECK(ro->tim[i] <= Tobs) would "
            f"abort. Use a start day <= {TOBS_DAYS - (last_end - mission_start):.1f}.")
    windows = []
    for i in range(N_SEASONS):
        yr, w = divmod(i, len(SEASON_PATTERN))
        offset, length = SEASON_PATTERN[w]
        start = mission_start + yr * YEAR_DAYS + offset
        end = start + length
        cadence = HIGH_CADENCE_DAYS if i in HIGH_CADENCE_SEASONS else LOW_CADENCE_DAYS
        windows.append((i, start, end, cadence))
    return windows


def main():
    ap = argparse.ArgumentParser(description="Generate RomanBaseline.dat (F146 GBTDS visits).")
    ap.add_argument("--mission-start", type=float, default=MISSION_START_DAY, metavar="DAYS",
                    help="day on the simulation clock (0 = first Rubin bulge visit, "
                         f"2026-04-11) where Roman season 0 begins. Default {MISSION_START_DAY:g}.")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    windows = build_season_windows(args.mission_start)

    rows = []  # [ID, RA, Dec, l, b, time, sig5, field, layout]
    for season_idx, start, end, cadence in windows:
        # Seasons alternate spring/autumn, and so does the roll.
        layout = season_idx % len(SEASON_PATTERN)
        for field, (l, b) in enumerate(FIELDS_BY_LAYOUT[layout]):
            ra, dec = galactic_to_radec(l, b)
            t = start
            while t <= end:
                rows.append([0, ra, dec, l, b, t, SIG5_PLACEHOLDER, field, layout])
                t += cadence

    rows = np.array(rows)

    # matchVisibleEpochs() and src/run/inputs.cpp assume the file is sorted by time.
    order = np.argsort(rows[:, 5], kind="stable")
    rows = rows[order]
    rows[:, 0] = np.arange(len(rows))  # renumber IDs after sorting

    with open(OUTPUT_PATH, "w") as f:
        f.write("#ID  RA  Dec  l  b  time  sig5  field  layout\n")
        for row in rows:
            f.write(
                f"{int(row[0])}  {row[1]:.6f}  {row[2]:.6f}  "
                f"{row[3]:.6f}  {row[4]:.6f}  {row[5]:.6f}  {row[6]:.3f}  "
                f"{int(row[7])}  {int(row[8])}\n"
            )

    print(f"Wrote {len(rows)} Roman F146 visits to {OUTPUT_PATH}")
    print(f"NlRoman = {len(rows)}   <-- run `python3 tools/sync_data_products.py` from the repo root to update config/data_products.h")
    print(f"Time span: {rows[:, 5].min():.1f} to {rows[:, 5].max():.1f} days "
          f"(mission start = {args.mission_start:g} d on the Rubin clock)")
    print(f"Fields: {len(FIELDS_BY_LAYOUT[0])} per roll, 2 rolls | Seasons: {N_SEASONS} "
          f"({len(HIGH_CADENCE_SEASONS)} high-cadence, {N_SEASONS - len(HIGH_CADENCE_SEASONS)} low-cadence)")


if __name__ == "__main__":
    main()
