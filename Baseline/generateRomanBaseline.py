"""
generateRomanBaseline.py

Generates RomanBaseline.dat: a synthetic per-visit observation log for Roman's
Galactic Bulge Time Domain Survey (GBTDS), F146 filter only. Output columns match
what Bulge_LSST.cpp's new Roman-baseline read block expects:

    ID  RA  Dec  l  b  time  sig5

(all reused from the existing BulgeBaseline.dat convention: `time` is in days
from a common survey t=0, RA/Dec in degrees, l/b Galactic in degrees.)

This is the Roman-side analog of readbaselineBulge.py, but generated from the
ROTAC 2025 GBTDS season/cadence design rather than queried from an LSST OpSim
database, since Roman has no equivalent public per-visit database to query —
its cadence is a fixed, known observing pattern.

-------------------------------------------------------------------------------
SOURCES (verify against these before treating any number below as final):
-------------------------------------------------------------------------------
- ROTAC 2025 Final Report and Recommendations
  https://assets.science.nasa.gov/content/dam/science/missions/rst/science/ROTAC-Report-20250424-v3.pdf
  (also arXiv:2505.10574): 6 fields; 6 high-cadence seasons (the first three and
  last three of 10 total bulge seasons over the 5-yr mission), each ~72 days,
  F146 every 12.1 minutes (67s exposures); 4 intervening low-cadence seasons.
- Roman GBTDS user documentation (roman-docs.stsci.edu/roman-community-defined-surveys/
  galactic-bulge-time-domain-survey): confirms 6 fields covering 1.7 deg^2 total,
  6 seasons "three early on, and three toward the end", each field observed every
  ~12 min in high-cadence seasons.
- Field centers: mtpenny/gbtds_optimizer, field_layouts/layout_40395.centers
  (M. Penny's GBTDS field-layout tool) — 5 contiguous fields at b=-1.2 deg plus
  1 Galactic Center field at (l,b)=(0,-0.125). This is a notional/community
  layout, NOT necessarily the final flight-adopted field centers — confirm
  against the ROTAC Fig. 4 field layout or later mission documentation.

-------------------------------------------------------------------------------
ASSUMPTIONS YOU SHOULD VALIDATE BEFORE TREATING THIS AS FINAL (flagged inline
with TODO(Ali) below too):
-------------------------------------------------------------------------------
1. LOW_CADENCE_DAYS = 3.0 — your own literature review cites "3-day F146/F213 in
   low-cadence seasons" per the ROTAC 2025 overguide, but the public Roman
   user-documentation site instead says low-cadence seasons use a "five-day
   cadence." These may refer to different design iterations. Pick the one that
   matches whichever ROTAC revision you're citing in the whitepaper.
2. SEASON_PERIOD_DAYS (~182.6 days between season starts, i.e. two visibility
   windows per year) is inferred from "visible ... for two 72-day stretches each
   spring and fall" (NASA GBTDS overview) combined with the ~111-day inter-season
   gap already in your lit review (365.2425/2 - 72 ~= 110.6 ~= 111). It is NOT an
   explicit ROTAC number — the real visibility windows are set by Roman's actual
   orbit/sun-angle constraints, not an exact half-year split. Treat this as a
   reasonable approximation, not a precise ephemeris.
3. MISSION_START_DAY — where Roman's 5-yr mission sits inside the simulation's
   10-yr Tobs window (LSST's survey length) is a modeling choice, not a physical
   constant. Defaults to 0 (Roman and LSST surveys start together); change this
   if your science case wants Roman offset relative to LSST's timeline.
4. SIG5_PLACEHOLDER — written as one fixed depth per visit. Roman is space-based
   so per-visit depth is far more uniform than LSST's (airmass/sky-brightness-
   dependent) depth, but this is still a placeholder, not a validated number.
   errRomanM() in Bulge_LSST.cpp currently ignores this column entirely (it only
   uses ro->mag/ro->err) — it's carried here for symmetry with BulgeBaseline.dat
   and in case you later want a per-visit-depth-dependent Roman error model.
5. No dithering, detector gaps, or per-field position angle are modeled — every
   visit to a field uses that field's exact center coordinates.
6. F087/F213 are NOT generated here — only F146 (filter index 6 in Bulge.h).
   Add a second generator (or a `filter` column + loop) if you extend the
   Fisher/light-curve code to use Roman's other bands.

-------------------------------------------------------------------------------
USAGE:
    python3 generateRomanBaseline.py
Writes ./Baseline/RomanBaseline.dat and prints the row count to set as
`NlRoman` in Bulge.h.
-------------------------------------------------------------------------------
"""

import os
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

# ============================================================================
# Configuration — every value here should be checked against the ROTAC 2025
# Final Report (arXiv:2505.10574) and your whitepaper before a real run.
# ============================================================================

OUTPUT_PATH = "./Baseline/RomanBaseline.dat"

YEAR_DAYS          = 365.2425
MISSION_YEARS      = 5
SEASONS_PER_YEAR   = 2                                   # bulge visible ~twice/year
N_SEASONS          = MISSION_YEARS * SEASONS_PER_YEAR    # 10 total bulge seasons
SEASON_LENGTH_DAYS = 72.0                                 # ROTAC 2025: six 72-day seasons
SEASON_PERIOD_DAYS = YEAR_DAYS / SEASONS_PER_YEAR         # ~182.6 days between season starts
                                                            # -> ~111-day gap between seasons
                                                            #    (see ASSUMPTION 2 above)

# High-cadence = first three and last three of the 10 seasons (ROTAC 2025)
HIGH_CADENCE_SEASONS = {0, 1, 2, N_SEASONS - 3, N_SEASONS - 2, N_SEASONS - 1}

HIGH_CADENCE_DAYS = 12.1 / (24.0 * 60.0)  # 12.1 minutes -> days (ROTAC 2025)
LOW_CADENCE_DAYS  = 3.0                    # TODO(Ali): confirm vs. "five-day" cited elsewhere

MISSION_START_DAY = 0.0  # where the 5-yr Roman mission sits inside the sim's Tobs window

# 6 GBTDS fields: 5 contiguous "bulge" fields (b=-1.2 deg) + 1 Galactic Center field.
# Source: mtpenny/gbtds_optimizer, field_layouts/layout_40395.centers (notional layout —
# confirm against ROTAC Fig. 4 / later mission docs before publication).
FIELDS_L_B = [
    (-0.417948, -1.200),
    (-0.008974, -1.200),
    ( 0.400000, -1.200),
    ( 0.808974, -1.200),
    ( 1.217948, -1.200),
    ( 0.000000, -0.125),  # Galactic Center field
]

SIG5_PLACEHOLDER = 24.0  # TODO(Ali): replace with a real per-visit depth model if needed


def galactic_to_radec(l_deg, b_deg):
    """Convert Galactic (l,b) [deg] to ICRS RA/Dec [deg], matching the convention
    already used in readbaselineBulge.py for BulgeBaseline.dat."""
    c = SkyCoord(l=l_deg * u.deg, b=b_deg * u.deg, frame="galactic")
    icrs = c.icrs
    return icrs.ra.deg, icrs.dec.deg


def build_season_windows():
    """Returns [(season_index, start_day, end_day, cadence_days), ...] for all
    N_SEASONS seasons, tagging each as high- or low-cadence per ROTAC 2025."""
    windows = []
    for i in range(N_SEASONS):
        start = MISSION_START_DAY + i * SEASON_PERIOD_DAYS
        end = start + SEASON_LENGTH_DAYS
        cadence = HIGH_CADENCE_DAYS if i in HIGH_CADENCE_SEASONS else LOW_CADENCE_DAYS
        windows.append((i, start, end, cadence))
    return windows


def main():
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    windows = build_season_windows()

    rows = []  # [ID, RA, Dec, l, b, time, sig5]
    for l, b in FIELDS_L_B:
        ra, dec = galactic_to_radec(l, b)
        for season_idx, start, end, cadence in windows:
            t = start
            while t <= end:
                rows.append([0, ra, dec, l, b, t, SIG5_PLACEHOLDER])
                t += cadence

    rows = np.array(rows)

    # Sort by time — matchVisibleEpochs()/main() in Bulge_LSST.cpp assume the
    # baseline file's `tim` column is pre-sorted, same convention as BulgeBaseline.dat.
    order = np.argsort(rows[:, 5])
    rows = rows[order]
    rows[:, 0] = np.arange(len(rows))  # renumber IDs after sorting

    with open(OUTPUT_PATH, "w") as f:
        f.write("#ID  RA  Dec  l  b  time  sig5\n")
        for row in rows:
            f.write(
                f"{int(row[0])}  {row[1]:.6f}  {row[2]:.6f}  "
                f"{row[3]:.6f}  {row[4]:.6f}  {row[5]:.6f}  {row[6]:.3f}\n"
            )

    print(f"Wrote {len(rows)} Roman F146 visits to {OUTPUT_PATH}")
    print(f"NlRoman = {len(rows)}   <-- set this constant in Bulge.h")
    print(f"Time span: {rows[:, 5].min():.1f} to {rows[:, 5].max():.1f} days")
    print(f"Fields: {len(FIELDS_L_B)} | Seasons: {N_SEASONS} "
          f"({len(HIGH_CADENCE_SEASONS)} high-cadence, {N_SEASONS - len(HIGH_CADENCE_SEASONS)} low-cadence)")


if __name__ == "__main__":
    main()
