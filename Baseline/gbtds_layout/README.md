# Adopted GBTDS field layout (vendored)

Inputs that define where Roman looks. Read by `Baseline/generateRomanBaseline.py` (field
centres), by `src/surveys/footprints.cpp` (readGbtdsLayout) at start-up (detector rectangles), and by
`analysis/gbtds_geometry.py`. Copied unchanged from M. Penny's GBTDS field-layout tool:

- repository: https://github.com/mtpenny/gbtds_optimizer
- commit: `7c2e5e5ec25d41183117769f51a491c4012c57c6`, accessed 2026-10-01

| file | upstream path | content |
|---|---|---|
| `gbtds_spring_2026.4.3.centers` | `field_layouts/` | six field centres (l, b) [deg], spring roll |
| `gbtds_autumn_2026.4.3.centers` | `field_layouts/` | the same six, autumn roll |
| `sca_layout_spring.txt` | repo root | 18 detector (SCA) outlines, 5 vertices each: `sca dl db` [deg] |
| `sca_layout_fall.txt` | repo root | the same, autumn (fall) roll |

**How they are used** (the convention of upstream's README and `plotFields.py`): a detector's
sky outline is its offsets ADDED to the field centre's (l, b) -- no cos(b) factor, which at
|b| <= 1.4 deg changes the l-extent by < 3e-4 of itself.

**Checked on import (2026-10-01):** 18 rectangles per file, each axis-aligned in (l, b),
0.12491 deg on a side; 0.28085 deg^2 per field, 1.685 deg^2 for six (STScI: "1.7 deg^2");
the fall layout is the spring layout rotated by 180 deg (exactly, to 1e-6 deg). Field outline:
l -0.280..+0.204, b -0.394..+0.394 about the centre (spring; mirrored in autumn).

The five contiguous fields sit at b = -1.400 and the Galactic-centre field at b = -0.221, in
both seasons; the two rolls shift the l centres (spring -0.318..1.318 and GC +0.055; autumn
-0.468..1.168 and GC -0.095). These agree with the centres measured from the Aladin view of
the real tiles in Deviation 59. The superseded notional layout was upstream's
`field_layouts/layout_40395.centers` (five fields at b = -1.2, GC at (0, -0.125)).
