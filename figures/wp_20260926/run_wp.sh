#!/usr/bin/env bash
# The whitepaper's Results figures, regenerated on the post-extinction-fix bulge run
# (runs/prod_bulge_20260924, commit a5028fe). Replaces the *_v3 set made from the pre-fix
# 2026-09-06 run. Each script's stdout (its printed numbers) goes to <name>.log here.
#
# f2/f3/f4/h5 read test5_detJ.dat: the detected rows of test5.dat, header kept
#   awk 'NR==1 || /^#/ {print; next} $40==1' test5.dat > test5_detJ.dat      (detJ = column 40)
# Valid because every script here selects detections before using a row, and a row's event
# weight depends only on that row and its sightline's nsim (from the map file), never on the
# other rows -- so dropping the 8.3M undetected rows changes no number, only the memory.
set -uo pipefail
cd "$(dirname "$(readlink -f "$0")")/../.."
R=runs/prod_bulge_20260924
o=figures/wp_20260926
py=.roman/bin/python
W="--map $R/files/MONTLMC/files/MapLMC5.dat --log $R/run.log"
P="--provenance $R/files/MONTLMC/files/run_provenance.txt"
E=$R/test5_detJ.dat
run() { local n=$1; shift; echo "$(date '+%T') $n"; $py "$@" > $o/$n.log 2>&1 || echo "  FAILED: $n (see $o/$n.log)"; }
run f2_tE   analysis/f2_gap_filling.py $E --param tE  -o $o/f2_gap_filling.png     $P $W
run f2_piE  analysis/f2_gap_filling.py $E --param piE -o $o/f2_gap_filling_piE.png $P $W
run f3      analysis/f3_characterization_map.py $E -o $o/f3_characterization_map.png $P $W
run f4      analysis/f4_fisher_precision.py $E -o $o/f4_fisher_precision.png $P $W
run f4_all  analysis/f4_fisher_precision.py $E --scope all -o $o/f4_fisher_precision_all.png $P $W
run h5      analysis/h5_astrometric_shift.py $E -o $o/h5_astrometric_shift.png $P $W
# h5_astrometry_summary.py reads a column extract with the old names ndwL/ndwR (h5_crosscheck.py);
# built from the header, not from awk indices.
$py - <<'PY'
import sys; sys.path.insert(0, "analysis"); import romanlib as R
cols = ("tE piE tetE u0 Ml Dl Ds Vt lon lat ndw_L ndw_R detL detR detJ okA_J okA_L okA_R sigpiE_J "
        "sigpiE_L sigpiE_R sigtetE_J sigtetE_L sigtetE_R relMl_J relMl_L relMl_R okB_J okB_L okB_R "
        "condB_J condB_L condB_R w_area nepL_pk nepR_pk magb_F146").split()
d = R.load_events("runs/prod_bulge_20260924/test5_detJ.dat", usecols=cols)[cols]
d.rename(columns={"ndw_L": "ndwL", "ndw_R": "ndwR"}).to_csv(
    "runs/prod_bulge_20260924/h5_extract.dat", sep=" ", index=False)
PY
run h5sum   analysis/h5_astrometry_summary.py $R/h5_extract.dat -o $o/h5_astrometry_summary.png $W
run h3      analysis/h3_satellite_parallax.py $R/h3_pair.dat --out-prefix $o/h3 $W
echo "$(date '+%T') all done"
