#!/usr/bin/env bash
# p6 then p7 on the post-extinction-fix bh/ns runs (Step R tables), each under memrun.
# Finished when p67.log ends with "all done".
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/../.."
o=figures/prod_20260925
runs="--run bh=runs/prod_bh_20260924 --run ns=runs/prod_ns_20260924"
for s in p6_synergy_resolution p7_forecast_figures; do
    echo "$(date '+%F %T') $s start"
    # -o is a filename PREFIX: p6_synergy.png, p7_precision.png, ... inside $o
    .roman/bin/python analysis/memrun.py analysis/$s.py $runs -o $o/${s%%_*} > $o/$s.log 2>&1
    echo "$(date '+%F %T') $s done: $(tail -1 $o/$s.log)"
done
echo "$(date '+%F %T') all done"
