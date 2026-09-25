#!/usr/bin/env bash
# yields_20260925.sh -- Step Y on the post-extinction-fix tables (runs/prod_*_20260924).
#
# ONE POPULATION AT A TIME. y1 holds every weightable draw in memory (bulge: ~12M rows), so the
# three are run in sequence, smallest first, each under analysis/memrun.py, which prints
# "[memrun] peak RSS" as the last line of its y1.log. An OOM kill is silent, so a population is
# finished only when its y1_yields.csv exists AND its y1.log ends with the [memrun] line.
#
#   figures/yield_20260925/<pop>/y1.log, y1_yields.{csv,md}    per population
#   figures/yield_20260925/y1_yields.csv                       the three concatenated
#   figures/yield_20260925/y3/                                 y3 on the concatenation
#
# Re-running skips any population already finished. Run from anywhere:
#   nohup runs/yields_20260925.sh > figures/yield_20260925/driver.log 2>&1 &

set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."
out=figures/yield_20260925
py=.roman/bin/python

for pop in bh ns bulge; do
    d=$out/$pop
    if [[ -f $d/y1_yields.csv ]] && tail -1 "$d/y1.log" | grep -q '^\[memrun\]'; then
        echo "$(date '+%F %T') $pop: already done"; continue
    fi
    mkdir -p "$d"
    echo "$(date '+%F %T') $pop: start"
    $py analysis/memrun.py analysis/y1_absolute_yield.py \
        --run "$pop=runs/prod_${pop}_20260924" -o "$d" > "$d/y1.log" 2>&1
    echo "$(date '+%F %T') $pop: done, $(tail -1 "$d/y1.log")"
done

{ head -1 "$out/bh/y1_yields.csv"
  for pop in bh ns bulge; do tail -n +2 "$out/$pop/y1_yields.csv"; done; } > "$out/y1_yields.csv"
$py analysis/y3_yield_vs_F.py "$out/y1_yields.csv" -o "$out/y3"
echo "$(date '+%F %T') all done"
