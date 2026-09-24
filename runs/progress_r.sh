#!/usr/bin/env bash
# progress_r.sh -- one-line-per-run progress of the post-extinction-fix re-runs (Step R).
#
# "entered / 1829" is a poor progress measure on its own: a sightline outside Roman's footprint
# costs 8-18x less than one inside it (PROGRESS.md 5d), and the scan does not visit them in a
# uniform mix. The better yardstick is the SAME population's previous run: same seed, same
# sightline order, same flags, differing only in the extinction law. So this prints, per run,
# the sightlines entered, the event table's size against the previous run's final size, and
# the CPU time spent. The size fraction is approximate -- the extinction fix changes how many
# stars are visible and so how many rows each sightline writes -- but it tracks the work.
#
#   runs/progress_r.sh            (from the repo root or from runs/)

cd "$(dirname "$(readlink -f "$0")")"
TOTAL=1829
declare -A NEW=( [bulge]=prod_bulge_20260924 [bh]=prod_bh_20260924 [ns]=prod_ns_20260924 )
declare -A TAB=( [bulge]=test5.dat [bh]=testbh.dat [ns]=testns.dat )
# Previous final table sizes (bytes). bulge: the v3 run, which had no --pair-satellite, so its
# fraction runs slightly high; bh/ns: the 2026-09-17 runs, identical flags.
declare -A OLD=( [bulge]=$(stat -c %s ../../roman_runs/2026-09-06_v3_h7/test5.dat 2>/dev/null || echo 0)
                 [bh]=$(stat -c %s prod_bh_20260917/testbh.dat)
                 [ns]=$(stat -c %s prod_ns_20260917/testns.dat) )

printf "%-6s %-8s %-10s %-9s %-11s %s\n" run state entered table vs_prev cpu
for pop in bulge bh ns; do
    d=${NEW[$pop]}
    s=$(./runctl.sh status "$d" 2>/dev/null)
    state=$(sed -n 's/^state: *\([^ (]*\).*/\1/p' <<<"$s")
    cpu=$(sed -n 's/.*cpu \([0-9:-]*\),.*/\1/p' <<<"$s")
    n=$(sed -n 's/^entered: *\([0-9]*\).*/\1/p' <<<"$s")
    sz=$(stat -c %s "$d/${TAB[$pop]}" 2>/dev/null || echo 0)
    frac="n/a"; [[ ${OLD[$pop]} -gt 0 ]] && frac=$(awk -v a="$sz" -v b="${OLD[$pop]}" 'BEGIN{printf "%.1f%%", 100*a/b}')
    printf "%-6s %-8s %4s/%-5s %-9s %-11s %s\n" "$pop" "${state:-?}" "$n" "$TOTAL" \
        "$(numfmt --to=iec "$sz")" "$frac" "${cpu:--}"
done
echo "disk free: $(df -h . | awk 'NR==2{print $4}')   mem available: $(free -m | awk '/^Mem/{print $7}') MB   load: $(cut -d' ' -f1-3 /proc/loadavg)"
