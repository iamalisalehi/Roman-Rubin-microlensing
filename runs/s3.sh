#!/usr/bin/env bash
# s3.sh -- Step S3: the short sample-event runs that feed the S2 plotter.
#
# WHY SHORT, AND WHY STOPPED BY HAND. A sample run exists to fill a handful of class quotas
# (samples/<pop>.spec), not to measure anything, so its per-sightline targets are small
# (--events 30 --lenses 10 --maxdraws 10000, against production's 300/50/5e4). The simulator
# has no "stop when the quotas are full" rule -- it scans all 1829 sightlines -- so once every
# class is filled the run is stopped with runctl. Every class needs Roman's footprint (the
# physics classes require ndw_R > 0, the others Roman epochs or a Roman detection), and in scan
# order the footprint sightlines sit at positions ~500-1100, so nothing fills before ~500.
#
# Roman coverage is IDENTICAL on every footprint sightline (50,401 epochs: one mission-wide
# visit list), so there is no "thinly covered" sightline to aim at. rubin_only and gap_filler
# are filled by TIME -- peaks near season edges, in the low-cadence seasons, or in the gaps --
# and are the classes most likely to fill last.
#
# The samples are illustrative. The sample run's RNG stream differs from production's (other
# targets), so these are not production events; nothing quantitative is read from them.
#
#   runs/s3.sh launch <bulge|bh|ns>   start that population's sample run -- refuses unless its
#                                     production run (runs/prod_<pop>_20260924) has finished
#   runs/s3.sh status                 per run: state, sightlines entered, class quotas filled
#   runs/s3.sh plot <bulge|bh|ns>     draw every dumped event with the S2 plotter

set -euo pipefail
here=$(dirname "$(readlink -f "$0")")
repo=$(dirname "$here")
cd "$here"

prod_done() {   # 0 if runs/prod_<pop>_20260924 is not running and entered all 1829 sightlines
    local s; s=$(./runctl.sh status "prod_$1_20260924")
    grep -q "^state: *not running" <<<"$s" && grep -q "^entered: *1829 " <<<"$s"
}

quotas() {      # "class kept/quota" per class, from the spec and the run's log
    local d=samples_$1_S3 spec
    spec=$d/samples/$1.spec
    awk '$1=="class"{print $2, $3}' "$spec" | while read -r c q; do
        k=0
        [[ -f $d/run.log ]] && k=$(grep -c "\[sample\] $c " "$d/run.log" || true)
        printf "%s %s/%s\n" "$c" "$k" "$q"
    done
}

case ${1:-} in
launch)
    pop=${2:?population}
    prod_done "$pop" || { echo "s3: prod_${pop}_20260924 has not finished -- not launching"; exit 1; }
    ./runctl.sh start "samples_${pop}_S3"
    ;;
status)
    for pop in bulge bh ns; do
        d=samples_${pop}_S3
        s=$(./runctl.sh status "$d")
        st=$(sed -n 's/^state: *\(.*\)/\1/p' <<<"$s" | cut -c1-40)
        n=$(sed -n 's/^entered: *\([0-9]*\).*/\1/p' <<<"$s")
        q=$(quotas "$pop" | tr '\n' ' ')
        full=$(quotas "$pop" | awk -F'[ /]' '$2<$3{o=1} END{print o?"":"  ALL FILLED -> runs/runctl.sh stop runs/'"$d"'"}')
        printf "%-6s %-40s entered %4s  %s%s\n" "$pop" "$st" "${n:-0}" "$q" "$full"
    done
    echo "disk used by samples: $(du -sh samples_*_S3/samples 2>/dev/null | awk '{printf "%s %s  ", $2, $1}')"
    ;;
plot)
    pop=${2:?population}
    d=$here/samples_${pop}_S3
    cd "$repo"
    .roman/bin/python analysis/s2_sample_lightcurves.py "$d/samples/$pop" \
        -o "figures/samples/s3_${pop}" --prov "$d/files/MONTLMC/files/run_provenance.txt"
    ;;
*)
    sed -n '2,25p' "$0"; exit 2
    ;;
esac
