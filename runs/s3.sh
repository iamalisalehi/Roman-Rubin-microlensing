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
# BATCHES. When a sample is not good enough, a further batch is drawn and the earlier events are
# KEPT (user, 2026-09-24). Batch "" is runs/samples_<pop>_S3; batch "b" is runs/samples_<pop>_S3b,
# and so on, each with its own spec (only the classes being redrawn, usually with tighter cuts)
# and `--start-index 500` so it starts at the footprint on a FRESH RNG stream -- the same flags
# from index 0 would redraw the same events. A batch's events keep their class folder but carry
# the batch letter in their name (astrometric_b001), so nothing is ever overwritten.
#
#   runs/s3.sh launch <pop> [batch]   start that sample run -- refuses unless the population's
#                                     production run (runs/prod_<pop>_20260924) has finished
#   runs/s3.sh status                 per sample run (every batch): state, entered, quotas filled
#   runs/s3.sh plot <pop> [batch]     draw every dumped event into figures/samples/<pop>/<class>/,
#                                     one self-contained folder per class (figures + data)

set -euo pipefail
here=$(dirname "$(readlink -f "$0")")
repo=$(dirname "$here")
cd "$here"

prod_done() {   # 0 if runs/prod_<pop>_20260924 is not running and entered all 1829 sightlines
    local s; s=$(./runctl.sh status "prod_$1_20260924")
    grep -q "^state: *not running" <<<"$s" && grep -q "^entered: *1829 " <<<"$s"
}

quotas() {      # "class kept/quota" per class, from the spec and the run's log; $2 = batch
    local d=samples_$1_S3${2:-} spec
    spec=$d/samples/$1.spec
    awk '$1=="class"{print $2, $3}' "$spec" | while read -r c q; do
        k=0
        [[ -f $d/run.log ]] && k=$(grep -c "\[sample\] $c " "$d/run.log" || true)
        printf "%s %s/%s\n" "$c" "$k" "$q"
    done
}

case ${1:-} in
launch)
    pop=${2:?population}; b=${3:-}
    prod_done "$pop" || { echo "s3: prod_${pop}_20260924 has not finished -- not launching"; exit 1; }
    d=samples_${pop}_S3$b
    ./runctl.sh start "$d"
    # runctl's `continue` resumes at run.startindex + sightlines entered; start wrote 0, which
    # is wrong for a batch launched with --start-index N.
    si=$(grep -o -- '--start-index [0-9]*' "$d/run.flags" | awk '{print $2}' || true)
    [[ -n $si ]] && echo "$si" > "$d/run.startindex"
    ;;
status)
    for d in samples_*_S3*/; do
        d=${d%/}; rest=${d#samples_}; pop=${rest%%_S3*}; b=${rest#*_S3}
        s=$(./runctl.sh status "$d")
        st=$(sed -n 's/^state: *\(.*\)/\1/p' <<<"$s" | cut -c1-40)
        n=$(sed -n 's/^entered: *\([0-9]*\).*/\1/p' <<<"$s")
        q=$(quotas "$pop" "$b" | tr '\n' ' ')
        full=$(quotas "$pop" "$b" | awk -F'[ /]' '$2<$3{o=1} END{print o?"":"  ALL FILLED -> runs/runctl.sh stop runs/'"$d"'"}')
        printf "%-8s %-40s entered %4s  %s%s\n" "$pop${b:+/$b}" "$st" "${n:-0}" "$q" "$full"
    done
    echo "disk used by samples: $(du -sh samples_*_S3/samples 2>/dev/null | awk '{printf "%s %s  ", $2, $1}')"
    ;;
plot)
    # One dedicated folder per population and class, each event self-contained in it:
    #   figures/samples/<pop>/run_provenance.txt
    #   figures/samples/<pop>/<class>/<class>_<id>.{pdf,png}              the figure
    #   figures/samples/<pop>/<class>/<class>_<id>_{epochs,model,params}.dat   its data
    # The .dat files are COPIED from the run directory, which stays the source of truth.
    # Safe to re-run while the sample run is still going: it redraws every event dumped so far.
    pop=${2:?population}; b=${3:-}
    d=$here/samples_${pop}_S3$b
    src=$d/samples/$pop
    prov=$d/files/MONTLMC/files/run_provenance.txt
    out=$repo/figures/samples/$pop
    cd "$repo"
    shopt -s nullglob
    events=("$src"/*_params.dat)
    [[ ${#events[@]} -gt 0 ]] || { echo "s3: no events dumped yet in $src"; exit 1; }
    mkdir -p "$out"
    cp "$prov" "$out/run_provenance${b:+_$b}.txt"
    for p in "${events[@]}"; do
        ev=$(basename "$p" _params.dat)          # e.g. roman_only_005
        cls=${ev%_*}                             # strip the 3-digit id: roman_only
        name=${cls}_$b${ev##*_}                  # batch b: roman_only_b005; batch "": unchanged
        mkdir -p "$out/$cls"
        for k in epochs model params; do cp "$src/${ev}_$k.dat" "$out/$cls/${name}_$k.dat"; done
        # Drawn from the copies, so the figure carries the batch-tagged name too.
        .roman/bin/python analysis/s2_sample_lightcurves.py "$out/$cls" --only "$name" \
            -o "$out/$cls" --prov "$prov" > /dev/null
        echo "s3: $pop/$cls/$name"
    done
    ;;
*)
    sed -n '2,25p' "$0"; exit 2
    ;;
esac
