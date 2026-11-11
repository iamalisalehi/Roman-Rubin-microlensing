#!/usr/bin/env bash
# Byte-identity and timing harness for refactors.
#
#   tests/identity_check.sh BASE_DIR NEW_DIR OUT_DIR [--quick] [--time N]
#
# BASE_DIR and NEW_DIR are built checkouts (each holds ./roman and ./fishertest). The same reference
# runs are made with each build, in separate run directories under OUT_DIR/{base,new}/<run>/, and every
# file a run wrote is compared. Exit 0 iff every compared file is identical.
#
# Reference runs (seed 42 throughout):
#   stub_776, stub_795        --stub, bulge, 10 events, 2 lenses, 1000 draws, one sightline
#   prod_bulge_795/_0/_1325   bulge, 300 events, 50 lenses, stride-roman 5, one sightline
#   prod_bh_795               the same for the BH population
#   pair_bulge_795            prod_bulge_795 plus --pair-satellite
#   dump_bulge_795            bulge 795, 40 events, 10 lenses, --dump-samples (samples/bulge.spec copy)
#   fisher_*                  fishertest in default, espl, --scale 1e-3, --eigen, --astro-variants,
#                             --finite-source, --lens-motion, --sweep, --sweep-astro modes
# --quick runs only the stubs and the fishertest modes.
# --time N then times N interleaved BASE/NEW pairs of prod_bulge_795 (CPU s = user + sys of the child,
# Python resource.getrusage(RUSAGE_CHILDREN)) and prints the median ratio NEW/BASE.
#
# Compared per run: the stdout log (log.txt), the exit status, every regular file the run wrote (the run
# directory's own files, files/MONTLMC/files/*, files/density/*, dumped samples) with cmp. Dropped before
# comparing: in run_provenance.txt AND in the stdout log (which echoes the provenance header) the lines
# starting "# git_commit" and "# built" (they differ between builds by design). The stdout log carries no
# wall-time, date or host lines (checked), so nothing else is filtered.
#
# Inputs are git-ignored and live only in the main checkout: env ROMAN_DATA (default below) is the data
# root; its inputs are symlinked into each run directory exactly as make_rundir in pipeline/pipeline.sh
# does. Two runs at a time at most.
set -euo pipefail

ROMAN_DATA=${ROMAN_DATA:-/home/ali/Documents/PhD/Offline_project/Roman}
PAR=2

usage() { echo "usage: $0 BASE_DIR NEW_DIR OUT_DIR [--quick] [--time N]" >&2; exit 2; }
[[ $# -ge 3 ]] || usage
BASE=$(realpath "$1"); NEW=$(realpath "$2"); mkdir -p "$3"; OUT=$(realpath "$3"); shift 3
QUICK=0; NTIME=0
while [[ $# -gt 0 ]]; do
    case $1 in
        --quick) QUICK=1 ;;
        --time)  [[ $# -ge 2 ]] || usage; NTIME=$2; shift ;;
        *) usage ;;
    esac
    shift
done
for b in "$BASE" "$NEW"; do
    [[ -x $b/roman && -x $b/fishertest ]] || { echo "ERROR: $b lacks a built roman and fishertest" >&2; exit 2; }
done
SPECDIR=$(realpath "$(dirname "$0")/../samples")
for f in Baseline CMD external files/ephemeris files/sigma_roman.txt files/ext/ext_tables.dat; do
    [[ -e $ROMAN_DATA/$f ]] || { echo "ERROR: $ROMAN_DATA/$f missing (set ROMAN_DATA)" >&2; exit 2; }
done

# A directory ./roman can run in: inputs linked, outputs real (make_rundir of pipeline/pipeline.sh).
make_rundir() {   # make_rundir DIR
    local d=$1
    rm -rf "$d"
    mkdir -p "$d/files/MONTLMC/files" "$d/files/density" "$d/files/ext"
    ln -sfn "$ROMAN_DATA/Baseline" "$d/Baseline"
    ln -sfn "$ROMAN_DATA/CMD" "$d/CMD"
    ln -sfn "$ROMAN_DATA/files/ext/ext_tables.dat" "$d/files/ext/ext_tables.dat"
    ln -sfn "$ROMAN_DATA/files/sigma_roman.txt" "$d/files/sigma_roman.txt"
    ln -sfn "$ROMAN_DATA/external" "$d/external"
    ln -sfn "$ROMAN_DATA/files/ephemeris" "$d/files/ephemeris"
}

# Reference runs: name|kind|args. kind roman -> ./roman ARGS; fisher -> ./fishertest ARGS.
STUB="--stub --population bulge --events 10 --lenses 2 --maxdraws 1000"
PROD="--seed 42 --events 300 --lenses 50 --stride-roman 5"
RUNS=(
    "stub_776|roman|$STUB --start-index 776 --end-index 777"
    "stub_795|roman|$STUB --start-index 795 --end-index 796"
)
FISH=(
    "fisher_default|fisher|"
    "fisher_espl|fisher|--source-model espl"
    "fisher_scale|fisher|--scale 1e-3"
    "fisher_eigen|fisher|--eigen"
    "fisher_astro_variants|fisher|--astro-variants"
    "fisher_finite_source|fisher|--finite-source"
    "fisher_lens_motion|fisher|--lens-motion"
    "fisher_sweep|fisher|--sweep"
    "fisher_sweep_astro|fisher|--sweep-astro"
)
if [[ $QUICK == 0 ]]; then
    RUNS+=(
        "prod_bulge_795|roman|--population bulge $PROD --start-index 795 --end-index 796"
        "prod_bulge_0|roman|--population bulge $PROD --start-index 0 --end-index 1"
        "prod_bulge_1325|roman|--population bulge $PROD --start-index 1325 --end-index 1326"
        "prod_bh_795|roman|--population bh $PROD --start-index 795 --end-index 796"
        "pair_bulge_795|roman|--population bulge $PROD --pair-satellite --start-index 795 --end-index 796"
        "dump_bulge_795|roman|--population bulge --seed 42 --events 40 --lenses 10 --stride-roman 5 --dump-samples bulge.spec --start-index 795 --end-index 796"
    )
fi
RUNS+=("${FISH[@]}")

# run_one BUILD NAME KIND ARGS... : in OUT/<label>/NAME, stdout to log.txt, exit status to exit_status.
run_one() {
    local label=$1 build=$2 name=$3 kind=$4; shift 4
    local d=$OUT/$label/$name rc=0
    make_rundir "$d"
    [[ $name == dump_* ]] && cp "$SPECDIR/bulge.spec" "$d/bulge.spec"   # its dir is relative: samples/bulge in the run dir
    if [[ $kind == roman ]]; then
        (cd "$d" && "$build/roman" "$@" >log.txt 2>log.err) || rc=$?
    else
        (cd "$d" && "$build/fishertest" "$@" >log.txt 2>log.err) || rc=$?
    fi
    echo "$rc" > "$d/exit_status"
}

# Apply the fixed filter to a file on stdout.
filtered() {   # filtered FILE
    case $(basename "$1") in
        run_provenance.txt|log.txt) grep -v -e '^# git_commit' -e '^# built' "$1" ;;
        *) cat "$1" ;;
    esac
}

NDIFF=0; NSAME=0
report() {   # report NAME REL STATUS
    printf '%-6s %s/%s\n' "$3" "$1" "$2"
    if [[ $3 == same ]]; then NSAME=$((NSAME + 1)); else NDIFF=$((NDIFF + 1)); fi
}

compare_run() {   # compare_run NAME
    local name=$1 b=$OUT/base/$1 n=$OUT/new/$1 rel
    local -a files
    mapfile -t files < <( { (cd "$b" && find . -type f -printf '%P\n'); (cd "$n" && find . -type f -printf '%P\n'); } | sort -u )
    for rel in "${files[@]}"; do
        [[ $rel == log.err ]] && continue   # stderr: compared only through exit_status
        if [[ -f $b/$rel && -f $n/$rel ]] && cmp -s <(filtered "$b/$rel") <(filtered "$n/$rel"); then
            report "$name" "$rel" same
        else
            report "$name" "$rel" DIFF
        fi
    done
}

# Run the reference runs, BASE and NEW of each, at most PAR at a time.
echo "== running $((${#RUNS[@]} * 2)) runs, $PAR at a time (output in $OUT) =="
rm -rf "$OUT/base" "$OUT/new"
active=0
for spec in "${RUNS[@]}"; do
    IFS='|' read -r name kind args <<<"$spec"
    for label in base new; do
        build=$BASE; [[ $label == new ]] && build=$NEW
        # shellcheck disable=SC2086
        run_one "$label" "$build" "$name" "$kind" $args &
        active=$((active + 1))
        if (( active >= PAR )); then wait -n || true; active=$((active - 1)); fi
    done
done
wait || true

echo "== comparing =="
for spec in "${RUNS[@]}"; do compare_run "${spec%%|*}"; done
echo "== identity: $NSAME same, $NDIFF DIFF =="

# --time N: N interleaved pairs of the bulge-795 production run, CPU seconds each.
timed() {   # timed BUILD DIR : prints user+sys of the roman child
    python3 - "$1" "$2" <<'PY'
import resource, subprocess, sys
build, d = sys.argv[1], sys.argv[2]
args = ("--population bulge --seed 42 --events 300 --lenses 50 --stride-roman 5 "
        "--start-index 795 --end-index 796").split()
with open(d + "/log.txt", "w") as out:
    subprocess.run([build + "/roman"] + args, cwd=d, stdout=out, stderr=subprocess.DEVNULL, check=True)
r = resource.getrusage(resource.RUSAGE_CHILDREN)
print(f"{r.ru_utime + r.ru_stime:.2f}")
PY
}
if (( NTIME > 0 )); then
    echo "== timing: $NTIME interleaved pairs of prod_bulge_795 (CPU s = user + sys) =="
    ratios=()
    for ((i = 1; i <= NTIME; i++)); do
        make_rundir "$OUT/time_base"; make_rundir "$OUT/time_new"
        tb=$(timed "$BASE" "$OUT/time_base"); tn=$(timed "$NEW" "$OUT/time_new")
        r=$(python3 -c "print(f'{$tn / $tb:.4f}')")
        ratios+=("$r")
        echo "pair $i: BASE $tb  NEW $tn  NEW/BASE $r"
    done
    med=$(python3 -c "import statistics, sys; print(f'{statistics.median(map(float, sys.argv[1:])):.4f}')" "${ratios[@]}")
    echo "median NEW/BASE: $med"
    rm -rf "$OUT/time_base" "$OUT/time_new"
fi

(( NDIFF == 0 ))
