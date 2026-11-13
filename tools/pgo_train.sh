#!/usr/bin/env bash
# Training runs for the profile-guided build (`make pgo`).
#   usage: tools/pgo_train.sh BINARY TRAIN_DIR
# BINARY is the instrumented roman (absolute path). Each run gets its own directory under TRAIN_DIR, with
# the inputs symlinked from the repo root as make_rundir in pipeline/pipeline.sh does; the repo's own
# files/ is never written to. Fixed seeds. Prints one line per run with its CPU seconds; exits non-zero if
# any run does.
set -euo pipefail

BIN=${1:?usage: pgo_train.sh BINARY TRAIN_DIR}
TRAIN=${2:?usage: pgo_train.sh BINARY TRAIN_DIR}
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
[[ -x $BIN ]] || { echo "pgo_train: $BIN is not executable" >&2; exit 1; }

rm -rf "$TRAIN"
mkdir -p "$TRAIN"

make_rundir() {
    local d=$1
    mkdir -p "$d/files/MONTLMC/files" "$d/files/density" "$d/files/ext"
    ln -sfn "$ROOT/Baseline" "$d/Baseline"
    ln -sfn "$ROOT/CMD" "$d/CMD"
    ln -sfn "$ROOT/files/ext/ext_tables.dat" "$d/files/ext/ext_tables.dat"
    ln -sfn "$ROOT/files/sigma_roman.txt" "$d/files/sigma_roman.txt"
    ln -sfn "$ROOT/external" "$d/external"
    ln -sfn "$ROOT/files/ephemeris" "$d/files/ephemeris"
}

run() {   # run NAME ARGS...
    local name=$1; shift
    local d=$TRAIN/$name
    make_rundir "$d"
    local TIMEFORMAT="pgo_train: $name  CPU user %U + sys %S s  (wall %R s)"
    # shellcheck disable=SC2030
    time ( cd "$d" && "$BIN" "$@" > run.log 2>&1 ) \
        || { echo "pgo_train: run $name failed; see $d/run.log" >&2; tail -15 "$d/run.log" >&2; exit 1; }
}

run stub  --stub --population bulge --events 10 --lenses 2 --maxdraws 1000 --start-index 10 --end-index 11
run bulge --seed 42 --events 300 --lenses 50 --stride-roman 5 --start-index 795 --end-index 796
run bh    --seed 42 --events 300 --lenses 50 --stride-roman 5 --start-index 795 --end-index 796 --population bh
echo "pgo_train: all runs finished"
