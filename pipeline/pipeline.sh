#!/usr/bin/env bash
# pipeline.sh -- run the whole forecast pipeline (downloads -> preprocessing -> ./roman -> analysis)
# on a fresh clone, on a laptop or on a Slurm cluster.
#
#   pipeline/pipeline.sh CONFIG STAGE [args]
#
#   check      what is present / missing / stale; changes nothing
#   setup      create the Python venv, install requirements.txt, build GSL in deps/gsl if needed (GSL=)
#   fetch      MIST bolometric-correction tables, the LSSTCam focal-plane map, the OpSim database
#              (if OPSIM_URL is set), the dust maps (only if EXT_TABLES=build)
#   prep       star lists, Rubin and Roman visit lists, extinction tables, config/data_products.h, ./roman
#   sim        run ./roman for every population, split into chunks of CHUNK_SIZE sightlines
#              (SCHEDULER=local: a pool of JOBS workers; SCHEDULER=slurm: array jobs + a merge/analyze job)
#   chunk POP K   run one chunk (internal: what both schedulers call)
#   merge      concatenate the chunks of each population into one run directory
#   analyze    the analysis scripts named in ANALYSES, logs next to their outputs
#   all        check, (setup), fetch, prep, sim, merge, analyze   (with slurm: stops after submitting)
#   status     chunks done / running / pending, CPU used, output sizes
#
# Every stage is idempotent: a step whose inputs (path, size, mtime) and options are unchanged since it
# last succeeded is skipped. FORCE=1 redoes it. See pipeline/README.md and pipeline/config.example.sh.

set -uo pipefail

SELF=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/$(basename "${BASH_SOURCE[0]}")
ROOT=$(dirname "$(dirname "$SELF")")
STARTDIR=$PWD
FORCE=${FORCE:-0}
export FORCE

die()  { echo "pipeline: ERROR: $*" >&2; exit 1; }
say()  { echo "$*"; }
warn() { echo "pipeline: warning: $*" >&2; }

usage() {
    sed -n '2,/^# last succeeded/p' "$SELF" | sed 's/^# \{0,1\}//' >&2
    exit 2
}

[[ $# -ge 2 ]] || usage
CONFIG_ARG=$1; STAGE=$2; shift 2
case $CONFIG_ARG in /*) CONFIG=$CONFIG_ARG ;; *) CONFIG=$STARTDIR/$CONFIG_ARG ;; esac
[[ -f $CONFIG ]] || die "config file not found: $CONFIG_ARG (copy pipeline/config.example.sh and edit it)"
case $ROOT in *[[:space:]]*) die "the repository path contains whitespace ($ROOT); move the clone" ;; esac
cd "$ROOT" || die "cannot cd to $ROOT"

# ------------------------------------------------------------------------------------------------
# Configuration: defaults, then the user's file, then validation.
# ------------------------------------------------------------------------------------------------
# shellcheck disable=SC1090
source "$CONFIG"

: "${BESANCON_CATALOGUE:=CMD/Besancon/bos10}"
: "${OPSIM_DB:=Baseline/baseline_v5.1.0_10yrs.db}"
: "${OPSIM_URL:=https://s3df.slac.stanford.edu/data/rubin/sim-data/sims_featureScheduler_runs5.1/baseline/baseline_v5.1.0_10yrs.db}"
: "${MIST_URL_BASE:=https://mist.science/BC_tables/v2}"
: "${EXT_TABLES:=}"
: "${CATALOGUE_DWARFS:=besancon}"
: "${CATALOGUE_FILL:=none}"
: "${ROMAN_MISSION_START:=306}"
: "${RUN_ROOT:=runs/pipeline_run}"
: "${POPULATIONS:=bulge bh ns}"
: "${STRIDE:=10}"
: "${STRIDE_ROMAN:=5}"
: "${EVENTS:=300}"
: "${LENSES:=50}"
: "${MAXDRAWS:=}"
: "${PAIR_SATELLITE:=1}"
: "${SEED:=42}"
: "${STUB:=0}"
: "${EXTRA_FLAGS:=}"
: "${CHUNK_SIZE:=20}"
: "${SCHEDULER:=local}"
: "${SLURM_OPTS:=}"
: "${SLURM_OPTS_FINISH:=$SLURM_OPTS}"
: "${SUBMIT:=1}"
: "${KEEP_CHUNKS:=1}"
: "${ANALYSES:=yields figures}"
: "${PYTHON:=.roman/bin/python}"
: "${GSL:=auto}"
: "${GSL_URL:=https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz}"
if [[ -z ${JOBS:-} ]]; then JOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1); fi

absp() { case $1 in /*) printf '%s\n' "$1" ;; *) printf '%s\n' "$ROOT/$1" ;; esac; }
P_BESANCON=$(absp "$BESANCON_CATALOGUE")
P_OPSIM=$(absp "$OPSIM_DB")
P_RUN=$(absp "$RUN_ROOT")
PY=$(absp "$PYTHON")
P_EXT=""; [[ -n $EXT_TABLES && $EXT_TABLES != build ]] && P_EXT=$(absp "$EXT_TABLES")

STAMPS=$ROOT/pipeline/.stamps          # prep/fetch stamps: one per data configuration, so one per clone
LOGS=$ROOT/pipeline/logs
RUN_STAMPS=$P_RUN/.stamps
BIN_SRC=$ROOT/roman

is_uint() { [[ $1 =~ ^[0-9]+$ ]]; }

validate_config() {
    local bad=0 e v a
    err() { echo "pipeline: config error: $*" >&2; bad=1; }
    for e in STRIDE STRIDE_ROMAN EVENTS LENSES CHUNK_SIZE JOBS; do
        { is_uint "${!e}" && (( ${!e} >= 1 )); } || err "$e must be a positive integer (got '${!e}')"
    done
    is_uint "$SEED" || err "SEED must be a non-negative integer (got '$SEED')"
    [[ -z $MAXDRAWS || $MAXDRAWS =~ ^[0-9]+(\.[0-9]+)?([eE][0-9]+)?$ ]] || err "MAXDRAWS must be a number or empty (got '$MAXDRAWS')"
    if is_uint "$STRIDE" && is_uint "$STRIDE_ROMAN" && (( STRIDE_ROMAN >= 1 )) && (( STRIDE % STRIDE_ROMAN != 0 )); then
        err "STRIDE_ROMAN ($STRIDE_ROMAN) must divide STRIDE ($STRIDE)"
    fi
    for v in PAIR_SATELLITE STUB SUBMIT KEEP_CHUNKS; do
        [[ ${!v} == 0 || ${!v} == 1 ]] || err "$v must be 0 or 1 (got '${!v}')"
    done
    case $GSL in auto|system|build) ;; *) err "GSL must be auto, system or build (got '$GSL')" ;; esac
    [[ $SCHEDULER == local || $SCHEDULER == slurm ]] || err "SCHEDULER must be local or slurm (got '$SCHEDULER')"
    case $CATALOGUE_DWARFS in besancon|empirical) ;; *) err "CATALOGUE_DWARFS must be besancon or empirical (got '$CATALOGUE_DWARFS')" ;; esac
    case $CATALOGUE_FILL in none|kroupa|koshimoto) ;; *) err "CATALOGUE_FILL must be none, kroupa or koshimoto (got '$CATALOGUE_FILL')" ;; esac
    [[ $ROMAN_MISSION_START =~ ^-?[0-9]+(\.[0-9]+)?$ ]] || err "ROMAN_MISSION_START must be a number of days (got '$ROMAN_MISSION_START')"
    [[ -n $POPULATIONS ]] || err "POPULATIONS is empty"
    for v in $POPULATIONS; do [[ $v =~ ^[A-Za-z0-9_-]+$ ]] || err "bad population name '$v'"; done
    for a in $ANALYSES; do
        case $a in yields|figures) ;; *) err "ANALYSES: unknown entry '$a' (yields, figures)" ;; esac
    done
    case $P_RUN$P_BESANCON$P_OPSIM in *[[:space:]]*) err "RUN_ROOT, BESANCON_CATALOGUE and OPSIM_DB must not contain whitespace" ;; esac
    [[ $SCHEDULER == slurm && -z $SLURM_OPTS ]] && warn "SCHEDULER=slurm with empty SLURM_OPTS: the jobs will get your cluster's defaults (time limit! memory!)"
    (( bad == 0 )) || exit 2
}

# The binary's own opinion of the population names (it lists the valid ones when it rejects one).
# --stride 0 makes a VALID population stop at the stride check, before it reads any data file.
check_populations() {   # check_populations BINARY
    local pop out
    for pop in $POPULATIONS; do
        out=$("$1" --population "$pop" --stride 0 2>&1 </dev/null)
        if [[ $out == *"unknown population"* ]]; then
            echo "$out" | head -20 >&2
            die "population '$pop' (POPULATIONS in $CONFIG_ARG) is not one the binary knows"
        fi
    done
}

# ------------------------------------------------------------------------------------------------
# Stamps: a step is current when its stamp equals (options + path/size/mtime of every input) and all
# its outputs exist.
# ------------------------------------------------------------------------------------------------
stat_sm() { stat -L -c '%s %Y' "$1" 2>/dev/null || stat -L -f '%z %m' "$1" 2>/dev/null || echo "? ?"; }
md5f() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }

file_sig() {
    local f=$1
    if [[ -d $f ]]; then
        local g
        while IFS= read -r g; do echo "$g $(stat_sm "$g")"; done < <(find "$f" -type f 2>/dev/null | LC_ALL=C sort)
    elif [[ -e $f ]]; then
        echo "$f $(stat_sm "$f")"
    else
        echo "$f MISSING"
    fi
}

stamp_text() {   # stamp_text OPTS INPUT...
    local opts=$1; shift
    echo "opts $opts"
    local f; for f in "$@"; do file_sig "$f"; done
}

all_exist() { local f; for f in "$@"; do [[ -e $f ]] || return 1; done; return 0; }

# step_state STAMP OPTS "INPUTS" "OUTPUTS"  ->  current | missing | stale
step_state() {
    local st=$1 opts=$2 ins=$3 outs=$4
    # shellcheck disable=SC2086
    all_exist $outs || { echo missing; return; }
    [[ -f $st ]] || { echo stale; return; }
    # shellcheck disable=SC2086
    [[ "$(cat "$st")" == "$(stamp_text "$opts" $ins)" ]] && echo current || echo stale
}

fmt_dur() { local s=$1; if (( s >= 3600 )); then printf '%dh%02dm' $((s/3600)) $((s%3600/60)); elif (( s >= 60 )); then printf '%dm%02ds' $((s/60)) $((s%60)); else printf '%ds' "$s"; fi; }

# run_step NAME STAMP LOG OPTS "INPUTS" "OUTPUTS" CMD...   (inputs/outputs: space-separated paths)
run_step() {
    local name=$1 st=$2 log=$3 opts=$4 ins=$5 outs=$6; shift 6
    local state
    state=$(step_state "$st" "$opts" "$ins" "$outs")
    if [[ $FORCE != 1 && $state == current ]]; then
        say "  skip   $name (up to date)"; return 0
    fi
    [[ $state == stale && -f $st ]] && say "  stale  $name (inputs or options changed since it last ran)"
    mkdir -p "$(dirname "$st")" "$(dirname "$log")"
    local want t0=$SECONDS
    # shellcheck disable=SC2086
    want=$(stamp_text "$opts" $ins)
    rm -f "$st"
    say "  run    $name   (log: ${log#"$ROOT"/})"
    if "$@" >"$log" 2>&1; then
        printf '%s\n' "$want" > "$st"
        say "  done   $name in $(fmt_dur $((SECONDS - t0)))"
        return 0
    fi
    say "  FAILED $name after $(fmt_dur $((SECONDS - t0))); last lines of ${log#"$ROOT"/}:"
    tail -12 "$log" | sed 's/^/         | /'
    return 1
}

# ------------------------------------------------------------------------------------------------
# The prep/fetch steps, as data, so `check` and `prep` cannot disagree about what they need.
# step_def NAME sets S_STAMP S_OPTS S_INS S_OUTS S_FN.
# ------------------------------------------------------------------------------------------------
COMPONENT_FILES="CMD/components/thin_disk.dat CMD/components/bulge.dat CMD/components/thick_disk.dat CMD/components/halo.dat CMD/components/provenance.txt"
EEM=CMD/empirical/EEM_dwarf_UBVIJHK_colors_Teff.txt
PREP_STEPS="catalogue lens_ml rubin_visits roman_visits fovmap extinction sync build"

step_def() {
    S_STAMP=$STAMPS/$1.stamp; S_OPTS=""; S_INS=""; S_OUTS=""; S_FN=""
    case $1 in
        fovmap)     S_FN=do_fovmap; S_INS=Baseline/lsstcam_fov/export_fov_map.py; S_OUTS=Baseline/lsstcam_fov/fov_map.txt ;;
        opsim)      S_FN=do_opsim;  S_OUTS=$P_OPSIM; S_STAMP="" ;;
        dustmaps)   S_FN=do_dustmaps; S_OUTS="dustmaps/decaps/decaps_mean.h5 dustmaps/marshall/marshall.h5"; S_STAMP="" ;;
        catalogue)
            S_FN=do_catalogue; S_OPTS="dwarfs=$CATALOGUE_DWARFS fill=$CATALOGUE_FILL"
            S_INS="$P_BESANCON CMD/BolometricCorrection.py CMD/Rubin CMD/Roman $EEM"; S_OUTS=$COMPONENT_FILES ;;
        lens_ml)
            S_FN=do_lens_ml; S_OPTS="dwarfs=$CATALOGUE_DWARFS"
            S_INS="$P_BESANCON CMD/lens_ml_table.py CMD/BolometricCorrection.py $EEM"; S_OUTS=CMD/components/lens_ml.dat ;;
        rubin_visits)
            S_FN=do_rubin_visits; S_OPTS="db=$P_OPSIM"
            S_INS="$P_OPSIM Baseline/readbaselineBulge.py analysis/gbtds_geometry.py Baseline/gbtds_layout config/parameters.h"
            S_OUTS=Baseline/BulgeBaseline.dat ;;
        roman_visits)
            S_FN=do_roman_visits; S_OPTS="mission_start=$ROMAN_MISSION_START"
            S_INS="Baseline/generateRomanBaseline.py Baseline/gbtds_layout config/parameters.h"; S_OUTS=Baseline/RomanBaseline.dat ;;
        extinction)
            S_FN=do_extinction; S_OPTS="ext=${EXT_TABLES:-UNSET}"
            if [[ $EXT_TABLES == build ]]; then S_INS="maps.py analysis/dustref.py analysis/gbtds_geometry.py"; else S_INS=$P_EXT; fi
            S_OUTS=files/ext/ext_tables.dat ;;
        sync)
            S_FN=do_sync; S_INS="$COMPONENT_FILES CMD/components/lens_ml.dat Baseline/BulgeBaseline.dat Baseline/RomanBaseline.dat files/sigma_roman.txt tools/sync_data_products.py"
            S_OUTS=config/data_products.h ;;
        build)      S_FN=do_build; S_OUTS=roman; S_STAMP="" ;;
        *) die "internal: unknown step $1" ;;
    esac
}

# State of a step for `check`. Steps without a stamp (opsim, dustmaps, build) look at their outputs.
state_of() {
    step_def "$1"
    if [[ $1 == build ]]; then
        [[ -x roman ]] || { echo missing; return; }
        make -q roman >/dev/null 2>&1 && echo current || echo stale
    elif [[ -z $S_STAMP ]]; then
        # shellcheck disable=SC2086
        all_exist $S_OUTS && echo current || echo missing
    else
        step_state "$S_STAMP" "$S_OPTS" "$S_INS" "$S_OUTS"
    fi
}

run_named_step() {   # run_named_step NAME
    step_def "$1"
    if [[ -z $S_STAMP ]]; then
        state=$(state_of "$1")
        if [[ $FORCE != 1 && $state == current && $1 != build ]]; then say "  skip   $1 (present)"; return 0; fi
        "$S_FN"; return
    fi
    run_step "$1" "$S_STAMP" "$LOGS/$1.log" "$S_OPTS" "$S_INS" "$S_OUTS" "$S_FN"
}

# ------------------------------------------------------------------------------------------------
# Step bodies
# ------------------------------------------------------------------------------------------------
need_py() { [[ -x $PY ]] || die "Python not found: $PY -- run: pipeline/pipeline.sh $CONFIG_ARG setup"; }

do_mist() {   # do_mist SYSTEM(LSST|Roman) DIR(Rubin|Roman)
    local sys=$1 dir=$2 tmp n
    tmp=$(mktemp -d "${TMPDIR:-/tmp}/mist.XXXXXX") || return 1
    curl -fL --retry 3 -sS -o "$tmp/$sys.txz" "$MIST_URL_BASE/$sys.txz" || { rm -rf "$tmp"; return 1; }
    mkdir -p "$tmp/x" && tar -xJf "$tmp/$sys.txz" -C "$tmp/x" || { rm -rf "$tmp"; return 1; }
    n=$(find "$tmp/x" -type f -name "feh*_afe*.$sys" | wc -l)
    if (( n != 75 )); then echo "expected 75 feh*_afe*.$sys files in the archive, found $n"; rm -rf "$tmp"; return 1; fi
    mkdir -p "CMD/$dir"
    find "$tmp/x" -type f -name "feh*_afe*.$sys" -exec cp {} "CMD/$dir/" \;
    rm -rf "$tmp"
    echo "$n files -> CMD/$dir/"
}
# The MIST tables are decided by the file count, not by a stamp (nothing else writes there).
mist_state() { local sys=$1 dir=$2 n; n=$(find "CMD/$dir" -maxdepth 1 -name "feh*_afe*.$sys" 2>/dev/null | wc -l); (( n == 75 )) && echo current || echo missing; }

do_fovmap() { need_py; "$PY" Baseline/lsstcam_fov/export_fov_map.py; }

do_opsim() {
    [[ -n $OPSIM_URL ]] || die "OpSim database $P_OPSIM is missing and OPSIM_URL is empty: download a baseline .db and set OPSIM_DB"
    mkdir -p "$(dirname "$P_OPSIM")"
    say "  download $OPSIM_URL (~0.8 GB) -> $P_OPSIM"
    curl -fL --retry 3 -C - -o "$P_OPSIM.part" "$OPSIM_URL" && mv "$P_OPSIM.part" "$P_OPSIM"
}

do_dustmaps() {
    need_py
    say "  downloading the DECaPS mean map (7 GB) and the Marshall map (5 MB) into dustmaps/ ..."
    mkdir -p "$ROOT/dustmaps/decaps" "$ROOT/dustmaps/marshall"
    "$PY" - <<PY
from dustmaps.config import config
config["data_dir"] = "$ROOT/dustmaps"
import dustmaps.decaps, dustmaps.marshall
dustmaps.decaps.fetch(mean_only=True, silence_warnings=True)
dustmaps.marshall.fetch()
PY
}

do_catalogue() {
    need_py
    (cd CMD && "$PY" BolometricCorrection.py --input "$P_BESANCON" --dwarfs "$CATALOGUE_DWARFS" \
        --fill-floor "$CATALOGUE_FILL" --out-dir components)
}
do_lens_ml() {
    need_py
    "$PY" CMD/lens_ml_table.py --input "$P_BESANCON" --dwarfs "$CATALOGUE_DWARFS" --out "$ROOT/CMD/components/lens_ml.dat"
}
do_rubin_visits()  { need_py; (cd Baseline && "$PY" readbaselineBulge.py --db "$P_OPSIM" --no-plots); }
do_roman_visits()  { need_py; "$PY" Baseline/generateRomanBaseline.py --mission-start "$ROMAN_MISSION_START"; }
do_extinction() {
    if [[ $EXT_TABLES == build ]]; then
        need_py
        all_exist dustmaps/decaps/decaps_mean.h5 dustmaps/marshall/marshall.h5 \
            || { echo "dust maps not downloaded: run the fetch stage with EXT_TABLES=build"; return 1; }
        "$PY" maps.py
    else
        mkdir -p files/ext
        if [[ $(cd "$(dirname "$P_EXT")" && pwd -P)/$(basename "$P_EXT") != "$ROOT/files/ext/ext_tables.dat" ]]; then
            cp "$P_EXT" files/ext/ext_tables.dat || return 1
            [[ -f $(dirname "$P_EXT")/ext_provenance.json ]] && cp "$(dirname "$P_EXT")/ext_provenance.json" files/ext/
        fi
        echo "extinction table: $(wc -c < files/ext/ext_tables.dat) bytes"
    fi
    return 0
}
do_sync() { need_py; "$PY" tools/sync_data_products.py; }
do_build() {
    local log=$LOGS/build.log
    mkdir -p "$LOGS"
    if [[ $FORCE != 1 ]] && make -q roman >/dev/null 2>&1; then say "  skip   build (./roman up to date)"; return 0; fi
    say "  run    build   (log: ${log#"$ROOT"/})"
    local t0=$SECONDS
    if make -j "$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)" roman >"$log" 2>&1; then
        say "  done   build in $(fmt_dur $((SECONDS - t0))); compiler warnings: $(grep -c 'warning:' "$log")"
        return 0
    fi
    say "  FAILED build; last lines of ${log#"$ROOT"/}:"; tail -15 "$log" | sed 's/^/         | /'
    return 1
}

# ------------------------------------------------------------------------------------------------
# Stages: check
# ------------------------------------------------------------------------------------------------
N_REQ=0; N_PROD=0
item() {   # item CLASS(req|prod) STATUS LABEL DETAIL
    printf '  %-8s %-26s %s\n' "$2" "$3" "$4"
    if [[ $2 != ok && $2 != current && $2 != note ]]; then
        if [[ $1 == req ]]; then N_REQ=$((N_REQ + 1)); else N_PROD=$((N_PROD + 1)); fi
    fi
}
have() { command -v "$1" >/dev/null 2>&1; }

check_tools() {
    say "Tools (install them yourself):"
    local t
    for t in g++ make awk curl tar xz; do
        if have "$t"; then item req ok "$t" "$(command -v "$t")"; else item req MISSING "$t" "not on PATH"; fi
    done
    if [[ $GSL != build ]] && gsl_system; then
        item req ok GSL "system headers found"
    elif gsl_local; then
        item req ok GSL "local build in deps/gsl"
    elif [[ $GSL == system ]]; then
        item req MISSING GSL "gsl/gsl_matrix.h not found: install libgsl-dev (apt) / gsl (brew, module load gsl), or set GSL=auto"
    else
        item prod missing GSL "not installed: setup builds a static copy in deps/gsl"
    fi
    if have python3 || [[ -x $PY ]]; then item req ok python3 "$(command -v python3 || echo "$PY")"; else item req MISSING python3 "needed to create the venv"; fi
    if [[ $SCHEDULER == slurm ]]; then
        if have sbatch; then item req ok sbatch "$(command -v sbatch)"; else item req MISSING sbatch "SCHEDULER=slurm but sbatch is not on PATH (use SUBMIT=0 to only write the scripts)"; fi
    fi
    have git || item req note git "not found: the binary will be named 'nogit' in RUN_ROOT/bin"
}

check_python() {
    say "Python environment (pipeline setup builds it):"
    if [[ ! -x $PY ]]; then item prod MISSING "$PYTHON" "no such interpreter: run the setup stage"; return; fi
    local miss
    miss=$("$PY" - <<'PY' 2>/dev/null
import importlib
bad = []
for m in ("numpy", "scipy", "pandas", "astropy", "matplotlib", "dustmaps"):
    try:
        importlib.import_module(m)
    except Exception:
        bad.append(m)
print(" ".join(bad))
PY
)
    if [[ -n $miss ]]; then item prod MISSING "$PYTHON" "missing modules: $miss -- run the setup stage"
    else item prod ok "$PYTHON" "$("$PY" -V 2>&1)"; fi
}

check_inputs() {
    say "Inputs you provide (or that fetch downloads):"
    if [[ -f $P_BESANCON ]]; then
        local py=python3; [[ -x $PY ]] && py=$PY
        if msg=$("$py" pipeline/validate_inputs.py besancon "$P_BESANCON" 2>&1); then item req ok Besancon "$BESANCON_CATALOGUE: $msg"
        else item req MISSING Besancon "$msg"; fi
    else
        item req MISSING Besancon "$BESANCON_CATALOGUE not found. Not downloadable by script (web form): see pipeline/README.md, CMD/Besancon/RERUN_bos10.md"
    fi
    if [[ -f $P_OPSIM ]]; then
        local py=python3; [[ -x $PY ]] && py=$PY
        if msg=$("$py" pipeline/validate_inputs.py opsim "$P_OPSIM" 2>&1); then item req ok "OpSim db" "$OPSIM_DB: $(echo "$msg" | head -1)"; echo "$msg" | tail -n +2 | sed 's/^/           /'
        else item req MISSING "OpSim db" "$msg"; fi
    elif [[ -n $OPSIM_URL ]]; then item prod MISSING "OpSim db" "$OPSIM_DB missing; fetch downloads it from OPSIM_URL"
    else item req MISSING "OpSim db" "$OPSIM_DB missing and OPSIM_URL empty"; fi
    if [[ -z $EXT_TABLES ]]; then item req MISSING "EXT_TABLES" "not set: path to a prebuilt ext_tables.dat, or 'build' (needs the 7 GB DECaPS map, ~35 min)"
    elif [[ $EXT_TABLES == build ]]; then item req ok EXT_TABLES "build (maps.py)"
    elif [[ -f $P_EXT ]]; then item req ok EXT_TABLES "$EXT_TABLES"
    else item req MISSING EXT_TABLES "$EXT_TABLES not found"; fi
    local f
    for f in files/sigma_roman.txt $EEM Baseline/gbtds_layout/sca_layout_spring.txt; do
        [[ -f $f ]] && item req ok "$(basename "$f")" "in the clone" || item req MISSING "$(basename "$f")" "$f should come with the git clone"
    done
}

check_products() {
    say "Built by the pipeline:"
    local s st pop
    item prod "$(mist_state LSST Rubin)" "MIST Rubin tables" "CMD/Rubin (75 files)"
    item prod "$(mist_state Roman Roman)" "MIST Roman tables" "CMD/Roman (75 files)"
    for s in catalogue lens_ml rubin_visits roman_visits fovmap extinction sync build; do
        st=$(state_of "$s")
        item prod "$st" "$s" "$(case $s in
            catalogue) echo "CMD/components/{thin_disk,bulge,thick_disk,halo}.dat" ;;
            lens_ml) echo CMD/components/lens_ml.dat ;;
            rubin_visits) echo Baseline/BulgeBaseline.dat ;;
            roman_visits) echo Baseline/RomanBaseline.dat ;;
            fovmap) echo Baseline/lsstcam_fov/fov_map.txt ;;
            extinction) echo files/ext/ext_tables.dat ;;
            sync) echo config/data_products.h ;;
            build) echo ./roman ;; esac)"
    done
    if [[ -x roman ]]; then
        for pop in $POPULATIONS; do
            if [[ $(./roman --population "$pop" --stride 0 2>&1 </dev/null) == *"unknown population"* ]]; then
                item req MISSING "population $pop" "./roman does not know it (./roman --population $pop lists the valid names)"
            fi
        done
    fi
}

stage_check() {
    say "Pipeline check ($CONFIG_ARG; populations: $POPULATIONS; scheduler: $SCHEDULER)"
    check_tools; say; check_inputs; say; check_python; say; check_products
    say
    if (( N_REQ + N_PROD == 0 )); then say "Everything is in place."; return 0; fi
    say "$N_REQ item(s) you must provide, $N_PROD item(s) the pipeline will build or refresh."
    return 1
}

# ------------------------------------------------------------------------------------------------
# Stages: setup, fetch, prep
# ------------------------------------------------------------------------------------------------
stage_setup() {
    if [[ ! -x $PY ]]; then
        case $PY in
            */bin/python|*/bin/python3) ;;
            *) die "PYTHON=$PYTHON does not exist and is not <venv>/bin/python, so setup cannot create it" ;;
        esac
        have python3 || die "python3 not found"
        say "  creating the venv $(dirname "$(dirname "$PY")")"
        python3 -m venv "$(dirname "$(dirname "$PY")")" || die "python3 -m venv failed (Debian/Ubuntu: apt install python3-venv)"
    fi
    run_step setup "$STAMPS/setup.stamp" "$LOGS/setup.log" "python=$PY" "requirements.txt" "$PY" \
        "$PY" -m pip install -r requirements.txt || return 1
    setup_gsl
}

# GSL: the system copy, or a static build in deps/gsl that the Makefile picks up by itself.
gsl_system() { echo '#include <gsl/gsl_matrix.h>' | g++ -E -x c++ - >/dev/null 2>&1; }
gsl_local()  { [[ -f $ROOT/deps/gsl/include/gsl/gsl_matrix.h && -f $ROOT/deps/gsl/lib/libgsl.a ]]; }
build_gsl() {
    local tmp
    tmp=$(mktemp -d "${TMPDIR:-/tmp}/gsl.XXXXXX") || return 1
    curl -fL --retry 3 -sS -o "$tmp/gsl.tar.gz" "$GSL_URL" \
        && tar -xzf "$tmp/gsl.tar.gz" -C "$tmp" \
        && (cd "$tmp"/gsl-*/ && ./configure --prefix="$ROOT/deps/gsl" --disable-shared --enable-static \
            && make -j"$JOBS" && make install)
    local rc=$?
    rm -rf "$tmp"
    return $rc
}
setup_gsl() {
    case $GSL in
        system) return 0 ;;
        auto)   gsl_system && return 0 ;;
        build)  ;;
        *) die "GSL=$GSL: expected auto, system or build" ;;
    esac
    have make && have curl || die "building GSL needs make and curl"
    run_step gsl "$STAMPS/gsl.stamp" "$LOGS/gsl.log" "url=$GSL_URL" "" "$ROOT/deps/gsl/lib/libgsl.a" build_gsl || return 1
}

stage_fetch() {
    local rc=0 sys dir
    mkdir -p "$LOGS"
    say "fetch:"
    for pair in "LSST Rubin" "Roman Roman"; do
        sys=${pair% *}; dir=${pair#* }
        if [[ $FORCE != 1 && $(mist_state "$sys" "$dir") == current ]]; then say "  skip   MIST $sys tables (present)"
        elif do_mist "$sys" "$dir" >"$LOGS/mist_$sys.log" 2>&1; then say "  done   MIST $sys tables -> CMD/$dir/"
        else say "  FAILED MIST $sys tables (log: pipeline/logs/mist_$sys.log)"; rc=1; fi
    done
    run_named_step fovmap || rc=1
    if [[ -f $P_OPSIM ]]; then say "  skip   OpSim database (present: $OPSIM_DB)"
    elif [[ -n $OPSIM_URL ]]; then do_opsim && say "  done   OpSim database" || { say "  FAILED OpSim download"; rc=1; }
    else say "  OpSim database $OPSIM_DB is missing and OPSIM_URL is empty: download one yourself"; rc=1; fi
    if [[ $EXT_TABLES == build ]]; then run_named_step dustmaps || rc=1; fi
    if [[ -f $P_BESANCON ]]; then say "  ok     Besancon catalogue present"
    else say "  NOTE   Besancon catalogue $BESANCON_CATALOGUE not found. It cannot be fetched by script (web form);"
         say "         the parameters to enter are in CMD/Besancon/RERUN_bos10.md (summary: pipeline/README.md)."; fi
    return $rc
}

stage_prep() {
    need_py
    [[ -f $P_BESANCON ]] || die "Besancon catalogue not found: $BESANCON_CATALOGUE (see pipeline/README.md)"
    [[ -f $P_OPSIM ]]    || die "OpSim database not found: $OPSIM_DB (run the fetch stage, or set OPSIM_DB)"
    [[ -n $EXT_TABLES ]] || die "EXT_TABLES is not set: give the path of a prebuilt ext_tables.dat, or 'build'"
    [[ $EXT_TABLES == build || -f $P_EXT ]] || die "EXT_TABLES file not found: $EXT_TABLES"
    "$PY" pipeline/validate_inputs.py besancon "$P_BESANCON" || exit 1
    "$PY" pipeline/validate_inputs.py opsim "$P_OPSIM" || exit 1
    [[ $(mist_state LSST Rubin) == current && $(mist_state Roman Roman) == current ]] \
        || die "MIST bolometric-correction tables missing: run the fetch stage"
    say "prep:"
    local s
    for s in $PREP_STEPS; do
        run_named_step "$s" || die "prep stopped at step '$s'"
    done
    say "prep done."
}

# ------------------------------------------------------------------------------------------------
# Stages: sim (plan, chunk, schedulers)
# ------------------------------------------------------------------------------------------------
build_flags() {   # build_flags POP
    local f="--population $1 --events $EVENTS --lenses $LENSES --stride $STRIDE --stride-roman $STRIDE_ROMAN --seed $SEED"
    [[ $PAIR_SATELLITE == 1 ]] && f="$f --pair-satellite"
    [[ -n $MAXDRAWS ]] && f="$f --maxdraws $MAXDRAWS"
    [[ $STUB == 1 ]] && f="$f --stub"
    [[ -n $EXTRA_FLAGS ]] && f="$f $EXTRA_FLAGS"
    echo "$f"
}

plan_get() { sed -n "s/^$2=//p" "$1" | head -1; }

# A directory ./roman can run in: the shared inputs linked, its own outputs real.
make_rundir() {   # make_rundir DIR
    local d=$1
    mkdir -p "$d/files/MONTLMC/files" "$d/files/density" "$d/files/ext"
    ln -sfn "$ROOT/Baseline" "$d/Baseline"
    ln -sfn "$ROOT/CMD" "$d/CMD"
    ln -sfn "$ROOT/files/ext/ext_tables.dat" "$d/files/ext/ext_tables.dat"
    ln -sfn "$ROOT/files/sigma_roman.txt" "$d/files/sigma_roman.txt"
}

git_desc() { git -C "$ROOT" describe --always --dirty 2>/dev/null || echo nogit; }

write_provenance() {
    local p=$P_RUN/pipeline_provenance.txt f
    [[ -f $p ]] && return 0
    {
        echo "# pipeline provenance -- written at the first 'sim' of this RUN_ROOT"
        echo "date          $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
        echo "host          $(hostname)"
        echo "user          ${USER:-?}"
        echo "git_describe  $(git_desc)"
        echo "git_branch    $(git -C "$ROOT" rev-parse --abbrev-ref HEAD 2>/dev/null || echo ?)"
        echo "binary        $BIN  md5 $(md5f "$BIN")"
        echo "config        $CONFIG_ARG (copied below)"
        echo "besancon      $P_BESANCON  $(stat_sm "$P_BESANCON")  (size bytes, mtime; not hashed)"
        echo "opsim_db      $P_OPSIM  $(stat_sm "$P_OPSIM")  (size bytes, mtime; not hashed)"
        echo "# md5 of the data products this run reads:"
        for f in CMD/components/thin_disk.dat CMD/components/bulge.dat CMD/components/thick_disk.dat CMD/components/halo.dat \
                 CMD/components/lens_ml.dat CMD/components/provenance.txt Baseline/BulgeBaseline.dat Baseline/RomanBaseline.dat \
                 Baseline/lsstcam_fov/fov_map.txt files/ext/ext_tables.dat files/sigma_roman.txt config/data_products.h; do
            echo "md5 $(md5f "$f")  $f"
        done
        echo "# ---- config file ----"
        cat "$CONFIG"
    } > "$p"
}

ranges() { tr ' ' '\n' | awk 'NF{ if (n && $1==prev+1) { prev=$1 } else { if (n) out=out (start==prev?start:start"-"prev) ","; start=$1; prev=$1; n=1 } } END{ if (n) out=out (start==prev?start:start"-"prev); print out }'; }

chunk_dir() { printf '%s/%s/chunks/%04d\n' "$P_RUN" "$1" "$2"; }
chunk_done() { [[ -f $(chunk_dir "$1" "$2")/DONE ]]; }

stage_sim() {
    need_py
    local missing=""
    for f in roman Baseline/BulgeBaseline.dat Baseline/RomanBaseline.dat Baseline/lsstcam_fov/fov_map.txt files/ext/ext_tables.dat \
             CMD/components/lens_ml.dat CMD/components/bulge.dat files/sigma_roman.txt; do
        [[ -e $f ]] || missing="$missing $f"
    done
    [[ -z $missing ]] || die "not ready to simulate, missing:$missing -- run the prep stage"
    local s
    for s in $PREP_STEPS; do
        [[ $(state_of "$s") == current ]] || warn "prep step '$s' is not up to date; running with what is on disk (rerun prep to refresh)"
    done
    mkdir -p "$P_RUN/bin" "$RUN_STAMPS" "$P_RUN/logs"
    # The binary is COPIED into RUN_ROOT once, so a later `make` cannot change a run in progress.
    BIN=$P_RUN/bin/roman.$(git_desc)
    if [[ ! -e $BIN ]]; then cp "$BIN_SRC" "$BIN"
    elif ! cmp -s "$BIN_SRC" "$BIN"; then
        die "$BIN exists but differs from ./roman. A run must use one binary: delete RUN_ROOT/bin and the chunks to restart, or use another RUN_ROOT"
    fi
    check_populations "$BIN"
    write_provenance
    cp "$CONFIG" "$P_RUN/config_used.sh"

    local pop flags n k nchunks plist=""
    for pop in $POPULATIONS; do
        mkdir -p "$P_RUN/$pop/chunks"
        flags=$(build_flags "$pop")
        local plan=$P_RUN/$pop/plan.txt
        if [[ -f $P_RUN/$pop/MERGED && ! -d $P_RUN/$pop/chunks ]]; then say "$pop: already merged (chunks removed); nothing to do"; continue; fi
        if [[ -f $plan ]]; then
            [[ $(plan_get "$plan" flags) == "$flags" && $(plan_get "$plan" chunk_size) == "$CHUNK_SIZE" ]] \
                || die "$P_RUN/$pop was planned with other flags or CHUNK_SIZE ($(plan_get "$plan" flags), chunks of $(plan_get "$plan" chunk_size)). Use a new RUN_ROOT, or delete $P_RUN/$pop to start it over."
            [[ $(plan_get "$plan" binary) == "$BIN" ]] || die "$P_RUN/$pop was started with binary $(plan_get "$plan" binary), now $BIN. Use a new RUN_ROOT."
        else
            say "$pop: dry run to count the sightlines ..."
            make_rundir "$P_RUN/$pop/dry"
            (cd "$P_RUN/$pop/dry" && "$BIN" $flags --dry-run > dry.log 2>&1) || { tail -15 "$P_RUN/$pop/dry/dry.log" >&2; die "dry run failed for $pop (log: $P_RUN/$pop/dry/dry.log)"; }
            n=$(sed -n 's/^ *total  *\([0-9][0-9]*\) sightlines.*/\1/p' "$P_RUN/$pop/dry/dry.log" | tail -1)
            is_uint "${n:-x}" && (( n > 0 )) || die "could not read the sightline count from $P_RUN/$pop/dry/dry.log"
            rm -rf "$P_RUN/$pop/dry"
            nchunks=$(( (n + CHUNK_SIZE - 1) / CHUNK_SIZE ))
            { echo "population=$pop"; echo "n_sightlines=$n"; echo "chunk_size=$CHUNK_SIZE"; echo "n_chunks=$nchunks"
              echo "flags=$flags"; echo "binary=$BIN"; } > "$plan"
        fi
        n=$(plan_get "$plan" n_sightlines); nchunks=$(plan_get "$plan" n_chunks)
        local pend=""        # comma-separated chunk numbers still to run
        for (( k = 0; k < nchunks; k++ )); do
            if [[ $FORCE == 1 ]] || ! chunk_done "$pop" "$k"; then pend="$pend,$k"; fi
        done
        pend=${pend#,}
        say "$pop: $n sightlines in $nchunks chunks of $CHUNK_SIZE; $(echo "${pend//,/ }" | wc -w) to run"
        plist="$plist $pop:$pend"
    done

    if [[ $SCHEDULER == local ]]; then sim_local "$plist"; else sim_slurm "$plist"; fi
}

sim_local() {
    local plist=$1 jobs="" maxk=0 e pop ks k
    for e in $plist; do pop=${e%%:*}; ks=${e#*:}; for k in ${ks//,/ }; do (( k > maxk )) && maxk=$k; done; done
    # chunk-major order, so populations finish together
    for (( k = 0; k <= maxk; k++ )); do
        for e in $plist; do
            pop=${e%%:*}; ks=",${e#*:},"
            [[ $ks == *",$k,"* ]] && jobs="$jobs$pop $k"$'\n'
        done
    done
    if [[ -z $jobs ]]; then say "sim: every chunk is already done."; return 0; fi
    say "sim: $(printf '%s' "$jobs" | grep -c .) chunks, $JOBS workers (each ./roman uses ~0.7 GB), nice 10"
    printf '%s' "$jobs" | xargs -n 2 -P "$JOBS" "$SELF" "$CONFIG" chunk
    local rc=$?
    local left=0 pop2 nch
    for pop2 in $POPULATIONS; do
        [[ -f $P_RUN/$pop2/plan.txt ]] || continue
        nch=$(plan_get "$P_RUN/$pop2/plan.txt" n_chunks)
        for (( k = 0; k < nch; k++ )); do chunk_done "$pop2" "$k" || left=$((left + 1)); done
    done
    if (( left > 0 )); then say "sim: $left chunk(s) did not finish (see chunks/*/run.log); rerun sim to redo them"; return 1; fi
    say "sim: all chunks done."
    return $rc
}

sim_slurm() {
    local plist=$1 sd=$P_RUN/slurm e pop ks jid ids="" opt
    have sbatch || [[ $SUBMIT == 0 ]] || die "sbatch not found (SUBMIT=0 writes the scripts without submitting)"
    mkdir -p "$sd"
    if [[ $SUBMIT == 1 ]] && have squeue && [[ $FORCE != 1 ]]; then
        local busy; busy=$(squeue -h -u "${USER:-$(id -un)}" -o '%j' 2>/dev/null | grep -c '^rm_' || true)
        (( busy > 0 )) && die "you already have $busy queued/running rm_* jobs: rerunning sim now would duplicate their chunks (FORCE=1 to override)"
    fi
    cp "$CONFIG" "$sd/config.sh"       # the jobs read this snapshot, not the file you may edit meanwhile
    for e in $plist; do
        pop=${e%%:*}; ks=${e#*:}
        [[ -n $ks ]] || { say "slurm: $pop has nothing to run"; continue; }
        {
            echo '#!/usr/bin/env bash'
            echo "#SBATCH --job-name=rm_$pop"
            echo "#SBATCH --output=$sd/${pop}_%A_%a.out"
            echo "#SBATCH --array=$(echo "${ks//,/ }" | ranges)"
            for opt in $SLURM_OPTS; do echo "#SBATCH $opt"; done
            echo "# One array task = one chunk of $CHUNK_SIZE sightlines of population $pop."
            echo 'set -euo pipefail'
            echo "exec \"$SELF\" \"$sd/config.sh\" chunk $pop \"\$SLURM_ARRAY_TASK_ID\""
        } > "$sd/$pop.sbatch"
        bash -n "$sd/$pop.sbatch" || die "generated $sd/$pop.sbatch has a syntax error"
    done
    {
        echo '#!/usr/bin/env bash'
        echo "#SBATCH --job-name=rm_finish"
        echo "#SBATCH --output=$sd/finish_%j.out"
        for opt in $SLURM_OPTS_FINISH; do echo "#SBATCH $opt"; done
        echo "# Merge every population's chunks, then run the analyses ($ANALYSES)."
        echo 'set -euo pipefail'
        echo "\"$SELF\" \"$sd/config.sh\" merge"
        echo "\"$SELF\" \"$sd/config.sh\" analyze"
    } > "$sd/finish.sbatch"
    bash -n "$sd/finish.sbatch" || die "generated $sd/finish.sbatch has a syntax error"
    say "slurm: wrote $(ls "$sd"/*.sbatch | wc -l) scripts in $sd/"
    if [[ $SUBMIT != 1 ]]; then
        say "slurm: SUBMIT=0, nothing submitted. To submit by hand:"
        for e in $plist; do pop=${e%%:*}; [[ -n ${e#*:} ]] && say "   sbatch --parsable $sd/$pop.sbatch"; done
        say "   sbatch --dependency=afterok:<the array job ids> --kill-on-invalid-dep=yes $sd/finish.sbatch"
        return 0
    fi
    for e in $plist; do
        pop=${e%%:*}; [[ -n ${e#*:} ]] || continue
        jid=$(sbatch --parsable "$sd/$pop.sbatch") || die "sbatch failed for $pop"
        jid=${jid%%;*}; ids="$ids:$jid"
        say "slurm: $pop -> job array $jid"
    done
    if [[ -n $ids ]]; then jid=$(sbatch --parsable --dependency="afterok$ids" --kill-on-invalid-dep=yes "$sd/finish.sbatch")
    else jid=$(sbatch --parsable "$sd/finish.sbatch"); fi
    say "slurm: merge+analyze -> job ${jid%%;*} (runs when the arrays all succeed; if a chunk fails, fix it, rerun sim, and submit finish.sbatch yourself)"
}

# One chunk, in its own run directory. Done = ./roman exited 0 (DONE is written then and only then).
stage_chunk() {
    local pop=${1:-} k=${2:-}
    [[ -n $pop && $k =~ ^[0-9]+$ ]] || die "usage: pipeline.sh CONFIG chunk POP K"
    local plan=$P_RUN/$pop/plan.txt
    [[ -f $plan ]] || die "$plan not found: run the sim stage first"
    local n c nch flags bin s e d
    n=$(plan_get "$plan" n_sightlines); c=$(plan_get "$plan" chunk_size); nch=$(plan_get "$plan" n_chunks)
    flags=$(plan_get "$plan" flags); bin=$(plan_get "$plan" binary)
    (( k < nch )) || die "chunk $k out of range (0..$((nch - 1)))"
    s=$(( k * c )); e=$(( s + c )); (( e > n )) && e=$n
    d=$(chunk_dir "$pop" "$k")
    if [[ $FORCE != 1 && -f $d/DONE ]]; then say "chunk $pop $k [$s,$e): already done"; return 0; fi
    rm -rf "$d"; make_rundir "$d"
    echo "$$ $(hostname)" > "$d/RUNNING"
    # shellcheck disable=SC2064
    trap "rm -f '$d/RUNNING'" EXIT        # expanded now on purpose: $d is local and gone when the trap fires
    say "chunk $pop $k [$s,$e) start $(date '+%H:%M:%S') on $(hostname)"
    local t0=$SECONDS rc nice_cmd=""
    [[ $SCHEDULER == local ]] && have nice && nice_cmd="nice -n 10"
    # shellcheck disable=SC2086
    ( cd "$d" && $nice_cmd "$bin" $flags --start-index "$s" --end-index "$e" > run.log 2>&1 )
    rc=$?
    local tf; tf=$(mktemp "${TMPDIR:-/tmp}/times.XXXXXX"); times > "$tf"
    local cpu; cpu=$(awk 'NR==2{ for (i=1;i<=2;i++){ split($i,a,/[ms]/); t+=a[1]*60+a[2] } printf "%.0f", t }' "$tf"); rm -f "$tf"
    local entered; entered=$(grep -c 'NEW STEP' "$d/run.log" 2>/dev/null || true)
    if (( rc != 0 )); then
        say "chunk $pop $k FAILED (exit $rc); see $d/run.log"; tail -5 "$d/run.log" | sed 's/^/   | /'
        echo "failed rc=$rc" > "$d/FAILED"; return 1
    fi
    (( entered == e - s )) || warn "chunk $pop $k: log shows $entered sightlines entered, expected $((e - s))"
    { echo "start=$s"; echo "end=$e"; echo "rc=0"; echo "wall_s=$((SECONDS - t0))"; echo "cpu_s=$cpu"
      echo "sightlines_entered=$entered"; echo "host=$(hostname)"; echo "finished=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"; } > "$d/DONE"
    say "chunk $pop $k [$s,$e) done in $(fmt_dur $((SECONDS - t0))) (cpu ${cpu}s)"
}

# ------------------------------------------------------------------------------------------------
# merge, analyze, status
# ------------------------------------------------------------------------------------------------
stage_merge() {
    need_py
    mkdir -p "$RUN_STAMPS"
    local pop d rc=0 ins f nch k
    for pop in $POPULATIONS; do
        d=$P_RUN/$pop
        [[ -f $d/plan.txt ]] || { say "merge: $pop: no plan.txt (sim has not run)"; rc=1; continue; }
        if [[ -f $d/MERGED && ! -d $d/chunks ]]; then say "  skip   merge $pop (merged; chunks removed)"; continue; fi
        ins=""; nch=$(plan_get "$d/plan.txt" n_chunks)
        for (( k = 0; k < nch; k++ )); do ins="$ins $(chunk_dir "$pop" "$k")/DONE"; done
        if [[ $FORCE != 1 ]] && [[ -f $d/MERGED ]] && [[ $(step_state "$RUN_STAMPS/merge_$pop.stamp" "chunks=$nch" "$ins" "$d/MERGED") == current ]]; then
            say "  skip   merge $pop (up to date)"; continue
        fi
        # shellcheck disable=SC2086
        if ! all_exist $ins; then
            "$PY" pipeline/merge_chunks.py "$d" || rc=1      # prints which chunks lack DONE
            continue
        fi
        if run_step "merge $pop" "$RUN_STAMPS/merge_$pop.stamp" "$P_RUN/logs/merge_$pop.log" "chunks=$nch" "$ins" "$d/MERGED" \
              "$PY" pipeline/merge_chunks.py "$d"; then
            sed 's/^/         /' "$P_RUN/logs/merge_$pop.log"
            if [[ $KEEP_CHUNKS == 0 ]]; then rm -rf "$d/chunks"; say "  removed $d/chunks (KEEP_CHUNKS=0)"; fi
        else rc=1; fi
    done
    return $rc
}

prov_get() { sed -n "s/^# *$2  *\([^ #]*\).*/\1/p" "$1" | head -1; }

# Detection-only copy of a population's table: header kept, rows with detJ == 1; the column is found
# from the header, not by a fixed index.
make_detj() {   # make_detj TABLE OUT
    awk 'NR==1 { print; h=$0; sub(/^#[ \t]*/, "", h); n=split(h, c, /[ \t]+/); for (i=1;i<=n;i++) if (c[i]=="detJ") j=i; if (!j) { print "no detJ column in the header" > "/dev/stderr"; exit 2 } next }
         /^#/ { print; next } $j == 1' "$1" > "$2.tmp" && mv "$2.tmp" "$2"
}

FAILED_ANALYSES=""
# analysis NAME LOG "INPUTS" CMD...   -- never aborts the others; records failures. AN_OUTS: further files
# (besides the log) the step must have produced to count as done.
analysis() {
    local name=$1 log=$2 ins=$3; shift 3
    if ! run_step "$name" "$RUN_STAMPS/an_$(echo "$name" | tr ' /' '__').stamp" "$log" "" "$ins" "$log ${AN_OUTS:-}" "$@"; then
        FAILED_ANALYSES="$FAILED_ANALYSES
   $name   (log: $log)"
    fi
}

stage_analyze() {
    need_py
    mkdir -p "$RUN_STAMPS"
    local a pop d tag out do_y=0 do_f=0 pops=""
    for a in $ANALYSES; do [[ $a == yields ]] && do_y=1; [[ $a == figures ]] && do_f=1; done
    for pop in $POPULATIONS; do
        [[ -f $P_RUN/$pop/MERGED ]] || die "$pop has not been merged: run the merge stage first"
        pops="$pops $pop"
    done
    local A=$P_RUN/analysis
    if (( do_y )); then
        say "analyze: yields"
        local csvs=""
        for pop in $pops; do
            d=$P_RUN/$pop; out=$A/yields/$pop; mkdir -p "$out"; tag=$(prov_get "$d/files/MONTLMC/files/run_provenance.txt" population_tag)
            AN_OUTS="$out/y1_yields.csv" analysis "y1 $pop" "$out/y1.log" "$d/test$tag.dat $d/files/MONTLMC/files/MapLMC$tag.dat" \
                "$PY" analysis/memrun.py analysis/y1_absolute_yield.py --run "$pop=$d" -o "$out"
            [[ -f $out/y1_yields.csv ]] && csvs="$csvs $out/y1_yields.csv"
        done
        if [[ -n $csvs ]]; then
            local all=$A/yields/y1_yields.csv first=1 c
            : > "$all.new"
            for c in $csvs; do if (( first )); then cat "$c" >> "$all.new"; first=0; else tail -n +2 "$c" >> "$all.new"; fi; done
            if cmp -s "$all.new" "$all"; then rm -f "$all.new"; else mv "$all.new" "$all"; fi    # keep the mtime if unchanged
            mkdir -p "$A/yields/y3"
            analysis "y3" "$A/yields/y3/y3.log" "$all" "$PY" analysis/y3_yield_vs_F.py "$all" -o "$A/yields/y3"
        fi
    fi
    if (( do_f )); then
        say "analyze: figures"
        local runs="" det W P o
        for pop in $pops; do
            d=$P_RUN/$pop; tag=$(prov_get "$d/files/MONTLMC/files/run_provenance.txt" population_tag)
            [[ -n $tag ]] || { FAILED_ANALYSES="$FAILED_ANALYSES
   $pop: no population_tag in run_provenance.txt"; continue; }
            det=$d/test${tag}_detJ.dat; o=$A/figures/$pop; mkdir -p "$o"
            W="--map $d/files/MONTLMC/files/MapLMC$tag.dat --log $d/run.log"
            P="--provenance $d/files/MONTLMC/files/run_provenance.txt"
            AN_OUTS="$det" analysis "detJ extract $pop" "$o/detJ_extract.log" "$d/test$tag.dat" make_detj "$d/test$tag.dat" "$det"
            [[ -f $det ]] || continue
            local In="$det $d/files/MONTLMC/files/MapLMC$tag.dat $d/run.log"
            # shellcheck disable=SC2086
            {
            # f1 counts the non-detected draws too (its N_events column), so it reads the FULL table
            analysis "f1 $pop"      "$o/f1.log"      "$d/test$tag.dat $d/files/MONTLMC/files/MapLMC$tag.dat $d/run.log" "$PY" analysis/f1_results_table.py "$d/test$tag.dat" --fields-only -o "$o/f1.csv" $P $W
            analysis "f2_tE $pop"   "$o/f2_tE.log"   "$In" "$PY" analysis/f2_gap_filling.py "$det" --param tE -o "$o/f2_gap_filling.png" $P $W
            analysis "f2_piE $pop"  "$o/f2_piE.log"  "$In" "$PY" analysis/f2_gap_filling.py "$det" --param piE -o "$o/f2_gap_filling_piE.png" $P $W
            analysis "f3 $pop"      "$o/f3.log"      "$In" "$PY" analysis/f3_characterization_map.py "$det" -o "$o/f3_characterization_map.png" $P $W
            analysis "f4 $pop"      "$o/f4.log"      "$In" "$PY" analysis/f4_fisher_precision.py "$det" -o "$o/f4_fisher_precision.png" $P $W
            analysis "f4_all $pop"  "$o/f4_all.log"  "$In" "$PY" analysis/f4_fisher_precision.py "$det" --scope all -o "$o/f4_fisher_precision_all.png" $P $W
            analysis "h5 $pop"      "$o/h5.log"      "$In" "$PY" analysis/h5_astrometric_shift.py "$det" -o "$o/h5_astrometric_shift.png" $P $W
            if [[ $PAIR_SATELLITE == 1 ]]; then
            analysis "h3 $pop"      "$o/h3.log"      "$d/h3_pair.dat $d/files/MONTLMC/files/MapLMC$tag.dat $d/run.log" "$PY" analysis/h3_satellite_parallax.py "$d/h3_pair.dat" --out-prefix "$o/h3" $W
            fi
            }
            runs="$runs --run $pop=$d"
        done
        if [[ -n $runs ]]; then
            local xin=""; for pop in $pops; do xin="$xin $P_RUN/$pop/MERGED"; done
            mkdir -p "$A/figures"
            # shellcheck disable=SC2086
            analysis "p6" "$A/figures/p6.log" "$xin" "$PY" analysis/memrun.py analysis/p6_synergy_resolution.py $runs --detections-only -o "$A/figures/p6"
            # shellcheck disable=SC2086
            analysis "p7" "$A/figures/p7.log" "$xin" "$PY" analysis/memrun.py analysis/p7_forecast_figures.py $runs -o "$A/figures/p7"
        fi
    fi
    if [[ -n $FAILED_ANALYSES ]]; then
        say "analyze: FAILED:$FAILED_ANALYSES"
        return 1
    fi
    say "analyze: all done; results in $A/"
}

stage_status() {
    [[ -d $P_RUN ]] || { say "status: $RUN_ROOT does not exist yet (nothing simulated)"; return 0; }
    say "Run root: $P_RUN"
    [[ -f $P_RUN/pipeline_provenance.txt ]] && sed -n 's/^\(git_describe\|date\) */  \1: /p' "$P_RUN/pipeline_provenance.txt"
    local pop plan n nch c k d done_n run_n pend_n fail_n cpu cpu_all=0 ent
    for pop in $POPULATIONS; do
        plan=$P_RUN/$pop/plan.txt
        if [[ ! -f $plan ]]; then say "$pop: not planned (sim has not run)"; continue; fi
        n=$(plan_get "$plan" n_sightlines); nch=$(plan_get "$plan" n_chunks); c=$(plan_get "$plan" chunk_size)
        done_n=0; run_n=0; pend_n=0; fail_n=0; cpu=0; ent=0
        for (( k = 0; k < nch; k++ )); do
            d=$(chunk_dir "$pop" "$k")
            if [[ -f $d/DONE ]]; then done_n=$((done_n + 1)); cpu=$((cpu + $(plan_get "$d/DONE" cpu_s)))
            elif [[ -f $d/RUNNING ]]; then
                read -r pid host < "$d/RUNNING"
                if [[ $host != "$(hostname)" ]] || kill -0 "$pid" 2>/dev/null; then
                    run_n=$((run_n + 1)); ent=$((ent + $(grep -c 'NEW STEP' "$d/run.log" 2>/dev/null || true)))
                else fail_n=$((fail_n + 1)); fi        # RUNNING left behind by a killed job
            elif [[ -f $d/FAILED ]]; then fail_n=$((fail_n + 1))
            else pend_n=$((pend_n + 1)); fi
        done
        cpu_all=$((cpu_all + cpu))
        printf '%s: %d sightlines, %d chunks of %d: %d done, %d running (%d sightlines entered), %d pending, %d failed/interrupted; CPU of finished chunks %s\n' \
            "$pop" "$n" "$nch" "$c" "$done_n" "$run_n" "$ent" "$pend_n" "$fail_n" "$(fmt_dur "$cpu")"
        [[ -d $P_RUN/$pop/chunks ]] && echo "    chunks on disk: $(du -sh "$P_RUN/$pop/chunks" | cut -f1)"
        if [[ -f $P_RUN/$pop/MERGED ]]; then
            echo "    merged: $(sed -n 's/^file \([^ ]*\) lines=\([0-9]*\) bytes=\([0-9]*\)/\1 (\2 lines, \3 B)/p' "$P_RUN/$pop/MERGED" | grep -E '^test|^h3' | tr '\n' ';')"
        fi
    done
    say "Total CPU of finished chunks: $(fmt_dur "$cpu_all")"
    if [[ -d $P_RUN/analysis ]]; then say "analysis outputs: $(find "$P_RUN/analysis" -type f \( -name '*.png' -o -name '*.csv' -o -name '*.md' \) | wc -l) files under $P_RUN/analysis"; fi
}

# ------------------------------------------------------------------------------------------------
# all
# ------------------------------------------------------------------------------------------------
stage_all() {
    N_REQ=0; N_PROD=0
    # Hard prerequisites only: what the pipeline cannot supply for itself.
    ( check_tools; check_inputs ) | tee "$ROOT/pipeline/.check.tmp"
    if grep -qE '^  MISSING' "$ROOT/pipeline/.check.tmp"; then rm -f "$ROOT/pipeline/.check.tmp"; die "prerequisites missing (above); fix them, then rerun (details: pipeline.sh $CONFIG_ARG check)"; fi
    rm -f "$ROOT/pipeline/.check.tmp"
    stage_setup  || die "setup failed"
    stage_fetch  || die "fetch failed"
    stage_prep
    stage_sim    || die "sim failed"
    [[ $SCHEDULER == slurm ]] && { say "all: jobs submitted (slurm); the finish job merges and analyses. Watch with: pipeline/pipeline.sh $CONFIG_ARG status"; return 0; }
    stage_merge  || die "merge failed"
    stage_analyze
}

# One brace group, so bash has parsed the whole dispatch (and `exit`) before running any stage: editing this
# file while a stage is running cannot corrupt the running script.
{
    validate_config
    case $STAGE in
        check)   stage_check ;;
        setup)   stage_setup ;;
        fetch)   stage_fetch ;;
        prep)    stage_prep ;;
        sim)     stage_sim ;;
        chunk)   stage_chunk "$@" ;;
        merge)   stage_merge ;;
        analyze) stage_analyze ;;
        all)     stage_all ;;
        status)  stage_status ;;
        *)       usage ;;
    esac
    exit $?
}
