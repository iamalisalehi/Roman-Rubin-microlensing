#!/usr/bin/env bash
# runctl.sh -- start, pause, resume and continue a production run of ./roman.
#
# WHY THIS EXISTS. A production scan is hours of CPU on a laptop that sleeps, throttles and
# gets closed. Three facts about this program decide how it can be interrupted:
#
#   1. Every output is APPEND-mode on a resume (--start-index > 0) and TRUNCATING otherwise.
#      Resuming without --start-index wipes the run. That has happened once already and cost
#      1,936,653 rows -- DEVIATIONS.md 34.
#   2. The per-sightline streams are flushed after every sightline (DEVIATIONS.md 42), so a
#      killed run keeps everything up to the last COMPLETED sightline. Before that fix a kill
#      lost the buffered tail silently.
#   3. The resume index is the number of sightlines the scan ENTERED -- `grep -c 'NEW STEP'`
#      on the log -- NOT the line count of MapLMC<tag>.dat, which undercounts because a
#      no-coverage or barren sightline is entered but writes no map row.
#
# So there are two ways to interrupt, and they are not equivalent:
#
#   pause/resume (SIGSTOP/SIGCONT) -- costs nothing, loses nothing, keeps the process and its
#       RNG stream exactly where they were. This is the one to use for "I need my cores back
#       for an hour". A stopped process still holds its memory; it does not survive a reboot.
#
#   stop/continue (SIGTERM, then a fresh process with --start-index) -- survives a reboot, but
#       the LAST sightline was mid-flight when the signal landed and its rows are already on
#       disk. Continuing past it leaves a half-sampled sightline whose event count no longer
#       represents the sky area it stands for, which quietly biases every w_area-weighted
#       yield. So `continue` TRUNCATES that sightline's rows first and redoes it whole. That
#       is what was done by hand for the v3 run's chunk boundary; here it is a script.
#
# USAGE
#   runctl.sh start  <dir>         launch (flags come from <dir>/run.flags)
#   runctl.sh status <dir>         progress, rate, ETA, output sizes
#   runctl.sh pause  <dir>         SIGSTOP
#   runctl.sh resume <dir>         SIGCONT
#   runctl.sh stop   <dir>         SIGTERM, leaving the run continuable
#   runctl.sh continue <dir> [--dry-run]   truncate the partial sightline and relaunch
#
# <dir> is a run directory: Baseline/, CMD/ and files/{density,ext,sigma*} symlinked to the
# repo, files/MONTLMC/files/ a REAL directory, and the binary copied in so a later rebuild
# cannot change what is running.

set -euo pipefail

die() { echo "runctl: $*" >&2; exit 1; }

cmd=${1:-}; dir=${2:-}
[[ -n $cmd && -n $dir ]] || die "usage: runctl.sh {start|status|pause|resume|stop|continue} <dir>"
dir=$(cd "$dir" && pwd) || die "no such directory: $dir"

pidfile=$dir/run.pid
logfile=$dir/run.log
flagfile=$dir/run.flags
idxfile=$dir/run.startindex     # the --start-index the CURRENT process was launched with

# The binary is whatever roman.* was copied into the run directory. Copied, not symlinked:
# a run must not change under a `make`.
binary=$(ls "$dir"/roman.* 2>/dev/null | head -1) || true

running_pid() {
    [[ -f $pidfile ]] || return 1
    local p; p=$(cat "$pidfile")
    [[ -n $p ]] && kill -0 "$p" 2>/dev/null && { echo "$p"; return 0; }
    return 1
}

proc_state() {  # one letter from ps: R/S running or sleeping, T stopped, Z zombie
    ps -o state= -p "$1" 2>/dev/null | tr -d ' ' | cut -c1
}

entered() {     # sightlines this run has ENTERED, across all chunks
    local base=0
    [[ -f $idxfile ]] && base=$(cat "$idxfile")
    local here=0
    [[ -f $logfile ]] && here=$(grep -c 'NEW STEP' "$logfile" || true)
    echo $(( base + here ))
}

case $cmd in

start)
    running_pid >/dev/null && die "already running (pid $(cat "$pidfile"))"
    [[ -n ${binary:-} && -x $binary ]] || die "no executable roman.* in $dir"
    [[ -f $flagfile ]] || die "no $flagfile -- write the run's flags there, one line"
    # A fresh start truncates the outputs by design. Refuse if a previous run left data,
    # rather than silently destroying it; `continue` is the way to extend a run.
    if [[ -s $dir/run.log ]]; then
        die "$logfile is not empty -- this directory holds a run already. Use 'continue', or move it aside."
    fi
    echo 0 > "$idxfile"
    flags=$(cat "$flagfile")
    cd "$dir"
    # setsid so the run survives the terminal that launched it; nice so it yields to the
    # user's interactive work rather than making the laptop unusable for hours.
    setsid nice -n 10 "$binary" $flags >> "$logfile" 2>&1 &
    echo $! > "$pidfile"
    echo "runctl: started pid $(cat "$pidfile") in $dir"
    echo "runctl: flags: $flags"
    ;;

status)
    n=$(entered)
    if p=$(running_pid); then
        st=$(proc_state "$p")
        case $st in
            T) state="PAUSED (SIGSTOP)" ;;
            Z) state="ZOMBIE" ;;
            *) state="running" ;;
        esac
        # Elapsed CPU, not wall clock: a paused or throttled run makes wall clock a lie.
        cpu=$(ps -o cputime= -p "$p" | tr -d ' ')
        rss=$(ps -o rss= -p "$p" | tr -d ' ')
        echo "state:     $state (pid $p, cpu $cpu, rss $(( rss / 1024 )) MB)"
    else
        echo "state:     not running"
    fi
    echo "entered:   $n sightlines"
    if [[ -f $logfile ]]; then
        # Rate over the whole log is the honest number to quote: per-sightline cost varies
        # 8-18x between footprint and outside strata, so a recent-window rate extrapolates
        # badly in both directions (PROGRESS.md 5d).
        start_epoch=$(stat -c %Y "$logfile")
        now=$(date +%s)
        echo "log:       $(du -h "$logfile" | cut -f1), last modified $(( (now - start_epoch) / 60 )) min ago"
        tail -3 "$logfile" | sed 's/^/           /'
    fi
    # A run still in its ~170 s startup has produced no .dat files at all, so the globs match
    # nothing. Without the guard the last [[ -f ]] is the command that sets the exit status and
    # `status` reports failure for a perfectly healthy run -- which would teach the reader to
    # ignore its exit code.
    shopt -s nullglob
    for f in "$dir"/test*.dat "$dir"/files/MONTLMC/files/*.dat; do
        printf "output:    %-28s %s\n" "$(basename "$f")" "$(du -h "$f" | cut -f1)"
    done
    shopt -u nullglob
    exit 0
    ;;

pause)
    p=$(running_pid) || die "not running"
    kill -STOP "$p"
    echo "runctl: paused pid $p -- nothing is lost, the process keeps its state and memory."
    echo "runctl: it will NOT survive a reboot. Use 'stop' for that."
    ;;

resume)
    p=$(running_pid) || die "not running"
    [[ $(proc_state "$p") == T ]] || die "pid $p is not stopped"
    kill -CONT "$p"
    echo "runctl: resumed pid $p"
    ;;

stop)
    p=$(running_pid) || die "not running"
    # A stopped process cannot act on SIGTERM; wake it first or it lingers.
    [[ $(proc_state "$p") == T ]] && kill -CONT "$p"
    kill -TERM "$p"
    sleep 2
    kill -0 "$p" 2>/dev/null && { sleep 5; kill -KILL "$p" 2>/dev/null || true; }
    echo "runctl: stopped. Entered $(entered) sightlines."
    echo "runctl: the last one was mid-flight -- 'continue' will truncate and redo it."
    ;;

continue)
    running_pid >/dev/null && die "still running -- stop it first"
    [[ -n ${binary:-} && -x $binary ]] || die "no executable roman.* in $dir"
    n=$(entered)
    [[ $n -gt 0 ]] || die "nothing entered yet -- use 'start'"
    resume_at=$(( n - 1 ))   # redo the sightline that was interrupted, whole
    here=$(dirname "$(readlink -f "$0")")
    args=(--dir "$dir" --resume-at "$resume_at")
    [[ ${3:-} == --dry-run ]] && args+=(--dry-run)
    python3 "$here/resume_truncate.py" "${args[@]}"
    [[ ${3:-} == --dry-run ]] && exit 0
    echo "$resume_at" > "$idxfile"
    : > "$logfile.$(date +%Y%m%d_%H%M%S).chunk"   # marker only; the log keeps accumulating
    flags=$(cat "$flagfile")
    cd "$dir"
    setsid nice -n 10 "$binary" $flags --start-index "$resume_at" >> "$logfile" 2>&1 &
    echo $! > "$pidfile"
    echo "runctl: continued at --start-index $resume_at, pid $(cat "$pidfile")"
    ;;

*)
    die "unknown command: $cmd"
    ;;
esac
