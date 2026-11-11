#!/usr/bin/env bash
# Smoke test of the pipeline's compressed-table path (COMPRESS_TABLES=1 -> .colxz per chunk, merged by
# concatenation), against COMPRESS_TABLES=0 (plain) on the stub patch, same seed.
#
#   tests/pipeline_colxz_smoke.sh                 # modes 1 and 0
#   WITH_GZIP=1 tests/pipeline_colxz_smoke.sh     # also gzip, compared with 0
#
# Environment: PYTHON (interpreter with numpy/pandas; default .roman/bin/python of this checkout),
# SMOKE_ROOT (scratch directory the test creates and may delete; default a fresh mktemp -d),
# KEEP=1 keeps it. Exit 0 = every check passed.
set -u
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
cd "$REPO" || exit 2
PY=${PYTHON:-$REPO/.roman/bin/python}
[[ -x $PY ]] || { echo "no python at $PY (set PYTHON)"; exit 2; }
SR=${SMOKE_ROOT:-$(mktemp -d "${TMPDIR:-/tmp}/colxz_smoke.XXXXXX")}
mkdir -p "$SR"
SR=$(cd "$SR" && pwd -P)
case $SR in /|"$HOME"|"$REPO") echo "refusing SMOKE_ROOT=$SR"; exit 2 ;; esac
PIPE="bash $REPO/pipeline/pipeline.sh"
CT="$PY $REPO/analysis/coltable.py"
FAILS=0; PASSES=0
ok()   { PASSES=$((PASSES + 1)); echo "ok   $*"; }
fail() { FAILS=$((FAILS + 1)); echo "FAIL $*"; }
check() { local msg=$1; shift; if "$@" >/dev/null 2>&1; then ok "$msg"; else fail "$msg"; fi; }

mk_config() {   # mk_config FILE ROOT MODE
    cp "$REPO/pipeline/config.example.sh" "$1"
    cat >> "$1" <<CFG
STUB=1; STRIDE=1; STRIDE_ROMAN=1; EVENTS=10; LENSES=2; MAXDRAWS=3000; CHUNK_SIZE=10
POPULATIONS=bulge
RUN_ROOT=$2
COMPRESS_TABLES=$3
PYTHON=$PY
ALLOW_DIRTY=1
OGLE_TABLE3=
JOBS=2
SEED=42
${EXTRA_CFG:-}
CFG
}

peak_sampler() {   # peak_sampler ROOT OUTFILE : max of du -sb every 5 s until killed
    local m=0 s
    while :; do
        s=$(du -sb "$1" 2>/dev/null | cut -f1); [[ -n $s ]] && (( s > m )) && { m=$s; echo "$m" > "$2"; }
        sleep 5
    done
}

run_mode() {   # run_mode NAME MODE : sim + merge in $SR/NAME, prints the disk peak
    local name=$1 mode=$2 root=$SR/$1 cfg=$SR/$1.sh t0=$SECONDS
    mk_config "$cfg" "$root" "$mode"
    echo 0 > "$SR/$name.peak"
    peak_sampler "$root" "$SR/$name.peak" & local sp=$!
    $PIPE "$cfg" sim   > "$SR/$name.sim.out" 2>&1; local rs=$?
    $PIPE "$cfg" merge > "$SR/$name.merge.out" 2>&1; local rm=$?
    kill $sp 2>/dev/null; wait $sp 2>/dev/null
    local pk; pk=$(du -sb "$root" | cut -f1); [[ $(cat "$SR/$name.peak") -gt $pk ]] && pk=$(cat "$SR/$name.peak")
    echo "$pk" > "$SR/$name.peak"
    echo "== mode $mode: sim exit $rs, merge exit $rm, $((SECONDS - t0)) s, disk peak (du -sb every 5 s) $pk B, final $(du -sb "$root" | cut -f1) B"
    [[ $rs == 0 && $rm == 0 ]] || { tail -20 "$SR/$name.sim.out" "$SR/$name.merge.out"; fail "mode $mode: sim/merge exit codes $rs/$rm"; return 1; }
    ok "mode $mode: sim and merge exit 0"
}

plain_tables()  { find "$1" -name 'test*.dat' ! -name '*_detJ.dat' -type f; }

# ---- 1. empty chunks: the merge of .colxz chunks decodes to what the plain merge writes ------------------
echo "== empty-chunk merge (merge_chunks.merge_file on plain and .colxz chunks)"
E=$SR/empty; rm -rf "${E:?}"; mkdir -p "$E"
if "$PY" - "$REPO" "$E" <<'PY'
import os, sys
repo, E = sys.argv[1:3]
sys.path.insert(0, os.path.join(repo, "pipeline")); sys.path.insert(0, os.path.join(repo, "analysis"))
import merge_chunks, coltable
H = "# a b c detJ\n"
rows = lambda s, n: "".join(f"{s+i} {0.5*(s+i)} {s-i} {(s+i)%2}\n" for i in range(n))
cases = {
  "normal":       [H + rows(0, 7), H + rows(10, 5), H + rows(20, 6)],
  "empty_first":  ["", H + rows(10, 5), H + rows(20, 6)],
  "empty_mid":    [H + rows(0, 7), "", H + rows(20, 6)],
  "hdronly_first":[H, H + rows(10, 5), H + rows(20, 6)],
  "hdronly_mid":  [H + rows(0, 7), H, H + rows(20, 6)],
  "all_hdronly":  [H, H, H],
  "all_empty":    ["", "", ""],
}
bad = 0
for name, chunks in cases.items():
    for kind in ("plain", "colxz"):
        d = os.path.join(E, name + "_" + kind)
        cds = []
        for k, txt in enumerate(chunks):
            cd = os.path.join(d, f"c{k}"); os.makedirs(cd); cds.append(cd)
            p = os.path.join(cd, "test5.dat")
            open(p, "w").write(txt)
            if kind == "colxz":
                coltable.encode(p, p + ".colxz"); os.remove(p)
        rel = "test5.dat" + (".colxz" if kind == "colxz" else "")
        merge_chunks.merge_file(rel, cds, d)
    plain_p = os.path.join(E, name + "_plain", "test5.dat")
    plain = open(plain_p, "rb").read()
    cx = os.path.join(E, name + "_colxz", "test5.dat.colxz")
    out = os.path.join(E, name + ".dec")
    with open(out, "wb") as fh:
        coltable.decode(cx, fh)
    dec = open(out, "rb").read()
    nrows = sum(b.nrows for b in coltable.iter_blocks(cx))
    want_rows = max(plain.count(b"\n") - 1, 0)
    good = dec == plain and nrows == want_rows
    print(("ok   " if good else "FAIL ") + f"empty-chunk case {name}: decoded {len(dec)} B vs plain merge {len(plain)} B, rows {nrows} (want {want_rows})")
    bad += not good

# refusals: both forms must refuse exactly as the plain merge does
def refuses(name, mk, kinds):
    global bad
    for kind in kinds:
        d = os.path.join(E, name + "_" + kind)
        cds = []
        for k in range(2):
            cd = os.path.join(d, f"c{k}"); os.makedirs(cd); cds.append(cd)
        mk(kind, cds)
        rel = "test5.dat" + (".colxz" if kind == "colxz" else "")
        try:
            merge_chunks.merge_file(rel, cds, d)
            good, msg = False, "merge accepted it"
        except SystemExit as e:
            good, msg = True, str(e)[:90]
        print(("ok   " if good else "FAIL ") + f"refusal {name} ({kind}): {msg}")
        bad += not good

def write(kind, cd, txt):
    p = os.path.join(cd, "test5.dat")
    open(p, "w").write(txt)
    if kind == "colxz":
        coltable.encode(p, p + ".colxz"); os.remove(p)

def cut_last(kind, cds):
    write(kind, cds[0], H + "1 2 3\n4 5")
    write(kind, cds[1], H + rows(10, 3))
refuses("cut_last_row", cut_last, ("plain", "colxz"))

def zero_byte(kind, cds):
    write(kind, cds[0], H + rows(0, 3)); write(kind, cds[1], H + rows(10, 3))
    open(os.path.join(cds[1], "test5.dat.colxz"), "wb").close()
refuses("zero_byte_colxz", zero_byte, ("colxz",))

def truncated(kind, cds):
    write(kind, cds[0], H + rows(0, 3)); write(kind, cds[1], H + rows(10, 300))
    q = os.path.join(cds[1], "test5.dat.colxz")
    b = open(q, "rb").read(); open(q, "wb").write(b[:len(b) // 2])
refuses("truncated_colxz", truncated, ("colxz",))
sys.exit(1 if bad else 0)
PY
then ok "empty-chunk cases"; else fail "empty-chunk cases"; fi

# ---- 2. the runs -------------------------------------------------------------------------------------------
run_mode m1 1
run_mode m0 0
if [[ ${WITH_GZIP:-0} == 1 ]]; then run_mode mg gzip; fi

T1=$(find "$SR/m1/bulge" -maxdepth 1 -name 'test*.dat.colxz' | head -1)
T0=$(find "$SR/m0/bulge" -maxdepth 1 -name 'test*.dat' ! -name '*_detJ.dat' | head -1)
[[ -n $T1 ]] && ok "mode 1: merged table is ${T1##*/}" || fail "mode 1: no merged test*.dat.colxz in $SR/m1/bulge"
[[ -n $T0 ]] && ok "mode 0: merged table is ${T0##*/}" || fail "mode 0: no merged plain test*.dat"

check "mode 1: no plain test*.dat in any chunk or the merged directory" test -z "$(plain_tables "$SR/m1/bulge")"
n_chunks=$(find "$SR/m1/bulge/chunks" -maxdepth 2 -name 'test*.dat.colxz' | wc -l)
check "mode 1: every chunk holds a .colxz ($n_chunks of 4)" test "$n_chunks" -eq 4
check "mode 1: no .tmp left" test -z "$(find "$SR/m1" -name '*.tmp')"
if [[ -n $T1 && -n $T0 ]]; then
    $CT decode "$T1" > "$SR/m1.dec"
    check "mode 1: merged .colxz decodes byte-identical (cmp) to the mode-0 merged table" cmp "$SR/m1.dec" "$T0"
    rows=$($CT rows "$T1"); lines=$(( $(wc -l < "$T0") - 1 ))
    check "mode 1: row count $rows equals the plain table's $lines" test "$rows" -eq "$lines"
    check "mode 1: header decoded once" test "$(grep -c '^#' "$SR/m1.dec")" -eq "$(grep -c '^#' "$T0")"
    for c in "$SR"/m1/bulge/chunks/*/; do
        k=${c%/}; k=${k##*/}
        check "mode 1: chunk $k decodes byte-identical to mode 0's chunk" cmp <($CT decode "$c"test*.dat.colxz) "$SR"/m0/bulge/chunks/$k/test*.dat
    done
fi
check "mode 1: other merged files equal mode 0's (h3_pair, MapLMC, LpLMC)" \
    bash -c "cd $SR/m0/bulge && for f in h3_pair.dat files/MONTLMC/files/MapLMC*.dat files/MONTLMC/files/LpLMC*.dat; do cmp \$f $SR/m1/bulge/\$f || exit 1; done"

if [[ ${WITH_GZIP:-0} == 1 ]]; then
    TG=$(find "$SR/mg/bulge" -maxdepth 1 -name 'test*.dat.gz' | head -1)
    [[ -n $TG ]] && ok "gzip: merged table is ${TG##*/}" || fail "gzip: no merged .gz"
    check "gzip: merged .gz decodes (gzip -dc) identical to plain" cmp <(gzip -dc "$TG") "$T0"
    check "gzip: chunks stay plain" test -n "$(plain_tables "$SR/mg/bulge/chunks")"
fi

# ---- 3. analyze: exit 0 on mode 1, detJ subsets identical to mode 0 -----------------------------------------
$PIPE "$SR/m1.sh" analyze > "$SR/m1.analyze.out" 2>&1; ra1=$?
$PIPE "$SR/m0.sh" analyze > "$SR/m0.analyze.out" 2>&1; ra0=$?
# On the stub patch the cross-population steps (y3, p6) fail for want of data, in plain mode too: what is
# asserted is that mode 1 behaves exactly like mode 0 (same exit code, same failed steps) and that every
# per-population step reads the .colxz and writes the same numbers.
failed_steps() { sed -n '/^analyze: FAILED:/,$p' "$1" | sed 's/(log:.*//' | sort; }
check "mode 1: analyze exit code equals mode 0's ($ra1 vs $ra0)" test "$ra1" -eq "$ra0"
check "mode 1: same failed analysis steps as mode 0" test "$(failed_steps "$SR/m1.analyze.out")" = "$(failed_steps "$SR/m0.analyze.out")"
check "mode 1: no per-population analysis step failed" bash -c "! grep -E '^   (y1|f[0-9]|h[35]|detJ|w2)' <(sed -n '/^analyze: FAILED:/,\$p' $SR/m1.analyze.out)"
check "mode 1: y1_yields.csv identical to mode 0's" cmp "$SR/m1/analysis/yields/y1_yields.csv" "$SR/m0/analysis/yields/y1_yields.csv"
check "mode 1: the f1 results table is identical to mode 0's" cmp "$SR/m1/analysis/figures/bulge/f1.csv" "$SR/m0/analysis/figures/bulge/f1.csv"
D1=$(ls "$SR"/m1/bulge/test*_detJ.dat 2>/dev/null | head -1); D0=$(ls "$SR"/m0/bulge/test*_detJ.dat 2>/dev/null | head -1)
if [[ -n $D1 && -n $D0 ]]; then check "detJ subsets identical (cmp), mode 1 vs 0" cmp "$D1" "$D0"; else fail "detJ subset missing (mode 1: '$D1', mode 0: '$D0')"; fi
check "detJ subset is a plain file in mode 1" test "$(head -c1 "${D1:-/dev/null}")" = "#"
if [[ ${WITH_GZIP:-0} == 1 ]]; then
    $PIPE "$SR/mg.sh" analyze > "$SR/mg.analyze.out" 2>&1
    DG=$(ls "$SR"/mg/bulge/test*_detJ.dat 2>/dev/null | head -1)
    check "gzip: detJ subset identical to mode 0's" cmp "$DG" "$D0"
fi

# ---- 4. a header mismatch between chunks makes the merge refuse ------------------------------------------------
echo "== header mismatch"
cp "$SR/m1.sh" "$SR/hm.sh"; sed -i "s|^RUN_ROOT=.*|RUN_ROOT=$SR/hm|" "$SR/hm.sh"
rm -rf "${SR:?}/hm"; cp -a "$SR/m1" "$SR/hm"
c1=$(ls "$SR"/hm/bulge/chunks/0001/test*.dat.colxz)
$CT decode "$c1" | sed '1s/$/ extracol/' > "$SR/hm.plain"
$CT encode "$SR/hm.plain" "$c1" >/dev/null
rm -f "${SR:?}/hm/bulge/MERGED" "$SR"/hm/bulge/test*.dat.colxz
FORCE=1 $PIPE "$SR/hm.sh" merge > "$SR/hm.merge.out" 2>&1; rhm=$?
check "merge refuses (exit != 0) when a chunk's header differs (rc=$rhm)" test "$rhm" -ne 0
check "the refusal names the header" grep -qi 'header\|leading' "$SR/hm.merge.out"
check "no merged table written after the refusal" test -z "$(ls "$SR"/hm/bulge/test*.dat.colxz 2>/dev/null)"
rm -rf "${SR:?}/hm"

# ---- 5. an interrupted chunk is re-run, the merge is identical --------------------------------------------------
echo "== interrupted chunk"
d2=$SR/m1/bulge/chunks/0002
$CT decode "$d2"/test*.dat.colxz > "$SR/plain2.tmp"
t2=$(basename "$(ls "$d2"/test*.dat.colxz)" .colxz)
rm -f "${d2:?}/DONE" "${d2:?}/$t2.colxz"; cp "$SR/plain2.tmp" "$d2/$t2"; echo junk > "$d2/$t2.colxz.tmp"
$PIPE "$SR/m1.sh" sim > "$SR/m1.resim.out" 2>&1; rr=$?
check "re-run of sim after the interruption exits 0 (rc=$rr)" test "$rr" -eq 0
check "chunk 2 has DONE again, a .colxz, no plain table and no .tmp" \
    bash -c "test -f $d2/DONE && test -f $d2/$t2.colxz && test ! -e $d2/$t2 && test ! -e $d2/$t2.colxz.tmp"
$PIPE "$SR/m1.sh" merge > "$SR/m1.remerge.out" 2>&1; rm_=$?
check "re-merge exits 0 (rc=$rm_)" test "$rm_" -eq 0
check "merged .colxz still decodes identical to mode 0" cmp <($CT decode "$T1") "$T0"
# a chunk that has DONE but holds both forms is refused, naming the chunk
cp "$SR/plain2.tmp" "$d2/$t2"
FORCE=1 $PIPE "$SR/m1.sh" merge > "$SR/m1.both.out" 2>&1; rb=$?
check "merge refuses a DONE chunk holding both a plain table and a .colxz (rc=$rb)" test "$rb" -ne 0
check "  ... and names chunk 0002" grep -q '0002' "$SR/m1.both.out"
rm -f "${d2:?}/$t2"
check "merged .colxz untouched by the refused merge" cmp <($CT decode "$T1") "$T0"

# ---- 6. KEEP_CHUNKS=0: a second sim after the merge starts no chunk ------------------------------------------------
echo "== sim after a KEEP_CHUNKS=0 merge"
EXTRA_CFG="KEEP_CHUNKS=0" run_mode mk 1
TK=$(find "$SR/mk/bulge" -maxdepth 1 -name 'test*.dat.colxz' | head -1)
check "KEEP_CHUNKS=0: chunks directory removed after the merge" test ! -d "$SR/mk/bulge/chunks"
cp "$TK" "$SR/mk.before"
$PIPE "$SR/mk.sh" sim > "$SR/mk.resim.out" 2>&1; rk=$?
check "second sim exits 0 (rc=$rk)" test "$rk" -eq 0
check "second sim prints a skip line" grep -q 'skip   sim bulge' "$SR/mk.resim.out"
check "second sim started no chunk and created no chunks/ directory" bash -c "test ! -d $SR/mk/bulge/chunks && ! grep -q 'start' $SR/mk.resim.out"
check "merged table untouched by the second sim" cmp "$TK" "$SR/mk.before"
check "KEEP_CHUNKS=0 merged table still decodes identical to mode 0" cmp <($CT decode "$TK") "$T0"

# ---- report -------------------------------------------------------------------------------------------------------
echo "== disk peak (du -sb RUN_ROOT, sampled every 5 s): mode 1 $(cat "$SR/m1.peak") B, mode 0 $(cat "$SR/m0.peak") B${WITH_GZIP:+, gzip $(cat "$SR/mg.peak" 2>/dev/null) B}"
for k in "$SR"/m1/bulge/chunks/*/DONE; do echo "   $(basename "$(dirname "$k")") DONE: $(tr '\n' ' ' < "$k")"; done
echo "== $PASSES passed, $FAILS failed (scratch: $SR)"
[[ ${KEEP:-0} == 1 || $FAILS -gt 0 ]] || rm -rf "${SR:?}"
exit $((FAILS > 0))
