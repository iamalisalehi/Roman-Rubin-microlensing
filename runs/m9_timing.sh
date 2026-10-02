#!/usr/bin/env bash
# M9 (2026-10-02): per-sightline CPU/wall timing on a stratified sample, production flags, for
# bulge / bh / ns in parallel. Sightlines are reproducible individually (Deviation 77), so each is
# run alone with --start-index i --end-index i+1. Output: runs/m9_timing_<pop>/times.csv
set -u
R=$(cd "$(dirname "$0")/.." && pwd)
IDX_F="$1"; IDX_O="$2"
run_pop() {
  pop=$1; D=$R/runs/m9_timing_$pop
  mkdir -p $D/files/MONTLMC/files && cd $D
  ln -sfn $R/Baseline Baseline; ln -sfn $R/CMD CMD
  for f in density ext sigmaA_LSST.txt sigma_roman.txt; do ln -sfn $R/files/$f files/$f; done
  cp $R/roman ./roman
  echo "pop,stratum,index,wall_s,user_s,sys_s" > times.csv
  for s in F O; do
    if [ $s = F ]; then idx="$IDX_F"; else idx="$IDX_O"; fi
    for i in $idx; do
      t0=$(date +%s%N)
      /usr/bin/env bash -c "TIMEFORMAT='%U %S'; time ./roman --population $pop --events 300 --lenses 50 --stride-roman 5 --pair-satellite --start-index $i --end-index $((i+1)) > log_$i.txt 2>&1" 2> t_$i.txt
      t1=$(date +%s%N)
      read u sy < <(tail -1 t_$i.txt)
      echo "$pop,$s,$i,$(( (t1 - t0) / 1000000000 )),$u,$sy" >> times.csv
    done
  done
  echo DONE >> times.csv
}
for pop in bulge bh ns; do run_pop $pop & done
wait
echo ALL DONE
