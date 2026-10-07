#!/bin/bash
# E5: projection floors at the rezone times (operator independent; see RunFloors in experiments/exp5_ale.cpp).
# Appends to results/exp5_floors.csv.  2 jobs at a time.
set -u
cd "$(dirname "$0")/.."
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
OUT=results/exp5_floors_parts
mkdir -p $OUT results/exp5_logs
JOBS=$(mktemp)
for nrez in 20 5 2; do
  for p in 3 2; do
    for eps in 0.3 0.6; do
      echo "fl_p${p}_q2_e${eps}_r${nrez} -N 8 -p $p -q 2 -eps $eps -K 100 -rezones $nrez" >> $JOBS
    done
  done
done
for nrez in 5 2; do echo "fl_p3_q3_e0.6_r${nrez} -N 8 -p 3 -q 3 -eps 0.6 -K 100 -rezones $nrez" >> $JOBS; done
run() { tag=$1; shift; [ -s results/exp5_floors_parts/$tag.csv ] && return
  mpirun -np 1 --oversubscribe bin/exp5_ale -floors "$@" -csv results/exp5_floors_parts/$tag.csv > results/exp5_logs/$tag.log 2>&1; }
export -f run
cat $JOBS | xargs -P 2 -L 1 bash -c 'run "$@"' _
rm -f $JOBS
python3 scripts/merge_csv.py results/exp5_floors.csv results/exp5_floors_parts/*.csv
