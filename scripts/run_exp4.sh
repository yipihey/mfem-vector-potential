#!/bin/bash
# E4: repeated remaps, ping-pong A (uniform) <-> B (f1 eps 0.3), N = 8, 1 MPI rank per run.
#   s1: field abc, no B0 (helicity), p=1,2,3, n=100, all 7 ops, diagnostics every 5 steps
#   s2: field abc WITH B0, p=2, n=100, all ops
#   s3: long run n=400, p=2, ops A_int,A_pt,B_int,B_l2c, diagnostics every 20 steps (no B0)
#   s3b: long run n=400, p=2, ops A_l2,B_l2 (added after s3 showed gauge growth of A_int), every 20 steps
# Usage: scripts/run_exp4.sh [parts...]  parts: s1_p1 s1_p2 s1_p3 s2 s3 s3b   (default: all, serial)
# Each part appends to results/exp4_parts/<part>.csv; scripts/merge_csv.py merges into results/exp4_repeat.csv
set -u
cd "$(dirname "$0")/.."
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
BIN=${BIN:-bin/exp4_repeat}
NP=${NP:-1}
mkdir -p results/exp4_parts results/exp4_logs
PARTS=${*:-s1_p1 s1_p2 s1_p3 s2 s3}

run() {  # part, args...
  local part=$1; shift
  local log=results/exp4_logs/$part.log
  echo "== [$part] $*" | tee -a "$log"
  local t0=$(date +%s)
  mpirun -np $NP --oversubscribe $BIN "$@" -tag "$part" -csv results/exp4_parts/$part.csv >> "$log" 2>&1
  local rc=$?
  echo "   rc=$rc  $(( $(date +%s) - t0 )) s" | tee -a "$log"
}

for part in $PARTS; do
  case $part in
    s1_p1) run $part -N 8 -p 1 -field abc -no-b0 -epsB 0.3 -n ${N1:-100} -every 5 -ops all ;;
    s1_p2) run $part -N 8 -p 2 -field abc -no-b0 -epsB 0.3 -n ${N1:-100} -every 5 -ops all ;;
    s1_p3) run $part -N 8 -p 3 -field abc -no-b0 -epsB 0.3 -n ${N1:-100} -every 5 -ops all ;;
    s2)    run $part -N 8 -p 2 -field abc -b0 -epsB 0.3 -n ${N1:-100} -every 5 -ops all ;;
    s3)    run $part -N 8 -p 2 -field abc -no-b0 -epsB 0.3 -n ${N3:-400} -every 20 -ops A_int,A_pt,B_int,B_l2c ;;
    s3b)   run $part -N 8 -p 2 -field abc -no-b0 -epsB 0.3 -n ${N3:-400} -every 20 -ops A_l2,B_l2 ;;
    *) echo "unknown part $part" ;;
  esac
done
