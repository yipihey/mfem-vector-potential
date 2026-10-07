#!/bin/bash
# E4-gauge: repeated remaps (ping-pong uniform <-> f1 eps 0.3, N=8, abc) with a gauge fix after every
# A-route transfer, plus the one-remap exp3 variant and a contemporaneous cost comparison.  1 MPI rank/run.
#   g_p2   : p=2, n=100, every 5, ops A_pt,A_int,A_l2, -gauge coulomb
#   g_p3   : p=3, same
#   g_long : p=2, n=400, every 20, ops A_int,A_l2, -gauge coulomb
#   g_jac  : p=2, n=100, A_int, -gauge jacobi:5 and jacobi:20           (experimental local variant)
#   g_b0   : p=2, n=100, A_int, -gauge coulomb, WITH the mean field B0
#   g_cost : n=20, p=2,3: B_l2c, and A_int/A_l2 with -gauge none and -gauge coulomb (same machine state, run serially)
#   e3     : exp3 one remap, N=8, eps 0.3 -> uniform, abc + B0, p=1..3, A_pt,A_int,A_l2, -gfix coulomb
# Usage: scripts/run_exp4_gauge.sh [parts...]   (default: all, serial; at most 2 concurrent invocations please)
# Parts append to results/exp4_gauge_parts/<part>.csv; scripts/merge_csv.py merges (see end of this script's usage
# notes in results/exp4_gauge_summary.md):  results/exp4_gauge.csv  (exp4 parts)  and  results/exp3_gfix.csv (e3).
set -u
cd "$(dirname "$0")/.."
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
NP=${NP:-1}
D=results/exp4_gauge_parts
mkdir -p $D results/exp4_gauge_logs
PARTS=${*:-g_p2 g_p3 g_long g_jac g_b0 g_cost e3}

run() {  # part, binary, args...
  local part=$1 bin=$2; shift 2
  local log=results/exp4_gauge_logs/$part.log
  echo "== [$part] $bin $*" | tee -a "$log"
  local t0=$(date +%s)
  mpirun -np $NP --oversubscribe bin/$bin "$@" -tag "$part" -csv $D/$part.csv >> "$log" 2>&1
  local rc=$?
  echo "   rc=$rc  $(( $(date +%s) - t0 )) s" | tee -a "$log"
}
E4="-N 8 -field abc -epsB 0.3"

for part in $PARTS; do
  case $part in
    g_p2)   run $part exp4_repeat $E4 -p 2 -no-b0 -n 100 -every 5 -ops A_pt,A_int,A_l2 -gauge coulomb ;;
    g_p3)   run $part exp4_repeat $E4 -p 3 -no-b0 -n 100 -every 5 -ops A_pt,A_int,A_l2 -gauge coulomb ;;
    g_long) run $part exp4_repeat $E4 -p 2 -no-b0 -n 400 -every 20 -ops A_int,A_l2 -gauge coulomb ;;
    g_jac)
      run $part exp4_repeat $E4 -p 2 -no-b0 -n 100 -every 5 -ops A_int -gauge jacobi:5
      run $part exp4_repeat $E4 -p 2 -no-b0 -n 100 -every 5 -ops A_int -gauge jacobi:20 ;;
    g_b0)   run $part exp4_repeat $E4 -p 2 -b0 -n 100 -every 5 -ops A_int -gauge coulomb ;;
    g_cost)
      for p in 2 3; do
        run $part exp4_repeat $E4 -p $p -no-b0 -n 20 -every 4 -ops B_l2c
        run $part exp4_repeat $E4 -p $p -no-b0 -n 20 -every 4 -ops A_int,A_l2 -gauge none
        run $part exp4_repeat $E4 -p $p -no-b0 -n 20 -every 4 -ops A_int,A_l2 -gauge coulomb
      done ;;
    e3)
      for p in 1 2 3; do
        run $part exp3_remap -N 8 -p $p -field abc -b0 -eps1 0.3 -m2 uniform -ops A_pt,A_int,A_l2 -gfix coulomb
      done ;;
    *) echo "unknown part $part" ;;
  esac
done
