#!/bin/bash
# E3: one remap M1 -> M2.  1 MPI rank per run (the machine is shared; two parts can run concurrently).
#   (a) M1 = f1 eps in {0.1,0.3,0.5} -> M2 uniform, p=1..4 (N=8) and p<=2 (N=16), abc + B0, all 7 ops
#   (b) M1 = f1 eps 0.3 -> M2 = f2 eps 0.15, p=1..3, N=8,16, mod2, no B0, all ops
#   (c) gauge sensitivity: (a) with eps 0.3, N=8, p=1..3, g=0.5, k=2pi and 6pi, A-route ops
#   (d) resolution change: M1 = f1 eps 0.3 N=8 -> M2 uniform N=12, 16, p=1..3, abc + B0
#   (e) quadrature study: nq in {2,3,4,6,8}, ops A_int,B_int, p=1..3, N=8, eps 0.3 -> uniform
# Usage: scripts/run_exp3.sh [parts...]      (default: a b c d e, serial)
#        e.g. scripts/run_exp3.sh a c e   and   scripts/run_exp3.sh b d   as two concurrent jobs
# Each part appends to results/exp3_parts/<part>.csv; scripts/merge_csv.py merges into results/exp3_remap.csv
set -u
cd "$(dirname "$0")/.."
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
BIN=${BIN:-bin/exp3_remap}
NP=${NP:-1}
mkdir -p results/exp3_parts results/exp3_logs
PARTS=${*:-a b c d e}

run() {  # part, args...
  local part=$1; shift
  local log=results/exp3_logs/$part.log
  echo "== [$part] $*" | tee -a "$log"
  local t0=$(date +%s)
  mpirun -np $NP --oversubscribe $BIN "$@" -tag "$part" -csv results/exp3_parts/$part.csv >> "$log" 2>&1
  local rc=$?
  echo "   rc=$rc  $(( $(date +%s) - t0 )) s" | tee -a "$log"
}

for part in $PARTS; do
  case $part in
    a)
      for eps in 0.1 0.3 0.5; do
        for p in 1 2 3 4; do run a -N 8 -p $p -field abc -b0 -eps1 $eps -m2 uniform -ops all; done
      done
      for eps in 0.1 0.3 0.5; do
        for p in 1 2; do run a -N 16 -p $p -field abc -b0 -eps1 $eps -m2 uniform -ops all; done
      done ;;
    b)
      for N in 8 16; do
        for p in 1 2 3; do run b -N $N -p $p -field mod2 -no-b0 -eps1 0.3 -m2 deformed:0.15:2 -ops all; done
      done ;;
    c)
      for gk in 2 6; do
        for p in 1 2 3; do run c -N 8 -p $p -field abc -b0 -eps1 0.3 -m2 uniform -gauge 0.5 -gk $gk -ops A_pt,A_int,A_l2; done
      done ;;
    d)
      for N2 in 12 16; do
        for p in 1 2 3; do run d -N 8 -p $p -field abc -b0 -eps1 0.3 -m2 fine:$N2 -ops all; done
      done ;;
    e)
      for p in 1 2 3; do
        for nq in 2 3 4 6 8; do run e -N 8 -p $p -field abc -b0 -eps1 0.3 -m2 uniform -nq $nq -ops A_int,B_int -no-gradfrac; done
      done ;;
    *) echo "unknown part $part" ;;
  esac
done
