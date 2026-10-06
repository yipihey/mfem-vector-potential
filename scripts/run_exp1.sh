#!/bin/bash
# E1: static convergence on the undeformed periodic box, 1 MPI rank.
# p=1..4, N=4,8,16 (+32 for p<=2), fields abc and mod2, pointwise and
# integrated projections; skip combinations with > 3e6 ND dofs.  Plus a short
# series with the mean field B0 (-b0) at N=4,8.  q=1: the undeformed box is
# affine, the geometric order has no influence on the results (only on cost).
# Usage: scripts/run_exp1.sh [out.csv]     (existing file is replaced)
set -u
cd "$(dirname "$0")/.."
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
OUT=${1:-results/exp1_static.csv}
mkdir -p results
rm -f "$OUT"
BIN=bin/exp1_static
nd_dofs() { python3 -c "
p=$1; N=$2
print(N**3*(3*p+6*p*(p-1)+3*p*(p-1)**2))"; }
run() {  # p N field proj [extra]
  local p=$1 N=$2 field=$3 proj=$4; shift 4
  local nd; nd=$(nd_dofs "$p" "$N")
  if [ "$nd" -gt 3000000 ]; then echo "skip p=$p N=$N (ND dofs $nd)"; return; fi
  echo "== exp1 p=$p N=$N field=$field proj=$proj $*"
  mpirun -np 1 $BIN -N "$N" -p "$p" -q 1 -field "$field" -proj "$proj" -csv "$OUT" "$@" > /tmp/exp1_last.log 2>&1 \
    || { echo "FAILED"; tail -5 /tmp/exp1_last.log; }
}
for p in 1 2 3 4; do
  Ns="4 8 16"; if [ "$p" -le 2 ]; then Ns="4 8 16 32"; fi
  for N in $Ns; do
    for field in abc mod2; do
      for proj in pt int; do run "$p" "$N" "$field" "$proj"; done
    done
  done
done
for p in 1 2 3 4; do
  for N in 4 8; do
    for proj in pt int; do run "$p" "$N" abc "$proj" -b0; done
  done
done
echo "done -> $OUT"
