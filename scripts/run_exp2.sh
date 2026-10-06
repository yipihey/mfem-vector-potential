#!/bin/bash
# E2: distorted periodic meshes, 1 MPI rank per run.
#  * main matrix (integrated projection, field abc, no B0):
#      p=1..4 at N=8 (and N=16 for p<=2); geometric order q in {1,2,p} (deduped);
#      f1: eps = 0,0.05,...,0.3,0.4,...,1.0 ; f2 and f3: eps = 0.1,0.2,0.3,0.4,0.6,0.8,1.0
#      (reduced lists for p=4: f1 0,0.1,0.2,0.3,0.4,0.6,0.8; f2,f3 0.2,0.4,0.8).
#      A series stops when a mesh has non-positive Jacobians (recorded as
#      valid=0) or min detJ < 0.15 (normalised) -- that is the "invalid eps".
#  * extras: pointwise-projection series (N=8, q=2, f1 and f3); abc + mean
#    field B0 series (N=8, q=2; slice flux, Pi_RT(B0) exactness, integrated and
#    pointwise B0 projection); mod2 (non-Beltrami, multi-mode) on f1 (p<=3).
#
# Usage:
#   scripts/run_exp2.sh [out.csv]            # everything; p-sets {1,2} and {3,4} run as 2 concurrent jobs
#   PLIST="3 4" STAGES="main extras" OUT=part.csv scripts/run_exp2.sh   # a subset, serial, APPENDS to OUT
# Env: PLIST (p values), STAGES ("main", "extras"), OUT, BIN (executable).
set -u
cd "$(dirname "$0")/.."
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
BIN=${BIN:-bin/exp2_distort}
mkdir -p results

series() {  # csv p N q variant proj field "eps list" [extra args]
  local csv=$1 p=$2 N=$3 q=$4 fv=$5 proj=$6 field=$7 epss=$8; shift 8
  local log; log=$(mktemp)
  for e in $epss; do
    echo "== exp2 p=$p N=$N q=$q f$fv eps=$e proj=$proj field=$field $*"
    mpirun -np 1 $BIN -N "$N" -p "$p" -q "$q" -fvariant "$fv" -eps "$e" -proj "$proj" \
        -field "$field" -csv "$csv" "$@" > "$log" 2>&1
    local rc=$?
    if [ $rc -ge 3 ] && [ $rc -le 4 ]; then echo "   stop series (rc=$rc: degenerate/invalid mesh)"; break; fi
    if [ $rc -ne 0 ]; then echo "   FAILED rc=$rc"; tail -5 "$log"; break; fi
  done
  rm -f "$log"
}

lists() {  # p -> sets F1 F23 FB (b0 f1) FB3 (b0 f3)
  if [ "$1" -ge 4 ]; then
    F1="0 0.1 0.2 0.3 0.4 0.6 0.8"; F23="0.2 0.4 0.8"; FB="0 0.2 0.4 0.6"; FB3="0.4"; FP="0 0.1 0.2 0.3 0.4 0.6"
  else
    F1="0 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"; F23="0.1 0.2 0.3 0.4 0.6 0.8 1.0"
    FB="0 0.1 0.2 0.3 0.4 0.5 0.6"; FB3="0.2 0.4 0.6"; FP="$F1"
  fi
}

run_plist() {  # csv stages p...
  local csv=$1 stages=$2; shift 2
  for p in "$@"; do
    lists "$p"
    if [[ " $stages " == *" main "* ]]; then
      local Ns="8"; if [ "$p" -le 2 ]; then Ns="8 16"; fi
      for N in $Ns; do
        local qs; qs=$(printf "%s\n" 1 2 "$p" | sort -nu | tr '\n' ' ')
        for q in $qs; do
          series "$csv" "$p" "$N" "$q" 1 int abc "$F1"
          series "$csv" "$p" "$N" "$q" 2 int abc "$F23"
          series "$csv" "$p" "$N" "$q" 3 int abc "$F23"
        done
      done
    fi
    if [[ " $stages " == *" extras "* ]]; then
      series "$csv" "$p" 8 2 1 pt abc "$FP"
      series "$csv" "$p" 8 2 3 pt abc "$F23"
      for proj in int pt; do
        series "$csv" "$p" 8 2 1 "$proj" abc "$FB" -b0
        series "$csv" "$p" 8 2 3 "$proj" abc "$FB3" -b0
      done
      if [ "$p" -le 3 ]; then series "$csv" "$p" 8 2 1 int mod2 "0 0.1 0.2 0.3 0.4 0.5 0.6"; fi
    fi
  done
}

if [ -n "${PLIST:-}" ]; then
  run_plist "${OUT:-results/exp2_distort.csv}" "${STAGES:-main extras}" $PLIST
else
  OUT=${1:-results/exp2_distort.csv}
  TMP1=results/.exp2_part_a.csv; TMP2=results/.exp2_part_b.csv
  rm -f "$TMP1" "$TMP2" "$OUT"
  run_plist "$TMP1" "main extras" 1 2 &
  run_plist "$TMP2" "main extras" 3 4 &
  wait
  cat "$TMP1" > "$OUT"; tail -n +2 "$TMP2" >> "$OUT"
  rm -f "$TMP1" "$TMP2"
  echo "done -> $OUT"
fi
