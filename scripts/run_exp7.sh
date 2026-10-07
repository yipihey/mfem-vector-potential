#!/bin/bash
# run_exp7.sh -- E7 timing matrix.  Usage:  scripts/run_exp7.sh STAGE [nowait]
#   stage 1: p=1 (N=16,32), p=2 N=8,16, p=3 N=8,16, p=4 N=8   (ranks 1/2/4 for p2N16, p3N16)
#   stage 2: p=2 N=24 (ranks 1,2,4)
#   stage 3: p=2 N=32 and p=4 N=12 (only if stage 2 shows it is feasible), 1 rank
# Every run is done twice (tags r1, r2; r2 skipped if r1 took > 300 s); the analysis keeps the min.
# The matrix starts only when no exp5_ale process has been seen for 2 consecutive checks 60 s apart
# (poll up to 90 min; afterwards runs with 1 rank only).
cd "$(dirname "$0")/.."
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 OMP_NUM_THREADS=1
STAGE=${1:-1}
CSV=${CSV:-results/exp7_perf.csv}
BCSV=${BCSV:-results/exp7_baseline.csv}
LOGD=${LOGD:-results/exp7_logs}
mkdir -p $LOGD

if [ "$2" != "nowait" ]; then
  idle=0; waited=0
  while [ $idle -lt 2 ] && [ $waited -lt 5400 ]; do
    if pgrep -f exp5_ale > /dev/null; then idle=0; else idle=$((idle+1)); fi
    [ $idle -lt 2 ] && sleep 60 && waited=$((waited+60))
  done
  if [ $idle -lt 2 ]; then
    echo "WARNING: exp5_ale still running after 90 min; running with 1 rank only" | tee -a $LOGD/_wait.log
    ONLY1=1
  else
    echo "idle after ${waited}s wait" | tee -a $LOGD/_wait.log
  fi
fi

OPS=("A_pt none" "A_int none" "A_int coulomb" "A_int jacobi:20" "A_l2 coulomb" "B_l2 none" "B_l2c none")

run1() {  # N p np op gauge tag
  local N=$1 p=$2 np=$3 op=$4 g=$5 tag=$6
  local log=$LOGD/N${N}_p${p}_np${np}_${op}_${g/:/}_${tag}.log
  local t0=$(date +%s)
  mpirun -np $np --oversubscribe bin/exp7_perf -N $N -p $p -op $op -gauge $g -tag $tag -csv $CSV > $log 2>&1
  local rc=$?
  local dt=$(( $(date +%s) - t0 ))
  echo "N=$N p=$p np=$np $op $g $tag rc=$rc wall=${dt}s  $(tail -1 $log | cut -c1-200)" | tee -a $LOGD/_progress.log
  LASTWALL=$dt
}

runcfg() {  # N p np
  local N=$1 p=$2 np=$3
  [ -n "$ONLY1" ] && [ $np -gt 1 ] && return
  for tag in r1 r2; do
    mpirun -np $np --oversubscribe bin/exp7_perf -N $N -p $p -op none -tag $tag -csv $BCSV >> $LOGD/_baseline.log 2>&1
  done
  for o in "${OPS[@]}"; do
    set -- $o
    run1 $N $p $np $1 $2 r1
    if [ $LASTWALL -le 300 ]; then run1 $N $p $np $1 $2 r2; fi
  done
}

case $STAGE in
 1)
  for c in "16 1" "32 1"; do set -- $c; runcfg $1 1 1; done
  runcfg 8 2 1; runcfg 8 3 1; runcfg 8 4 1
  for np in 1 2 4; do runcfg 16 2 $np; done
  for np in 1 2 4; do runcfg 16 3 $np; done
  ;;
 2) for np in 1 2 4; do runcfg 24 2 $np; done ;;
 3) runcfg 32 2 1; runcfg 12 4 1 ;;
 0) runcfg 4 2 1; runcfg 4 2 2; runcfg 4 1 1 ;;   # smoke test (use CSV=/tmp/.. LOGD=/tmp/..)
esac
echo "stage $STAGE done" | tee -a $LOGD/_progress.log
