#!/bin/bash
# E5: smooth ALE cycle with rezones, N = 8, K = 100 steps, 1 MPI rank per job, JOBS (default 2) jobs at once.
#   lag      : Lagrangian references (rezones 0; all ops identical so A_int only), p=2,3 x eps 0.3,0.6 (q=2),
#              plus p=3 with q=3 (separates the q=2 geometry error) and p=2 eps=0.3 with b0
#   main20/main5/main2/main1 : p=2,3 x eps 0.3,0.6 x rezones 20/5/2/1 x {A_int none, A_int coulomb, A_l2 coulomb, B_int, B_l2c, B_l2}
#   apt      : A_pt (gauge none), p=2, eps 0.3, rezones 1,5,20
#   b0       : mean field B0 on, p=2, eps 0.3, rezones 5, ops B_int, B_l2c, B_l2
#   q3       : p=3 with q=3, eps=0.6, 5 rezones, ops B_l2c, A_int+Coulomb, A_l2+Coulomb, B_l2 (geometry control)
#   gpre     : A_int/A_l2 with the Coulomb gauge applied on the deformed source mesh before the transfer (-gauge coulomb_pre)
#   extra    : A_l2 without gauge fixing (p=2, eps 0.3, rezones 1,5,20)
# (rezones=1 is a rezone at t=T only, i.e. a uniform->uniform identity transfer; rezones=2 is the smallest case with a
#  genuine transfer from a deformed mesh: one at T/2, one at T.)
# Thinned for the time budget: r=1 only for p=2, eps=0.3; B_int p=3 r=20 only for eps=0.3; p=3 logged every 2nd step (every 5th at r=20).
# Usage: scripts/run_exp5.sh [groups...]   (default: lag main20 b0 main5 main2 main1 apt extra)
# Each job appends to results/exp5_parts/<tag>.csv, log in results/exp5_logs/<tag>.log;
# at the end the parts are merged into results/exp5_ale.csv.
set -u
cd "$(dirname "$0")/.."
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
BIN=${BIN:-bin/exp5_ale}
JOBS=${JOBS:-2}
K=${K:-100}
EVERY3=${EVERY3:-2}   # logging stride for p=3 (time budget: diagnostics cost ~1.3 s/step at p=3)
mkdir -p results/exp5_parts results/exp5_logs
GROUPS_=${*:-lag main20 b0 main5 main2 main1 apt extra}
JOBLIST=$(mktemp)

add() {  # tag, args...
  local tag=$1; shift
  [ -s results/exp5_parts/$tag.csv ] && return   # already done (delete the part file to rerun)
  echo "$tag $*" >> $JOBLIST
}

for g in $GROUPS_; do
  case $g in
    main20|main5|main2|main1)
      nrez=${g#main}
      for p in 3 2; do
        for eps in 0.3 0.6; do
          [ $nrez -eq 1 ] && { [ $p -ne 2 ] || [ $eps != 0.3 ]; } && continue          # r=1 (identity transfer): p=2, eps=0.3 only
          for cfg in "B_int none" "B_l2c none" "B_l2 none" "A_int none" "A_int coulomb" "A_l2 coulomb"; do
            set -- $cfg
            [ $nrez -eq 20 ] && [ $p -eq 3 ] && [ $eps = 0.6 ] && [ $1 = B_int ] && continue   # thinned: B_int p=3 r=20 (17 min/run) only for eps=0.3
            ev=1; [ $p -ge 3 ] && ev=$EVERY3     # p=3: log every EVERY3-th step (rezone steps and t=T always logged)
            [ $nrez -eq 20 ] && [ $p -ge 3 ] && ev=5
            add main_p${p}_e${eps}_r${nrez}_$1_$2 -N 8 -p $p -q 2 -eps $eps -K $K -rezones $nrez -op $1 -gauge $2 -every $ev -tag main
          done
        done
      done ;;
    lag)
      for p in 2 3; do for eps in 0.3 0.6; do
        add lag_p${p}_e${eps} -N 8 -p $p -q 2 -eps $eps -K $K -rezones 0 -op A_int -tag lag
      done; done
      for eps in 0.3 0.6; do
        add lagq3_p3_e${eps} -N 8 -p 3 -q 3 -eps $eps -K $K -rezones 0 -op A_int -every $EVERY3 -tag lag_q3
      done
      add lagb0_p2_e0.3 -N 8 -p 2 -q 2 -eps 0.3 -K $K -rezones 0 -op A_int -b0 -tag lag_b0 ;;
    apt)
      for nrez in 20 5 1; do
        add apt_p2_e0.3_r${nrez} -N 8 -p 2 -q 2 -eps 0.3 -K $K -rezones $nrez -op A_pt -gauge none -tag apt
      done ;;
    b0)
      for op in B_int B_l2c B_l2; do
        add b0_p2_e0.3_r5_$op -N 8 -p 2 -q 2 -eps 0.3 -K $K -rezones 5 -op $op -b0 -tag b0
      done ;;
    q3)   # geometry control: p=3 with q=3 (the q=2 mesh interpolation error of the material map dominates at p=3), eps 0.6, r=5
      for cfg in "B_l2c none" "A_int coulomb" "A_l2 coulomb" "B_l2 none"; do
        set -- $cfg
        add q3_p3_e0.6_r5_$1_$2 -N 8 -p 3 -q 3 -eps 0.6 -K $K -rezones 5 -op $1 -gauge $2 -every 5 -tag q3
      done ;;
    gpre)   # extra: Coulomb gauge applied on the deformed source mesh BEFORE the transfer
      for nrez in 20 5 2; do
        for eps in 0.3 0.6; do
          for op in A_int A_l2; do
            add gpre_p2_e${eps}_r${nrez}_${op} -N 8 -p 2 -q 2 -eps $eps -K $K -rezones $nrez -op $op -gauge coulomb_pre -tag gpre
          done
        done
      done
      for nrez in 20 5; do
        for op in A_int A_l2; do
          add gpre_p3_e0.6_r${nrez}_${op} -N 8 -p 3 -q 2 -eps 0.6 -K $K -rezones $nrez -op $op -gauge coulomb_pre -every 5 -tag gpre
        done
      done ;;
    extra)
      for nrez in 20 5 1; do
        add extra_p2_e0.3_r${nrez}_A_l2_none -N 8 -p 2 -q 2 -eps 0.3 -K $K -rezones $nrez -op A_l2 -gauge none -tag extra
      done ;;
    *) echo "unknown group $g" ;;
  esac
done

runjob() {
  tag=$1; shift
  log=results/exp5_logs/$tag.log
  rm -f results/exp5_parts/$tag.csv.tmp
  t0=$(date +%s)
  mpirun -np 1 --oversubscribe $BIN "$@" -csv results/exp5_parts/$tag.csv.tmp > $log 2>&1
  rc=$?
  if [ $rc -eq 0 ]; then mv results/exp5_parts/$tag.csv.tmp results/exp5_parts/$tag.csv; fi
  echo "[$tag] rc=$rc $(( $(date +%s) - t0 )) s" | tee -a results/exp5_logs/_summary.log
}
export -f runjob
export BIN
echo "$(wc -l < $JOBLIST) jobs, $JOBS at a time"
cat $JOBLIST | xargs -P $JOBS -L 1 bash -c 'runjob "$@"' _
rm -f $JOBLIST
python3 scripts/merge_csv.py results/exp5_ale.csv results/exp5_parts/*.csv
