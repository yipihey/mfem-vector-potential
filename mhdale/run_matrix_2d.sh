#!/bin/bash
# E6 side note: 2-D TG (scalar H1 potential A), Disk-4x4-quad, -rs 2, -ok 2 -ot 1 -or 1 -oa 3
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
cd /home/user/mfem-vector-potential/external/MHD-ALE || exit 1
TF=${TF:-0.3}
run() { # name, extra opts
  NAME=$1; shift
  OD=output/E6_2d/$NAME; mkdir -p $OD
  CMD="mpirun --oversubscribe -np 2 ./MHD -p 0 -dim 2 -m ./mesh/Disk-4x4-quad.mesh -rs 2 -tf $TF -ok 2 -ot 1 -or 1 -oa 3 -s 3 -cfl 0.5 -g1 0.85 -g2 1.02 -mst 1 -rmv 1 -rme 1 -rmr 3 -cgt 1e-12 -bpt 1 -ptf plot_time/plot_time_TaylorGreen.dat -od $OD $@"
  echo "$CMD" > $OD/cmd.txt
  T0=$(date +%s.%N); $CMD > $OD/MHD.out 2> $OD/MHD.err; T1=$(date +%s.%N)
  echo "$(echo "$T1 - $T0" | bc -l) s wall" > $OD/walltime.txt
  echo "$NAME done: $(cat $OD/walltime.txt)"
}
run noale_rs2_tf$TF -no-ale
for I in 10 2; do for A in 0 1; do run fsri${I}_rma${A}_rs2_tf$TF -ale -fsr -fsri $I -rma $A; done; done
