#!/bin/bash
# Usage: run_one.sh <name> <rs> <tf> <orders: hi|lo> <ale-opts...>
# Runs one 3-D Taylor-Green MHD-ALE job on 2 MPI ranks; output in external/MHD-ALE/output/E6/<name>
# hi orders: -ok 2 -ot 1 -or 1 -oa 3     lo orders: -ok 1 -ot 0 -or 0 -oa 2
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
EXE=${EXE:-./MHD}
cd /home/user/mfem-vector-potential/external/MHD-ALE || exit 1
NAME=$1; RS=$2; TF=$3; ORD=$4; shift 4
if [ "$ORD" = hi ]; then ORDS="-ok 2 -ot 1 -or 1 -oa 3"; else ORDS="-ok 1 -ot 0 -or 0 -oa 2"; fi
OD=output/E6/$NAME; mkdir -p $OD
CMD="mpirun --oversubscribe -np 2 $EXE -p 0 -dim 3 -m ./mesh/Cube-4x4x4-hex.mesh -rs $RS -tf $TF $ORDS -s 3 -cfl 0.5 -g1 0.85 -g2 1.02 -mst 1 -rmv 1 -rme 1 -rmr 3 -cgt 1e-12 -bpt 1 -ptf plot_time/plot_time_TaylorGreen.dat -od $OD $@"
echo "$CMD" > $OD/cmd.txt
T0=$(date +%s.%N)
$CMD > $OD/MHD.out 2> $OD/MHD.err
T1=$(date +%s.%N)
echo "$(echo "$T1 - $T0" | bc -l) s wall" > $OD/walltime.txt
echo "$NAME done: $(cat $OD/walltime.txt)"
