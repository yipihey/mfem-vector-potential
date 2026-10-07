#!/bin/bash
# tf = 0.25 set: the fsri-2 -rma 2 run costs ~100 s per remap on the shared machine, so tf=0.5 (14 remaps) was aborted
R=/home/user/mfem-vector-potential/mhdale/run_one.sh
$R fsri2_rma0_rs1_tf0.25 1 0.25 hi -ale -fsr -fsri 2 -rma 0
$R noale_rs1_tf0.25 1 0.25 hi -no-ale
$R fsri2_rma2_rs1_tf0.25 1 0.25 hi -ale -fsr -fsri 2 -rma 2
