#!/bin/bash
# E6 main matrix: 3-D TG, -rs 1, orders -ok 2 -ot 1 -or 1 -oa 3, tf=0.5
R=/home/user/mfem-vector-potential/mhdale/run_one.sh
TF=${TF:-0.5}
$R noale_rs1_tf$TF 1 $TF hi -no-ale
$R fsri10_rma0_rs1_tf$TF 1 $TF hi -ale -fsr -fsri 10 -rma 0
$R fsri10_rma2_rs1_tf$TF 1 $TF hi -ale -fsr -fsri 10 -rma 2
$R fsri2_rma0_rs1_tf$TF 1 $TF hi -ale -fsr -fsri 2 -rma 0
$R fsri2_rma2_rs1_tf$TF 1 $TF hi -ale -fsr -fsri 2 -rma 2
