#!/bin/bash
# Runs bin/selftest (FD curl check; D_h C_h = 0 and C_h G_h = 0 on eps=0.25
# deformed meshes for p=1..4; constant-B0 exactness + slice fluxes; integrated
# projection commutation; helicity/energy) on 1 rank (and optionally more:
# NP="1 2 4" scripts/selftest.sh) and prints PASS/FAIL.
cd "$(dirname "$0")/.."
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
make -s bin/selftest || { echo "FAIL (build)"; exit 1; }
rc_all=0
for np in ${NP:-1}; do
  echo "=== selftest, $np MPI rank(s) ==="
  mpirun -np "$np" --oversubscribe bin/selftest > /tmp/vp_selftest_$np.log 2>&1
  rc=$?
  grep -E "^(FAIL|INFO)|ALL PASS|SOME FAIL" /tmp/vp_selftest_$np.log
  grep -c "^PASS" /tmp/vp_selftest_$np.log | sed 's/^/PASS checks: /'
  if [ $rc -ne 0 ]; then rc_all=1; fi
done
if [ $rc_all -eq 0 ]; then echo "SELFTEST: PASS"; else echo "SELFTEST: FAIL (see /tmp/vp_selftest_*.log)"; fi
exit $rc_all
