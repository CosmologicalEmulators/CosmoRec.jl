#!/bin/bash
# Paired CMB certification, SERIAL. Native = ORIGINAL tools/CAMB-cosmorec + original libCosmoRec (unmodified; same CAMB source and make
# recipe as the bridge build). Julia = bridge CAMB + CosmoRec.jl. Same runner/preset (../camb/camb_forward_c36.py, preset worker =
# PRODUCTION settings, lmax 9500), same python, same points, OMP 1, same tables root.
# EXIT CONTRACT: 0 = every point's lensed TT AND EE passed the pointwise 5e-4 gate at every ell 2..9500 (gate.py);
# nonzero = failed or broken (see *_native.exit, julia_points.exit, *_gate.exit, gate.log). ALL_DONE is written ONLY on success.
# usage: run_cmb_pair.sh [points...]   (default: fid cornerA cornerB)
set -uo pipefail
C36=/home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk36_camb_replacement
D=$C36/bridge/cmb_pair
ORIG=/home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/tools/CAMB-cosmorec
PY=/home/marcobonici/miniconda3/bin/python
J=/home/marcobonici/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia
export OMP_NUM_THREADS=1 COSMOREC_TABLES_ROOT=/home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/tools/CosmoRec
PTS=${@:-fid cornerA cornerB}
rm -f $D/ALL_DONE $D/ALL_FAILED
{
  echo "campaign start $(date -Is) host $(hostname) points: $PTS"
  sha256sum $C36/bridge/camb_bridge/camb/camblib.so $ORIG/camb/camblib.so $C36/bridge/shim/cosmorecjl_shim.c
  $J --version; $PY --version 2>&1
} > $D/run_header.txt 2>&1
rc=0
for p in $PTS; do
  echo "[native $p] start $(date -Is)"
  if PYTHONPATH=$ORIG /usr/bin/time -v $PY $C36/camb/camb_forward_c36.py $D/$p.json $D/${p}_native native worker > $D/${p}_native.log 2>&1; then
    echo 0 > $D/${p}_native.exit
  else
    g=$?; echo $g > $D/${p}_native.exit; echo "NATIVE RUN FAILED: $p exit $g (see ${p}_native.log)"; rc=1
  fi
done
echo "[julia] start $(date -Is)"
if /usr/bin/time -v $J --threads=1 --project=$C36/bridge/julia_env $D/run_julia_points.jl $PTS > $D/julia_points.log 2>&1; then
  echo 0 > $D/julia_points.exit
else
  g=$?; echo $g > $D/julia_points.exit; echo "JULIA RUN FAILED exit $g (see julia_points.log)"; rc=1
fi
for p in $PTS; do
  if [ -f $D/${p}_julia_cls.txt ] && [ -f $D/${p}_native_cls.txt ]; then
    $PY $C36/camb/cl_compare.py $D ${p}_julia ${p}_native "$p Julia-bridge vs native" >> $D/compare.log 2>&1
    echo "cl_compare exit $? (diagnostic only)" >> $D/compare.log
    if $PY $D/gate.py $D ${p}_julia ${p}_native >> $D/gate.log 2>&1; then
      echo 0 > $D/${p}_gate.exit
    else
      g=$?; echo $g > $D/${p}_gate.exit; echo "GATE FAIL: $p exit $g (see gate.log)"; rc=1
    fi
  else
    echo 2 > $D/${p}_gate.exit; echo "GATE BROKEN: $p missing outputs"; rc=1
  fi
done
echo "campaign end $(date -Is) rc=$rc" >> $D/run_header.txt
if [ $rc -eq 0 ]; then echo done > $D/ALL_DONE; echo "CAMPAIGN PASS: $PTS"; else echo fail > $D/ALL_FAILED; echo "CAMPAIGN FAIL"; fi
exit $rc
