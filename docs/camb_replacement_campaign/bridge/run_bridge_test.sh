#!/bin/bash
# Focused bridge test, one serial process. Log bridge_test.log, exit code in bridge_test.exit.
B=/home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk36_camb_replacement/bridge
export OMP_NUM_THREADS=1 COSMOREC_TABLES_ROOT=/home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/tools/CosmoRec
/usr/bin/time -v /home/marcobonici/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia --threads=1 --project=$B/julia_env $B/bridge_test.jl > $B/bridge_test.log 2>&1
echo $? > $B/bridge_test.exit
