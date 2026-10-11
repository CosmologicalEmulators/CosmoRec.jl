#!/bin/bash
# Builds the CosmoRec.jl bridge CAMB (scratch). Pattern of ../camb/build_camb_scratch.sh: rsync of the ORIGINAL tools/CAMB-cosmorec (read-only, no .git)
# -> bridge/camb_bridge (NO source patch); COSMOREC_PATH = bridge/shim_lib holding libCosmoRec.a = the C shim only + a no-op `lib` Makefile;
# GSL_LINK empty (nothing of native CosmoRec/GSL is linked). Same make target/flags as the scratch build otherwise. Log: bridge/build_camb_bridge.log.
set -eu
B=/home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk36_camb_replacement/bridge
ORIG=/home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/tools/CAMB-cosmorec
(cd $ORIG && git rev-parse HEAD; git -C $ORIG status --short | wc -l) > $B/orig_camb_state_before.txt
mkdir -p $B/shim_lib
gcc -O2 -fPIC -c $B/shim/cosmorecjl_shim.c -o $B/shim_lib/cosmorecjl_shim.o
rm -f $B/shim_lib/libCosmoRec.a; ar rcs $B/shim_lib/libCosmoRec.a $B/shim_lib/cosmorecjl_shim.o
printf 'lib:\n\t@echo "no-op: CosmoRec.jl bridge shim libCosmoRec.a"\ntidy:\n\t@true\n' > $B/shim_lib/Makefile
rsync -a --delete --exclude .git $ORIG/ $B/camb_bridge/
rm -f $B/camb_bridge/camb/camblib.so
(cd $B/camb_bridge/fortran && make python RECOMBINATION_FILES="recfast cosmorec" COSMOREC_PATH=$B/shim_lib/ GSL_LINK="") > $B/build_camb_bridge.log 2>&1
SO=$B/camb_bridge/camb/camblib.so
sha256sum $B/shim/cosmorecjl_shim.c $B/shim_lib/libCosmoRec.a $SO $ORIG/camb/camblib.so > $B/build_hashes.txt
{ echo "## nm -D (shim exports)"; nm -D --defined-only $SO | grep -i "cosmorec" ;
  echo "## native CosmoRec symbols in camblib.so (expect none)"; nm -C $SO | grep -E "Cosmos::|Xe_frac|CosmoRec\(int|compute_DPesc|gsl_" | wc -l ;
  echo "## ldd (expect no libgsl)"; ldd $SO; } > $B/link_proof.txt
(cd $ORIG && git rev-parse HEAD; git -C $ORIG status --short | wc -l) > $B/orig_camb_state_after.txt
cmp $B/orig_camb_state_before.txt $B/orig_camb_state_after.txt && echo "original CAMB tree unchanged"
echo "build done: $SO"
