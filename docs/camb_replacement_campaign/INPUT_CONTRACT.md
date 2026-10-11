# CAMB -> recombination-provider input contract: bitwise proof (2026-10-11)

**Question (Marco):** are the outputs CAMB passes INTO CosmoRec.jl the same values native CosmoRec
receives, or does the bridge route feed it something different?

**Answer: bitwise identical, proven directly.** On both input paths (explicit Y_He, and Y_He
derived from CAMB's BBN table), every value CAMB hands the provider is byte-for-byte the same on
the two routes. Additionally, isolating the provider swap behind byte-identical CAMB binaries
shows the pure CosmoRec.jl-vs-CosmoRec spectra difference is sub-ppm (<= 9.5e-7 at fid and
deg_bnd1), i.e. the campaign-level ~1e-6..1.2e-5 numbers vs the original native camblib contain a
separate, bounded camblib-rebuild component (<= 2.6e-6 on TT/EE) that is not attributable to the
recombination swap.

## Mechanism

The bridge camblib is the original tools/CAMB-cosmorec source (clean tree, commit fa3f097)
rebuilt with the same make recipe, swapping ONLY the libCosmoRec archive for the shim
(bridge/shim/cosmorecjl_shim.c). The shim was extended (2026-10-11) with an input-contract mode:

- `COSMORECJL_DUMP=<file>`: before dispatching, the shim writes every CAMB->provider input as one
  raw dump. The SAME C code writes the dump on both routes, so the two files are comparable with
  `cmp` (no cross-language layout risk).
- `COSMORECJL_NATIVE_LIB=<libCosmoRec.so>`: the shim forwards the call, all pointers unchanged,
  to the ORIGINAL native CosmoRec via dlopen. The .so is built from the unmodified original
  tools/CosmoRec/libCosmoRec.a (commit 0867690, the CosmoRec.jl port source) with
  -Wl,--whole-archive against the same libgsl.so.28 (tools/camb-cosmorec-env/lib) the original
  camblib resolves. Without the env var the shim behaves exactly as before (Julia entry, abort
  on unregistered call / nonzero status; the native mode keeps the same abort discipline).

So both routes run behind the **byte-identical camblib binary**; the only difference is which
library answers `cosmorec_calc_h_cpp_` (CosmoRec.h:153-158; called from cosmorec.f90:211-213).

## Dump layout (240253 bytes at fid/deg_bnd1)

magic "CRJLDUMP1" (9 B) | int32 ver=1 | int32 runmode | float64 runpars[21] | float64
omega_c, omega_b, omega_k, num_nu, h0, t_cmb, y_he (7) | int32 nz | float64 z_Hz[nz] | float64
Hz[nz] | int32 len | float64 za_in[len] | int32 label.  (nz = len = 10000; zrec is passed twice:
as the H-table grid and as the output grid, cosmorec.f90:211-213; Hz is CAMB's own background
solver output.)

## Results

| check | fid (explicit Y_He) | deg_bnd1 (BBN-derived Y_He) |
|---|---|---|
| input dump cmp (native vs julia route) | **bitwise identical** (240253 B) | **bitwise identical** |
| pure swap gate: julia vs native-behind-bridge-camblib, lensed TT / EE | 4.641e-7 / 4.415e-7 | PASS (gate) |
| ... unlensed TT / EE (worst) | 4.893e-7 / 9.537e-7 | <= 8.7e-7 (TE diag) |
| TE/BB/PP diagnostics | <= 3.6e-7 | <= 8.7e-7 |

Decomposition at fid (max relative, lensed TT/EE):
- julia vs original-native camblib (campaign numbers): ~3e-6 level
- camblib-rebuild effect (native-behind-bridge-camblib vs original-native camblib, same CosmoRec
  objects): 2.2e-6 / 2.6e-6 lensed, 1.5e-6 / 2.0e-6 unlensed TT/EE; larger *relative* values only
  on near-zero columns (unlensed BB/TE zero crossings, abs <= 3.3e-8). Same source (clean trees,
  fa3f097 + 0867690), same-day toolchain: a build-level artifact, bounded 200x under the 5e-4
  budget, NOT attributable to the recombination swap.
- pure provider swap (identical camblib): 4.4e-7 / 4.6e-7 lensed.

Both routes' bridge input logs (`fid_julia_bridge.txt` etc.), per-point metas (params +
cosmorec_cfg identical) and the unregistered-call abort remain the secondary evidence; this dump
comparison is the direct proof.

## Artifacts (scratch, not in repo)

- bridge/natlog/libCosmoRec.so (from tools/CosmoRec/libCosmoRec.a + original libgsl.so.28)
- bridge/cmb_pair/input_dump_{native,julia}.bin (fid), input_dump_{native,julia}_bbn.bin (deg_bnd1)
- bridge/cmb_pair/{fid,deg_bnd1}_natlog_cls.txt + logs; gate output in the run logs
