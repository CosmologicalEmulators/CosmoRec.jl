# Frozen intended implementation: CosmoRec.jl as CAMB's recombination provider

Marco, 2026-10-11 (verbatim intent): "When I say one to one, I do not mean a floating-point equivalence. The goal is to be able to get the
same input from CAMB as CosmoRec and provide an output for CAMB such that, if we replace CosmoRec with CosmoRec.jl, we see no actual
difference (meaning that the CMB TT and EE up to l=9500 have a difference lower than 0.05%). This should be the final goal!"

This file freezes WHAT is being validated so that campaign results are attributable. It is the operative record for the bridge route.

## The intended implementation (frozen)

Process model (one process, same thread, `OMP_NUM_THREADS=1` — the production worker setting):

    Julia -> PyCall -> python camb -> ctypes camblib.so (Fortran) -> C shim `cosmorec_calc_h_cpp_` -> Julia `@cfunction` -> CosmoRec.jl

1. **Entry**: `cosmorecjl_bridge.jl` `camb_cosmorec_entry` receives exactly CAMB's arguments (Omega_c, Omega_b, Omega_k, N_eff, H0,
   T_CMB, Y_He, z_Hz[10000], Hz[10000], nz, za[m], len, label, runmode, runpars[21]) — CAMB's OWN derived Y_He and N_eff are consumed,
   never reconstructed. It writes Xe/Tm back into CAMB's arrays. Any exception -> status 1 -> the shim ABORTS (no native fallback,
   no zero output).
2. **Science**: `CosmoRec.cosmorec_calc_h` (chunk35 standalone entry): own closure constants, own preliminary history, runmode-0 full
   history. No native-produced initialization data enters a new cosmology's calculation.
3. **Solver callbacks** (verbatim extract of the package's own test callbacks, `production_callbacks.jl`, sha256 a420aede...):
   - preliminary history: `solve_rodas_dense` — Rodas5P dense, reltol 1e-12, abstol [1e-16, 1e-16, 1e-12];
   - feedback passes: `solve_phase5` — **Rodas4P, reltol 1e-10, a1 1e-18, aex 1e-12, tstops = znodes** (the measured production default,
     CMB Tier A verified vs native-history CAMB at l = 2..10000);
   - tail: `solve_tail5` — Rodas5P, reltol 1e-12, abstol [1e-16, 1e-16, 1e-12].
4. **Tables**: user-provided original CosmoRec data tree via `COSMOREC_TABLES_ROOT` (read-only; never copied, never redistributed),
   loaded once per process by `load_runmode0_tables`.
5. **GIL**: Python ctypes releases the GIL around the CAMB call; the callback re-acquires it for its whole duration
   (`PyGILState_Ensure/Release`). Without this, PyCall finalizers (`pydecref`) segfault in `PyObject_Free` (proved: bridge_test_nogil.log).

## Explicitly NOT required for the replacement

- Bitwise/ULP history parity with native, exact native Gear replication, native spline arithmetic. The `RecfastGear.jl` native-parity
  work is PARKED as a diagnostic; it is not on the replacement path and must not be made the default for this campaign.
- The six old chunk10 internal-population gates are diagnostics, not forward-replacement blockers.
- AD/gradient correctness is a SEPARATE scope (own tests; failures stay visible, never hidden).

## Acceptance gates (in force)

- **Focused**: `bridge_test.jl` — real CAMB->shim->Julia->CAMB roundtrip, counters, native-unlinked proof, unregistered-call abort,
  repeated A-B-A calls in one process, bitwise determinism for repeated identical inputs. Must pass before any campaign.
- **Certification**: per cosmology, LENSED TT and EE pointwise `max |C_Julia/C_native - 1| < 5e-4` at EVERY ell 2..9500 (`gate.py`;
  exit contract 0/1/2). TE/BB/PP and unlensed spectra are diagnostics. Ladder: fid + cornerA + cornerB first, then the 36-point matrix.

## Provenance (this freeze)

- Original trees read-only: `tools/CosmoRec` (0867690), `tools/CAMB-cosmorec` (state recorded in orig_camb_state_before/after.txt).
- Bridge build: `build_camb_bridge.sh` — rsync of original CAMB, NO source patch, `COSMOREC_PATH` = shim-only `libCosmoRec.a`,
  `GSL_LINK` empty. `link_proof.txt`: 0 native CosmoRec/GSL symbols in camblib.so, no libgsl in ldd.
- Hashes (`build_hashes.txt`): shim cosmorecjl_shim.c ed9ec23f..., shim libCosmoRec.a f1f62d2a..., bridge camblib.so 598ddf24...,
  original camblib.so e70a6976....
- Julia env: `bridge/julia_env` (Manifest sha fd13c348...), Julia 1.12.6 (`~/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia`).
- Measured bridge behavior (lmax-2000 feasibility run): Julia recombination ~38.8 s/call (lmax-independent), CAMB wall ~40.6 s;
  CAMB re-tabulated x_e vs returned arrays 1.38e-8 max rel (CAMB's own thermo-grid interpolation), T_b 9.3e-12.

## Change discipline

Any change to entry chain, callbacks, tolerances, tables, or build invalidates prior campaign results: bump this file with the new
hashes and re-run the ladder from the focused gate. No gate/reference weakening. No commits/pushes of the working trees.
