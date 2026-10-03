# Chunk 5d acceptance: output assembly and the returned history

## Gates
- `test/chunk5d_output_native.jl`: **13/13** (focused). Full root `Pkg.test()` (`COSMOREC_NATIVE_DATA_DIR` set; tests registered through 5d): **49407/49407**, exit 0, 22m35.7s, `chunk5/pkgtest5d.log` (49394 + 13).
- Native comparison: native rows through the Julia spline/fallback code reproduce the native returned `Xe_arr`, `Te_arr` on all 10000 CAMB-grid nodes to 2.2e-16 (assembly is exact); fallbacks (z >= 3000, z = 0) 2.1e-16; the Julia pass inside the stored range differs from the native by Xe 1.6e-7, Te 3.0e-8 (native solver tolerance level).
- Benchmarks: assembly 0.78 ms; complete single pass 798 ms.
## Accepted scope
Native `output_CosmoRec` rows and `return_solution_to_calling_program` (GSL natural cubic splines of ln Xe, ln Te with native end nudge, preliminary-history fallback outside the stored range), plus the complete single pass `recombination_history`. This is native iteration 0 (`runmode 1`), NOT the production CAMB history (two further diffusion iterations).
## Limitations
One cosmology/table; native solver tolerance level agreement; machine-specific data path (D1).
## Still not a full CosmoRec model
Gradients (5e), diffusion/feedback (phase 6+), runmode 0 parity with `native_camb_cosmorec_thermo.txt`.
