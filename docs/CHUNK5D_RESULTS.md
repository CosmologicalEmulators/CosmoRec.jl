# Chunk 5d results: output assembly and the returned history

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al.; cite Chluba & Thomas 2010, Chluba & Sunyaev 2006, Rubino-Martin et al. 2008. No native source, library or table was changed or bundled; the tests require `COSMOREC_NATIVE_DATA_DIR` (D1).

## Scope
`recombination_output_rows`, `return_solution_to_grid`, `recombination_history` (`src/RecombinationODE.jl`): the native `output_CosmoRec` rows (z, Xe, Te: 3000 ODE nodes followed by tail rows 2..200 = 3199 rows) and `return_solution_to_calling_program` (`CosmoRec.cpp:545-597`): GSL natural cubic splines (4b) of `ln Xe` and `ln Te` against ascending z, evaluated with the native end nudge for `zmin < z < zmax` (stored range 0.001 < z < 3000); at or beyond the ends the Recfast preliminary history `Xe_Seager(z)`, `Te(z)` (4b accessors). `recombination_history` runs the whole single pass (pass + tail + assembly) for a CAMB-like grid.

## Native oracle
The native `Xe_arr`, `Te_arr` returned on the CAMB 10000-node grid (z = 1e4 -> 0) by the original `cosmorec_calc_h_cpp_` with `runmode = 1` (fixture `native_ode_pass_grid.txt`), plus the native stored rows.

## Results (focused 13/13)
- Stored rows: ODE `z` bitwise, tail `z` to 1e-13 (0.0 observed); 3199 distinct descending nodes.
- Spline stage in isolation (NATIVE rows in, Julia splines/fallback out): Xe 2.2e-16, Te 2.0e-16 relative to the native returned arrays (all 10000 grid nodes, inside and outside the stored range): the assembly is an exact port.
- Fallbacks (7001 grid nodes: z >= 3000 and z = 0): Xe 2.1e-16, Te 0.0.
- Inside (2999 grid nodes): Julia pass rows vs native: Xe max 1.6e-7 (z = 916), Te max 3.0e-8 (z = 1.0001); medians 2.9e-12 and 0.0. Gate 5e-6.
- The native single-pass history (iteration 0) is NOT the production CAMB history (two further diffusion iterations): gate model D2.

## Benchmarks (`benchmark/chunk5cd_benchmarks.jl`)
`return_solution_to_grid` (3199 rows -> 10000 nodes) 0.78 ms (783 KiB); complete single pass (pass + tail + assembly) 798 ms (904 MiB).

## Files
`src/RecombinationODE.jl`, `src/CosmoRec.jl`, `test/chunk5d_output_native.jl`, `test/runtests.jl`, `test/fixtures/native_ode_pass_grid.txt`, this file, `docs/ACCEPTANCE_CHUNK5D.md`. No Git operation.
