# Chunk 5c results: the Recfast tail (z = 50 -> 0.001)

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al.; cite Chluba & Thomas 2010, Chluba & Sunyaev 2006, Rubino-Martin et al. 2008. No native source, library or table was changed or bundled; the tests require `COSMOREC_NATIVE_DATA_DIR` (D1).

## Scope
`recfast_tail`, `recfast_tail_inputs`, `init_xarr_log`, `recfast_tail_rhs!` (`src/RecombinationODE.jl`) port `compute_Recfast_part` (`CosmoRec.cpp:278-324`), `Cosmos::recombine_using_Recfast_system` (`Cosmos.cpp:795-840`) and `Xe_frac_rescaled` (`Recfast++.cpp:336-441`) on top of the 4c Recfast RHS: inputs from the last ODE node (`Xe_Hi = 1 - X1s`, `Xe_Hei = max(1e-30, 0)`, `Xei = Xe = xp`, `TMi = TCMB rho`, `dXei = -dX1s/dz`), 200 logarithmic nodes `init_xarr_log(zi, 0.001, 200)` (accumulated like the C++), rescaling factor `ff = dXei / f_xp(zi, y0)` applied to the `xp` equation (native `fcn_rescaled`), `Xe = Xe_H + Xe_He`. The Recfast system uses the loaded Cosmos H(z) (4b), hence the reproduced native low-z Hubble bug (z < 1: +32% at z = 0.01) acts on the tail exactly as natively. The stiff solve is the caller's (tests: Rodas5P, analytic Jacobian, 1e-12).

## Native oracle
The tail rows of the Chunk 5b native pass (`output_CosmoRec` rows 3001..3199, fixture `native_ode_pass_output.txt`) and the native inputs of `compute_Recfast_part` (fixture `native_ode_pass_summary.txt`: `RECFAST_INPUT` = `ze`, `Xe_Hi`, `Xei`, `rhoi`, `dXe_dz` computed by the harness through the ORIGINAL `fcn_effective(ze, Level_I)` with `flag_He = 0` at the last native state, and the final `Level_I.X`). The native Recfast tail solver (rel 1e-8, abs (1e-10, 1e-8, 1e-10)) is the error scale.

## Results (focused 17/17)
- Tail grid: bitwise equal to the native rows (0.0); last node 0.001.
- Native inputs from the native final state through the Julia mapping: `Xe_Hi`, `Xei`, `TMi` to 1e-14; `dXei` differs by 2.1e-4 relative: `dX1s/dt` at z = 50 is a difference of large rates (1.8e-20 net), conditioning-limited (gate 1e-3). Its downstream effect is negligible (end-to-end below).
- Tail from the NATIVE inputs vs the native rows (isolates the solver): Xe 1.03 native-tolerance units (relative 5.9e-5, dominated by the native absolute tolerance 1e-8 at Xe ~ 2e-4 near z = 0.001), Te 166 units (relative 1.8e-6 at z = 2.38, a global error accumulated by the native solver; gate 500 units). Rescaling factor ff = 1.0152648; `ff f_xp(zi, y0)` reproduces the native derivative to 1e-12.
- End to end (tail from the Julia pass): Xe 1.4e-7, Te 3.0e-8 relative to the native rows.
- Limitations: native Recfast tail error is much larger than the Julia solver error; the `dXei` input is cancellation-limited; one cosmology.

## Benchmarks (`benchmark/chunk5cd_benchmarks.jl`, `chunk5/bench5cd.log`)
`recfast_tail` (200 nodes, Rodas5P 1e-12) median 2.76 ms (647 KiB; first sample includes compilation).

## Files
`src/RecombinationODE.jl`, `src/CosmoRec.jl`, `test/chunk5_helpers.jl`, `test/chunk5c_recfast_tail_native.jl`, `test/runtests.jl`, `benchmark/chunk5cd_benchmarks.jl`, this file, `docs/ACCEPTANCE_CHUNK5C.md`. No Git operation.
