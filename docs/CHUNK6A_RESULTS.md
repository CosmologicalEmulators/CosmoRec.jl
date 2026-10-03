# Chunk 6a results: HI diffusion-PDE coefficients (populations, pd, Dnem)

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al. Use must be acknowledged; cite Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39). No native source or data was changed.

## Scope
`src/HIPDECoefficients.jl`: `load_hi_bitot_table` (the file stores ln Bitot), `get_rates_all` (A, B, Bitot, all R), `HIPopulationSplines` (native `compute_Xi_HI_splines`: last stored row dropped; ln Xe, ln Xi, rho-1), `hi_rp_rm_pd` (`calc_HI_Rp_Rm_pd_all`), `hi_pde_coefficients` (native `set_up_splines_for_HI_pd_Rp_effective` on `init_xarr(ze/1.0001, zs*1.0001, nz)`, `nz` = number of stored rows with z >= ze/1.0001, the native trailing-zero behaviour kept), `hi_pd`, `hi_Dnem`. Production range zs = 2500, ze = 500.

## Native oracle
`chunk6/native_capture/harness6a.cpp` after one ORIGINAL runmode-1 pass. `cosmorec_calc_h_cpp_` clears the loaded CAMB H on return (CosmoRec.cpp:1029), but production evaluates the coefficients INSIDE the call with the loaded H. The harness therefore re-arms `cosmos.init_Hubble` with the same table before the coefficients and calls `clear_Hubble` before it returns. The latter fixed a "corrupted double-linked list" abort at exit. Two fresh-process outputs are byte-identical.
REJECTED evidence, preserved and not used to choose tolerances: the first capture evaluated the coefficients with the analytic H after `clear_Hubble`; Julia's background then disagreed (bg:H 2.96e-4). It is kept in `chunk6/native_capture/out/failed_run{1,2}_cleared_hubble` and `chunk6/failed_fixture_cleared_hubble_native_hi_pde_coefficients.txt`.
Fixture: `test/fixtures/native_hi_pde_coefficients.txt` (sha256 98aa7316...b855, tested).

## Methodology for Dnem
Dnem is a DIFFERENCE between a population-ratio term and exp(-x), so near equilibrium it cancels. Errors are measured against the source-derived residual scale: `Ni/N1s/(2l+1) + e^{-x}` for s and d states, and `(nL + e^{-x}) |1 + (1/pd - 1) PS|` for np. The raw relative and absolute errors are reported as well. Before choosing any tolerance, two genuine defects were fixed: a double log of Bitot, and the coefficient-grid size (native nz counts the rows above the range too).

## Results (full suite run, `chunk6/pkgtest6a.log`)
- Native populations as input: populations 2.1e-16, get_rates_all 2.2e-16, Rp/Rm and pd 1.8e-16, background 0 (bitwise), pd spline 3.4e-16. Dnem spline: scaled 2.2e-16, raw relative 9.7e-9, absolute 1.1e-23.
- End-to-end from the Julia 5b pass: pd 3.4e-16. Dnem: scaled 2.25e-7, raw relative 0.186. The raw number is a cancellation artefact. The scaled number tracks the excited-population difference of the Julia pass, 2.23e-7, i.e. the native solver tolerance. The focused run gave 2.18e-7 / 2.17e-7 (raw 0.244). That is a run-to-run solver-level variation of the Julia stiff solve; the gate scaled < 10 x pop_excited holds in both runs.
- AD (`test/chunk6a_pde_coefficients_ad.jl`, 24000 population inputs):
  - The ForwardDiff directional derivative agrees with a 256-bit central difference to 5.0e-15.
  - The prepared Mooncake VJP agrees with the ForwardDiff gradient to 7.0e-16, and two independent preparations give identical results.
- Full root `Pkg.test()`: **50275/50275**, exit 0, 32m33.1s (= 49444 + 825 + 6).

## Benchmarks (`benchmark/chunk6a_benchmarks.jl`, `chunk6/bench6a.log`)
| path | median |
|---|---|
| population splines (3000 rows) | 811 us |
| get_rates_all | 624 ns |
| pd/Dnem coefficient splines | 4.45 ms |
| hi_Dnem evaluation | 24 ns |
| ForwardDiff directional (24000 inputs) | 7.48 ms |
| Mooncake prepare_gradient (cold, includes compilation) | 1.01 s |
| Mooncake prepared hot gradient | 54.8 ms |

## Limitations
Derivatives are not claimed across spline cells (C2 only), the discrete coefficient-grid size, rate-table stencils or the range ends. One cosmology. The Dnem raw relative error is not meaningful near equilibrium (see methodology).
