# Chunk 4b results: explicit Cosmos accessors and GSL-compatible splines

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al. Use must be acknowledged; cite Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39), Rubino-Martin, Chluba & Sunyaev 2008. No native source or data was changed.

## Scope (and nothing more)
`src/CosmosAccessors.jl` (included from `src/CosmoRec.jl`, exports added):
- `NaturalCubicSpline`/`natural_cubic_spline`: exact coefficient algorithm of `gsl_interp_cspline` (natural, tridiagonal LDL solve, `b,c,d` per interval); `spline_eval` (bisection interval, `SplineDomainError` outside the node range = native GSL error handler abort, no extrapolation); `spline_eval_native` (native `calc_spline_JC` 1e-14 relative end nudge).
- `CosmosConstants`, `RecfastSplines`, `recfast_splines`: native `Cosmos::calc_coeff_X_spline` from the 6000-node history in NATIVE order (z descending): splines of `ln Xe`, `dXe`, `ln Xe_H`, `dXe_H`, `rho = TM/TCMB`, `ln max(Xe_He, 1e-20)` on ascending z. `Xe_H = n_e,H/n_H` and `Xe_He = n_e,He/n_H` are IONIZED populations, not neutral fractions; the native accessors derive `X1s = 1 - exp(spline)` and `XHeI1s = fHe - exp(spline)` from them.
- `hubble_table`/`cosmos_H`: native `Hubble::init` (spline of `ln H` vs `ln z`, valid for `zmin < z < zmax`) with the native analytic form outside. The native bug (node 0 keeps `ln H = 0` because `lHz` is filled from `k = 1`) is REPRODUCED, not repaired; the effect (H wrong for z <~ 1) is confirmed against the native outputs.
- Accessors: `cosmos_TCMB, Nb, NH, fHe, sigT, rho_g_1_cm3, Xe_Seager, Xe_b, dXe_dz, Ne, Ntot, kappa_cool, Te_Tg, Te, X1s, Xp, XHeI1s, XHeII1s, NHeI, NHeII, NHeIII`, `SahaBoltz_HI1s/HeI1s/HeII1s/HeII/HeIII` (including the `Ttoz` round trip of the native code), the `z >= 3500` (`z_saha`) and `z >= zsRe` branches, and `saha_inputs_at` (feeds the 4a initializer).
NOT implemented / not claimed: Recfast++ itself (Stage 4c), the cosmology closure (`Omega_L`, `Nb0`, `rho_g_gr`, `H0` are explicit inputs), any ODE, the sampled-HeI switch.

## Native fixtures (text, plain values, hashes tested)
- `test/fixtures/native_recfast_history.txt` (749574 B, SHA-256 `f893af2d8a1f52c0278e422a5490616acb347081a10c09981de0b33020e97716`): the 6000 native history nodes (z descending 25000 -> 0; columns z, Xe_H, Xe_He, Xe, dXe/dz, dXe_H/dz, TM[K]), verbatim from the corrected `chunk4/native_capture/out4b/history1.txt` (byte-identical to `out/history1.txt`).
- `test/fixtures/native_hubble_input.txt` (418321 B, SHA-256 `bcfa199b78d756f433202f0997bf4315761b28e6193b363745fdbc25c1f2ed99`): the CAMB 2.0.4 H(z) table passed to the native library (an INPUT, not a native output).
- `test/fixtures/native_cosmos_accessors.txt` (388645 B, SHA-256 `b34bffd9234a2fa334533ce5defd3af06239e0a56c97587374fce39e44622ee9`): COS2 constants, 658 ACC records (21 accessors each, including the exact history nodes and node mid-points, the ends, 3500, zsRe, low z), 6 SBT records. Produced by the external read-only harness `harness4b.cpp` (hash in the fixture header, `chunk4/native_capture/hashes.txt`) on the production batch route; two fresh-process outputs byte-identical. No history file beyond these is bundled.

## Results
Focused native suite `test/chunk4b_cosmos_native.jl`: **14097/14097**. Maximum relative errors against the native values over 658 redshifts: H 1.8e-16, TCMB/Te/Te_Tg/sigT 0 (bitwise), Nb 4.2e-16, NH 4.2e-16, Xe_Seager 1.6e-16, Xe_b 1.8e-16, Xp 1.1e-16, X1s 8.1e-15, XHeI1s 8.2e-16, XHeII1s 7.0e-16, NHeI 1.0e-15, NHeII 1.0e-15, dXe/dz 1.6e-14, kappa_cool 4.3e-16, Ne 4.3e-16, Ntot 4.5e-16, rho_g 2.1e-16, all five Saha-Boltzmann functions 0 (bitwise); NHeIII 6.0e-7 is the cancellation `(fHe - XHeI1s - XHeII1s) NH` at z >= 3500 (same operations as native; tolerance 1e-5). 4a fed by 4b (`saha_inputs_at` -> `saha_initial_state`) reproduces the native 4a vectors at 6 redshifts to 4.6e-16. Spline reproduction of the nodes, the independently checked z grid (strictly decreasing, 25000 -> 0), the z_saha switch sides, the zsRe closed forms, the endpoint nudge, domain errors and invalid constructors are tested.
Focused AD suite `test/chunk4b_cosmos_ad.jl`: **16/16**: ForwardDiff z-derivatives of the accessor/Saha-input map vs 256-bit central differences at z = 300, 800, 1500, 2999.5, 4000, 8000: 1.3e-11 (z = 300), 2.6e-13, 4.9e-16, 4.5e-16, 9.4e-16, 4.8e-19; directional derivative over all 24001 inputs (history node values and z) vs 256-bit central difference 3.8e-15; prepared Mooncake VJP (two independent preparations, changed inputs, 24001 inputs) vs the ForwardDiff gradient 4.2e-16, `gA == gB`.
Full root `Pkg.test()` (`chunk4b/pkgtest.log`, 16:59:11 to 17:21:32 EDT): **37814/37814**, `Testing CosmoRec tests passed`; 37814 = 23701 (through 4a) + 14097 + 16. The supervisor independently reran the focused native (14097/14097) and AD (16/16) suites.

## Limitations
- Piecewise definitions: the map is C2 across spline nodes (derivative continuity across a node 1.6e-11 relative), but NOT differentiable across the `z_saha = 3500` switch (X1s jumps ~1e-3 relative, its z-derivative by 1.7e-4), `zsRe`, the end nudge, the loaded-H edges `zmin/zmax` or the `Xe_He` floor. Outside the spline range every accessor throws.
- The native low-z Hubble bug is reproduced (H(1e-3) is ~10x the table, +32% at z = 0.01, -29% at 0.5, <1e-4 above z ~ 1); it does not affect z >= 100 (<1.4e-8) and is not a Julia defect.
- One cosmology and one table; `Omega_L`, `Nb0`, `rho_g_gr` are supplied (native-printed), not derived.
- AD tests use interior z only; z = 300 shows the largest finite-difference gap (1.3e-11) from exp-conditioned Saha terms.

## Benchmarks (`benchmark/chunk4b_benchmarks.jl`, `chunk4b/bench.log`, BenchmarkTools)
| path | median | memory |
|---|---|---|
| six splines from 6000 nodes | 592 us | 3.58 MiB, 234 allocs |
| Hubble table (10000 nodes) | 206 us | 938 KiB, 36 allocs |
| `cosmos_Xe_Seager` | 43 ns | 0 B |
| `cosmos_H` (loaded table) | 61 ns | 0 B |
| `saha_inputs_at` | 163 ns | 0 B |
| inputs + initialize + pack | 392 ns | 336 B, 4 allocs |
| ForwardDiff gradient, 24001 inputs (1 sample, includes compile) | 17.7 s | 45.6 GiB, 72% GC |
| Mooncake `prepare_gradient` (3 samples; first includes compilation) | 50.8 ms | 113.7 MiB |
| Mooncake prepared hot gradient | 10.2 ms | 8.71 MiB |

## Files
`src/CosmosAccessors.jl`, `src/CosmoRec.jl`, `test/chunk4b_helpers.jl`, `test/chunk4b_cosmos_native.jl`, `test/chunk4b_cosmos_ad.jl`, `test/runtests.jl`, the three fixtures above, `benchmark/chunk4b_benchmarks.jl`, this file, `docs/ACCEPTANCE_CHUNK4B.md`. No Git operation was performed.
