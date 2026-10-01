# Chunk 3e results: native `DPesc_coh` fallback of the H-I absorber

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al. Use must be acknowledged; cite Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39), Rubino-Martin, Chluba & Sunyaev 2008. Data tables are not bundled.

## Scope
Only the explicit coherent-scattering integral used when `DP_interpol_S/T` leave the pre-tabulated DP table, plus its plumbing into `hi_absorption_rhs!` / `fcn_effective!`. No ODE, no background, no helium feedback/diffusion.
Route: pure Julia (`src/DPescCoh.jl`); no GSL; the SciML/AD route is the ordinary Julia code path (ForwardDiff, Mooncake). No claim about the complete CosmoRec model.

## Native evidence (hard gate)
- External read-only harness (outside any git repo): `.../chunk3e/native_capture/harness3e.cpp`, linked to the original `libCosmoRec.a`, `libRecfast++.a`, GSL; startup through the original chain; the ODE is never run. It calls the ORIGINAL `call_DP_Singlet/Triplet`, `DPesc_appr_I_sym`, `DP_interpol_S/T`, and the original `Voigtprofile_Dawson` / `Photoionization_cross_section_SH` objects.
- 198 queries (9 redshift states x {S,T} x 11 perturbations of tau_S, eta_c, T; scalar inputs printed with each result). At z = 1500, 1200, 800 the native `DP_interpol_S/T` returned `call_DP_*` bit-for-bit for 64 queries (10 of 11 at z=1500; the other states are inside the table). Three fresh-process runs byte-identical.
- Fixture `test/fixtures/native_dpesc_coh.txt` (plain values, SHA-256 in the header and in `test/chunk3e_helpers.jl`), containing hashes of harness, binary, ini, libraries and read sources.
- Caveat: unlike the Chunk 3d harness, the config switches (Diffusion etc.) are not asserted at runtime here; the DPesc math does not depend on them, and the end-to-end check below uses the Chunk 3d fixture states.

Exact algorithm: `DPesc_coh` -> `P=p_ij(pd tau)+Dpij`, `Pesc=pd P/(1-(1-pd)P)`, returns `Pesc-p_ij(tau)` (no fcorr in the fallback). `Dpij=DPesc_appr_I_sym`: `epsabs = P_d*1e-5` first, then `|r| 1e-6`, seven intervals [0,1,2,4,10,30,100,1e4] of `Inner(x)+Inner(-x)`; analytic inner integral; Voigt (Mihalas expansion, Dawson, wing switch |x|=30); `Integrate_using_Patterson_adaptive` (1,3,...,255 nested rule, stop `|dr|<=max(epsabs,|r|1e-6)`, `DC_sumprod` pairwise summation; its refine branch is dead since the rule routine always returns 0).

## Julia implementation
- Patterson nodes/weights are generated from first principles (BigFloat Stieltjes extension + interpolatory weights), not copied. Offline read-only comparison with the native tables: node differences <= 5.6e-17, relative weight differences <= 5.1e-16. Tests check exactness degree, interleaving and weight sum.
- Dawson/erf rational coefficients are mechanically transcribed from `routines.cpp`; atomic line data come from the fixture (`nu21, Gamma, AM` of `nP_S/T_profile(2)`); H 1s cross section shape uses the closed-form Gaunt factor (equals the native evaluation to ~1e-15, tested at 15 frequencies).
- Fallback semantics inside `dp_correction`: T outside the table returns the fallback directly (unscaled Tg, abs values, no fcorr); eta outside a sheet -> the fallback is that sheet's DP; tau outside -> every eta row = fallback then the four eta weights are applied (as native). Opt in with `fcn_effective(...; dp_fallback = tr -> dpesc_fallback(NATIVE_DPESC, tr))`; without it `DPTableDomainError` is still thrown (no silent fallback).

## Parity (all captured queries; focused test 1512 tests)
| quantity | max rel. error |
|---|---|
| singlet Dpij / correction | 2.6e-14 / 1.9e-15 |
| triplet Dpij / correction | 1.7e-10 / 8.8e-10 |
| Voigt phi / xi_Int | 2.2e-16 / 2.1e-16 |
| `fcn_effective` z=1500,1200,800, absorber on, all 15 components | 4.6e-16 |
Triplet error is conditioning: Dpij ~1e-7 is a cancellation of O(1) terms; at the worst row native is 4.0e-11 from the 256-bit evaluation of the same formula, Float64 Julia 2.1e-10. Tolerances in the tests: S 1e-12, T 1e-9 (Dpij) / 5e-9 (correction), fcn 1e-13. The (on - off) absorber part also matches native (rtol 1e-9).

Integration refinement (z1500 singlet, forced Patterson level on all 7 intervals vs native Dpij): level 5 (31 pts) 1.5e-13, level 6 6.0e-10, level 7 and 8 5.96e-10. The native value therefore differs from the fully converged integral by 6e-10 (its stop criterion), and the Julia port reproduces the native stopped value, not the converged one (as required). Orders used at z=1200: (4,4,4,4,4,3,5) = (15,15,15,15,15,7,31 points).

## AD (rules pinned to the primal-selected orders; 281 tests)
x = (tau_S, eta_c, T, pd), singlet / triplet / combined, at z=1500,1200,800.
- BigFloat(256) ForwardDiff vs 256-bit central difference: 0 (exact AD formulas).
- Float64 ForwardDiff vs 256-bit FD: singlet <= 1.9e-14; triplet and combined <= 3.3e-9 (same conditioning as above).
- Directional derivative: <= 2.7e-15. Prepared Mooncake VJP (two independent preparations, reused at 3 changed parameter points, 3 seeds each) vs ForwardDiff J'w: <= 1.2e-14.
- Adaptive (primal-selected) and pinned-order ForwardDiff Jacobians are identical.
- Limitation: differentiation is through the fixed final rule; across a Patterson order switch the map is not differentiable (the Jacobian changes by 3.8e-10 relative when sub-interval 1 is lowered one level at z1200), nor across the |x|=30 wing switch, sigma(nu) edge/cut branches, or DP table seams. No global differentiability claim. No `DPTableDomainError` counts as gradient success.

## Tests / runtime
- Focused: `test/chunk3e_dpesc_native.jl` (1512 pass), `test/chunk3e_dpesc_ad.jl` (281 pass).
- One full root `Pkg.test()`: 22656/22656 pass (18m14s; log `.../chunk3e/pkgtest.log`). An earlier accidental full run, launched with its output discarded, was killed by me before completion and is not counted as a pass.
- Benchmark (`benchmark/chunk3e_benchmarks.jl`, evals=1, z1200 singlet): primal integral 27 us (2.9 KiB), full `dpesc_coh` 28 us, ForwardDiff Jacobian (1x4) 52 us, Mooncake prepare (cold) 1.5 ms, prepared hot gradient 333 us. Patterson rule generation happens at precompilation (BigFloat) and is paid once.

## Failed approaches / notes
- Initial comparison scripts used wrong fixture columns (reported ~100% errors); fixed, not a code defect.
- Copying native node tables was avoided on purpose.
- Not done: ODE, helium feedback/diffusion, background, `epsabs`-induced order switches under AD.
