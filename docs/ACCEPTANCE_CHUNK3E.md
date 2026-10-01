# Chunk 3e acceptance: native coherent-DP absorber fallback

## Independent gates

The supervisor reran both focused suite entry points in the instantiated local Julia test environment:

- `test/chunk3e_dpesc_native.jl`: **1512/1512** assertions; direct original native `call_DP_Singlet/Triplet` parity, integral refinement levels, Voigt/x_i intermediate values, Patterson-rule checks and end-to-end `fcn_effective` parity for the absorber-on states at z=1500/1200/800.
- `test/chunk3e_dpesc_ad.jl`: **281/281**; 256-bit ForwardDiff vs 256-bit finite difference, float64-conditioned deviations, JVP, and prepared Mooncake VJPs against ForwardDiff. Worst listed VJP discrepancy `1.2e-14`; no gradients claimed across discrete adaptive-order or wing/domain switches.
- Worker root `Pkg.test()` **22,656/22,656**, exit0, log `chunk3e/pkgtest.log`. The full suite includes 20,863 from prior chunks + 1512 + 281; supervisor independently verified both additions but did not rerun the full root suite.

## Accepted scope

`src/DPescCoh.jl` ports the original native `DPesc_coh` coherent-scattering integral used when native `DP_interpol_S/T` data falls outside its stored temperature/eta/tau grid. It ports the actual nested Patterson integration, analytic inner integral, Mihalas Voigt/Dawson representation and native stopping conditions/order selection; it does not simply extrapolate the DP table. Fixture `test/fixtures/native_dpesc_coh.txt` captures 198 direct native queries (nine states, singlet/triplet, input perturbations) and provenance/hash data. Harness/objects are external; original CosmoRec source, static library and user native data were not modified.

The fallback is opt-in in `fcn_effective(...; dp_fallback=...)`; no callback still throws `DPTableDomainError` instead of silently switching physics. Table-call plumbing reproduces the native fallback behavior for out-of-table T, eta sheets, and tau rows. End-to-end composition with the saved native `fcn_effective` vectors at previously failing z=1500/1200/800 now matches with maximum relative error `4.6e-16`.

## Numerical/AD limitations

- Singlet `Dpij` max relative error `2.6e-14`; triplet `1.7e-10`, and triplet final DP correction max `8.8e-10`. The latter is cancellation of O(1) terms leaving ~1e-7. The native adaptive integration stops at its own specified tolerance; Julia intentionally matches that native stopped value instead of replacing it with a more-converged answer. Refinement/error details are in `docs/CHUNK3E_RESULTS.md`.
- AD differentiates the selected fixed final quadrature rule. Derivatives may jump at Patterson-order-selection, |x|=30 wing, cross-section, or DP-table branch seams. Do not infer global differentiability; future solver work must test the actual path/state distribution and preserve the exact fallback/table switch semantics.

## Still not a full CosmoRec model

No initial conditions/Recfast history/Saha state generator, no helium state switch, no cosmology/H(z) module, no SciML recombination ODE, no diffusion/feedback PDE or 3ODE/2PDE correction iterations, no full trajectory or CAMB-spectrum validation. `DPescCoh` closes one default absorber-domain hole and makes previous native fixtures executable; it does **not** prove the entire trajectory or gradients thereof. Continue with native initialization/switch/background fixtures and the approved SciML solver pipeline. No commit.
