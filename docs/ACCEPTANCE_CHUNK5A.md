# Chunk 5a acceptance: ODE right-hand side

## Gates
- `test/chunk5a_ode_rhs_native.jl`: **1966/1966** (focused); `test/chunk5a_ode_rhs_ad.jl`: **35/35**.
- Full root `Pkg.test()` with `COSMOREC_NATIVE_DATA_DIR` set: **49335/49335**, `Testing CosmoRec tests passed`, exit 0, 22m06.8s, log `chunk5/pkgtest5a.log` (47334 previous + 1966 + 35).
- Errors: dy/dz 9.5e-15 (`flag_He = 1`), 2.5e-16 (`flag_He = 0`); AD: 256-bit 0.0, Float64 3.6e-16, VJP 2.4e-16. Benchmarks: RHS 51.8 us / 0.84 us, ForwardDiff Jacobian 56 us, Mooncake prepared gradient 0.76 ms.
## Accepted scope
Packed 12-/7-state `dy/dz` from the 3d RHS with the explicit full production tables (read from the user-specified native directory, never bundled; tests fail loudly if unset), including the new `flag_He = 0` branch and the DP fallback. Native Jacobian not ported.
## Limitations
Piecewise table/DP/absorber branches (interior points only); explicit-state (not trajectory) 7-state oracle; machine-specific data path by decision D1.
## Still not a full CosmoRec model
No diffusion/feedback (phase 6+), no CAMB-level (runmode 0) parity.
