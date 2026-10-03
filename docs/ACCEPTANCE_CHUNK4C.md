# Chunk 4c acceptance: preliminary Recfast++ history (production terms)

## Gates
- `test/chunk4c_recfast_native.jl`: **8475/8475** (focused; log `chunk4c/focused_final.log`): 900 direct native `evaluate_Recfast_System` calls, 36 Saha/`Te_QS` probes, 7 `H_z` probes, the 6000-node grid, the Saha segments and the ODE nodes against the native history.
- `test/chunk4c_recfast_ad.jl`: **29/29**: ForwardDiff through the solve vs finite differences, JVP, 256-bit check of the closed-form Saha node, prepared Mooncake VJPs (independent preparations, changed parameters).
- Full root `Pkg.test()`: **46318/46318**, `Testing CosmoRec tests passed`, exit 0, 20m19.5s, log `chunk4c/pkgtest.log` (37814 previous + 8475 + 29). This count is for 4c only; the 4b result is separate.
- Benchmarks (BenchmarkTools): RHS 429 ns, grid 316 us, full 6000-node history 18.4 ms, ForwardDiff Jacobian 39 ms, Mooncake prepared hot gradient 720 ms (cold prepare separate: first sample includes compilation).

## Accepted scope and measured errors
Production Recfast++ terms only (F = 1.14, A2s1s = 8.2206 as observed natively; no DM/PMF/reionization/CF; factors 1). Native RHS reproduced bitwise for f3, f4 and to 3.6e-16 / 3.8e-16 (relative to the two cancelling terms) for f1, f2 (raw worst 2.0e-6 at equilibrium nodes); `Te_QS`, Saha helpers bitwise; grid 2.6e-15; Saha-segment values bitwise except Xe on segment 2 (2.5e-11, ulp amplification through d + sqrt(d^2 + ...)). The ODE part agrees with the native Gear solver to its own tolerance: at most 0.20 / 1.04 / 48.9 native-tolerance units for Xe_He / Xe_H / TM (Xe 1.04); derivatives compared via the solver's dense output (3.3e-2 of the peak, median 3e-7). The native solver is a reference, not reproduced step by step; Rodas5P results are never written into a fixture.
Gradients: ForwardDiff vs finite differences 1.6e-7 (best step), JVP 3.3e-9, 256-bit check of the Saha node 0.0 (AD formulas) / 1.9e-10 (Float64 conditioning); Mooncake VJP 2.6e-10 (random projection) up to 4.7e-7 (tiny Xe_He block).

## Limitations to carry forward
The Saha-segment breaks and node grid are discrete (a 1e-12 relative parameter change moved the ODE start node by one); no derivative of grid/solver steps is claimed. The reverse route needs an explicit Jacobian/time gradient and a plain-function RHS; QNDF/FBDF adjoints were rejected (wrong TM gradient). No 256-bit ODE solve. The low-z native Hubble bug and the `A2s1s = 8.2206` oddity are reproduced transparently, not fixed. Ill-conditioned f1/f2 at equilibrium nodes make native `dXe` (solver derivative) not reproducible from the RHS.

## Still not a full CosmoRec model
No sampled-HeI switch (4d, separate), no CosmoRec ODE/PDE, no diffusion/feedback, no history or CAMB spectrum parity.
