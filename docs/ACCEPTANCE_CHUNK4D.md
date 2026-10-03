# Chunk 4d acceptance: sampled-HeI state switch (discrete reset and packing)

## Gates
- `test/chunk4d_helium_switch_native.jl`: **984/984** (focused): 80 native switch points (40 switched), decision, post-state, `HeI_Atoms.Xi`, pre/post packed vectors and lengths bitwise equal; decision logic, `flag_He`, dimension errors.
- `test/chunk4d_helium_switch_ad.jl`: **32/32**: frozen-branch reset/pack operator, ForwardDiff vs 256-bit central differences (0.0), prepared Mooncake VJPs (independent preparations, changed inputs, 0.0), explicit no-event-sensitivity tests.
- Full root `Pkg.test()`: **47334/47334**, `Testing CosmoRec tests passed`, exit 0, 20m53.9s, log `chunk4d/pkgtest.log` (46318 previous + 984 + 32).
- Benchmarks (BenchmarkTools): decision 14 ns, decision + reset + packing 33 ns, ForwardDiff Jacobian 0.84 us, Mooncake prepared hot gradient 0.69 us (cold prepare separate).

## Accepted scope
Native decision `(zs < 200 || fHe - Xi(0) <= 1e-7) && flag_He`, reset (He 1s = fHe, other helium slots 1e-300, `Xe`/H/`rho` unchanged) and the 7-entry `flag_He = 0` packing, as a discrete operator on explicit state with a frozen branch decided on primal values. The reset depends on `fHe`. The harness executes the inline native statements on the original objects (the driver code is not callable) and then the original `copy_LI_to_ysol`; this provenance is stated in the fixture header.

## Limitations to carry forward
No ODE, no solver restart, no event location in z, and NO continuous event sensitivity (the decision is piecewise constant and the state dimension changes); derivatives hold only while the primal branch is unchanged. Only the production atoms and the 1 -> 0 transition are covered.

## Still not a full CosmoRec model
Phase 5 (the SciML recombination ODE of the full model with the 4a/4b/4c/4d pieces), diffusion/feedback, iteration loops, and history/spectrum parity are untouched.
