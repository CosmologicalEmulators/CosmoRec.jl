# Chunk 1a acceptance: limited linear stiff-solve probe

The supervisor read all changed probe/test/benchmark/diagnostic code and independently reran the repaired validation on 2026-09-30 EDT. This accepts only the scoped linear analytic problem, not nonlinear recombination physics, arbitrary callback support, or a production CosmoRec.jl implementation.

## Independent evidence

- `julia --project=/home/marcobonici/Desktop/work/CosmologicalEmulators/CosmoRec.jl -e 'using Pkg; Pkg.test()'`: exit 0, **61/61 assertions**, test execution 4m24.4s.
- Fresh-process `benchmark/chunk1a_route_diagnostic.jl throw`: GaussAdjoint and QuadratureAdjoint both reach the deliberate MooncakeVJP sentinel through `adjoint_sensitivity_backpass`; the default sensealg path does not. The sentinel is isolated to that process, not installed package source or package tests.
- `benchmark/chunk1a_benchmarks.jl` in the persistent benchmark environment: exit 0. Median hot primal / hot prepared reverse: Rodas5P with solver-internal AutoFiniteDiff **0.246 ms / 7.778 ms**; QNDF **0.280 ms / 6.546 ms**. Allocation/memory and warmed preparation definitions are in the logs. Cold compilation was not measured.
- Analytic-oracle tests exercise parameter-dependent initial values and rates, early/interior/late observations, tolerance refinement, non-cancelling projections, and independent preparation reuse across changed parameters.

Supervisor logs are retained externally under
`/home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk1a_repair01/`:
`supervisor_repaired_pkgtest.log`, `supervisor_route.log`, `supervisor_benchmark.log`.

## Limits and retained findings

Default Rodas5P's solver-internal AD fails in the backward adjoint solve's time-gradient computation. AutoFiniteDiff works for this linear toy; nonlinear adequacy is not proven. QNDF derivatives converge against the analytic oracle, but finite-difference numerical references can be noise-dominated and are diagnostic-only for that candidate. Test gates no longer accept a 100%-error finite-difference ceiling.

The package still has no recombination implementation or native physics fixtures. Native-original TEXT fixtures and unconditional tests must precede each later translation. The native fixture contract is preliminary: actual capture grids/matching intervals and source provenance must be checked when the capture chunk starts; no native endpoint correction is authorized by this acceptance.

## Next authorized worker scope

Chunk 1b only: a small radiation-shaped, method-of-lines/structured-linear-solve compatibility probe with a fixed nonuniform frequency-like mesh, active time/parameter-dependent coefficients and boundary inputs, independent manufactured reference, and full ForwardDiff/Mooncake solve-stage checks. Keep it test-only and stop for independent supervision before any radiation physics or native-table translation.
