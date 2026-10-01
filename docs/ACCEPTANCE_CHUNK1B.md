# Chunk 1b acceptance: small linear structured PDE probe

The supervisor inspected the helper, nonuniform stencil construction, sparse Jacobian assembly, numerical/AD tests, diagnostics and benchmarks, then independently reran validation on 2026-09-30 EDT / 2026-10-01 UTC.

## Independent evidence

- Full `Pkg.test()`: **495/495 assertions**, exit 0, test execution 14m57.4s. Includes all 61 Chunk1a and 434 Chunk1b tests unconditionally.
- Fresh-process structured route diagnostic: ForwardDiff uses a 7x7 `SparseMatrixCSC{Dual}` Sparspak factorization. GaussAdjoint and QuadratureAdjoint each reach the MooncakeVJP hook during preparation and evaluation, and use 7x7 Float64 sparse factorizations. No installed package source is edited; process-local instrumentation is outside package tests.
- Independent bounded BenchmarkTools run: exit 0. Mesh_a (7 unknowns), median primal3.799ms / warmed prep8.246ms / hot prepared reverse98.142ms. Mesh_b (10 unknowns), primal6.519ms / prep11.156ms / reverse155.143ms. Allocation counts reproduce the worker measurements. These are unoptimized small-probe timings, not production predictions.
- Test coverage separately checks polynomial derivative stencils and boundary lifting, coefficient/source/initial/boundary maps, operator Jacobian, complete solve derivatives against manufactured sensitivities, tolerance refinement and independent cache reuse.

Independent logs are in the external analysis directory `chunk1b_resume01/`: `supervisor_pkgtest.log`, `supervisor_route_count.log`, `supervisor_benchmarks.log`.

## Accepted scope and limitations

Accepted route: `Rodas5P(autodiff=AutoFiniteDiff(), linsolve=SparspakFactorization())`, analytic sparse primal Jacobian, continuous SciMLSensitivity adjoint using MooncakeVJP, on the two tiny NONUNIFORM manufactured linear-in-state PDE problems. The polynomial source is derived from independent continuum expressions, while map/operator tests separately check dependencies that cancel on the manufactured trajectory.

This is not the CosmoRec radiation PDE, not the native characteristic-traced boundary, not a frequency-integral correction test, and not a nonlinear population solve. Mesh scaling and nonlinear adequacy remain unproven. QNDF reverse remains a retained LIMITED diagnostic, not an accepted accuracy path. Quadrature is checked at the reported default tolerance point, not claimed to have the same full refinement study as Gauss.

The first four focused assertion failures and the walltime-interrupted package run remain in logs. The boundary-support assertions were corrected to the actual five-point stencil; the default-tolerance Quadrature ceiling was recalibrated, while the Gauss refinement gate remains tight. Passing a gate does not establish a general solver-tolerance bound on all possible gradients.

## Next scope

A small Chunk1c nonlinear population-shaped compatibility check is warranted before native population physics: an analytically soluble quadratic loss RHS with active rate and initial-state dependencies, verified through ForwardDiff and prepared Mooncake. This is a targeted extension of the existing compatibility probe, not permission for a new framework, a nonlinear radiation model, or a large third demonstration. Stop for review afterward. Original-code text fixtures precede real physics translation.
