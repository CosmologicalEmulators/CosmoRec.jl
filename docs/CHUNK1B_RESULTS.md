# Chunk 1b results: radiation-shaped method-of-lines structured-solve probe

Toy only (linear in u). No CosmoRec physics, no native fixtures (contract: `docs/fixture_contract.md`).
Checkout `/home/marcobonici/Desktop/work/CosmologicalEmulators/CosmoRec.jl`, branch `develop`, HEAD `80e3c6a`, nothing staged/committed. Julia 1.12.6.
Versions (`benchmark` env; `chunk1b_resume01/versions.txt`): SciMLSensitivity 7.119.12, SciMLBase 3.57.0, OrdinaryDiffEqRosenbrock 2.7.5, OrdinaryDiffEqBDF 2.4.12,
Mooncake 0.5.61, DifferentiationInterface 0.7.21, ForwardDiff 1.4.6, LinearSolve 5.18.2 (new), Sparspak 0.3.15 (new), BenchmarkTools 1.8.0.

## Files (Chunk 1b)
test/chunk1b_helpers.jl, test/chunk1b_radiation_shaped_probe.jl, test/runtests.jl (adds the include; 1a include kept), Project.toml (extras/targets: + LinearSolve, Sparspak, SparseArrays; still no `[deps]`),
benchmark/{chunk1b_benchmarks.jl, chunk1b_route_diagnostic.jl, chunk1b_sweep_diagnostic.jl, chunk1b_candidate_probe.jl}, benchmark/Project.toml (+LinearSolve, Sparspak, SparseArrays), this file. `src/` still an empty stub.
Logs: `cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk1b/` (focused_test_run1.log = FIRST FAILED run, kept; focused_test_run2.log, route_diagnostic_{count,throw}.log, sweep_diagnostic.log, candidate_probe.log) and `.../chunk1b_resume01/` (pkg_test_full.log, benchmarks.log, versions.txt). The older `chunk1b/pkg_test.log` is the walltime-interrupted run (SIGTERM, not a failure).

## Problem
`u_s = kappa(s) u_xx + v(s) u_x - gamma(s) u + src(s,x)`, theta=[k0,v0,g0,A,B,c]=[2,3,1,1,.5,.4].
kappa=k0(1+0.1 sin s), v=v0(1+0.5 s), gamma=g0(1+0.2 cos s) (time dependent, positive). r(theta)=c+0.3k0+0.2v0+0.1g0.
Manufactured solution (cubic in x): `u = A(1+x+x^2) e^{-r s} + B x^3 e^{-s}`; Dirichlet `u(0,s)=A e^{-rs}`, `u(1,s)=3A e^{-rs}+B e^{-s}`; `u0 = A(1+x+x^2)+B x^3`.
Source `src = u_s - kappa u_xx - v u_x + gamma u` from hand-written CONTINUUM derivatives (`ux_exact, uxx_exact, us_exact`), never from the matrices.
Meshes: mesh_a 9 points nonuniform `[0,.07,.17,.3,.45,.62,.8,.9,1]` (7 interior unknowns); mesh_b 12 points, power-law clustered (10 interior unknowns). Observed at s=0.05, 0.3, 1.0 (saveat; saveat points).
Discretisation: Fornberg weights on five-node stencils (centered inside, clamped one-sided at the ends), exact to degree 4; Dirichlet values enter via lifting vectors (boundary node appears only in rows 1-2 and last two rows). Interior blocks A1, A2 and the Jacobian share a fixed union sparse pattern (bandwidth 3).

## Structured route (REAL sparse solve, no densification)
ODEFunction with analytic `jac` writing `J.nzval .= k.*nz2 .+ v.*nz1 .- g.*nzI` on a `SparseMatrixCSC` `jac_prototype`; solver `Rodas5P(autodiff=AutoFiniteDiff(), linsolve=SparspakFactorization())`.
Test C2 inspects an initialised integrator: `f.jac_prototype`, `integ.cache.W`, `linsolve.A` are `SparseMatrixCSC{T}` for T=Float64 and for `ForwardDiff.Dual`; `linsolve.alg isa SparspakFactorization`; cacheval type is Sparspak `SparseSolver`; nnz < size^2/1.4.
Route diagnostics (`route_diagnostic_count.log`): ForwardDiff solve builds Sparspak caches for `SparseMatrixCSC{Dual{...,6}}` (7x7); Mooncake runs build them for `SparseMatrixCSC{Float64}` 7x7. Dual support is therefore exercised through assembly AND factorization. Dense W + Sparspak fails (`getcolptr(::Matrix)`): a misconfiguration, not used; default linsolve with sparse input picks UMFPACK/KLU-type routes (works, diagnostic only).
MooncakeVJP route (process-local sentinel, `route_diagnostic_throw.log`): Gauss and Quadrature reach `get_paramjac_config(::MooncakeLoaded, ::MooncakeVJP,...)` via `adjoint_sensitivity_backpass`; count run: hook 2 calls in prepare, 4 after one gradient, both adjoints. No overrides in Pkg.test.

## Tests
Focused: `julia --project=<repo>/benchmark -e 'using Test; include("<repo>/test/chunk1b_radiation_shaped_probe.jl")'` -> **434/434** (5m40s; run2). Full unconditional: `julia --project=<repo> -e 'using Pkg; Pkg.test()'` -> **495/495** (61 Chunk 1a + 434), exit 0, 11m14s.
Sections: A primitives on both meshes (Fornberg vs textbook uniform stencils, topology/bandwidth, degree 0-4 exactness at 1e-11/1e-10, degree 5 must FAIL, lifting vs analytic cubic, signs/types); B coefficient/boundary/init/source maps (ForwardDiff vs hand-derived gradients, source vs central FD, continuum identity) and assembled operator derivatives wrt u, boundary values (directional), k, v, g;
C primal solve (retcode, 4 tolerances, componentwise in atol+rtol|u| units: Rodas5P < 1 (observed <=0.12), QNDF < 5 (observed <=1.9), plus proto/default-linsolve/dense variants); C2 structure evidence; D ForwardDiff full solve Jacobian vs analytic (col-scaled, worst index logged, Rodas5P < 10 reltol, QNDF < 100 reltol, monotone, all six columns > 1e-2); E Mooncake map projections vs ForwardDiff transpose (3 seeds, one preparation reused over 3 thetas, < 1e-11); F Mooncake solve stage; G cache reuse.
Calibration data (`sweep_diagnostic.log`, per mesh): Rodas5P ForwardDiff col-scaled error 7e-7/5.9e-7 (tol 1e-5) to 3.9e-12/6.7e-12 (1e-11); QNDF 1.2e-4 to 2e-10/4.2e-10. Mooncake-Gauss projected-gradient error, Rodas5P: 3.1e-8 -> 4.0e-13 (mesh a), 2.7e-7 -> 5.4e-13 (mesh b); in the full test, worst over 3 thetas at reltol 1e-7/1e-9/1e-11: mesh_a 1.1e-9/2.3e-11/3.3e-13, mesh_b 1.3e-9/2.2e-11/2.2e-13 (gate < 0.5 reltol, monotone, final < 1e-11).
Quadrature(MooncakeVJP) at the default tolerances (abstol 1e-10, reltol 1e-9): worst over 3 thetas 2.6e-11 (a), 2.8e-11 (b).
Preparation: one `prepare_gradient` per objective/tolerance reused over changed same-shape theta; G: two independent preparations, interleaved, repeat bitwise-equal, return to nominal bitwise-equal, each < 1e-10 vs analytic.

### The four original failures (run1) and their causes - no real gradient failure was relaxed
1-2. `L2L/L2R` support assertion (mesh_a, mesh_b): TEST BUG. I asserted the boundary node enters rows 1:4; the correct centered/clamped topology (already asserted elsewhere: row 3 starts at node 1) gives rows 1:2 and N-2:N-1 (`[1,2]`,`[6,7]` / `[9,10]` printed). Operator values were already passing. Assertion corrected to the actual stencil, still an exact index equality.
3-4. Quadrature gradient `< 1e-11` (2.59e-11 mesh_a, 2.83e-11 mesh_b): CALIBRATION error. That call uses the default solver tolerances abstol 1e-10/reltol 1e-9, so the ceiling 1e-11 was 0.01 reltol, tighter than the solver can promise; observed error is 0.03 reltol, comparable to Gauss at the same tolerance (2.3e-11). Ceiling changed to 1e-10 (= 0.1 reltol, ~4x above observation). Tight-tolerance Gauss results (3e-13) show there is no systematic gradient error.

## Benchmarks (`benchmark/chunk1b_benchmarks.jl`, run in `<repo>/benchmark`, 1 thread, Rodas5P(FD)+Sparspak+analytic jac, Gauss+MooncakeVJP, abstol 1e-10/reltol 1e-9, warmed, interpolated args, evals=1 for prepare/hot)
| | mesh_a primal | mesh_a prepare | mesh_a hot reverse | mesh_b primal | mesh_b prepare | mesh_b hot reverse |
|---|---|---|---|---|---|---|
| median | 3.02 ms | 8.07 ms | 100.7 ms | 6.78 ms | 11.2 ms | 156.4 ms |
| memory / allocs | 1.9 MiB / 38047 | 3.4 MiB / 50730 | 19.5 MiB / 301886 | 2.6 MiB / 43096 | 4.5 MiB / 56668 | 27.4 MiB / 361964 |
(benchmarks.log headings read "N=9, 8 unknowns"/"N=12, 11 unknowns": a mislabel since corrected in the script; the meshes are 9 points/7 unknowns and 12 points/10 unknowns, timings unaffected.)
Hot reverse / primal 33x (a), 23x (b). Toy 7-10 unknown system, one thread; NOT production performance. Prepare is warmed (no compile time). Large allocation counts in the reverse pass are accepted as unoptimised.

## Limitations / rejected paths
- Default solver-internal autodiff under Mooncake with the state-carrying RHS struct fails at `get_pf`/`prepare_pullback_cache`: "Differentiating through a FunctionWrappersWrapper whose wrapped function itself carries differentiable state" (Rodas5P and QNDF, `candidate_probe.log`). Workaround: `autodiff=AutoFiniteDiff()` (also needed on this time-dependent problem; tested). Not tried: restructuring the RHS to remove carried state.
- QNDF default linsolve under Mooncake: `Ptr{Mooncake.NoTangent}` / `@from_chainrules` error. QNDF(AutoFiniteDiff)+Sparspak runs, but its reverse gradient converges slowly (6.6e-5..6.6e-7 mesh a, 2.8e-4..5.7e-6 mesh b; full test 5.0e-5/1.1e-5/3.0e-6 at reltol 1e-7/1e-9/1e-11). **LIMITED / NOT ACCEPTED**: only a monotone-refinement diagnostic is tested.
- Linear problem only: finite-difference solver Jacobian is exact up to rounding here; nothing is shown for nonlinear state-dependent RHS.
- Manufactured-solution cancellation: the source is built to cancel the discretised operator action so that coefficient (k0,v0,g0) sensitivity enters the solve only through r(theta) and the source; coefficient/boundary/init/source derivatives are therefore tested separately at map/operator level (section B, E), not only through the solve.
- Dense W with Sparspak is a misconfiguration; analytic sparse Jacobian assembled on a fixed union pattern (not generic sparsity detection).
- Mesh sizes tiny (7-10 unknowns); Sparspak scaling to the real PDE size untested.

## Next step
Chunk 1c candidate: same harness with a nonlinear-in-u term (e.g. quadratic source) and sparse Jacobian from AD/analytic, to test whether solver-internal AutoFiniteDiff remains adequate for the adjoint and to try a stateless RHS closure with default autodiff.

Work stopped after Chunk 1b; no physics translation started. Awaiting supervisor verification.
