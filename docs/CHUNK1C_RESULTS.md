# Chunk 1c results: nonlinear-state stiff ODE compatibility probe

Scope: ONE two-state, four-parameter quadratic-loss toy. Nothing here establishes general nonlinear compatibility, a nonlinear PDE, or anything about CosmoRec physics. No native fixtures. Branch `develop`, HEAD `80e3c6a`, nothing staged/committed; Julia 1.12.6; versions as in `docs/CHUNK1B_RESULTS.md` (SciMLSensitivity 7.119.12, SciMLBase 3.57.0, OrdinaryDiffEqRosenbrock 2.7.5, OrdinaryDiffEqBDF 2.4.12, Mooncake 0.5.61, DifferentiationInterface 0.7.21, ForwardDiff 1.4.6, BenchmarkTools 1.8.0; no dependency change in Chunk 1c).

## Files
`test/chunk1c_helpers.jl`, `test/chunk1c_nonlinear_probe.jl`, `test/runtests.jl` (adds the include), `benchmark/chunk1c_benchmarks.jl`, `benchmark/chunk1c_hook_blocker.jl` (reduced reproducer), this file. Chunk 1a/1b files untouched.
Logs: `cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk1c/` (focused_test_run1.log = first run with 3 failures, kept; focused_test_run2.log; calibration_sweep.log; benchmarks.log) and `.../chunk1c_resume01/` (pkg_test_full.log/.rc, hook_blocker_{j,tg,jtg}.log). `chunk1c/` has no usable full-test log: the first Pkg.test invocation was cut by the walltime and is neither a pass nor a failure.

## Problem and oracle
theta=[a,c,x0,y0]=[1000,1,1,1]; `x' = -a x^2`, `y' = -c y`, s in [0,0.5]. Solving `dx/x^2 = -a ds`: `x(s)=x0/(1+a x0 s)`, `y(s)=y0 e^{-cs}`. The RHS is nonlinear in x; Jacobian `diag(-2 a x, -c)` is state dependent; the RHS is autonomous (df/ds = 0). No forcing; the nonlinearity is not cancelled.
Observation at saveat s = 0.0005 (early), 0.02 (interior), 0.5 (late), flat `[x,y]` x 3 = 6 outputs. Column maxima of the analytic Jacobian: [2.2e-4, 0.30, 0.44, 1.0] (a, c, x0, y0), so all four parameters are resolvable (dx/da at s=0.5 is -x0^2 s/(1+a x0 s)^2 ~ -2e-6, non-zero and compared with a column-scaled metric). u0 and p are built from theta inside every objective; `solve` uses `saveat`; outer arithmetic is `dot(w, observation)` with seeded w.

## Tests (focused 108/108, full 603/603)
Focused: `julia --project=<repo>/benchmark -e 'using Test; include("<repo>/test/chunk1c_nonlinear_probe.jl")'` -> **108/108**, 4m00s (run2). First run: 105 pass, 3 fail, all exact `==` on Float64 values (0.7^2 = 0.48999999999999994); changed to `isapprox` with rtol 1e-15; not a numerical problem.
Full unconditional: `julia --project=<repo> -e 'using Pkg; Pkg.test()'` -> **603/603** (61 Chunk 1a + 434 Chunk 1b + 108 Chunk 1c), exit 0, 12m54s test time / 13m03s wall (`chunk1c_resume01/pkg_test_full.log`).
Coverage: RHS values, mutation exercised (rhs! writes into du; hooks write into J/dT; no global caches), ForwardDiff state/parameter derivatives vs hand values, analytic Jacobian vs ForwardDiff, state dependence, zero time gradient of the ORIGINAL rhs, oracle ODE residual and hand derivatives; Mooncake VJPs of 4 seeded RHS projections vs ForwardDiff (one preparation reused over 3 (u,p) points, < 1e-13); primal; ForwardDiff solve Jacobian; central FD; Mooncake+Gauss; Quadrature; two independent preparations.

## Numerical results (abstol/reltol sweep reltol = 1e-5, 1e-7, 1e-9, 1e-11; abstol = reltol/10)
Solver-internal autodiff is a separate choice from the outer AD. "Rodas5P(AutoFiniteDiff)" = `Rodas5P(autodiff=AutoFiniteDiff())`, default linear solver (dense 2x2); "Rodas5P default" = `Rodas5P()`; "QNDF" = default `QNDF()`.

| quantity | Rodas5P(AutoFiniteDiff) | Rodas5P default | QNDF |
|---|---|---|---|
| primal componentwise error / (atol+rtol\|u\|) | 0.03-0.11 | 0.03-0.09 | 1.3, 3.5, 12.6, 29.9 (grows as tol tightens) |
| ForwardDiff(solve) col-scaled error vs analytic | 4.5e-6, 3.3e-8, 3.4e-10, 5.5e-12 | 4.5e-6, 3.3e-8, 3.4e-10, 5.7e-12 | 2.4e-4, 5.0e-7, 1.4e-8, 3.6e-10 |
| Mooncake+Gauss(MooncakeVJP), max proj. gradient error over 3 theta | 2.3e-6, 1.8e-8, 8.9e-11, 5.7e-13 | FAILS (see below) | 4.6e-6, 6.4e-8, 4.6e-9, 3.1e-11 |
| Mooncake+Quadrature(MooncakeVJP) | 8.9e-11 (1e-9), 5.7e-13 (1e-11) | - | not tested |
Gates: Rodas primal < 0.5; ForwardDiff < 1.0 x reltol (Rodas; observed <= 0.45) and < 50 x reltol (QNDF; observed <= 36); Gauss < 0.5 x reltol (Rodas; observed <= 0.23) and < 10 x reltol (QNDF; observed <= 4.6); monotone refinement for all. The QNDF gates are calibrated from the observation, not a derived bound. Per-theta ForwardDiff projected-gradient vs analytic also tested at 3 thetas. Thetas: nominal, [800,1.5,1.3,0.7], [1500,0.6,0.8,1.4].
Central finite differences of the Rodas5P(AutoFiniteDiff) solve (tol 1e-12/1e-11), col-scaled vs analytic: h=1e-2: 1.85e-5, 1e-3: 1.85e-7, 1e-4: 1.9e-9, 1e-5: 2.0e-9. Truncation O(h^2) (x100 per decade) until about h=1e-4, then a noise/floor ~2e-9 (h=1e-6: 6e-9 in calibration). Gate: >50x drop per decade for h=1e-2..1e-4 and error(1e-4) < 1e-8. FD is secondary and is not claimed to be better than ForwardDiff.
QNDF is classified LIMITED / not accepted as an accuracy path: its primal error is not controlled to the requested tolerance (up to ~30x) and ForwardDiff error is up to ~36x reltol; gradients still converge monotonically and agree with the oracle. (QNDF with `AutoFiniteDiff` gives the same ForwardDiff errors and similar Gauss errors, 3.9e-6, 1.1e-7, 1.2e-9, 3.1e-11, in `calibration_sweep.log`; not part of the package tests.)

## Prepared cache reuse
For each objective and tolerance one `prepare_gradient` is reused over 3 changed same-shape theta. Separately two independent preparations (w1, w2; Rodas5P(AutoFiniteDiff), Gauss+MooncakeVJP, default tolerances abstol 1e-12/reltol 1e-11) were called interleaved over nominal, theta2, theta3, nominal: repeated call bitwise equal, each < 1e-10 vs analytic, return to nominal bitwise equal to the first result.

## Hooks, solver-internal autodiff and what the backward solver needs
- Accepted configuration: `Rodas5P(autodiff=AutoFiniteDiff())`, `GaussAdjoint(autojacvec=MooncakeVJP())`, outer `AutoMooncake(; config=nothing)` through DifferentiationInterface. Solver-internal AutoFiniteDiff (a Jacobian/time-gradient policy of the Rosenbrock integrator) is not the outer AD.
- `ODEFunction(...; jac=nl_jac!)` supplies the analytic state Jacobian `diag(-2 a x, -c)` of the ORIGINAL (autonomous) RHS to the forward integrator; `tgrad=nl_tgrad!` supplies its time gradient, which is exactly zero because that RHS is autonomous. Tested: `jac` hook primal (tolerance sweep), ForwardDiff (col-scaled < 1e-10), Mooncake+Gauss (proj. error < 1e-11, exploratory 3.9e-14). With a Dual-valued `p` the hooks are Dual-compatible.
- The BACKWARD adjoint ODE is constructed by SciMLSensitivity; the original hooks are NOT evidence that it has a Jacobian or time-gradient. Its time dependence is real (the adjoint RHS depends on the stored forward trajectory), so a zero time-gradient must never be supplied for it, and none is. What the backward solve needs still comes from solver-internal autodiff: with default `Rodas5P()` (with or without the `jac` hook) the Mooncake+Gauss gradient fails with "First call to automatic differentiation for time gradient" (same error class as Chunk 1a; tested in the package tests via the message), and `AutoFiniteDiff` is the workaround. This inference (hooks do not reach the backward solve) rests on the unchanged failure with the `jac` hook present; the SciMLSensitivity source was not edited or further audited.
- Default `QNDF()` runs under Mooncake on this problem (no linsolve/state-carrying wrapper in this toy, unlike the Chunk 1b PDE).
- `AutoFiniteDiff` adequacy here relies on a 2x2 smooth RHS; it is not shown for stiffer, larger or non-smooth nonlinear RHS.

## Retained failures and limits
1. Julia compiler segfault, NOT FIXED, only avoided: `ODEFunction(rhs!; jac=..., tgrad=...)` with BOTH hooks, Rodas5P(AutoFiniteDiff), Gauss+MooncakeVJP, during `prepare_gradient` -> SIGSEGV 11 in `decay_derived` (cgutils.cpp) from `emit_ccall`, reached from Mooncake `_ssa_to_ids` (stack in `chunk1c_resume01/hook_blocker_jtg.log`). Reduced reproducer `benchmark/chunk1c_hook_blocker.jl`: mode `j` (jac only) works, `tg` (tgrad only) works, `jtg` (both) segfaults (exit 139). Not isolated further (Julia vs Mooncake vs SciMLBase); the accepted route does not set both hooks, and the package tests do not include that combination.
2. An earlier, different segfault in the same stack signature came from my own helper: `Val(hooks)` of a RUNTIME Symbol inside the differentiated objective (type-unstable dispatch). Fixed in the helper by taking `hooks::Val` as argument. Observed with the same compiler/Mooncake stack, so unstable dispatch inside a Mooncake-prepared objective should be avoided (not further characterised).
3. Default `Rodas5P()` fails under Mooncake+Gauss (time-gradient error above); retained as a tested expected failure.
4. QNDF limited (above). QuadratureAdjoint is tested only at the two tightest tolerances with Rodas5P(AutoFiniteDiff).
5. Scalar-parameter quadratic toy only; no PDE coupling, no state-dependent coefficients beyond a x^2, no large systems or sparse path in this chunk (the sparse structured path was Chunk 1b, linear only). No claim about the CosmoRec population equations.

## Benchmarks (`benchmark/chunk1c_benchmarks.jl`, `<repo>/benchmark` env, 1 thread, `chunk1c/benchmarks.log`)
Rodas5P(AutoFiniteDiff), Gauss+MooncakeVJP, abstol 1e-12/reltol 1e-11, 3 saveat times; interpolated args; prepare and hot use evals=1; prepare is WARMED (no compile time measured; no cold-compilation claim).
| | median | min | mean | memory | allocs | samples |
|---|---|---|---|---|---|---|
| primal | 0.657 ms | 0.645 ms | 0.666 ms | 13.2 KiB | 246 | 100 |
| prepare_gradient (warmed) | 1.950 ms | 1.799 ms | 30.008 ms (outlier(s) in mean) | 1495.4 KiB | 19689 | 15 |
| hot prepared reverse | 10.036 ms | 8.680 ms | 10.201 ms | 5291.5 KiB | 97393 | 60 |
Hot reverse / primal 15.3. Toy two-state timings, not production performance.

Work stopped after Chunk 1c; no native physics translation started. Awaiting supervisor verification.
