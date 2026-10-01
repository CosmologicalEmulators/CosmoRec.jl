# Chunk 1a results (repaired): stiff SciML solve <-> AD compatibility probe

Toy problem only (`y1'=-a y1+b y2`, `y2'=-c y2`, theta=[a,c,b,y10,y20], nominal [1000,1,50,1,1]);
reference is the closed form. No CosmoRec physics, no native fixtures (contract: `docs/fixture_contract.md`;
native original-code text fixtures + unit tests remain mandatory before each later physics stage).

Checkout `/home/marcobonici/Desktop/work/CosmologicalEmulators/CosmoRec.jl`, branch `develop`, HEAD `80e3c6a`
(nothing committed/staged). Julia 1.12.6. Resolved: SciMLSensitivity 7.119.12, SciMLBase 3.57.0,
DiffEqBase 7.21.3, OrdinaryDiffEqCore 4.18.1, OrdinaryDiffEqDifferentiation 3.12.4,
OrdinaryDiffEqRosenbrock 2.7.5, OrdinaryDiffEqBDF 2.4.12, Mooncake 0.5.61, DifferentiationInterface 0.7.21,
ADTypes 1.24.0, ForwardDiff 1.4.6, BenchmarkTools 1.8.0, ChainRulesCore 1.26.1.
Test env versions: `chunk1a_sonnet55/pkgtest_resolved_versions.txt`; benchmark env: `benchmark/Manifest.toml` (untracked/ignored; seeded from the same resolution).

## Files
Project.toml (test deps only via `[extras]`/`[targets]`, `julia = "1.12"` = the only version tested; no `[deps]`),
src/CosmoRec.jl (empty stub), test/runtests.jl, test/chunk1a_helpers.jl (shared toy), test/chunk1a_stiff_ad_probe.jl,
docs/fixture_contract.md, this file, benchmark/{Project.toml, chunk1a_benchmarks.jl, chunk1a_route_diagnostic.jl,
chunk1a_sweep_diagnostic.jl, chunk1a_reverse_sweep_diagnostic.jl, chunk1a_rodas_default_blocker.jl}.
`docs/REVIEW_CHUNK1A.md` is supervisor-owned and untouched. `Project.toml.bak` removed (generated backup); the stale root Manifest removed (ignored, generated).

## 1. Which derivative route runs (runtime evidence)
`julia --project=<repo>/benchmark benchmark/chunk1a_route_diagnostic.jl {none|count|throw}` (logs `chunk1a_sonnet55/route_diagnostic_*.log`).
Process-local: the script redefines `SciMLSensitivity.get_paramjac_config(::MooncakeLoaded, ::MooncakeVJP, ...)`
(resolved source `SciMLSensitivityMooncakeExt.jl:53-62`, body re-implemented verbatim for `count`, throwing a sentinel for `throw`); no installed source edited, no override in `Pkg.test()`.
- `count`: hook called 2x in `prepare_gradient`, 2 more per `gradient`, for both GaussAdjoint and QuadratureAdjoint (`MooncakeVJP`); results unchanged. With `sensealg=nothing` hook calls = 0.
- `throw`: sentinel reached during preparation; stack (relevant frames): `get_paramjac_config <- _adjoint_sensitivities <- adjoint_sensitivities <- adjoint_sensitivity_backpass` (Gauss) and `... adjointdiffcache <- ODEQuadratureAdjointSensitivityFunction <- ODEAdjointProblem ...` (Quadrature).
- Conclusion: with `sensealg=GaussAdjoint/QuadratureAdjoint(autojacvec=MooncakeVJP())` the SciMLSensitivity continuous-adjoint backpass runs and builds Mooncake pullback caches (RHS-VJP via Mooncake). The keyword is not decorative. The default-`sensealg` path (`nothing`) does not reach the MooncakeVJP hook; its route was not characterised here and is not used.
- Rodas5P blocker is consistent with this: with default Rodas5P the failure `FirstAutodiffTgradError` occurs inside the adjoint's backward `solve` (Rosenbrock `calc_tderivative!`, called from `_adjoint_sensitivities`), because the adjoint RHS contains the Mooncake VJP. Fix = solver-internal `autodiff=AutoFiniteDiff()` (`rodas_default_blocker.log`). **Caveat:** the toy RHS is linear, so a finite-difference Jacobian is exact to rounding; adequacy for the nonlinear physics RHS is untested and must be re-checked (or an analytic/AD-compatible time gradient supplied) in later chunks.

## 2. Tests (all against the analytic oracle)
Commands: `julia --project=<repo> -e 'using Pkg; Pkg.test()'` -> **61/61 pass, exit 0** (4m32s; `pkgtest.log`).
Focused: `julia --project=<repo>/benchmark -e 'using Test; include("test/chunk1a_stiff_ad_probe.jl")'` -> 61/61 (2m50s; `focused_test.log`).
(Earlier 27-assertion version superseded: its FD gate with ceiling 1.0 was removed.)

Error metrics: Jacobian: column-scaled absolute error `|dJ_ij|/max_i|J_ref,ij|` (no relative floor; true-zero entries such as dy2/dy... are judged by absolute deviation vs their column scale). Gradient of a seeded projection `w'g`: `max|g-ref|/max|ref|`.

### ForwardDiff through public `solve` vs analytic Jacobian (col-scaled), tolerances (abstol,reltol)
| tol | Rodas5P | QNDF |
|---|---|---|
| 1e-6/1e-5 | 8.4e-5 | 1.6e-4 |
| 1e-8/1e-7 | 1.1e-7 | 2.6e-6 |
| 1e-10/1e-9 | 9.0e-11 | 1.4e-7 |
| 1e-12/1e-11 | 2.0e-11 | 3.4e-9 |
Gates: monotone decrease, `err < 50*reltol` (Rodas5P) / `500*reltol` (QNDF; calibrated from this sweep), final < 1e-9 / 1e-7. Worst entry throughout: row 1 (y1 early), column a (dy1/da, J=-3.0777e-4). Full per-column/worst-element/FD data for tolerances down to 1e-13: `sweep_diagnostic.log`.

### Central FD vs analytic (tight tol 1e-12/1e-11), col-scaled
- Rodas5P: h=1e-2: 1.6e-5, 1e-3: 1.6e-7, 1e-4: 1.6e-9 (truncation O(h^2), factor 100/decade), h=1e-5: ~1.7e-9 (noise floor begins), h=1e-6: 2.8e-8. Gate: errors fall >20x per decade for h=1e-2..1e-4 and error(1e-4) < 1e-8. **Earlier version's "min over a nonconvergent sweep" metric removed.**
- QNDF: h=1e-2 1.6e-5 (truncation), 1e-3 8.4e-7, 1e-4 1.9e-5, 1e-5 1.7e-5, 1e-6 1.0e-4; at looser tolerances the noise is far larger (e.g. abstol 1e-8: h=1e-5 gives 0.12, h=1e-6 gives 5.7 col-scaled). The QNDF numerical output is not smooth in theta at FD step scales (step/order selection), so FD of the solve is noise-dominated (error ~ tol/h). **This is an FD-reference failure, not a ForwardDiff failure** (ForwardDiff vs analytic converges, table above). The earlier 0.28 relative error was this noise amplified by the 1e-8 relative floor on near-zero entries. QNDF FD is therefore a retained diagnostic, not a gate; it is not accepted as a reference.

### Mooncake + continuous adjoint vs analytic gradient (Gauss, `max over 3 thetas` of projected-gradient error)
| tol | Rodas5P(autodiff=AutoFiniteDiff) | QNDF |
|---|---|---|
| 1e-8/1e-7 | 6.5e-11 | 1.2e-7 |
| 1e-10/1e-9 | 2.0e-12 | 5.4e-9 |
| 1e-12/1e-11 | 4.2e-14 | 6.3e-10 |
QuadratureAdjoint at 1e-12/1e-11: Rodas5P 5e-14, QNDF 4.5e-11 (test). Full sweep incl. 1e-6 and both adjoints: `reverse_sweep_diagnostic.log` (all converge, 1.7e-5 worst at 1e-6 for QNDF). thetas: nominal, [1200,0.8,40,1.3,0.7], [900,1.5,60,0.6,1.4]; objective builds theta-dependent u0/p inside, calls `solve` with `saveat` (interpolated observation at s=0.02), dot with seeded w.
Gates: error < 1e3*reltol each, monotone refinement, final < 1e-10 (Rodas5P) / 1e-8 (QNDF).

### Prepared cache reuse
Two independently prepared objectives (w1,w2), one preparation each reused over theta sequence nominal, theta2, theta3, nominal: repeated call bitwise equal (`g1 == g1b`), each gradient < 1e-10 vs analytic, return to nominal reproduces first gradient bitwise.

## 3. Candidate classification
- Rodas5P (+AutoFiniteDiff solver autodiff for reverse): accepted for ForwardDiff and Mooncake+Gauss/Quadrature (toy only; nonlinear caveat above). Default-autodiff Rodas5P: rejected under Mooncake reverse (blocker).
- QNDF: ForwardDiff vs analytic converges but ~20-100x less accurate than Rodas5P at equal tolerance; Mooncake Gauss/Quadrature converge to analytic (6e-10 / 5e-11 at tight tol). FD cross-check unusable (noise). Status: limited candidate, passes the analytic gates above; no claim beyond this toy.
- `MooncakeAdjoint` (discrete, explicit-only per docs): not used.

## 4. Benchmarks (BenchmarkTools, Julia 1.12.6, 1 thread; `benchmark.log`, exit 0)
Command: `julia --project=<repo>/benchmark benchmark/chunk1a_benchmarks.jl`. Configuration = tested: tight tol, Gauss+MooncakeVJP, Rodas5P(autodiff=AutoFiniteDiff) and QNDF. `prepare` is warmed (compile excluded) and NOT cold compilation; compile-inclusive first call not measured. Prepared reverse and prepare use `evals=1`.
| | primal | prepare_gradient (warmed) | hot prepared reverse | reverse/primal |
|---|---|---|---|---|
| Rodas5P(FD) | median 0.253 ms, 13.2 KiB, 246 allocs | median 1.12 ms (mean 23.6 ms, outliers), 727.6 KiB, 9878 allocs | median 8.08 ms, 3831.8 KiB, 74174 allocs | 31.9 |
| QNDF | 0.274 ms, 15.4 KiB, 272 | 1.67 ms, 1477.4 KiB, 20636 | 6.27 ms, 2888.6 KiB, 62445 | 22.9 |
Warm `prepare_gradient` is cheaper than one gradient because Mooncake memoises compiled rules; do not read it as total set-up cost.

## 5. Retained failures / open items
- Terminated earlier benchmark attempts (SIGTERM traces from the killed session, one orphaned process I killed with SIGKILL) preserved as `chunk1a_sonnet55/benchmark_terminated_attempt.log` and `/tmp/benchmark_run1.log`; not an algorithm bug.
- Rodas5P default-autodiff blocker (above); linear-RHS caveat for the AutoFiniteDiff workaround.
- QNDF FD noise (above). Hook-route evidence is for MooncakeVJP config only; the default `sensealg=nothing` Mooncake route not analysed.
- Future physics stages need native text fixtures and tolerance calibration; none exist.

STOPPED AFTER CHUNK 1a; no Chunk 1b or physics translation started.
