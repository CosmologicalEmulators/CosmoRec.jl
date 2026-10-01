# Chunk 1a supervising review: not yet accepted

The supervisor independently ran `julia --project=/home/marcobonici/Desktop/work/CosmologicalEmulators/CosmoRec.jl -e 'using Pkg; Pkg.test()'`: all 27 assertions passed, exit 0, approximately three minutes of test execution. This verifies the current assertions, not that the requested scientific gates are adequate.

## P1: Do not turn a failed gradient check into success by widening its threshold

`test/chunk1a_stiff_ad_probe.jl:187-243` accepts QNDF finite-difference consistency with `fd_ceiling["QNDF"]=1.0`. Reported best error is approximately 0.284 under the current elementwise relative/floor metric. The near-zero-entry floor and observable/parameter scaling must be examined before interpreting that number, but accepting up to 100% disagreement is not a meaningful gate.

Report absolute discrepancies, the worst row/column and derivative values, scaled parameter-direction errors, and scientifically interpretable per-output or per-column errors. Compare the analytic expression, ForwardDiff of the numerical solve, and finite differences at several step sizes and solver tolerances. Show the finite-difference truncation/noise regime rather than merely taking the smallest error in a nonconvergent sweep. A tiny true derivative should be judged with an explicit absolute criterion, not a misleading relative criterion or a denominator chosen to make tests green.

QNDF can be an explicitly rejected/limited candidate in the diagnostics if evidence supports that decision. Its failed gate must not be presented as an accepted AD path. At least one stiff-solver path must pass meaningful analytic and AD checks before this chunk is accepted. Do not claim convergence from only one loose/tight pair when the requested sensitivity comparison has not been swept.

## P1: Prove which solve sensitivity path actually runs

The current Mooncake result agrees with ForwardDiff to approximately 1e-13, yet the default Rodas5P preparation failed in the solver-internal ForwardDiff time-gradient path and was replaced with AutoFiniteDiff. Neither fact alone proves whether the selected GaussAdjoint/MooncakeVJP rule is active or whether the numerical solver was directly traced.

Identify the resolved DiffEqBase/SciMLSensitivity/Mooncake extension and selected solve rule. Obtain runtime evidence that the continuous adjoint and MooncakeVJP construction/execution occur, not simply that `sensealg` was passed as a keyword. A separate fresh-process, process-local diagnostic can instrument or deliberately make the narrow SciMLSensitivity MooncakeVJP configuration path raise a sentinel, then check whether preparation/evaluation reaches it. Do not edit installed package source. Retain the diagnostic and its output; do not leave global method overrides in the normal package tests.

Correct the misleading "same discrete solve"/"Discrete-vs-discrete" annotations at lines 309-314. If the real path is direct discrete differentiation, report it honestly; that may still satisfy the human's Mooncake-through-stiff-solve requirement, but it does not establish the proposed continuous-adjoint route. Do not invent a failure or a success for either route.

Parameter-dependent initial states and coefficient construction must remain in the objective. Repeated prepared evaluation at changed inputs must be checked against the analytic oracle, not only against another AD path that could reproduce the same numerical artifact. Add a real reverse-error tolerance/refinement sweep if testing continuous sensitivities.

## P1: Benchmark the tested configuration and finish the handoff

`benchmark/chunk1a_benchmarks.jl` uses default Rodas5P for the Mooncake objective, while the accepted test configuration uses `MOONCAKE_ALG=Rodas5P(autodiff=AutoFiniteDiff())`. Match them. Factor the shared toy definition into a small helper rather than including the entire testset in the benchmark. Use BenchmarkTools interpolation, explicit sample limits and evals=1 where prepared mutable state is reused. Report hot primal, prepared-cache construction, hot prepared reverse, allocations and limitations. A warmed BenchmarkTools construction measurement is not cold compilation and must not be called "cold-ish" compilation.

Previous benchmark logs contain interrupted/terminated Julia traces; preserve the failures and distinguish them from an algorithm bug. The last Claude response "The benchmark is running in the background. I'll resume once it finishes" is not a final chunk report. Do not leave an unmanaged background process as a success report. Either wait for a bounded benchmark to finish and inspect its result, or explicitly report a blocker with logs and exit status.

## P2: Dependency and artifact hygiene

The module currently contains no runtime code, yet Project.toml duplicates test/benchmark packages in both deps and extras, including Test and BenchmarkTools. Use proper extras/targets for test-only dependencies. Put the standalone benchmark in a persistent local benchmark project, developed against the package if needed, instead of creating a /tmp environment or adding test-only dependencies to the package runtime. Record resolved versions used for each result.

Read and remove the accidental Project.toml.bak if it contains only generated backup data; do not commit generated backup/config artifacts. Honor the actual tested Julia/dependency compatibility rather than asserting unsupported compatibility from the skeleton. Keep the repository on develop, no commits/pushes/worktrees. No edits to unrelated projects or installed Julia packages.

## Required repair order and stop gate

1. Resolve the sensitivity-route evidence first. If there is a real unsupported-route blocker, reduce/report it before further expansion.
2. Repair analytical/scaled gradient acceptance tests and candidate status.
3. Fix test dependency/helper organization and run focused tests, then Pkg.test().
4. Run the matching bounded BenchmarkTools measurements in a persistent local environment.
5. Write an actual `docs/CHUNK1A_RESULTS.md` containing source/versions, exact commands, errors and convergence tables, accepted/rejected paths, reuse checks, benchmarks, and retained failures. No unverifiable claim that every candidate passed.

STOP after repaired Chunk 1a. No radiation prototype, no CosmoRec physics translation, no native fixture generation yet. Every later translated stage must have original-code text fixtures and unconditional tests; the existing toy analytic oracle does not replace those future fixtures.
