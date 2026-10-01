# Chunk 1c acceptance: nonlinear scalar recombination-shaped ODE

The supervising agent independently ran the complete unconditional suite on 2026-10-01 EDT: `julia --project=<repo> -e 'using Pkg; Pkg.test()'` -> exit 0, **603/603 tests**, 12m46.2s test execution. This completes only the approved numerical-AD toolchain probes. The package `src/` remains an empty stub: there is still NO translated CosmoRec physics.

## Independent evidence

- Chunk1a stiff analytic linear ODE and Chunk1b radiation-shaped sparse PDE gates remain unconditional; the full 603-result run includes both.
- Chunk1c focused repaired suite passed108/108; first run with three exact-Float64 assertion failures remains in the logs. Repairs changed only floating equality to tight isapprox checks.
- The inspected probe uses `x'=-a*x^2`, `y'=-c*y`, with independent oracle `x=x0/(1+a*x0*s)`, `y=y0 exp(-c*s)`, theta `[a,c,x0,y0]`. The nonlinear state Jacobian is explicitly checked. Initial state and parameters are built inside the differentiated public solve objective; early/interpolated/late outputs, seeded projections, tolerance refinement, cache reuse and independent preparations are tested.
- Claimed results reproduced in the independent full-suite output: ForwardDiff/Rodas5P(AutoFiniteDiff) and Mooncake/Gauss errors converge against analytic derivatives with solver tolerance. Solver-internal finite differencing is distinct from MooncakeVJP outer sensitivity. QNDF is limited; do not read its permissive factor as evidence of equal accuracy to Rodas.
- Reproduction of the **combined original-ODE `jac+tgrad` segfault is an intentionally retained failure** (`benchmark/chunk1c_hook_blocker.jl` and logs), not a supported route and not part of the passing package suite. Exact trigger attribution among Julia/Mooncake/SciML remains unresolved. The implementation avoids the combination. Don't claim the compiler crash was fixed.
- Independently verified source edits are confined to the expected `test/`, `benchmark/`, and results-document locations. The original model sources and `src/` remain untouched. No commits/staging/pushes.

## Limits

This proves one smooth scalar quadratic-loss system. It is not a broad stiff nonlinear solver certification, not a coupled multilevel recombination RHS, and not a nonlinear radiation-PDE result. In particular `AutoFiniteDiff` adequacy for CosmoRec's table/interpolation/coupled RHS remains unknown. No native-code fixture exists yet.

Original worker/full-suite data and preserved segfault/first-test logs are under the local analysis root `chunk1c/` and `chunk1c_resume01/`. The independent 603-test log is `chunk1c_resume01/supervisor_pkgtest.log`; worker physics/results report is `docs/CHUNK1C_RESULTS.md`.

## Stop point

Chunk 1c is accepted **only for its analytic toy scope**. Before any native-physics translation, the next task must freeze small plain-text fixtures from the pinned CosmoRec/CAMB sources, resolve the accuracy/build configuration and preserve the known H(z) endpoint behavior as an explicit baseline. Then translate and validate only the first physics component against those fixtures. No chunk4, ODE, PDE, or package implementation may claim native-CosmoRec fidelity before that reference gate.
