# SciML manual and current ODE performance/type-stability audit

Read-only assessment, 2026-10-03. No production callback, algorithm, tolerance,
grid, dependency, or physics was changed. The working tree includes the previously
validated fallback-specialization patch. Diagnostics are outside git under
`cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk17/`.

## Official manuals consulted

- [Code optimization](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/faster_ode_example/):
  reduce RHS allocations and redundant work; use an in-place RHS; consider static
  arrays for small systems; improve Jacobian construction; select methods by benchmarks.
- [Specialization levels](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/):
  `AutoSpecialize` trades runtime for latency; `FullSpecialize` is recommended when
  hot runtime matters. Installed-version behavior must also be checked directly.
- [Common solver options](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/):
  `saveat` already defaults `save_everystep` and `dense` to false. `calck=false`
  cannot be used with our required interpolated `saveat` outputs. `tstops` can
  address known discontinuities; it changes numerical stepping and needs validation.
- [Solver recommendations](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/):
  Rosenbrock methods are appropriate for small stiff systems. A different method
  must be compared at achieved output accuracy, not just identical tolerance numbers.

## Fresh type checks on the optimized working tree

JET 0.12.1 and `code_warntype`, Julia 1.12.6, existing isolated tools environment.
Output: `chunk17/static/{run.log,jet_report.txt,code_warntype_rhs.txt}`; exit 0.

| Concrete call | JET optimization reports |
|---|---:|
| 12-state recombination RHS | 0 |
| 7-state RHS, with/without feedback | 0 |
| time-gradient callback | 0 |
| PDE coefficients, step, integrals and setup kernels | 0 |
| allocating callback `jac5!` | 865 |
| public history, CosmoRec-module filter | 25 |

The explicitly flagged 7-/12-state RHS wrappers and assembled RHS return
`Vector{Float64}` in `code_warntype`. No `Any` or `Core.Box` appears in those
printed optimized call paths. This is evidence for the tested concrete calls,
not proof of all parameter/type/branch combinations.

The Jacobian reports include runtime chunk/configuration selection. ForwardDiff
has an inner function barrier that restores concrete Dual types; these reports
must not be interpreted as 865 executed errors or all RHS calls being dynamic.
The 25 public-history reports include low-frequency `Any` containers and setup/
fallback paths. See the earlier profiling analysis for measured hotness.

Thus the hot RHS is concretely inferred, but the entire orchestration/callback
path is not fully type-grounded. Type stability also does not imply allocation
freedom: the validated 12-state RHS still allocates 1440 B/20 objects per call;
the 7-state RHS allocates 1040 B/14 objects.

## Which ODE is most expensive?

`chunk17/ode_work.jl` records actual solver blocks from the public three-pass,
two-PDE calculation. The recording callback reproduces public Xe/Tm bitwise.
Then BenchmarkTools replays groups with their original block inputs, warm,
five samples, evals=1, pinned to P-core CPU 8. These are block-solve costs,
not full-pass assembly or a before/after speedup claim.

| Pass | System | Blocks | RHS evaluations | Jacobians | Group median |
|---|---|---:|---:|---:|---:|
| 0 | 12-state H+He | 27 | 180126 | 22416 | 0.627 s |
| 0 | 7-state H only | 9 | 53106 | 6587 | 0.117 s |
| 1 | 12-state H+He | 27 | 176766 | 21966 | 0.669 s |
| 1 | 7-state H only | 9 | 53058 | 6575 | 0.123 s |
| 2 | 12-state H+He | 27 | 178606 | 22242 | 0.677 s |
| 2 | 7-state H only | 9 | 53010 | 6560 | 0.123 s |

- H+He accounts for about 77% of solver-counted RHS evaluations and 85% of
  these measured block-solve times. It is about 5.4 times the H-only segment's
  total cost. That is not a claim that each 12-state call is 5.4 times slower.
- The three passes have similar cost; there is no single dominant iteration.
- The four H+He blocks covering roughly z=1869 down to 1672 account for
  45–46% of each pass's RHS evaluations. The first H-only block after the
  switch (z~1681 to 1484) contributes another 15%.
- About 29000 Jacobians are constructed per pass, with fewer than 200 rejected
  steps. The work is predominantly accepted small steps, not repeated failure.
- The Recfast tail is millisecond-scale in the earlier stage benchmarks;
  it is not the relevant optimization target. Nor are the small LU systems.

## Actual saving and specialization, not assumed defaults

The current function is declared `AutoSpecialize`. Inspecting an initialized
Rodas5P integrator shows an actual `FunctionWrappersWrapper` RHS on this installed
version. Its options are `dense=false`, `save_everystep=false`, `calck=true`.
We already avoid storing every accepted step. Disabling output interpolation
would not preserve the requested history and is not a free optimization.

## Bounded candidate probes (not installed changes)

`chunk17/callback_probe.jl` tests two possibilities. Completed corrected output:
`chunk17/work/callback_retry/callback_probe.txt`, exit 0. The first attempt
incorrectly treated the stored pre-reset switch state as a 7-state state; that
diagnostic failure is retained in `work/callback_run.log`. The retry applies the
exact public helium reset/packing before testing the H-only block.

### Reuse a prepared, fixed-chunk state-Jacobian configuration

At representative real states, using `ForwardDiff.jacobian!` and a reusable
`JacobianConfig` gives bitwise-identical matrices to the current callback:

| State size | Current median | Cached median | Current bytes | Cached bytes |
|---|---:|---:|---:|---:|
| 12 | 4.03 us | 3.10 us | 13200 | 7168 |
| 7 | 2.09 us | 1.32 us | 8064 | 4416 |

Each uses 1000 BenchmarkTools samples, evals=1, compilation excluded. This saves
about 23%/37% of the measured callback time and 45% of its allocation bytes.
It does not imply the same reduction in full prediction time. The output function
inside this probe still allocates RHS temporaries.

This is a primal state-Jacobian test with p and z held fixed, NOT proof that a
shared cache remains differentiable through an entire solve. A real implementation
needs per-solve/per-chain caches with the correct element/tag types, updated p/z,
and ForwardDiff/Mooncake validation before adoption.

### Explicit `FullSpecialize`

On the actual first H+He block and first post-reset H-only block, explicit
`ODEFunction`/`ODEProblem` FullSpecialize gives bitwise-identical saved states,
RHS counts and accepted-step counts:

| Block | Auto median | Full median | RHS calls / accepted steps, both |
|---|---:|---:|---|
| 12-state, z=3000 to 2951 | 20.34 ms | 18.60 ms | 7010 / 865 |
| 7-state, z=1681 to 1484 | 69.90 ms | 61.69 ms | 35018 / 4377 |

Five samples per configuration, same process/core/inputs. These roughly 9–12%
block-level improvements are promising exploratory timings, not a thermally
matched full-history benchmark. More specialization can increase compilation.
This would change the caller-provided solver configuration, not a library default.

## Recommended next steps

1. Validate an explicit FullSpecialize callback over complete histories and
   derivative paths; benchmark full predictions with the existing thermal protocol.
2. Prepare/reuse state-Jacobian buffers/configuration per solve, preserving AD
   element types and independent mutable workspaces for concurrent chains.
3. Remove remaining small RHS temporaries with an AD-safe workspace (or benchmark
   an appropriate small-static-array design). StaticArrays would be a new dependency
   and is not justified merely because the manual shows a speedup on another ODE.
4. Only then study alternative stiff solvers/restart layout or tolerance work-precision
   curves. Compare achieved Xe/Tm/Cl accuracy; do not loosen the target or remove
   the original clamp just to run faster. Sparse/Krylov machinery is not indicated
   for these 7/12-state systems with currently small linear-algebra cost.

The checks and prototypes above do not modify production source/tests/settings,
do not validate true reverse through Rodas5P, and do not implement full cosmological
initialization. No commit or push was made for this assessment.
