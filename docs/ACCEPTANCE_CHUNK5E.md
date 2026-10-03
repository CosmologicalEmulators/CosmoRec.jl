# Chunk 5e acceptance: derivatives of the complete single pass

## Gates
- `test/chunk5e_history_ad.jl`: **37/37** (focused, twice: `focused5e_run5.log`, `focused5e_run6.log`).
- Full root `Pkg.test()` with `COSMOREC_NATIVE_DATA_DIR` set (all Phase 5 tests registered): **49444/49444**, `Testing CosmoRec tests passed`, exit 0, 29m41.7s, `chunk5/pkgtest5e.log` (49407 + 37).
- Errors: ForwardDiff vs finite differences (elasticity) 1.2e-5 best step; JVP 3.3e-8; prepared Mooncake VJP through the Rodas5P steps vs ForwardDiff 1.3e-9 (random projection) to 1.3e-7 (pass Xe block). Benchmarks: primal pipeline 715 ms, ForwardDiff Jacobian 1.39 s, Mooncake cold prepare 4.2 s, hot gradient 1.10 s.
## Accepted scope
Gradients of the complete single pass with respect to `[F, A2s1s, hscale, nbscale]` with the switch node, node grid, initial state and preliminary history frozen; forward mode through `solve`, reverse mode by Mooncake directly through the solver steps.
## Limitations
Continuous adjoints are unstable for this system (documented failure); no event-location or grid derivatives; Yp/h100/T0/table derivatives not included (they enter through frozen preliminary parts); solver-tolerance-level agreement; machine-specific data path (D1).

## CORRECTION (2026-10-02, Phase 10 investigation)
The reverse-mode route was described above as "Mooncake directly through the Rodas5P steps (discretize-then-differentiate), no sensealg". That claim is NOT established. With no explicit `sensealg`, the installed SciMLSensitivity 7.119.12 (`src/concrete_solve.jl:252-258`) selects `ForwardDiffSensitivity()` whenever `length(u0) + length(p) <= 100`. The 5e problems (12 or 7 states, 2 or 5 parameters) fall in that branch. The accepted 5e gradients were therefore most likely computed by forward-mode sensitivities inside the solve rule, not by reverse differentiation of the steps. The numerical agreement gates of 5e are unaffected. The route description is wrong and is superseded by the Phase 10 reverse-mode investigation (`docs/STATUS_PHASES7_10.md`). For large parameter vectors the same default selects GaussAdjoint (`:367-389`), the continuous adjoint that failed in 5e.

Precise route (confirmed in `chunk10/probe_mc_sensealg.log`: the default and explicit ForwardDiffSensitivity give bitwise-identical gradients on a 2-parameter block): the 5e reverse-mode gradient is OUTER Mooncake (the prepared gradient of the pipeline) combined with INNER ForwardDiffSensitivity (forward-mode sensitivities computed inside the SciMLSensitivity solve rule). It is not Mooncake reverse through the Rodas5P steps. In the tested versions (LinearSolve 5.18.2, Mooncake 0.5.61, SciMLSensitivity 7.119.12, SciMLBase 3.57.0), true reverse through the steps fails in the LinearSolve cached `solve!(cache)` Mooncake rule (`chunk10/mre_linearsolve/`). This is a version-specific limitation of that rule, not a general limitation of Julia AD.
