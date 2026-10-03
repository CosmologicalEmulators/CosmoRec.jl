# Chunk 5e results: derivatives of the complete single pass

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al.; cite Chluba & Thomas 2010, Chluba & Sunyaev 2006, Rubino-Martin et al. 2008. No native source, library or table was changed or bundled; the tests require `COSMOREC_NATIVE_DATA_DIR` (D1).

## Parameter set and frozen parts
p = `[F, A2s1s, hscale, nbscale]`: the Recfast fudge factor and H 2s-1s rate (enter the Recfast tail), and multiplicative scales of H(z) and of the baryon density (NH in the pass, Omega_b in Recfast) applied consistently in the pass and the tail (`hscale`, `nbscale` keywords of `recombination_rhs!`, `ode_background`, `recombination_pass`, `recfast_tail_inputs`, `recombination_history`). FROZEN at their primal values (no derivative claimed): the node grid, the helium-switch node (1342; discrete, no event-location derivative), the initial Saha state and the 4b/4c preliminary history (so Yp, h100, T0 and the rate tables are not differentiated). Structural check: F and A2s1s cancel in the rescaled Recfast `xp` equation (`ff ∝ 1/F`) and do not act above z = 50 (tested).
`recombination_pass(...; k_switch, nodes)` solves/reports only selected nodes (plus the switch and end nodes) for the frozen branch; it reproduces the full-grid pass at those nodes to 5e-8.

## Routes and results (focused 37/37)
- ForwardDiff through the complete pipeline (pass + switch + 7-state + tail + spline assembly to 9 grid redshifts; Rodas5P, explicit Jacobian/time gradient, primal error control, reltol 1e-12) vs central finite differences (h = 1e-2, 1e-3, 1e-4), measured as an ELASTICITY error `|J - J_fd| |p| / |f|` (column-scaled errors are meaningless here because the F column is ~0): best 1.2e-5, worst step 1.2e-4 (finite differences are limited by the solver tolerance); gate 1e-4. JVP along a random direction vs symmetric differences: 3.3e-8 (best step).
- Prepared Mooncake VJP (two independent preparations, two parameter points, four weight blocks) DIRECTLY THROUGH THE Rodas5P STEPS (discretize-then-differentiate: `ODEProblem{true, FullSpecialize}`, plain-function RHS with the model as a global constant, explicit ForwardDiff Jacobian and time gradient, no `sensealg`), on the selected-node pipeline, vs the ForwardDiff projection `J'w`: random projection 1.3e-9, pass Te block 2.1e-9, tail block 2.6e-9, pass Xe block 1.3e-7 (Mooncake route at reltol 1e-10 vs ForwardDiff at 1e-12); `gA == gB`. Gate 1e-6.
- Discreteness check: a 1e-3 change of H(z) did not move the switch node (0 nodes); the switch is still a discrete operator (Chunk 4d).

## Failed approaches (logs kept in `chunk5/`)
1. Column-scaled FD comparison: errors 5e6 (the F column is ~0) -> replaced by the elasticity metric (`focused5e_run1.log`).
2. Continuous adjoints (GaussAdjoint / QuadratureAdjoint with MooncakeVJP, the Chunk 4c route): first `RateTableDomainError(Te = 0)` because SciMLSensitivity probes the RHS with a zero state (fixed by a test-only rho <= 0 guard), then gradients wrong by O(1): the backward adjoint solves abort with `dt` below epsilon and states growing to 1e59 (adjoint Jacobian entries up to 3e21; `diag5e.log`, `focused5e_run2-4.log`). Not accepted.
3. Direct Mooncake through Rodas5P with explicit Jacobian/tgrad matched ForwardDiff bitwise on a single 12-state solve (`diag5e2.log`) and is the accepted reverse route (`focused5e_run5/6.log`, 37/37 twice).

## Benchmarks (`benchmark/chunk5e_benchmarks.jl`, `chunk5/bench5e.log`)
Primal complete pipeline 715 ms (884 MiB); ForwardDiff Jacobian of the complete pipeline (18 x 4) 1.39 s; selected-node primal 673 ms, its ForwardDiff Jacobian 1.32 s; Mooncake `prepare_gradient` cold 4.23 s (one sample, includes compilation); prepared hot gradient 1.10 s (1.35 GiB).

## Limitations
No derivative of the switch location, node grid, solver step sequence (the Mooncake route differentiates the realized steps), initial state or preliminary history; the reverse route needs the model as a module-level constant and a plain-function RHS; gradients agree to the solver tolerance level only; one cosmology; machine-specific data path (D1).

## Files
`src/RecombinationODE.jl` (scales, node subsets), `src/CosmoRec.jl`, `test/chunk5e_helpers.jl`, `test/chunk5e_history_ad.jl`, `test/runtests.jl`, `benchmark/chunk5e_benchmarks.jl`, this file, `docs/ACCEPTANCE_CHUNK5E.md`. No Git operation.

## CORRECTION (2026-10-02, Phase 10 investigation)
The reverse-mode route was described above as "Mooncake directly through the Rodas5P steps (discretize-then-differentiate), no sensealg". That claim is NOT established. With no explicit `sensealg`, the installed SciMLSensitivity 7.119.12 (`src/concrete_solve.jl:252-258`) selects `ForwardDiffSensitivity()` whenever `length(u0) + length(p) <= 100`. The 5e problems (12 or 7 states, 2 or 5 parameters) fall in that branch. The accepted 5e gradients were therefore most likely computed by forward-mode sensitivities inside the solve rule, not by reverse differentiation of the steps. The numerical agreement gates of 5e are unaffected. The route description is wrong and is superseded by the Phase 10 reverse-mode investigation (`docs/STATUS_PHASES7_10.md`). For large parameter vectors the same default selects GaussAdjoint (`:367-389`), the continuous adjoint that failed in 5e.

Precise route (confirmed in `chunk10/probe_mc_sensealg.log`: the default and explicit ForwardDiffSensitivity give bitwise-identical gradients on a 2-parameter block): the 5e reverse-mode gradient is OUTER Mooncake (the prepared gradient of the pipeline) combined with INNER ForwardDiffSensitivity (forward-mode sensitivities computed inside the SciMLSensitivity solve rule). It is not Mooncake reverse through the Rodas5P steps. In the tested versions (LinearSolve 5.18.2, Mooncake 0.5.61, SciMLSensitivity 7.119.12, SciMLBase 3.57.0), true reverse through the steps fails in the LinearSolve cached `solve!(cache)` Mooncake rule (`chunk10/mre_linearsolve/`). This is a version-specific limitation of that rule, not a general limitation of Julia AD.
