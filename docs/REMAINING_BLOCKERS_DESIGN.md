# Remaining blockers: design options (approval-dependent; NOTHING here is implemented or accepted)

**Partial approval, 2026-10-02 08:34:** Marco approved changing ONLY the final-history maximum pointwise relative Xe error gate from 5e-7 to **1e-5 (10 ppm)**. The last measured 1.71e-6 is within this gate. This does not approve changing population/PDE/gradient gates, the numerical formulation, solver dependencies, or background-model policy. The 5e-7 discussions below describe the historical stricter Xe requirement; the other production blockers and proposed implementation choices remain unresolved.

Status at handoff (2026-10-02, about 07:00): the latest full suite (`chunk10/pkgtest7_10_run2.log`) gave 50809 passed, 9 failed, 0 errored, 0 broken (rc 1, 47m40s). Phases 7-10 are NOT accepted. The supervisor's independent stage reruns (7a-9a, local AD) are scoped kernel validation only. All gates are unchanged. See `docs/STATUS_PHASES7_10.md` for the evidence. This document proposes routes for the three remaining blocker classes. Each route lists its source basis, an MRE/repro, the acceptance tests it would need, its costs, and a ready-to-use implementation prompt. Every route needs Marco's decision before any numerical formulation, gate or scope changes.

Terminology used throughout:
- **Frozen 128-bit diagnostic map**: the Julia composition with the ODE step sequences (109 blocks) and the Patterson stopping levels recorded at p0 and replayed, evaluated in BigFloat(128) (`chunk10/diag10b_big*.log`). It is smooth by construction and is NOT the production map.
- **Adaptive public map**: `recombination_history_diffusion` with adaptive Rodas5P and native-rule quadrature decisions, in Float64. This is what users call.
- **Native baseline**: the ORIGINAL `cosmorec_calc_h_cpp_(runmode = 0)` (Gear/BDF ODE solver, `Development/ODE_PDE_Solver/ODE_solver_Rec.cpp`, Float64), fixture `native_cosmorec_runmode0.txt`.

---------------------------------------------------------------------------------------------------------------------------------------
## B1. Precision-stable runmode-0 parity (10a final Xe 1.07e-6 / 1.71e-6 vs gate 5e-7; pass-1 populations; DF feedback range)

### Evidence (what is established, and what is not)
- Established, conditioning of the composition:
  - 128-bit vs Float64 primal on the FROZEN map: 4.4e-6.
  - Codegen alone (`--check-bounds=yes`): final Xe moves by about 6e-7, reproduced bitwise in a focused run.
  - Smooth native-tolerance perturbations of the pass-0 populations: final Xe moves by 1.1-1.2e-6.
  - Dnem cancellation factor kappa up to 1.9e8; PDE outputs with 59-90 of 199 outputs containing unconverged quadrature sub-integrals.
- Established, wiring: with native PDE outputs injected, Julia reproduces the native final Xe to 1.9e-7. The ODE, feedback, tail and assembly are therefore consistent with native.
- NOT established:
  - how far the ADAPTIVE public map is from its own exact-arithmetic value;
  - how far the NATIVE baseline is from its exact-arithmetic value. The native code is Float64 Gear/BDF, with its own round-off and its own discretization;
  - hence whether a 5e-7 Julia-vs-native agreement is achievable by ANY Float64 implementation.

### Step 0 (diagnostic only, no formulation change): compare the REAL public map at higher precision
Goal: measure, without changing any numerical formulation, the precision floor of (a) the adaptive public Julia map and (b) the native map.
1. Julia adaptive public map in BigFloat(128): run `recombination_history_diffusion` with BigFloat state, parameters and tolerances. The public code is type-generic; Rodas5P, the splines, the PDE and Patterson all accept BigFloat. Keep the tolerances numerically identical (reltol 1e-12 etc.), so the step selection follows the same rules but is driven by higher-precision error estimates.
   - Compare BigFloat-public vs Float64-public: this is the Float64 round-off of the public map, not of the frozen map.
   - Compare BigFloat-public at reltol 1e-12 vs 1e-14: this is the discretization floor.
   - Cost estimate from the frozen run: about 320 s per composed evaluation at 128 bit with recorded steps. An adaptive 128-bit run may take 2-5x more steps, so about 10-30 min per run and about 2-3 GiB. Needed runs: 2 tolerances x 1 point; save Xe/Te on the 10000-node CAMB grid as text.
2. Native precision floor, with read-only native code and harness-only changes:
   - (a) Repeat `harness10a` with different but equally valid native settings, e.g. ODE tolerances x0.1 through the harness's parameter object, if exposed; otherwise document that it is not exposed. Measure the native self-spread.
   - (b) Compiler-flag spread: build the harness against the existing libCosmoRec.a (cannot change) vs a SEPARATE local rebuild of the ORIGINAL sources in a scratch copy outside the native tree, with `-O0` vs `-O2` vs `-ffloat-store`. Measure the native codegen spread. The scratch-copy rebuild must leave the native tree untouched; Marco to approve.
3. Acceptance of Step 0 (diagnostic): the report states three numbers: e_J = |public64 - public128|, e_N = native self-spread, and d = |public128 - native|.
   - If d <= 5e-7 while e_J and e_N exceed it, the 5e-7 gate measures Float64 noise, and Marco decides whether it stays.
   - If d > 5e-7 with small e_J, there is a genuine modelling/solver difference to find. Next candidates are the Gear/BDF vs Rodas5P discretization and the DF z-grid.

### Options after Step 0 (each needs Marco's approval)
- O1 (no formulation change; scope of the gate): relate the runmode-0 gate to the measured native self-spread e_N, e.g. |Julia - native| <= max(5e-7, 2 e_N). This changes acceptance and is NOT to be done without Marco.
- O2 (higher-precision kernels, formulation-preserving): compute the hypersensitive pieces in Double64 (DoubleFloats.jl) or BigFloat inside Float64 interfaces, namely the Dnem/DnLj differences (kappa ~ 1e8), the PDE march, the integrals, and then round to Float64 at the stage boundary.
  - Source basis: the cancellations are at `HI_pd_Rp_splines_effective.cpp` (Dnem) and `Solve_PDEs_integrals.cpp` (DnLj, DRtot, the sign-changing integrand).
  - Acceptance: 10a unchanged gates; the stage outputs move toward the 128-bit public map by >= 10x; benchmarks reported.
  - Cost: Double64 is about 5-20x slower in the PDE stage (73 ms now, so about 0.4-1.5 s per stage); with 2 stages per run, +1-3 s per runmode-0 call.
  - Risk: it does not remove the native's own noise. If e_N > 5e-7, agreement with native can still fail.
- O3 (native-solver replication): port the native Gear/BDF solver (orders 1-5, Newton with BiCG/Gaussian elimination, the native step control) as an optional solver callback, so the passes reproduce the native discretization rather than Rodas5P.
  - Source: `ODE_solver_Rec.cpp`, about 1500 lines.
  - Acceptance: pass-0 nodewise vs native within the native tolerance units (<= 1 for X1s/rho) with a smaller spread than Rodas5P; 10a gates.
  - Cost: large (a chunk of its own); AD through a hand-written BDF needs its own tests; Mooncake reverse would then use non-cached `\`, which also helps B2.
- O4 (quadrature): replace the Patterson stopping decisions by fixed high-order rules. This CHANGES the native algorithm, so it is a formulation change; only with Marco. It removes the unconverged-piece variability but departs from native parity.

### Next implementation prompt (Step 0, diagnostic only)
"Phase 10 Step 0 (diagnostic, no formulation or gate change): evaluate the ADAPTIVE public `recombination_history_diffusion` in BigFloat(128) at reltol 1e-12 and 1e-14 (same rules), save Xe/Te on GRID5 as text, and report e_J = |Float64 - 128| and the discretization spread. Separately, measure the native self-spread with harness-only variations (no native-tree modification; any rebuild only in a scratch copy, and only with approval). Report e_J, e_N and d; do not change any test gate."

---------------------------------------------------------------------------------------------------------------------------------------
## B2. True reverse mode through the production ODE solver (Rodas5P)

### Evidence
- Installed versions: LinearSolve 5.18.2, Mooncake 0.5.61, SciMLSensitivity 7.119.12, SciMLBase 3.57.0.
- `ext/LinearSolveMooncakeExt.jl:117-118` errors in the reverse rule of `solve!(cache)` for every consumption pattern of the 2x2 MRE: the returned `sol.u` (the pattern the error message recommends), `cache.u`, a zero-cotangent primal decision, and explicit LUFactorization. Non-cached `solve(LinearProblem(A, b))` gives the exact gradient. Files: `chunk10/mre_linearsolve/mre.jl`, `mre.log`.
- The supervisor showed `DI.prepare_pullback` also fails, because Mooncake 0.5.61 `interface.jl:669-670` resets the pullback with zero rdata (`mre_pullback.jl`, `supervisor_pullback.log`). The comment in LinearSolveMooncakeExt.jl:37-41 is stale for this version.
- Consequences:
  - `MooncakeAdjoint()` and `SensitivityADPassThrough()` error;
  - the default sensealg picks ForwardDiffSensitivity for <= 100 parameters (this is what 5e used) and GaussAdjoint above (zero/unstable).
- Working today: outer Mooncake + inner ForwardDiffSensitivity. Through a checkpointed chain on the full composition it is at 1.9e-6 vs ForwardDiff (FAIL at 1e-6) and costs 2548 s, 15.8 GiB. The cost of ForwardDiffSensitivity scales with length(u0) + length(p), here 7-12 + 802.

### Isolated validation workflow (no shared-depot change, no issue/PR submission)
1. Create a separate, throwaway Julia environment with its OWN depot (`JULIA_DEPOT_PATH=<scratch>/depot:<shared, read-only>`), e.g. under `local_analysis/.../chunk10/linsolve_isolated/`. `Pkg.develop` a LOCAL CLONE of LinearSolve into the scratch directory and pin Mooncake 0.5.61. The shared `~/.julia/packages` is never edited.
2. Reproduce `mre.jl` in that environment unchanged; it must fail identically (the baseline).
3. Diagnose in the clone only, with candidate root-cause checks:
   - (a) whether `sol.dx.data.u` actually receives the cotangent. Instrument the rule in the clone to print the norms of `∂u`, `λ` and `tu` for MRE case 1. Because the `sol` returned by `zero_fcodual(solve!(cache))` aliases `cache.u`, the cotangent may accumulate on the CACHE's fdata (`fdata(_cache.dx).fields.u`) instead of `sol.dx`. In that case `∂u` read from `sol.dx` is zero and line 117 fires.
   - (b) compare with the Mooncake 0.5.61 rules API: whether the rule must read the cotangent from both `sol.dx` and `_cache.dx.fields.u`.
4. Candidate fix (in the clone only): read `∂u = sol.dx.data.u .+ fdata(_cache.dx).fields.u` (and zero the latter), and remove the error for a legitimately zero cotangent, i.e. return zero ∂A, ∂b.
   - Validate on the MRE (cases 1-5, gradients vs ForwardDiff and analytic) and on a 2-parameter Rodas5P block with `MooncakeAdjoint()` vs ForwardDiff (dot test 1e-8).
   - Then validate the packed-feedback block and the full composition (1e-6 gate).
5. Upstream path: prepare a patch file and an MRE text locally (`chunk10/linsolve_isolated/PATCH.diff`, `ISSUE_DRAFT.md`). Do NOT submit. Marco decides on submission.
6. Alternative without touching LinearSolve:
   - a test-only Rosenbrock/Rodas integrator written with non-cached `\` (Mooncake-native), replaying the recorded Rodas5P step sequence;
   - or O3's native Gear/BDF port.
   Either changes the solver route used for gradients and needs approval. Acceptance: its primal matches the public Rodas5P map to <= 1e-10 on the recorded steps, and Mooncake vs ForwardDiff passes 1e-8 per block.
- Costs:
  - Isolated environment and diagnosis: about 1-2 h, plus 10-20 min precompile per environment.
  - Each full-composition reverse validation: about 45 min and 16 GiB (one job at a time).
  - Cost of a working true-reverse route: UNVERIFIED. The often-quoted 2-5x one primal (about 10-60 s here) is an UNMEASURED heuristic, not an expected or guaranteed cost; Mooncake tape size and cached-solve overheads on this problem are unknown. Any acceptance must report measured cold `prepare_gradient`, first gradient and prepared hot gradient (BenchmarkTools, evals = 1) and peak RSS, against today's measured hybrid baseline of 2548 s and 15.8 GiB.

### Next implementation prompt
"B2 isolated validation: in a scratch environment with its own depot, develop a local clone of LinearSolve 5.18.2 with Mooncake 0.5.61 pinned. Reproduce chunk10/mre_linearsolve/mre.jl, instrument the solve!(cache) rrule in the CLONE ONLY to locate where the cotangent of sol.u lands, and test a candidate fix on MRE cases 1-5, then on a 2-parameter Rodas5P block with MooncakeAdjoint vs ForwardDiff. Do not modify the shared depot, do not submit issues or PRs; write PATCH.diff and ISSUE_DRAFT.md locally."

---------------------------------------------------------------------------------------------------------------------------------------
## B3. Remaining work toward the ORIGINAL goal: whole cosmological initialization and background derivatives

This is NOT an optional scope expansion. The original goal is a fully differentiable CosmoRec, and the current derivatives cover only four parameters: two Recfast-tail parameters and two uniform background proxies. They are therefore INCOMPLETE with respect to that goal, independently of the B1/B2 blockers. The current result must not be described as cosmological sensitivities.

### What the current derivatives cover
p = [F, A2s1s (Recfast tail only), hscale, nbscale]. hscale and nbscale multiply H(z) and N_H(z) uniformly. They are proxies, not cosmological parameters.

### What is frozen (no derivative path today)
- **Preliminary Recfast history**: the Cosmos accessors (`CosmosAccessors.jl`) are built from the NATIVE 6000-node history fixture (`native_recfast_history.txt`): Xe_Seager, rho = TM/TCMB, Xe_H, Xe_He splines. They enter the Saha initial state, the output fallback outside the stored range, and the accessors used in the passes.
- **Saha initial state** at zstart (`SahaInit.jl`, fed by the accessors).
- **Loaded Hubble table** (the CAMB H(z) input, `hubble_table`), including the reproduced native low-z bug.
- **Cosmology closure constants**: `CosmosConstants` (Y_p, T_CMB0, Nb0, H0, Omega_k, Omega_L, Omega_m, Omega_rel, fac_mHemH, zsRe) are native-printed inputs, not derived from (omega_b, omega_c, h, N_nu, T_CMB, Y_He).
- **Discrete structure**: the helium switch node (k_switch = 1342, z = 1680.91), the node grids, the PDE frequency and coefficient grids, the Patterson decisions, the polint_JC stencils, the rate-table stencils, the feedback range ends (500, 2000), and the Recfast tail grid.

### Design for full cosmological sensitivities
1. Parameter vector theta = (omega_b, omega_c, h, N_nu, T_CMB, Y_He, plus F, A2s1s) as the public input.
2. Port the cosmology closure (`Development/Cosmology/Cosmos.cpp`: H(z) without a loaded table, N_H, fHe, Omega_rel from N_nu and T_CMB), with native-parity fixtures from harness4b-style probes at 3-5 cosmologies (native read-only).
3. Make the preliminary history differentiable: run the Recfast++ system (already ported, `RecfastPP.jl`, chunk 4c) inside the pipeline instead of loading the native history. Native-parity test vs `native_recfast_history.txt` at the fiducial point, then vs harness runs at other cosmologies.
4. H(z) POLICY (a decision for Marco; never swap models silently). In production, H(z) is an EXTERNALLY SUPPLIED CAMB table H(z; theta), which the native code splines (including its low-z node-0 bug). The real choice is HOW that table is differentiated:
   - (a) keep the external CAMB table as the model and take dH/dtheta from a differentiable source. A differentiable background solver, or CAMB-side derivatives supplied as a table, are passed in with the table, and the spline is differentiated w.r.t. its node values (already AD-transparent). Native parity stays exact at the fiducial point.
   - (b) the analytic Cosmos closure (`Cosmos.cpp` H(z) without a loaded table). This is a DIFFERENT model from production: its H(z) differs from the CAMB table (dark energy, massive neutrinos, numerical details). Using it changes the primal and needs explicit approval, with its own parity tests vs a native run without a loaded table.
   Option (a) is the production-faithful route. Option (b) is not a drop-in replacement.
5. Discrete structure: keep it frozen at the nominal values, with explicit checks that it is unchanged under the perturbation (as for k_switch), and report the derivative as "frozen-branch".
- Acceptance: native parity per cosmology (existing gates); derivative tests vs 128-bit FD on the frozen map and vs ForwardDiff for the reverse route (1e-6); a statement of every frozen dependency.
- Cost: cosmology closure and fixtures about 1 chunk; the Recfast++-in-pipeline about 1 chunk plus a re-run of the 5b-10a fixtures at the fiducial point (expected unchanged); per-cosmology native harness runs about 1 s each.

### Next implementation prompt
"B3 (after Marco's H(z) policy decision): port the native Cosmos closure constants (N_H, fHe, Omega_rel from N_nu/T_CMB; H(z) per the DECIDED policy: (a) the external CAMB table with supplied dH/dtheta node derivatives, or (b) the analytic closure as an approved model change with its own native parity run without a loaded table). Use native fixtures at 3-5 cosmologies (read-only harness). Then replace the loaded native preliminary history by the in-pipeline Recfast++ solve (RecfastPP.jl) and verify parity at the fiducial point, without changing any production gate. Keep all discrete structure frozen; report derivatives w.r.t. (omega_b, omega_c, h, N_nu, T_CMB, Y_He) as frozen-branch, tested vs 128-bit FD on the frozen map, with cold/prepared/hot benchmarks."

---------------------------------------------------------------------------------------------------------------------------------------
## Summary of decisions needed from Marco
1. B1: approve Step 0 (diagnostic only). After it: O1 (gate relative to the measured native spread), O2 (Double64 kernels), O3 (native Gear/BDF port) or O4 (quadrature change, a departure from native).
2. B2: approve the isolated LinearSolve clone validation; decide later on an upstream submission; or approve a test-only non-cached integrator route.
3. B3 (remaining work toward the original goal, not optional): approve the plan for true cosmological parameters (closure constants, in-pipeline preliminary history), and decide the H(z) policy: differentiate the external CAMB table, option (a), or switch to the analytic closure, option (b), which is a model change.
