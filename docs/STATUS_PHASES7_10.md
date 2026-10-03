# SCIENTIFIC FORWARD VALIDATION (Boltzmann-facing, 2026-10-03) — separate from B2/B3 and from the test-gate status below

Details are in `docs/BOLTZMANN_FORWARD_VALIDATION.md`.

**Method.** Native CosmoRec v3.0b and the Julia port (injected Rodas5P test callback 1e-12/1e-18/1e-14; no library default) receive identical inputs. Their complete X_e/T_m histories are injected into one CAMB 2.0.4 build through a diagnostic scratch hook; the original CAMB, CosmoRec and env are untouched. Points: the fiducial plus ω_b + 0.001, ω_c + 0.005, H0 = 70 and Y_He = 0.24, each with native-provenance initialization from parameterized native harnesses that are bitwise-validated at the fiducial.

**History criteria (unchanged).**
- X_e ≤ 1e-5: PASS at 5/5 (1.45–2.88e-6).
- T_m < 1e-7: PASS at 4/5; **FAIL at ω_b + 0.001 (1.358e-7)**.

**C_ℓ** (ℓ 2–3000; unlensed TT/EE/TE, lensed TT/EE/BB/TE, φφ):
- port error max fractional ≤ 1.27e-6 (TE correlation-normalized ≤ 5.7e-7);
- ≤ 6.7e-5 σ_CV at any ℓ (full-sky cosmic variance; diagonal per-spectrum sums ≤ 2.4e-6, not a joint χ²);
- the native's own CAMB-path floor is 0.1–1e-6, and the authors' 0.1% scale is about 800× larger.

**Not demonstrated:**
- Julia-only cosmological initialization (B3);
- CLASS or other settings;
- parameter-level inference;
- gradients/true reverse (B2).

Internal pass-1 and PDE parity are reported separately. Pass-1 X1s fails its gate at all points, as at the fiducial. `oc` has a single-node 2s/2p outlier at the z = 500 feedback cut.

# UPDATE 2026-10-02 afternoon (a1 fix, Marco-approved) — read first; supersedes the counts below

**What changed.** The runmode-0 test/reference callback configuration used for these native comparisons went from (reltol, a1, aex) = (1e-12, 1e-16, 1e-14) to (1e-12, **1e-18**, 1e-14) in every injected test callback carrying that triple: `SOLVE10`, `TOL5E`, `SOLVE5B/5C/5D` and the 6a end-to-end pass. Reverse-route, tail and benchmark configurations are unchanged. No physics, native, dependency or gate change, and no library default: the public functions take a required caller-supplied solve callback, and `Project.toml` has no solver dependency. Downstream callers still choose their own solver and tolerances, and the benchmark costs refer to the injected test callbacks. Details and the full table are in `docs/A1_FIX_VALIDATION.md`. The attribution evidence is in `docs/ORIGINAL_VS_PORT_ACCURACY_INVESTIGATION.md` §5.6.

**Latest full suite** (`chunk12/pkgtest_a1_run1.log`; unconditional root Pkg.test, `--check-bounds=yes`): **50812 passed, 6 failed, 0 errored, 0 broken; rc 1; 55m51s; peak 16.4 GiB. NOT green; Phases 7–10 NOT accepted.**

**Resolved:**
- 10a:104 PDE0/PDE1, by the a1 fix: 1.075e-3 / 1.039e-3 → 3.27e-5 / 5.12e-5 (≤ 1e-3).
- 10a:118 final Xe, by Marco's approved 1e-5 gate; the a1 fix also lowered it from 1.71e-6 to 1.45e-6.

**Still failing:**
| gate | old | new | bound | note |
|---|---|---|---|---|
| 10a:83 | 1.655e-4 | 1.655e-4 (unchanged) | ≤ 1e-4 | native-injected PDE1 2g3s, scaled metric |
| 10a:119, pass-1 X1s | 3.25 units | 5.21 units | < 1 | |
| 10a:120, pass-1 2s/2p | 15.4 / 18.3 | 15.95 / 16.18 units | < 10 | |
| 10b:42 | 2.37e-3 | 3.04e-3 | < 1e-4 | Float64 FD reference noise |
| 10b:56 | 3.51e-5 | 2.72e-5 | < 1e-5 | Float64 FD reference noise |
| 10d:46 | 1.00093549296115e-6 | same (bitwise) | < 1e-6 | |

The original native code does not meet the pass-1 gates against itself (investigation §5.9).

**Hybrid reverse probes** (outer Mooncake + inner ForwardDiffSensitivity; NOT true reverse). With their own tolerances they are unchanged by the edit. Labelled diagnostics with the reverse-route a1 also set to 1e-18 shrink the mismatch without passing:
- reduced: 1.47e-5 → 6.95e-6;
- full: 1.9e-6 → 1.3599712114113554e-6 (gate 1e-6: FAIL).

**Cost** (matched BenchmarkTools, default flags):
- primal pass, history and runmode-0: 1.31–1.35×;
- ForwardDiff Jacobians (5e / 10b): 1.23× / 1.20×;
- memory: about 1.3×.

**Not fixed by a1:**
- the native-baseline pass-1 precision (gate policy is Marco's decision);
- B2: true reverse through Rodas5P, blocked by the LinearSolve 5.18.2 / Mooncake 0.5.61 cached-solve rule;
- B3: whole cosmological initialization, unimplemented, with the H(z) policy pending;
- the 10b FD reference design;
- the 10d Float64 AD agreement.

**Remaining consistency item:** the standalone benchmark scripts still configure a1 = 1e-16, deliberately unchanged. The matched old/new benchmark of the affected routes is already done (`docs/A1_FIX_VALIDATION.md` §4); aligning the scripts is Marco's decision.

**Task status:** the approved a1 validation task is complete, with a PARTIAL fix. Not all issues are fixed.

**Gradient-gate review** (`docs/FINITE_DIFFERENCE_ACCURACY_REVIEW.md`, read-only; no gate changed):
- 10b:42/56 are dominated by the FD reference map. Central differences sample across clamp-induced non-smooth residual features of the parameter-response map (`eps_A_effective`, `src/RateTable.jl:165`). The clamp-on/off FD test-reference effect is 3.04e-3 → 1.44e-4 (elasticity residual; counterfactual, still > 1e-4), plus a roughness of about 4e-7. This is not evidence of a missing gradient.
- 10d: a 128-bit exact-map reference puts both Mooncake and ForwardDiff within 1.05e-6 of exact (at the floor; no AD-bug signature).
- The non-smooth parameter-response residuals (≈ 3e-5 neighbour increments of the cubic-detrended relative response; native and Julia scans are different maps) are NOT an Xe error and not comparable with the fiducial cross-code Xe gate (1.45e-6, PASS). Their C_ℓ impact is unmeasured; they are a possible sensitivity risk for gradient use.

# MORNING SUMMARY (2026-10-02, updated ~07:15) — historical; counts superseded by the update above
**User-approved update (2026-10-02, 08:34):** the final-history maximum pointwise relative electron-fraction error gate is now **1e-5 (10 ppm)**, replacing 5e-7. The previously measured 1.71e-6 is within this gate. Only this Xe gate changed; population, PDE-correction, temperature, and gradient gates are unchanged. The full-suite counts below were obtained BEFORE this change; a fresh full result is pending. Historical comparisons against 5e-7 below are retained as evidence, not the current Xe acceptance policy.

**Focused verification after that change:** `chunk10/focused10a_xe1e5_checkbounds.log`, run with `--check-bounds=yes`, reports **38 passed, 5 failed, 0 errored, 0 broken**. The final Xe assertion now passes at **1.7057284588522985e-6**. The remaining failures are the injected-population PDE correction (line 82), two full-path PDE correction checks (line 103), and the pass-1 population checks (lines 118–119). No new full-suite count is claimed.

**Latest full suite** (`chunk10/pkgtest7_10_run2.log`): **50809 passed, 9 failed, 0 errored, 0 broken; rc 1; 47m40s. NOT green; Phases 7-10 NOT accepted.**
**Full hybrid reverse result**: outer Mooncake (checkpointed chain) + inner ForwardDiffSensitivity on the full 3-ODE / 2-PDE / tail objective is at 1.9e-6 vs ForwardDiff, against a 1e-6 gate: FAIL. It is not true reverse through Rodas5P.
Design options for the remaining blockers: `docs/REMAINING_BLOCKERS_DESIGN.md`.
**Independent supervisor rerun** of the stage-level native tests (`chunk10/supervisor_chunk{7,8,9}*.log`, rc 0): 7a 99/99, 7b 115/115, 7c 47/47, 8a 35/35, 9a 161/161. The local 7-10 AD (20/20) was also rerun independently earlier. This is SCOPED kernel/oracle validation of the individual stages under the existing fixed bounds. It is NOT production acceptance. The caveats remain: conditioning (kappa up to 1.9e8; Float64 vs 128-bit 4.4e-6 on the frozen composed map), the 8a bound is a stopping-rule ESTIMATE and 59-90 of 199 outputs contain unconverged sub-integrals, and all 10a/10b/10d blockers stand. The global result is the fully registered suite: 50809 passed / 9 failed.

- **Accepted:** up to Phase 6a (full suite 50275/50275 at that point).
- **Implemented, passing their focused native tests:** 7a, 7b, 7c, 8a, 9a, and the local 7-10 AD (20/20). All are registered in runtests.jl.
- **NOT accepted (blocking gates, unchanged, failing):**
  1. Remaining production parity: pass-1 populations outside their native-tolerance gates; the native-injected PASS1 2g 3s at 1.66e-4 > 1e-4; the full-path DF feedback range at 1.08e-3 > 1e-3 under `--check-bounds=yes`. Final Xe 1.07e-6 (default flags) / 1.71e-6 (`--check-bounds=yes`) is within the user-approved 1e-5 gate.
  2. Composed derivative (10b): ForwardDiff vs FD fails (FD is noise-limited). 128-bit FD on the frozen map: Float64 ForwardDiff within 5.6e-6 (hscale) / 4.9e-6 (nbscale).
  3. Composed reverse mode: the full 3-ODE / 2-PDE / tail objective via checkpointed Mooncake (inner ForwardDiffSensitivity) is at 1.9e-6 vs ForwardDiff (gate 1e-6), FAIL. 10d background-link regression 15/16 (1.0009e-6 at one point).
  4. True reverse through the Rodas5P steps is BLOCKED by the cached `solve!(cache)` rule of LinearSolve 5.18.2 / Mooncake 0.5.61 (MRE in chunk10/mre_linearsolve). This is version-specific, not a general Julia-AD limit.
- **Defects found and fixed tonight:**
  - erased background derivatives (a value/type fast path in hi_diffusion_stage);
  - Float(Dual) promotion bugs;
  - the 5e route claim (it actually ran outer Mooncake + inner ForwardDiffSensitivity; corrected in the 5e docs).
- **Key evidence (conditioning):**
  - The composed Float64 output differs from 128-bit by 4.4e-6 on the frozen map.
  - Bounds-check codegen alone moves the final Xe by ~6e-7.
  - Native-tolerance population perturbations move the final Xe by ~1.1e-6.
  - The native C++ output is subject to the same conditioning.
- **Latest full suite** (`chunk10/pkgtest7_10_run2.log`, after the background fix, with 10d registered; 47m40s, rc 1): **50809 passed, 9 FAILED**. That is the 8 earlier failures (10a:82, 10a:103 x2, 10a:117-119, 10b:42, 10b:56) plus 10d:46. The 10a numbers are bitwise identical to run 1, as expected, because the fix leaves the primal unchanged. NOT green.
- **Decision for Marco:** whether the runmode-0 gates should be related to this measured double-precision conditioning, or whether a higher-precision or native-solver-replicating route is required. No gate was changed.

# Status of Phases 7-10 (HI radiation PDE, correction integrals, feedback, runmode-0 iteration) — NOT ACCEPTED

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al. No native source or data was changed. Nothing was staged, committed or pushed.

The production (runmode 0) path is **not accepted**. The full root suite was last green at Phase 6a (50275/50275, `chunk6/pkgtest6a.log`). All Phase 7-10 tests (7a, 7b, 7c, 8a, 9a, the local 7-10 AD, 10a, 10b) are now registered unconditionally. The first fully registered run, `chunk10/pkgtest7_10_run1.log` (48m55s, rc 1), gave **50794 passed, 8 FAILED**:
- 10a:82 (native-injected PASS1 2g 3s at 1.66e-4 > 1e-4);
- 10a:103 twice (full-path 2g 3s/3d in the feedback range at 1.08e-3 / 1.04e-3 > 1e-3);
- 10a:117-119 (final Xe 1.71e-6 > 5e-7; pass-1 X1s 3.2 units > 1; 2s/2p 15.4 / 18.3 units > 10);
- 10b:42 (composed FD elasticity 2.4e-3 > 1e-4);
- 10b:56 (JVP 3.5e-5 > 1e-5).
10c (the composed Mooncake objective) is not registered yet; it is blocked as explained below. All logs are under `cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk{7,8,9,10}/`, and failed runs are preserved there.

## Implemented (src/)
- `HIPDEProfiles.jl` (7a): A_SH, HI_Transition_Data, matrix elements; two-photon/Raman profiles from the ORIGINAL `two-photon-data` tables, which are read in place and never bundled. The tests locate them next to COSMOREC_NATIVE_DATA_DIR and fail loudly if they are missing. Also the PDE frequency grid and the arm_PDE_solver ratios. Production configuration only (nShells 3, nS_2gamma 3, nS_Raman 2).
- `HIPDEDefine.jl` (7b): def_PDE_Lyn_and_2s1s.
- `HIPDESolver.jl` (7c): Lagrange O2 stencils, Step_PDE_O2t, polint_JC lower boundary, the production stepping loop.
- `HIPDEIntegrals.jl` (8a): DI1_2s, DF_2gamma 3s/3d, DF_Raman 2s. It also provides a diagnostic count of unconverged quadrature pieces and `PattersonLevels` record/replay.
- `HIDiffusionFeedback.jl` (9a): setup_DF_interpol_data / interpolate_DF, plus the corrections in `fcn_effective!` through `RecombinationModel.diffusion`.
- `RecombinationDiffusion.jl` (10a): `recombination_history_diffusion`, the 3-ODE / 2-PDE iteration. `ScaledBackground` applies the hscale/nbscale background to the PDE stage too.

## Native oracles (external harnesses, two fresh runs byte-identical each)
harness7a (setup and coefficients), harness7c (ORIGINAL PDE run plus a stepping-loop replica that is bitwise equal to it), harness8a (native pd/Dnem spline knots), harness9a (feedback; the first capture is REJECTED because it set the DI1 flag after the spline setup), harness10a (the ORIGINAL runmode-0 call plus a CosmoRec() replica whose xe/tb are bitwise equal to it).

## Focused results
- 7a 99/99: grid bitwise; atomic data <= 1.7e-14; profiles and ratios <= 6e-16.
- 7b 115/115: coefficients <= 1e-15 of each array's maximum. D is 4.8e-14 and Dnem_eff 2.9e-13 because of the native cancellation in the Dp/Dn correction.
- 7c 48/48 (fixed bounds): final spectrum 2.8e-15. The raw early-step core differences of 2.7e-9 are under the FIXED bound 4 eps kappa_max = 1.8e-7, with kappa_max = 2e8 asserted on the fixed native input history (measured 1.88e8). The Dnem projection is a diagnostic only.
- 8a (fixed bound 1e-4 of the residual scale): passes. The scale derivation and its caveats are below.
- 9a 161/161: interpolation 1.3e-16. The correction term is at 0.69 rounding units of the native F_on - F_off difference.
- Local AD on dynamic inputs (`test/chunk7_10_pde_ad.jl`, 20/20):
  - Directional vs 256-bit FD: 7b 1.6e-11, 7c 2.1e-11, 8a (3-step march + integrals) 1.5e-10, 9a 1.6e-15.
  - Prepared Mooncake VJP vs ForwardDiff: <= 2.2e-15.
  - THIS IS LOCAL STAGE AD, NOT A COMPOSED RUNMODE-0 GRADIENT.

## UNMET / BLOCKING
1. **10a full Julia path vs the ORIGINAL runmode-0 call:**
   - Final Xe 1.07e-6 against the gate 5e-7.
   - Pass-1 X1s at 3.6 native-tolerance units (gate 1); 2s/2p at 15/18 units (gate 10).
   - Te 3.8e-8 passes.
   - The gates are ordinary blocking assertions. The 10a test fails and is not registered.
2. **10a native-injected diagnostic, native populations of PASS1 into the Julia PDE:** the 2g 3s output is at 1.66e-4 of the residual scale, against the fixed 1e-4 bound. FAIL.
3. **10b composed runmode-0 derivative (ForwardDiff vs FD):** best elasticity error 2.5e-3 to 2.9e-3 against the 5e gate 1e-4, and JVP 8.5e-5 against 1e-5. FAIL. Not resolved by freezing the Patterson decisions, so quadrature stopping is NOT a proven cause. Also, ForwardDiff through the frozen-replay map differs from ForwardDiff through the free map (norm still to be recorded).
4. **No prepared-Mooncake gradient of the composed runmode-0 objective exists yet.**
5. **No public-path benchmarks have been run for 7-10, and there is no full Pkg.test including 7-10.**

## Evidence so far (diagnostics, not conclusions)
- Native PDE outputs injected into Julia: pass 1 matches native PASS1 at solver level (X1s 1.4e-8, Xe 2.0e-7), and the final Xe matches the original call to 1.9e-7. The ODE, feedback, tail and assembly wiring reproduces native when given native PDE outputs.
- Native populations injected into the Julia PDE (PASS0): DF errors scaled <= 2.5e-5 (rel_to_max 3.6e-7). With the Julia pass-0 populations, which differ from native at solver level (~5e-7 in excited states), the DF errors in the feedback range are 1.4e-4 to 3e-4 of the maximum. The largest DF differences (55% for 2g, 15% for Raman) are at z >= 2000, where the native interpolate_DF returns 0.
- Near-equilibrium cancellation:
  - The Dnem cancellation factor kappa is up to 1.9e8.
  - Smooth perturbations of the native pass-0 populations at the native solver tolerance (rtol 1e-7 for rho/X1s, 1e-6 for excited states) move the final runmode-0 Xe by 1.07-1.22e-6 (`chunk10/diag10a_cond.log`). Uncorrelated node noise gives an unphysical 6e-4 to 8e-4 (`diag10a_cond_run1_uncorrelated_noise.log`).
- Quadrature: 2g 3s/3d and Raman have unconverged sub-integrals (never meeting the native stopping rule) at 90/78/59 of 199 outputs. There the native returns the 255-point value without any error guarantee.
- The residual scale for the resonance-split integrals is w(sum_i |r_i| + |DRtot|). Derivation: each sub-integral i stops when |r_k - r_{k-1}| <= max(epsabs, 1e-5 |r_i|), with epsabs = 1e-5 |running sum| <= 1e-5 sum_j |r_j|. For converged pieces the stopping-rule error ESTIMATE is therefore <= 9e-5 sum_j |r_j| over the <= 9 pieces. This is an estimate, not a bound, and it does not apply to unconverged pieces. The raw and absolute errors are recorded alongside in every test.
- Composed derivative, stage by stage (`chunk10/diag10b_stages.log`):
  - pass0: Dual vs Float primal 1.3e-11, derivative vs FD 1.6e-6. Fine.
  - Stage-0 DF: Dual vs Float primal differs by 65% of the global maximum, i.e. the hypersensitive z >= 2000 outputs.
  - pass1 / final: the FD error grows as h decreases (final: 1.8e-3, 3.5e-3, 7e-2 for h = 1e-3, 1e-4, 1e-5), which indicates ~1e-6-level non-smoothness of the composed map in p. Its source is not proven yet. A first candidate was p-dependent adaptive ODE step sequences, so they were frozen too (`chunk10/diag10b_steps.log`). Each ODE block's accepted steps were recorded at p0 (109 blocks) and replayed non-adaptively through tstops, together with the frozen Patterson decisions. Results:
  - Frozen replay vs adaptive at p0: 6.4e-7 (the replay is not bitwise; to be understood).
  - ||J_adaptive - J_frozen||, elasticity units: 1.3e-6.
  - FD of the frozen map STILL does not converge: elasticity error 3.4e-3, 6.1e-3, 3.9e-3, 2.5e-2, 0.38 for h = 1e-2 ... 1e-6.

  So neither the quadrature decisions nor the step sequences explain it. The frozen Float64 map itself shows output irregularity of about 4e-7 relative in the final Xe. The leading hypothesis is Float64 round-off amplified by the hypersensitive PDE/integral stage (the 8a ulp-response is up to 9e-5 of the residual scale). A 128-bit BigFloat evaluation of the frozen map, with its central difference, is running to test this (`chunk10/diag10b_big.log`).

## Decisions that need Marco
- Any change of the runmode-0 acceptance gates (e.g. relating them to the native output's own conditioning) or of the native-parity scope.

## Update: extended-precision test of the round-off hypothesis (`chunk10/diag10b_big.log`)
This covers ONE direction (hscale) of the FROZEN diagnostic map, i.e. ODE steps and Patterson decisions recorded at p0 and replayed:
- The 128-bit BigFloat central differences converge (h = 1e-6 and 1e-8 agree to 10 digits).
- The Float64 ForwardDiff derivative agrees with that extended-precision derivative to 5.6e-6 of |f0|.
- The Float64 primal differs from the 128-bit evaluation of the same map by 4.4e-6 relative.

This SUPPORTS the hypothesis that double-precision round-off, amplified by the PDE/integral stage, sets an accuracy floor of ~4e-6 on the composed output. That floor is larger than the 5e-7 parity gate. It is NOT yet shown for nbscale, F or A2s1s, for the adaptive public map, or as an explanation of the Julia-native difference. The gates are unchanged. Next higher-precision runs will save the full output/FD vectors as text.

## Update: reverse mode (Mooncake) through the composition
- One full PDE stage, prepared Mooncake gradient: dot test vs ForwardDiff 3.3e-13; cold prepare 169 s, hot 1.7 s, 12.1 GiB RSS (`chunk10/probe_mc_stage.log`).
- REDUCED composition (2 ODE passes, 1 PDE stage, no tail): FAILED with VJP relative error 1.0 (`probe_mc_composed_reduced.log`; rc 0 is not a pass).
- Isolation (`probe_mc_links_run1_ii_zero_gradient.log`, `probe_mc_isolate.log`):
  - The feedback RHS alone passes (dot test <= 2e-16), and so does a 2-parameter ODE block (2.8e-9).
  - A packed-feedback ODE block gives a ZERO Mooncake gradient.
- Cause, from the source: the default sensealg is GaussAdjoint for more than 100 parameters. Explicit `MooncakeAdjoint()` vs the default is under test (`probe_mc_sensealg.log`).

## Update: reverse-mode route through the production ODE solver (Rodas5P), evidence-based
1. Sensitivity dispatch, verified in the installed sources: SciMLSensitivity 7.119.12, SciMLBase 3.57.0, `concrete_solve.jl:252-258, 367-389`. Without an explicit sensealg it selects ForwardDiffSensitivity for length(u0) + length(p) <= 100, otherwise GaussAdjoint. `chunk10/probe_mc_sensealg.log` confirms it: the default and explicit ForwardDiffSensitivity give identical numbers on a 2-parameter block. The 5e gradients were therefore forward-mode sensitivities inside the rule, NOT reverse through the steps; the 5e docs carry a CORRECTION.
2. Packed-feedback ODE block (3178 parameters), same probe:
   - default (GaussAdjoint): zero gradient, FAIL;
   - explicit ForwardDiffSensitivity: fixed feedback w.r.t. hscale/nbscale PASS (1.0e-9), but packed coefficients as inputs FAIL (dot test 1.8e-2), arbitration by FD running (`probe_mc_sensealg2.log`);
   - MooncakeAdjoint(): ERROR inside the LinearSolve rule.
3. LinearSolve 5.18.2 / Mooncake 0.5.61 MRE (`chunk10/mre_linearsolve/mre.jl`, `mre.log`). The rule for `solve!(cache)` (ext/LinearSolveMooncakeExt.jl:117-118) errors for EVERY consumption pattern: the returned `sol.u` (the recommended one), `cache.u`, a zero-cotangent primal decision, and explicit LUFactorization. Non-cached `solve(LinearProblem(A, b))` returns the exact gradient. OrdinaryDiffEq's Rosenbrock methods use cached LinearSolve, so Mooncake reverse through the Rodas5P steps is NOT available with the installed versions. The shared package installation is not patched.
4. Consequence: a prepared Mooncake gradient of the composed 3-ODE / 2-PDE objective with the production solver is BLOCKED in these versions, unless a supported route is found:
   - ForwardDiffSensitivity costs about length(u0) + length(p) dual partials per solve and its packed dot test is not yet verified;
   - a test-only integrator using non-cached linear solves would change the solver route and needs Marco's approval as a scope change.
   The Mooncake gradient through a full PDE stage works (3.3e-13).
5. `chunk10/probe_mc_sensealg2.log`:
   - SensitivityADPassThrough() errors in the same LinearSolve cached-solve rule, even with 2 parameters.
   - FD arbitration of the packed-coefficient direction: FD 3.39e-4 at h = 1e-4 and 3.41e-4 at 1e-5; FDS-in-Mooncake 3.32e-4 (about 2% off); direct ForwardDiff 3.90e-4 (about 15% off). NEITHER is verified in that direction. It perturbs the spline coefficients independently, the sensitivity is forced by a knot-wise non-smooth forcing, and error control is on the primal only.
   - The composed derivative in the natural parameter directions is checked separately: 10b ForwardDiff vs 128-bit FD on the frozen map gave 5.6e-6 for hscale.
   - Next: the REDUCED composition (2 ODE, 1 PDE, no tail) with explicit ForwardDiffSensitivity and the feedback passed as 796 DF node values, Mooncake vs ForwardDiff on the identical map, with a hard assertion (`probe_mc_composed_reduced_fds.jl`).
6. Preparation API: the LinearSolveMooncakeExt comment (lines 37-41) recommends pullback preparation instead of gradient preparation. The supervisor tested `DI.prepare_pullback` / `value_and_pullback` on the 2x2 sol.u MRE (`mre_linearsolve/mre_pullback.jl`, `supervisor_pullback.log`), and it also fails, during preparation: Mooncake 0.5.61 `interface.jl:669-670` resets the pullback with zero rdata. That comment is stale for this version, and changing the preparation API does NOT fix true reverse mode.
7. Scope of the limitation: it is specific to the tested LinearSolve 5.18.2 / Mooncake 0.5.61 cached-solve rule. Non-cached `solve(LinearProblem)` is correct (MRE case 4), and Mooncake through the full PDE stage is correct (3.3e-13). It is NOT a general limitation of Julia AD.

## Update: code-generation sensitivity (independent evidence of round-off amplification)
The full suite (`Pkg.test` runs with `--check-bounds=yes`) and the focused run (default flags) give different 10a numbers for the same code. A focused run with `--check-bounds=yes` (`chunk10/focused10a_run4_checkbounds_yes.log`) reproduces the suite BITWISE:

| metric | `--check-bounds=yes` | default flags |
|---|---|---|
| final Xe vs native | 1.7057e-6 | 1.0747e-6 |
| pass-0 2p | 2.39 units | 3.22 units |
| 2g 3s, feedback range | 1.08e-3 | 1.4e-4 |

Bounds checks change code generation (SIMD and summation order) at the ulp level, and the composition amplifies this to ~1e-6 in the final Xe. This is consistent with the 128-bit primal deviation of 4.4e-6 on the frozen map. It is evidence about conditioning; it does not change any gate.
- nbscale direction, frozen map, 128 bit (`chunk10/diag10b_big_nbscale.log`, vectors in `diag10b_big_nbscale_vectors.txt`):
  - the 128-bit central differences converge (h = 1e-6 and 1e-8 agree to 6 digits);
  - Float64 ForwardDiff vs 128-bit FD: 4.87e-6 at worst (Xe at z = 1000), with Te outputs <= 6.8e-8;
  - Float64 vs 128-bit primal: 4.4e-6.
  So hscale and nbscale are now checked on the frozen map; F and A2s1s (tail-only parameters) and the adaptive public map are not.
8. REDUCED composition (2 ODE passes, 1 PDE stage, NO tail; NOT the production objective) with explicit ForwardDiffSensitivity in every solve and the feedback passed as 796 DF node values (`chunk10/probe_mc_composed_reduced_fds.log`, vectors in `_vectors.txt`):
   - the prepared Mooncake gradient now tracks ForwardDiff on the identical map: hscale -3.17544 vs -3.17342, nbscale 3.33734 vs 3.33526;
   - relative error 6.2e-4, against 1.0 under the default GaussAdjoint;
   - the inactive F and A2s1s gradients are exactly 0;
   - it FAILS the 1e-6 gate;
   - cost: cold prepare 652 s, gradient 265 s, 15.3 GiB.
   Hypothesis under test: discretization error of the per-node forward sensitivities (primal-only error control), as in probe 2. Convergence check with reltol 1e-12 for both routes is running (`probe_mc_composed_reduced_fds_tol12.log`).
   - Convergence check at reltol 1e-12 (`probe_mc_composed_reduced_fds_tol12.log`): g = [-3.17549, 3.33735] and J'w = [-3.17342, 3.33528], i.e. 6.19e-4. UNCHANGED from 1e-10, so the discrepancy is SYSTEMATIC, not discretization. Every `_primal` in src/ was audited: all are decisions or fixed redshifts, none strips a p-dependent value inside the reduced objective. Next: arbitrate one feedback block along the physical direction dDF/dhscale with 128-bit FD.
9. Background-scale link of the PDE stage (fixed native rows, parameters hscale and nbscale only):
   - DEFECT FIXED: `hi_diffusion_stage` used the unscaled accessors when hscale == nbscale == 1.0 were Float64 values. That preserves the primal, but it erased the explicit background derivatives for Mooncake, which sees plain Float64 parameters. The stage now always builds `ScaledBackground`, and the type promotion still includes h/nb. Focused regression `test/chunk10d_background_link_ad.jl` (exact unity and nearby values; Float/Dual primal bitwise; non-zero derivative; Mooncake VJP vs ForwardDiff < 1e-6, with per-block-max weights). It runs when memory frees.
   - Link values with 1/|f0| weights (`probe_link_pde.log`, `probe_link_pde_big.log`, vectors in `probe_link_pde_big_vectors.txt`):

     | route | J'w | relative error vs the 128-bit derivative |
     |---|---|---|
     | 128-bit ForwardDiff on the frozen map (exact reference) | [-0.3407908528, +0.3407908528] (exactly antisymmetric) | — |
     | Float64 ForwardDiff | | 4.55e-4 |
     | Float64 Mooncake | | 8.87e-5 |

   - ulp-level input perturbations move Float64 ForwardDiff by ~4e-4 (`probe_link_pde_ulp.log`).
   - This is evidence of strong conditioning: both Float64 routes are limited at 1e-4 to 5e-4 on this functional. It is NOT a passing 1e-10 link gate and NOT a root-cause proof. Gates unchanged.
   - Earlier evidence: an adaptive-solve 128-bit FD arbitration of one feedback block was INVALID (p-dependent steps; `probe_block_physdir.log`). On that block ForwardDiff and Mooncake+FDS agree to 3.6e-8.
10. Reduced composition rerun AFTER the background-scale fix (`probe_mc_composed_reduced_fds_bgfix.log`, vectors in `_bgfix_vectors.txt`):
    - g = [-3.1733722, 3.3352699] vs ForwardDiff J'w = [-3.1734212, 3.3352569], i.e. 1.47e-5 (was 6.2e-4 before the fix); F/A2s1s exactly 0;
    - still FAILS the unchanged 1e-6 gate.
    The fix removed the dominant systematic error, which confirms that the erased background derivative was a real defect. The remaining 1.5e-5 is comparable to the measured Float64 derivative conditioning (5e-6 composed vs 128-bit; 4.6e-4 / 8.9e-5 at the stage link). That is a plausible explanation, not a proof.
11. `test/chunk10d_background_link_ad.jl`, registered (`chunk10/focused10d_run1.log`): 15/16, ONE FAILURE.
    - Exact unity: the background derivative is non-zero (it was zero on the old fast path); Float/Dual primal bitwise at all 4 points.
    - Mooncake vs ForwardDiff: 3.0e-7 at [1, 1], 3.5e-7 at [1.001, 1], 5.2e-7 at [1, 0.999], **1.0009e-6 at [0.998, 1.002] > 1e-6** (FAIL; the gate is not broadened).
12. FULL composed objective, reverse mode by checkpointed chain rule (one prepared Mooncake gradient per link, inner ForwardDiffSensitivity), vs ForwardDiff on the identical map, running: `probe_mc_composed_full_chain.jl/.log`. This is NOT true reverse through the Rodas5P steps (blocked by the cached-solve rule).
13. FULL composed objective, reverse mode by checkpointed chain rule (`chunk10/probe_mc_composed_full_chain.log`, vectors in `_vectors.txt`; 2548 s, 15.8 GiB peak):
    - Mooncake chain g = [-7.4e-15, 4.0e-21, -2.9800520, 3.0995585] vs ForwardDiff J'w on the identical map = [6.5e-17, -2.5e-19, -2.9800579, 3.0995636];
    - relative error **1.9e-6 > 1e-6: FAIL**.
    - F and A2s1s are ~0 on both routes (F cancels in the rescaled Recfast tail, as in 5e).
    - The link VJPs took 1348, 103, 875, 39 and 74 s (L4, L3, L2, L1, L0).
    For context only, not a pass: Float64 ForwardDiff on the frozen composed map is itself 5e-6 from the 128-bit derivative. Inner solves use ForwardDiffSensitivity; this is NOT true reverse through the Rodas5P steps.
