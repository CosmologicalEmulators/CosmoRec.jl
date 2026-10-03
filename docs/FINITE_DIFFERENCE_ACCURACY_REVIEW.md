# Review: finite-difference and gradient accuracy gates (10b, 10d, hybrid) — read-only assessment, 2026-10-02

No production code, test, gate, dependency or native source was changed. All probes are scratch diagnostics in `A/chunk13/`, where `A = cmbcheb_test/local_analysis/cosmorec_differentiability_20260930`.

**Configuration.** Probes use the current runmode-0 test/reference callback (`TOL5E` = reltol 1e-12, a1 1e-18, aex 1e-14; injected Rodas5P; the library has no solver default) and `--check-bounds=yes`, as Pkg.test does. Exceptions are labelled.

**Latest suite** (`A/chunk12/pkgtest_a1_run1.log`): 50812 passed, 6 failed. The gradient-type failures are 10b:42 3.042e-3 (bound 1e-4), 10b:56 2.720e-5 (1e-5), 10d:46 1.000935e-6 (1e-6), and the hybrid probe 1.3599712e-6 (1e-6).

## 1. What each gate measures (units; not fractional spectrum errors)

| gate | definition (source) | units |
|---|---|---|
| 10b:42 | `min_h max_{i,j} \|J_ij − FD_ij(h)\| \|p_j\| / \|f0_i\|` (`test/chunk10b_runmode0_ad.jl:33-42`, `elast_err`, `fd_jac` in `test/chunk5e_helpers.jl`). h ∈ {1e-2, 1e-3, 1e-4} relative; central FD of `full_out10_frozen` (Patterson decisions replayed from P0, **ODE steps adaptive**, re-decided at each p); outputs Xe, Te at the 9 `ZG5E` redshifts; p = [F, A2s1s, hscale, nbscale] | error in the elasticity d ln f / d ln p (dimensionless) |
| 10b:56 | `min_h max_i \|J v − FD_v(h)\|_i / \|f0_i\|` along v = randn(Xoshiro(10), 4) .* \|p\| (`:47-56`) | elasticity along a random direction (each parameter moves by h·randn_j relative) |
| 10d:46 | `max_j \|g_j − (J'w)_j\| / max_j \|(J'w)_j\|`; g = prepared Mooncake VJP, J'w = ForwardDiff; w = randn/blockmax/n (block-max weights); PDE stage on fixed native rows, no ODE; Patterson decisions primal (`test/chunk10d_background_link_ad.jl:28-46`) | relative error of a 2-vector of weighted-sum derivatives |
| hybrid (probe) | same form; outer Mooncake + inner ForwardDiffSensitivity vs ForwardDiff on the identical composed map | same |

## 2. 10b: the FD reference map, not ForwardDiff, dominates the failure

### 2.1 Measured facts (`A/chunk13/fd_noise_probe.{jl,log}`, `fd_analysis.{py,log}`, `fd_secants.log`; all primal vectors saved)

- **Same primal, different schedules.** At P0 the adaptive map A (= the 10b FD reference) and its fixed-step replay S give **bitwise-identical** outputs. `ForwardDiff.jacobian(full_out10) == ForwardDiff.jacobian(A)` bitwise.
- **ForwardDiff is stable across schedules.** ForwardDiff of A vs of S differs by 7.9e-7 (hscale) and 2.4e-6 (nbscale) in elasticity.
- **128-bit check, from an earlier configuration** (`chunk10/diag10b_big*.log`, a1 = 1e-16, before the a1 edit):
  - map: fixed steps and fixed Patterson levels;
  - 128-bit central FD at h = 1e-6 and 1e-8 agree to about 2e-8;
  - Float64 ForwardDiff differs from it by 5.6e-6 (hscale) and 4.9e-6 (nbscale) in elasticity.
- **What that 128-bit result is a reference for.** It references the *local smooth piece* of the fixed-step map only. Data and physical constants are Float64 and enter BigFloat exactly. Primal-decided branches (`eps_A_effective`, stencils) are taken on the 128-bit primal, and their agreement with the Float64 branches was not verified. It is NOT a reference for the adaptive map, the native code, or the current configuration.

**FD error vs h** (elasticity, max over the 18 outputs; the worst entry is always an Xe output at z = 800–1000):

| map, parameter | h = 1e-2 | 5e-3 | 2e-3 | 1e-3 | 5e-4 | 2e-4 | 1e-4 | 1e-5 |
|---|---|---|---|---|---|---|---|---|
| A, hscale | 2.8e-3 | 4.4e-3 | 1.7e-4 | 1.7e-4 | 2.3e-4 | 1.3e-3 | 2.1e-3 | 2.2e-2 |
| A, nbscale | 4.2e-3 | 5.5e-3 | 1.1e-2 | **3.04e-3** | 2.8e-4 | 1.1e-3 | 8.6e-3 | 3.7e-2 |
| S, hscale | 2.7e-3 | 4.8e-3 | 1.2e-4 | 8.5e-4 | 6.1e-4 | 1.6e-3 | 3.7e-3 | 1.5e-2 |
| S, nbscale | 3.9e-3 | 5.3e-3 | 1.1e-2 | 1.6e-4 | 1.9e-3 | 4.9e-4 | 2.6e-3 | 1.7e-2 |

At the test h set, the max over parameters, minimized over h, is 3.04e-3 (nbscale, h = 1e-3, Xe(z = 800)). This reproduces 10b:42 exactly. F and A2s1s columns: elasticities ≤ 5e-17; FD errors ≤ 3e-8.

### 2.2 Mechanism

1. **Model discontinuity: the `eps_A_effective` clamp** (`src/RateTable.jl:165` = native `get_effective_rates.HI.cpp:320`; the A coefficient jumps by about 1e-4 when \|e^fxy − 1\| crosses 1e-4).
   - **Metric.** Δ(d) = (f(d) − f0 − J·d·p)/\|f0\| is the deviation of the sampled 10b map from ForwardDiff's local tangent, relative, at fixed output Xe(z = 800); the numbers below are increments of Δ between adjacent samples.
     - Δ is a **non-smooth/discrete residual of the numerical parameter-response map**. It is NOT an error in Xe, NOT a native-vs-Julia discrepancy, and NOT comparable with the fiducial cross-code Xe gate (1e-5; measured 1.45e-6, PASS).
   - On the 10b map, these residual increments are (`jump_epsA.log`):
     - nbscale: +1.2e-5 at d ∈ (1.1, 1.2)e-3, +3.4e-5 at (1.6, 1.7)e-3, and ≈ 7e-6 between −5e-4 and −1e-3;
     - hscale (cubic-detrended response, `hscan_compare.log`): ≈ 3.2e-5 at d ∈ (3.50, 3.75)e-3.
   - Removing the clamp **in-process only** (counterfactual; `jump_epsA_probe.jl noclamp`) removes these discrete residual features. Δ stays ≲ 1e-6 over the scanned interval, and the ForwardDiff elasticity is unchanged (−0.8134893 vs −0.8134870).
   - At the worst 10b:42 entry (nbscale, h = 1e-3, Xe(z = 800)), FD − ForwardDiff (an elasticity residual of the test reference, not an Xe error) goes from −3.0416e-3 to **+1.443e-4** with the clamp removed (21× smaller; `jump_fd_check.log`; the supervisor reproduced it independently). At h = 2e-3 it is −1.05e-4, with a max over outputs of 2.1e-4.
   - 1.443e-4 is still above 1e-4, and the counterfactual is not the original model. **No actual gate passes.**
   - The fixed-step map S shows the same residual features, so they are not caused by the ODE step controller.
2. **Residual primal roughness** σ ≈ 3.5–5.2e-7 (rms, relative; cubic-fit residual over 41 points at 1e-6 spacing; worst at Xe(z = 800–1000); median over outputs ~1e-8).
   - It is the same with the clamp removed (4.2e-7) and in the fixed-step map (4.3e-7), so it is neither the clamp nor the step sequence.
   - **Origin unresolved:** Float64 rounding amplified by the PDE-stage conditioning (κ up to 1.9e8), or other primal-decided branches such as the rate-table/polint stencils. The `14/41 distinct values` in the nbscale scan is Te(z = 2500), with elasticity 6e-10 (about 4 ulps per step). That is quantization, not noise.
   - Independent noise contributes about 0.71 σ/h to a central FD: about 3e-4 at h = 1e-3 and 3e-3 at h = 1e-4. That matches the measured small-h errors.
3. **Adaptive-step residual features** (≈ 5e-6 in Δ; e.g. A vs S at nbscale d = −1e-3: Δ ≈ 7e-6 vs 2e-6) exist but are second-order here.
4. **Truncation:** the cubic/quintic-fit estimate (e3 ≈ 170–280) overpredicts the h = 1e-2 error by 10×, because the fit absorbed the discrete residual features. It is **PROVISIONAL and not used**.

The "noise floor ~3.4e-4" fit in `fd_analysis.log` is likewise **PROVISIONAL**. It is a model, not a measured floor, and its residual includes discrete features and possibly higher-order terms.

### 2.3 Classification

- **10b:42** is an FD-reference / map-semantics failure. The reference samples across clamp-induced non-smooth features of the map at h ≈ 1e-3 and is roughness-dominated below that. The clamp-on/off FD result (3.04e-3 → 1.44e-4) is a test-reference effect, not a model-accuracy pass. It is **not evidence that ForwardDiff misses the gradient**: on the local smooth piece, ForwardDiff agrees with the 128-bit FD of the fixed-step map to about 5e-6 (old configuration).
- **What ForwardDiff does give:** the derivative of the local smooth branch (the almost-everywhere derivative). It ignores the clamp's discrete contributions by construction. Whether the public sensitivity should account for them is a **contract question** (§6), not an AD bug.
- **10b:56:** the same reference map. The best measured value, 2.72e-5, is the size of the noise term at h = 1e-2 (0.71 σ / 1e-2 ≈ 3e-5). The clamp's share along the random direction was not isolated. FD-diagnostic failure; not decomposed further.

## 3. 10d: an exact 128-bit reference on the same map (`A/chunk13/d10_bigfloat_reference.{jl,log,txt}`)

Setup: the exact 10d functional (rows, weights and points). The Float64 Patterson decisions are recorded at each p and replayed in 128-bit ForwardDiff. Replayed Float64 ForwardDiff == logged test values, bitwise.

| p | Mooncake vs 128-bit | ForwardDiff vs 128-bit | test metric \|g − r\|/max\|r\| | 128-bit invariant h·∂h + nb·∂nb |
|---|---|---|---|---|
| [1, 1] | 2.2e-7 | 2.4e-7 | 3.0e-7 | 0 |
| [1.001, 1] | 5.0e-7 | 5.0e-7 | 3.5e-7 | 0 |
| [1, 0.999] | 5.4e-7 | **1.05e-6** | 5.2e-7 | 0 |
| [0.998, 1.002] | 8.8e-7 | 1.4e-7 | **1.0009e-6** | 0 |

- Both Float64 AD routes are accurate to 1.4e-7 – 1.05e-6 against the exact derivative.
- The gate's own reference, ForwardDiff, exceeds 1e-6 at one point. The failing difference is Mooncake 8.8e-7 + ForwardDiff 1.4e-7.
- **Classification:** the Float64 accuracy floor of two correct AD routes on a cancellation-dominated stage. There is no AD-backend bug signature. A 1e-6 bound on the *mutual difference* is at the per-route floor.
- The exact ratio invariant gives a free, reference-free accuracy check. The Float64 routes violate it by 1.6e-8 – 9.6e-7 (`10d_invariant.log`).

## 4. Hybrid (outer Mooncake + inner ForwardDiffSensitivity), 1.3599712e-6

- It compares two Float64 derivative routes on the identical map; there is no exact reference for the full chain.
- The measured route floors are about 1e-6 for the stage (§3) and about 5e-6 composed (fixed-step, old configuration). 1.36e-6 is consistent with them, but causality is **unresolved**.
- It is separate from the B2 infrastructure blocker: true reverse mode through Rodas5P hits the reproduced LinearSolve 5.18.2 / Mooncake 0.5.61 cached-solve rule failure.

## 5. Native vs Julia around the eps_A threshold

- **State space, at the threshold** (`chunk11/rhs_scan`, investigation §5.1–5.4):
  - The native and Julia right-hand sides are identical to ≤ 1.1e-13, including the same 3d jump at z = 938.9617.
  - The native's pass-0 output error there (−2.82e-9 X1s) is **native solver sensitivity**: its error controller never rejects, and one accepted order-5 step crosses the jump. The native step cap removes it, converging to the Julia value within 7.4e-13.
- **Parameter space** (`A/chunk13/native_hscan.sh`, `hscan_compare.{py,log}`; H(z) → H(z)(1 + d), 29 points, all runs rc 0 and complete):
  - **Metric** (`hscan_compare.py`): Y(d) = Xe(d)/Xe(0) − 1 at a fixed output z; R(d) = Y minus a least-squares cubic in d; the reported number is the largest neighbour increment of R. It is a **non-smooth residual feature of each code's sampled parameter-response curve**, NOT a step in Xe at fixed cosmology and NOT a native-vs-Julia Xe error. It cannot be compared with the 1e-5 fiducial cross-code Xe gate (measured 1.45e-6, PASS).
  - **The maps are not the same.** The native scan rescales the native H table, including the native's own initialization (preliminary history, Saha state). Julia's `hscale` multiplies H inside the solver, with the initial state and loaded closures frozen. Only the presence and location of non-smooth features can be compared, not smooth slopes or values.
  - Largest residual increment at Xe(z ≈ 800): native 3.15e-5 (original) and 3.19e-5 (step cap 1e-4); Julia 3.21e-5. All at d ∈ (3.50, 3.75)e-3. At Xe(z ≈ 500): 6.1 / 6.0e-6 native vs 6.2e-6 Julia, same interval.
  - Reading: the feature is unchanged by native step capping, so it is not native step-size sensitivity, and it disappears in Julia with the clamp removed. This is **consistent with** a shared non-smoothness of the original model (the clamp) appearing in both numerical maps. The native-side attribution to the clamp is inferred, not run: a native clamp-free run would be a model change.
  - Mismatch: at z ≈ 1000 the native residual feature is larger (1.1e-5 vs 1.4e-6), which is consistent with the maps not being identical.

## 6. Reasonableness of the bounds, and what to do instead (proposals only; no gate changes)

- **1e-4 elasticity (10b:42).** Inherited from the single-pass 5e test, where FD_best is 1.15e-5. For the composed map with this reference it is below the measured FD floor: ≥ 3e-3 on the actual model at the test h, and 1.4e-4 at best on the smooth counterfactual. **Unsupported** as an AD-correctness criterion with this reference.
- **1e-5 JVP (10b:56).** Below the measured roughness term (about 3e-5 at h = 1e-2). Larger h runs into the clamp-induced residual features. **Unsupported** with this reference.
- **1e-6 VJP (10d) and the hybrid 1e-6.** Not derived; at the measured per-route Float64 accuracy (up to 1.05e-6 against exact). As written they cannot separate a bug from Float64 noise.
- **Studies instead of loosening:**
  - Same-map high-precision references: record and replay ALL primal branches (steps, Patterson, and per-evaluation clamp/stencil decisions), then compare ForwardDiff with a 128-bit derivative of that exact piece. The bound would be derived from the measured Float64-vs-128-bit primal error times the conditioning.
  - Exact invariants, such as the 10d ratio invariant.
  - FD step-convergence curves over ≥ 10 values of h, with the branch decisions frozen.
  - Peak- or block-normalized denominators where outputs approach zero.
- **The model has a hard threshold** (the clamp). Options, each needing Marco's decision:
  - (a) Event handling: locate the crossings and restart. The map stays piecewise smooth, with well-defined discontinuities.
  - (b) Domain restriction: not viable, because crossings occur throughout the history.
  - (c) Smoothing the clamp: a model departure from native parity.
  - (d) A different public sensitivity contract. For example, the derivative of the local smooth branch (ForwardDiff semantics) with a documented bound on the non-smooth residual of the parameter response (measured: neighbour residual increments ≈ 3e-5 of the relative response Y at Xe(z = 800), for 2.5e-4 parameter spacing near d ≈ 3.5e-3, metric as in §5), or a scale-Δ smoothed sensitivity.
- **Three separate questions:**
  - Numerical correctness: ForwardDiff is accurate to about 1e-6 (stage) and about 5e-6 (composed, fixed-step) against exact references on the same piece.
  - Native equivalence: at fixed cosmology the right-hand sides match to ≤ 1.1e-13, and the fiducial cross-code final Xe is 1.45e-6, which passes the approved 1e-5 gate. The parameter-response scans of the two codes (different maps) show non-smooth residual features at the same interval.
  - Scientific impact: NOT measured. The non-smooth parameter-response residuals are a possible sensitivity risk for gradient-based use, but their C_ℓ / observable impact is unquantified. They are not an Xe accuracy failure.

## 7. Scope: registered tests vs scratch probes; native independence

- **Registered tests** (actual values; no gate changed): 10b:42 3.042e-3, 10b:56 2.720e-5, 10d:46 1.000935e-6. All still FAIL.
- **Scratch probes** (diagnostic only, never acceptance evidence):
  - the clamp-removed (`noclamp`) runs are a counterfactual model, not a fix;
  - the 128-bit references cover only the same frozen-decision piece;
  - the noise-floor fit is a model.
- **Native code.** The native code is an independent implementation, and nothing here changes it. The native H-scan comparison is a benchmark/accuracy concern about the shared model's non-smooth parameter response: both codes' (non-identical) response scans show non-smooth residual features at the same parameter interval. It is not evidence of a native defect, and not a validation of the Julia derivative.

## 8. Remaining genuine showstoppers and uncertainty

- **Showstoppers** (not resolved by this review):
  - B2, true reverse mode through Rodas5P (the LinearSolve/Mooncake cached-solve rule);
  - B3, whole cosmological initialization (not implemented; H(z) policy pending);
  - the gate policy for 10a:119/120, where the native fails against itself;
  - FD/VJP gates that are not derived from a valid same-map reference.
- **Possible risk, not a showstopper or a failed gate:** clamp-induced non-smooth residuals of the parameter-response map. For gradient-based inference they matter; their C_ℓ-level impact is unmeasured and needs a dedicated sensitivity study.
- **Uncertain:**
  - the origin of the σ ≈ 4e-7 roughness (rounding vs other primal branches);
  - the truncation level;
  - the clamp's share in 10b:56;
  - the full-chain hybrid causality;
  - the native-side clamp attribution (inferred, not run);
  - whether the old-configuration 128-bit result (a1 = 1e-16) carries over unchanged to the current callback.
- **Not claimed:**
  - that the FD error is pure roundoff;
  - that the noise floor is 3.4e-4;
  - that removing the clamp is a fix;
  - that 10d is a serious gradient failure;
  - that any error is scientifically harmless.
