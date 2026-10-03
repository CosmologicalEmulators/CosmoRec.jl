# Original CosmoRec vs CosmoRec.jl: attribution of the remaining runmode-0 discrepancies

> **CORRECTION (2026-10-02, after the matched-flag ablation, §5.6).** An earlier draft blamed 10a:103 on the excited-state tolerance `aex = 1e-14`, and compared runs with different compiler flags. The supervisor flagged both points, and the ablation **refutes** the aex attribution.
> - **Cause.** The factor that decides the gate is **a1**, the absolute tolerance of rho, X1s and He1s.
> - **a1 = 1e-18 passes everywhere.** For all 6 reltol/aex combinations, under both `--check-bounds=yes` and default flags, 10a:103 passes at 5.12e-5 against the original. That value is the original's own step error.
> - **a1 = 1e-16 (production) fails or varies.** Every case lies between 1.4e-4 and 1.46e-3, is not monotone, and depends on the flags.
> - **aex has no consistent one-factor effect.**
>
> **Other findings.**
> - The production configuration under `--check-bounds=yes` reproduces the test values exactly.
> - Runs with the same flags are bitwise reproducible across processes.
> - The excited-state budget weakness (aex above every excited population) remains a source-derived observation; it is not measured to affect any gate.
>
> **Superseded wording.** §1, §5.6, §7, §8 and §9 have been rewritten. The earlier mixed sweep is kept in §5.6 only as a record of the earlier reasoning.
>
> **Unaffected.** The source-verified native findings are not touched by this correction: the native solver never rejects a step on its error estimate (§4.1, §5.2), the `eps_A_effective` pass-0 step and its step-cap convergence (§5.1–5.4), the pass-1 forcing and the native self-failure of 10a:118/119 (§5.7–5.9), and the 10b reference noise (§5.11).

> **Terminology (interface scope).** In this document, "production configuration" or "production" tolerances means the runmode-0 test/reference callback configuration used for these native comparisons: the Rodas5P callbacks `SOLVE10` / `SOLVEP5E` / `SOLVE5B-5D` injected by the tests. It is NOT a built-in library setting. `recombination_pass` (`src/RecombinationODE.jl:148`), `recombination_history` (`:319`) and `recombination_history_diffusion` (`src/RecombinationDiffusion.jl:62`) take a required caller-supplied solve callback, and `Project.toml` has no ODE-solver dependency or default. Downstream callers choose their own solver and tolerances.

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al. Use must be acknowledged; cite Chluba & Thomas 2010 (MNRAS 412, 748). The original source, library and data tree was used read-only (hashes verified at the end of the investigation). No production Julia code, test, gate, dependency or native file was changed. No Git operation was performed.

Date: 2026-10-02. Paths are relative to `A = cmbcheb_test/local_analysis/cosmorec_differentiability_20260930` unless they start with `src/`, `test/` or `docs/`.

## 1. Question and short answer

The question was whether the remaining original-vs-Julia differences are:
- porting mistakes;
- accuracy limits of the original solvers;
- limits of the Julia solvers;
- harness or configuration mistakes;
- or a mixture.

**Short answer: a mixture, with no porting mistake found in the physics of any stage.**

| Failing gate (latest run, `chunk10/pkgtest7_10_run2.log`, `chunk10/focused10a_xe1e5_checkbounds.log`) | Value / bound | Classification | Status of the evidence |
|---|---|---|---|
| 10a:103 (×2), full-path PDE0/PDE1 2g3s, max\|Δ\|/max in 500<z<2000 | 1.075e-3, 1.039e-3 / 1e-3 | **Julia solver configuration (our mistake): `a1 = 1e-16`.** Above z ≈ 1813, X1s (7e-10 at z = 3000) has an absolute-dominated Rodas weight, an effective relative tolerance of up to 1.3e-7. The PDE stage is hypersensitive to X1s. Matched-flag one-factor result: a1 1e-16 → 1e-18 makes 103 pass in all 12 reltol/aex/flag cases (5.12e-5). aex has no consistent effect. | Confirmed by the matched-flag ablation (§5.6) |
| 10a:118, pass-1 X1s | 3.25 / 1 unit (1e-7 rel) | **Original-formulation conditioning, together with the original's own numerical noise.** The pass-1 solve is forced by PDE0 errors, and PDE0 is hypersensitive to pass-0 populations. The Julia pass-1 integration is faithful. The **native fails this gate against itself.** | Confirmed mechanism (§5.7–5.9). The exact-arithmetic conditioning floor is NOT established. |
| 10a:119, pass-1 2s/2p | 15.4, 18.3 / 10 units (1e-5 rel) | Same as 10a:118 | Same as 10a:118 |
| 10a:82, injected Julia PDE stage vs native, `:scaled` metric | 1.66e-4 / 1e-4 | **Stage implementation at the round-off/conditioning level, not a port error.** At identical inputs the difference is 3.6e-7 of peak. The native's own 1e-15-input response reaches 1.6e-4 in the same metric. | Confirmed comparison (§5.10). The origin of the 3.6e-7 is a hypothesis (evaluation order) |
| 10a:117, final Xe | 1.71e-6 / 1e-5 (approved) | **Passes.** It lies at the native's own reproducibility: 1.0–1.8e-6 under 1e-15 H noise, and 1.6e-6 against the step-converged native. | Confirmed |
| 10b:42, 10b:56, composed Jacobian vs Float64 FD | 2.4e-3 / 1e-4; 3.5e-5 / 1e-5 | **Test design.** The Float64 FD reference is limited by primal noise, and a 128-bit BigFloat FD agrees with ForwardDiff to 5.6e-6. No native involvement. | Confirmed on the frozen-decision map only (§5.11) |
| 10d:46, Mooncake VJP vs ForwardDiff, PDE stage | 1.0009e-6 / 1e-6 | Julia-internal agreement of two Float64 AD modes; not an original-vs-port question. | Not investigated further here |

**New finding about the original, from instrumented native evidence (§4.1, §5.2–5.4):**
- The original Gear/BDF error estimate **never rejects a step**. The rejection test at `ODE_solver_Rec.cpp:1244` cannot be reached.
- Steps are accepted with local error estimates up to 217× the tolerance in pass 0, and up to 2.8e5× in passes 1–2.
- At an original-model discontinuity (`get_effective_rates.HI.cpp:320`, the `eps_A_effective` switch), one accepted order-5 step leaves a 2.8e-9 X1s error that **no tolerance setting removes**. Capping the native step size removes it and makes the native converge to the Julia solution, which agrees to 7.4e-13.

## 2. Error definitions (common to every table; `chunk11/compare/cmp.py`)

- pointwise = max\|a−b\|/\|b\|; peak = max\|a−b\|/max\|b\|; zmax = z of the pointwise maximum.
- Comparison ranges:
  - PDE outputs: 500 < z < 2000 (the feedback range; native `interpolate_DF` returns 0 outside it).
  - FINAL: z < 3000 (this document's `cmp.py` comparisons only; the 10a:117/118 test assertion itself takes the maximum over all 10000 nodes of `ZREC10`).
  - Pass rows: all 3000 nodes.
- "Units" (10a gates) = \|Δ\|/(rtol\|b\| + atol), with the native solver tolerances: X1s/rho 1e-7/1e-12, 2s/2p 1e-6/1e-80.
- `:scaled` (10a:82) = max over all z of \|Δ\| divided by the Julia residual scale `o.scale` (columns 6–9 of the sweep PDE files).
- Near-zero references: the 2g3s/2g3d/R2s outputs change sign in the feedback range. Their pointwise errors (1e-2 to 2e-1 in every comparison, native-vs-native included) are dominated by points near zero. Where the two disagree, this document uses the peak-normalized value as the observable-level number.

## 3. Controls and provenance

- **Scratch copies** (`chunk11/native_scratch/`). There are two.
  - `CosmoRec_unmodified/` is an rsync of the original tree without objects, temp or outputs. Its `Rec_database` is a symlink, used read-only. All 155 sources are hash-identical (`source_hashes_original.txt` = `source_hashes_copy.txt`).
  - `CosmoRec_diag/` is a copy with isolated numeric-control and logging changes only, activated by environment variables:
    - v1 (`diag.diff`, used by `sweep/`): `COSMOREC_DIAG_TOL_SCALE`, `COSMOREC_DIAG_MONITOR_ALL`, `COSMOREC_DIAG_STRICT`;
    - v2 (`diag_v2.diff`, 49 changed lines in `Modules/main.CosmoRec.cpp` and `Development/ODE_PDE_Solver/ODE_solver_Rec.cpp`): adds `COSMOREC_DIAG_TRACE` (logs every step attempt; changes no numbers) and `COSMOREC_DIAG_DZZ_MAX` / `COSMOREC_DIAG_DZZ_ZWIN` (caps Dz/z).
    - No physical equation was touched.
  - Build: `make lib` with `GSL_INC_PATH/GSL_LIB_PATH = tools/camb-cosmorec-env`, flags `-Wall -pedantic -O2 -fPIC`. The harness is linked by `build_harness.sh` (the same flags as the chunk10 harness).
  - Hashes: `hashes.txt` (v1) and `hashes_v2.txt` (v2: diag_v2.diff 96caf3bf…, lib 14b06266…, harness 52bb75d6…, the v1 harness kept as `harness10a_diag.pre_trace` 5ba3852b…).
- **Bitwise baselines.** Each probe10a.txt below is byte-identical to the original-library capture `chunk10/native_capture/out/run1/probe10a.txt`. The supervisor independently confirmed the first two.
  - `out_unmodified/`: unmodified scratch build;
  - `out_diag_default/`: diag v1 with the environment unset;
  - `out_trace_off/`: diag v2 without the trace;
  - `out_trace_on/`: diag v2 **with the trace on**;
  - `out_dzcap_default/`: diag v2 with the cap unset.
  - In every one of them the replica Xe/Tb is bitwise equal to the original.
- **Completeness of compared runs.** Every native case used in a comparison has 15427 lines and per-tag counts identical to the original capture (checked by `run_dzcap.sh` / `status.txt`, and for `sweep/` and `hpert/` in this session).
  - The v2 runs record rc = 0 and exit status 0 (`dzcap/*/status.txt`, `time.txt`).
  - For `sweep/` and `hpert/` the exit status was not recorded at the time, but their outputs are complete.
  - Failed or stalled runs are never compared. They are kept as negative results (§6).
  - `native_scratch/run_monall.log` records a failed launcher (a missing `run_monall.sh`). The monall outputs used here come from the separate sweep loop (`sweep/monall_*`, complete).
- **Common inputs.** Every native and Julia run in this document has the same:
  - CAMB H table, `chunk10/native_capture/camb_H_input.txt`;
  - parameters (`theta5()`), flags (runmode 0, `Diffusion_correction` 1, `DI1_2s_correction_on` 1, `Diff_iteration_max` 2, nS_2γ 3, nS_R 2, induced_flag 2);
  - grid (`GRID5`, 10000 nodes; PDE zs 2500, ze 500), and helium switch (k_switch = 1342, z = 1680.9103034346122, asserted in 10a and printed by every sweep).
  - The two earlier harness defects (the wrapper clearing the loaded H; the DI1 flag ordering) were fixed before chunk 10 and are not present in these captures.
- **Julia side.**
  - Runs use the public adaptive map `recombination_history_diffusion(RM5, theta5(), D, solver, solve_tail5, GRID5[:,1])` (`chunk11/julia_sweep/julia_public_sweep.jl`), with Rodas5P and `internalnorm = primal_norm`.
  - Tolerance parameters: reltol, `a1` (rho, X1s, He1s) and `aex` (excited states). The production test configuration is reltol 1e-12, a1 1e-16, aex 1e-14 (`test/chunk10a_runmode0_native.jl:27`).
  - Default compiler flags; the 10a test log was run with `--check-bounds=yes`.

## 4. Source reading (original, line numbers of the unmodified file)

### 4.1 Gear/BDF solver, `Development/ODE_PDE_Solver/ODE_solver_Rec.cpp`

- **Error estimate** (`ODE_Solver_estimate_error_and_next_step_size`, 1123–1250).
  - Predictor (extrapolation) minus corrector.
  - Only components with `rel_vector[l] != 0` are checked (1161).
  - Per component (`ODE_Solver_error_check`, 1095–1116): `err_aim = rtol|y|`, or `atol` if that is larger. `fac = err_aim/|Δy|`; the minimum over components is used.
  - `power = 1/order` (1157).
- **Step-size update** (1192–1197): `Dz_z_new = Dz_z·min(5, fac^power)`, clipped to [1e-10, Dz_z_max]; down_thresh 0.7, up_thresh 1.3.
- **Rejection is unreachable.**
  - In the `Dz_z*down_thresh > Dz_z_new` branch (1219), `Dz_z = Dz_z_new` is set first (1221). After `n_down >= 3` it is divided further by 2^(n_down−1), or reset to `Dz_z_min`.
  - Only then is `if(Dz_z*0.5 > Dz_z_new) return 1;` tested (1244). With `Dz_z ≤ Dz_z_new` and both positive, it is never true.
  - All other paths return 0. The error estimate therefore only shrinks the *next* step and never redoes the current one.
  - **Instrumented confirmation** (§5.2): 0 rejections in 18672 step attempts.
- **Redo only on Newton failure.**
  - Newton (1276–1395) re-evaluates the Jacobian at every iteration (1322).
  - The linear solve is Gaussian elimination for neq ≤ 50 (360).
  - Convergence is checked only on monitored components: \|ΔY\| < max(rtol\|y\|, atol) (1353). `tolJac = tolSol/2` (1276) with `tolSol = minrtol` (666).
  - After 25 loops it returns 1 (1377). `Solve_history` then halves the step and retries (1472–1478).
  - Measured: at most 2 Newton loops per step, and no Newton failure, in the default run.
- **Steps are clipped to the output node** (1456). Each `Solve_history(zs, zend)` call (1407–1530; default `Dz_in` 1e-8, `Dz_max = 0.5 zs`, 1536) integrates one output interval, so a step never exceeds the node spacing. `Dz_z_last` carries over between intervals (1523).
- **Tolerances** (`Modules/main.CosmoRec.cpp:16–36`):
  - rho and X1s: 1e-7/1e-12;
  - 2s and 2p: 1e-6/1e-80;
  - He1s: 1e-7/1e-12; He 2¹S and 2¹P: 1e-6/1e-80;
  - 3s, 3p, 3d, He 2³S and 2³P1: unmonitored (rtol = 0). They enter neither the error estimate nor the Newton convergence test.

### 4.2 Effective rates, `Rec_database/Effective_Rates.HI/get_effective_rates.HI.cpp`

- `const double eps_A_effective = 1.0e-4` (71); `if(fabs(exp(fxy)-1.0) <= eps_A_effective) fxy = 0.0;` (320 for A, and the same at 416).
  - This is a **discontinuity in the original model**: when the interpolated log-correction crosses \|e^fxy − 1\| = 1e-4, A_nl jumps by ≈1e-4 relative.
  - Julia reproduces it deliberately (`src/RateTable.jl:165`).
- Stencil: `locate`/`stencil_start` (`src/RateTable.jl:95–144`) mirror native 233–285. Four-point Lagrange interpolation is continuous in value but not in slope when the stencil moves.

### 4.3 PDE stage, quadrature, feedback

Previously established line by line (`docs/CHUNK7_9_RESULTS.md`):
- grid and profiles: bitwise or ≤ 1.7e-14;
- def_PDE coefficients: ≤ 2.9e-13 (the Dnem cancellation);
- the march: final spectrum 2.8e-15 of maximum (`Solve_PDEs.cpp:786–975`, `Step_PDE_O2t`);
- the Patterson correction integrals: ≤ 1e-4 of the residual scale; 59–90 of 199 outputs contain unconverged sub-integrals in **both** codes, because the stopping decisions are reproduced;
- `interpolate_DF` / `fcn_effective` with feedback: 3.5e-13.

The Dnem/DnLj cancellation reaches κ = 1.88e8 on the fixed native history (7c). This conditioning belongs to the original formulation and dominates everything downstream (§5.5–5.9).

## 5. Experiments and results

### 5.1 Right-hand side at identical states (`chunk11/rhs_scan/`)

- **Config:** 2001 linearly interpolated flag_He = 0 states between native pass-0 nodes 2095 (z = 939.2297) and 2096 (z = 938.2461).
  - Native: the ORIGINAL `fcn_effective(…, -1)`, which is the call the solver uses (`ODE_solver_Rec.cpp:400`), linked to the original library (`harness_rhs_scan.cpp`, `build.sh`).
  - Julia: `recombination_rhs(z, y, RM5; flag_He = false)` (`julia_scan.jl`).
- **Result** (`scan_compare.log`, `scan_jumps.log`):
  - rho, X1s and the n = 3 RHS agree to ≤ 1.1e-13 relative everywhere.
  - 2s/2p agree to ≤ 1e-13 except at node 2096 itself. There the 2s/2p derivative cancels to 5e-17 (against ~1e-7 at neighbouring points), giving 2.5e-3 relative but 4e-20 absolute: a near-zero reference.
  - **Both codes show the same jump** in the 3d RHS at z = 938.9617 (first difference 6.1e-8 against 3.1e-9 typical).
- **Cause** (`epsA_check.jl/.log/.txt`): fxy(3d) crosses 1e-4 (9.9995e-5) between z = 938.96170 and 938.96120, so the `eps_A_effective` branch switches. The stencil (lx, ly) = (385, 35) does not change. This is a common, original-model discontinuity.

### 5.2 Native step trace (`native_scratch/out_trace_on/trace.tsv`, `trace_summary.log`; output bitwise = original)

The harness runs the original and its replica, so the trace has 12 ODE segments: three passes × (He on, He off) × 2. The replica segments are identical to the original's.

| pass, segment | steps | rejected by error estimate | accepted with error > tol | max error / tol (component, z) |
|---|---|---|---|---|
| 0, He on (3000→1680.9) | 1413 | 0 | 75 | 36.8 (He 2¹P, 1691.7) |
| 0, He off (1680.9→50) | 1695 | 0 | 35 | 217 (2p, 936.28) |
| 1, He off | 1700 | 0 | 40 | 2.7e5 (2p, 497.57) |
| 2, He off | 1701 | 0 | 133 | 2.8e5 (2p, 497.57) |

Newton: 1 loop for 17790 attempts, 2 loops for 882, no failure.

**Pass 0 around the anomaly:**
- From 939.230 to 938.246 the native takes a single order-5 step over the whole output interval (Dz/z = 0.5 cap; the node spacing limits the step).
- The 2p error estimate goes from fac = 3.5e3 on the previous step to 0.0289, i.e. 35× the tolerance, because the step straddles the 3d `eps_A` jump. The step is accepted.
- The next steps (Dz_z 0.25 → 0.09 → 0.03, still node-limited) carry 141–217× tolerance, because the BDF history now contains the jump.

**Passes 1–2:** the worst steps sit at z = 497.57, just below the feedback cut at z = 500, where `interpolate_DF` drops to 0. This is a hypothesis from the location only and was not separately tested.

### 5.3 Why the pass-0 X1s step is tolerance-independent (`compare/node2096_julia_vs_native.log`)

Pass-0 X1s minus the resolved Julia solution (rt1e-14/aex1e-20), absolute, at native nodes 2095 / 2096 / 2097:

| run | node 2095 | node 2096 | node 2097 |
|---|---|---|---|
| native original, s1e-1, s1e-2, s1e-3, monall_s1e-3 | +3.1e-13 | −2.82e-9 | −3.13e-9 |
| Julia rt1e-10 | −8.3e-12 | −1.2e-10 | … |
| Julia rt1e-12 | −4.8e-14 | −3.3e-12 | … |
| Julia rt1e-13 | −7.4e-15 | +7.4e-14 | … |

- Native: identical to three digits for every tolerance scale. Scaling the tolerances only shrinks the *next* step, which is still clipped to the node, so the same accepted step is taken.
- Julia: its own error at this switch shrinks as the tolerance tightens, because Rodas5P rejects and refines.

### 5.4 Native step-size cap, native-only convergence (`native_scratch/dzcap/`, `run_dzcap.sh`, `compare/node2096_dzcap.log`)

Windowed cap of Dz/z for output intervals starting at 930 ≤ zs ≤ 945 (all runs complete, rc 0):

| cap | X1s at node 2096 minus resolved Julia | native self-change at node 2096 |
|---|---|---|
| none (original) | −2.82e-9 | |
| 1e-4 | +3.9e-10 | 1e-4 vs 1e-5: 3.6e-10 |
| 1e-5 | +3.0e-11 | 1e-5 vs 1e-6: 2.9e-11 |
| 1e-6 | **+7.4e-13** | |

The native converges at about first order, as expected across a jump, to the value the Julia adaptive map already produces. **The −2.8e-9 step is a native ODE-solver accuracy limitation at a common model discontinuity, not a port error.**

### 5.5 Native global step caps vs resolved Julia (`compare/dzcap_cmp_1e-5.log`, `dzcap_cmp_1e-6.log`; caps 1e-3 to 1e-6, wall 1.2 s / 5 s / 45 s / 7 min 42 s)

Peak-normalized values; pointwise values are in the logs. J14 = Julia rt1e-14/a1 1e-18/aex1e-20; J13 = rt1e-13/1e-17/1e-20.

| quantity | J13 vs J14 | native orig vs J14 | cap1e-6 vs J14 | native cap1e-5 vs cap1e-6 | native orig vs cap1e-6 |
|---|---|---|---|---|---|
| PASS0 X1s | 7.6e-13 | 3.1e-9 | 4.2e-12 | 4.3e-11 | 3.1e-9 |
| PASS0 2s | 4.8e-11 | 1.4e-8 | 1.8e-11 | 1.9e-10 | 1.4e-8 |
| PDE0 DI1 | 2.3e-6 | 3.8e-6 | 1.5e-7 | 2.1e-9 | 3.8e-6 |
| PDE0 2g3s | 5.3e-5 | 1.2e-5 | 5.5e-6 | 1.1e-6 | 1.2e-5 |
| PDE0 R2s | 2.4e-5 | 3.3e-5 | 1.9e-6 | 9.6e-7 | 3.3e-5 |
| PASS1 X1s | 2.8e-7 | 1.4e-7 | 1.3e-7 | 2.0e-7 | 1.9e-7 |
| PASS1 2s | 1.0e-6 | 5.0e-7 | 5.5e-7 | 9.4e-7 | 6.9e-7 |
| FINAL Xe (pointwise) | 1.1e-6 | 1.5e-6 | 8.4e-7 | 1.4e-6 | 1.6e-6 |

- **Pass 0 and PDE0:** the native converges under step capping, toward the resolved Julia result. The original's pass-0 error (3e-9 X1s) and PDE0 error (1e-5 to 3e-5 of peak) are native step-size errors.
- **Converged native vs J14 PDE0:** the remaining 5.5e-6 (2g3s) is inside Julia's own 1e-13→1e-14 change. The two are consistent, but neither is shown converged below ~1e-6 of peak.
- **Pass 1 and FINAL** do not converge in either code at ~1–3e-7 (X1s peak). Native-only self-changes are the same size as the cross-code differences (§5.7–5.9).
- The native tolerance sweep (`compare/native_selfconv_s.log`; tightening s → 1e-3) and monitor-all (`native_monall.log`, up to 1.8e-4 of peak in PDE0 R2s) are consistent with this. As the supervisor noted, they are numerical-control variants and do not prove original error by themselves. The step-cap convergence above is what establishes the native pass-0/PDE0 error, together with agreement between the two independent solvers.

### 5.6 Gate 10a:103 vs Julia configuration: matched-flag one-factor ablation (`chunk11/ablation/`, `ablation_cmp.log`)

**Config.**
- Driver: `run_ablation.sh`, which calls the same public-map driver `julia_public_sweep.jl`.
- Environment: the same Julia 1.12.6 and project environment (`chunk2/test_env`) as every other Julia run here.
- Factors: reltol {1e-12, 1e-13, 1e-14} × a1 {1e-16, 1e-18} × aex {1e-14, 1e-20}.
- Flags: run once with `--check-bounds=yes` (as Pkg.test and the focused 10a run) and once with default flags.
- Repeats: one fresh-process repeat of the production configuration per flag setting.
- Completion: all 4 processes rc 0; 12 + 12 + 1 + 1 cases complete (8 vector files each); k_switch = 1342 in every pass.

**Gate metrics** are computed exactly as in the test (`compare/ablation_cmp.py`), against the native original, which is bitwise the 10a fixture capture.

**Reproduction and determinism.**
- The production configuration under `--check-bounds=yes` reproduces the 10a test values to every printed digit: 103 = 1.08e-3 / 1.04e-3, 118 = 3.25, 119 = 15.4 / 18.3, 117 = 1.71e-6. The environment therefore matches the test, and the earlier 1.075e-3 vs 1.41e-4 gap is a **compiler-flag effect, not an environment effect**.
- With the same flags, fresh-process repeats are **bitwise** identical. The default-flag run is bitwise identical to the earlier sweep `rt1.0e-12`.
- The determinism block of `ablation_cmp.log` also prints `DIFF max|.|/max` scalars over whole PASS/PDE matrices. **Only its "bitwise" verdicts are meaningful.** Those scalars include the z column (different units, dominating the max) and are not accuracy metrics.
- Per-quantity flag effect, z excluded (`ablation/flag_perquantity.log`, from `compare/ablation_flag_perquantity.py`, cmp.py definitions); identical configuration, check-bounds vs default:

| quantity | rt1e-12 / a1 1e-16 / aex1e-14 (production) | rt1e-12 / a1 1e-18 / aex1e-14 | rt1e-14 / a1 1e-18 / aex1e-20 |
|---|---|---|---|
| PASS0 X1s pointwise (zmax) | 5.5e-9 (z = 2467.8) | 3.7e-10 (2620.3) | 8.5e-11 (2567.2) |
| PASS0 X1s peak | 5.8e-11 | 2.8e-11 | 3.7e-13 |
| PASS0 2s pointwise / peak | 7.4e-7 / 3.6e-10 | 6.3e-7 / 4.4e-10 | 3.3e-8 / 1.6e-11 |
| PDE0 2g3s peak (feedback range) | **1.19e-3** | 4.2e-6 | 1.1e-5 |
| PDE0 R2s peak | 2.7e-4 | 8.6e-6 | 1.3e-6 |
| PDE1 2g3s peak | 1.15e-3 | 4.9e-6 | 1.4e-5 |
| PASS1 X1s pointwise / peak | 4.3e-7 / 1.7e-7 | 3.4e-7 / 1.8e-7 | 4.6e-7 / 1.3e-7 |
| PASS1 2s pointwise | 1.2e-5 | 1.5e-5 | 5.7e-6 |
| FINAL Xe pointwise / peak | 1.9e-6 / 1.5e-7 | 1.0e-6 / 1.5e-7 | 1.0e-6 / 2.1e-7 |

- **Reading of the table.**
  - With a1 = 1e-16, the largest flag-induced pass-0 X1s change sits at z ≈ 2468, inside the region where a1 dominates the X1s weight. The PDE stage turns it into a 1.2e-3-of-peak change in 2g3s.
  - With a1 = 1e-18 that PDE change falls to ~5e-6.
  - Pass-1 X1s (1.3–1.8e-7 peak) and final Xe (~1e-6 pointwise) move by the same amount in every configuration. That matches the native's own 1e-15-input reproducibility (§5.9), so tightening a1 does not remove it.

**max103** = the maximum over the eight 10a:103 metrics (DI1, 2g3s, 2g3d, R2s × PDE0, PDE1); bound 1e-3. Values are against the native original; values against the step-converged native cap1e-6 are in parentheses.

| reltol, aex | a1 = 1e-16, check-bounds | a1 = 1e-16, default | a1 = 1e-18, check-bounds | a1 = 1e-18, default |
|---|---|---|---|---|
| 1e-12, 1e-14 (production a1) | **1.08e-3** (1.07e-3) | 3.04e-4 (3.03e-4) | 5.12e-5 (1.68e-5) | 5.12e-5 (1.19e-5) |
| 1e-12, 1e-20 | 1.46e-3 (1.46e-3) | 4.52e-4 (4.52e-4) | 5.12e-5 (3.93e-5) | 5.12e-5 (2.31e-5) |
| 1e-13, 1e-14 | 3.57e-4 (3.57e-4) | 4.80e-4 (4.78e-4) | 5.12e-5 (7.21e-6) | 5.12e-5 (7.60e-6) |
| 1e-13, 1e-20 | 4.93e-4 (4.95e-4) | 2.86e-4 (2.87e-4) | 5.12e-5 (3.33e-6) | 5.12e-5 (6.40e-6) |
| 1e-14, 1e-14 | 3.02e-4 (3.01e-4) | 2.35e-4 (2.34e-4) | 5.12e-5 (1.72e-5) | 5.12e-5 (4.79e-6) |
| 1e-14, 1e-20 | 2.35e-4 (2.35e-4) | 2.47e-4 (2.47e-4) | 5.12e-5 (6.41e-6) | 5.12e-5 (7.13e-6) |

**One-factor effects, matched flags.**
- **a1, 1e-16 → 1e-18 (decisive).** It improves max103 in **all 12** matched pairs: from 2.35e-4 to 1.46e-3, down to 5.12e-5.
  - With a1 = 1e-18, the value against the original is the same 5.12e-5 in all 12 cases. It is the PDE1 R2s difference, and equals the original's own step error (native original vs cap1e-6 PDE1 R2s: 5.1e-5 peak, §5.5).
  - Against the step-converged native, a1 = 1e-18 gives 3.3e-6 to 3.9e-5.
- **aex, 1e-14 → 1e-20 (no consistent effect).**
  - At a1 = 1e-16 it gets worse in 4 of 6 pairs (e.g. rt1e-12: 1.08e-3 → 1.46e-3 with check-bounds; 3.04e-4 → 4.52e-4 with default flags) and better in 2.
  - At a1 = 1e-18, against cap1e-6, the changes are within a factor of 2.5 in both directions.
- **reltol (no monotone effect at a1 = 1e-16).** With check-bounds: 1.08e-3 / 3.57e-4 / 3.02e-4. With default flags: 3.04e-4 / 4.80e-4 / 2.35e-4. At a1 = 1e-18, against the original, it has no effect.
- **Compiler flag.** At a1 = 1e-16 the flag changes the production value 3.6× (1.08e-3 vs 3.04e-4). At a1 = 1e-18 it changes nothing against the original, and at most 2× against cap1e-6.

**Other gates in the same 24 runs** (against the original):
- 10a:117 passes in all of them (9.4e-7 to 3.0e-6).
- 10a:118 fails in all (X1s 2.4 to 6.3 units).
- 10a:119 fails in all (2s 11 to 23, 2p 12 to 29 units).

This is consistent with §5.7–5.9: no tolerance setting fixes the pass-1 gates, because the original itself does not meet them.

**Mechanism.**
- **Source-derived and measured** (`ablation/a1_budget.log`):
  - The Rodas weight of X1s is a1 + reltol\|X1s\|, and X1s is the neutral fraction: 7.4e-10 at z = 3000, 2.6e-8 at the PDE start z = 2500, 6e-6 at z = 2000.
  - With a1 = 1e-16 and reltol 1e-12 the weight is absolute-dominated for all z > 1813. The effective relative tolerance is up to 1.3e-7, and 3.8e-9 at z = 2500.
  - With a1 = 1e-18 these become z > 2141, up to 1.3e-9, and 3.8e-11.
  - The PDE stage is hypersensitive to X1s: 1e-9 white noise on X1s moves 2g3s by 1.4e-2 of peak (`mechanism/pde_response.log`).
- **The same applies to He1s**, which is controlled by a1 and small at high z. It was not separately ablated.
- **The excited-state budget.** aex = 1e-14 lies above the excited populations (2s ~1e-15 to 2e-14; 3s ~1e-17 to 1e-29), and Julia–Julia excited-state differences are near-white noise (`compare/pairwise.log`). That is a source-derived weakness. **It was not measured to affect 10a:103**, and this document no longer claims it does.
- aex ≤ 1e-24 fails (Rodas5P `Unstable` at z0 = 3000; `julia_sweep/sweep_aex.log`).

**Cost.** These are single wall times from the ablation logs, not BenchmarkTools, and exclude the first case of each process, which includes compilation.

| config | default flags | check-bounds |
|---|---|---|
| production (rt1e-12 / a1 1e-16 / aex1e-14) | 2.7 s (earlier sweep, bitwise-identical results) | warm not measured (70.9 s and 73.4 s include compilation) |
| minimal passing change (rt1e-12 / a1 1e-18 / aex1e-14) | 3.4 s (≈1.3× production) | 4.6 s |
| rt1e-13 / a1 1e-18 | 5.6–5.8 s | 8.0 s |
| rt1e-14 / a1 1e-18 | 8.6–8.7 s | 11.9 s |

**Superseded earlier sweep** (`compare/gate103_by_config.log`): default flags, with a1 changing jointly with reltol in the aex1e-20 rows (1e-16 / 1e-16 / 1e-17 / 1e-18 for rt1e-11 to 1e-14). It is kept as a record. Its "monotone convergence with aex = 1e-20" was the **a1** change, and its "7× compile-flag" comparison mixed a check-bounds test value with a default-flag sweep.

### 5.7 Pass 1 is driven by PDE0 (diagnostic cross-injection; `mechanism/df_inject/`, `compare.log`)

Diagnostic only: Julia pass 1 (public `recombination_pass`, rt1e-14/aex1e-20) is fed native PDE0 outputs.
- **Common PDE0 input:** Julia pass 1 vs native pass 1 (step-capped native) differ by 1.8e-10 in X1s peak. Fed the ORIGINAL native PDE0, Julia vs the original pass 1 differ by 4.2e-9 (the original's own step error). **The Julia pass-1 integration is faithful.**
- **Forced response:** feeding the PDE0 of the capped native with 1e-15 H noise, instead of without it, changes Julia pass 1 by exactly the native end-to-end response, digit for digit: X1s 2.32e-7 peak, 5.93e-7 pointwise at z = 1318.9; 2s 1.43e-5 pointwise. The same holds for the profile from z = 1700 to 800.
  - That PDE0 input change is only DI1 1.3e-9, 2g3s 4.9e-7, 2g3d 3.5e-7 and R2s 1.5e-6 of peak.
- **Shape matters more than size** (`pass1_conditioning.log`): a smooth uniform 1e-6 scaling of each PDE0 output moves pass-1 X1s by only 2e-10 to 5e-9 (total feedback effect: X1s 7.6e-3 peak). The ulp-induced change of similar size moves it about 50× more. So pass 1 responds strongly to the non-smooth part of PDE0 errors.

### 5.8 No unstable mode in pass 1 (`mechanism/pass1_slow_mode.jl/.log/.txt`)

- **Method:** slow-manifold rate λ_eff = J_ss − J_sf J_ff⁻¹ J_fs for X1s, with rho and the excited states eliminated, on the diagonally scaled Jacobian. It is evaluated at the resolved Julia pass-1 states, with and without feedback.
  - Float64 vs 256-bit agree to 15 digits.
  - The raw-eigenvalue attempt (`pass1_eigen.log`) is a NEGATIVE result: the scaling is too poor (~1e18 per unit z).
- **Result:** λ_eff > 0 (stable toward decreasing z) everywhere in 900–1680; feedback changes it by < 4%.
- So the exponential growth of the pass-1 response, from 1e-15 at z = 1700 to 2e-7 at z = 1300, is **forced** by the PDE0 input differences. It is not an instability of the equations.

### 5.9 Native reproducibility under 1e-15 input noise (`compare/native_hpert.log`, `hpert_dzcap.log`)

H(z) is multiplied by (1 + 1e-15 noise), with two seeds, in the original and in the step-capped (1e-5) native. All runs are complete.

| quantity | original, seed1 / seed2 | cap1e-5, seed1 / seed2 | 10a gate |
|---|---|---|---|
| PASS0 X1s peak | 8.7e-12 / 8.8e-12 | 8.3e-13 / 3.1e-12 | — |
| PDE0 2g3s peak | 6.9e-6 / 1.0e-5 | 4.9e-7 / 6.7e-7 | — |
| PDE0 R2s peak | 1.7e-5 / 2.5e-5 | 1.5e-6 / 1.6e-6 | — |
| **PASS1 X1s pointwise** | 3.4e-7 / 4.2e-7 | 5.9e-7 / 4.2e-7 | < 1e-7 (1 unit) |
| **PASS1 2s pointwise** | 7.9e-6 / 9.4e-6 | 1.4e-5 / 2.6e-5 | < 1e-5 (10 units) |
| **PASS1 2p pointwise** | 1.42e-5 / 1.04e-5 | 2.6e-5 / 3.2e-5 | < 1e-5 |
| FINAL Xe pointwise | 1.1e-6 / 1.0e-6 | 1.6e-6 / 1.8e-6 | ≤ 1e-5 |

- Step capping removes about 10× of the original's PDE0 sensitivity, which was step-sequence noise. The pass-1 response stays.
- **The native, compared with itself, fails the 10a:118 X1s gate in all four cases and the 2p gate in all four.**
- **Not established:** whether exact arithmetic would show the same pass-1 floor. Under 1e-15 noise the capped native's pass-0 response (8e-13) is still about 1000× the naive δH/H expectation, so residual solver noise is likely still present. The claim is scoped to the measured Float64 native reproducibility at these settings, not to an intrinsic floor.

### 5.10 Gate 10a:82 in native-only terms (`compare/gate82_scaled_native_selfref.log`, `mechanism/pde_response.log`)

- **At identical native rows**, the Julia PDE stage vs the native stage differ by 2–3.6e-7 of peak. The Julia stage's response to a native input change (s1 → s1e-3 rows) equals the native's to 3 digits (1.23e-5 / 1.08e-5 / 3.27e-5).
- **In the 10a:82 `:scaled` metric**, the native compared with itself under 1e-15 H noise reaches 2.8e-5 to 1.6e-4 (PDE1 2g3s, seed 2), and 5.1e-4 for the capped native. Original vs cap1e-6 gives 1.4e-4.
- This is an end-to-end response, not a stage-level one at identical input, so it is a scale reference, not a like-for-like bound. It shows that 1e-4 of the residual scale lies at the native's own Float64 reproducibility.
- The 3.6e-7 stage difference at identical input is attributed, as a **hypothesis**, to different floating-point evaluation order and libm, amplified by the κ ~ 2e8 cancellation. That is consistent with the fixed 7c bound 4 eps κ = 1.8e-7 on raw core differences. It was not proven line by line beyond the 7a–9a stage tests.

### 5.11 Composed derivative references (Julia only; `mechanism/jacobian_config.log`, `chunk10/diag10b_big.log`, `diag10b_big_nbscale.log`, `probe_link_pde_big.log`)

- Float64 FD vs ForwardDiff elasticity: 4.2e-3 / 2.5e-3 / 2.9e-3 at h = 1e-2 / 1e-3 / 1e-4 (10b configuration), and 4.2e-3 / 3.0e-3 / 2.0e-3 (resolved configuration). There is no h² scaling.
- Float64 primal vs 128-bit BigFloat (`setprecision(BigFloat, 128)`, `chunk10/diag10b_big.jl`) on the frozen-decision map: 4.4e-6 relative.
- 128-bit FD converges (h = 1e-6 and 1e-8 agree) and agrees with Float64 ForwardDiff to 5.6e-6 (direction), and 4.9e-6 per output (nbscale).
- So the 10b gates compare against a reference that is noise-limited at about 4.4e-6/h. This is a test-design issue, scoped to the frozen-decision map. It says nothing about the native.

## 6. Negative results (kept; never used as comparison cases)

- `sweep/s1e-4`: native tolerances × 1e-4 ran over 11 min CPU without finishing pass 0 and was stopped. Hypothesis: the Newton test with tolJac = rtol/2 ≈ 5e-12 cannot be reached in Float64.
- `sweep/strict_s1`, `strict_s1e-2`, `strict_monall_s1e-3` (v1 `COSMOREC_DIAG_STRICT`, rejecting when fac < 1): stalled (> 14.5 min) and were stopped. The original controller is not designed for rejection (its estimate mixes orders and history), so this variant is not a usable reference.
- `native_scratch/run_monall.log`: failed launcher (missing script; my `pkill -f` pattern also killed my own shell). The monall cases were produced by the sweep loop instead.
- Julia aex ≤ 1e-24 (all reltols tried): `Unstable` at z0 = 3000.
- `mechanism/pass1_eigen.log`: raw-eigenvalue stability analysis, unusable because of scaling (see §5.8).

Each stopped native case has a `NEGATIVE_RESULT.txt` in its directory.

## 7. Confirmed vs hypotheses

**Confirmed (native and Julia evidence, common inputs):**
1. The original error controller never rejects steps (source plus trace).
2. The pass-0 X1s step at z ≈ 938.7 is a native step error at the original-model `eps_A_effective` discontinuity. The step-capped native converges to the Julia solution, to 7.4e-13.
3. The native pass-0 and PDE0 errors in the original run are step-size errors; they converge under native step capping.
4. 10a:103 is caused by the Julia `a1 = 1e-16` absolute tolerance of rho/X1s/He1s. This is a matched-flag one-factor result: a1 = 1e-18 passes in all 12 reltol/aex/flag cases. aex has no consistent effect, and the earlier aex attribution is refuted (§5.6).
5. Julia pass-1 integration and feedback are faithful at common PDE0 input (1.8e-10).
6. Pass-1 sensitivity is forced by PDE0 differences; the X1s slow mode is stable.
7. The native fails the pass-1 gates (10a:118/119) against itself under 1e-15 H noise and against its own step-converged run.
8. The 10b FD reference is noise-limited (frozen map).

**Hypotheses (not proven):**
- The pass-1/2 worst steps at z = 497.57 come from the feedback cut at z = 500 (location only).
- The 3.6e-7 stage-level PDE difference comes from evaluation order and libm.
- The pass-1 floor of ~2e-7 X1s is intrinsic in exact arithmetic (residual solver noise is likely; see §5.9).
- `s1e-4` stalls in the Newton test.

**Not tested** (common to both codes, since the Julia stage reproduces the native stage): continuum convergence of the PDE grid and time step, and of the Patterson quadrature. This is model discretization error and does not contribute to original-vs-port differences.

## 8. Our own mistakes (candid)

1. **Test/reference callback tolerances** (not a library default; `SOLVE10`, `test/chunk10a_runmode0_native.jl:27`; `abstol5`, `test/chunk5_helpers.jl:53–60`).
   - **The mistake:** `a1 = 1e-16` gives X1s (and He1s) an effective relative tolerance of up to 1.3e-7 at z > 1813. That causes 10a:103, with an outcome that depends on the compiler flag.
   - `aex = 1e-14` is also a weak budget for the excited states. That is source-derived, but it was not measured to affect any gate.
2. **My earlier misattribution** of 10a:103 to aex, made from a sweep where a1 changed together with reltol and from a mixed-flag comparison. The supervisor caught it, and the matched ablation refuted it.
3. **An earlier wrong claim.** In this investigation I first described the native as "rejecting only when the error exceeds ~0.5^(−order) × tol". That was wrong: it never rejects. The supervisor flagged that claim as uninstrumented, and the trace corrected it.
4. Earlier (already fixed): the ScaledBackground fast path; Float(Dual) promotions; the wrong 5e route claim; harness defects (the wrapper clearing the loaded H; DI1 flag ordering).
5. **Test design.**
   - The 10b Float64 FD reference is noise-limited.
   - The pass-1 gates were set from the native solver tolerances, not from the native's measured reproducibility.
   - 10a:82 uses a scale-normalized metric whose denominator is small where the outputs cancel.
6. The `pkill -f` pattern that killed my own shell (process hygiene, no data loss).

No porting mistake in the physics was found:
- RHS ≤ 1.1e-13 on the probed path;
- PDE stage 2–3.6e-7 of peak at identical input, with an identical response to input changes;
- pass-1 integration at common input 1.8e-10;
- an identical discontinuity at the same z in both codes.

The search was not exhaustive: for example, the He-era (z > 1680) RHS was not z-scanned the same way.

## 9. Recommended next fixes (proposals only; nothing implemented; gate and scope changes need Marco's approval)

1. **Julia solver configuration** (no gate change). Set a1 to 1e-18 in the runmode-0 test/reference callback configuration (`SOLVE10`, and the corresponding 5d configuration if one is shared). The library has no solver default to change; downstream callers choose their own. Applied 2026-10-02, see `docs/A1_FIX_VALIDATION.md`.
   - Measured, matched flags: 10a:103 goes 1.08e-3 → 5.12e-5 (check-bounds) and 3.04e-4 → 5.12e-5 (default). That equals the original's own step error, and the result becomes flag-insensitive.
   - Cost: about 1.3× (3.4 s vs 2.7 s, single default-flag timings), to be confirmed with BenchmarkTools.
   - Tightening reltol to 1e-13 improves agreement with the step-converged native (≤ 7.6e-6) at about 2.1× cost. That is optional, because the gate is already met at 1e-12.
   - Changing aex is not supported by the measurements.
   - Check the 5d/10b/10c AD paths at the new a1. This is a test-callback change that needs approval (approved and applied 2026-10-02).
2. **Pass-1 gates 10a:118/119** (Marco's decision). The native cannot meet them against itself. Options:
   - (a) keep them failing as documented;
   - (b) derive fixed bounds from the native's own measured reproducibility envelope. That would be the max over {ulp seed 1, ulp seed 2, cap1e-6} of native-vs-native: X1s 5.9e-7 pointwise, 2s 2.6e-5, 2p 3.2e-5, from native-only runs and so independent of Julia. A pre-declared factor would be added, never self-widening;
   - (c) peak-normalized metrics (X1s 1.4e-7 to 2.2e-7 peak for Julia vs native, against a native self envelope of 1.5e-7 to 2.3e-7).
3. **10a:82** (Marco's decision): the same choice. The native-only envelope in that metric reaches 1.6e-4 (original) and 5.1e-4 (capped). A peak-normalized stage gate at identical input (measured 3.6e-7) would test the port more sharply.
4. **10b:** replace the Float64 FD reference with a higher-precision FD on the frozen-decision map (BigFloat at 320 s per evaluation, or Double64), or derive the bound from the measured primal noise / h. This is a test-design change that needs approval.
5. **10d:** treat it as a Julia-internal AD-consistency gate. Its 1e-6 bound sits at Float64 noise for this hypersensitive stage (`probe_link_pde_big.log`: Float64 ForwardDiff vs 128-bit 4.6e-4, Mooncake 8.9e-5 for that probe's weighting).
6. **Upstream note** (outward-facing; only if Marco wants it): the unreachable rejection at `ODE_solver_Rec.cpp:1244` and the `eps_A_effective` discontinuity could be reported to the CosmoRec authors. They limit the original's pass-0 accuracy to ~3e-9 in X1s and its PDE0 accuracy to ~1e-5 of peak. Their effect on final Xe was not separated from the pass-1 floor, which is ~1e-6 pointwise. Whether that matters physically was not assessed here.

## 10. File index (all under `A/chunk11/`)

- **native_scratch/:**
  - copies, `build_harness.sh`, `diag.diff`, `diag_v2.diff`, `hashes*.txt`, `source_hashes_*.txt`, `ODE_solver_Rec.cpp.diag_v1_copy`;
  - baselines: `out_unmodified/`, `out_diag_default/`, `out_trace_off/`, `out_trace_on/` (`trace.tsv`, `trace_summary.log`), `out_dzcap_default/`;
  - runs: `sweep/`, `hpert/`, `dzcap/{win,glob}_*/` (with `env.txt`, `status.txt`, `time.txt`, `stdout.txt`), `run_dzcap.sh`.
- **rhs_scan/:** `states.txt`, `harness_rhs_scan.cpp`, `build.sh`, `native_scan.txt`, `julia_scan.jl/.txt`, `scan_compare.log`, `scan_jumps.log`, `epsA_check.jl/.txt/.log`.
- **julia_sweep/:** `julia_public_sweep.jl`, per-tag `PASS*/PDE*/FINAL` vectors, `sweep_*.log`.
- **compare/:** `cmp.py`, `dzcap_cmp.py`, and the logs cited above (`node2096_*`, `dzcap_cmp_*`, `hpert_dzcap.log`, `pass1_profile.log`, `gate103_by_config.log`, `gate82_scaled_native_selfref.log`, `native_*`, `julia_*`, `pairwise.log`, `pde_range_pops.log`, `pass0_*`).
- **ablation/:** `run_ablation.sh`, `cb_yes/`, `cb_default/`, `cb_yes_repeat/`, `cb_default_repeat/` (sweep-layout vectors), `*.log`, `status.txt`, `julia_version.txt`, `ablation_cmp.log` (from `compare/ablation_cmp.py`; its matrix-peak DIFF scalars are not accuracy metrics), `flag_perquantity.log` (from `compare/ablation_flag_perquantity.py`), `a1_budget.log`, `supervisor_comparison.log` (independent rerun by the supervisor).
- **mechanism/:** `pde_response.*`, `jacobian_config.*`, `pass1_eigen.*` (negative), `pass1_slow_mode.*`, `df_inject/` (native PDE0/PASS1 extracts, `julia_pass1_from_native_pde0.jl`, `pass1_conditioning.jl`, vectors and logs).
