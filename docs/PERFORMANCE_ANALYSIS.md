# Performance analysis of the runmode-0 forward prediction (diagnosis only; no optimization applied) — 2026-10-02

Notes on scope:
- **Not committed**, per instruction. All profiles, logs and tables are text files in `A/chunk15/`, where `A = cmbcheb_test/local_analysis/cosmorec_differentiability_20260930`.
- No source, test, gate, dependency or tolerance change.
- The analysis tool (JET 0.12.1) lives in an isolated env, `A/chunk15/profenv`. It was installed offline with a private depot first in `JULIA_DEPOT_PATH`; nothing was downloaded, and the shared depot was not written.

**Workload.**
- `recombination_history_diffusion` at the fiducial cosmology, outputs X_e and T_m on the 10000-node grid.
- Caller-supplied Rodas5P callback: reltol 1e-12, a1 1e-18, aex 1e-14. This is the test/reference callback, not a library default.
- Julia 1.12.6, one thread, default flags.

## Current status (2026-10-03): candidate 1 main-reviewed, not committed

**Main-reviewed: correctness/AD/numerical identity and allocation reduction validated; wall time: ONE main-verified, hardware-conditioned matched pair under the predeclared pinned-frequency/thermal/ps-proxy protocol: a matched thermally throttled laptop measurement, not a universal certified algorithmic speedup, statistical CI or unthrottled peak (FULL median ×0.766, 23.39 % less elapsed, 1.305× throughput).** This covers candidate 1 only; it is not a project-wide green status. Main's independent focused 5a native run, including the allocation regression: 1968/1968 pass, exit 0 (`A/chunk16/main_5a_verification.log`).

The diagnosis below (§1–§5) is preserved unchanged as the pre-patch record, including the candidate-1 hypothesis and its unmeasured "~5–20 %" planning estimate.
Full results: `docs/PERFORMANCE_SPECIALIZATION_RESULTS.md`; artifacts in `A/chunk16/`.

- **Patch:** 6 function signatures (2 in `FcnEffective.jl`, 4 in `HIAbsorption.jl`) across 7 changed source lines (the `hi_absorption_rhs!` signature spans 2 lines). `dp_fallback::F = nothing … where {F}` in `fcn_effective!` and `fcn_effective`; `fallback::F = nothing … where {F}` in `dp_correction`, `hi_abs_singlet`, `hi_abs_triplet` and `hi_absorption_rhs!`. No arithmetic, physics, grid, tolerance or clamp change.
- **Numerics:** bitwise identical to the HEAD copy:
  - probe values (RHS incl. off-table fallback states, no-fallback RHS, Jacobians, d/dz, mixed Dual, `fcn_effective`);
  - the public X_e/T_m histories at the fiducial cosmology and the 4 chunk14 cosmologies.

  So no CAMB rerun is needed, and the Boltzmann validation carries over unchanged, including the pre-existing `ob` T_m FAIL.
- **AD:** prepared DI-Mooncake vs ForwardDiff ≤ 3.4e-16, bitwise equal to baseline. Focused 3c/3d/3e/5a native + AD suites all pass (5a native 1968 = 1966 + 2 new).
- **Runtime proof:** the fallback slot of the compiled keyword bodies, read by argument name, was `Function` (`fcn_effective!`, `hi_absorption_rhs!`, `hi_abs_singlet`, `hi_abs_triplet`) or `Any` (`dp_correction`) on the RHS path at baseline. After the patch it is concrete everywhere: the factory closure or `var"#42#43"{Nothing, DPescModel, Bool}`.
- **Allocations (measured):**

  | | memory | allocation events |
  |---|---|---|
  | full forward call | ×0.631 (3.97 → 2.50 GB) | ×0.323 (56.3 M → 18.2 M) |
  | 12-state RHS | 3296 → 1440 B | 77 → 20 |
  | 7-state RHS | unchanged | unchanged |
- **Regression test:** new unconditional allocation-equality test in `test/chunk5a_ode_rhs_native.jl`. It fails on the HEAD copy (`3296 == 1440`) and passes with the patch.
- **Full suite:** root `Pkg.test` exit 1, 50814 passed / 6 failed. That is baseline 50812/6 plus 2; the 6 failures are the known ones at the same sites with identical evaluated values.
- **Wall time (rerun 2026-10-03, main-verified single pair):** warm-then-cool in-process acquisition, both processes pinned to CPU 8, predeclared rule unchanged: FULL median 3085.5 → 2363.9 ms (×0.766; 23.39 % less elapsed, 1.305× throughput), min ×0.759, non-overlapping samples, outputs byte-identical. Matched thermally throttled laptop measurement (both heated 57/59 → 100/97 °C during the batch); not a certified algorithmic speedup, CI or peak. Limitations: one pair; ps lifetime-average foreign criterion; one-sided interval contention log. The first acquisition (historical) established no timing. Details: `docs/PERFORMANCE_SPECIALIZATION_RESULTS.md` §6b.
- **Remaining:** candidates 2–7 untouched; the composed-map ForwardDiff Jacobian benchmark was not run.

## ODE callback optimizations (2026-10-03, approved; one step at a time)

Details: `docs/ODE_CALLBACK_OPTIMIZATIONS.md`. Baseline for each step: the current validated working tree (incl. the `F` patch), not `c5a12b9`.
- **Step 1** — explicit `ODEFunction{true, SciMLBase.FullSpecialize}` / `ODEProblem{true, SciMLBase.FullSpecialize}` in the caller reference callback (test helpers only; the library has no solver default).
  - Runtime `FunctionWrappersWrapper` removed; per-solve solver statistics identical.
  - Public X_e/T_m (5 cosmologies) and the 10b 4-parameter ForwardDiff Jacobian byte-identical; focused 5a–5e and 10a unchanged (10a's 3 known failures identical).
  - New regression test: fails on baseline, passes now.
  - One accepted matched thermally throttled pair: FULL median 2325.9 → 2163.9 ms (×0.930, 6.97 % less elapsed), allocations unchanged.
  - Main review passed.
- **Step 2** — per-solve prepared state-Jacobian cache (`PreparedJac5`, live p/z, fixed chunk, counted oracle fallback) in the caller callback.
  - Public X_e/T_m (5 cosmologies), the 10b Jacobian and all solver stats byte-identical / identical; 0 fallbacks over the public call and under Dual solves.
  - New local J and dJ (y, z, hscale) and prepared-Mooncake-through-cache tests pass; regression tests fail on Step 1.
  - One accepted matched pair: FULL median 2136.6 → 2046.2 ms (×0.958); memory ×0.767. One outlier sample per side, so the ranges overlap.
  - Main-accepted; callback batch root `Pkg.test` 50876 passed / the same 6 known failures.
- **Step 3, increment 1** — private primal RHS workspace (`src/RHSWorkspace.jl`, shared private cores; public allocating path unchanged as the oracle) via the per-solve `RHSWS5` callback.
  - Public outputs, stats and the 10b Jacobian byte-identical; 0 fallbacks.
  - Per-call 0 B (was 1440/1040 B).
  - One accepted matched pair: FULL median 2085.8 → 1886.8 ms (×0.905); memory ×0.512; allocations ×0.258. One outlier sample per side, so the ranges overlap.
  - Main-accepted: a single hardware-conditioned pair, with the GC-outlier caveat.
- **Step 3, increment 2** — buffered state Jacobian (`WSJac5`, family tag `Tag{BufRHS5,V}`, default check, live p/z) in the caller callback.
  - Outputs, stats and the 10b Jacobian byte-identical; 0 fallbacks.
  - Jacobian kernel 0 B.
  - One accepted matched pair: FULL median 1813.5 → 1710.6 ms (×0.943); memory ×0.509; allocations ×0.624. Non-overlapping ranges (one base GC outlier).
  - The off-table Float64 second-order Mooncake/ForwardDiff difference is pre-existing and precision-sensitive (design doc §8); the gate there is no-new-regression parity.
  - Main-accepted.
- **Batch closed (2026-10-03):**
  - final root `Pkg.test` exit 1, **51126 passed / the same 6 known failures** (identical values), 0 errored or broken;
  - whole prediction ≈ 1.711 s median (hardware-conditioned matched pairs), 501 MB / 2.80 M allocations, primal RHS and Jacobian 0 B per call;
  - the Float64 off-table higher-derivative limit is unchanged; not a fully green project; all code uncommitted.
- **CMB tolerance frontier (2026-10-04, research only; `docs/ODE_TOLERANCE_CMB_FRONTIER.md`):**
  - opt-in `reltol 1e-11 / a1 1e-12 / aex 1e-14` passes STRICT (lensed TT/EE/TE Gaussian ΔChi² ≤ 0.01 and per-column criteria) on five cosmologies through ell 9000, with 6.06× fewer RHS evaluations;
  - the observed ~4.4× wall-time pair is UNMATCHED (not certified);
  - a12 is the conservative alternative;
  - nothing adopted.
- **CMB 1e-4 relative target through ell 10000 (2026-10-04, research only; `docs/ODE_TOLERANCE_CMB_1E4_L10000.md`):**
  - PRACTICAL tier (positive spectra pointwise ≤ 1e-4, zero-safe TE) met on five cosmologies under both flags by `reltol 1e-9 / a1 1e-10 / aex 1e-14` (~15.3× fewer RHS evaluations; at the boundary, 8.3e-5/8.7e-5) and, with an 8e-5 guard band, by `1e-9 / 1e-11` (~11.6×);
  - the literal raw-TE goal is met by none of the screened profiles;
  - timing ENVIRONMENT BLOCKED;
  - nothing adopted.

## 1. Measured costs (BenchmarkTools; clean isolated run, `A/chunk15/bench_clean/stage_bench.txt`, Julia rc 0)

Every benchmark is warm, `$`-interpolated, evals = 1. Stage inputs are produced by the same calls as the public function, and their outputs reproduce it bitwise.

| item | n | min | median | memory | allocations |
|---|---|---|---|---|---|
| **full forward call** | 5 | 3.527 s | **3.575 s** | 3783.8 MiB | 56,314,010 |
| ODE pass 0 (no feedback) | 7 | 1.098 s | 1.102 s | 1188 MiB | 18.86 M |
| ODE pass 1 (feedback 0) | 7 | 1.079 s | 1.094 s | 1172 MiB | 18.52 M |
| ODE pass 2 (feedback 1) | 7 | 1.087 s | 1.091 s | 1184 MiB | 18.71 M |
| PDE stage (whole, `hi_diffusion_stage`) | 20 | 78.5 ms | 86.6 ms | 118 MiB | 102,719 |
| …of which `hi_pde_corrections` (march + 199 × 4 integrals) | 20 | 74.1 ms | 83.6 ms | 112 MiB | 71,506 |
| …of which `hi_pde_march` (200 steps) | 20 | 47.3 ms | 48.1 ms | 1.3 MiB | 2,380 |
| …of which `hi_pde_coefficients` | 20 | 4.7 ms | 4.9 ms | 4.0 MiB | 30,845 |
| Recfast tail | 50 | 4.4 ms | 7.2 ms | 11.4 MiB | 76,229 |
| output rows + grid assembly | 50 | 0.82 ms | 0.99 ms | 0.8 MiB | 96 |

- **Shares** (stage medians vs the full median): the three ODE passes are ≈ 3.29 s (**≈ 92%**), the two PDE stages ≈ 0.17 s (≈ 5%), and the tail plus assembly ≈ 0.2%.
- **Memory and allocations are cumulative per call**, not peak RSS (peak process RSS was 1.7 GiB). They are identical to the earlier matched benchmark (`A/chunk12`: 3784 MiB).
- **Timing variance.** The earlier median was 3.244 s. Memory and allocation counts are identical and no flags or code changed, so the 10% timing difference is treated as machine-state variance, **not** a regression.
- **Discarded runs.** An earlier in-driver benchmark (median 6.63 s, `run1/perf_driver.log`) overlapped a foreign precompilation burst from another session's AbstractCosmologicalEmulators benchmark env (load average 14), and is discarded as evidence. The `@elapsed` stage timings in `run1/perf_driver.log` are exploratory instrumentation only.

**Kernels** (same clean run; states from pass 0):

| kernel | min | median | bytes | allocations |
|---|---|---|---|---|
| `recombination_ode!` 12-state (He on, z = 2705.9) | 2.01 μs | 2.80 μs | 3296 | **77** |
| `recombination_ode!` 7-state (He off, z = 1033.7) | 0.78 μs | 0.92 μs | 1040 | 14 |
| `recombination_ode!` 7-state with feedback | 0.83 μs | 1.00 μs | 1040 | 14 |
| `jac5!` (caller callback: unprepared `ForwardDiff.jacobian`) 12-state | 4.49 μs | 5.89 μs | 17936 | 87 |
| `jac5!` 7-state | 1.80 μs | 2.30 μs | 8064 | 24 |
| `tgrad5!` (caller callback) 7-state | 1.14 μs | 1.31 μs | 1888 | 18 |
| `hi_pde_rhs_coefficients!` | 170.5 μs | 176.2 μs | 0 | 0 |
| `hi_pde_step!` | 216.5 μs | 223.9 μs | 0 | 0 |
| `hi_pde_integrals` (one output) | 103.6 μs | 196.6 μs | 577392 | 319 |

**Solver work per pass** (`run1/solver_stats_summary.txt`; stats callback bitwise-identical to the public call):
- About 232,000 RHS evaluations, 29,000 Jacobians, 29,000 W/LU factorizations, 29,000 accepted and about 150 rejected steps.
- The 12-state segment makes 77% of the RHS calls.
- The four 12-state blocks from z 1869 to 1672 (the helium-recombination tail, 4 of 36 blocks) make **≈ 46% of all RHS calls**. The first 7-state block after the helium switch makes 15%.

**Cold costs, reported separately** (`run1/perf_driver.log`):
- 10.5 s for package load plus test-helper includes. That includes native table loading and fixture parsing, so it is not library-only.
- 51.5 s for the first forward call, including compilation.

## 2. Where CPU time goes (sampling profile; `run1/cpu_profile.jls`, postprocessed without recollection)

- **Collection:** 3 warm calls, delay 0.5 ms. 26,305 samples lie inside the forward call and form the denominator.
- **Excluded:** about 26,300 idle listener/other-thread samples (the all-task "50% utilization" view) are not counted as work.
- **Caveat:** the collection window partly overlapped the foreign precompilation, so relative shares are used, not absolute times.

**Non-overlapping context partition** (`run1/cpu_context_split.txt`; first matching rule):

| context | share |
|---|---|
| plain Float64 RHS (`recombination_ode!` via the solver) | **52.8%** |
| `jac5!` (caller Jacobian: RHS on Duals) | 19.0% |
| solver internals (Rodas5P stepping, interpolation, LU) | 12.6% |
| `tgrad5!` (caller time-gradient: RHS on Duals) | 9.5% |
| PDE stages | 5.3% |
| other (tail, assembly) | 0.8% |

**Non-overlapping leaf-first partition by component** (`run1/cpu_postprocess.txt` §2):

| component | share |
|---|---|
| H-I absorber of HeI photons (`HIAbsorption.jl` / `DPescCoh.jl`) | **25.5%** |
| GC | **16.2%** |
| rate-table interpolation (`RateTable.jl`: `get_rates`, `_setup`, `_A`, `log_qnl_qe`) | 15.7% |
| Rodas5P / OrdinaryDiffEq | 12.7% |
| other RHS | 8.4% |
| ForwardDiff machinery | 7.7% |
| helium rates (`HeRHS.jl`) | 6.7% |
| cosmos splines | 3.6% |
| PDE (coefficients 1.2, step 0.6, integrals 0.1, setup 0.1) | 2.0% |
| LU | 0.8% |

Inclusive (overlapping) shares are listed separately in `cpu_postprocess.txt` §1 and are **not** summed; for example, `hi_absorption_rhs!` is 34.3% inclusive.

**Hottest source lines (inclusive):**
- `FcnEffective.jl:93` → `hi_absorption_rhs!` (`HIAbsorption.jl:371`), 34%;
- `hydrogen_rhs!` (`RHS.jl:121`), 18.7%;
- `get_rates` (`RateTable.jl:215`), 12.1%;
- `HIAbsorption.jl:378` / `:385` (the `hi_abs_singlet` / `hi_abs_triplet` calls), 14.2% / 12.0%;
- `dp_lookup` (`HIAbsorption.jl:291`), 7.7%;
- `helium_base_rhs!` (`HeRHS.jl:200`), 8.7%.

**Dispatch and GC leaves by context** (% of forward samples):
- GC: 7.0 in the plain RHS, 5.6 in `jac5!`, 1.2 in `tgrad5!`, 0.4 in PDE.
- Dynamic dispatch / method-cache lookup (including `int64hash`): 3.6 in the plain RHS, 1.0 in `jac5!`, 0.5 in `tgrad5!`.
- The dispatch callers are almost exclusively `hi_absorption_rhs!` lines 378–391.

## 3. Allocations (Profile.Allocs, sample_rate 1e-4, one warm call; `run1/allocs_profile.txt`)

**Totals:** 5,750 sampled events, estimating ≈ 57.5 M events and ≈ 3.81 GiB per call. BenchmarkTools measured 56.3 M allocations and 3783.8 MiB (= 3.695 GiB) per call.
- The two approximately agree, within sampling error: events ≈ 2% higher and bytes ≈ 3% higher in the estimate. They are **not** identical.
- Sampling is **uniform per event, not byte-weighted**, so these are statistical estimates.

**By site (count share):**
- `HIAbsorption.jl:378`: 24.2%;
- `HIAbsorption.jl:385`: 20.1%;
- lines 379–391 together: about another 25%;
- `get_helium_rates` (`HeRHS.jl:130`): 7.0%;
- `fcn_effective!` (`FcnEffective.jl:78`, `:84`): 3.5% and 3.2%;
- `ode_unpack` (`RecombinationODE.jl:60`): 3.3%;
- `recombination_rhs!` (`RecombinationODE.jl:84`): 3.1%;
- `get_rates` (`RateTable.jl:215-217`): 8.3% together;
- `jac5!`: 1.1% of events but 12.4% of sampled bytes.

**By type:**
- boxed `Float64`: 41.5%;
- `Memory{Float64}`: 13.0%;
- `Vector{Float64}`: 10.2%;
- **heap-allocated immutable structs** `HIAbsConstants`, `HIAbsLine`, `DPTable`, `BitotSeries`, `FcorrSpline`: about 15% together;
- keyword tuples `@NamedTuple{f_t, fallback::<closure>, line::HIAbsLine, c::HIAbsConstants}` and `Tuple{Float64, <closure>, HIAbsLine, HIAbsConstants}`.

## 4. Static analysis, and why it disagrees with the profile on the absorber

**JET 0.12.1 `@report_opt` at real states** (`static/jet_report.txt`):
- **0 reports** for `recombination_ode!` (12-state, 7-state, 7-state with feedback), `tgrad5!`, `hi_pde_rhs_coefficients!`, `hi_pde_step!`, `hi_pde_integrals` and `hi_pde_coefficients`.
- `@report_call` on the 12-state RHS: 0.
- `@code_warntype` (via explicit `flag_He = false/true` wrappers, so the analyzed path is the executed one) infers `Vector{Float64}` (`static/code_warntype_rhs.txt`).

**`jac5!` has 865 reports.** The caller's unprepared `ForwardDiff.jacobian(closure, u)` selects the chunk size at runtime, giving `JacobianConfig(f, x, ::Chunk)` and `Array{Dual{…,_A}}` with an uninferred `_A` at the outer level.
- **This does not prove that every Dual RHS evaluation is dynamically dispatched.** ForwardDiff passes the `JacobianConfig{Tag,T,N}` through an inner function barrier, which can restore concrete specialization at runtime.
- The captured runtime instances show exactly that: concrete `Dual{Tag{jac5!…}, Float64, 7}` and `Dual{…, Float64, 12}` specializations of the RHS chain exist (28 each, `static/specializations.txt`).
- The `jac5!` reports are therefore attributed to **outer config/closure construction and dispatch, plus its allocations** (87 per call).
- The dispatch measured inside the Dual RHS (≈ 1% of samples) is at least partly the **independent** absorber-fallback despecialization below, which applies to Dual element types too. The split between the two causes was not measured.

**The forward call has 25 reports, all low-frequency** (the profile puts them at ≈ 0.05% together):
- block bookkeeping in `recombination_pass` held in `Any` containers (57 blocks per pass);
- the per-pass `Any[]` loop of `recombination_history_diffusion`;
- `get_rates_all` → `lagrange4(x::Any)` in the once-per-stage coefficient setup;
- `NaturalCubicSpline{T} where T` in the PDE integrals;
- the rarely executed `dpesc_coh` / `sobolev_p` off-table fallback.

**Contradiction resolved: runtime despecialization** (`static/specializations.txt`). After a plain 12-state RHS call, the runtime holds the instance

`#fcn_effective!#45 :: (…, Function, Bool, Nothing, typeof(fcn_effective!), Vector{Float64}, Float64, Vector{Float64}, EffectiveBackground{…}, EffectiveModel)`

so the keyword `dp_fallback` is **despecialized to `Function`**. Several `#dp_correction#20` instances have `fallback::Any`. These coexist with concretely typed instances.
- **The mechanism** is Julia's documented heuristic: a method is not specialized on a `Function` argument that is only passed through, not called. `dp_fallback` / `fallback` is only forwarded along `fcn_effective!` → `hi_absorption_rhs!` → `hi_abs_singlet` / `hi_abs_triplet` → `dp_correction`.
- **The consequences match every measurement:**
  - each forwarding call is a runtime dispatch with heap-boxed arguments (the boxed structs and keyword tuples in §3);
  - the returned `dS` / `dT` are `Any`, so `g[…] += dS` (lines 379–391) boxes `Float64`;
  - this explains the 77 allocations of the 12-state RHS, against 14 for the 7-state RHS where the absorber is off.
- JET and `code_warntype` analyze with fully concrete argument types, so neither sees this. **"0 JET reports" therefore does not mean "no runtime dispatch"** for this path.

## 5. Candidates (proposals only; implement one at a time after review, each with a focused numerics + ForwardDiff/Mooncake check and a matched BenchmarkTools before/after)

| # | candidate | where | evidence | rough ceiling (Amdahl, of the 3.57 s call) | risk |
|---|---|---|---|---|---|
| 1 | Force specialization of the pass-through fallback (`fallback::F … where {F}` / `dp_fallback::F`) in `fcn_effective!`, `hi_absorption_rhs!`, `hi_abs_singlet`, `hi_abs_triplet`, `dp_correction` | `src/FcnEffective.jl:68`, `src/HIAbsorption.jl:347-371` | §2–§4: dispatch/method-cache leaves 3.6 + 1.0 + 0.5% of forward samples; about 70% of allocation events at these lines; GC 16% of samples | **Heuristic planning estimate, NOT an upper bound: ~5–20%.** The estimate assumes the removed dispatch leaves (≈ 5%) plus a share of GC roughly proportional to the removed events, plus some boxed-arithmetic overhead. GC cost is not linear in event count, bytes, generations or roots, so the event-share × GC product is not rigorous. The boxed-arithmetic share is unknown. The only measured fact is the affected CPU: the absorber is 34% inclusive and 25.5% of the leaf partition. The real gain must be measured | Low numerical risk: the arithmetic is unchanged, so bitwise-identical outputs are **expected** (to be verified); slightly more compilation. AD: expected unaffected mathematically, but **ForwardDiff and Mooncake must be verified before acceptance** (not measured). The keyword syntax `fallback::F = nothing … where {F}` was verified to work in an isolated Julia 1.12 probe (review); no source patch exists yet |
| 2 | Remove per-RHS-call heap temporaries (`ode_unpack` X vector, `g`, `dXHe` zeros, `RateResult` vectors in `get_rates`, `get_helium_rates`) via stack (`StaticArrays`/tuples) or type-generic caches | `RecombinationODE.jl:60,84`, `FcnEffective.jl:78,84`, `RateTable.jl:215-217`, `HeRHS.jl:130` | §3 (about 25% of events); 14 allocations even in the 7-state RHS | the remaining GC share after #1 (a few %) plus allocation cost | **Moderate for AD:** caches must be element-type-generic (Duals, Mooncake tangents) and free of aliasing in reverse mode; StaticArrays would be a dependency change and needs approval |
| 3 | Caller-side Jacobian: prepared `JacobianConfig` with a fixed `Chunk{N}` (or DifferentiationInterface prep) instead of `ForwardDiff.jacobian(closure, u)` per call; same for `tgrad5!` | `test/chunk5_helpers.jl:44-50` (the **injected test callback**; the library has no default) | `jac5!` 19% + `tgrad5!` 9.5% of samples; 865 JET reports from outer config construction (the inner Dual RHS is concretely specialized behind ForwardDiff's function barrier, §4); 87 allocations per call | removes the outer config construction, dispatch and allocations only; most of the ≈ 28% is the Dual RHS work itself, which remains (affected by #1 instead) | A test/benchmark-callback change only (identical Jacobian expected); downstream callers choose their own |
| 4 | `ODEProblem{true, SciMLBase.FullSpecialize}` in the callback | caller | FunctionWrappers ≈ 2.4% self | ≈ 2% | Caller-side; longer compilation |
| 5 | Solver work concentration: 46% of RHS calls in z 1869 → 1672, 15% in the first post-switch block; 36 integrator restarts per pass (block = 50/200 nodes) | `recombination_pass` | §1 stats | potentially large | **Numerics change** (step sequences, tolerances, block layout); outside this diagnosis; needs a separate accuracy-preserving study and approval |
| 6 | PDE integrals allocate 577 KB per output (80 MB per stage) | `HIPDEIntegrals.jl` | §1 kernels | ≤ 1% | Low; minor |
| 7 | Cold start: about 51 s first-call compilation | package | §1 | latency only | A PrecompileTools workload; separate concern |

**Not addressed:** gradient-route performance. The ForwardDiff Jacobian of the composed map, 6.5 s at a1 = 1e-18 per `A/chunk12`, is *expected* to benefit from #1, since the absorber-fallback despecialization also applies to Dual element types; this is not measured. Mooncake/true-reverse (B2) is blocked by the backend rule, which is separate from primal performance. The 60-minute hybrid experiment was not repeated.

## 6. Commands, runs and terminal states

| step | command / script | Julia exit | output |
|---|---|---|---|
| tool env | `JULIA_DEPOT_PATH=A/chunk15/depot:~/.julia JULIA_PKG_OFFLINE=true julia --project=A/chunk15/profenv -e 'Pkg.add("JET")'` | 0 | JET 0.12.1 |
| profile + bench driver (first attempt) | `julia --project=A/chunk2/test_env A/chunk15/perf_driver.jl A/chunk15/run1` | **1** (`UndefVarError` `cur`, soft scope, `perf_driver.jl:125`, after the profile was saved; the allocation section never ran; log preserved in `run1/time.txt`) | `run1/perf_driver.log`, `cpu_profile.jls`, `cpu_flat*.txt`, `cpu_tree.txt`, `cpu_by_thread.txt` |
| profile postprocess (no recollection) | `profile_postprocess.jl`, `profile_split.jl` on `run1/cpu_profile.jls` | 0 | `run1/cpu_postprocess.txt`, `run1/cpu_context_split.txt` |
| allocation profile | `allocs_profile.jl` | 0 | `run1/allocs_profile.txt` |
| static analysis | `JULIA_LOAD_PATH="@:A/chunk15/profenv:@stdlib" julia --project=A/chunk2/test_env static_jet.jl A/chunk15/static` | 0 | `static/jet_report.txt`, `static/code_warntype_rhs.txt` |
| runtime specializations | `static/specializations.jl` | 0 | `static/specializations.txt` |
| clean BenchmarkTools (stages + kernels; quiet machine, load ≈ 2–3) | `stage_bench.jl A/chunk15/bench_clean` | 0 | `bench_clean/stage_bench.txt` |
