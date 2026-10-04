# Candidate 1: forced specialization of the pass-through fallback callable — results (2026-10-02)

Status: **main-reviewed: correctness/AD/numerical identity and allocation reduction validated; wall time: ONE main-verified, hardware-conditioned matched pair under the predeclared pinned-frequency/thermal/ps-proxy protocol: a matched thermally throttled laptop measurement, not a universal certified algorithmic speedup, statistical CI or unthrottled peak — FULL public call median 3085.53 → 2363.86 ms (23.39 % less elapsed time, 1.305× throughput), min 24.1 % lower (§6b).** Earlier acquisition (§6) established no timing. Implemented in the working tree, not committed. Main's independent focused 5a native run incl. the allocation regression: 1968/1968 pass, exit 0 (`A/chunk16/main_5a_verification.log`). Artifacts are text files in `A/chunk16/`,
`A = cmbcheb_test/local_analysis/cosmorec_differentiability_20260930`. Diagnosis and hypotheses: `docs/PERFORMANCE_ANALYSIS.md` (§2–§5, candidate 1).

## 1. The patch (6 function signatures — 2 in FcnEffective.jl, 4 in HIAbsorption.jl — across 7 changed source lines, since the `hi_absorption_rhs!` signature spans 2 lines; 2 files; no arithmetic, physics, grid, tolerance or clamp change)

The fallback callable was only passed through these functions, so Julia's no-specialization heuristic for `Function` arguments compiled them with
`dp_fallback::Function` / `fallback::Any` (runtime dynamic dispatch, boxed keyword tuples). A keyword type parameter forces specialization:

| function | file | change |
|---|---|---|
| `fcn_effective!`, `fcn_effective` | `src/FcnEffective.jl` | `dp_fallback = nothing` → `dp_fallback::F = nothing … where {F}` |
| `dp_correction`, `hi_abs_singlet`, `hi_abs_triplet`, `hi_absorption_rhs!` | `src/HIAbsorption.jl` | `fallback = nothing` → `fallback::F = nothing … where {F}` |

Public keywords, `nothing` defaults, the factory closure/callables, and the Float64/Dual/BigFloat and off-table fallback paths are unchanged.
Diff: `A/chunk16/patch.diff`. Baseline source state (HEAD `c5a12b9`, sha256 of both files): `A/chunk16/baseline/source_state.txt`.

**Regression test added** (unconditional, normal suite): `test/chunk5a_ode_rhs_native.jl`, testset
"dp_fallback callable is specialized (no extra allocation from the pass-through)". At the in-table 12-state point z = 2500 (fallback never invoked),
it warms both models through one wrapper `alloc5a!(du, r, m) = @allocated recombination_rhs!(du, r.z, r.y, m; flag_He = true)`. It then requires the RHS with the
factory fallback (`RM5`) to equal the no-fallback model's RHS exactly, and to allocate exactly as many bytes.
- Against the HEAD copy: **FAILS**, `Evaluated: 3296 == 1440`, 1967 passed / 1 failed, exit 1 (`A/chunk16/probe2/regress/base_5a_newtest.log`).
- Patched: **PASSES**, 1968/1968.

## 2. How the baseline was obtained without touching the working tree

`git archive HEAD` → `A/chunk16/baseline_pkg/CosmoRec.jl`, with sha256 equal to the saved baseline. It has its own env (the copied test env with only the CosmoRec path changed)
and its own depot first in `JULIA_DEPOT_PATH`, so its compile caches did not go to `~/.julia`. Probe runs print `pathof(CosmoRec)`, which confirms the copy was loaded,
and the test helpers come from the copy. No stash, reset or worktree.

## 3. Numerics: bitwise

| check | result |
|---|---|
| probe values, `%.17g` (`spec_check.jl`): RHS at all pattern-0 5a rows incl. off-table fallback states z = 1500/1200/1000/800 and inactive-absorber states; RHS with `nothing` fallback; ForwardDiff Jacobians; d/dz; mixed Float64-z/Dual-state RHS; `fcn_effective` on all 3d states with and without fallback; DI-Mooncake and ForwardDiff gradients | `cmp` **identical** base vs after (`A/chunk16/probe2/{base,after}/spec_values.txt`) |
| no-fallback throws | only `DPTableDomainError`, at the same 4 off-table states, in both. The probe catch is typed and rethrows anything else, so a new `MethodError` would abort the probe instead of looking like a NaN |
| public X_e/T_m history, fiducial + ob, oc, h70, yhe (chunk14 drivers, same inputs) | all 5 `inj_julia.txt` and the native-injected companions **byte-identical** to `A/chunk14` |
| consequence | no CAMB/Boltzmann rerun needed; the Boltzmann validation (`docs/BOLTZMANN_FORWARD_VALIDATION.md`) carries over unchanged, including the pre-existing `ob` T_m 1.358e-7 FAIL against the < 1e-7 criterion. No gate widened |

## 4. AD

- Prepared DI-Mooncake vs ForwardDiff gradients of w'RHS: max relative difference 1.311e-16 (12-state, z = 2500), 3.390e-16 (12-state off-table fallback, z = 1500), 1.319e-16 (7-state, z = 400). These are bitwise identical to baseline.
- Mixed Dual cases (Dual state with Float64 z, and Dual z) are bitwise identical to baseline.
- Focused suites, each in a fresh process with `--check-bounds=yes` (`A/chunk16/{baseline,after}/focused/`):

| suite | baseline | patched |
|---|---|---|
| 3c H-I absorber native / AD | 8870 / 462 pass | 8870 / 462 pass |
| 3d `fcn_effective` native / AD | 566 / 2203 pass | 566 / 2203 pass |
| 3e DPesc fallback native / AD | 1512 / 281 pass | 1512 / 281 pass |
| 5a ODE RHS native / AD (ForwardDiff vs 256-bit FD, prepared Mooncake) | 1966 / 35 pass | **1968** (+2 new) / 35 pass |

## 5. Runtime proof that specialization changed (structural, not a string filter)

`spec_check.jl` reads the keyword-body argument named `fallback`/`dp_fallback` via `Base.method_argnames` and prints exactly that slot's type from every
compiled `MethodInstance` after the probe workload (`A/chunk16/probe2/{base,after}/spec_runtime.txt`). Counts are numbers of instances:

| keyword body | baseline slot types | patched slot types |
|---|---|---|
| `fcn_effective!` | **`Function` ×8 (1 Float64-path, 7 Dual-path)**, factory closure ×12, `Nothing`, `typeof(fb3e)` | factory closure `var"#RecombinationModel##0#RecombinationModel##1"` ×12 (5 Float64, 7 Dual), `Nothing`, `typeof(fb3e)` |
| `hi_absorption_rhs!`, `hi_abs_singlet`, `hi_abs_triplet` (each) | **`Function` ×8** + concrete as above | concrete only (same set) |
| `dp_correction` | **`Any` ×8 (1 Float64, 7 Dual)**, `var"#42#43"{Nothing, DPescModel, Bool}` ×14, `Nothing` | `var"#42#43"{Nothing, DPescModel, Bool}` ×14 (6 Float64, 8 Dual), `Nothing` |

No fallback slot is `Function` or `Any` after the patch. (The concrete instances that already existed at baseline come from direct test calls with concrete types; the
`Function`/`Any` instances are the ones the ODE RHS path actually ran through.)

**Allocations per warm RHS call** (`@allocated`, `spec_runtime.txt`): 12-state z = 2500: **3296 → 1440 B**. 7-state z = 400: 1040 → 1040 B. The 7-state call
has helium off; it is unchanged, consistent with the despecialized pass-through being on the helium/H-I absorber path, which the probe shows was the only
`Function`/`Any` slot site.

## 6. Allocations and timing — HISTORICAL first acquisition (unmatched timing; superseded for timing by §6b; allocation figures still valid)

**Protocol.**
- `A/chunk16/guarded_bench.jl` (r1) and `guarded_bench_v2.jl` (r2) cover the full forward call, ODE pass 0, the RHS kernels (12- and 7-state) and the `jac5!` callback Jacobians.
- Inputs and the callback (Rodas5P, reltol 1e-12 / a1 1e-18 / aex 1e-14) are the same as chunk15. Settings: `$` interpolation, evals = 1, FULL ≥ 5 samples, pass 7, kernels 5000. Compilation is excluded.
- Runs are interleaved base r1 → after r1 → base r2 → after r2, each in a fresh process and each preceded by a read-only gate, `wait_quiet[_v2].sh`. Every gate passed QUIET: 120 s at ≤ 75 °C, zero throttle delta, no foreign process > 50 % CPU.
- Each batch logs at start and end:
  - package and hottest-core temperature;
  - mean `scaling_cur_freq`;
  - **deltas** of the package/core throttle counters (never the cumulative values);
  - load;
  - the foreign > 50 %-CPU set (v2 only; v1 recorded it at the end only).
- Pairs are assessed mechanically by `bench_match.py` against `bench_matching_rule.txt`, which was fixed before any timed data was taken. MATCHED requires all of:
  - the same foreign set;
  - |ΔT_pkg| ≤ 5 °C at the start and at the end;
  - mean frequency within 10 %;
  - throttle rates both zero, or within a factor of 1.5.

  Report: `A/chunk16/gbench/match_report.txt`; raw: `gbench/guarded_bench.txt`. No threshold was changed after the data.

**Allocations — measured, independent of thermals** (identical in both rounds):

| batch | memory base → patched | allocation events base → patched |
|---|---|---|
| FULL `recombination_history_diffusion` | 3967610928 → 2503627696 B (**×0.631**) | 56314010 → 18195488 (**×0.323**) |
| ODE pass 0 | 1246156528 → 753675376 B (×0.605) | 18858956 → 6036350 (×0.320) |
| `recombination_ode!` 12-state | 3296 → 1440 B (×0.437) | 77 → 20 (×0.260) |
| `recombination_ode!` 7-state | 1040 → 1040 B | 14 → 14 |
| `jac5!` 12-state | 17936 → 13200 B (×0.736) | 87 → 30 (×0.345) |
| `jac5!` 7-state | 8064 → 8064 B | 24 → 24 |

(Allocation bytes and event counts are not RSS, and the event fraction is not a GC-time budget.)

**Timing (historical, first acquisition) — NOT established under the predeclared rule (environment-limited); see §6b for the main-verified rerun.**
- Every gate passed cool (58–59 °C). Within the ≈ 1 min of compilation and warm-up before the first timed batch, however, a single Julia thread drove the package to 88–100 °C. All long batches (FULL, pass 0) therefore ran in sustained thermal throttling (package throttle rate 77.7–92.6 events/s, all-CPU mean frequency 1.17–2.07 GHz at the batch ends). The ~1 s kernel batches ran in transients.
- Results of the 12 pairs:
  - r1 FULL: **CONDITIONAL**. Criteria 2–4 are met: T 100/100 → 100/98 °C, frequency within 10 %, pkg rate 90.2 vs 77.7 /s. The foreign set at batch start was not observed (v1). Under these conditions it measured min 3093.1 → 2312.6 ms (×0.748), median 3157.6 → 2333.2 ms (×0.739). This is a **matched steady-state throttled throughput** observation, not unthrottled peak performance, and it is conditional on criterion 1.
  - r2 FULL: **UNMATCHED**. ΔT at start 88 vs 95 °C > 5 °C; end frequency 1172 vs 1466 MHz > 10 %. Its timings (×0.742 min, ×0.750 median) are not used.
  - All r1/r2 pass-0 and kernel pairs: **UNMATCHED** (temperature, frequency or throttle-rate criteria); timings not used.
- No speedup figure is claimed. The unmatched ratios point the same way as the CONDITIONAL one, but that is not evidence under the rule. A valid timing needs a controlled host (cooling, or a sustained-load-stable state matched on both sides); changing the governor, fans or jobs requires user approval and was not done.
- Limitation: the logged frequency is the mean over all logical CPUs, not the frequency of the core that ran the benchmark.

**Process notes.**
- v1 counted its own short-lived `ps` child as a foreign 100 % process (a false gate reset in `gate_base_1.log` at 00:02:19), and recorded the foreign set only at batch end. Both were fixed in v2: the own process group is excluded and the set is recorded at start and end. Thresholds are unchanged.
- The v1 chain was stopped during gate r2, with no timed batch running; its log is kept as `gate_base_2_v1_stopped.log`.
- The older un-logged baseline numbers (`A/chunk16/baseline/stage_bench.txt`, FULL median 3509.99 ms) have no thermal record and are not used for ratios.

## 6b. Timing rerun 2026-10-03 (acquisition fix), `A/chunk16/rerun_20261003T023027/` — MAIN-VERIFIED (single pair, scoped)

Requested explicitly by Marco. No source/test/gate/physics/tolerance/dependency change. Earlier data (§6) preserved. The protocol and rule were written before any data, in
`PROTOCOL.txt`; a later addendum narrowed the scope without changing any threshold.

**Acquisition change.**
- Before: the gate ran before a fresh Julia process, so load/JIT/warm-up heated the CPU to 88–100 °C before timing.
- Now (`guarded_bench_v3.jl`; after-r1 used `v3b`, which adds logging only): each process
  1. loads, sets up, and compiles and warms both the target and the BenchmarkTools thunk;
  2. then waits for an in-process QUIET gate, with the same thresholds as before (120 s at ≤ 75 °C, zero throttle-counter delta, no ps-listed foreign process > 50 % CPU);
  3. then times immediately, with no restart in between.
- Both processes were pinned with `taskset -c 8` to the same logical CPU 8 (P-core `core_id` 16, sibling CPU 9, max 4.8 GHz; P-cores cpu0–11 with HT, E-cores cpu12–19).
- The frequency criterion uses the pinned CPU's `scaling_cur_freq`. A 0.5 s helper sampler (own process group) logged package temperature, CPU 8 frequency and throttle counters during each batch.

**Scope (supervisor narrowing).**
- Primary only: the public `recombination_history_diffusion` (fiducial, 10000-node X_e/T_m, Rodas5P reltol 1e-12 / a1 1e-18 / aex 1e-14, default flags, 1 thread), 5 samples, evals = 1, `$`-interpolated.
- base-r1 was stopped while in the gate of the optional 12-state batch, with no timed batch running. Secondaries were skipped. after-r1 ran FULL only. Round 2 was not needed.
- Total time: launch 02:31 to finish 02:55.

**The pair (base = HEAD `c5a12b9` copy, after = candidate 1).** Assessment by the unchanged rule, `match_r1.txt`:

| | base-r1 | after-r1 | criterion | met |
|---|---|---|---|---|
| gate before the batch | QUIET (120 s) | QUIET (120 s) | — | — |
| foreign > 50 % (ps) start / end | none / none | none / none | identical | yes |
| package temperature start → end | 57 → 100 °C | 59 → 97 °C | ≤ 5 °C | yes (2, 3) |
| pinned CPU 8 frequency start → end | 4800 → 4700 MHz | 4800 → 4787 MHz | ≤ 10 % | yes |
| throttle rate, package / core counters | 75.3 / 150.6 /s | 71.6 / 143.3 /s | ratio ≤ 1.5 | yes (1.05) |
| during batch (sampler): package temperature min/median/max | 57 / 97 / 100 °C (n = 39) | 59 / 94 / 98 °C (n = 31) | descriptive | — |
| during batch: CPU 8 frequency min/median/max | 4053 / 4791 / 4800 MHz | 4269 / 4800 / 4800 MHz | descriptive | — |
| all-CPU mean frequency start → end | 1529 → 1682 MHz | 1394 → 1230 MHz | reported, not a criterion | would fail (end) |

**Verdict: ACCEPTED by the unchanged rule; main-verified (results.txt, match_r1.txt, thresholds, pinned CPU 8, output cmp = 0).** Label: matched thermally throttled laptop measurement. Both batches started cool (57/59 °C) and heated during the batch to 100/97 °C, i.e. NOT a strictly steady thermal state at all times; they entered package throttling at near-identical rates, while the pinned core
stayed near its 4.8 GHz maximum. This is not an unthrottled-peak measurement.

| FULL public call | base-r1 | after-r1 | after / base |
|---|---|---|---|
| all 5 samples [ms] | 3085.53, 3074.50, 3081.84, 3144.65, 3128.20 | 2363.86, 2332.82, 2351.69, 2370.56, 2367.05 | — |
| min [ms] | 3074.50 | 2332.82 | 0.759 (24.1 % lower) |
| median [ms] | 3085.53 | 2363.86 | 0.766 (23.39 % less elapsed; 1.305× throughput) |
| GC per sample [ms] | 307.4–340.5 (median 332.9) | 180.3–204.3 (median 197.1) | — |
| non-GC median [ms] | 2778.1 | 2162.8 | 0.779 |
| memory | 3967610928 B (3.695 GiB) | 2503627696 B (2.332 GiB) | 0.631 |
| allocation events | 56314010 | 18195488 | 0.323 |

- **Margins.** The sample ranges do not overlap. The pairwise extremes bound the ratio between 0.742 (min after / max base) and 0.771 (max after / min base). This is the reliable statement for this single matched pair; it is not a population confidence interval.
- **Output check.** The public X_e/T_m output of one warm call in each process is byte-identical (`cmp` of `output_base-r1.txt` vs `output_after-r1.txt`).
- **Cold costs (separate, not warm samples).** Process start to first public result: base 84.1 s, after 45.7 s. These are not comparable: base loads from its own compile-cache depot with a different package path, and the CPU state differed. They are reported only for completeness.

**Limitations.**
- One accepted pair on one host.
- The foreign-load criterion uses `ps` `%CPU`, which is a lifetime average, not current load. after-r1 additionally logged measured contention from `/proc/stat` deltas minus its own process: 0.42–0.49 foreign busy cores in the last gate intervals, 0.52 over the batch. base-r1 (v3) did not log this, so the contention comparison is asymmetric.
- `scaling_cur_freq` is kernel-reported, not APERF/MPERF.
- The 12-state/7-state RHS secondaries were not timed in this rerun; §6 allocation figures apply.

## 7. Full suite (root `Pkg.test`, run once, unconditional)

`COSMOREC_NATIVE_DATA_DIR=… julia --project=. -e 'using Pkg; Pkg.test()'`, 22:59–23:55, log `A/chunk16/fullsuite/pkgtest.log`, exit **1**:
**50814 passed, 6 failed, 0 errored, 0 broken**. Baseline (`A/chunk12/pkgtest_a1_run1.log`): 50812 passed, 6 failed.
- The +2 passes are the new regression assertions.
- The 6 failures are the same known pre-existing ones, at the same sites (`chunk10a_runmode0_native.jl:83,119,120`, `chunk10b_runmode0_ad.jl:42,56`,
`chunk10d_background_link_ad.jl:46`), with **identical printed `Evaluated` values**.
- No new failure; no check was relaxed.

## 8. Remaining limits

- Wall time: ONE main-verified, hardware-conditioned matched pair under the predeclared pinned-frequency/thermal/ps-proxy protocol: a matched thermally throttled laptop measurement, not a universal certified algorithmic speedup, statistical CI or unthrottled peak (§6b; median ×0.766). §6 (first acquisition) established none and is kept as historical record. Limitations: one pair, one host; ps lifetime-average foreign criterion; interval contention logged on the after side only. The allocation reduction, bitwise numerics, AD agreement and the runtime slot types are established.
- 7-state RHS still allocates 1040 B/call and the 12-state 1440 B/call (candidate 2, per-call temporaries); not addressed here.
- The ForwardDiff Jacobian of the composed map (4-parameter, chunk12-style) was not re-benchmarked.
