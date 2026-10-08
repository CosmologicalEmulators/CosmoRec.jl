# ODE caller-callback optimizations (approved 2026-10-03; one step at a time, each behind main review)

Scope: the **caller-supplied** Rodas5P reference callback in `test/chunk5_helpers.jl` (`ODEFUN5`, `solve_phase5`). The library has **no solver default**:
callers that build their own `ODEProblem` get these benefits only if they make the same choice. Plan and diagnostics: `docs/SCIML_ODE_PERFORMANCE_REVIEW.md`.
Update 2026-10-07: the `solve_phase5` defaults are now Rodas4P / reltol 1e-10 / a1 1e-18 / aex 1e-12 / `tstops = znodes` (measured -50% forward config, CMB Tier A verified); the legacy tight config (Rodas5P / 1e-12 / 1e-18 / 1e-14, `tstops = false`) is pinned explicitly at all fixture-gated tests and benchmark call sites, and the step numbers below refer to it. The library itself still has no solver default.
Artifacts are text files outside git in `A/chunk18/`, where `A = cmbcheb_test/local_analysis/cosmorec_differentiability_20260930`. Nothing is committed.

**Baseline for every step** = the CURRENT validated working tree, not `c5a12b9`:
- the uncommitted `F`-specialization patch (`src/FcnEffective.jl`, `src/HIAbsorption.jl`) and its 5a regression test;
- the untracked analysis docs;
- HEAD and `origin/develop` = `c5a12b9`, fetched fresh.

It was snapshotted (sha256, working-tree diff, Julia 1.12.6, package versions) to `A/chunk18/baseline_snapshot/`, without git stash or reset. Step gains are incremental
relative to this baseline, so the earlier ×0.766 of the `F` patch is not counted again. Baseline evidence:
- full `Pkg.test` 50814 passed / the same 6 known failures;
- public X_e/T_m at the fiducial and 4 nearby cosmologies byte-identical to chunk14/chunk16;
- the `ob` T_m 1.358e-7 miss of the < 1e-7 target is pre-existing and unchanged.

## Step 1 — explicit `FullSpecialize` in the reference callback — MAIN-REVIEWED (pass)

**Change** (test helpers only; no src, algorithm, tolerance, `saveat`, norm, physics, clamp or grid change; the Recfast-tail callback is unchanged):

```julia
const ODEFUN5 = ODEFunction{true, SciMLBase.FullSpecialize}(recombination_ode!; jac = jac5!, tgrad = tgrad5!)
prob_phase5(u0, z0, znodes, p) = ODEProblem{true, SciMLBase.FullSpecialize}(ODEFUN5, u0, (z0, znodes[end]), p)   # used by solve_phase5
```

All runmode-0 test/benchmark routes (`SOLVE5B/5C/5D/10`, `solve_phase5` users) go through this one callback. `test/chunk5e_helpers.jl` and
`test/chunk10_reverse_helpers.jl` already used `FullSpecialize`.

Order note: the probe, focused suites and histories below ran with `ODEProblem{true, SciMLBase.FullSpecialize}(ODEFUN5, …)` written inline in `solve_phase5`. It was then moved, unchanged, into the accessor `prob_phase5` for the regression test. After the move, 5b (46/46) and the timing processes' public outputs (byte-identical to baseline) were rerun.

**Runtime proof** (`A/chunk18/step1_probe.jl` → `step1/probe_probe.txt`, exit 0). Same process, same `src`; the old Auto callback is reproduced verbatim:
- **Wrapper:** with Auto, an initialized Rodas5P integrator has `integ.f.f::FunctionWrappersWrapper{…}`. With FullSpecialize, `integ.f.f::typeof(recombination_ode!)`.
- **Bitwise:** the installed `solve_phase5` equals the probe's FullSpecialize solve byte for byte (pass 0).

**Solver work: identical, per solve, across all 108 solves of the public three-pass call:**

| system | solves | RHS evaluations (nf) | Jacobians | factorizations (nw) | accepted | rejected |
|---|---|---|---|---|---|---|
| 12-state | 81 | 535498 | 66624 | 66917 | 66624 | 293 |
| 7-state | 27 | 159174 | 19722 | 19890 | 19722 | 168 |

**Numerics and AD.**
- Public X_e/T_m (fiducial): Auto vs Full byte-identical.
- ForwardDiff 4-parameter Jacobian of the 10b runmode-0 route (18×4, `P0_5E`): Auto vs Full byte-identical, and byte-identical to `A/chunk12/bench_a1_jacobian10b_NEW_a1_1e-18.txt`.
- Histories (chunk14 drivers, fiducial + ob, oc, h70, yhe): all `inj_julia.txt` and native-injected files byte-identical to the baseline (`A/chunk16/after/hist`), so no CAMB rerun is needed.

**Focused tests** (fresh processes, `--check-bounds=yes`, `A/chunk18/step1/focused/`; pass/fail read from the logs):

| test | result | baseline |
|---|---|---|
| 5a native | 1968 pass | 1968 |
| 5a AD (ForwardDiff vs 256-bit FD, prepared Mooncake) | 35 pass | 35 |
| 5b | 42 pass | 42 |
| 5c | 17 pass | 17 |
| 5d | 13 pass | 13 |
| 5e (gradients of the complete single pass) | 37 pass | 37 |
| 10a | 40 pass / 3 fail (exit 1) | 40 / 3 |

The 3 10a failures are the known ones: same sites, identical `Evaluated` values 1.655e-4, 5.206, 3.042e-3.
(The helper `run_focused.sh` logged `rc=0` for every file because `$(date)` was expanded before `$?`. Results were therefore read from the logs; the fixed helper is `A/chunk18/run_focused_v2.sh`.)

**Regression test** (normal suite, unconditional): `test/chunk5b_ode_pass_native.jl`, testset "reference callback is FullSpecialize (no FunctionWrappersWrapper
around the RHS)". It checks `ODEFUN5`, the problem built by `prob_phase5` (the one `solve_phase5` uses), and a live integrator's `integ.f.f === recombination_ode!`.
- Baseline construction (snapshot helpers with the same accessor building the old Auto problem): 42 pass / **4 fail**, including `integ.f.f` = `FunctionWrappersWrapper{…}`.
- New: **46/46 pass** (`A/chunk18/step1/regress/`).
- A first version also passed on baseline for its integrator assertions, because it built its own FullSpecialize problem. It was corrected to inspect `prob_phase5`; both logs are kept.

**Timing** (`A/chunk18/step1/bench/`; predeclared `BENCH_PROTOCOL.txt`; driver `guarded_bench_v3b.jl`, unchanged, sha256 recorded).
- Baseline = same `src` with the snapshot test helpers (Auto); after = the new helpers.
- FULL public `recombination_history_diffusion`: fiducial, 10000-node X_e/T_m, Rodas5P reltol 1e-12 / a1 1e-18 / aex 1e-14, 1 thread, 5 samples, evals = 1.
- Each process loads, compiles and warms the target and the BenchmarkTools thunk, then waits for an in-process QUIET gate, then times.
- Both processes pinned to CPU 8 (P-core). Launch 10:08, done 10:14; round 2 not needed.

Assessment by the unchanged rule (`match_r1.txt`, via `A/chunk18/match_v3.py`, which reproduces the main-verified rerun verdict exactly):

| | baseline (Auto) | FullSpecialize | criterion | met |
|---|---|---|---|---|
| gate | QUIET | QUIET | — | — |
| foreign > 50 % (ps) start / end | none / none | none / none | identical | yes |
| interval foreign busy cores over batch (/proc/stat, own process excluded) | 0.45 | 0.52 | reported | — |
| package T start → end | 59 → 98 °C | 59 → 94 °C | ≤ 5 °C | yes |
| pinned CPU 8 frequency start → end | 4800 → 4800 MHz | 4800 → 4700 MHz | ≤ 10 % | yes |
| throttle rate, package / core | 75.1 / 150.3 /s | 71.3 / 142.6 /s | ratio ≤ 1.5 | yes |
| all-CPU mean frequency start → end | 1482 → 1243 MHz | 1352 → 856 MHz | reported, not a criterion | would fail |

**Verdict: ACCEPTED — matched thermally throttled laptop measurement.** Both batches started cool and heated during the batch; this is not unthrottled peak, a CI or a certified algorithmic speedup.

| FULL public call | baseline (Auto) | FullSpecialize | after / baseline |
|---|---|---|---|
| samples [ms] | 2325.91, 2295.66, 2317.03, 2328.76, 2330.26 | 2163.89, 2170.89, 2137.08, 2136.76, 2195.63 | — |
| min [ms] | 2295.66 | 2136.76 | 0.931 (6.92 % lower) |
| median [ms] | 2325.91 | 2163.89 | 0.930 (6.97 % less elapsed; 1.075× throughput) |
| pairwise extremes | — | — | 0.917 … 0.956 (non-overlapping samples) |
| GC median [ms] | 196.6 | 193.4 | — |
| non-GC median [ms] | 2130.9 | 1977.5 | 0.928 |
| memory | 2503627696 B | 2503506160 B | 1.0000 |
| allocation events | 18195488 | 18197000 | 1.0001 |

- **Output check:** X_e/T_m of one warm call per process are byte-identical, and identical to the `F`-patch rerun output.
- **Cross-check:** this baseline median (2325.9 ms) is within 1.6 % of the main-verified `F`-patch after-run (2363.9 ms, same code, earlier session).
- **Cold:** process start to first result 45.5 s (baseline) vs 42.6 s (FullSpecialize). Single observations, reported separately; the extra FullSpecialize compilation is not visible at this resolution.
- **Interpretation:** the gain is time, not memory. The wrapper's dispatch overhead on 694672 RHS calls per prediction is removed, while allocations stay the same. This is consistent with the block-level probe (9–12 % on single blocks, `docs/SCIML_ODE_PERFORMANCE_REVIEW.md`); the whole prediction gains less because it includes the PDE stages, tail and assembly.

**Not yet run for Step 1:** the full root `Pkg.test`. It is planned once at the end of the accepted callback batch, per instruction; the focused suites above cover the touched routes.

**Status:** main review passed (main's 5b verification 46/46, rc 0, `A/chunk18/step1/main_5b_verification.log`); Step 2 started from the Step-1 working tree.

## Step 2 — per-solve prepared state Jacobian — MAIN-ACCEPTED (callback batch closed by the root `Pkg.test`)

**Baseline:** the Step-1 working tree, snapshotted to `A/chunk18/step2_baseline_snapshot/`. Gains are incremental relative to Step 1.

**Change** (`test/chunk5_helpers.jl` only):

```julia
mutable struct RHSAt5{P,Z}; p::P; z::Z; end                     # RHS at the LIVE (p, z) of the current call
(r::RHSAt5)(du, u) = (recombination_ode!(du, u, r.p, r.z); nothing)
struct PreparedJac5{R<:RHSAt5,C,Y}; f::R; cfg::C; y::Y; nfallback::Base.RefValue{Int}; end
PreparedJac5(u0, p, z0)   # JacobianConfig(f, similar(u0), u0, Chunk{length(u0)}()) from the solve's ACTUAL types
(j::PreparedJac5)(Jm, u, p, z)   # p, z, u types match the prepared ones: j.f.p = p; j.f.z = z; ForwardDiff.jacobian!(Jm, j.f, j.y, u, j.cfg)
                                 # otherwise: nfallback += 1; jac5!(Jm, u, p, z)   (unprepared oracle; never a conversion)
odefun5(u0, p, z0) = ODEFunction{true, SciMLBase.FullSpecialize}(recombination_ode!; jac = PreparedJac5(u0, p, z0), tgrad = tgrad5!)
prob_phase5(u0, z0, znodes, p) = ODEProblem{true, SciMLBase.FullSpecialize}(odefun5(u0, p, z0), u0, (z0, znodes[end]), p)
```

- One cache per problem/solve, so there is no global workspace and nothing is shared across solves, blocks or chains.
- The live `p` and `z` of every call are written in before differentiating; nothing physical is frozen.
- The chunk equals ForwardDiff's default for n ≤ 12 (12 or 7), so the arithmetic is unchanged.
- Unchanged: `jac5!`, `tgrad5!` (not rewritten in this step), `ODEFUN5` (kept as the unprepared reference), Rodas5P, tolerances, grids, physics, and the Mooncake reverse-route callbacks (`mkfun` in `chunk5e_helpers.jl`, which do not use `jac5!`).
- Not in scope: the B2 cached-LinearSolve true-reverse bug is untouched.

**Probe** (`A/chunk18/step2_probe.jl` → `step2/probe_probe.txt`, exit 0; reference = the Step-1 callback verbatim):

| check | result |
|---|---|
| prepared caches per public call | 108 (one per solve), **0** oracle fallbacks |
| public X_e/T_m | byte-identical to Step 1 |
| per-solve solver stats (nf, njacs, nw, naccept, nreject) | identical for all 108 solves |
| ForwardDiff 4-parameter 10b Jacobian (18×4, Dual solves) | byte-identical to the Step-1 file; 108 Dual-typed caches, **0** fallbacks |
| jac kernel 12-state (indicative, 1000 samples, not thermally matched) | 4.01 → 2.99 µs; 13200 → 5776 B; 30 → 21 allocations |
| jac kernel 7-state (same) | 2.04 → 1.41 µs; 8064 → 3888 B; 24 → 15 allocations |

The remaining kernel allocations are RHS temporaries evaluated on Duals (Step 3 target).

**Local Jacobian correctness: new unconditional regression tests** in `test/chunk5a_ode_rhs_ad.jl`.
- States: the real 5a states, 12-state z = 2500, 12-state off-table DP fallback z = 1500, 7-state z = 400.
- The prepared cache vs the fresh unprepared oracle:
  - primal J, bitwise;
  - **d J/d y, d J/d z, d J/d hscale** by ForwardDiff (the cache is built from Dual types, as a Dual solve does), bitwise;
  - BigFloat state, bitwise;
  - one cache with alternating live (y, p, z), and two independent caches interleaved, each equal to the oracle (no stale p/z);
  - a type-mismatched call (Dual z into a Float64-z cache) takes the counted oracle fallback, with derivative equal to the oracle's.
- `test/chunk5b_ode_pass_native.jl`: each problem gets its own cache; a real 12-state block (40 nodes) solves **bitwise** like the unprepared callback, with the same nf/njacs/naccept and 0 fallbacks.
- Regression proof: on the Step-1 helpers both files fail (`UndefVarError: PreparedJac5`; 5b 47 pass / 2 errored, 5a AD 35 pass / 6 errored, `step2/regress/`). With Step 2: 5b **49/49**, 5a AD **74/74**.
- The alternating-call tests show **independent cache state**. Actual concurrent (multi-threaded) chain execution was **not** tested and is not claimed.

**Prepared DI-Mooncake through the new cache** (main's coverage gate; `A/chunk18/step2_mooncake_probe.jl` → `step2/mooncake_probe.txt`, exit 0).
- Scalar projection s = Σ W ⊙ J (fixed random W), prepared `AutoMooncake` gradients (prepare once, two gradients), vs the ForwardDiff gradient through the oracle.
- Three paths: unprepared `jac5!`; a cache built inside the function; one **prebuilt** Float64 cache whose p/z mutation and reused config Mooncake differentiates through.

| case | d/d u | d/d z | d/d hscale |
|---|---|---|---|
| 7-state z = 400 (oracle / fresh / prebuilt) | 1.845e-16 (all three) | 6.400e-15 (all three) | 0 (all three) |
| 12-state z = 2500 (oracle / fresh / prebuilt) | 1.809e-14 (all three) | 3.314e-13 (all three) | 1.349e-15 (all three) |

- Values are relative errors vs ForwardDiff. Every path's repeated gradient is identical, with 0 cache fallbacks. The new cache adds no Mooncake difference: in each case its error equals the unprepared oracle's.
- Attempts 1–3 of this probe failed on **harness** bugs, kept as `mooncake_*_attempt{1,2,3}*`:
  1. a Float64 state was passed with a Dual z;
  2. `hscale`/`nbscale` were given different types;
  3. a parse error from an inline comment.

  None of these was a Mooncake or cache failure.
- This covers Mooncake through the Jacobian callback locally. It is not reverse mode through a full Rodas5P solve: the Mooncake solve routes use their own callbacks.

**Focused suites** (`A/chunk18/run_focused_v2.sh`, real exit codes; `step2/focused/`):

| test | result |
|---|---|
| 5a AD | 74 pass (35 + 39 new) |
| 5b | 49 pass (46 + 3 new) |
| 5a native | 1968 |
| 5c | 17 |
| 5d | 13 |
| 5e (complete single-pass gradients, incl. prepared Mooncake reverse route) | 37 |
| 10a | 40 pass / 3 fail, rc 1; the known failures with identical `Evaluated` values |

**Histories:** fiducial + ob, oc, h70, yhe, all files byte-identical to the baseline, so no CAMB rerun is needed.

**Timing** (`step2/BENCH_PROTOCOL.txt`, predeclared; base = Step-1 snapshot helpers, after = Step-2 helpers; same `src`; FULL public call, 5 samples; warm then in-process QUIET gate; pinned CPU 8).
Rule assessment (`step2/bench/match_r1.txt`): **ACCEPTED**.
- foreign none/none;
- package T 57/59 → 96/100 °C;
- CPU 8 4800/4800 → 4800/4800 MHz;
- throttle rates 73.7/77.6 and 147.4/155.2 /s;
- interval foreign busy cores 0.39/0.48;
- the all-CPU mean frequency would fail at both ends (1735/1488 → 1608/1367); it is reported, not a criterion.

| FULL public call | Step 1 | Step 2 | Step 2 / Step 1 |
|---|---|---|---|
| samples [ms] | 2137.09, 2119.18, 2136.56, 2116.24, **2486.26** | 2030.01, 2034.02, 2046.16, **2360.09**, 2097.49 | — |
| min [ms] | 2116.24 | 2030.01 | 0.959 (4.07 % lower) |
| median [ms] | 2136.56 | 2046.16 | 0.958 (4.23 % less elapsed; 1.044× throughput) |
| GC median [ms] | 192.6 | 156.2 | — |
| memory | 2503506160 B | 1920093792 B | **0.767** |
| allocation events | 18197000 | 17427122 | 0.958 |

- Matched thermally throttled laptop measurement.
- **Margin caveat:** each side has one outlier sample (bold), so the ranges overlap; the pairwise extremes are 0.817 … 1.115. Min and median are robust to the single outliers, but this pair is less decisive than Step 1's. The predeclared protocol allows a round 2 only for a rejected pair, so none was run.
- Outputs of both processes are byte-identical, and identical to Step 1's.

**Status:** main review positive (main noted: GC outliers overlap; a single pair, not a precise CI). Remaining: main's focused 5b/5a run, then the new normal-suite test `test/chunk5f_prepared_jac_mooncake.jl` (prepared DI-Mooncake through the cache: 7- and 12-state; oracle d/du, cache d/du, cache d/dz with live mutation; repeat gradient; 1e-12 vs ForwardDiff), then the root `Pkg.test` once for the callback batch. Step 3 (core RHS workspace) is **design only** (`docs/ODE_CORE_WORKSPACE_DESIGN.md`) after this review; implementation waits for the design gate and Marco's choice if a public contract changes.

## Callback batch (Steps 1–2): normal-suite CI and root `Pkg.test`

- **5f** (`test/chunk5f_prepared_jac_mooncake.jl`, new, in `runtests.jl`): **16/16**, rc 0 (focused, Mooncake 0.5.61, `A/chunk18/step2/ci5f/`); 16/16 in the root suite (Mooncake **0.5.62**).
- **Main's independent focused runs:** 5b 49/49 and 5a AD 74/74, exit 0.
- **Root `Pkg.test`, run once** (`COSMOREC_NATIVE_DATA_DIR=… julia --project=. -e 'using Pkg; Pkg.test()'`; 11:55–13:01; `A/chunk18/fullsuite_callback_batch/pkgtest.log`): exit **1**, **50876 passed, 6 failed**, 0 errored.
  - 50876 = 50814 (previous) + 39 (5a AD) + 7 (5b) + 16 (5f).
  - The 6 failures are the known ones at the same sites (`chunk10a_runmode0_native.jl:83,119,120`, `chunk10b_runmode0_ad.jl:42,56`, `chunk10d_background_link_ad.jl:46`) with **identical** `Evaluated` values. No new failure; no check relaxed.

## Step 3 — core RHS workspace

Design: `docs/ODE_CORE_WORKSPACE_DESIGN.md` (v2; main approved increment 1 = **private, primal-only** workspace, conditional on the batch above).
Increment-1 baseline = the Step-2 working tree, snapshot `A/chunk19/primal_baseline_snapshot/` (timing reference: FULL median 2046.2 ms, 1.920 GB, 17.43 M allocations).

### Increment 1 — private primal RHS workspace — MAIN-ACCEPTED; root `Pkg.test` 51056 passed / the same 6 known failures

**Allocation inventory first** (`A/chunk19/alloc_inventory.jl` → `alloc_inventory_baseline.txt`; `Profile.Allocs` sample rate 1, warmed single calls,
actual source path). The per-call objects are exactly:
- `X` (`ode_unpack`), `g` (`recombination_rhs!`);
- `dXH`, `dXHe` (`fcn_effective!`; `dXHe` even with helium off);
- H rates `A`, `B`, `R` (`get_rates`);
- He rates `A`, `B`, `R` (`get_helium_rates`).

This applies to the 12-state, the 7-state and the 7-state + feedback calls, Float64 and outer-Dual. Finding: in the outer-Dual case at z = 2706 the H-rate `A`
is **Float64**. In the detailed-balance branch (`db`, Tg/T0 − 1 > 2000) `_A = exp(log_qnl_qe(Tg))` depends on Tg only, so `A`'s element type is
branch-dependent. The workspace therefore keeps **two A buffers** (Tg-typed for `db`, promote(Tg, Te)-typed otherwise), selected by the same `db`.

**Source change** (private, additive; no export, public struct, dependency, StaticArrays or global):
- **Shared private cores, so the physics is not duplicated:**
  - `_fcn_effective_core!` (the former body of `fcn_effective!`, same statements and order, `dp_fallback::F … where {F}` kept);
  - `_ode_unpack!` and `_pack_rhs!` (the former bodies' assignments);
  - `_fill_hrates!` (the `get_rates!` loop) and `_fill_herates!` (the `get_helium_rates` loop).

  The public `fcn_effective!`, `ode_unpack`, `recombination_rhs!`, `get_rates!` and `get_helium_rates` call these with the **same allocations as
  before**, and remain the oracle. The comprehension-based `get_rates` is untouched.
- **New `src/RHSWorkspace.jl`:**
  - `_RHSWorkspace{TM,TY,TZ,TH,TN,TF,TX,TG,TA,TAdb,TB,TR,THe}` and `_rhs_workspace(y0, z0, rm; hscale, nbscale)`.
  - `_recombination_rhs_ws!(f, z, y, rm, ws; …) -> Bool`: one compatibility check per call (model type, argument types, feedback element type,
    table dimensions vs buffer sizes). A mismatch runs the public path and returns `false`; it never converts.
  - **Contract:** a workspace holds element types and buffer sizes only, no model values.
  - Buffer element types are derived from the actual dependencies: X from (y, fHe); g as in `recombination_rhs!` incl. the feedback type; H A/A_db/B/R and He
    from the same return expressions **at a zero interpolant** (no table value interpolated at a fictitious stencil; main QA item).
- **Test callback** (`test/chunk5_helpers.jl`): `RHSWS5`, a per-solve primal workspace functor with a counted fallback, as the ODE `f`. `PreparedJac5`
  (Jacobian) and `tgrad5!` still use the allocating `recombination_ode!` (Jacobian workspace = increment 2, tgrad later).

**Evidence** (`A/chunk19/inc1/`, re-run after main's QA edits in `inc1_postQA/`):
- **Probe** (`step3_probe.jl`), unchanged after QA:
  - public X_e/T_m byte-identical to the Step-2 callback;
  - all 108 per-solve stats identical;
  - **0** fallbacks (RHS workspace + Jacobian cache) over the public call and under the Dual solves of the 18×4 10b Jacobian;
  - the 10b Jacobian byte-identical to Step 2;
  - full-call allocations (one sample): **1920099088 B / 17427135 → 984247552 B / 4492087**.
- **Per-call allocations** (`@allocated`, warmed): workspace path **0 B**, vs the public path 1440 B (12-state) and 1040 B (7-state), unchanged.
- **New `test/chunk5g_rhs_workspace.jl`** (normal suite, unconditional; **179/179** post-QA). The private path vs the public allocating oracle:
  - bitwise at all 132 5a native rows (all patterns incl. the off-table DP fallback) and with the diffusion-feedback model;
  - bitwise for outer Dual state + Dual hscale, Dual hscale/nbscale, Dual z, Dual state, and BigFloat (12-state above and below the `db` switch, 7-state);
  - rate-buffer element types equal to `get_rates`/`get_helium_rates` output types in the branch taken;
  - type-mismatched call → public path (`false`);
  - ForwardDiff d/dy, d/dz, d/dhscale through the private path equal to the public ones;
  - prepared DI-Mooncake with **one workspace reused by two RHS evaluations in one taped function** (the first buffer overwritten by the second), and with a **fresh workspace constructed inside the objective** (main QA item), both within 1e-12 of ForwardDiff through the allocating path, with identical repeated gradients;
  - two independent workspaces alternating;
  - **model contract** (main QA item): a same-type model with a different value (`spin_forbidden` flipped) uses the workspace and matches the public path for that model; a different-type model (with diffusion) falls back;
  - per-call allocations below the public path.
- **5b** (`RHSWS5` is the solver's RHS; a 12-state block bitwise vs the allocating callback; 0 fallbacks, fresh workspace per problem): **50/50**.
- **Focused suites** (real exit codes):

  | test | result | baseline |
  |---|---|---|
  | 5f | 16 | 16 |
  | 5a native | 1968 | 1968 |
  | 5a AD | 74 | 74 |
  | refactored units: 2 native / 2 AD | 3760 / 1309 | identical |
  | 3b native / invariants / AD | 1783 / 95 / 320 | identical |
  | 3d native / AD | 566 / 2203 | identical |
  | 5c / 5d / 5e | 17 / 13 / 37 | 17 / 13 / 37 |
  | 10a | 40 / 3 fail, rc 1 | the same 3 failures, identical values |

- **Histories:** fiducial + ob, oc, h70, yhe, all files byte-identical to the baseline, so no CAMB rerun is needed.

**Timing** (`inc1/BENCH_PROTOCOL.txt`, predeclared). Baseline = the Step-2 tree as a **separate package copy** (`A/chunk19/primal_baseline_snapshot`, own env with the CosmoRec
path changed, own depot first; src and test both differ); after = the working tree. FULL public call, 5 samples, warm then in-process QUIET gate,
pinned CPU 8. Rule assessment (`inc1_postQA/bench/match_r1.txt`): **ACCEPTED**.
- foreign none/none;
- package T 56/57 → 98/98 °C;
- CPU 8 4800/4774 → 4800/4800 MHz;
- throttle rates 76.0/83.0 and 151.9/166.1 /s;
- interval foreign busy cores 0.41/0.50;
- the all-CPU mean frequency would fail (reported, not a criterion).

| FULL public call | Step 2 | increment 1 | ratio |
|---|---|---|---|
| samples [ms] | 2060.84, 2073.25, **2416.92**, 2111.94, 2085.84 | 1886.84, 1861.65, 1869.05, **2214.30**, 1901.76 | — |
| min [ms] | 2060.84 | 1861.65 | 0.903 (9.67 % lower) |
| median [ms] | 2085.84 | 1886.84 | 0.905 (9.54 % less elapsed; 1.106× throughput) |
| GC median [ms] | 173.0 | 79.5 | — |
| memory | 1920093792 B | 983754096 B | **0.512** |
| allocation events | 17427122 | 4491534 | **0.258** |

- **Main-verified:** a single hardware-conditioned matched pair (thermally throttled laptop, pinned P-core) gives ≈ 9.5 % less elapsed. This is **not** a global
  certification, CI or peak-performance claim; the memory (−48.8 %) and allocation-event (−74.2 %) reductions are measured independently of thermals.
- **Caveat:** one GC-heavy outlier sample per side (bold), so the ranges overlap; the pairwise extremes are 0.770 … 1.075.
- Outputs of both processes are byte-identical, and identical to Step 2's.
- The baseline median (2085.8 ms) is within 1.9 % of the earlier Step-2 after-run (2046.2 ms).
- **Cold costs not comparable:** the baseline package copy precompiled into its fresh private depot (897 s to first result vs 53 s).

**Root `Pkg.test` for increment 1.**
- **First attempt** (`A/chunk19/fullsuite_inc1/pkgtest.log`, started 14:42): terminated by **SIGTERM at 14:54**, with no summary and no exit file. The cause is unknown; the session's processes were lost and there was no host reboot. The log is preserved.
- **Retry** (`A/chunk19/fullsuite_inc1_retry_20261003T174133/`): a detached runner (own session, PID 2949850), the same source hashes as the interrupted attempt (checked), start hashes recorded, 17:41–18:45:34.
  - Exit **1**: **51056 passed, 6 failed**, 0 errored, 0 broken (63m52.8s).
  - 51056 = 50876 + 179 (5g) + 1 (5b).
  - The 6 failures are the known ones, at the same sites with **identical** expressions and `Evaluated` values vs the callback-batch run (main also compared all 6).
  - 5g 179/179, 5b 50/50, 5f 16/16. Mooncake 0.5.62 (root environment).

**Status: increment 1 main-accepted** (main verified: 0 B warm per call for 12- and 7-state; all 5 histories and the paired benchmark outputs byte-identical; ×0.905 median as a single hardware-conditioned pair with the GC-outlier caveat). Increment 2 (the Jacobian workspace with the family tag) needs its own gate, with the
CosmoRec-level tag checks. tgrad later.

### Increment 2 — buffered state Jacobian (family tag) — MAIN-ACCEPTED (main's 5h 70/70 + 5b 50/50, exit 0; final suite 51126 passed / the same 6 known failures)

Design and proof gate: `docs/ODE_CORE_WORKSPACE_DESIGN.md` §8–§9. Main passed the proof gate after:
- componentwise legacy parity of all Mooncake vectors;
- BigFloat 128/256 new == legacy bitwise;
- the off-table Float64 Mooncake-vs-ForwardDiff difference shown to be **pre-existing** and precision-sensitive in both engines (retained, not claimed fixed).

Baseline = the increment-1 tree (snapshot `A/chunk20/inc2_baseline_snapshot/`); `src` is unchanged by this increment.

**Change** (`test/chunk5_helpers.jl` only):
- `BufRHS5`: functor writing through a Dual-typed private `_rhs_workspace`, with the LIVE `p`/`z` written per call and a counted fallback.
- `bufrhs5_tag`: family tag `Tag{BufRHS5,V}`; passes the default check; `ForwardDiff.tagcount` is called explicitly (internal-API dependency, isolated here).
- `WSJac5` via the factory `prepare_wsjac5`: per solve, fixed chunk, default `jacobian!` check. On a type mismatch it calls the counted oracle `jac5!`, never a conversion.
- `odefun5` uses it. `PreparedJac5` and `jac5!` stay as oracles; `tgrad5!` is unchanged.

**Evidence** (`A/chunk20/inc2/run_20261003T204839/`, detached runner, chain exit 0):
- **Probe** (`chunk20/inc2_probe.jl`, reference = the increment-1 callback verbatim):
  - 108 caches, **0** fallbacks;
  - public X_e/T_m byte-identical; all per-solve stats identical;
  - 10b 18×4 Jacobian (Dual solves) byte-identical to increment 1, 0 fallbacks;
  - full-call allocations (one sample): **983759392 B / 4491547 → 501099568 B / 2801905**;
  - Jacobian kernel 0 B (was 5776 / 3888 B), indicative.
- **New `test/chunk5h_jacobian_workspace.jl`** (normal suite, unconditional): 52/52 in the chain, then **70/70** after main's coverage addition (`chunk20/inc2/5h_livepz_*`, rc 0; source hashes recorded; helpers unchanged). It covers:
  - values, ForwardDiff dJ/du, dJ/dz, dJ/dhscale and d²J/dz² bitwise vs `jac5!` at 12-state z = 2500, off-table z = 1500 and 7-state z = 400;
  - two independent caches alternating live (y, p, z);
  - a counted type-mismatch fallback (Dual z into a Float64 cache → oracle, equal derivative, count 1);
  - tag ordering, the default `checktag`, and the tag type `Tag{BufRHS5,Float64}`;
  - prepared Mooncake, one tape with two Jacobian evaluations through one prebuilt cache: in-table within 1e-12 of ForwardDiff (unchanged gate); off-table z = 1500 **bitwise parity with the legacy oracle's Mooncake vector** (labelled no-new-regression, not accuracy);
  - **added (main's coverage gap):** prepared Mooncake **d/dz and d/dhscale** of a J projection through one prebuilt cache at in-table 12/7-state, re-using the preparation at changed z/hscale (live mutation, no stale parameters), within 1e-12 of ForwardDiff through the oracle, repeated gradients identical, 0 fallbacks;
  - BigFloat(128) directional second derivative (k = 10, z = 1500), new == legacy;
  - allocations per Jacobian call: **0 B** vs 5776 B (12-state) / 3888 B (7-state).
- **5b** (`WSJac5` is the solver's Jacobian; block bitwise vs the allocating callback; 0 fallbacks): 50/50.
- **Focused suites** (real exit codes): 5g 179, 5f 16, 5a native 1968, 5a AD 74, 5c 17, 5d 13, 5e 37; **10a** 40 / 3 fail, rc 1, failures identical.
- **Histories:** fiducial + 4 cosmologies, all files byte-identical, so no CAMB rerun is needed.

**Timing** (`chunk20/inc2/BENCH_PROTOCOL.txt`, predeclared). Base = the same `src` with increment-1 helpers; after = increment-2 helpers. FULL public call,
5 samples, warm then in-process QUIET gate, pinned CPU 8. Rule (`bench/match_r1.txt`): **ACCEPTED**.
- foreign none/none;
- package T 52/52 → 93/96 °C;
- CPU 8 4800/4800 → 4800/4800 MHz;
- throttle rates 36.7/46.8 and 73.3/93.5 /s;
- interval foreign busy cores 0.30/0.29;
- the all-CPU mean frequency would fail (reported, not a criterion).

| FULL public call | increment 1 | increment 2 | ratio |
|---|---|---|---|
| samples [ms] | 1814.13, 1805.19, 1813.47, **2113.66**, 1807.78 | 1710.57, 1713.41, 1710.42, 1707.91, 1713.72 | — |
| min [ms] | 1805.19 | 1707.91 | 0.946 (5.39 % lower) |
| median [ms] | 1813.47 | 1710.57 | 0.943 (5.67 % less elapsed; 1.060× throughput) |
| GC median [ms] | 72.9 | 33.7 | — |
| memory | 983754096 B | 500719296 B | 0.509 |
| allocation events | 4491534 | 2801460 | 0.624 |

- Single hardware-conditioned matched pair (throttled laptop, pinned P-core); not unthrottled peak, CI or certification.
- One GC outlier in the base (2113.7 ms). The ranges do **not** overlap (pairwise extremes 0.808 … 0.949).
- Outputs of both processes are byte-identical, and identical to increment 1's.

**Not run for increment 2:** the root `Pkg.test`. Main releases it once after main's independent 5h/native runs.

**Status: increment 2 main-accepted; batch closed by the final root `Pkg.test`** (see Batch summary). Time-gradient workspace: not started.

## Batch summary (F patch, Steps 1–2, Step 3 increments 1–2) — all UNCOMMITTED

Each row is **incremental vs the previous accepted state**: no double credit. Every row is ONE main-verified, thermally logged, matched pair on a throttled
laptop (pinned P-core CPU 8; warm then in-process cool gate; unchanged matching rule). It is **not** an unthrottled peak, CI or universal certification. The
memory and allocation-event figures are measured independently of thermals. FULL = public `recombination_history_diffusion`, fiducial, 10000-node X_e/T_m,
caller Rodas5P callback (reltol 1e-12, a1 1e-18, aex 1e-14).

| step | change | FULL median (ms) | ratio | memory | allocation events | caveat |
|---|---|---|---|---|---|---|
| F patch | force specialization of the pass-through DP fallback (`src`, 6 signatures) | 3085.5 → 2363.9 | 0.766 | 3.695 → 2.332 GiB (×0.631) | 56.31 M → 18.20 M (×0.323) | samples non-overlapping |
| Step 1 | explicit `FullSpecialize` callback (test helpers) | 2325.9 → 2163.9 | 0.930 | unchanged | unchanged | samples non-overlapping |
| Step 2 | per-solve prepared state Jacobian `PreparedJac5` (test helpers) | 2136.6 → 2046.2 | 0.958 | ×0.767 | ×0.958 | one outlier per side; ranges overlap |
| Step 3 inc 1 | private primal RHS workspace (`src/RHSWorkspace.jl`, shared private cores; callback `RHSWS5`) | 2085.8 → 1886.8 | 0.905 | 1.920 → 0.984 GB (×0.512) | 17.43 M → 4.49 M (×0.258) | one GC outlier per side; ranges overlap |
| Step 3 inc 2 | buffered state Jacobian `WSJac5` (family tag; test helpers) | 1813.5 → 1710.6 | 0.943 | 0.984 → 0.501 GB (×0.509) | 4.49 M → 2.80 M (×0.624) | one base GC outlier; ranges non-overlapping |

- Successive baselines were remeasured in separate sessions, so the rows are **not** a single end-to-end ratio; they should not be multiplied as if they were.
- Per call (warm): RHS 0 B (was 1440 B 12-state / 1040 B 7-state); Jacobian 0 B (was 13200/8064 B before Step 2, 5776/3888 B after).
- Cold costs (first call) are reported per step only where comparable. Package-copy baselines with fresh private depots are not comparable.

**Numerics at every step:**
- public X_e/T_m at the fiducial and 4 nearby cosmologies byte-identical, so no CAMB rerun was needed;
- per-solve solver stats identical;
- 10b 4-parameter ForwardDiff Jacobian byte-identical;
- the known `ob` T_m 1.358e-7 miss and the 6 known suite failures unchanged.

**AD scope and limits:**
- ForwardDiff values and nested derivatives through every new cache are bitwise equal to the allocating oracles.
- Prepared DI-Mooncake through the RHS workspace and the Jacobian caches agrees with ForwardDiff within 1e-12 at the in-table states tested.
- At the off-table DP-fallback absorber state (12-state z = 1500), Float64 **second-order** derivatives (derivatives of J) are precision-sensitive in **both** engines. The disagreement is pre-existing, with the same vectors bitwise in the legacy oracle; BigFloat references differ from Float64 ForwardDiff by 235 % on one component, and Mooncake has the wrong sign there. The gate there is legacy **parity** (no new regression). Float64 off-table Hessian accuracy is **not** claimed, nor corrected. The culprit operation is not identified.
- True reverse mode through Rodas5P (B2) and the time-gradient workspace are **out of scope**: `tgrad5!` is still the allocating version.

**Process notes kept on record** (all negative results are preserved in `A/chunk16`–`A/chunk20`):
- the helper `rc` bug (`run_focused.sh` v1);
- harness bugs in several probes: a Float64/Dual buffer mix, constructor overwrites, a misplaced `end`, and an invalid parse check;
- one SIGTERM-interrupted `Pkg.test`, re-run;
- one regression test that initially passed on the baseline, corrected.

**Final root `Pkg.test` (closing the batch)** — `A/chunk20/fullsuite_final_20261003T221434/`. It was gated on main's independent 5h 70/70 + 5b 50/50 (exit 0) and run detached (runner PID 2988849, Julia PID 2990095), 22:23:04–23:29:13.
- **Verdict:** exit **1**, **51126 passed, 6 failed, 0 errored, 0 broken** (65m56.5s).
- 51126 = 51056 (increment-1 suite) + 70 (5h).
- The 6 failures are the known ones, at the same sites with **identical** expressions and `Evaluated` values (checked here and by main).
- Main verified that the sha256 of all 84 recorded `src`/`test` files was unchanged through the suite.
- Root Mooncake 0.5.62.

**State at the end of the batch:**
- whole public prediction ≈ **1.711 s** median (last matched pair, hardware-conditioned);
- **501 MB / 2.80 M** allocation events per call;
- primal RHS and state Jacobian **0 B** per warm call.

Unchanged: the 6 known suite failures, the `ob` T_m miss, and the Float64 off-table higher-derivative limit. This is **not** a fully green project. All code is **UNCOMMITTED**. The batch is complete; there is no further implementation (time-gradient workspace, B2/B3 and new optimizations are not started).
