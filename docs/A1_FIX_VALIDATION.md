# a1 fix validation (2026-10-02): runmode-0 test/reference callback configuration a1 = 1e-16 → 1e-18

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al. No native source, library or data was modified. No physics formulation, dependency or acceptance gate was changed. No Git operation was performed. Paths are relative to `A = cmbcheb_test/local_analysis/cosmorec_differentiability_20260930` unless they start with `src/`, `test/` or `docs/`.

**Status: the approved a1 validation task is COMPLETE, with a PARTIAL fix; not all issues are fixed.** **Verdict: the a1 fix resolves 10a:103 (both stages) and nothing else. Labelled diagnostics show the same tightening in the separate reverse route shrinks the hybrid mismatches (2.1× reduced, 1.4× full) without meeting their gates. The suite is NOT green: 50812 passed, 6 failed, 0 errored, 0 broken; rc 1. Phases 7–10 remain NOT accepted.**

## 1. Change (approved by Marco; minimal, a1 only)

**Scope (interface):** CosmoRec.jl has no built-in solver configuration. `recombination_pass(rm, solve_ode; …)` (`src/RecombinationODE.jl:148`), `recombination_history(rm, θ, solve_pass, solve_tail, zgrid; …)` (`:319`) and `recombination_history_diffusion(rm, θ, d, solve_pass, solve_tail, zgrid; …)` (`src/RecombinationDiffusion.jl:62`) all take a REQUIRED caller-supplied solve callback, and `Project.toml` has no ODE-solver dependency or default. This change edits only the injected Rodas5P callbacks of the test suite (the runmode-0 test/reference callback configuration used for these native comparisons). No library source or API default changed. **Downstream callers still choose their own solver and tolerances**, and nothing here guarantees their accuracy or cost.

The runmode-0 test/reference callback configuration used for these native comparisons is (reltol, a1, aex) = (1e-12, 1e-16, 1e-14) → (1e-12, **1e-18**, 1e-14). reltol, aex, the algorithm (Rodas5P) and `internalnorm` are unchanged. Every test callback carrying that triple was changed; nothing else was touched.

| file:line | constant | used by |
|---|---|---|
| `test/chunk10a_runmode0_native.jl:28` | `SOLVE10` | 10a runmode-0 public map, injected passes |
| `test/chunk5e_helpers.jl:17` | `TOL5E` → `SOLVEP5E` | 5e forward route, 10b composed ForwardDiff map, `PASS0_5E`/`K5E`/`SEL5E` |
| `test/chunk5b_ode_pass_native.jl:12` | `SOLVE5B` | 5b single pass |
| `test/chunk5c_recfast_tail_native.jl:12` | `SOLVE5C` | 5c pass + tail |
| `test/chunk5d_output_native.jl:11` | `SOLVE5D` | 5d single-pass history |
| `test/chunk6a_pde_coefficients_native.jl:80` | inline | 6a end-to-end coefficients from a Julia pass |

**Deliberately NOT changed**, because they carry different existing tolerances:
- the 5b tighter reference pass (1e-13/1e-17/1e-16);
- the reverse routes (5e `:rev` and 10c/`chunk10_reverse_helpers.jl`: reltol 1e-10, a1 1e-14, aex 1e-12);
- the Recfast tail solver;
- the `abstol5` default.

The benchmark scripts (`benchmark/chunk5b_benchmarks.jl`, `chunk5cd_benchmarks.jl`, `chunk7_10_benchmarks.jl`) still benchmark a1 = 1e-16. Their logged results belong to that configuration, and the matched comparison is in §4. Updating them is a decision for Marco.

The exact configuration before and after is recorded in `A/chunk12/config_before.txt` and `config_after.txt` (sha256 of every file, the lines, and the unified diff); pre-edit copies are in `A/chunk12/before_copies/`.

**Why a1** (source-derived, `docs/ORIGINAL_VS_PORT_ACCURACY_INVESTIGATION.md` §5.6): the Rodas weight of X1s is a1 + reltol·|X1s|, and X1s is the neutral fraction (7.4e-10 at z = 3000). With a1 = 1e-16 that weight is absolute-dominated for z > 1813, an effective relative tolerance of up to 1.3e-7. The PDE stage is hypersensitive to X1s. The matched-flag ablation found a1, not aex, decisive for 10a:103.

## 2. Every previously failing gate (old = `chunk10/pkgtest7_10_run2.log`, new = `A/chunk12/pkgtest_a1_run1.log`; both full root Pkg.test, `--check-bounds=yes` as Pkg.test sets it)

Line numbers in 10a moved by +1 because of the added comment (old 82/103/117/118/119 = new 83/104/118/119/120). Every new value below is identical in the focused run (`A/chunk12/focused_a1_*.log`) and in the full suite.

| gate (new line) | quantity, normalizer | bound | old value | new value | old → new | what changed it |
|---|---|---|---|---|---|---|
| 10a:83 | native-injected populations → Julia PDE1, 2g3s, max over all z of \|Δ\| / Julia residual scale | ≤ 1e-4 | 1.6550050940187906e-4 | 1.6550050940187906e-4 (bitwise) | FAIL → FAIL | nothing: native-injected input, no ODE solve |
| 10a:104, PDE0 | full path, max\|Δ\|/max\|native\| in 500<z<2000, each of DI1/2g3s/2g3d/R2s | ≤ 1e-3 | 1.0750e-3 (2g3s) | max of 4 = 3.27e-5 (R2s); 2g3s 1.23e-5 | **FAIL → PASS** | a1 |
| 10a:104, PDE1 | same, stage 2 | ≤ 1e-3 | 1.0393e-3 (2g3s) | max of 4 = 5.12e-5 (R2s); 2g3s 1.87e-5 | **FAIL → PASS** | a1 |
| 10a:118 | final Xe, `maximum(abs.(Xe .- Xn) ./ Xn)` over ALL 10000 nodes of the tested output grid `ZREC10` (= `GRID5[:,1]`), no z cut | ≤ 1e-5 (Marco-approved; was 5e-7) | 1.7057e-6 | 1.4494e-6 | FAIL (at 5e-7) → **PASS** | the PASS comes from the approved gate; a1 lowered the value by 15% |
| 10a:119 | pass-1 X1s, max \|Δ\|/(1e-7\|native\| + 1e-12) | < 1 | 3.246 | **5.206** | FAIL → FAIL (worse) | a1 moves it within the native's own pass-1 spread |
| 10a:120 | pass-1 2s, 2p, max \|Δ\|/(1e-6\|native\|) | < 10 | 15.38, 18.32 | 15.95, 16.18 | FAIL → FAIL | same |
| 10b:42 | composed ForwardDiff vs Float64 central FD on the frozen-decision map, best elasticity error over h = 1e-2, 1e-3, 1e-4 | < 1e-4 | 2.372e-3 | 3.042e-3 | FAIL → FAIL | the FD reference is noise-limited (test design) |
| 10b:56 | ForwardDiff directional derivative vs Float64 central FD, best over h, max \|Δ\|/\|f0\| | < 1e-5 | 3.506e-5 | 2.720e-5 | FAIL → FAIL | same |
| 10d:46 | prepared Mooncake VJP vs ForwardDiff, PDE stage on fixed native rows, worst of 4 points, block-max-weighted | < 1e-6 | 1.00093549296115e-6 | 1.00093549296115e-6 (bitwise) | FAIL → FAIL | nothing: fixed native rows, no ODE solve |

**Hybrid-gradient gates.** These are probe-level, not registered tests. Every route is outer prepared Mooncake + INNER ForwardDiffSensitivity in every ODE solve, i.e. NOT true reverse through Rodas5P. The route uses its own reverse-route tolerances (reltol 1e-10, a1 1e-14, aex 1e-12), which the test-callback edit does NOT touch. Its map depends on the edit only through `K5E` = 1342 and `SEL5E`, both unchanged (printed in the diagnostics). So the prior results apply unchanged after the callback edit.

A *separately labelled* diagnostic changed the reverse-route a1 to 1e-18, with everything else equal (default flags, `chunk2/test_env`, as the prior probes; no source change since the prior run except the TOL5E line):

| probe | metric, bound | prior (reverse a1 1e-14) | diagnostic (reverse a1 1e-18, labelled) | result |
|---|---|---|---|---|
| REDUCED composition (2 ODE, 1 PDE, no tail) | max\|g − J'w\|/max\|J'w\|, < 1e-6 | 1.469e-5 (`chunk10/probe_mc_composed_reduced_fds_bgfix.log`) | 6.945e-6 (`A/chunk12/probe_reduced_hybrid_reverse_a1_1e-18.log`, vectors `_vectors.txt`; rc 3, 1329 s, 15.3 GiB) | FAIL → FAIL (2.1× smaller) |
| FULL composed objective (3 ODE, 2 PDE, tail; checkpointed chain) | same, < 1e-6 | 1.9e-6 (`chunk10/probe_mc_composed_full_chain.log`) | 1.3599712114113554e-6 (`A/chunk12/probe_full_hybrid_reverse_a1_1e-18.log`, vectors `_vectors.txt`; rc 3, 3618 s, 16.0 GiB) | FAIL → FAIL (1.4× smaller) |

The labelled reverse-route a1 reduces both hybrid mismatches but meets neither gate. The remaining 1.4e-6 is the same size as the measured Float64 derivative conditioning: Float64 ForwardDiff on the frozen composed map is itself 5e-6 from the 128-bit derivative. That is consistent but not proven, and neither the gate nor the reverse-route tolerance was changed in any test.

## 3. Affected chunks that passed before: still passing (focused, `--check-bounds=yes`; full suite identical)

| chunk | old | new | notable observables, old → new |
|---|---|---|---|
| 5b single pass | 42/42 | 42/42 | X1s units 0.1139 → 0.1139; 2p units 2.386 → 2.597; Xe rel 1.730e-7 → 1.714e-7 |
| 5c tail | 17/17 | 17/17 | e2e Xe 1.384e-7 → 1.500e-7; e2e Te 3.05e-8 → 3.09e-8 |
| 5d history | 13/13 | 13/13 | inside Xe 1.643e-7 → 1.627e-7; Te 2.99e-8 → 3.02e-8 |
| 5e AD | 37/37 | 37/37 | FD_best 9.53e-6 → 1.15e-5; JVP_dir 3.45e-8 → 2.71e-8; VJP:all 1.14e-9 → 1.21e-9; k_switch shift 0 |
| 6a coefficients | 825/825 | 825/825 | e2e Dnem_scaled 2.25e-7 → 2.20e-7; e2e Dnem_raw 0.186 → 0.0042 (reported only, near-zero reference) |

## 4. Matched benchmarks (BenchmarkTools, default flags, same process, compilation excluded, evals = 1; `A/chunk12/a1_fix_benchmarks.jl`, `bench_a1.log`, `bench_a1_table.txt`)

| public route | samples | a1 1e-16 median | a1 1e-18 median | ratio | memory, old → new |
|---|---|---|---|---|---|
| `recombination_pass` (primal) | 10 | 0.747 s | 1.012 s | 1.35 | 904 → 1188 MiB |
| `recombination_history` 5d, 10000 nodes (primal) | 10 | 0.757 s | 1.022 s | 1.35 | 906 → 1191 MiB |
| `recombination_history_diffusion` 10a (primal) | 5 | 2.481 s | 3.244 s | 1.31 | 2926 → 3784 MiB |
| `ForwardDiff.jacobian`, 5e pipeline (18 × 4) | 3 | 1.591 s | 1.959 s | 1.23 | 2290 → 3001 MiB |
| `ForwardDiff.jacobian`, 10b runmode-0 (18 × 4) | 3 | 5.407 s | 6.487 s | 1.20 | 7681 → 9823 MiB |

- Prepared hot gradients: none of the prepared Mooncake routes uses the edited configuration (5e `:rev`, 10c and the hybrid probes use the reverse-route tolerances; 10d uses fixed rows). They were therefore not re-benchmarked. 10d's value is bitwise unchanged, and the hybrid map is unchanged as argued in §2.
- These timings are for the INJECTED test callbacks: Rodas5P, reltol 1e-12, aex 1e-14, with a1 1e-16 vs 1e-18. For example, 3.24 s is the 10a runmode-0 primal with the a1 1e-18 callback. They are not an API default or a performance guarantee; a caller's own solver choice determines its cost. The earlier benchmark logs (old a1) were not overwritten.
- No regression beyond the expected cost of tighter error control: about 1.2–1.35× time and memory.
- Saved: the derivative vectors (`bench_a1_jacobian{5e,10b}_{OLD_a1_1e-16,NEW_a1_1e-18}.txt`) and their comparison (`jacobian_old_vs_new.log`).
  - Physical columns (hscale, nbscale) change by 1.7e-10 / 1.5e-11 (5e) and 2.8e-6 / 1.1e-8 (10b).
  - The F / A2s1s columns are ~1e-18, effectively zero because F cancels in the rescaled Recfast tail; their relative change is a near-zero-reference artifact.

## 5. Full suite (`A/chunk12/pkgtest_a1_run1.log`)

- **Command:** `COSMOREC_NATIVE_DATA_DIR=<native Rec_database> julia --project=. -e 'using Pkg; Pkg.test()'` from the repository root. It is unconditional: no skips, no `@test_broken`, no environment gates. This is the same command as `pkgtest7_10_run2`.
- **Result:** **50812 passed, 6 failed, 0 errored, 0 broken; exit status 1**. START 12:37:32, END 13:33:23 (Test Summary 54 m 56.8 s; wall 55 m 51 s); peak RSS 16.4 GiB.
- **Failures:** 10a:83, 10a:119, 10a:120, 10b:42, 10b:56, 10d:46 (values in §2).
- **Before:** 50809 passed, 9 failed. The 3 resolved failures are 10a:104 ×2 (a1) and 10a:118 (the approved Xe gate).
- **Log hygiene:** in both logs, stderr `Info` messages from the chunk-1 probes interleave with some `Chunk10a` observable lines (e.g. `pde1:R2s:feedback_range:rel_to_max`, printed later at line 1307). The values are identical to the focused run.

## 6. NOT fixed by this change (stated separately)

1. **Native-baseline precision of pass 1 (10a:119/120).** The original itself does not meet these gates. Under 1e-15 H noise the native compared with itself gives X1s 3.4–5.9e-7 and 2p 1.0–3.2e-5 pointwise (investigation §5.9). No Julia tolerance setting meets them: all 24 ablation runs failed. The a1 fix makes X1s worse (3.25 → 5.21 units), which is consistent with that spread. This is a gate-policy decision for Marco; no gate was changed.
2. **10a:83:** a stage-level comparison at identical native input, at the native's own Float64 reproducibility in that metric (native vs itself under 1e-15 noise: up to 1.6e-4). Unaffected by any solver setting.
3. **10b:42/56:** a Float64 FD reference limited by primal noise. A **128-bit** BigFloat FD (`setprecision(BigFloat, 128)` in `chunk10/diag10b_big.jl` and `diag10b_big_nbscale.jl`) agrees with Float64 ForwardDiff to 5.6e-6 (direction) and 4.9e-6 (nbscale) on the frozen map (`chunk10/diag10b_big.log`, `diag10b_big_nbscale.log`). The 256-bit FD checks of the local 7–10 kernel AD tests are a different test. This is test design; no gate was changed.
4. **10d:46:** agreement between two Float64 AD modes on the hypersensitive PDE stage (fixed rows), bitwise unchanged.
5. **Hybrid reverse gates** (outer Mooncake + inner ForwardDiffSensitivity): both still fail, even with the labelled reverse-route a1 1e-18. Reduced: 1.47e-5 → 6.95e-6. Full: 1.9e-6 → 1.36e-6.
6. **True reverse mode through Rodas5P (B2)** stays BLOCKED by the cached `solve!(cache)` rule of LinearSolve 5.18.2 / Mooncake 0.5.61 (independently reproduced MRE in `chunk10/mre_linearsolve`). A solver tolerance does not touch that defect.
7. **Whole cosmological initialization (B3)** is still unimplemented, and the derivatives cover only four parameters. The preliminary Recfast history comes from the native fixture, and the Saha initial state, the loaded Hubble table and the closure constants are frozen. The H(z) policy is Marco's decision (`docs/REMAINING_BLOCKERS_DESIGN.md` B3).
8. **Remaining consistency item:** the standalone benchmark scripts (`benchmark/chunk5b_benchmarks.jl`, `chunk5cd_benchmarks.jl`, `chunk7_10_benchmarks.jl`) still configure a1 = 1e-16. They were deliberately not changed. The matched old/new comparison of the affected routes is already done in §4. Aligning these scripts is Marco's decision.

## 7. Artifacts (`A/chunk12/`)

- **Configuration:** `config_before.txt`, `config_after.txt`, `before_copies/`.
- **Focused runs:** `run_focused.sh`, `focused_status.txt` (rc, wall time and peak RSS per file), `focused_a1_*.log` and `.time` (failed logs preserved).
- **Full suite:** `pkgtest_a1_run1.log`.
- **Benchmarks:** `a1_fix_benchmarks.jl`, `bench_a1.log`, `bench_a1.time`, `bench_a1_table.txt`, Jacobian vectors, `jacobian_old_vs_new.log`.
- **Hybrid diagnostics:** `probe_reduced_hybrid_reverse_a1_1e-18.{jl,log,time}`, `_vectors.txt`, `_primal.txt`, `_jacobian.txt`; `probe_full_hybrid_reverse_a1_1e-18.{jl,log,time}`, `_vectors.txt`.
- **Ablation evidence:** `A/chunk11/ablation/`.
