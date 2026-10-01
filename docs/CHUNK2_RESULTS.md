# Chunk 2: native HI effective-rate lookup (`get_rates`) in pure Julia

Scope: **hydrogen effective-rate table lookup only** (default 3-shell H, `nS_effective = 500`, 5 resolved states, `neqres = 5`).
No population RHS, no ODE, no helium, no Chunk 3. Nothing is committed or staged.

## NOTICE and citations (CosmoRec licence)

The fixtures contain outputs and table values produced by **CosmoRec** (J. Chluba et al., v3.0b, 2017). Its README licence requires that
use is acknowledged in publications and that the following are cited: Chluba & Thomas 2010, MNRAS 412, 748; Chluba & Sunyaev 2006, A&A 446, 39.
Also considered by the authors: Chluba, Vasil & Dursi 2010, MNRAS 407, 599; Switzer & Hirata 2008, PRD 77, 083006; Grin & Hirata 2010,
PRD 81, 083005; Ali-Haimoud & Hirata 2010, PRD 82, 063521; Rubino-Martin et al. 2010, MNRAS 403, 439. Bug reports to the CosmoRec author.
"No guarantee for the correctness of the outputs is given" (upstream). This repo contains only a small subset of table values (72 of 500
rows) and the returned numbers of the unmodified native `get_rates()`; the native source and binaries are untouched.

## Native capture (gate A)

* Harness outside all git repos: `.../cosmorec_differentiability_20260930/chunk2/native_capture/` (`harness.cpp`, `build_cmd.txt`, `compile.log`).
  Links the existing `libCosmoRec.a`, `libRecfast++.a`, static GSL 2.8. Uses `Gas_of_Atoms(nS=3,Z=1,Np=1.0,Qlines_on=false,Rec_flag=0,mflag=-2)`,
  `nS_effective=500`, table `Rec_database/Effective_Rates.HI/Effective_Rate_Tables.nS_3/`. No native edit, no table copy, no `make`.
  The harness uses a layout copy of the native `res_state_Data` struct (native lines 42-66) to read the knots for the table dump.
* CosmoRec sha `086769055f61ae0c244a53dd381ee65b624d0ac3`; both native repos `git status` clean after the capture.
* Two fresh processes: byte-identical outputs; reversed query order: byte-identical per query.
* Fixtures (plain text, NOTICE/provenance headers incl. library/table hashes, flags, GSL version, constructor, state order, units statement, query contract):

| file | bytes | sha256 |
|---|---|---|
| `test/fixtures/native_hrates_table_window.txt` | 315734 | `907ba3f82db4f0328ea595c93d7ac524d7cea5491b135c9865c50fc725309209` |
| `test/fixtures/native_hrates_get_rates.txt` | 36535 | `9297b02fa093c1da324b17368fabde5483cbe9ca53dc0e027dadc93181f5fc4d` (records `b4d069e3...45c65`) |
| `test/fixtures/native_hrates_domain_probes.txt` | 1631 | `eef02bece59830617d53c1eefe84eaeb09cc76d0e6bb4c35c281d8ddd3fccbd5` |

The immutable chunk-1d fixture `native_camb_cosmorec_thermo.txt` (sha `394bf61f...6c96`) is unchanged. The table window holds 72 contiguous-run rows
(union of the 4-row stencils of all queries), so `locate` on the window reproduces the native stencils. 47 queries x 5 states = 235 records.
Units of A, B, R were not independently verified here.

## Native source trace (`Rec_database/Effective_Rates.HI/get_effective_rates.HI.cpp`)

| lines | meaning | Julia (`src/RateTable.jl`) |
|---|---|---|
| 71, 72, 300 | `eps_A_effective=1e-4`, `z_detailed_balance=2000`, `Tg/2.725-1.0 > z` -> `A = exp(log_qnl_qe(Tg))` | `_setup`, `_A` (`db`) |
| 233-241 | exact 4-point Lagrange weights | `lagrange4` |
| 243-249 | `log_qnl_qe_func` | `log_qnl_qe` |
| 254-267 | load checks, `logTg=log(Tg)`, `logrho=log(Te/Tg)`, out-of-table message + `exit(0)` | `_setup`, `stencil_start` |
| 269-292 | `locate_JC` (clamped), stencil `lx = j-1` if `j>=1` else `0`, weights | `locate`, `stencil_start` |
| 285 | out-of-bounds exit | `RateTableDomainError(:stencil_exceeds_table)` |
| 305-321 | A: bicubic of `A - B - log_qnl_qe(knot)`, clamp `|exp(fxy)-1|<=eps -> 0`, `exp(fxy+log_qnl_qe(Tg))` | `_A` |
| B, R (to ~335) | B = exp(1D Lagrange), `R[m][i<=m]=0` else exp(1D Lagrange) | `_B`, `_R` |

`get_rates_all` (line 350) is not substituted. API: `AtomicRateTable` (explicit, typed, no global cache or implicit file load), `get_rates`, `get_rates!`,
`load_rate_table(dir, 500, states)` and `read_native_rate_file(io)` (explicit loader, not used by tests against the real production directory).

## Semantics

* Coordinates `log(Tg)` and `log(Te/Tg)`; interior: bicubic for A, cubic for B and R; level order and R orientation as native (row = resolved state, column = lower-triangle zeros).
* Cell clamp: cells 0 and 1 share stencil start 0 (1-based: 1); a stencil needing a row/column beyond the table is an error.
* Native failure modes (documented, **not** patched, Julia throws `RateTableDomainError` instead of `exit`/UB): out of table -> message + `exit(0)`;
  Tg in [9430.68, 9540.22) -> SIGSEGV (bounds check admits an OOB read); Te/Tg in [1.075, 1.1) -> OOB read of the A row (garbage); Tg at the max -> `exit(0)`.
* Native NaN: Tg below about 55.6 K makes `exp(h_kb*nuion/Tg)` overflow for the n=2 states (A = NaN); the Julia port reproduces the NaN pattern.
* Branch decisions use **primal values**. Finding: `ForwardDiff.Dual` comparisons break value ties by the partials (`Dual(2000.0,1,0) > 2000.0 == true`), which put the AD
  branch on the DB side at the exact tie `Tg/2.725-1 == 2000` while the primal was on the table side. Fixed with `_primal` in `locate`, bounds, DB and clamp tests.
* Derivatives are piecewise: dA/dTe = 0 in the DB branch and in the clamp; at stencil seams and branch switches they jump (reported below, never compared across).

## Results (focused tests, persistent env `chunk2/test_env`, absolute `--project`)

| gate | file | tests | result |
|---|---|---|---|
| B/C | `test/chunk2_rate_table_native.jl` | 3760 | all pass |
| D | `test/chunk2_rate_table_ad.jl` | 1309 | all pass |

Native comparison (B): 1645 values (47 queries x 5 states x A, B, R0..R4), NaN pattern equal, zeros exact, `rtol = 1e-14` required; **worst observed relative difference 1.59e-16**; no tolerance tuning.
Independent properties (C): Newton divided-difference oracle on native table values (interior, rtol 1e-12/1e-9), closed-form detailed-balance expression on both sides of the switch
(smallest double with `Tg/2.725-1 > 2000` is 5452.7250000000013), eps clamp both sides (11 clamped, 24 unclamped, 19 unclamped within 3e-4 of the threshold), exact-knot Lagrange weights `(0,1,0,0)`,
exact polynomial reproduction on a synthetic table (66 random points), query-order independence, domain errors for all 10 native probes plus non-finite inputs, loader parse test.

AD (D):

| comparison | worst relative error |
|---|---|
| ForwardDiff Jacobian vs hand-derived local-polynomial derivatives (788 entries, interior + low/high boundary stencils + DB) | 5.7e-13 |
| Mooncake VJP vs ForwardDiff J'w, A | 1.9e-13 |
| ... B | 2.7e-13 |
| ... R | 2.0e-13 |
| ... A+B+R | 9.3e-13 |
| Mooncake, one prepared cache reused over 18 changed same-shape (query, weight) pairs | 1.7e-13 |
| Mooncake, independent preparations (same 18 pairs) | 1.7e-13 |

Seams / branches (reported): d/dTg jump across the Tg=100-knot stencil seam is 1.06 (relative, max over finite outputs; the stencil changes), 2.5e-6 at the Tg=250 knot seam;
each side matches its own analytic local polynomial. DB switch: dA/dTe below is nonzero, above is exactly 0. Clamp: ratio 0.99984 is smooth (dA/dTe = -2.9e-19 for the first state),
0.99986 is clamped (dA/dTe = 0).

**Known limitation (tested, not hidden):** Mooncake at the native low-Tg overflow point (Tg = 31 K) returns NaN for every gradient component, even when the objective uses only finite outputs,
because the reverse pass multiplies a zero cotangent by the infinite derivative of the overflowed (unused) A. ForwardDiff stays finite. Reverse mode is therefore only valid for Tg above about 55.6 K.

## Benchmarks (BenchmarkTools, scalar lookup only, after warm-up; `benchmark/chunk2_benchmarks.jl`)

2000 samples x 10 evals, one thread, table window fixture. Table branch: `get_rates` 520 ns median (464 B, 6 allocs), `get_rates!` 470 ns (0 alloc).
DB branch: `get_rates` 322 ns (464 B, 6 allocs), `get_rates!` 316 ns (0 alloc). The machine was shared with unrelated jobs; treat as indicative.

## Thermal status and process log

* Focused tests and the benchmark were run only when `sensors` Package id 0 was below 85 C (58-67 C).
* **Full suite subsequently independently passed:** after the human explicitly authorized ignoring throttling, the supervisor ran `julia --project=<repo> -e 'using Pkg; Pkg.test()'`, exit 0, **5969/5969**, 13m36s; authoritative log `local_analysis/cosmorec_differentiability_20260930/chunk2/final_full_pkgtest.log`.
* Earlier attempts, including a worker PTY attempt, hit CPU package thermal conditions and were stopped; their logs remain (including `supervisor_pkgtest_pty.log` and the first `supervisor_pkgtest.log`). The later completed full run supersedes those attempts. No CMB-lite processes were terminated.
* Failed attempts: harness dynamic GSL link (libgsl.so.28 not found; fixed by static link); first domain exit-code capture reported a wrong rc (rewritten, confirmed by direct runs 139/0);
  first test run had my own mistaken expectations (reversed one-argument `occursin`, wrong clamp-label assumption, synthetic table crossing the clamp, 35 vs 20 packed outputs) - fixed in the tests, not by loosening tolerances;
  a thermal gate bug in my shell helper (parsed `crit = +100.0` instead of the reading) which caused a spurious 7-minute wait; the `Dual` tie-break finding above.

## How to run

```
julia --project=<chunk2/test_env> -e 'using Test; include("test/chunk2_rate_table_native.jl")'
julia --project=<chunk2/test_env> -e 'using Test; include("test/chunk2_rate_table_ad.jl")'
julia --project=<chunk2/test_env> benchmark/chunk2_benchmarks.jl
```
Both test files are also included from `test/runtests.jl`. STOPPED after Chunk 2.
