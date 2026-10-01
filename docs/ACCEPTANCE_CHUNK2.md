# Chunk 2 review: native HI effective-rate lookup

## Scope / result

Focused native-rate capture/loader/interpolation work is **source-verified; rate-stage focused tests are independently accepted**. This is only the CosmoRec `get_rates` lookup for H I 3-shell / effective 500-state tables. It is not a population RHS, a recombination history, helium, or a production-quality differentiability claim. `src/RateTable.jl` is the first physics translation. Do not start Chunk 3 until the project-wide test gate is completed in a safe thermal window.

## Source-backed native references

Original read-only sources:

- `/home/marcobonici/Desktop/work/CosmologicalEmulators/cmbcheb_test/tools/CosmoRec`, SHA `086769055f61ae0c244a53dd381ee65b624d0ac3`.
- `Rec_database/Effective_Rates.HI/get_effective_rates.HI.{h,cpp}`. Native loader uses `load_rates(path,500,Gas_of_Atoms(3,1,1.0,false,0,-2))`; query function is the actual exported `get_rates(Tg,Te,A,B,R)` at `cpp:254-347`.

A small read-only C++ harness compiled OUTSIDE git and linked the existing untouched `libCosmoRec.a`, `libRecfast++.a`, and GSL2.8. It captures:

- Two native text fixtures: a compact window over the rate table's four-point interpolation stencils and direct native `get_rates` outputs at47 labeled queries. Both independently recaptured outputs compare byte-for-byte; original library and table hashes, state ordering, parameters, compile/link flags, version and CosmoRec citation notice accompany the fixtures.
- A third native domain/edge-probe fixture. Original OOB behavior, including native SIGSEGV/undefined reads, is recorded; Julia maps invalid inputs to a typed exception rather than intentionally reproducing memory corruption.

Fixture paths: `test/fixtures/native_hrates_table_window.txt`, `native_hrates_get_rates.txt`, and `native_hrates_domain_probes.txt`. Combined fixture payload is under400KiB, plain text; no full 338MB database or original code copied. The original C++ repositories remain untouched. Attribution and required paper citations are in `docs/CHUNK2_RESULTS.md` and the fixture NOTICE.

The native fixture reproducer is under analysis `chunk2/native_capture/`. The exported original `load_rates/get_rates` symbols, their expected atom/table setup, source lines and separate direct-C++ results are there. A native capture failure could not be substituted with Python or Julia synthetic values.

## Julia implementation reviewed

`src/RateTable.jl` reproduces native log(Tγ), log(Tₑ/Tγ) local 4-point Lagrange interpolation: bicubic A in transformed log space; 1D B/R; exact resolved-level/index ordering; detailed-balance branch at `Tg/2.725-1>2000`; `eps_A_effective=1e-4` clamp. It uses explicit `AtomicRateTable` input (no implicit loading/global caches) and exposes allocating plus in-place get_rates paths. The native data arrays remain caller-supplied; production table loading against the complete original database is not a checked-in dependency yet.

Dual-aware branch decisions use primal values, verified at the detailed-balance tie. Clamp/stencil seams are genuinely piecewise; tests evaluate derivatives separately on each side, not as if the function were globally smooth. Native output agreement is independently gated to `rtol=1e-14`; observed worst across1,645 returned numbers is `1.59e-16` with exact zero/NaN-pattern checks.

## Independent focused verification

In persistent `analysis/.../chunk2/test_env`:

- Native values/domain/interpolation tests: **3760/3760**, reproduced independently.
- ForwardDiff and prepared Mooncake tests: **1309/1309**, reproduced independently. Worst Mooncake-vs-ForwardDiff Jᵀw relative errors `9.3e-13`; one preparation is reused on18 changed queries and separate preparations are checked.
- BenchmarkTools scalar benchmarks independently reproduced: allocating lookup ~520ns interior/~322ns detailed-balance; in-place ~467ns/~316ns with zero allocations. These are table-sized scalar measurements only.

Worker logs/data: `analysis/.../chunk2/native_rates.log`, `rate_ad.log`, `benchmarks.log`, `native_fixture/`. Supervisor logs: `analysis/.../chunk2/supervisor_native_rates.log`, `supervisor_rate_ad.log`, `supervisor_benchmarks.log`.

## Unresolved gradients

Mooncake gives NaN gradients below ~55.6K where the original native rate computation overflows an unused A branch and the reverse pass encounters zero cotangent times infinity. ForwardDiff remains finite. The cosmological default-history limiter/matching behavior below native Recfast handoff still requires a full-physics design; do not hide or smooth this native-function pathology.

**Full-suite status supersedes the earlier thermal pauses:** after the user explicitly authorized continuing regardless of throttling, the supervisor ran the complete root suite in `analysis/.../chunk2/final_full_pkgtest.log`: **5969/5969 passing, exit0, 13m36s** (all1a/1b/1c/1d plus both Chunk2 suites). Earlier attempts and SIGTERM/SIGKILL thermal logs are preserved; they were interrupted runs, not code failures. A later Julia/Pkg precompile again reached the laptop's thermal limit under unrelated CPU tasks, so don't mislabel the current result: the completed full suite is the valid evidence. No unrelated processes were killed.

## Stop point

Chunk2 is accepted for its direct C++ reference harness, compact native rate fixtures, independently verified native/ForwardDiff/Mooncake focused tests and the subsequent **5969/5969 full Pkg.test**. QNDF reverse is still limited, low-T Mooncake-overflow behavior is documented. Next is a separately fixture-gated ODE-kernel step; this acceptance does not validate a recombination history or authorize omission of any active default RHS channels.
