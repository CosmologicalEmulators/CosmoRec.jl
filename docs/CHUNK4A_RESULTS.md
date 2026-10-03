# Chunk 4a results: explicit-input Saha initialization and ODE-vector packing

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al. Use must be acknowledged; cite Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39), Rubino-Martin, Chluba & Sunyaev 2008. No native source or data was changed; no table or history file is bundled.

## Scope (and nothing more)
Implemented in `src/SahaInit.jl` (included from `src/CosmoRec.jl`, exports added):
- `SahaConstants`/`NATIVE_SAHA_CONSTANTS`, `SahaLevels`/`NATIVE_SAHA_HYDROGEN`/`NATIVE_SAHA_HELIUM` (explicit level data: statistical-weight numerator, `Eion` [erg], `mu_red`, `ME_scale`).
- `saha_lte_hydrogen`, `saha_lte_helium`: native `Ni_NeNc_LTE` (H `g/2/gc lambdac^3 (2 pi kb_mec2 T ME mu)^-3/2 exp(Eion/kB/T)`, He `gw/4 lambdac^3 (2 pi kb_mec2 mu ME T)^-3/2 exp(Eion/kB/T)`, same factor ordering as the C++).
- `SahaInputs(fHe, Xe_Seager, Xp_raw, NH, Te, Tg, XHeII, XHeI1s)`: every value the native routines read from the `Cosmos` object is an explicit argument (`XHeII = NHeII/NH`, `XHeI1s = NHeI/NH`, `Te = cosmos.Te(z)`, `Tg = TCMB(z)`); `z` is not an input.
- `saha_initial_state` = native `Set_Hydrogen_Levels_to_Saha` then `Set_HeI_Levels_to_Saha` into `Level_I.X` (layout `[Xe | H 1s,2s,2p,3s,3p,3d | HeI x7 | rho]`, length 15): H: `X[1]=min(Xe_Seager,1+fHe)`, `Xp=min(Xp_raw,1)`, `X_i = Xe Xp NH f_i(Te)`, ground override `1-Xp`, `rho=1`; He (does not touch `X[1]`): unclipped `Xe_Seager`, `X_i = Xe XHeII NH f_i(Tg)`, ground override `XHeI1s`.
- `pack_ysol`, `NATIVE_RESHI=(1,2,3,4,5)`, `NATIVE_RESHE=(1,2,3,5)`: native `copy_LI_to_ysol` for `flag_He >= 1` only: `[rho, H 1s, resolved H, He 1s, resolved He]` (length 12).
NOT implemented / not claimed: Recfast++ history or its `Xe_H`/`Xe_He` (these are IONIZED populations `n_e,H/n_H`, `n_e,He/n_H`, not neutral fractions X1s/XHeI1s), Cosmos accessors/splines and their z>=3500 branches, H(z), cosmology, the sampled-HeI state switch (`flag_He = 0` branch of the packing), any ODE/PDE. The 4a function takes accessor values as inputs and does not validate them.

## Native fixture provenance
`test/fixtures/native_saha_init.txt` (107 lines, SHA-256 `d61d10b689e572e6848bd5361503e16825d9ecb32f79a41164b3e8cd886e7656`, tested). Mechanically extracted (`grep`, no edits) from `.../chunk4/native_capture/out/probe1.txt` (SHA-256 `14450f2b90fc1f1b70d7e31c93fb078c80a77ecea782d580d7da053389c2ea2c`), produced by the external read-only harness `harness4.cpp` (SHA-256 `65699f725de4504de6b2cc46343ce7a05668a24e2fdce7b084b6767b1c758888`) calling the ORIGINAL `Set_Hydrogen_Levels_to_Saha`, `Set_HeI_Levels_to_Saha`, `copy_LI_to_ysol` on the CAMB-adapter batch route (thermo-fixture cosmology `H0=67.36, ombh2=0.02237, omch2=0.12, Mnu=0.06, omk=0, TCMB=2.7255, nnu=3.046, YHe=0.24568275240335613`, H(z) from CAMB as an input); two fresh-process outputs byte-identical; `CosmoRec()` never run. Source/library/harness hashes: `chunk4/native_capture/hashes.txt`. Records: 12 redshifts (3600, 3500, 3400.5, 3400, 3399, 3000, 2800, 2500, 2000, 1500, 1200, 800): inputs (`SAHAIN`), native LTE factors (`LTEH`, `LTEHE`), native `Level_I.X` (`XLI`), packed vector (`YSOL`), level data/constants, layout, resolved-level maps, `fHe`, and `SP` (raw `Xe_Seager` at 9 of the 12 z).
Caveat: `SAHAIN` prints the clipped `Xe`; the clip is inactive at all 12 states (tested; raw equals clipped bitwise where `SP` has it), so the unclipped value used by helium equals it.

## Results
Focused native suite `test/chunk4a_saha_native.jl`: **765/765**. Maximum relative errors against the native values: full `X` vector and packed `y` 5.06e-16, H LTE 4.16e-16, He LTE 4.41e-16 (tolerance 1e-14). Bitwise equalities (tested): Xe slot, `rho == 1`, H and He ground overrides, native-LTE products (`Xe Xp NH LTE`) for every level of every state, and `pack_ysol(native X) == native YSOL` exactly for all 12 states.
Focused AD suite `test/chunk4a_saha_ad.jl`: **280/280** (z = 3000 base point only, 2000, 1500, 1200, 800; state vector and packed y; 8 inputs):
- 256-bit ForwardDiff vs 256-bit central difference: 0.0. Float64 ForwardDiff vs 256-bit central difference: <= 1.16e-15 (tolerance 1e-12).
- Directional derivative vs `J v`: <= 2.0e-28 relative (tolerance 1e-10). Prepared Mooncake VJPs (two independent preparations, reused at changed inputs, 3 seeds each) vs ForwardDiff `J'w`: <= 2.2e-16 (tolerance 1e-9), `gA == gB`.
Full root `Pkg.test()` (`chunk4a/pkgtest.log`, 16:13:36 to 16:30:55 EDT): `CosmoRec.jl | 23701 23701 17m12.9s`, `Testing CosmoRec tests passed`, `rc=0`; 23701 = 22656 (previous) + 765 + 280. The supervisor independently reran both focused suites (765/765, 280/280; logs `supervisor_native.log`, `supervisor_ad.log`) with identical error maxima. Tests are unconditional (no environment gates); type inference of `saha_initial_state` and `pack_ysol` is checked with `@inferred`.
An earlier focused AD run (`chunk4a/ad_run1.log`, tolerance 1e-11 before it was tightened to 1e-12) passed 280/280 as well; no failed run occurred.

## Branches and limitations
- Branch/override behavior is reproduced, not smoothed: `Xe` clip at `1+fHe` (H only), `Xp` clip at 1 (H only), ground-state overrides (H 1s `=1-Xp`, He 1s `=XHeI1s`, whose derivatives w.r.t. the Saha inputs are zero/-1/1 by construction, tested analytically), `rho` fixed to 1 (zero derivative). In a clipped region the derivative with respect to the clipped input is zero; at the clip edge (`Xp_raw` near 1: at z = 3000 it is 7.4e-10 below the clip) no derivative is claimed, so the AD tests use interior points only (z = 3000 at its base point without perturbed points).
- The clip and override cases have no native numerical oracle (native states are all interior): they are tested with synthetic explicit inputs against the same formulas.
- No domain checks (as native): NaN/negative/zero temperatures propagate; `exp(Eion/(kB T))` can overflow.
- Level data/constants for the production 3-shell H and 2-shell HeI atom only; other shell counts require new data. `pack_ysol` does not implement the helium-off packing.
- Parametric coverage: one cosmology, 12 redshifts. The four unperturbed chunk3d `IN` rows were matched bit-for-bit by the C++ harness in standalone mode at design time, but are not re-tested in Julia (see `docs/CHUNK4_DESIGN.md` section 11).

## Benchmark (`benchmark/chunk4a_benchmarks.jl`, `chunk4a/bench.log`, BenchmarkTools, evals = 1, z = 1500 state)
| path | median | memory |
|---|---|---|
| `saha_initial_state` (15-vector) | 244 ns | 176 B, 2 allocs |
| `pack_ysol` (12-vector) | 21 ns | 160 B, 2 allocs |
| initialize + pack | 250 ns | 336 B, 4 allocs |
| ForwardDiff Jacobian (12 x 8) | 1.22 us | 4.09 KiB, 9 allocs |
| Mooncake `prepare_gradient` (3 samples, first includes compilation; median not a cold time) | 62.5 us | 28.8 KiB, 350 allocs |
| Mooncake prepared hot gradient | 2.49 us | 1.42 KiB, 25 allocs |
Fixture parsing is outside the timed region.

## Files
`src/SahaInit.jl`, `src/CosmoRec.jl` (include/export), `test/chunk4a_helpers.jl`, `test/chunk4a_saha_native.jl`, `test/chunk4a_saha_ad.jl`, `test/runtests.jl` (includes), `test/fixtures/native_saha_init.txt`, `benchmark/chunk4a_benchmarks.jl`, `docs/CHUNK4_DESIGN.md` (section 11 amendments), this file, `docs/ACCEPTANCE_CHUNK4A.md`. No Git operation was performed.
