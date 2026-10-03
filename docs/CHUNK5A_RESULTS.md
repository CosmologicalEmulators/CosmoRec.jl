# Chunk 5a results: the packed-state ODE right-hand side dy/dz (12- and 7-state)

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al.; cite Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39), Rubino-Martin, Chluba & Sunyaev 2008. No native source, library or table was changed or bundled.

## Data policy (decision D1, resolved by the user)
The Phase 5 tests REQUIRE the environment variable `COSMOREC_NATIVE_DATA_DIR` = the native `Rec_database` directory (`.../cmbcheb_test/tools/CosmoRec/Rec_database`; read-only, never copied). `test/chunk5_helpers.jl` raises an error (the suite fails loudly, never skips) if it is unset, not a directory, lacks `Effective_Rates.HI`, `Effective_Rates.HeI` or `Pesc_Data`, or if `f.corr.dat` (`<dir>/../Development/Recombination/Data.fcorr/f.corr.dat`) is missing. The loader `load_native_effective_model` reads: HI `Effective_Rate_Tables.nS_3/Rates_n*_l*.nS_500.dat`, HeI `Effective_Rate_Tables.HeI.res_2/HeI_Rates_*.nS_30.dat` (also the Bitot column), `Pesc_Data/DP_Coll_Data.31.fac_50.neff_30.dat`, `f.corr.dat`. Run the tests with `COSMOREC_NATIVE_DATA_DIR=... julia --project=. -e 'using Pkg; Pkg.test()'`.

## Scope
`src/RecombinationODE.jl` (included from `src/CosmoRec.jl`, exports added):
- `RecombinationModel(eff, cosmos[, dp_fallback])`: 3d `EffectiveModel` (full tables), 4b `CosmosAccessors` (background and loaded Hubble), native DP fallback (3e) for off-table absorber queries.
- State (native `ysol`): 12 entries `[rho, X(H 1s), X(H 2s,2p,3s,3p,3d), X(He 1^1S), X(He 2^1S, 2^1P, 2^3S, 2^3P1)]` for `flag_He = 1`; 7 entries (hydrogen only) for `flag_He = 0`. `ode_unpack` is the native `copy_ysol_to_LI` (X[1] = Xe recomputed from the ground states; the two UNRESOLVED triplet levels are placeholders, tested not to enter `g`).
- `recombination_rhs!(f, z, y, rm; flag_He)`: native `fcn_effective(int*, ...)` (col < 0): `f = g/dz_dt`, `dz_dt = -H(z)(1+z)`, `g` = the 3d `fcn_effective!`, now with a `flag_He` keyword (`compute_fractions` with `XHeII = 0`, `xe = xp`, no helium equations, no HI absorber of the helium lines) - the `flag_He = 0` branch is new and verified against the native.
- Side fix found by AD tests: `RateResult` is parametrized by three element types (the table value `A` depends on `Te` while `B`, `R` depend only on `Tg`).
NOT captured/claimed: the native hand-written Jacobian `jac_effective` (informational in the design; the port uses an exact AD Jacobian), the native solver.

## Native fixture
`test/fixtures/native_ode_rhs.txt` (SHA-256 `f9b9e03dd14acc91267d97ece7599fc8c49b223508e5e221557535737ae09906`): 132 direct calls of the ORIGINAL exported `fcn_effective(int*, double*, double*, double*, int)` at explicit states (60 `flag_He = 1` states: the original Saha initialization at z = 3000, 2800, 2500, 2200, 2000, 1800, 1500, 1200, 1000, 800 times 6 perturbation patterns; 72 `flag_He = 0` states at z = 1500 ... 50 times 6 patterns), with the native background (`Tg`, `NH`, `Hz`, `Xp`) printed per state. External harness `chunk5/native_capture/harness5a.cpp` (built with the Chunk 3/4 flags against the original libraries; native source untouched), two fresh-process outputs byte-identical (0.09 s).

## Results (focused: native 1966/1966, AD 35/35)
- Background from the 4b accessors vs the native background: TCMB, NH 1e-14; H 1e-13; Xp 1e-12.
- dy/dz, all 132 states, component-wise relative error (floor 1e-12 of the largest component): `flag_He = 1` 9.5e-15, `flag_He = 0` 2.5e-16 (gate 1e-12). Off-table DP fallback queries are exercised (z <= 1500).
- AD of the RHS in the state (12 x 12 and 7 x 7) at z = 2500, 1500, 1200 (12-state) and 400, 100 (7-state): 256-bit ForwardDiff vs 256-bit central differences 0.0 (<= 2e-51); Float64 ForwardDiff vs 256-bit 3.6e-16 worst; prepared Mooncake VJP (two independent preparations, changed states) vs ForwardDiff `J'w` 2.4e-16 worst (gates 1e-12).
- Limitations: table stencil seams, the DP table domain, the Patterson order and the absorber switch z = 3400 are piecewise (interior points only). The 7-state states are explicit states, not a trajectory.

## Benchmarks (`benchmark/chunk5a_benchmarks.jl`, `chunk5/bench5a.log`)
Table load + model construction 2.4 s (one cold measurement incl. compilation). RHS `flag_He = 1` (z = 1500, DP fallback) 51.8 us (8.4 KiB); `flag_He = 0` 0.84 us; ForwardDiff Jacobian 12 x 12 55.9 us; Mooncake `prepare_gradient` 482 ms (2 samples, includes compilation); prepared hot gradient 763 us.

## Files
`src/RecombinationODE.jl`, `src/FcnEffective.jl` (flag_He), `src/RateTable.jl` (RateResult), `src/CosmoRec.jl`, `test/chunk5_helpers.jl`, `test/chunk5a_ode_rhs_native.jl`, `test/chunk5a_ode_rhs_ad.jl`, `test/runtests.jl`, `test/fixtures/native_ode_rhs.txt`, `benchmark/chunk5a_benchmarks.jl`, this file, `docs/ACCEPTANCE_CHUNK5A.md`. No Git operation.
