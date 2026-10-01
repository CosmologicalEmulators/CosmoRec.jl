# Chunk 3d results: assembled first-pass pointwise default H/He RHS (`fcn_effective`)

**NOTICE (original-use acknowledgment).** The reference values are outputs of the original CosmoRec v3.0b (J. Chluba et al.; native SHA
`086769055f61ae0c244a53dd381ee65b624d0ac3`, clean and read-only). Use must be acknowledged; cite Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006
(A&A 446, 39) and Rubino-Martin, Chluba & Sunyaev 2008; also considered: Chluba, Vasil & Dursi 2010; Switzer & Hirata 2008; Grin & Hirata 2010;
Ali-Haimoud & Hirata 2010. Bugs: Jens@Chluba.de. The fixture is a small text subset of derived values; no native table or source is bundled.

## Scope (what is and is not claimed)
Claimed: ONE pointwise evaluation of the production-default first-pass `fcn_effective(z, Data_Level_I&)` derivative `g = dX/dt` (15 equations), assembled in Julia from the
accepted H (3a), HeI-base (3b) and H-I-absorption (3c) blocks plus the matter-temperature equation, against the ORIGINAL exported `fcn_effective` on 12 explicit physical states.

NOT claimed / not ported: ODE integration or any recombination history; initial conditions; HeI diffusion/PDE; HeI radiative feedback (`HeISTfeedback`, `HeI_Feedback`);
H Ly-alpha diffusion; exotic sources (DM annihilation/decay, magnetic fields); the Cosmology module (`Tg = TCMB(z)`, `NH(z)`, `H(z)`, `fHe` are explicit inputs, `EffectiveBackground`);
the native explicit `DPesc` integral (`DPTableDomainError` when the absorber query leaves the DP table: all states with z <= 1500 in this fixture); Recfast completion; CAMB likelihood/spectra.
`src/CosmoRec.jl` is **not** a finished CosmoRec.

## Native oracle gate (done before any Julia composition)
Harness and logs live OUTSIDE all git repositories: `cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk3d/native_capture/` (`harness3d.cpp`, `build.sh`, `hashes.txt`, `table_hashes.txt` with 84 table files).
- Chain: original `read_entries_from_parameter_file` (ini mirrors the CAMB batch runmode 0 defaults) -> `set_startup_data_CR` -> `set_array_dimensions_and_parameters` -> `allocate_memory` ->
  `Set_Hydrogen_Levels_to_Saha(z)`, `Set_HeI_Levels_to_Saha(z)`, explicit multiplicative scaling of the 1s / excited populations, `X[neq-1] = rho` -> ORIGINAL `fcn_effective` (three calls per state: absorber on, on again, absorber off).
- The hidden `Parameters_form_initfile` is reached by forward declaration + `extern`, with source-verified layout copies (`global_variables.cpp:121-143`, `191-205`) and a runtime layout check (mismatch aborts, exit 3). No native source/binary was modified; `CosmoRec()` is never called.
- Runtime-proved configuration (`CF` records): `nShells = 3`, `nS_effective = 500`, `nShellsHeI = 2`, `nS_effective_HeI = 30`, `HI_absorption` remapped to 1 with `_HI_abs_appr_flag = 1`, `spin_forbidden = 1`, `HeI_Feedback = HeISTfeedback = 0`,
  `Diffusion_correction_is_on = Diffusion_correction_HeI_is_on = 0`, `DM_annihilation = DM_decay = magnetic_fields = 0`, `f_t = f_b = 1`, `flag_He = 1`, `fHe = 0.08170746361996017`.
- Layout (`LY`): `neq = 15`, `nHIeq = 6`, `nHeIeq = 7`, `index_HI = 1`, `index_HeI = 7`: `X[0] = Xe` (ignored by `fcn_effective`; zeroed in `g` and only touched by the absorber), `X[1..6]` = H 1s,2s,2p,3s,3p,3d, `X[7..13]` = HeI 1^1S0,2^1S0,2^1P1,2^3S1,2^3P0,2^3P1,2^3P2, `X[14] = rho = Tm/Tg`. `g` is d/dt in 1/s (cosmic time) per H nucleus.
- `compute_fractions` (flag_He = 1): `Xp = 1 - X[iHI]`, `XHeII = fHe - X[iHeI]`, `XHeIII = 0`, `Xe = Xp + XHeII` (the stored `X[0]` is not used). Implemented as `effective_fractions` and checked against the native `FR` records.
- 12 states: z = 3400.5 / 3400 / 3399 (absorber threshold both sides), 3000, 2600, 2200, 1800, 1500, 1200, 800 (several states carry perturbed rho and perturbed 1s/excited populations: z3000_pert, z2200, z1800_pert, z1200). All inputs are inside the native H table domain and the HeI Bitot/rate domains.
- Repeatability: fresh-process repeat byte-identical (`run1.txt` = `run2.txt`); reversed state order identical per label (60/60 comparisons); `on_cached` (second call, same z and rho) equals `on` for all states (native static cache hygiene); no NaN/inf in the final run.
- **Preserved failure** (`preserved_invalid_run/`, with `NOTE.txt`): the first attempt scaled the proton fraction so that `X1s = 1 - fXp Xp_saha` became negative (-0.002) at two states; native returned +-inf in g[1], g[3], g[5] with the absorber on or off. These are invalid physical inputs, not an oracle; the final states use multiplicative scaling of the populations, which stays positive.
- Absorber exercised: exactly zero above `zcrit = 3400` (z = 3400.5), active at z <= 3400 on native slots 0, 1, 7, 9, 12 only; native contribution magnitude (max over slots) between 5e-21 (z = 800) and 5e-12 (z = 1500).
- **Native-oracle limitation:** for z <= 1500 here the DP query is outside the pre-tabulated DP table, so native used the unported explicit `DPesc` integral. Julia throws `DPTableDomainError` for those states with the absorber on (tested) and is compared with the native absorber-off vector there.
- Late addition: the harness also emits the full H `lgTg` grid (`HLG`, fixture tag `lgTgfull`), needed so the compact table windows keep the true grid positions (the earlier H window fixture has gaps). Native outputs were regenerated byte-identically except for that new record (diff checked).

Fixture: `test/fixtures/native_fcn_effective.txt` (plain text, ~400 KB, sha256 `9f21e5e451a9eee5c83fc9df460bdf2f311cf363e9a7b5b8fb6f7c347cfedc34`, verified by a test). Records: `CF` config, `LY` layout, `BK` background, `IN` inputs, `FR` fractions, `GO` outputs (modes on / on_cached / off),
H-table stencil windows (`lgrho`, `lgTgfull`, `lgTg`, `B`, `R`, `A`), HeI grid and windows (`TH`, `L`, `W`), DP table/Bitot windows (`DT`, `DI`, `DH`, `DS`, `DR`) and the 40 `f.corr.dat` knots (`F`). Tables are NaN-filled outside the captured windows, so any query outside them fails loudly.

## Julia composition (`src/FcnEffective.jl`)
`fcn_effective!(g, z, X, bg::EffectiveBackground, m::EffectiveModel)` / `fcn_effective(...)`. Explicit workspace: `EffectiveModel` holds the H table, HeI table, DP table, Bitot series, fcorr spline, atom/constants and the two switches (`spin_forbidden`, `hi_absorption`); nothing is cached or loaded globally.
Native order: zero `g`; `g[14]` assigned by the matter-temperature equation; H block accumulated into `g[1:6]`; HeI base block into `g[7:13]`; absorber last (native `ODEdef_CosmoRec.cpp:284-321`). Production switches (diffusion, feedback, exotic OFF) are documented in the file header. `z` is only used through its primal for the absorber branch (`z <= 3400`).

## Results
- `test/chunk3d_fcn_native.jl`: **566 / 566**. All 15 components x 12 states x {absorber on, absorber off} match native with `rtol = 1e-13`; **worst relative error 8.3e-15** (state z3000, mode on, component g[0], the absorber-only electron slot); typical 1e-16. Also: fixture sha256, runtime-config keys, native `on_cached == on`, ground-state electron fraction (`Xe = Xp + XHeII`, `X[0]` scrambling has no effect), absorber threshold both sides and at `prevfloat`/`nextfloat(3400)`, exact zero of the non-absorber slots, absorber transfer invariants (electron created = H 1s destroyed, HeI population conserved), population-conservation of the two-photon/Lyman/interlevel transfer terms (H and HeI), `d rho/dt = -H` at rho = 1, `DPTableDomainError` out of the DP table.
- `test/chunk3d_fcn_ad.jl`: **2203 / 2203** (about 85 s, including compilation). Parameter vector `p = (X[15], Tg, NH, Hz, fHe)`; components rho, hydrogen, helium, absorber, full (absorber off and on), per native state.
  - ForwardDiff run in 256-bit arithmetic equals the 256-bit central difference (step 1e-20 relative) to < 1e-13 after rounding to Float64 (AD formulas exact).
  - Float64 ForwardDiff vs the 256-bit central difference: worst column-scaled error **1.7e-8** (below_zcrit absorber). It is not an AD error: the Float64 *value* itself differs from the 256-bit value by ~3e-12 because of near-cancelling terms in the absorber, and the bound used is 1e-7. 
  - Directional derivative vs `J v`: worst 9.6e-14.
  - Two independent prepared Mooncake preparations (reused at varied points, 3 weight vectors each) vs ForwardDiff `J'w`: worst relative error **2.9e-15**, identical to each other.
  - Full composition on = off + absorber at the Jacobian level (9 states); above `zcrit` the absorber Jacobian is exactly zero; the `X[0]` column is exactly zero.
- **No differentiability across switches is claimed:** the absorber is a branch in `z` at `zcrit = 3400` (values jump from 0 to the full contribution; `z` is never differentiated), the table stencils are piecewise cubic (derivatives are one-sided at node crossings), and the DP-domain edge throws. The AD points are interior, with perturbations of 1e-5 relative.
- `benchmark/chunk3d_benchmarks.jl` (BenchmarkTools, `evals = 1`, laptop, thermally limited; state z2600): `fcn_effective!` 1.14 us median (16 allocs, 1.06 KiB; absorber off 0.81 us); `ForwardDiff.jacobian` (15 x 19) 7.5 us; Mooncake hot prepared gradient 148 us (first `prepare_gradient` ~0.8 ms median, 0.75 s first sample with compilation).
- Full `Pkg.test()` (single process, unconditional, log `chunk3d/pkgtest.log` in the external analysis directory): **20863 / 20863 passed**, exit 0, 16m42s test time (18094 prior + 566 native + 2203 AD). Run under the human's thermal-throttling override.

## Owner decisions / open items
- Out-of-DP-table absorber queries (z <= 1500 in this fixture, and generally wherever native falls back to the explicit `DPesc` integral) are not ported. The absorber is inactive only above z = 3400, so the DP-table coverage of an actual trajectory must be established (or the integral ported) before any ODE step.
- The Cosmology background (`TCMB`, `NH`, `H`, `fHe`) is still an explicit input; the xp/xe relation is the native `compute_fractions` and was not reinterpreted.
- Float64 AD of the triplet absorber term carries ~2e-8 relative conditioning error against the exact derivative (value-level rounding, not an AD defect).
