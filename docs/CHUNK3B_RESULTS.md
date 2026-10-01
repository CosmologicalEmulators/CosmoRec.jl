# Chunk 3b results: default neutral-helium effective-population base RHS terms

**NOTICE (original-use acknowledgment).** The reference values come from outputs of the original CosmoRec v3.0b (J. Chluba et al.; native SHA
`086769055f61ae0c244a53dd381ee65b624d0ac3`, clean and read-only). Use must be acknowledged; cite Chluba & Thomas 2010 (MNRAS 412, 748) and
Chluba & Sunyaev 2006 (A&A 446, 39); also considered: Chluba, Vasil & Dursi 2010; Switzer & Hirata 2008; Grin & Hirata 2010;
Ali-Haimoud & Hirata 2010; Rubino-Martin et al. 2010. Bugs: Jens@Chluba.de. The fixture is a small text subset of derived values.

## Scope (what is and is not claimed)
Accepted for review: ONLY the helium terms of `fcn_HeI_effective` (`Modules/ODEdef_CosmoRec.cpp:217-250`) that need neither the HeI diffusion PDE,
helium feedback, nor H-I absorption. This is not a helium ODE, not the coupled H+He RHS and not default correctness.
**Excluded, recorded as later dependencies:**
- H-I absorption of HeI photons (`ODEdef_CosmoRec.cpp:~284-322`, Ly-a and the intercombination line; `parameters.CR.HI_absorption==1 && z<=zcrit_HI(=3400)*f_t`). **Chunk3c MUST add its native reference and Julia tests before any coupled first-pass H+He RHS claim.** The Chunk3a hydrogen harness (`flag_He=0`) is not evidence for it.
- HeI diffusion correction (`Diffusion_correction_HeI_is_on`, off in the first pass, set at `CosmoRec.cpp:439`), HeI feedback (`HeISTfeedback=0` in production), population-to-radiation coupling, initialization/Saha, the He switch, Recfast handoff, `rescale_rates`/`f_t`, the ODE itself.
- `read_DP_Data` and the Pesc/DP tables (needed only by absorption).

## Production settings checked against source
Runmode-0 defaults (`CosmoRec.cpp:836-838`): `nShellsHeI=2`, `spin_forbidden=1`, `HeI_Feedback=0`, `HI_absorption=2` (CAMB maps it to 1), `Diffusion_flag=1`.
Constructor defaults (`global_variables.cpp:159-163`) are not overridden: `nS_effective_HeI=30`, `set_HeI_Quadrupole_lines=1`.
`setup_Helium_atom` (`Modules/HeI_routines.cpp:70-140`): `Gas_of_HeI_Atoms(2, 10, 10, 10, -2)`, path `Effective_Rate_Tables.HeI.res_2/`, `load_rates(path, 30, ...)`.
The quadrupole/intercombination additions (`nQ=nTS=10`) do not add transitions to the n = 2 loop: the harness lists exactly 3 channels.

## Line-level native -> Julia mapping (`src/HeRHS.jl`)
| native | Julia |
|---|---|
| `get_rates_HeI` (`get_effective_rates.HeI.cpp:236-305`): 1-D, `lx = locate_JC row` (no HI `lx-1` shift), 4-point Lagrange in log Tg, `A = exp(log(qnl/qe))` explicit detailed balance, `B,R = exp(sum a_k table)`, `R[m][i<=m]=0` | `get_helium_rates`, `helium_A`, `helium_stencil_start` |
| `ODE_HeI_effective::evaluate_effective_Rci_Ric_terms` (`ODE_effective.cpp:~459-470`), `Nc = NH (fHe - Xi(0))` | `helium_continuum!` |
| `evaluate_effective_Rij_terms` (`:~498-520`), `nuij = nu_ion(j) - nu_ion(i)` | `helium_interlevel!` |
| `ODE_effective::evaluate_2s_two_photon_decay(Xi(0), Xi(1), Sing.Level(2,0).Get_Dnu_1s2(), const_HeI_A2s_1s)` | `helium_two_photon!` |
| loop `k=1..nres-1`, `Get_Trans_Data(ik,1,0,0,0)`, `w = Get_gw(ik)/T.gwp`, skip `Get_Level_index(2,1,1,1)` iff `spin_forbidden==0`, skip `A21==0` | `helium_lyman!` (reuses `lyman_channel!` of Chunk3a with `w`) |
| call order of `fcn_HeI_effective` | `helium_base_rhs!` (two methods: fixed rates / explicit table) |

Level mapping (native `Xi` index; Julia `X` index = +1): 0 `(1,0,0,0)`, 1 `(2,0,0,0)`, 2 `(2,1,0,1)`, 3 `(2,0,1,1)`, 4 `(2,1,1,0)`, 5 `(2,1,1,1)`, 6 `(2,1,1,2)`; 7 levels total.
Resolved states (table order m -> native index): 0->1, 1->2, 2->3, 3->5. Levels 4 and 6 (2^3P0, 2^3P2) are never touched by these terms; 0 (1^1S0) is not a table state.
Channels: k=1 -> index 2 (2^1P, A21=1.7989e9, w=3); k=2 -> index 3 (2^3S, A21=1.2724e-4, w=3); k=3 -> index 5 (2^3P1 intercombination, A21=177.58, w=3).
Units: Tg K, NH cm^-3, Hz s^-1, X = N/NH, dX/dt per second (cosmic time), Ai cm^3, Bi s^-1, Rij s^-1. Accumulating `+=` exactly like native.
`HeliumRateTable` is explicit (no cache, no default load): `lgTg`, natural-log `B[i,m]`, `R[i,m,j]` plus `gw`, `nu_ion`, `mu_red`; `A` is never tabulated (native enforces detailed balance). `load_helium_rate_table(dir)` reads the four native files on request. The hydrogen `AtomicRateTable` is not reused.

## Native fixture (FIRST gate)
- Out-of-repo, only in `cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk3b/native_capture/`: `harness3b.cpp` (sha256 `43f339efebf6816ff64cd5f878e9f067af8677d1958f7fde40183627c6fed14c`), binary `1eebf50c...bc3`, `cases.txt` (`3c9fcae1...2e32df2`; 40 cases), `build_cmd.txt`, `hashes_native.txt`, `assemble_fixture.py`, `gen_cases.py`.
- Built against the original static `libCosmoRec.a` (`2ac075e8...b69b783`), `libRecfast++.a` (`441a3eaa...a5fc`), GSL 2.8 static (`libgsl.a d675b920...7227a`, `libgslcblas.a ea9213cb...c6c`), g++ 13.3.0, `-Wall -pedantic -O2 -fPIC`, with `-IDevelopment/Helium.v1.1` first (the Makefile's `HELIUM_DIR`, whose `HeI_Atom.h` matches the library; the first attempt with `Helium/` failed to compile and was corrected, not worked around).
- Calls the original `load_rates(.., Gas_of_HeI_Atoms&)`, `get_rates_HeI`, `ODE_HeI_effective::*`, `ODE_effective::*` from the library. Layout copy only of `res_state_Data_HeI` (read-only, to print stencil rows). The base RHS is those same native calls in `fcn_HeI_effective` order (`fcn_HeI_effective` itself is not called because it needs the file-local `parameters` global).
- Two fresh-process runs are byte-identical (`run1.txt == run2.txt`, sha256 `9575abdf...b23172`); header repeat identical. Native checkout `git status` is clean (0 entries); no output or source was written there.
- Package fixture `test/fixtures/native_herhs_components.txt`: 110684 bytes, text only, **sha256 `eaef3d847eebe5e0c6fdc12d5b59cb3e1420c68ce84a1b550113d69b1f05e56b`**; contains the NOTICE, hashes of the 4 table files, level map, units, 40 cases (physical z = 200-3300 states, non-equilibrium states, unresolved-level populations non-zero, `Xe=0`, all-neutral, almost-ionized, zero excited, the lowest allowed Tg, and nodes 1,2,3,100,250,400,495 each at -1e-9, exact, +1e-9, the highest allowed stencil, a mid interval), the exact 500-point `lgTg` grid and the stencil windows (232 rows). The hydrogen fixture is unchanged (sha256 `03aaffdf...e4b`).

## Tests (all unconditional; no skips)
| file | result |
|---|---|
| `test/chunk3b_rhs_native.jl` (native comparison, sum rules, spin_forbidden, domain errors) | 1783/1783 |
| same file, hand-derived invariants | 95/95 |
| `test/chunk3b_rhs_ad.jl` (ForwardDiff + Mooncake) | 320/320 |
| focused re-run of Chunk2 native (3760) and Chunk3a native (179+93) | unchanged, pass |
Both Chunk3b files are in `test/runtests.jl`. The full `Pkg.test` was not re-run in this chunk (27 min / 98C in Chunk3a); only the focused suites were.

Native comparison (gate 1): rates `A,B,R` from the sparse native table vs `RT` at every probe within 2e-13 relative; each component and the combined RHS agree with the accumulate-from-sentinel native vectors to 1e-14 (fixed native rates) and 3e-13 (table-driven); increments (after minus sentinel) agree within `4 eps*1.25e-3 + 2e-12 max|inc|` (observed worst 1.1e-9 of max(|inc|, 1e-10), a sentinel-subtraction floor, not a model error). Exact structural checks: resolved-only writes of Rci/Ric, 2^3P0/2^3P2 untouched, accumulation signs, antisymmetry, sum rules.
Invariants (gate 2), written from the equations independent of the helpers: Saha/Boltzmann equilibrium annihilates every term and the combined RHS; transfer terms sum to zero; total helium changes only through the continuum; hand-evaluated 2-photon and Lyman values; node values equal tabulated values; near-continuity across seams; first/last allowed stencils; `A` finite at both table ends.

AD (gates 3-4), point `phys_z2300` (Tg 6300 K) with populations scaled by `[1,1.1,0.9,1.2,0.8,1.05,0.95]` (off equilibrium):
256-bit central differences vs ForwardDiff Jacobian, worst per-column-scaled error: 2PH 2e-16, RCI 2e-16, RIJ 2e-15, ALL_RATES 2e-15, TABLE 2e-13, COSMO 3e-13, LY 2e-11. Directional derivatives agree within 2e-13. Mooncake prepared VJP vs ForwardDiff J'w (3 seeds x 4 points, preparation reused at perturbed values, two independent preparations that must be identical): worst 4e-15 (LY), others <= 6e-16. Analytic Jacobian entries are asserted for 2PH and RCI. Cosmological-scalar derivatives (Tcmb, omega_b, Omega_m, h through a test-defined background; the background module is not part of this chunk) are nonzero and consistent.
Seams (one-sided only, never across): ForwardDiff vs 256-bit difference at 1e-6 from nodes 100/250, the lowest and the highest allowed stencil; worst 3.3e-7 relative at the lowest-Tg end, Mooncake VJP 2e-16. At an exact node the primal selects the upper stencil (documented, tested).

## Documented numerical limitations (not tuned away)
1. **Native stencil read past the table.** `lx + 3 > N` exits, but `lx + 3 == N` (the last 3 intervals above `lgTg[497]`, Tg > ~9430 K) makes native read `BiVec[N]` out of range (undefined behaviour). Julia throws `RateTableDomainError(:stencil_exceeds_table)` for `j+3 > N` (1-based), i.e. it is stricter than native for those 2 points above the last safe stencil; the highest valid probe is `top_lx496_minus`. Out of table (`Tg<547.7 K` or `>9540 K`) native `exit(0)`s; Julia throws.
2. **Sobolev `p = (1-exp(-tau))/tau` is evaluated without `expm1`** (as native). For the 2^3S line `tau ~ 1e-7`, so relative roundoff `eps/tau ~ 2e-9` is present in native values and in all AD derivatives; and `tau = 0` (both populations 0, or `w X1s == Xnp`) gives NaN in native and Julia. The LY Jacobian test therefore uses `fdtol = 1e-10`, not 1e-12.
3. **Seam derivative accuracy.** `d log rate/d logTg = sum a_k' f_k` with `|f_k|` up to ~70 and `a_k' ~ 175` has Float64 roundoff of ~1e-12 of the rate; at the lowest-Tg end where the derivative is small this is up to 3.3e-7 relative. The seam tolerance is 1e-6. Points closer than 1e-6 (relative in Tg) to a node were not used for AD because `log(Tg)` rounding (1e-15) is then a 1e-9 relative error in the weights.
4. No temperature switch exists in the helium rate lookup (single table, Tg = Te); the only piecewise behaviour is the stencil start and the `spin_forbidden` channel skip.
5. Equilibrium populations overflow below ~2000 K (Saha with the ground state); those states are outside the model domain and not tested. The sparse test table has finite rows only in captured windows, so perturbation tests keep Tg inside one native interval (asserted).

## Benchmark (BenchmarkTools, explicit args, `evals=1`, 2000 samples; CPU ~67C at start)
`benchmark/chunk3b_benchmarks.jl`, log in `chunk3b/benchmark.log`: `get_helium_rates` median 278 ns (400 B, 6 allocs); `helium_base_rhs!` from fixed rates 123 ns median, 0 allocations; through the table 393 ns median (400 B, 6 allocs). No retune or campaign.

## Files
`src/HeRHS.jl` (sha256 `104841bf...6c2c`), `src/CosmoRec.jl` (include + exports), `test/chunk3b_helpers.jl`, `test/chunk3b_rhs_native.jl`, `test/chunk3b_rhs_ad.jl`, `test/fixtures/native_herhs_components.txt`, `test/runtests.jl`, `benchmark/chunk3b_benchmarks.jl`, this file. No `docs/ACCEPTANCE*.md`, accepted fixtures or rate tables were edited; no commits, staging, pushes or tags.

## Remaining before a coupled first-pass H+He RHS
Chunk3c (H-I absorption of HeI Ly-a and intercombination photons, with `read_DP_Data`/Pesc tables and the electron/H-1s coupling signs), then the HeI diffusion and feedback decisions, the Saha initialization and He switch (`Set_HeI_Levels_to_Saha`), the background and the ODE.
