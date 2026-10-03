# Chunk 4 design: initialization slice (design + native reference gate only)

Status: DESIGN AND NATIVE-ORACLE GATE ONLY. No Julia source, tests, fixtures, ODE, PDE or Git action were added in this turn. Everything below is a source trace of the
ORIGINAL CosmoRec v3.0b (SHA `086769055f61ae0c244a53dd381ee65b624d0ac3`, (c) J. Chluba et al.; cite Chluba & Thomas 2010, Chluba & Sunyaev 2006, Rubino-Martin et al. 2008)
and of the pinned CAMB-cosmorec checkout (`fa3f097343fbbe427cc04b4f5f0041c22c6ec764`), plus a new external native capture. Paths: `N4 = .../cosmorec_differentiability_20260930/chunk4/native_capture/`
(outside every Git repo), `CR = .../cmbcheb_test/tools/CosmoRec`.

## 0. What this design does NOT claim
No Recfast++ port, no Julia Cosmology/H(z), no ODE or SciML recombination solve, no sampled-HeI state switch, no diffusion/feedback, no full history/spectrum parity, no claim that the
final CAMB thermodynamics are reproduced. Nothing here is implemented. The oracle values below are native outputs for ONE cosmology (the thermo-fixture cosmology) and the default `runmode 0`.

## 1. Production path (what the oracle must mimic)
CAMB `fortran/cosmorec.f90` (`TCosmoRec_init`): `Nz=10000`, grid `z_i = 1e4 - (i-1) 1e4/(Nz-1)` (descending), `Hz(i) = 1/dtauda * (1+z)^2 / MPC_in_sec` [1/s] with `MPC_in_sec = Mpc/c`,
`Mpc = 3.085677581e22` (`fortran/constants.f90:44-46`), `runpars = -1` except `(fdm, accuracy=0, 0, A2s1s=0, B0=0, nB=-1, PMF=2, 1, 1, 0, 6, 0,0,0,0)`, then
`CosmoRec_calc_h_cpp(runmode=0, runpars, OmegaC, OmegaB, OmegaK, N_eff, H0, tcmb, Yhe, zrec, Hz, Nz, zrec, xrec, tmrec, Nz, label)` -> `cosmorec_calc_h_cpp_` (`CR/CosmoRec.cpp:1013`):
`cosmos.init_Hubble(z_Hz, Hz, nz)` FIRST, then the batch `CosmoRec(...)` (`CosmoRec.cpp:729-950`).
Batch values for `runmode 0, accuracy 0` (all verified at runtime by the harness): `nz=3000, zstart=3000, zend=0`; `Yp=YHe`, `T0=tcmb`, `Omega_m=omegac+omegab`, `Omega_L=0` input (closed flat inside `Cosmos::init`; observed 0.68613569155460208),
`Omega_k`, `h100=H0/100`, `Nnu=N_eff`; `Rec.F=1.14`, `Rec.f_ann=0`, `B0=0`; CR: `nShells=3, nS_effective=500, nShellsHeI=2, HI_absorption=2` (remapped to 1 with `_HI_abs_appr_flag=1` in `set_startup_data_CR`, `setup_functions.cpp:~405`),
`spin_forbidden=1, HeI_Feedback=0, Diffusion_flag=1, induced_flag=2, nS_2gamma=3, nS_Raman=2`, `Diff_iteration_max=2`. `set_startup_data_CR` builds the cosmology with `include_CF=0` (`setup_functions.cpp:~531`), sets up atoms.
Observed oddity (not interpreted): after startup `Rec.A2s1s` reads 8.2206 although `runpars[3]=0` (replaced by default in `check_and_set_params`, `setup_functions.cpp:303`).
Cosmology of the oracle = thermo-fixture cosmology: `H0=67.36, ombh2=0.02237, omch2=0.12, Mnu=0.06 (1 massive), omk=0, TCMB=2.7255, nnu=3.046, YHe=0.24568275240335613`; H(z) input table generated with CAMB 2.0.4 (`camb_H_input.py`; an INPUT, never a native output).

## 2. Stage A: preliminary Recfast++ history and the Cosmos accessors
Source map
- `Cosmos::init`/`set_initials` (`Cosmos.cpp:~213`, `:287-293`): `n_Xe=6000`, `zsRe=2.5e4`; `init_splines -> calc_coeff_X_spline(6000, 2.5e4, 0)` (`Cosmos.cpp:393-480`).
- It calls `recombination_history(M, zs, 0, zarr, Xe_H, Xe_He, Xe, dXe, dX_H, TM)` (`Cosmos.cpp:345-390`) = Recfast++ `Xe_frac` (`Recfast++.cpp:119`, RHS `evaluate_Recfast_System`, `evalode.Recfast.cpp:99`; own stiff solver `ODE_solver.Recfast.cpp`, 1011 lines; Recfast++ is ~2.8k lines) with the Hubble function `Hfunction -> cosmos.H`.
- Recfast's `Xe_H` is `n_e,H / n_H` (hydrogen ionized fraction), NOT neutral `X1s`; `Xe_He` is `n_e,He / n_H`, NOT neutral `XHeI1s`. The original spline inputs (GSL natural cubic, `routines.cpp:347-413`, evaluated in z with `x` nudged by 1e-14 at the ends, NO extrapolation: an out-of-range query goes to the GSL error handler) are `ln Xe`, `dXe` (linear), `ln Xe_H`, `dXe_H` (linear), `rho = TM/TCMB`, `ln max(Xe_He,1e-20)`. The low-z accessors derive `X1s = 1-exp(spline(ln Xe_H))` and `XHeI1s = fHe-exp(spline(ln Xe_He))`; both switch to direct Saha formulas for `z>=3500`.
- Accessors (`Cosmos.h:355-433`, `Cosmos.cpp:749-781`): `Xe_Seager(z)` = `exp(spline)` for `z<zsRe` else `(1-Yp/2(2-1/fac_mHemH))/(1-Yp)`; `X1s(z)` = `1-exp(spline)` for `z<3500` else `SahaBoltz_HI1s(TCMB)`; `Xp = 1-X1s`; `XHeI1s = fHe - exp(spline)` (<3500) else `SahaBoltz_HeI1s`; `XHeII1s = exp(spline)` (<3500) else `SahaBoltz_HeII1s`; `NHeI/NHeII = X*NH`; `Te_Tg` spline (or `1/(1+1/kappa_cool)` for z>=zsRe); `Te = TCMB*Te_Tg`; `NH=(1-Yp)Nb`; `TCMB=T0(1+z)`.
- Hubble (`Cosmos.cpp:567`, `Hubble::init`): inside `(1e-10, 1e4)` a natural cubic spline of `ln H` vs `ln z` through the CAMB table; outside the analytic `H0 sqrt(OL + (1+z)^2 (Ok + (1+z)(Om + Orel(1+z))))`.
  NATIVE HAZARD (reproduced): `Hubble::init` never assigns `lHz[0]` (the `z=0 -> 1e-10` node keeps `ln H = 0`), so `H(z)` is wrong at low z: `H(1e-3)` is 10.7x the table, `H(0.01)` +32 %, `H(0.5)` -29 %, `H(1)` -6e-5, `H(10)` -7e-7, `H(100)` -1.4e-8, `H(3000)` -1.5e-9 (relative to a log-log linear interpolation of the table; the last values mix in the cubic-vs-linear difference). Irrelevant at the initialization redshift, but any later Julia `H(z)` claim must state it.
Oracle (new): see section 6, files `probe1.txt` (records HZ, SP, COS, CFG) and `history1.txt` (all 6000 nodes x 7 columns).
Design position: Stage A reproduction requires the Recfast++ ODE and its solver step control; it is NOT part of the first slice. The first slice consumes the preliminary history/spline outputs as explicit inputs.

## 3. Stage B: per-level initialization at a supplied z (the bounded slice)
Source: `Set_Hydrogen_Levels_to_Saha` (`CR/Modules/HI_routines.cpp:150-175`), `Set_HeI_Levels_to_Saha` (`CR/Modules/HeI_routines.cpp:140-171`), `Gas_of_Atoms::Xi_Saha/Ni_NeNc_LTE` (`Atom.cpp:1346-1359`),
`Electron_Level_HeI_{Singlet,Triplet}::Xi_Saha/Ni_NeNc_LTE` (`HeI_Atom.cpp:960-972, :1332-1344`), packing `copy_LI_to_ysol` (`ODEdef_CosmoRec.cpp:87-107`), call order `main.CosmoRec.cpp:111-112`.
Equations (native names; X = N/N_H):
- H, executed first: `X[0] = min(Xe_Seager(z), 1+fHe)`; `Xe = X[0]`; `Xp = min(1 - X1s(z), 1)`; for `i=0..nHIeq-1`: `X[1+i] = Xe Xp NH f_i(Te)`, `Te = TCMB(z) Te_Tg(z)`,
  `f_i(T) = (g_i/2) lambda_c^3 (2 pi kb_mec2 T ME_i mu_i)^(-3/2) exp(Eion_i/(kB T))`, `g_i = 2(2 l_i + 1)`; then OVERRIDE `X[1] = 1 - Xp`; `X[neq-1] = 1` (rho = Te/Tg is set to exactly 1 regardless of the spline value).
- He, executed second (never touches `X[0]`): `Xe = Xe_Seager(z)` (UNclipped), `XHeII = NHeII(z)/NH(z)`, `TM = TCMB(z)` (not Te); `X[iHe+i] = Xe XHeII NH f^He_i(TM)`,
  `f^He_i(T) = (gw_i/4) lambda_c^3 (2 pi kb_mec2 mu_i ME_i T)^(-3/2) exp(Eion_i/(kB T))`; then OVERRIDE `X[iHe] = NHeI(z)/NH(z)`. Same formula for singlet, j-resolved triplet and non-j triplet levels (here: none non-j).
- Constants/atomic data (all printed by the oracle, records CONST/LVLH/LVLHE): `lambda_c=2.4263102175e-10`, `kb_mec2=kB/(me c^2)=1.6863720498875446e-10`, `kB=1.3806504e-16`, per-level `Eion_ergs`, `mu_red` (H 0.99945567942448077, He 0.99986292543644184), `ME_scale=1`.
Layout (default production): `neq=15`, `nHIeq=6` (1s,2s,2p,3s,3p,3d), `nHeIeq=7` (singlets 1^1S,2^1S,2^1P gw=1,1,3; triplets 2^3S,2^3P0,2^3P1,2^3P2 gw=3,1,3,5; `indexT=3`, `indexT_no_j=7=nl`), `index_HI=1`, `index_HeI=7`; `X = [Xe | H(6) | HeI(7) | rho]`.
Packed ODE vector (`nres = 2 + 5 + 1 + 4 = 12`, `flag_He=1`): `y[0]=rho`, `y[1]=X[1]` (H 1s), `y[2+k]=X[1+get_HI_index(k)]` with `RESHI = 1..5`, `y[7]=X[iHe]`, `y[8+k]=X[iHe+get_HeI_index(k)]` with `RESHE = 1,2,3,5`. `X[0]` is not in `y`.
Domain/branch behavior to preserve: no checks at all in the Saha functions (NaN/negative inputs propagate); `Xe` clip `1+fHe`; `Xp` clip 1; ground-state overrides; the `z<3500` / `z>=3500` accessor switch changes which formula feeds `Xp, NHeI, NHeII` (continuous to ~1e-3 relative across 3500); `z>=zsRe` closed forms.
Existing-fixture linkage: the accepted `native_fcn_effective.txt` `IN` rows of the four unperturbed states (3400.5, 3400, 3399, 3000) equal this capture's Saha output run in the same standalone-ini mode BIT-FOR-BIT (max relative difference 0.0), which verifies the capture definition; those states are for the standalone `Yp=0.245`, analytic-H cosmology and differ from production by 3e-4 (Xe) to 3e-3 (HeI 1s) at z=3000, so they are NOT the production initializer oracle.
Design-equation verification (scratch script `chunk4/verify_design_equations.jl`, 12 z values, native printed scalars as input): products `X_i = Xe Xp NH f_i` / ground overrides reproduce the native `X` with relative error 0.0 (given the printed `f`); `f_i` recomputed from the constants/levels above matches native to 4.2e-16 (H) and 4.4e-16 (He); packed `y` identical (0.0).

## 4. Stage C (separate, later chunk): sampled-HeI-state switch
`main.CosmoRec.cpp:163-197`: at each step, if `zs<200 || fHe - HeI_Atoms.Xi(0) <= Xi_HeI_switch` (`Xi_HeI_switch=1e-7`, `global_variables.cpp:37`, recorded in the oracle) and `flag_He==1`: `flag_He=0`, `Level_I.X[index_HeI]=fHe`, other He entries `1e-300`, `neq = 2 + nres_HI = 7`, solver re-initialized. It needs an ODE trajectory (`HeI_Atoms.Xi(0)` is the evolving ground population) and changes the state dimension, i.e. a discrete event; it cannot be an oracle without running the solver and is excluded here. (Not to be confused with the Chunk 3 `HeI_Feedback/diffusion` switches.)

## 5. Fixture gap analysis
| output needed | existing committed fixture | verdict |
|---|---|---|
| A: `Xe_Seager`, `X1s`, `Xp`, `XHeI1s`, `XHeII1s`, `NHeI`, `NHeII`, `Te_Tg`, `H(z)` at z | none (`native_camb_cosmorec_thermo.txt` is CAMB-resampled FINAL `x_e,T_b,H` only; no internals, header says so) | insufficient -> new oracle (`SP`, `HZ`) |
| A: Recfast++ preliminary arrays | none | new oracle `history1.txt` (6000 x 7) |
| B: Saha inputs (Xe, Xp, NH, Te, XHeII, ...) | `native_fcn_effective.txt` keeps only `B` (Tg,NH,H,fHe), perturbed `IN`, `F`; not the unscaled inputs, LTE factors, or production cosmology | insufficient -> new oracle (`SAHAIN`, `LTEH`, `LTEHE`) |
| B: per-level X at z | only 4 unperturbed standalone-cosmology states (z>=3000) | cross-check only; new oracle `XLI` (12 z values 800..3600, production cosmology) |
| B: packed ODE vector | none | new oracle `YSOL` |
| B: level data/constants | none | new oracle `CONST`, `LVLH`, `LVLHE`, `LI`, `RESHI`, `RESHE` |
| C | none | out of scope (needs a trajectory) |
| end-to-end final `x_e, T_b` | `native_camb_cosmorec_thermo.txt` (245 rows) | valid ONLY as a far-future gate |
The existing fixtures remain valid for what they document (RHS parity at explicit state/background); no re-capture of their contents is proposed.

## 6. New native oracle (external, read-only, not in Git)
`N4/harness4.cpp` (+`build.sh`, same compile line as chunk3d/3e, linked to the original `libCosmoRec.a`, `libRecfast++.a`, GSL; `#define private public` is used in the harness only to read `Cosmos::n_Xe/zsRe/loaded_Hz` and call the private `Cosmos::recombination_history`; no native source/library touched).
It performs the front half of the batch entry with the parameter assignments of section 1 plus `init_Hubble` with `camb_H_input.txt` (sha in `hashes.txt`), `set_startup_data_CR`, `set_array_dimensions_and_parameters`, `allocate_memory`; `CosmoRec()` (ODE) is never called. The 6000-row Recfast history output's fields `Xe_H` and `Xe_He` are ionized-electron abundances (`n_e,H/n_H` and `n_e,He/n_H` respectively), not neutral-atom populations. It records:
CFG/COS/LI/RESHI/RESHE/HEIDX/HLVL (configuration, derived cosmology, layout), CONST/LVLH/LVLHE (constants, level data), HZ (17 probes of `cosmos.H` incl. low z and the 1e4 edge), SP (21 z values 25000..1 of all accessors), SAHAIN/LTEH/LTEHE/XLI/YSOL (12 z: 3600, 3500, 3400.5, 3400, 3399, 3000, 2800, 2500, 2000, 1500, 1200, 800, full inputs and outputs), and `history1.txt` (the arrays fed to the splines, from the same call).
Reproducibility: two fresh-process runs byte-identical (`probe1==probe2`, `history1==history2`; sizes 25148 B and 748714 B; ~0.1 s per run). Hashes of harness, binary, build script, H input, outputs, libraries and ALL read sources (CosmoRec, CAMB, repo fixtures) in `N4/hashes.txt`.
Optional `ini` mode (`harness4 H out hist <ini>`) reproduces the standalone chunk3d startup solely for the bit-for-bit cross-check above.
Key native values (production cosmology): z=3000: `Xe_Seager=1.0819197658689494`, `X1s=7.4189054899e-10`, `XHeI1s=8.9556506e-05` (He starts as HeII), `NH=5118.6917333354`, `Te=8179.2234489`, `TCMB=8179.2255`; z=1500: `Xe=0.96146409307`, `X1s=0.0447558886`, `XHeI1s=0.0757893415`.

## 7. Proposed staging (smallest safe next step)
- 4a (recommended next, small, fixture-first): explicit-input Julia `saha_levels` / `pack_ysol` for H and HeI with the atomic constants and level tables as explicit data, validated bit-for-bit on the production-cosmology oracle (inputs `Xe, Xp, NH, Te, Tg, XHeII, XHeI` taken from the oracle) AND the four chunk3d `IN` rows (their inputs also required: re-capture from `harness4 ... ini` is already available). Pure algebra: no solver, no spline.
- 4b: Cosmos accessors from tabulated history nodes: natural cubic spline in z (and ln X variants), 1e-14 end nudging, z>=3500 Saha branches, closed forms for z>=zsRe, and the loaded-H cubic spline in log-log including the unset `lHz[0]` node. Oracle: `history1.txt` + `SP`/`HZ`. Still no ODE.
- 4c: Recfast++ preliminary history (own stiff solver, ~2.8k lines; DM/B-field/variation-of-constants/reionization modules disabled by production). Large; needs its own design and an oracle on solver-state outputs; parity of the 6000 nodes to the original solver tolerance is not guaranteed by any SciML solver, so it must be reviewed before starting.
- 4d: the Stage C helium switch inside the later ODE chunk.
Integration-test sequence: 4a unit -> 4a + `copy_LI_to_ysol` packing -> 4b accessors vs `SP` -> 4a fed by 4b vs `XLI` -> (later) 4c history vs `history1.txt` -> (much later) ODE vs `native_camb_cosmorec_thermo.txt`.

## 8. AD and boundary tests (for 4a/4b, after approval)
Differentiable inputs: `Xe, Xp, NH, Te, Tg, XHeII, XHeI` (4a); history node values and `z` (4b, through the cubic spline coefficients). Tests: ForwardDiff vs 256-bit central differences (as in Chunk 3), directional derivatives, Mooncake prepared VJP reused across parameter changes with two independent preparations.
Branch/boundary cases NOT to be called differentiable: the `Xe` clip `1+fHe`, `Xp` clip 1, ground overrides (they zero the Saha derivative of that entry), z=3500 accessor switch, z=zsRe switch, end-nudge branches of the splines, the spline out-of-range abort (must become an explicit error; never silently extrapolate), unset `lHz[0]` region (z < ~1). `exp(Eion/(kB T))` overflow/underflow for extreme T is native behavior (no guard) and must be tested as such.

## 9. Solver/background dependencies
4a needs no solver and no background (all scalars explicit). 4b needs only the tabulated nodes. Anything that produces those nodes (Recfast++, `H(z)`, Omega_L closure, `rho_g`, `Nb0`) is Stage 4c. `fac_mHemH`, `const_*` physical constants must be taken from `physical_consts.h` (values printed in the oracle where used).

## 10. Limitations / open items for review
1. The oracle covers ONE cosmology; parameter dependence is to be probed only after 4a/4b are accepted.
2. The H(z) table is a Python/CAMB-generated input; it equals what the adapter passes up to CAMB's own arithmetic (`h_of_z/MPC_in_sec`), not byte-verified against the Fortran `dtauda` expression.
3. The `lHz[0]` hazard is reproduced but not yet traced to an effect on any recombination output; it is expected to be negligible at z>=100 (<1.4e-8 relative `H`).
4. `Rec.A2s1s=8.2206` after startup is observed, not explained here (affects Recfast++/rates, not Stage B).
5. Committed fixture for 4a/4b would be derived from `probe1.txt`; the 6000-row history (749 kB) is proposed as a decimated subset (every node in z in [2900, 3600] plus every 20th) unless the reviewer wants it whole.
6. Stage C and everything after are intentionally untouched.

## 11. Amendments after the 4a implementation (see `docs/CHUNK4A_RESULTS.md`)
- API as built: `SahaInputs(fHe, Xe_Seager, Xp_raw, NH, Te, Tg, XHeII, XHeI1s)` (the last two are the native quotients `NHeII/NH`, `NHeI/NH` as printed by the oracle), `SahaLevels`/`SahaConstants` explicit data, `saha_initial_state`, `pack_ysol` (flag_He >= 1 branch only); no new atom abstraction beyond `SahaLevels` (the existing `HydrogenAtom`/`HeliumAtom` structs lack `Eion`).
- The 4a Julia tests use the production-cosmology fixture only. The four unperturbed chunk3d `IN` rows were compared bit-for-bit with the harness in section 3 (done at design time in C++), but are NOT re-tested in Julia because their scalar inputs are not part of the committed chunk3d fixture; the "AND the four chunk3d IN rows" clause of section 7 is therefore withdrawn for 4a.
- `SAHAIN` prints the CLIPPED `Xe`; the clip is inactive at all 12 native states (tested), so the unclipped value needed by helium equals it. The clip/branch edges are tested with synthetic explicit inputs, not a native oracle.
