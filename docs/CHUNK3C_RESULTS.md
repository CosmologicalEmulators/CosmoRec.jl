# Chunk 3c results: default H-I absorption of HeI photons (standalone absorption gate)

**NOTICE (original-use acknowledgment).** The reference values come from outputs of the original CosmoRec v3.0b (J. Chluba et al.; native SHA
`086769055f61ae0c244a53dd381ee65b624d0ac3`, clean and read-only). Use must be acknowledged; cite Chluba & Thomas 2010 (MNRAS 412, 748),
Chluba & Sunyaev 2006 (A&A 446, 39) and Rubino-Martin, Chluba & Sunyaev 2008; also considered: Chluba, Vasil & Dursi 2010; Switzer & Hirata 2008;
Grin & Hirata 2010; Ali-Haimoud & Hirata 2010. Bugs: Jens@Chluba.de. The fixture is a small text subset of derived values; no table is bundled.

## Scope (what is and is not claimed)
Accepted for review: ONLY the H-I absorption block of `fcn_HeI_effective` (`Modules/ODEdef_CosmoRec.cpp:284-321`) as a standalone, pure Julia gate.
This is **not** a coupled H+He RHS, not initialization, not an ODE solve and not default-run correctness. No "complete first-pass H+He RHS" is claimed;
that needs a total H/He fixture exercising the production `fcn_effective` flow.

Excluded / not ported: the native explicit `DPesc_coh` integral used when a query leaves the pre-tabulated DP table (Julia throws
`DPTableDomainError`); HeI diffusion (the block is off while `Diffusion_correction_HeI_is_on`; exposed only as a flag); feedback; hydrogen RHS coupling.

## Active defaults (exact)
- CAMB/CosmoRec runmode 0 sets `HI_absorption = 2`; `setup_functions.cpp:411` remaps it to `1` with `_HI_abs_appr_flag = 1`.
- Branch: `HI_absorption==1 && z <= zcrit_HI(3400) f_t && !Diffusion_correction_HeI_is_on`. Intercombination additionally `z <= zcrit_HI_Int(3400) f_t && spin_forbidden==1` (default 1).
- The approximation flag applies `fcorr_Loaded(Tg)` (natural cubic spline of `f.corr.dat`, argument `Tg/2.725-1`, clamped) to the singlet `Pesc` only. The triplet has `pd = 1`, no fcorr.
- `f_a = f_b = f_t = 1` in the fixture; non-unit `f_t`, `f_b` are implemented but NOT validated natively.
- Fixture table: `DP_Coll_Data.31.fac_50.neff_30.dat`, 150 sheets, `lnT` descending (about 4120-9404 K), 31 x 31 per sheet and channel.

## Native-to-Julia source map
| Native | Julia (`src/HIAbsorption.jl`) |
|---|---|
| `ODEdef_CosmoRec.cpp:284-321` branch, slots and signs | `hi_absorption_rhs!` |
| `ODEdef_Xi.HeI.cpp:377-475` `evaluate_HI_abs_HeI`, `_Intercombination` | `hi_abs_singlet`, `hi_abs_triplet`, `line_channel` |
| `DPesc_HI_abs_tabulation.cpp:34-50` `pd` + `get_Bitot_HeI` (`get_effective_rates.HeI.cpp:440-486`) | `pd_singlet`, `helium_Bitot`, `BitotSeries` |
| `DP_interpol_S/T` (`:74-...`), `_eta`/`_tau` cubic, T Lagrange | `dp_correction`, `dp_lookup`, `sheet_dp`, `dp_sheet_start` |
| `read_DP_Data` (`:686-745`) | `read_native_dp_table` (explicit, no cache) |
| `fcorr_Loaded` (`ODEdef_Xi.HeI.cpp:216-292`), GSL cspline | `FcorrSpline`, `read_native_fcorr`, `fcorr` |
| `Sobolev.cpp` `p_ij`, `tau_S_gw`, `exp_nu` | `sobolev_p`, `tau_S_gw`, `exp_nu` |

Formulas: `d = Dpij A21/(efac-1) (Xi efac - Xj gw/gwp)`, `gw=3`, `gwp=1`, `efac = exp(min(700, h_kb Dnu/Tg))`,
`tauS = A lambda^3/(8 pi H) |Ni (Nj/Ni gwi/gwj - 1)|`, `eta = cl NH XH1s sig_c / H`. Signs: HeI 1s `+d`, 2^1P (singlet) or 2^3P1 (triplet) `-d`,
electron `+d`, H 1s `-d`; all updates accumulate (`+=`). `Dpij = Pesc - p_ij(tauS)`, `P = p_ij(pd tauS) + DP`, `Pesc = fc pd P/(1-(1-pd)P)`,
table coordinates `(ln T f_t, ln eta, ln(tauS pd))`.

## Fixture provenance
`test/fixtures/native_hiabs_components.txt`: 394973 B, sha256 `f9c669f772431d0aa8b43ce7101587ccec84481c38b295760c56a6594f2d902b` (text, 17 digits).
Produced by an external, read-only C++ harness (`.../chunk3c/native_capture/harness3c.cpp`, sha256 `2b69cf44829674172d574c079d30211e1862a237412db548e9e96ed3cff7fdcf`;
binary `0d97816e5647449bf23b1593ed5410706aa274d14da3079aa23e2a91b2c95291`; build cmd `21482c01c216414296be75efcb0402dc4baac10547bbff474223580443a1fa5e`;
`cases.txt` `aaf81bb428ef8d05ef1bdd7d815b52d6580f8b9bd965258965be5ece6383116d`) linked against the unmodified static libraries
(`libCosmoRec.a` `2ac075e8...b783`, `libRecfast++.a` `441a3eaa...a5fc`). Source/table hashes (in the fixture header): DP table `ad0d8271...ee77f85`,
`f.corr.dat` `b25e3324...be269`, `ODEdef_CosmoRec.cpp` `6ad6db99...e3bd`, `DPesc_HI_abs_tabulation.cpp` `ab64f7ef...19ea`, `ODEdef_Xi.HeI.cpp` `27b18598...f3f`,
`routines.cpp` `2d71d672...3bd7`, `Sobolev.cpp` `894c978a...a7f8`. Two fresh native processes gave byte-identical output (`run1.txt == run2.txt`).
Native checkout left clean; no native edits. The native-edit-free harness could call the exported absorber interpolation functions directly.

Content: 63 cases (history points over the table range, zcrit at -1e-6/0/+1e-6/+20/+40, nonequilibrium and equilibrium populations, XH1s and X0 scalings,
T-node, node+-, mid-cell and 0.0001/0.9999 cell fractions, four OUT_* off-table cases); per case the inputs, intermediates (`eta`, `tauS`, `efac`, `pd`,
Bitot, fcorr), the 4-sheet stencil window of DP, the Bitot stencil, per-channel native outputs (`DP_interpol_S/T` and `dXe`) for both approximation flags,
and 504 combined-vector records (flag x spin x {zero-start, sentinel-start}). Inputs are not declared outputs.
**Limitation of the combined records:** the `K` vectors are a transcription of the branch/sign statements of lines 284-321 using native callee outputs
(`fcn_HeI_effective`'s `parameters` global is file-local, so the function itself could not be called); the callees (interpolators, pd, fcorr) are the real native code.
Cases with `inTable = 0` (9 fully off-table, 4 with only the triplet off) were not evaluated natively (native would run `DPesc_coh`); their `V` fields are NaN and the
Julia tests assert a `DPTableDomainError` instead. `fallback.txt` (informational native explicit-integral values) is not part of the package.

## Numerical precision (Julia vs native fixture, worst relative error over all cases)
| Quantity | Max rel. error | Quantity | Max rel. error |
|---|---|---|---|
| fcorr | 2.1e-16 | `DP_interpol_S`, flag 1 / flag 0 | 2.6e-13 / 6.5e-16 |
| eta, efac, Bitot, pd | 0 (bitwise) | `dXe_S`, flag 1 / flag 0 | 5.3e-13 / 3.3e-14 |
| tauS_S | 3.5e-16 | `DP_interpol_T`, `dXe_T` | 3.4e-16, 3.6e-16 |

Test tolerance is 1e-9 relative for DP and dXe (observed: 5e-13; the singlet flag-1 error comes from the cancellation `Pesc - p_ij` times `fcorr`), and 1e-12 to 1e-14 for intermediates.
Combined vectors: relative 1e-9 of the increment, plus an absolute floor of `4 eps` times the sentinel (1.25e-3) for the sentinel-start records. Nothing was tuned.

## Tests
- `test/chunk3c_hiabs_native.jl`: 8870 tests, unconditional, no skipped branches. Every fixture output compared elementwise (intermediates, per channel for both flags, all 504 combined vectors),
  output-index mapping (slots 1, 2, 3, 5, 8 and untouched slots 4, 6, 7, 9; relocated indices), signed transfer conservation, zcrit at -1e-6/0/+1e-6, diffusion and spin switches,
  accumulation, p_ij branches, `exp_nu` cap, a synthetic cubic table (reproduction at nodes and off-node, edge clamps, all domain errors incl. NaN), Bitot and fcorr limits.
  The sparse windowed tables are NaN outside the captured stencils so any wrong index fails.
- `test/chunk3c_hiabs_ad.jl`: 462 tests (see below). Included from `test/runtests.jl` after the 3b files.
- Full `Pkg.test()` (single process, log `chunk3c/pkgtest.log`): **18094 / 18094 passed**, 16m50s (8762 prior + 9332 new). The laptop was throttling (98 C at the end; the human approved continuing). Focused runs: native 1.7-3.1 s, AD about 50 s.

## AD tests and tolerances
ForwardDiff Jacobians vs 256-bit central differences (h = 1e-20 relative, per-column scaling), independent directional derivative, and Mooncake VJPs
(`prepare_gradient(obj, AutoMooncake(; config=nothing), x, Constant(w))`; each prep reused at 4 points x 3 seeds, two independent preps compared, `gA == gB`) vs ForwardDiff `J' w`.
Differentiated inputs: `Tg, NH, Hz, XH1s, X[1..7]`, plus a test-defined (Tcmb, omega_b, h) scaling of the captured background. Table values/axes are constants.

| Check | Worst observed | Tolerance |
|---|---|---|
| FD (256-bit) singlet A/B/C | 1.2e-13 / 1.8e-12 / 1.5e-14 | 1e-8 |
| FD triplet A/B/C | 2.0e-9 / 3.1e-14 / 1.8e-14 | 1e-8 |
| FD combined A/B/C | 2.0e-9 / 4.7e-13 / 1.5e-14 | 1e-8 |
| FD cosmological scalars | 2.0e-9 | 1e-7 |
| directional derivative | 1.3e-14 | 1e-12 |
| Mooncake VJP vs J'w (all) | <= 4.2e-16 | 1e-9 |
| FD, Tg column, seam points | 3.2e-14 | 1e-6 |

The 2e-9 (case A) is the Float64 roundoff of `p_ij = (1 - e^-tau)/tau` at small tau (native form, no `expm1`), as in 3b. Tolerances are not tuned; 1e-8 was chosen before running.
Independent invariants: rows 4, 6, 7, 9 identically zero; electron row equals HeI-1s row and minus the H-1s row; HeI column sums vanish; columns of X[2], X[4], X[5], X[7] vanish;
electron production equals the sum of the 2P losses; closed-form `dp_correction` and `d/dtau` for a DP = 0 table; derivative through a DP linear in `ln eta`; `pd = 1, fc = 1` gives exactly 0.

Branches (one-sided only; smoothness across them is NOT claimed): z vs zcrit (active side differentiable, inactive side exactly zero gradient for ForwardDiff and Mooncake;
jump = the whole value, 2.8e-12 in dXe); T-sheet stencil seam (points +-1e-6 off a node, each differentiated one-sided; the Tg-derivative differs across the seam by a relative **8.8**,
i.e. the interpolant is C0 but not C1 across sheet nodes on this table); table cell seams in eta/tau; p_ij branches at tau 1e-10 and 500; native clamps; spin switch.
Domain errors (`DPTableDomainError`) are asserted at eta/tau/T edges and NaN.

## Benchmarks (BenchmarkTools, `benchmark/chunk3c_benchmarks.jl`, explicit args, `evals = 1`; log `chunk3c/bench_run2.log`)
| Operation | Median | Alloc |
|---|---|---|
| fcorr | 34 ns | 0 |
| helium_Bitot | 40 ns | 0 |
| dp_lookup (singlet) | 121 ns | 0 |
| hi_abs_singlet / triplet | 187 / 147 ns | 0 |
| hi_absorption_rhs! (both) | 317 ns | 0 |
| ForwardDiff.jacobian 9x11 | 1.77 us | 4.3 KiB |
| Mooncake prepare_gradient (warmed; first sample compiles) | 254 us | 587 KiB |
| Mooncake hot gradient (prepared) | 12.4 us | 1.75 KiB |
Timings were taken on a throttling laptop and are indicative only.

## Known limitations
1. The out-of-table fallback `DPesc_coh` is not ported. The table is centred on a reference history (`HI.HeI.1s.dat`) with +-50x bands in tau and eta, so histories far from it
   (large parameter shifts) will throw `DPTableDomainError`; for AD over cosmological parameters this limits the usable range. This is the largest gap for production use.
2. Native `exit`s at Tg above about 9430 K inside `pd` (HeI table limit); Julia throws `DPTableDomainError`.
3. The combined `K` vectors are a transcription of lines 284-321 with native callees, not a call to `fcn_HeI_effective`.
4. `f_t`, `f_b != 1` are implemented but not natively validated. Only one DP table (`neff_30`) was captured; the `fac_50` / `31` header variants were not.
5. Only stencil windows are bundled; correctness of the full-table reader `read_native_dp_table` and `read_native_fcorr`/`read_native_helium_bitot` is checked structurally (small synthetic tests and the windows), not against a full native table dump.
6. Non-smooth points: z=zcrit, T-sheet and eta/tau cell seams (derivative jump measured above), p_ij tau branches, native clamps.
7. Harness header contains two stray lines (`HDRIDENT`, `0`) in the hash list; harmless to the parser.

## Budget
USD about 5.5 of 12 spent at the time of writing (about 6.5 remaining). Stopping here for the main-agent source/AD audit; no initialization or ODE work was started. No commits or stages were made.
