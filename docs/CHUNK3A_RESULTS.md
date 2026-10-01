# Chunk 3a: smooth hydrogen effective-rate / Compton / population RHS components in pure Julia

Scope: the right-hand-side **components** of the default 3-shell hydrogen effective model, evaluated at a given state. This is **not** the coupled complete
default CosmoRec ODE and nothing here solves or differentiates an ODE. Nothing is committed or staged.

## NOTICE and citations (CosmoRec licence)

The fixture contains outputs of **CosmoRec** (J. Chluba et al., v3.0b, 2017). Its README licence requires that use is acknowledged in publications and that
Chluba & Thomas 2010, MNRAS 412, 748 and Chluba & Sunyaev 2006, A&A 446, 39 are cited. Also considered by the authors: Chluba, Vasil & Dursi 2010,
MNRAS 407, 599; Switzer & Hirata 2008, PRD 77, 083006; Grin & Hirata 2010, PRD 81, 083005; Ali-Haimoud & Hirata 2010, PRD 82, 063521;
Rubino-Martin et al. 2010, MNRAS 403, 439. No guarantee for the correctness of the outputs is given upstream. The native sources and binaries are untouched
(`git status` of the native checkout is clean, SHA `086769055f61ae0c244a53dd381ee65b624d0ac3`). Only the small previously accepted rate-table window is used;
no further table data were copied.

## Reference-first gate: native oracle (done before any Julia physics)

All needed functions are exported by `libCosmoRec.a`, so **no access blocker**. An external harness (outside git; not compiled in this repo) calls the original
functions directly and does not re-implement any equation:

* `ODE_effective::evaluate_TM`, `evaluate_2s_two_photon_decay`, `evaluate_Ly_n_channel` (n = 2, 3)
* `ODE_HI_effective::evaluate_effective_Rci_Ric_terms`, `evaluate_effective_Rij_terms`
* `fcn_HI_effective` (native combined hydrogen part; `Diffusion_correction_is_on = 0`, quadrupole lines off)
* native `get_rates` (rates fed to every component) on the production atom `Gas_of_Atoms(3,1,1.0,false,0,-2)`, `nS_effective = 500`.

Harness and logs: `.../cosmorec_differentiability_20260930/chunk3a/native_capture/` (`harness3a.cpp` sha256 `08839f08...b4e1`, `build_cmd.txt`, `compile.log`, `header.txt`, `cases.txt`, `run1.txt`/`run2.txt`).
Compiler: g++ 13.3.0, `-Wall -pedantic -O2 -fPIC`, include dirs = all `Development/*` and `Modules` header dirs + `Rec_database/Effective_Rates.HI` + `camb-cosmorec-env/include`,
static link `libCosmoRec.a libRecfast++.a libgsl.a libgslcblas.a -lm` (GSL 2.8). Library hashes are in the fixture header.
Two fresh-process runs were byte-identical (`cmp`). Each component is called with nonzero **sentinel** "before" buffers, so accumulation (`+=`) versus assignment (TM) and the
before/after values per component are in the fixture, not just the sum.

Fixture `test/fixtures/native_hrhs_components.txt`: 41670 bytes, sha256 `03aaffdfdc877f1463c147fda5de7ab54acede90ff3c33abd08178e845ae2e4b`.
Rows: `G` (constants and atomic data read from native objects), `C` (case inputs: z, Tg, Te, rho, Xe, fHe, Xp, NH, Hz, X[1s,2s,2p,3s,3p,3d]), `RT` (native rates A, B, R per resolved state),
`K` (component, before[6], after[6]; TM scalar). 16 cases: 7 model-domain states (z = 200..2100, Te/Tg from 0.95 to 1.02, Tg = 548..5726 K) and 9 states at the accepted rate-window query
points (Tg from 200 to 5452.75 K, both sides of the detailed-balance switch at 5452.725 K), so the table-driven path also has a native oracle.
Populations and background (Hz, NH from a flat LCDM-like formula) are **harness-chosen inputs recorded per row**, not a solved history.

## Source line map (`CosmoRec` v3.0b)

| piece | native | Julia (`src/RHS.jl`) |
|---|---|---|
| Compton + adiabatic, `drho/dt = 8/3 sigT c rho_g(Tg) Xe/(1+Xe+fHe) (1-rho) - H rho`, `rho_g_fac = pi^2/15 kB (kB/(hbar c))^3/c^2/me` | `ODE_effective.cpp:85-98`, `:43` | `matter_temperature_rate`, `rho_g_fac` (assigns) |
| 2s-1s, `A2s1s (X2s - X1s e^{-h_kb nu21/Tg})`, `Lambda_ind = 1` | `ODE_effective.cpp:131-162` | `two_photon!` |
| Lyman Sobolev, `tau = A l^3/(8 pi H)(3 X1s NH - Xnp NH)`, `p = (1-e^-tau)/tau`, w = 3, n = 2..3 | `ODE_effective.cpp` (`evaluate_Ly_n_channel`), `:52-63` | `lyman_channel!`, `lyman!` |
| Rci/Ric, `B_m (A_m Xe Np - X_m)`, `Np = NH Xp` | `ODE_effective.cpp` `evaluate_effective_Rci_Ric_terms` | `continuum!` |
| Rij, i<j, `R (X_i - gw_i/gw_j e^{-h_kb nu_ul(n_i,n_j)/Tg} X_j)` | `ODE_effective.cpp:390-410` | `interlevel!` |
| combined hydrogen part | `ODEdef_CosmoRec.cpp:140-210` (`fcn_HI_effective`) | `hydrogen_rhs!` (two methods: given rates, or explicit `AtomicRateTable`) |
| ground-state reconstruction `xp = 1 - X1s`, `xHeII = fHe - XHeI1s`, `xe = xp + xHeII` | VERIFIED_ANALYSIS (ground populations) | `proton_fraction`, `electron_fraction` |

Level order (native `Get_Level_index = l + n(n-1)/2`): 1s, 2s, 2p, 3s, 3p, 3d; resolved states `(2,0),(2,1),(3,0),(3,1),(3,2)` map to `X[m+1]`.
Note the native sign convention: `nu_ul(n_i, n_j)` with `n_i < n_j` is negative (the source comment "also the sign is important"), so `exp(-xij) > 1` for 2->3 pairs and `xij = 0` within a shell.
The RHS time variable is cosmic time (1/s) exactly as native; division by `dz/dt` is outside these terms.

Atomic data (`gw`, `nu_ul`, `Dnu`, `A21`, `lambda21`, `A2s1s`, constants) are Float64 literals equal to the values read from the native objects; the fixture test checks equality with its `G` rows.
`sigT` is the native global `const_sigT` as seen by the harness before any `CosmoRec.cpp` initialisation (6.6524585583988891e-25 cm^2).

## Results

### Native comparison (`test/chunk3a_rhs_native.jl`, `rtol = 1e-13`, `atol = 0`)

Worst observed relative difference over all 16 captured states (sentinel buffers included): **0.0 (bit-identical)** for TM, 2PH, LY2, LY3, RCI, RIJ and FCN_HI;
the table-driven combined path (`hydrogen_rhs!(..., table, ...)` vs native `get_rates` + `fcn_HI_effective` at 9 states) is bit-identical except one state at **2.4e-16**.
The Julia pieces applied sequentially equal the native combined call. Native rates `RT` equal Julia `get_rates` to 1e-13 at the window queries. No tolerance tuning.

### Algebraic invariants (independent of the fixture)

Opposite transfers (2PH, Lyman, Rij) and population conservation; the sum of the combined RHS equals the sum of the continuum terms only; the exact state order; the ground-state
x_e relation (and that it is **not** `1 - sum(all populations)`); thermal sign and units (against an independent `a_rad` closed form, rtol 1e-12); detailed-balance equilibria give zero net flow.

### AD (`test/chunk3a_rhs_ad.jl`; point `tbl_int_a`: Tg = 3000 K, Te/Tg = 0.9, interior of the accepted table window; populations multiplied by (1, 1.1, 0.9, 1.2, 0.8, 1.05) so no output cancels exactly)

Each mapping is a fresh-output wrapper of the RHS function. Worst observed errors (relative to the largest entry of the column / gradient):

| component (inputs) | ForwardDiff J vs 256-bit central difference | directional derivative vs J v | prepared Mooncake VJP vs ForwardDiff J'w |
|---|---|---|---|
| TM (rho, Tg, Xe, fHe, Hz) | 2.1e-16 | 0 | 2.6e-24 |
| 2PH (Tg, X) | 4.1e-16 | 1.6e-16 | 4.4e-15 |
| LY (Tg, NH, Hz, X) | 5.2e-16 | 1.2e-15 | 9.6e-16 |
| RCI (Xe, Np, X, A, B) | 1.7e-16 | 7.7e-19 | 2.5e-16 |
| RIJ (Tg, X, R upper) | 4.3e-16 | 3.2e-16 | 3.5e-16 |
| combined from rates (Tg, Xe, Xp, NH, Hz, X, A, B, R) | 4.3e-16 | 4.4e-16 | 2.9e-16 |
| combined through the table (Tg, rho, Xp, fHe, NH, Hz, X; i.e. rate query coordinates Tg, Te = rho Tg) | 2.9e-12 | 7.3e-14 | 2.1e-15 |
| combined, cosmological scalars (Tcmb, omega_b, Omega_m, h) via a **test-defined** background, X, rho, Xp, fHe | 1.8e-11 | 6.0e-13 | 3.2e-16 |
| detailed-balance switch, one-sided (just below 5452.725 K and first above) | - | - | 4.2e-16 |

The finite-difference tolerance is 1e-12 for algebraic components and 1e-10 for the two table-driven mappings (observed 2.9e-12, 1.8e-11); VJP agreement is tested at 1e-9. The output is multi-output for each seed:
3 random weight vectors w at the base point plus 3 varied points (relative perturbation ~1e-3 of every input) per component, i.e. 12 seeded projections each.
**Cache reuse:** one `prepare_gradient` per mapping is reused at all 4 points and 3 weights (the weights enter as `Constant`), a second independent preparation is made for the same objective and the
two give identical gradients (`==`). The test background is a local formula (not the production cosmology) used only to show differentiability w.r.t. active cosmological scalars.
Derivatives are piecewise through the rate lookup (seams, DB switch, eps clamp): the switch was probed one-sided only and nothing was compared across it.

### Known limitation (tested, not hidden)

Tg below about 55.6 K (z < ~19): native `A` of the n = 2 states overflows to NaN. At Tg = 31 K the primal RHS is NaN in the 2s and 2p components only; the ForwardDiff Jacobian has NaN/Inf entries and
Mooncake's gradient of a finite output (1s) is NaN in several entries. I did **not** verify whether the native ODE callback can reach that range (the multilevel integration is expected to stop at the Recfast handoff
well above it, but that is unverified here), so the tests use points above it.

### Gate

Focused: native 272/272, AD 323/323 (`julia --project=<chunk2/test_env>`). Full `Pkg.test()` (PTY, no skipped or error-ignored paths, no env gate): **6564/6564 passed, exit 0**, 27m09s wall (5969 previous + 272 + 323), log `.../chunk3a/final_full_pkgtest.log`.

## Explicit gaps (not translated here)

Helium (HeI and HeII equations, `fcn_HeI_effective`, `XHeI` states); H absorption (`HI_absorption = 2`); diffusion correction (`Diffusion_flag = 1` in the full default) and DI1/DI2 corrections;
quadrupole lines; `rescale_rates` (`runmode > 0`) and the `f_t` scaling of Tg before the rate lookup; initialisation/Saha; Recfast handoff; the PDE; the `df_d*` native Jacobian functions; the
background/cosmology module (Hz, NH are inputs); switches; the ODE solve and its differentiation; free-electron balance (`Xe`, `Xp` are inputs of the components, and the ground-state reconstruction is provided separately).
The combined function is therefore **not** the complete default ODE.

## Files

`src/RHS.jl` (exported through `src/CosmoRec.jl`), `test/chunk3a_helpers.jl`, `test/chunk3a_rhs_native.jl`, `test/chunk3a_rhs_ad.jl`, `test/runtests.jl`, `test/fixtures/native_hrhs_components.txt`, this document.

## Thermal notes

Heavy JIT/Pkg.test drives the package to 98-100C within seconds on this machine. Our own `Pkg.test` was paused/resumed (SIGSTOP/SIGCONT on its own process tree only, by `chunk3a/thermal_watchdog.sh`,
pause at >= 92C, resume at <= 72C; log `chunk3a/thermal_watchdog.log`, 339 pause/resume cycles, which is why the wall time is long). No other process was touched.
