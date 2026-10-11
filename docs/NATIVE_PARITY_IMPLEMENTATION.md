# Native parity implementation: one-to-one port of ORIGINAL CosmoRec 0867690 (CAMB runmode 0, defaults)

Direction (Marco, 2026-10-11): reproduce the original algorithm itself, not the same equations solved by different integrators. Bitwise where feasible; otherwise an explicit ulp-level difference with a proven reason. There are no new references or acceptance thresholds. The SciML callbacks remain as explicit diagnostics, not as the faithful default.

## Breakdown (native source → Julia), status

| # | native component | source | Julia | status |
|---|---|---|---|---|
| 1 | Cosmos closure constants | `Cosmos.cpp:207-335,1141`, `Cosmos.h:326-333` | `StandaloneInit.jl` `cosmos_constants` | **bitwise** at 10 cosmologies (chunk35 fixture) |
| 2 | Recfast++ node grid / Saha segments | `Recfast++.cpp:119-273` | `RecfastPP.jl` `recfast_grid` | **bitwise** (fid all 6000 nodes, h50 subsample) after the pivot-redshift order fix: native `zc *= 2.726/T0` (`Recfast++.cpp:151`) |
| 3 | Recfast++ ODE solver: Gear/BDF orders 1–4, variable step, predictor/corrector coefficients, Newton with 4-point numerical Jacobian, Gaussian elimination (neq = 4), error/step controller, node driver | `Development/Cosmology/Recfast++/src/ODE_solver.Recfast.cpp:276-1011`; driver `Recfast++.cpp:284-325` | **new** `RecfastGear.jl` (`RFGearSolver`, `rf_gear_solve_history!`, `recfast_gear_solve`) | **fid: whole preliminary history bitwise** (Xe_H, Xe_He, Xe, dXe, dXe_H, TM; 6000 nodes). **h50: all states bitwise**; dXe/dXe_H differ at 25 of 1201 subsample nodes (z < 1.6, <= 99 ulp), traced to item 6b |
| 4 | Recfast++ RHS `evaluate_Recfast_System` + helpers (Boltzmann, SahaBoltz, alphaH/He, betaH/He; NH, f_He, TCMB) | `src/evalode.Recfast.cpp:52-220`, `cosmology.Recfast.cpp:51-120` | **new** `evaluate_recfast_system_native` (`RecfastGear.jl`), Float64, exact native operation order, glibc exp/pow; the generic AD form `recfast_rhs` is kept for Dual types | **bitwise** vs native at all 6000 nodes x 4 components, fid and h50 (native Hz fed back; `native_parity/rhs_parity.log`) |
| 5 | Recfast++ tail `Xe_frac_rescaled` (`fcn_rescaled`) | `Recfast++.cpp:336-441` | `RecombinationODE.jl` `recfast_tail` (SciML callback) | to do: drive with `RFGearSolver` (same solver; `fcn_rescaled` multiplies f[1] by ff) |
| 6b | Hubble table spline (`Hubble::init/H`, GSL cspline + `calc_spline_JC`) | `Cosmos.cpp:104-193`, `routines.cpp` | `hubble_table`/`cosmos_H` (glibc log/exp now on Float64) | fid bitwise at all 6000 probe z; **h50: 131/6000 H values differ (0.01 < z < 1.62, <= 127 ulp)** near the native node-0 ln H = 0 spike: GSL tridiagonal-solve / coefficient arithmetic not yet matched. NEXT |
| 6 | Cosmos splines of the preliminary history | `Cosmos.cpp:390-483`; GSL cspline | `CosmosAccessors.jl` | already GSL-algorithm port; re-verify bitwise once 3–4 are bitwise |
| 7 | Multilevel ODE solver (CosmoRec main passes): Gear with BiCGSTAB/Gauss, tolerances `main.CosmoRec.cpp:16-36` | `Development/ODE_PDE_Solver/ODE_solver_Rec.cpp` (~1500 lines) | currently SciML Rodas callbacks | to do (largest item): port like item 3; capture native step/Jacobian intermediates via scratch harness |
| 8 | Helium switch, Saha init, pass node grid | `main.CosmoRec.cpp:93-268`, `HI_routines.cpp:150-176`, `HeI_routines.cpp:140-177` | `HeliumSwitch.jl`, `SahaInit.jl`, `recombination_pass` | ported; re-verify bitwise once 7 exists |
| 9 | PDE stages (march `Step_PDE_O2t`, integrals, Patterson) | `PDE_solver.cpp:386-622`, `Solve_PDEs*.cpp` | `HIPDE*.jl` | ported (step assembly audited equal, p3_next); bitwise verification against a native per-step capture still to do |
| 10 | Output assembly / return to CAMB grid | `CosmoRec.cpp:545-597` | `return_solution_to_grid` | ported |

## Milestone 1 (this chunk): what exists

- **`src/RecfastGear.jl`** (new): faithful port of item 3. Native data flow is kept:
  - `Snptr[0..5]` pointer rotation;
  - the `Snew`/`Sguess` swap;
  - `ODE_solver_jac` 4-point formula with `Dyj = y*1e-12`;
  - Gaussian elimination (row-major, in place);
  - `ODE_Solver_error_check` on components 0–2 only;
  - `n_down` logic; order increase after `order+1` accepted steps with cap 4;
  - `Dz_z_max = 10/zs`; `Dz_z_last` carried between node calls;
  - the 4th state `t_cos`, initialised to `1/2/H(zs)`.

  `pow` goes through the C library (`_cpow`). `recfast_gear_solve` is the `recfast_history` callback (it returns the solver's `Sz.dy`, i.e. the RHS at the accepted state, as native).
- **`src/RecfastPP.jl:169`**: the pivot-redshift operation order now matches native (translation bug: Julia computed `4800*2.726/T0`, native computes `4800*(2.726/T0)`; 1 ulp at the break node, which shifted every ODE node).
- **`test/native_parity_recfast_gear.jl`** (new, unconditional, in `runtests.jl`). Exact facts only:
  - fid grid bitwise (6000 nodes);
  - Saha segments bitwise;
  - Gear ODE states bitwise from the ODE start through node 3485;
  - h50 grid and Saha bitwise.
  Remaining differences are printed. Results: fid Xe bitwise 4248/6000, max rel 6.9e-7 (Rodas callback previously 6.1e-5); TM max 3.6e-8 (previously 5e-7); h50 Xe max 3.3e-7.

## Milestone 2 (this chunk)

- `evaluate_recfast_system_native` (native order) drives the Gear port.
- C-library `exp`/`log`/`pow` on the Float64 paths that native evaluates with libm: `hubble_table`/`cosmos_H` (`CosmosAccessors.jl`); `rf_NH`, `rf_SahaBoltz_HeIII/HeII`, `rf_Te_QS` (`RecfastPP.jl`, via `_npow`/`_nexp`; generic AD fallbacks unchanged).
- Native capture: `C36/native_parity/native_rhs_probe.cpp` (scratch: the p2 probe without the z < 50 filter) gives all-node, all-component RHS records for fid and h50. The probe's own history equals the fixture bitwise.
- Test: `test/native_parity_recfast_gear.jl` now asserts the whole fid preliminary history bitwise, including the derivative output, and all h50 states. The h50 dXe residual is `@test_broken`.

## Next chunks (in order)

0. Hubble spline (item 6b): port the GSL natural-cspline coefficient computation (`gsl_interp_cspline` / `solve_tridiag`) and `calc_spline_JC` evaluation in exact native arithmetic. Test H bitwise at all 6000 h50 probe z, then the h50 dXe residual.

1. **RHS operation-order port (item 4).**
   - Native capture: `C36/p2/p2_native_probe.cpp` already prints RFRL records (native f at native states, z < 50). Extend it to every node and all 4 components. Optional: a scratch-harness dump of the Gear step sequence (z, order, Dz_z, Newton loops), to localise any remaining divergence.
   - Test: f bitwise at all native states; then the preliminary history bitwise (fid, h50).
2. Tail (item 5) on `RFGearSolver`; test against the native output rows below z = 50 (`native_cosmorec_runmode0.txt` NAT_OUTPUT).
3. Make `preliminary_history`'s faithful default the Gear callback, once 1 is bitwise. SciML stays as an explicit diagnostic.
4. Multilevel solver `ODE_solver_Rec.cpp` (item 7), with a native capture of the pass-0 step sequence via the scratch harness (C36 scratch build), then the three-pass composition.
5. Broaden to the 10 cosmologies (`C36/run_native.sh` captures + `verify_reconstruction.py`); then the CAMB integration, the ℓ ≤ 9500 spectra and the matched benchmark (not before).

## Limitations

- No full replacement is claimed. The six original gates are unchanged, and their status will be re-measured only after items 7 and 9 are native-exact.
- Bitwise parity assumes the same glibc `libm` (`exp`, `pow`, `log`) as the native build. Any library-level difference will be reported per function.
