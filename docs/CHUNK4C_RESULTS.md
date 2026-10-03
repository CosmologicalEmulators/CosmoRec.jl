# Chunk 4c results: preliminary Recfast++ history (production terms) and sensitivities

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al. Use must be acknowledged; cite Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39), Rubino-Martin, Chluba & Sunyaev 2008. No native source or data was changed.

## Scope (and nothing more)
`src/RecfastPP.jl` (included from `src/CosmoRec.jl`, exports added) ports `Xe_frac` of Recfast++ for the PRODUCTION configuration read back from the native objects (fixture record RFV): fudge `F = 1.14`, `A2s1s = 8.2206`
(observed native value after startup although the adapter passes 0; not interpreted), no DM annihilation/decay, no magnetic heating, no reionization, no Chluba-Thomas correction function (`include_CF = 0`),
`eval_ion_Tr = 0`, varying-constant factors 1 (`aS = mS = 1`, `pS = 0`, `RM = 0`). Contents: the SI constants (`RecfastConstants`), `recfast_rhs` (native `evaluate_Recfast_System`, 4 equations), `recfast_rhs3!` (3-state SciML form),
`rf_Te_QS`, the Recfast++ `rf_SahaBoltz_HeIII/HeII`, `recfast_grid` (the three Saha segments with their data-dependent breaks, the linear-then-logarithmic ODE grid, 6000 nodes z = 25000 -> 0) and `recfast_history`
(Saha segments in closed form + ODE part through a caller-supplied stiff solve). The native Hubble function is the loaded Cosmos H(z) (4b `cosmos_H`, including the reproduced low-z native bug, which is irrelevant for z >= 100).
The package itself has NO ODE dependency (`[deps]` is empty by design): the SciML `Rodas5P` solves live in `test/chunk4c_helpers.jl` and `benchmark/chunk4c_benchmarks.jl`.
Recfast `Xe_H = n_e,H/n_H` and `Xe_He = n_e,He/n_H` are IONIZED populations, not neutral fractions. NOT implemented: the native Gear solver itself (step-by-step reproduction), DM/PMF/reionization/CF terms, the full CosmoRec ODE.

## Native fixtures
- `test/fixtures/native_recfast_history.txt` (4b; the 6000-node native history, SHA-256 `f893af2d...7716`).
- `test/fixtures/native_recfast_rhs.txt` (243929 B + header, SHA-256 `81be03e49f5270c76121bddd75846294d8739c592bd970e98cb21476de350ea6`): RFC/RFV/CFG constants and switches, RFR 900 direct calls of the ORIGINAL `evaluate_Recfast_System` (300 history nodes x {native state, two perturbed states}; the native `Hz`, `NH`, `TCMB` are printed with each result), RFS 36 `Te_QS`/`SahaBoltz_HeIII`/`SahaBoltz_HeII` probes, RFH 7 `H_z` probes. Verbatim records of `chunk4/native_capture/out4c/probe1.txt` (SHA-256 `c3caaffe60a9081836cf3d62302f02eaae67f9606d0b462599aa1132504af94f`) written by `harness4c.cpp` (SHA-256 `c3ee76a5150822d6f0980f9193d0b94e0151bf55987d5a4d59c5139ed411cea1`); two fresh-process outputs byte-identical; the native ODE is not re-run (the history file is the native `Xe_frac` output).

## Equation / grid audit (what was checked against the source)
Native `evaluate_Recfast_System` (evalode.Recfast.cpp:99-258), rate functions (:40-98), `NH` (cosmology.Recfast.cpp), `Te_QS` and the Saha helpers (Recfast++.cpp:56-95), `Xe_frac` segment loops and grid (Recfast++.cpp:119-335), constants (constants.Recfast.h). Reproduced exactly: the loop/break semantics (a Saha segment's last node is overwritten by the first node of the next, `j` is not incremented on a break), the thresholds `5e-4`, `1e-5`, `100/1933/100` node counts, `zc` scaled by `2.726/T0`, the 800-switch of the grid.

## Results (focused: native 8475/8475, AD 29/29)
- 900 direct native RHS probes: f3 (matter temperature) and f4 (cosmic time) bitwise (0.0); f1, f2 (He and H ionization) agree to 3.6e-16 / 3.8e-16 relative to the magnitude of their two cancelling terms. The RAW relative error reaches 5.2e-7 (f1) / 2.0e-6 (f2) at 28 native equilibrium nodes where the two terms cancel to ~1e-10 (conditioning, not a defect). Native `Hz`, `NH` reproduced to 1e-14; `TCMB` bitwise.
- `Te_QS`, `SahaBoltz_HeIII`, `SahaBoltz_HeII` at 36 probes: bitwise (0.0). H(z) probes (z >= 1) to 1e-13.
- Node grid: all 6000 nodes within 2.6e-15 relative; segment boundaries equal (ODE start node 2033). The grid is a DISCRETE operator: a 1e-12 relative change of the parameters moves the Saha-segment break by one node (2033 -> 2034), so no derivative of the grid is meaningful.
- Saha segments (2032 nodes): Xe_H, Xe_He, dXe, dXe_H bitwise; TM 2.8e-16; Xe on segment 2: 2.5e-11 raw, equal to the amplification of 1 ulp of the native libm pow/exp through `d + sqrt(d^2 + ...)` with |d| up to 4e5 (test gate: 8 eps max(|d|, 1)); segment 1 and the ionized plateau bitwise.
- ODE part (3967 nodes) vs the native solver (variable-order Gear, rel tol 1e-8, abs tol (1e-10, 1e-8, 1e-10)): Julia reference = Rodas5P reltol 1e-12 (self-consistency vs reltol 1e-10: 3.5e-9). Differences in NATIVE TOLERANCE UNITS `|Delta|/(1e-8 |y| + atol)`: Xe_He max 0.20 (worst z = 1363), Xe_H max 1.04 (z = 0, where xp = 1.7e-4 and the native absolute tolerance 1e-8 dominates), Xe max 1.04, TM max 48.9 (z = 7.4); medians 2e-4, 2e-7, 5e-8. Raw relative: TM 5.0e-7, Xe 6.1e-5, Xe_He 1.2e-2 (z -> 0, values ~5e-14 far below the native absolute tolerance). The native solver error is thus about its own tolerance; this is the measured native error, not a defect of the port. Gates: 1.0 / 2.0 / 100 / 2.0 units (Xe_He / Xe_H / TM / Xe).
- dXe, dXe_H: the native values are the solver's own derivative (`Sz.dy`), NOT the RHS at the node: on the stiff equilibrium branch `f(z, y_num)` differs from the native `dy` by up to 4e-4 (the RHS is ill-conditioned there). The port therefore compares the derivative of the Rodas5P dense output (`recfast_history` accepts a `(states, derivatives)` callback result): max |Delta|/max|ref| = 3.3e-2, median 2.8e-7 (gate 0.1). With a plain callback the RHS at the node is used (documented as ill-conditioned).
- AD-driven reformulations (equal to the native expressions to rounding; tested by the 900 RHS probes): `rf_inv_boltzmann` (`min(1e150, exp(+E/kT))` instead of `1/max(1e-300, exp(-E/kT))`; native clamp <= 1e300, effect only on a ratio that is insensitive to it) and `rf_inv_sahaboltz` (`1/S` without forming `S ~ 1e295`), which removed second-order Inf/Inf NaNs in the solver's nested Jacobian.

## AD through the history solve (focused AD suite, 29/29)
theta = [F, A2s1s, Yp, Omega_b, h100, T0]; outputs = (Xe, TM at Saha node 1000) + (Xe_He, Xe_H, TM at 9 ODE nodes, z = 2500 ... 5). Node grid and segment breaks are frozen primal values (discrete operator).
- ForwardDiff through `solve` (Rodas5P, reltol 1e-12, error control on primal values; the default norm produces NaN partials at an exactly zero error estimate) vs central finite differences (h = 1e-3, 1e-4, 1e-5): col-scaled error 1.6e-7 (h = 1e-3, best), up to 1.0e-5 at the other steps (noise limited by the solver tolerance); directional derivative (JVP) vs symmetric difference 3.3e-9. 256-bit comparison is possible only for the closed-form Saha node: 256-bit ForwardDiff vs 256-bit central difference 0.0; Float64 ForwardDiff vs 256-bit 1.9e-10 (conditioning of `d + sqrt(d^2 + ...)`). A 256-bit ODE solve was not attempted.
- Prepared Mooncake VJP (two independent preparations, changed parameters; `GaussAdjoint(autojacvec = MooncakeVJP())` with `Rodas5P(autodiff = AutoFiniteDiff())` and an EXPLICIT analytic Jacobian/time-gradient (ForwardDiff of the plain RHS), reltol 1e-10) vs the ForwardDiff projection `J'w`: random-weight projection 2.6e-10; Saha nodes 9.9e-17; TM block 1.5e-10; Xe_H block 1.3e-7; Xe_He block 4.7e-7 (tiny gradients, limited by the adjoint tolerance). Without the explicit Jacobian the finite-difference Rodas forward solve is unstable; QNDF/FBDF adjoints ran but gave wrong TM gradients (error 1.0): not accepted. A closure RHS (captured data) is rejected by Mooncake; the RHS must be a plain function.
- Limitations: derivatives of the Saha-segment breaks, of the node locations and of solver steps are NOT claimed; T0 enters only through physics at frozen nodes; the native Gear solver is not differentiated (it is only a reference).

## Benchmarks (`benchmark/chunk4c_benchmarks.jl`, `chunk4c/bench.log`, BenchmarkTools)
| path | median | memory |
|---|---|---|
| one RHS evaluation | 429 ns | 0 B |
| node grid (6000 nodes) | 316 us | 47 KiB |
| full 6000-node history, Rodas5P 1e-10, dense output | 18.4 ms | 2.9 MiB |
| full history, saveat nodes | 20.1 ms | 1.3 MiB |
| selected-node history (9 ODE nodes) | 13.8 ms | 204 KiB |
| ForwardDiff Jacobian (27 x 6) through the solve | 39.3 ms | 1.15 MiB |
| Mooncake prepare_gradient (2 samples, first includes compilation: min 22.7 ms, max 459 ms) | 241 ms | 17.4 MiB |
| Mooncake prepared hot gradient | 720 ms | 40.8 MiB |

## Files
`src/RecfastPP.jl`, `src/CosmoRec.jl`, `test/chunk4c_helpers.jl`, `test/chunk4c_recfast_native.jl`, `test/chunk4c_recfast_ad.jl`, `test/runtests.jl`, `test/fixtures/native_recfast_rhs.txt`, `benchmark/chunk4c_benchmarks.jl`, this file, `docs/ACCEPTANCE_CHUNK4C.md`. Failed approaches kept in `chunk4c/` logs and diag scripts. No Git operation.
