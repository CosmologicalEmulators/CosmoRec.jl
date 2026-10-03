# Chunk 4 results (phase 4 complete): initialization slice

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al.; cite Chluba & Thomas 2010, Chluba & Sunyaev 2006, Rubino-Martin et al. 2008. No native source, library or table was changed; all new native outputs come from external read-only harnesses (`.../chunk4/native_capture/`, hashes in `hashes.txt` and in the fixture headers); the Hubble table is a CAMB-generated INPUT. Nothing is committed.
Design: `docs/CHUNK4_DESIGN.md`. Recfast `Xe_H = n_e,H/n_H` and `Xe_He = n_e,He/n_H` are IONIZED populations, not neutral X1s/XHeI1s fractions (corrected labels preserved).

| sub-chunk | source | fixtures | focused (native / AD) | full root `Pkg.test()` | docs |
|---|---|---|---|---|---|
| 4a Saha per-level initialization + `copy_LI_to_ysol` packing | `src/SahaInit.jl` | `native_saha_init.txt` | 765 / 280 | 23701/23701 | `CHUNK4A_RESULTS.md`, `ACCEPTANCE_CHUNK4A.md` |
| 4b Cosmos accessors, GSL natural cubic splines, loaded Hubble | `src/CosmosAccessors.jl` | `native_recfast_history.txt` (749 kB), `native_hubble_input.txt`, `native_cosmos_accessors.txt` | 14097 / 16 | 37814/37814 | `CHUNK4B_RESULTS.md`, `ACCEPTANCE_CHUNK4B.md` |
| 4c preliminary Recfast++ history (production terms) | `src/RecfastPP.jl` | `native_recfast_rhs.txt` (+4b history) | 8475 / 29 | 46318/46318 | `CHUNK4C_RESULTS.md`, `ACCEPTANCE_CHUNK4C.md` |
| 4d sampled-HeI state switch | `src/HeliumSwitch.jl` (+ packing branch) | `native_helium_switch.txt` | 984 / 32 | 47334/47334 | `CHUNK4D_RESULTS.md`, `ACCEPTANCE_CHUNK4D.md` |
Each full run was a fresh `Pkg.test()` (logs `chunk4a/`, `chunk4b/`, `chunk4c/`, `chunk4d/pkgtest.log`); every count is for that stage only. The earlier `Pkg.test()` started accidentally with discarded output was killed and is not counted.

## Worst errors against native outputs
4a: 5.1e-16 (X, packed y), LTE 4.4e-16. 4b: 1.6e-14 for every accessor except `NHeIII` 6.0e-7 (native cancellation); Saha-Boltzmann bitwise. 4c: direct native RHS f3, f4 bitwise and f1, f2 3.8e-16 relative to the cancelling terms (raw 2.0e-6 at equilibrium nodes); node grid 2.6e-15; Saha segments bitwise (Xe on segment 2: 2.5e-11, ulp amplification); ODE nodes within 0.20 / 1.04 / 48.9 native-tolerance units (Xe_He / Xe_H / TM), i.e. the native Gear solver's own error; derivatives 3.3e-2 of the peak (median 3e-7). 4d: bitwise (80 states).
Gradients: 4a FD 1.2e-15, VJP 1.2e-14 (relative to the 256-bit central differences / ForwardDiff); 4b z-derivatives 1.3e-11 worst, 24001-input directional 3.8e-15, VJP 4.2e-16; 4c ForwardDiff vs finite differences 1.6e-7 (best step), JVP 3.3e-9, Mooncake VJP 2.6e-10 (random projection) to 4.7e-7 (tiny Xe_He block); 4d 0.0.

## Limitations carried into phase 5 (all documented in the stage results)
- 4a: clip and override branches tested only with synthetic inputs; derivatives not claimed at clip edges.
- 4b: piecewise definitions (z_saha = 3500, zsRe, end nudge, loaded-H edges); low-z native Hubble bug (`lHz[0] = 0`) reproduced, not repaired; `Omega_L`, `Nb0`, `rho_g_gr` are inputs; one cosmology.
- 4c: discrete Saha-segment grid (a 1e-12 parameter change moved the ODE start node by one); native `dXe` is the solver derivative (compared via dense output); no step-by-step native solver reproduction; no 256-bit ODE solve; the reverse route needs a plain-function RHS with analytic Jacobian/time gradient and error control on primal values (QNDF/FBDF adjoints rejected); two AD-friendly reformulations (`rf_inv_boltzmann` with clamp 1e150 instead of 1e300, `rf_inv_sahaboltz`) equal native to rounding; the package has no ODE dependency (the SciML solves live in tests/benchmarks).
- 4d: frozen branch, no event-time sensitivity, production atoms only; the inline native driver statements were executed by the harness (not callable).
- Not included: DM/PMF/reionization/CF Recfast terms, diffusion/feedback, the full CosmoRec ODE, history or CAMB spectrum parity.

## Remaining (phase 5-10, not started)
Phase 5: the full recombination ODE assembling 4a-4d with `fcn_effective` (Chunk 3d) and the Hubble/background closure; then diffusion/feedback and iteration loops, derivative validation of the whole history, and the CAMB-level comparison with `native_camb_cosmorec_thermo.txt`. No phase-5 work was done.
