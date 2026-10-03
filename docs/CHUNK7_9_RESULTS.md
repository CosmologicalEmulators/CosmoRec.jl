# Chunks 7a-7c, 8a, 9a results: HI radiation PDE, correction integrals, feedback (stage-level; the runmode-0 composition is NOT accepted, see docs/STATUS_PHASES7_10.md)

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al. Use must be acknowledged; cite Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39). No native source or data was changed. No Git operation was performed.

## Scope
The production configuration of `compute_DPesc_with_diffusion_equation_effective` (nShells 3, nS_2gamma 3, nS_Raman 2, zs 2500, ze 500, dz_out 10, theta 0.999 for z > zs - 100 and 0.55 otherwise, Diff_iteration_max 2, induced_flag 2), and the feedback of `Modules/Diffusion_correction.cpp`. Other configurations throw (`hi_pde_setup`, the matrix elements beyond n = 3).

## Native oracles (external, read-only harnesses; two fresh-process outputs byte-identical for each)
| harness | content | fixture |
|---|---|---|
| harness7a | PDE setup: grid, ratios, atomic data; def_PDE coefficients at 9 z | native_hi_pde_setup.txt, native_hi_pde_def.txt |
| harness7c | ORIGINAL PDE run, plus a stepping-loop replica (bitwise equal) with per-step spectra | native_hi_pde_march.txt, native_hi_pde_integrals.txt |
| harness8a | the native pd/Dnem spline knots (read from the GSL objects) | native_hi_pde_coefficient_nodes.txt |
| harness9a | ORIGINAL interpolate_DF and fcn_effective with the correction on/off | native_hi_diffusion_feedback.txt |

Rejected evidence: the first harness9a capture switched DI1_2s_correction_on on after setup_DF_interpol_data, so the DI1 spline was never built (`chunk9/native_capture/out/failed_run*`).

## Results (focused; logs in chunk7, chunk8, chunk9)
- **7a, 99/99** (`focused7a_run2`):
  - grid bitwise; A_SH / Gamma_np / A_npks <= 1.7e-14 (libm log10/pow ulps in the Storey-Hummer recursion); profiles <= 5.1e-16; ratios <= 6e-16.
  - The native y = 0 and y = 1 probes are 0/0 (NaN) in the ORIGINAL code; NaN patterns are compared.
  - The first run failed on a wrong expectation (nS_2gamma = 4 is capped by nShells = 3 natively); it is preserved.
- **7b, 115/115:**
  - background, pd and aV bitwise; A, B, C <= 1e-15 of each array's maximum;
  - D 4.8e-14 and Dnem_eff 2.9e-13, from the native cancelling Dp/Dn correction.
- **7c, 48/48** (`focused7c_run3`):
  - final spectrum 2.8e-15 of its maximum; polint_JC lower boundary bitwise.
  - The early-step raw core differences of 2.7e-9 come from the Dnem cancellation. They lie under the FIXED bound 4 eps kappa_max = 1.8e-7, where kappa_max = 2e8 is asserted on the fixed native input history (measured 1.88e8 at z = 2499.5).
  - The projection on the Dnem response directions (residual 1.1e-13) is a diagnostic only.
  - The first run failed on placeholder gates; it is preserved.
- **8a** (`focused8a_run2`):
  - Single outputs from the native spectra and the full Julia stage both pass the FIXED bound of 1e-4 of the residual scale. The derivation and its caveats (an estimate, not a bound; 59-90 of 199 outputs contain unconverged sub-integrals) are in STATUS_PHASES7_10.md. Raw and absolute errors are recorded alongside.
  - Diagnostics: with the native knots injected the 2g errors are unchanged (2.4e-5). Julia's own ulp-level input response is up to 9e-5. The full-pipeline differences therefore sit at the conditioning level of these near-equilibrium outputs.
- **9a, 161/161:** interpolation 1.3e-16; F with the correction on/off 3.5e-13; the correction term at 0.69 rounding units of the native F_on - F_off difference.
- **Local AD on dynamic inputs** (`test/chunk7_10_pde_ad.jl`, 20/20; populations -> 7b, 7c, 8a with a 3-step march, and DF values -> 9a RHS):
  - ForwardDiff directional vs 256-bit central differences: <= 1.5e-10.
  - Prepared Mooncake VJPs vs ForwardDiff: <= 2.2e-15.
  - This is LOCAL stage AD only.

## Not accepted (see STATUS_PHASES7_10.md)
The runmode-0 composition (10a parity, 10b composed derivative, a composed reverse-mode gradient), the public-path benchmarks of 7-10, and a full Pkg.test including 7-10.

## Benchmarks (public path; `benchmark/chunk7_10_benchmarks.jl`, `chunk10/bench7_10.log`; BenchmarkTools, evals = 1)
Measured while one single-threaded diagnostic probe was running on another of the 20 cores.
| path | median | memory |
|---|---|---|
| 7a two-photon/Raman tables (load + splines) | 1.91 ms | 1.46 MiB |
| 7a PDE setup (grid + ratios, 2353 points) | 0.88 ms | 2.16 MiB |
| 7b def_PDE coefficients at one z | 172 us | 0 B |
| 7c Lagrange O2 weights | 255 us | 184 KiB |
| 7c one Step_PDE_O2t | 210 us | 0 B |
| 7c full march 2500 -> 500 (200 steps) | 43.0 ms | 1.31 MiB |
| 8a integrals at one output | 229 us | 564 KiB |
| 8a full PDE stage (march + 199 x 4 integrals) | 73.2 ms | 112 MiB |
| 8a ForwardDiff directional derivative of the stage (24000 inputs, 1 partial) | 154 ms | 194 MiB |
| 9a hi_diffusion_feedback | 12.9 us | 75 KiB |
| 9a recombination_rhs with feedback (7-state) | 875 ns | 1.1 KiB |
| 9a one ODE pass with feedback (Rodas5P, tight tolerances) | 765 ms | 882 MiB |
| 10a one diffusion stage from stored rows | 85.8 ms | 118 MiB |
| 10a full runmode-0 iteration (3 passes, 2 stages, tail, 10000-node assembly) | 2.53 s | 2.86 GiB |

Prepared Mooncake gradient through one full PDE stage (`chunk10/probe_mc_stage.log`, single runs, not BenchmarkTools): prepare 169 s (cold, includes compilation), first gradient 2.25 s, hot gradient 1.71 s, 12.1 GiB peak RSS.
