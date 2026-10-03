# Chunk 5b acceptance: one recombination pass

## Gates
- `test/chunk5b_ode_pass_native.jl`: **42/42** (focused).
- Full root `Pkg.test()` with `COSMOREC_NATIVE_DATA_DIR` set (tests registered through 5b): **49377/49377**, exit 0, 22m15.9s, `chunk5/pkgtest5b.log` (49335 + 42).
- Native comparison (one original pass, `runmode 1`, two byte-identical fresh runs): node grid bitwise; helium-switch node 1342 and `z = 1680.9103034346122` bitwise equal to the native `HeI_was_switched_off_at_z`; nodewise differences in native-tolerance units: X1s 0.114, rho 0.008, 2s 0.50, 2p 3.22 (monitored); Xe 1.7e-7, Te 8.0e-10 relative; unmonitored 3s/3p/3d <= 6.4e-7 relative.
- Benchmarks: full pass 758 ms (Rodas5P reltol 1e-12), 78 ms at 1e-8.
## Accepted scope
The ODE pass of the native `Xe_frac_effective_rates` without diffusion with the SciML solve supplied by the caller, the frozen-branch variant, observables per node. The Julia solution converges (tolerance study in the results) and agrees with the native within the native solver's own error.
## Limitations
Discrete switch node (no event-location derivative); excited levels are slaved/unmonitored (loose absolute tolerance by necessity); no native helium populations available for `runmode 1` (recovered through Xe); one cosmology; machine-specific data path (D1).
## Still not a full CosmoRec model
Recfast tail and assembly (5c/5d), gradients (5e), diffusion (phase 6+), runmode 0 parity.
