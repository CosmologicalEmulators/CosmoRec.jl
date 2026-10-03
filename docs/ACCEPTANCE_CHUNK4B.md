# Chunk 4b acceptance: explicit Cosmos accessors and splines

## Independent gates
- `test/chunk4b_cosmos_native.jl`: **14097/14097** (supervisor rerun): 658 native redshifts x 21 accessors, 6 Saha-Boltzmann temperatures, history/Hubble/accessor fixture hashes, spline node reproduction, branch and domain tests, 4a fed by 4b. Max relative errors: <= 1.6e-14 for every accessor except `NHeIII` (6.0e-7, native cancellation `fHe - XHeI1s - XHeII1s` at z >= 3500, same operations as native).
- `test/chunk4b_cosmos_ad.jl`: **16/16** (supervisor rerun): ForwardDiff vs 256-bit central differences (<= 1.3e-11), directional derivative over 24001 inputs 3.8e-15, prepared Mooncake VJP (independent preparations, changed inputs) 4.2e-16.
- Full root `Pkg.test()`: **37814/37814**, `Testing CosmoRec tests passed`, log `chunk4b/pkgtest.log` (23701 previous + 14097 + 16).
- Benchmarks (BenchmarkTools): spline set 592 us, `cosmos_Xe_Seager` 43 ns, `cosmos_H` 61 ns, Mooncake prepared hot gradient 10.2 ms (24001 inputs).

## Accepted scope
GSL-natural-cubic-spline port, the native Cosmos accessors and Saha-Boltzmann branches, the loaded-Hubble spline including the reproduced native `lHz[0]` bug, all from EXPLICIT history nodes, constants and Hubble table (text fixtures with tested hashes; the 749 kB history is plain text, no binaries). Recfast `Xe_H`/`Xe_He` are ionized populations, not neutral fractions. Out-of-range queries throw (native: GSL abort).

## Limitations to carry forward
Derivatives are not claimed across the 3500 and zsRe switches, spline ends or loaded-H edges; the low-z H bug is native behaviour; `Omega_L`, `Nb0`, `rho_g_gr` are inputs; one cosmology.

## Still not a full CosmoRec model
No Recfast++ integration (4c), no sampled-HeI switch (4d), no SciML recombination ODE, no diffusion/feedback PDE or iteration loop, no history or CAMB spectrum parity.
