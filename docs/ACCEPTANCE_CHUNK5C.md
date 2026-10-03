# Chunk 5c acceptance: Recfast tail

## Gates
- `test/chunk5c_recfast_tail_native.jl`: **17/17** (focused). Full root `Pkg.test()` (`COSMOREC_NATIVE_DATA_DIR` set; tests registered through 5c): **49394/49394**, exit 0, 21m52.7s, `chunk5/pkgtest5c.log` (49377 + 17).
- Native comparison: tail grid bitwise; tail from the native inputs vs the native rows: Xe 1.03 and Te 166 native-tolerance units (relative 5.9e-5 and 1.8e-6; the native Recfast solver's own error); end-to-end from the Julia pass: Xe 1.4e-7, Te 3.0e-8 relative; rescaling factor reproduces the native derivative to 1e-12.
- Benchmark: `recfast_tail` 2.76 ms.
## Accepted scope
`compute_Recfast_part` / `Xe_frac_rescaled` with the 4c RHS, the loaded H(z) including the reproduced native low-z Hubble bug, caller-supplied stiff solve.
## Limitations
`dXei` input conditioning-limited (2.1e-4); the native tail error (up to 1.8e-6 on Te) dominates every comparison; one cosmology; machine-specific data path (D1).
## Still not a full CosmoRec model
Assembly/grid output (5d), gradients (5e), diffusion (phase 6+).
