# Chunk 4a acceptance: explicit-input Saha initialization and ysol packing

## Independent gates
- `test/chunk4a_saha_native.jl`: **765/765** (supervisor rerun, `chunk4a/supervisor_native.log`); direct comparison with the original native `Set_Hydrogen_Levels_to_Saha`, `Set_HeI_Levels_to_Saha` and `copy_LI_to_ysol` outputs at 12 redshifts. Max relative errors: full X / packed y 5.06e-16, H LTE 4.16e-16, He LTE 4.41e-16 (tolerance 1e-14); packing from the native X reproduces the native y exactly; ground overrides, Xe slot and `rho = 1` are bitwise.
- `test/chunk4a_saha_ad.jl`: **280/280** (supervisor rerun, `chunk4a/supervisor_ad.log`); 256-bit ForwardDiff vs 256-bit finite differences 0.0, Float64 vs 256-bit FD <= 1.2e-15, directional <= 2e-28, prepared Mooncake VJPs (independent preparations, changed inputs) <= 2.2e-16.
- Full root `Pkg.test()`: **23701/23701**, exit 0, `Testing CosmoRec tests passed`, 17m12.9s, log `chunk4a/pkgtest.log` (22656 previous + 765 + 280).
- Benchmark numbers (BenchmarkTools, evals = 1): initialization 244 ns, packing 21 ns, both 250 ns, ForwardDiff Jacobian 1.22 us, Mooncake prepared hot gradient 2.49 us.

## Accepted scope
`src/SahaInit.jl` reproduces the native per-level Saha state initialization for the default 3-shell H / 2-shell HeI atoms and the native `copy_LI_to_ysol` packing (`flag_He >= 1`), with ALL cosmology-derived scalars, level data and constants supplied explicitly. Fixture `test/fixtures/native_saha_init.txt` (SHA-256 `d61d10b6...e7656`) is a mechanical extraction of a byte-reproducible original-library capture on the production CAMB-adapter route (provenance in the file header, `docs/CHUNK4A_RESULTS.md`, and the external `chunk4/native_capture/hashes.txt`). No native source or data was modified; no history file is bundled.

## Limitations to carry forward
- Clip (`Xe`, `Xp`) and override semantics are reproduced and unit-tested with synthetic inputs (no native oracle for clipped states); derivatives are zero in clipped regions and not claimed at clip edges; AD tests use interior points (z = 3000 base point only).
- One cosmology, 12 redshifts; no domain checks (native behaviour); level data fixed to the production atoms.
- The chunk3d `IN`-row cross-check of the design stage is not repeated in Julia.
- Recfast's `Xe_H` and `Xe_He` are ionized populations (`n_e,H/n_H`, `n_e,He/n_H`), not neutral X1s/XHeI1s fractions; 4a never consumes them directly.

## Still not a full CosmoRec model
No Cosmos accessors/splines (4b), no Recfast++ history (4c), no H(z)/cosmology, no sampled-HeI state switch, no SciML recombination ODE, no diffusion/feedback PDE or iteration loop, no history or CAMB spectrum parity.
