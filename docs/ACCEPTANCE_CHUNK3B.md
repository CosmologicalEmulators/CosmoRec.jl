# Chunk 3b acceptance: neutral-helium base RHS components

Chunk3b is accepted only for the scoped HeI population-rate components, not a complete default helium or H+He recombination RHS.

## Independent validation

Supervisor focused runs reproduced `test/chunk3b_rhs_native.jl`: **1,783/1,783 native assertions plus 95/95 independent physical invariants**, and `test/chunk3b_rhs_ad.jl`: **320/320 ForwardDiff/Mooncake assertions**. The combined native He components match the compact external C++ fixture. Independent compiled fixtures/logs are in the surrounding analysis `chunk3b/native_capture/`; native CosmoRec source remains clean. The fixture is text-only with source/table hashes, component state mapping, units and CosmoRec attribution.

The full unconditional package suite (chunks1a–1d,2,3a,3b) was then independently run by supervisor exec684: `julia --project=<repo> -e 'using Pkg; Pkg.test()'` -> **8,762/8,762**, exit0, 14m49.3s. The test log is external: `cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk3b_resume01/supervisor_pkgtest.log`.

## Scope / explicit exclusions

This chunk ports only the `fcn_HeI_effective` base continuum/photoionization, interlevel, two-photon and native HeI Lyman channels, with the production two-shell rate data and atom-state ordering. It includes focused ForwardDiff and Mooncake checks at query/interpolation map and RHS level; it does not solve an ODE.

**The active default H-I absorption of HeI photons remains unimplemented and uncaptured.** Production default `HI_absorption=2` is remapped to `1`; below `zcrit_HI=3400` this invokes H-I line absorption using HeI populations, DP/Pesc interpolation data, and Ly-alpha plus intercombination channels. Chunk3a's H-only fixture is not an oracle for this path. Chunk3b does not complete or validate the default coupled H+He RHS. Chunk3c must first capture/implement/test that path and its application to H/He/electron state derivatives.

Other pending components include helium diffusion/feedback, coupled initialization, redshift/background, sampled helium-switch state projection, and both ODE/PDE iteration paths. Do not write documentation describing the whole default RHS or a full CosmoRec history as implemented. No commits or pushes; user code is still uncommitted on develop.
