# Chunk 3c acceptance: default helium-photon absorption channels

## Independent component-level verification

The supervisor independently reran:

- `test/chunk3c_hiabs_native.jl`: **8,870 assertions passed**. This checks direct original-library absorber/Pesc/correction outputs, both approximation flags, both lines, all stored native query/window samples, domain behavior, accumulation slots, and independent sign/conservation invariants.
- `test/chunk3c_hiabs_ad.jl`: **462 assertions passed**. ForwardDiff Jacobians and 256-bit finite-difference references, one-sided switch/seam checks, prepared/reused Mooncake VJPs against ForwardDiff Jᵀw.

Original native source is read-only/unchanged; no native tables are modified. The focused package report and test artifacts are in `local_analysis/cosmorec_differentiability_20260930/chunk3c/`.

Worker also reports full `Pkg.test()`: **18,094/18,094**, 16m50.6s, `chunk3c/pkgtest.log`. The count includes previous suites. This has not been independently rerun by the supervisor, but the full log and exit status are retained. It ran under the human's explicit thermal-throttling override.

## What is actually validated

Only the table-bounded *standalone native helium-photon-to-neutral-hydrogen absorption function and its Julia port*: singlet Ly-alpha and intercombination outputs, `DP_interpol_S/T`, corrected escape probabilities, and their mapped increment signs/indices for valid table queries. The compact original fixture is `test/fixtures/native_hiabs_components.txt` (17-digit text, original sources/library/DP+f.corr hashes and attribution). Focused native and AD comparisons pass at their documented tolerances.

The native `fcn_HeI_effective` wrapper owns a file-local parameter object, so it could not be called directly from the stand-alone harness. **Fixture K rows for its combined state updates are transcribed branch/sign logic from source lines284–321 applied to actual native channel outputs**, not direct `fcn_HeI_effective` recordings. Channel quantities themselves are direct native function calls; test sign/conservation invariants validate the separate mapping.

## Important production limits / incomplete default physics

- Explicit-integral native `DPesc_coh` fallback is not ported. Julia throws a typed domain exception off the captured DP table. Cosmology/parameter-shifted trajectories outside the tabulated domain are not supported by this component yet.
- Only `neff=30` / `fac_50` DP data are represented. Some table’s native OOB reads are replaced with safe errors. `f_t`/`f_b != 1` not validated.
- Table-switch, T-sheet, eta/tau-cell seams, and opacity underflow are piecewise or nonsmooth; one-sided behavior is tested. Do not differentiate across the jumps as if smooth.
- The accepted HeI base chunk3b and H base chunk3a were tested separately. No **single total native `fcn_effective` state derivative fixture** has verified their composition. No coupled ODE solve, no source-driven fixture for the whole H/He RHS, no integration Jacobian/adjoint at CosmoRec system scale.
- Still unimplemented: the fallback integral, helium diffusion, helium feedback, Saha/Recfast initialization, sampled helium-state removal, production background and output interpolation, low-z matching and complete finite ODE→PDE→ODE→PDE→ODE chain.

Thus Chunk3c is accepted for its **isolated, bounded absorption-function gate only**. The complete first-pass H+He population RHS and CosmoRec.jl full-fidelity solver remain future work; do not label this component as an end-to-end default-model result.
