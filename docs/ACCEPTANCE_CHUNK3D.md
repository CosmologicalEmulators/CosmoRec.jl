# Chunk 3d acceptance: native-fixtured pointwise `fcn_effective`

## Independent numerical/AD verification

The supervisor independently ran the native comparison and AD suites through the persistent local test environment:

- `test/chunk3d_fcn_native.jl`: **566/566** assertions. It checks the 15-slot native `Data_Level_I.g` layout against calls to the original static-library `fcn_effective` for 12 explicit point states, with the default first-pass non-diffusion flags; maximum relative match **8.35e-15** in the supplied DP-table domain. An absorber-off fixture is also captured for explicitly out-of-DP-domain states.
- `test/chunk3d_fcn_ad.jl`: **2203/2203** assertions. It checks composed ForwardDiff/BigFloat finite-difference and directional derivatives plus prepared Mooncake VJPs. At low `Float64` roundoff through near-canceling absorber terms, FD error is ~1.67e-8; this is distinguished from the ~2.91e-15 reverse-VJP agreement.
- Worker reports the full unconditional suite **20,863/20,863**, exit0, in `chunk3d/pkgtest.log`. The supervisor did not run another 20-minute full suite at this boundary; both newly changed focused suites have been independently verified.

Native `fcn_effective(double,Data_Level_I&)` is actually exported by the original library, called directly by an external C++ harness (not a manual transcription), initialized through original startup and atomic/table functions, on 12 state vectors with controlled source-native Saha plus explicit perturbations. The harness confirms absorber switch and native caches/repeated calls, and both absorptive and nonabsorptive output vectors. Fixture `test/fixtures/native_fcn_effective.txt`; original source/library hashes, parameter flags, state order, units and repeatability appear in `docs/CHUNK3D_RESULTS.md`. Native CosmoRec/CAMB checkout sources were not changed.

## Correct accepted scope

`src/FcnEffective.jl` is a Julia composition of the already separately validated H, HeI-base and absorption blocks, plus the thermal-rate assignment and native ground-state fraction reconstruction. Inputs/background/tables are explicit. This is a **pointwise RHS mapping** in the cosmic-time derivative convention `dX/dt`, on finite supplied state vectors. It is NOT an initial-state generator, an ODE solve or a complete inference result.

## Hard limitation: incomplete native absorber domain

Some valid native reference states (`z=1500,1200,800` in this fixture) require the original explicit `DPesc_coh` integral because a DP/Pesc query leaves its stored table. That integral has not been ported. The native fixture stores the true absorber-on output, but the Julia test correctly expects `DPTableDomainError` on those rows and compares the Julia/base H+He output to the native absorber-OFF output. Thus the code is **not yet full default `fcn_effective` over the native physical domain**. Do not describe the coupled RHS as completely native-equivalent until either DPesc_coh is implemented with fixtures/AD gates or the physically admissible model/supported parameter domain is established not to hit this region. Other unimplemented physics are listed in `docs/CHUNK3D_RESULTS.md`: diffusion/feedback, initialization, helium switch, background/ODE/PDE, low-z match.

## Stop point

Chunk3d is accepted for the tested native pointwise composed RHS within its explicit DP-table support. Chunk3d does not establish that its absorber remains supported along a full recombination history. The next code stage cannot be called a complete ODE until it has source-backed Recfast/Saha initialization and helium-switch behavior plus a tested policy or fixture for explicit-DPesc fallback. No commits/staging/pushes.
