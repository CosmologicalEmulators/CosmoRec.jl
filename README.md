# CosmoRec.jl

A Julia translation of CosmoRec for computing cosmological recombination histories.

## Solver configurations

The library is solver-agnostic: callers supply the ODE solve callback. The
`solve_phase5` helper in `test/chunk5_helpers.jl` uses the adopted configuration
below by default; it is not an exported package API.

| Setting | Adopted helper default | Experimental, opt-in |
|---|---|---|
| Solver | `Rodas4P()` | `Rodas4P()` |
| Relative tolerance | `1e-10` | `1e-5` |
| Ground-state absolute tolerance (`a1`) | `1e-18` | `1e-18` |
| Excited-state absolute tolerance (`aex`) | `1e-12` | `1e-2 * 1e-5` |
| Forced step endpoints | `tstops = znodes` | `tstops = znodes` |
| ODE block size | 50 | 50 |

With the helper, select the experimental configuration explicitly with
`solve_phase5(args...; reltol = 1e-5, a1 = 1e-18, aex = 1e-2 * 1e-5, tstops = true)`.
It does **not** change the default. Fixture-based regression tests retain an
explicit tight reference configuration: `Rodas5P()`, `reltol = 1e-12`,
`a1 = 1e-18`, `aex = 1e-14`, and no forced step endpoints.

## Measured forward performance

Full runmode-0 histories (three ODE passes and two diffusion PDE stages), measured
on an Intel Core i7-13700H with Julia 1.12.6, one pinned CPU thread, and 10,000
output redshifts. Times are hot medians: compilation and warm-up are excluded.

| Implementation / configuration | Median time | Speedup over Julia tight reference |
|---|---:|---:|
| Julia tight reference | 1.689 s | 1.00× |
| **Julia adopted, `reltol = 1e-10`** | **0.856 s** | **1.97×** |
| **Julia experimental, `reltol = 1e-5`** | **0.330 s** | **5.12×** |
| Original C++ CosmoRec, hot public API | 0.377 s | — |
| Original C++ CosmoRec, core plus output interpolation | 0.357 s | — |

The experimental configuration is about 2.60× faster than the adopted one here.
The native public API includes per-call setup excluded from the Julia timings;
the native core measurement is closer, but the timing boundaries are not exactly
identical. Algorithms and tolerances also differ. Native outputs match the
reference fixtures bitwise. Measurements used five Julia samples and ten native
samples per bracket, cooling after warm-up; no throttle-counter growth was
recorded during measured samples. These are hardware-specific results, not
universal speed guarantees.

### Accuracy and experimental limitations

At the benchmark cosmology, the maximum fractional changes in **lensed** CAMB TT
and EE relative to native-CosmoRec-history CAMB over multipoles 2–10,000 were:

| Configuration | TT maximum error | EE maximum error |
|---|---:|---:|
| Adopted, `reltol = 1e-10` | 0.0000114% | 0.0000235% |
| Experimental, `reltol = 1e-5` | 0.000103% | 0.000114% |

Both meet the study's 0.05% budget on TT and EE. The experimental configuration
also passes the established ForwardDiff smoke comparison against the tight
reference (elasticity error 1.87e-6, threshold 1e-4). However, that check covers a
nine-node, frozen-branch map at **one parameter point**, not dense or domain-wide
gradient validation. Its internal diffusion-feedback diagnostic differs from
the native reference by about 0.8%, above the earlier internal accuracy gates.
Broader validation is required before adopting this experimental setting.
