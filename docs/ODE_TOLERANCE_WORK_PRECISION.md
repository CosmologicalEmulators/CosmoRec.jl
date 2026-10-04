# ODE tolerance work–precision study (2026-10-04)

Read-only experiments on the validated optimized tree. **No reference/default,
physics, clamp, grid, algorithm, dependency or CI gate was changed.** All data
and scripts are text files outside git in `A/chunk21/`, where
`A = cmbcheb_test/local_analysis/cosmorec_differentiability_20260930`.

## Reference and interpretation

The caller-supplied reference uses Rodas5P, FullSpecialize and the validated
primal/state-Jacobian workspaces: `reltol=1e-12`, ground/rho absolute tolerance
`a1=1e-18`, excited absolute tolerance `aex=1e-14`. These are **local numerical
error controls**, not twelve-digit physical accuracy of every state.

CosmoRec's authors quote approximately 0.1% hydrogen and 0.1–0.2% helium output
accuracy. No explicit matter-temperature accuracy claim was found. See
`COSMOREC_PUBLISHED_ACCURACY_REVIEW.md`; faithful native-code parity and
physical-model accuracy must not be conflated.

## Screening and robustness

`screen_20261004T0033/` and `bracket_20261004T0040/` contain 19 profile runs,
including repeated tight controls. They record all 10000-node histories,
native/tight errors, switch nodes, actual solver work and five-sample
BenchmarkTools timing/allocation/thermal data. The tight screening output is
byte-identical to the previously accepted optimized benchmark output.

`multicosmo_20261004T0048/` and `multicosmo_boundsyes_20261004T0057/` reuse
native-provenance fiducial/ob/oc/h70/yhe inputs, without native regeneration.
The tight debug-flag histories reproduce the previously certified four varied
cosmology histories bitwise. Normal/debug differences of order 1e-6 in Xe are
reported separately, and profile errors always use a matched-flags reference.

| Profile | reltol | a1 | Fiducial screen median | RHS evaluations |
|---|---:|---:|---:|---:|
| tight | 1e-12 | 1e-18 | 1.80–1.82 s | 694672 |
| a16 | 1e-12 | 1e-16 | 1.45 s | 550376 |
| a15 | 1e-12 | 1e-15 | 1.09 s | 408392 |
| r3e12_a16 | 3e-12 | 1e-16 | 1.19 s | 437448 |

**These timings are screening observations, not matched speedup certification.**
The work reductions are hardware independent: about 21% fewer RHS calls for
a16 and 41% fewer for a15. Switch indices are unchanged in these runs.

The coarse relative-tolerance settings `1e-10`, `1e-9`, `1e-8` are faster but
exceed the unchanged global Xe<=1e-5 study budget. The seemingly good fiducial
`reltol=1e-11` profiles also fail on other cosmologies or compiler flags. All
rejected data remain recorded. Error is not monotone in these coupled solves.

### History results across both flags and all five cosmologies

| Profile | Worst Xe/native | Worst additional Tm vs tight, normal flags |
|---|---:|---:|
| a16 | 2.92e-6 | 1.96e-8 |
| a15 | 2.36e-6 | 1.86e-8 |
| r3e12_a16 | 3.52e-6 | 3.06e-8 |

All three preserve the Xe budget and the conservative additional-temperature
budget. **They do not all preserve every strict native-temperature verdict.**
The existing Tm/native<1e-7 criterion already misses at `ob` in the certified
tight debug run (1.358e-7). With a16, only that existing point misses
(1.332e-7, debug). With a15, `oc` additionally misses very marginally:
1.00518e-7 under debug flags and 1.00561e-7 under normal flags. With
r3e12_a16, additional temperature misses also occur. These are **reported
failures**, not hidden passes or modified CI tolerances.

The small PDE corrections are more sensitive: a15 changes one correction by
up to ~0.8% relative to that correction's peak, while final history and CMB
errors remain much smaller. Existing strict internal-parity tests stay intact.

## Parameter-gradient comparison

`gradients_20261004T0102/` compares the existing four-parameter composed route
`[F,A2s1s,hscale,nbscale]`, with the tight switch node frozen, on nine Xe and
nine Tm outputs. Direct ForwardDiff-to-ForwardDiff comparison is used, not
adaptive-ODE finite differences. All Jacobians are finite.

Maximum dimensionless elasticity-Jacobian difference versus tight:
a16 **3.08e-6**, a15 **2.94e-6**, r3e12_a16 **3.81e-6**. The meaningful
hscale/nbscale columns differ by a few ppm. F/A2s1s columns are effectively
zero on this route (max absolute reference values ~1.5e-18 and 5.2e-21);
their large relative percentages are meaningless and their scaled effects
are negligible. This is not certification of derivatives through complete
cosmological initialization or of reverse mode through Rodas5P.

## CMB observable validation

`cmb_20261004T0107/` uses the existing isolated CAMB 2.0.4 injection runner,
CosmoRec recombination and Mead2020 lensing/nonlinear settings, with original
native outputs reused and all cosmology/configuration metadata checked.
Fresh tight and candidate spectra were generated for all five cosmologies.
All eight stored columns are tested over **ell=2:3000**: unlensed TT/EE/TE,
lensed TT/EE/BB/TE and phiphi. TE uses correlation-normalized error, not a
divergent fractional error at its zeros. CVL means full sky, unit multipole
bins, no noise or beam suppression; diagonal sums are not a joint likelihood.

| Profile | Max incremental spectrum error vs tight | Max CVL-normalized residual | Max spectrum error vs native |
|---|---:|---:|---:|
| a16 | 1.17e-6 | 6.40e-5 | 9.68e-7 |
| a15 | 1.26e-6 | 6.87e-5 | 7.93e-7 |
| r3e12_a16 | 1.93e-6 | 1.05e-4 | 1.29e-6 |

All meet the prospective study budgets (incremental spectrum error<=1e-5,
max CVL residual<=1e-3). This does not certify the entire cosmological domain
or multipoles beyond 3000. Improved agreement with native in one metric does
not establish improved physical-model accuracy.

## Matched runtime confirmation

The bounded attempt `paired_a15_20261004T0116/` reused the validated guarded
benchmark through numeric tolerance substitution only, with original and
configured driver hashes, a fixed P-core, compilation/warm-up before an
in-process cool gate, and unchanged acceptance thresholds.

**Both timed batches were BLOCKED**, despite the runner exiting zero. Another
Julia computation was active and the package remained near 100 C.
No wall-time speedup is certified from that attempt; no other jobs or host
settings were changed. Screening times suggest ~20% less elapsed for a16 and
~40% for a15, but these remain indicative. Fewer solver evaluations are proven.

## Recommendation, without changing the validated reference

- **Conservative candidate:** `reltol=1e-12, a1=1e-16, aex=1e-14`.
  It keeps all previously passing strict native-temperature points passing
  across the five tested cosmologies/flags, with fewer RHS evaluations and
  ppm-level CMB changes. The existing `ob` temperature miss remains explicit.
- **More aggressive, conditional candidate:** `a1=1e-15`, same reltol/aex.
  Its CMB and incremental history errors are similarly small and it reduces
  work more, but it introduces the marginal `oc` native-temperature miss.
  Using it requires an explicit decision that observable-level accuracy is
  the acceptance target rather than every old internal/native-temperature gate.

The tight benchmark/reference callback and all helper defaults and CI tests are unchanged. No fast setting has
been adopted. Confirm matched timing on a quiet host before quoting a precise
speed gain; extend the domain/ell validation before broader production use.
