# Native CosmoRec text-fixture contract (specification only)

**Status: this is a specification for future physics chunks. No native
CosmoRec fixtures are generated, required, or authorized by Chunk 1a.**
Chunk 1a's own numerical reference is the closed-form analytic solution of a
toy ODE (`test/chunk1a_stiff_ad_probe.jl`), not a native-code fixture.

This document fixes, in advance, what a later physics chunk must capture from
a native CosmoRec run (built from the pinned source revisions in
`VERIFIED_ANALYSIS.md`) before any Julia translation of that chunk's physics
begins, so that native vs. Julia comparisons are apples-to-apples and
reproducible. It does not authorize running the native code, modifying it, or
copying its atomic-data tables (table redistribution/licensing is an open
prerequisite, per the Chunk 0 design and its review, §"License, native
quirks and type sketches").

## 1. Provenance that MUST be recorded with every fixture

- CosmoRec git revision (expect `086769055f61ae0c244a53dd381ee65b624d0ac3`
  per `VERIFIED_ANALYSIS.md`, or the then-current pinned revision) and
  `CAMB-cosmorec` revision (expect `fa3f097343fbbe427cc04b4f5f0041c22c6ec764`).
- Exact build flags/Makefile variant used (compiler, optimization level,
  OpenMP on/off) -- OpenMP reductions can change floating-point rounding
  order and hence bit-for-bit (though not scientifically) reproducibility.
- The exact `runmode`, `accuracy` (`runpars[1]`), and every `runpars[]` value
  actually passed, matching `VERIFIED_ANALYSIS.md` ("Default calculation and
  active modes"): `runmode=0`, `accuracy=0`, 3 H shells, 500-shell effective
  rates, 2 He shells, 2-photon through n=3, Raman through n=2, induced mode 2,
  `Diff_iteration_min=0`, `Diff_iteration_max=2`, HI-absorption flag 2
  (remapped internally to flag 1 + approximate-treatment flag).
- Cosmological input values and units exactly as CAMB's adapter passes them
  (`Omega_b h^2`-derived `Omega_b`, `Omega_c`, `Omega_k`, `N_eff`, `H0` in
  km/s/Mpc, `T_cmb` in K, `Y_p` dimensionless mass fraction) -- see
  `B/fortran/cosmorec.f90`.
- The supplied `H(z)` array exactly as built by the CAMB adapter: **10,000
  uniformly spaced redshifts from 10,000 down to 0** (not log-spaced -- this
  corrects an error in an earlier draft report; see `VERIFIED_ANALYSIS.md`,
  "Default calculation and active modes") together with the corresponding
  `H(z)` values and units (1/s or equivalent, as constructed in
  `cosmorec.f90`).
- Whether the known native `Hubble` endpoint zero-fill
  (`C/Development/Cosmology/Cosmos.cpp:104-145`, documented as untested-impact
  in `VERIFIED_ANALYSIS.md`, "Untested native endpoint discrepancy") was left
  as-is (baseline fixture) or patched (corrected fixture) for this capture --
  **both must be captured as separately labeled fixture sets**; the baseline
  (unpatched) capture must exist first, per the Chunk 0 review ("Freeze the
  existing native behavior before any endpoint correction").

## 2. Every intermediate/output array to capture (not just final Xe/Te)

For the eventual full three-ODE/two-PDE composition
(`VERIFIED_ANALYSIS.md`, "Default calculation and active modes"), each of
the following is a **separate** named fixture, because each is a separate
future chunk's acceptance target (per the Chunk 0 review's corrected
dependency order):

1. Saha/Recfast preliminary initialization: `X_e`, `X_{H,1s}`, `X_{He,1s}`,
   `T_m` (or `rho=T_m/T_gamma`) at `z_start` and on a diagnostic grid
   approaching it from above.
2. Effective-rate table values: `A_i`, `B_i`, `R_{ij}` at a grid of
   `(T_gamma, T_m/T_gamma)` probe points spanning the interpolation domain,
   including points adjacent to the documented `eps`-clamp
   (`abs(exp(fxy)-1)<=1e-4`) and the detailed-balance temperature switch
   (`C/Rec_database/Effective_Rates.HI/get_effective_rates.HI.cpp:259-321`).
3. First ODE pass (no diffusion correction): full state vector
   `[rho, X_{H,1s}, X_{H,i...}, X_{He,1s}, X_{He,i...}]` at every native
   output redshift, plus the derived `X_e`.
4. Hydrogen radiation-PDE distortion `u(x,z)` (or a reduced summary:
   destruction probabilities `p_d(n,z)`, emission departures
   `Delta n_{em}(n,z)`) on the native nonuniform `x`-grid, at the native
   fixed `Delta z = 10` (theta = 0.999 burn-in then 0.55) output redshifts.
5. Correction integrals `DF^{(2gamma)}(n,z)`, `DF^{(Raman)}(n,z)`,
   `DI_{1,2s}(z)` (2s-1s path) as actually produced by
   `Solve_PDEs_integrals.cpp`.
6. Second and third ODE passes (with corrections applied), same layout as
   item 3.
7. Low-z Recfast completion segment: state at the matching redshift
   (value AND the rescaled derivative used for matching, per
   `C/Modules/main.CosmoRec.cpp:249-258`) and the completed history down to
   `z=0`.
8. Final interpolated public output: `X_e(z)`, `T_m(z)` on whatever output
   grid the calling convention requests (CAMB requests values at its own
   grid via spline evaluation of the native output; capture both the
   native internal output grid AND at least one off-grid query to exercise
   that final interpolation).

## 3. Grids, columns, and precision

- Every fixture is **plain text** (CSV or whitespace-delimited), never
  binary, per repository policy.
- Header row names every column explicitly (e.g. `z X_e X_H1s X_He1s rho`),
  states units, and states the exact redshift direction/ordering used
  (ascending or descending `z`; native output order must be recorded, not
  assumed).
- Numeric precision: full `%.17g`-equivalent (round-trippable double
  precision text), not a truncated display format -- elementwise agreement
  checks are only as good as the saved precision.
- Record the exact native `z`-grid used for each fixture (e.g. the batch
  `nz=3000` linear grid from `zstart=3000` to `zend=0`, or the PDE's own
  `Delta z=10` redshift sequence) -- do not assume one fixture's grid matches
  another's.

## 4. Parameter perturbations for derivative reference

- For each differentiated input (`Omega_b`, `Y_p`, `T_0`, and selected
  `H(z)` samples -- see `VERIFIED_ANALYSIS.md` for which inputs have a direct
  dependence), capture native output at the nominal value and at a geometric
  sweep of step sizes (suggested starting sweep: relative steps
  `{1e-2, 1e-3, 1e-4, 1e-5, 1e-6}`, both signs, to support central
  differences and a step-size-convergence plot) -- **this sweep is itself a
  calibration experiment, not a pre-decided tolerance**: the step at which
  truncation error and native floating-point/solver noise cross over must be
  read off the convergence plot per quantity, not assumed in advance.
- Record the native solver's own tolerances/switches active during each
  capture (ODE Newton/step tolerances, PDE `abstol`-equivalents if exposed,
  any verbosity/debug output that reveals step counts) so that later
  Julia-side tolerance calibration has a native point of comparison.

## 5. Determinism and comparison policy

- Repeat each native capture at least twice with identical inputs and diff
  byte-for-byte (or exactly numerically) before trusting it as a reference;
  non-determinism (e.g. from the global/static caches noted in
  `VERIFIED_ANALYSIS.md`) must be caught here, not discovered later as a
  spurious Julia-vs-native mismatch.
- Comparison policy for a future Julia port against a fixture is
  **elementwise relative/absolute agreement at a calibrated tolerance**, not
  a single aggregate norm; report the full elementwise distribution
  (max/median/worst-index), and treat the documented native clamps/switches
  (§2 item 2) as regions where exact equality is the target but a
  derivative-continuity test is explicitly **not** expected to pass smoothly
  -- that discontinuity must be reported, not hidden by a tolerance chosen
  to paper over it.
- Gradient/derivative comparisons against a fixture-derived finite-difference
  reference follow the same convergence-sweep logic as Chunk 1a's own
  internal checks (see `test/chunk1a_stiff_ad_probe.jl`): ForwardDiff through
  a discretized Julia reproduction is expected to agree tightly with finite
  differences of that *same* discretization, and only to agree with a
  continuous-adjoint (Mooncake+SciMLSensitivity) sensitivity up to a gap that
  shrinks under tolerance/discretization refinement -- never assume one-shot
  machine-precision equality between a continuous adjoint and a discrete
  derivative.
