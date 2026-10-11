# CosmoRec.jl forward CAMB replacement plan

Prepared by Claude Opus, high reasoning, 2026-10-11, session
`1452d25c-4cb0-4204-9368-77638dbc06a0`. Controller transcription of
the executable plan. Implementation resumes on Opus, medium reasoning.

## FINAL GOAL (corrected by Marco, 2026-10-11; supersedes the bitwise "one-to-one" reading)

"One to one" means the **CAMB contract**, not floating-point identity:

- CosmoRec.jl takes the same inputs from CAMB as CosmoRec (`cosmorec.f90:165-213`: Omega_c, Omega_b, Omega_k, `N_eff()`, H0, tcmb, Yhe, CAMB's z/H(z) table).
- It returns Xe/Tm to CAMB.
- Swapping CosmoRec for CosmoRec.jl must produce no actual difference.

**Acceptance gate:** pointwise `|C_l^Julia / C_l^native - 1| < 5e-4` in BOTH TT and EE at EVERY 2 <= l <= 9500.
- Report lensed (the production output) and unlensed separately.
- Gate verdict is per point, with the max and its argmax l.
- Conditions:
  - identical production CAMB settings (worker preset: CAMB 2.0.4 CosmoRec runmode 0, Mead2020 `NonLinear_both`, `DoLateRadTruncation` False, `lens_potential_accuracy` 8, `lens_output_margin` 2050, kmax 10, `k_per_logint` 130, AccuracyBoost 1, lSampleBoost 1, lAccuracyBoost 1.2, `min_l_logl_sampling` 6000, PPF, TCMB 2.7255, lmax 9500);
  - same inputs, generation domains and neutrino/YHe conventions as the workers (`C36/p4/WORKER_AUDIT.md`);
  - paired native vs Julia with identical CAMB numerics.
- Order: fid + 2 corners first, then the relevant domains.
- TE/BB/PP are reported as diagnostics; no threshold is invented for them.
- Native reference outputs are saved as text.

## EXECUTION STATUS AND REMAINING PHASES (2026-10-11, main agent; for the implementation model)

Phases 1-3 are DONE and PASSED. Frozen implementation, hardening, exit contracts, and the first full-l certification are complete:

- Frozen implementation record: `cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk36_camb_replacement/bridge/IMPLEMENTATION_FREEZE.md`.
- 3-point certification (fid + both domain corners, worker preset, lmax 9500): **PASS at every ell 2..9500.** Worst number:
  cornerB lensed EE 7.8e-6 (64x under 5e-4). Numbers and evidence pointers:
  `.../chunk36_camb_replacement/bridge/cmb_pair/CAMPAIGN_RESULTS_3POINT.md`.
- Focused bridge gate: 16/16 (`bridge_test_run4.log`), incl. repeated A-B-A calls in one process and bitwise determinism.
- Measured: warm Julia recombination ~1 s/cosmology (first call in a process ~39 s compile); native CAMB ~15 s/cosmology.

### Phase A - 36-point certification matrix (next)

Extend the paired campaign from 3 to 36 points covering the three production worker domains (`C36/p4/WORKER_AUDIT.md`):
3 representative baselines, 18 interior points (6 per domain, spread over H0, omega_b, omega_c, Mnu, w0+wa<-0.5 boundary interior,
ns, tau, YHe/Alens), 15 boundary points (domain corners and edges, incl. BBN-YHe cases with `YHe: null`, mnu 0 and max,
degenerate hierarchy with `neutrino_hierarchy: "degenerate"`, `standard_neutrino_neff: 3.044`, Alens 0.7 and 1.4).
Reuse `camb_forward_c36.py` preset `worker` and `run_cmb_pair.sh` (exit contracts already hard: 0/1/2). Generate point JSONs in
`bridge/cmb_pair/`, run serial, one heavy process, report per-point gate lines. DO NOT modify the gate or thresholds.

### Phase B - operational reliability

Repeated same-point runs in one process (A-B-A at worker preset), memory retention across >= 50 calls (bounded capture is in place),
missing-tables and bad-input failure behavior (must abort loudly, never fall back), and the interaction of the shim's fail-fast
abort with the distributed generation worker pattern (128 Julia workers, one CAMB each): document how a worker failure surfaces.

### Phase C - packaging into the repo

Move the validated bridge from the scratch analysis tree into maintained locations: shim source + build script + bridge setup +
production callbacks (or a package-side equivalent) with pinned env, explicit `COSMOREC_TABLES_ROOT` (user-provided, never
redistributed), native-vs-Julia build selector that does not touch scientific parameter construction. Then run the existing
50-cosmology generation smoke test end-to-end (temporary outputs outside git repos).

### Phase D - matched benchmark (only after correctness acceptance)

Native vs Julia under identical physics/settings/hardware: warm recombination cost (BenchmarkTools, prepared solvers), startup
and table-load cost, per-worker throughput, total CAMB wall. Serial, pinned-core, thermally stabilized (laptop throttling caveat).
Exclude diagnostic capture from production timings. Report startup separately; never compare native cold vs Julia warm.

### Standing constraints

Original trees read-only (`tools/CosmoRec` 0867690, `tools/CAMB-cosmorec`); no commits/pushes/tags; no gate or reference weakening;
AD correctness is a separate scope and stays visible; one heavy process at a time on Alcanter; Claude = Sonnet, medium reasoning,
bounded chunks, no automatic retries on quota.

Then a matched native/Julia benchmark (only after correctness acceptance).

**Scope separation:**
- Internal population/ULP gates, native bitwise parity (Gear port, spline order) and the six old gates are DIAGNOSTICS, not forward-replacement blockers.
- The Gear work is kept for review, but it is not expanded and not made the default.
- AD failures are reported separately and not declared resolved.
- The existing production callbacks are acceptable if the CMB gate passes: preliminary Rodas5P dense 1e-12; passes `solve_phase5` = Rodas4P 1e-10/1e-18/1e-12 with tstops; tail `solve_tail5`.

**Delivery proof** (the diagnostic paired-history injection is NOT delivery):
- CAMB -> CosmoRec.jl -> CAMB in the production process model (Julia -> PyCall -> camb -> camblib.so -> shim -> Julia), OMP 1;
- Julia call counter equal to the shim counter;
- CAMB using the returned arrays;
- native CosmoRec neither linked nor invoked;
- no fallback on failure.

Implementation (scratch, C36/bridge):
- `shim/cosmorecjl_shim.c` replaces `libCosmoRec.a`;
- `build_camb_bridge.sh`: original CAMB rsync, no source patch, no GSL;
- `cosmorecjl_bridge.jl`: the `@cfunction` entry over `cosmorec_calc_h`, which re-acquires the GIL;
- `bridge_test.jl`: focused proof;
- `cmb_pair/run_cmb_pair.sh`: certification launcher.

## Baseline and invariants

Main checkout `develop`, HEAD `171da92`; preserve existing chunk35 edits.
Original CosmoRec `0867690` and original CAMB sources are read-only. Scratch
copies may be patched with recorded diffs. No staging, commits, pushes,
releases, resets or worktrees. No unrelated emulator-zoo edits.

Analysis root: `cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/`.
New artifacts: `chunk36_camb_replacement/` (C36 below).

Focused tests first; unconditional suite inclusion; native fixtures as text.
One heavy process at a time. Engineer windows <=15 minutes, checkpoints each
window. Full suites run under controller-owned systemd/tmux execution, outside
Claude's process tree, not under a Claude timeout or merely nohup.
No acceptance gate relaxation or passing claims for documented failure envelopes.

## Open decisions / genuine blockers

* Native README says all rights reserved and permits use with acknowledgment
  and citations; no explicit redistribution/modification grant was found.
  Table bundling and public licensing of translated code require clarification
  from the rights holder. Local user-provisioned data can be used meanwhile.
* If preliminary-history differences are proven native solver error, faithful
  Recfast++ Gear port versus converged-native scientific reference requires a
  decision. The latter changes the parity reference, not automatically approved.
* Production high-ell acceptance needs explicit approval. Proposed STRICT tier
  from existing tolerance docs: frac <=1e-4, TE correlation metric <=1e-4,
  maximum CV residual <=0.01, diagonal sum <=0.01. This proposal is not an
  approved replacement for the original history gates. Native-floor ratio is
  a diagnostic only.
* CAMB bridge failure semantics and Julia foreign-OpenMP-thread safety need
  design evidence. Do not silently swallow errors or change CAMB threading.
* Keep the native H-table node-0 behavior unless a change is explicitly approved.

## 1. Current implementation validation

Last complete baseline: chunk33 suite 59665 passed / 6 failed / 59671,
about 160 minutes. Failed sites: chunk10a lines 83,119,120; chunk10b lines
42,56; chunk10d line46. Read baseline evaluated values, not just counts.
Chunk35 full suites were terminated and have no accepted tally.

1. Re-run focused chunk35 and save count/current-file evidence.
2. Copy verified chunk33 `full_suite/run.sh` (change OUT only) and run.jl to
   C36/full_suite. It uses >=20GiB admission, timeout18000, Julia threads1,
   heap hint8G, Pkg.test, and terminal outcome files. Controller launches it.
3. Acceptance: exactly prior failure sites/values, zero errors, chunk35 zero
   failures, expected pass count baseline plus new tests. Diagnose other results.
4. Save TALLY.md and append accurate section9 to chunk35 REPORT.md.

Resources: focused ~1min/1.3GiB; full suite ~160min/12GiB. Terminal files:
suite.done (exit_code/finished), suite.outcome, suite.log. Exit1 may represent
the expected six failures, not success; exit124 is timeout.

## 2. Preliminary Recfast++ discrepancies

1. Compare RHS at native node states, including z<50 lo/h50; compare grid,
   Saha breakpoints and interpolation. RHS/grid scaled gates 1e-12.
2. Julia convergence Rodas5P reltol1e-10,1e-12,1e-14 and 128-bit reference
   if supported; report nodewise low-z spread, diagnostic target1e-9 relative.
3. Scratch native build: reduce Recfast absolute tolerances by1e-6 and relative
   tolerance to1e-12, modifying only recorded scratch solver constants.
   Compare original/tight native and converged Julia. Native-error hypothesis
   supported only if tightened native moves toward Julia (target1e-8 relative)
   and native shift explains observed 3.9--9.2e-5 differences. Otherwise bisect
   per-node/state/step discrepancy to translation bug.
4. Only after evidence, decide faithful 3+1-state Gear port (~1000 native lines,
   8--12 windows) versus changed converged reference. Preserve original gates.

Artifacts C36/p2; ~10min heavy diagnostics/2GiB before possible Gear port.

## 3. Main history and prior six failures

1. Stage localization fid/ob/oc/lo/h50/neff via native pass rows into Julia PDE
   and native PDE outputs into Julia passes. Save first failing stage and
   native-tolerance units, C36/p3/stage_table.md.
2. Native self-spread: 1ulp input perturbations, recorded scratch build options;
   compare Float64 vs128-bit Julia public map where supported. For AD failures,
   finite-difference step sweep versus ForwardDiff on the frozen map. Conditioning
   claim requires measured self-spread, not rhetoric.
3. Low-ell localization via hybrid native/Julia histories split at z50,200,500
   through identical CAMB; capture reionization tau-to-z_re readout.
4. Fix proven translation errors with unconditional native-text-fixture tests.
   Conditioning-limited failing gates remain failures pending human decision.

~6--10 windows, ~40min serial heavy work/2GiB.

## 4. Production domain / high-ell certification

Read current three CAMB workers and actual prior ranges, then record a finite
configuration matrix C36/p4/config_matrix.tsv (~72 points target, not an
exponential corner expansion). Include representative permitted corners,
16 interiors, Mnu0/max, w0wa support, BBN CP%Yhe mass fraction, explicit-YHe,
Alens and each campaign's neutrino convention.

Production settings: CAMB2.0.4, CosmoRec, Mead2020, NonLinear_both,
DoLateRadTruncation=false, lens_potential_accuracy8, lens_output_margin2050,
kmax10, k_per_logint130, AccuracyBoost1, lSampleBoost1, lAccuracyBoost1.2,
min_l_logl_sampling6000, PPF, lmax9500, TCMB2.7255, omk0, Reion.set_tau.
yhe: nnu3.046/one massive/explicit helium; degenerate: standard3.044/three
degenerate or all massless for Mnu0/BBN; alens: fixed mass/one massive/BBN.

Use CAMB-internal dump inputs to remove Python-grid reference floor; paired
identical CAMB native-history vs Julia-history runs first, integrated bridge
later. Store dense references as text. Validate ALL lensed/unlensed spectra
TT/TE/EE/BB/PP through9500. TE full Knox variance
(TE^2+TT*EE)/(2ell+1), plus correlation metric; no fractional singularities.
Pointwise maxima, CV residuals and diagonal sums. Audit PP normalization.
Acceptance threshold unresolved above; no expensive training or generation.

~60min serial heavy/1.5GiB. Do not treat injection certification as swap-in proof.

## 5. Reproducible provisioning

Bridge project with explicit pinned solver dependencies/callbacks, user-supplied
table root and per-table path/hash manifest. No implicit root fallback/download.
Keep package solver policy explicit. License audit docs with exact native terms,
citations and GSL linkage considerations. Bundling blocked without permission.
~4 windows, lightweight work.

## 6. Real CAMB swap-in

Source-backed proposal (verify feasibility before final design): link-compatible
C shim replacing libCosmoRec, exact 17-argument `cosmorec_calc_h_cpp_` symbol,
registered Julia callback, no native fallback. Scratch CAMB built against shim,
not native lib. Python driver initializes persistent Julia and tables explicitly;
Julia unsafe_wrap boundary validates inputs, writes outputs, contains exceptions.
Callback status/error semantics and foreign OpenMP-thread adoption must be proven
before use; native CAMB invokes Init in OMP sections with maxthreads>1.

Evidence gates: symbol/link checks show no native Cosmos/CosmoRec/Xe_frac/GSL;
Julia call counter matches initialization including get_zstar; output arrays
match injection path bitwise for identical inputs; 50 repeated cosmologies equal
fresh-process results; OpenMP1 vs8 safely agrees; error path declared/tested;
200-call RSS/lifecycle check. No claim of swap-in if only diagnostic injection.
Then rerun production-domain/high-ell matrix actual bridge vs original CAMB.
~10--14 windows/~60min heavy/3GiB.

## 7. Matched benchmarks AFTER correctness acceptance

Native in-process C++ loop (A7-style) and Julia BenchmarkTools with interpolated
inputs, evals1, >=30 samples; same physics and analysis precision. Single-thread,
pinned core selected on actual machine, thermal stabilized, >=10 alternating
matched pairs, no concurrent heavy workloads. Measure warm full solve, setup,
first-result cold cost, 200-call persistent throughput, full integrated CAMB.
Report medians/spread/allocations and timed-output correctness. Separate solver
work/RHS counts from translation cost. Prior .85s single sample is not evidence.
Artifacts C36/p7; ~90min heavy/2GiB after all gates.

## Controller order

Phase1 runner/focused tests -> controller suite (no other heavy job) -> tally;
lightweight planning/provisioning audit can proceed during suite. Then phases2/3
diagnostics, decision gate if needed, integration/provisioning design, actual
replacement, production certification and benchmarks. Parent independently
reviews scientific changes; Claude self-reports are not scientific acceptance.
