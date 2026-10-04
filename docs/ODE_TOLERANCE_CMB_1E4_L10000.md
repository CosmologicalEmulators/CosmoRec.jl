# Loosest ODE tolerances with CMB relative error ~1e-4 through ell 10000 (2026-10-04) — research only, nothing adopted

**Request (Marco, 15:06 EDT):** find the loosest caller Rodas5P settings that keep the CMB **relative** error at about 1e-4 for the same input
parameters, through ell 10000.

- Primary target: the maximum pointwise fraction over ell 2:10000. The earlier Gaussian ΔChi² ≤ 0.01 criterion
  (`docs/ODE_TOLERANCE_CMB_FRONTIER.md`) is reported as a diagnostic only and does **not** veto a fraction-compliant profile.
- No source, test, helper, default, gate, math, clamp, grid, algorithm, dependency or API change. Nothing adopted or committed.
- Data, scripts and the protocol: `A/chunk23/` (`A = cmbcheb_test/local_analysis/cosmorec_differentiability_20260930`). `PROTOCOL.md` was written before any
  chunk23 data, with dated addenda before each data set; the tiers never changed. `summary23.{tsv,txt}` is the aggregated table.

## Setup

- **Varied:** reltol / a1 (rho/ground-state abstol) / aex (excited abstol) of the validated optimized callback; reference tight 1e-12 / 1e-18 / 1e-14.
- **Models:** the fiducial and native-provenance ob, oc, h70, yhe (same input parameters per case), 10000-node grids. Normal flags; the finalists are repeated with
  `--check-bounds=yes` against a tight run under the same flags.
- **CMB:** isolated CAMB 2.0.4, CosmoRec injection, Mead2020, chunk14 accuracy settings, via a local copy of `camb_forward.py` with LMAX 10000 /
  lmax_calc 12000 (the only change).
  - Every output was checked to have **10001 rows (ell 0..10000)** and metadata `lmax_out = 10000`, `lmax_calc = 12000`, and to be finite.
  - The 3k native spectra are not compared with these.
- Reused chunk22 histories for already-solved profiles; 10 new profile solves (cap 15–20).

## Floors (kept separate)

- **CAMB repeat:** two injections of the same tight history are bitwise identical (0). This is a deterministic repeat, **not** proof of accuracy.
- **Solver convergence of the reference** (tight vs tighter = 1e-13 / 1e-19 / 1e-15, 1.36 M RHS evaluations), ell ≤ 10000:
  - positive-spectrum fraction 1.07e-6, TE corr 8.3e-7;
  - **raw TE fraction 4.3e-4**, at ell ~1258 next to a TE zero crossing.
  - Candidate errors (~1e-5 to 1e-4) are 10–100× above the positive-spectrum floor.

## Tiers (predeclared) and the literal-tier result

- **LITERAL:** max |dC|/|C| ≤ 1e-4 for all TT/EE/BB/phiphi columns **and** raw |dTE|/|TE| ≤ 1e-4 wherever TE ≠ 0. **None of the screened non-trivial
  profiles passes**; the tight reference vs the tighter solve also gives 4.3e-4.
  - The raw TE fraction is mathematically defined wherever TE ≠ 0 (undefined only at exact zeros), but **ill-conditioned** and a poor physical metric
    next to TE's 11 sign changes, where |TE|/sqrt(TT EE) ~ 1e-4–1e-3.
  - The tight-vs-tighter value shows the limits of the reference and of the raw metric near zero. It is **not** proof that no tolerance could pass
    (tight vs itself is trivially zero).
  - What can be stated: "every spectrum, literally, within 1e-4" holds for **none of the screened profiles**.
- **PRACTICAL (zero-safe TE):** positive columns ≤ 1e-4 and TE corr = |dTE|/sqrt(TT EE) ≤ 1e-4.
- **GUARD-BAND:** PRACTICAL at 8e-5.
- Acceptance at a tier requires all five cosmologies through ell 10000.

## Results (max over the eight columns; worst over the cosmologies tested; work = mean RHS evaluations, tight averages ~694 k)

| profile (reltol / a1; aex 1e-14) | cosmologies | RHS ÷ vs tight | worst positive frac (where) | worst TE corr | worst raw TE frac | PRACTICAL | GUARD |
|---|---|---|---|---|---|---|---|
| r8_a10 (1e-8 / 1e-10) | fid | 18.5 | 1.88e-4 | 7.6e-5 | 4.7e-2 | **fail** | fail |
| r3e9_a10 (3e-9 / 1e-10) | 5 | 16.6 | 1.08e-4 (oc) | 3.9e-5 | 1.9e-1 | **fail** (oc) | fail |
| r7 (1e-7 / 1e-18) | fid | 15.9 | 1.90e-4 | 1.0e-4 | 7.2e-2 | **fail** | fail |
| **r9_a10 (1e-9 / 1e-10)** | 5 | **15.3** | **8.33e-5 (oc)** | 4.3e-5 | 2.0e-1 | **pass** | fail (oc) |
| — same, `--check-bounds=yes` | 5 | 15.4 | 8.68e-5 (oc) | 3.4e-5 | 2.4e-1 | **pass** | fail (oc) |
| r8_a12 (1e-8 / 1e-12) | 5 | 12.6 | 2.48e-4 (ob) | 8.4e-5 | 6.4e-1 | **fail** (ob, oc, yhe) | fail |
| r3e8 (3e-8 / 1e-18) | fid | 12.5 | 8.01e-5 | 2.9e-5 | 2.3e-2 | pass | fail |
| **r9_a11 (1e-9 / 1e-11)** | 5 | **11.6** | **6.69e-5 (ob)** | 3.6e-5 | 1.2e-1 | **pass** | **pass** |
| — same, `--check-bounds=yes` | 5 | 11.6 | 5.60e-5 (ob) | 3.2e-5 | 5.4e-2 | **pass** | **pass** |
| r2e8 (2e-8 / 1e-18) | fid | 11.5 | 4.78e-4 | 1.7e-4 | 2.2e-1 | **fail** | fail |
| r3e9_a12, r9_a12, r10_a12, r9_a14, r9, r8 | fid | 6.8–10.4 | ≤ 8.5e-5 | ≤ 3.3e-5 | ≤ 3.4e-2 | pass | pass (r8 fails GUARD) |
| r11_a12 / a12 (earlier recommendations) | fid | 6.1 / 5.4 | 1.4e-6 / 1.5e-6 | ≤ 4.7e-7 | ~1e-3 | pass | pass |

**Observations:**
- aex relaxation (r9_a10 with aex 1e-12) changes the work by ~1 %, so it is not a lever.
- The error is **non-monotone and cosmology dependent:** r3e8 passes but r2e8 fails badly; r8_a12 and r3e9_a10 pass on the fiducial but fail on other cosmologies.
- The worst pointwise errors of the finalists sit in the **high-ell damping tail**, mostly unlensed TT/EE at ell ~5000–9200.
  - The earlier ell ≤ 3000 definition missed this range.
  - The ell ≤ 9000 runs covered most of it, but with different profiles. The additional ell 9001–10000 range can also move the boundary.
- **Health:** all solves finite with successful retcodes; minimum state components positive; helium switch nodes unchanged; clamp unchanged.

**Diagnostics (not acceptance):**
- Gaussian lensed TT/EE/TE ΔChi² (3×3 CVL, full-sky approximation, not an experiment likelihood), worst over the cosmologies:
  - r9_a10: 7.8e-2 (normal) / 0.135 (check-bounds), mostly from ell 6001–10000;
  - r9_a11: 6.2e-2 / 5.1e-2.
- No profile exceeds ΔChi² = 1 (no WARNING). ΔChi² 0.08–0.14 corresponds to sqrt(ΔChi²) ≈ 0.28–0.37 CVL distance units in an ideal full-sky
  experiment to ell 10000. This is a significance caveat, not a veto of the pointwise user target.

**History-level** (reported separately; not criteria for this target): every profile faster than r11_a12 exceeds the old Xe ≤ 1e-5 history budget
(e.g. Xe vs tight up to 3.3e-4 for r9_a10 at oc). These are ~1e-4-level recombination-history differences whose CMB effect is the table above.

## TE Knox residuals (explicit user preference, 2026-10-04) — analysis by Claude, `A/chunk23/te_knox_20261004T161627/`

Marco prefers judging TE by its **full Knox variance** (with the TT and EE contributions) rather than by dTE/TE, which explodes at TE's zero crossings.
Tiers were added to `PROTOCOL.md` before computing. No new CosmoRec/CAMB run: all existing ell-10000 spectra were reused, after checking 10001 rows,
`lmax_out`/`lmax_calc` = 10000/12000 and finiteness.

**Definition.**
- General form, uncorrelated beam-deconvolved noise: sigma_TE² = [C_TE² + (C_TT+N_TT)(C_EE+N_EE)] / [(2ell+1) f_sky Δell].
- **Here:** f_sky = 1, Δell = 1, N = 0 (cosmic-variance limit), reference = tight. **No survey noise or beam is assumed.**
- Residual r_TE(ell) = (TE_cand − TE_tight)/sigma_TE, for lensed and unlensed TE, ell 2..10000.
- Three different quantities, never equated: **Knox σ units**; TE corr = |dTE|/sqrt(TT EE); the raw fraction |dTE|/|TE|. In particular, a 1e-4 *fraction* target is not a 1e-4 *σ* target.
- **Independent check:** a loop-based recomputation in plain Python floats (different evaluation order) agrees with the vectorized code to **4.3e-16** relative over all 45 profile/case curves (`independent_check.txt`, including the worst ells).
- Main's earlier values are confirmed: r9_a10 0.004544 σ (oc, unlensed TE, ell 9023); r9_a11 0.003843 σ (ob, unlensed TE, ell 9009).

**Tiers (user choice, not a default):** max |r_TE| ≤ 0.001 / 0.005 / 0.01 σ (both TE columns), combined with the **unchanged** positive-spectrum criterion
max |dC|/|C| ≤ 1e-4 (TT/EE/BB/phiphi, unlensed + lensed), on all five cosmologies (and both flags where run). The raw TE fraction is kept only as a historical diagnostic, not a veto.

| five-cosmology profile (aex 1e-14) | RHS ÷ vs tight | worst positive frac (both flags) | worst TE Knox (both flags) | ≤ 0.001 σ | ≤ 0.005 σ | ≤ 0.01 σ |
|---|---|---|---|---|---|---|
| r3e9_a10 (3e-9 / 1e-10), normal flags only | 16.6 | 1.08e-4 (**fails** positive) | 3.46e-3 | no | no (positive) | no (positive) |
| **r9_a10 (1e-9 / 1e-10)** | **15.3** | 8.68e-5 | **4.54e-3** (oc uTE ell 9023, normal) | no | **yes** | **yes** |
| r8_a12 (1e-8 / 1e-12), normal flags only | 12.6 | 2.48e-4 (**fails** positive) | 7.96e-3 | no | no (positive) | no (positive) |
| **r9_a11 (1e-9 / 1e-11)** | **11.6** | 6.69e-5 | **3.84e-3** (ob uTE ell 9009, normal) | no | **yes** | **yes** |

**Fastest tested candidate per TE budget:**
- **≤ 0.005 σ and ≤ 0.01 σ: `reltol 1e-9, a1 1e-10, aex 1e-14`** (~15.3× fewer RHS evaluations; positive spectra 8.7e-5, no margin vs the 8e-5 guard band). The guard-banded alternative `a1 1e-11` (~11.6×) has a worst TE of 3.8e-3 σ.
- **≤ 0.001 σ: not established within this study.** No five-cosmology ell-10000 profile qualifies. On the fiducial only, several slower profiles have max TE Knox below 0.001 σ (r3e9_a12 8.9e-4, 10.4×; r9 7.6e-4, 5.7×; r11_a12 4.6e-5 and a12 5.1e-5); r10_a12 (1.08e-3) and r9_a12 (1.31e-3) do not. None of them was run on all five cosmologies at ell 10000, so no 0.001 σ recommendation is made.
- The solver floor (tighter vs tight) is 8.7e-5 σ, and the CAMB repeat is 0.

**Per-case detail for the finalists** (max |r| and its ell; Σr² over ell, the diagonal TE-only sum; full rows incl. RMS, cumulative ranges, TE corr and raw fraction in `te_knox_cases.tsv`; curves in `curves/<flags>_<case>_<profile>.txt`):

| profile, flags | fid | ob | oc | h70 | yhe |
|---|---|---|---|---|---|
| r9_a10, normal: uTE max (ell) / Σr² | 1.73e-3 (8945) / 6.1e-3 | 9.0e-4 (7659) / 9.1e-4 | **4.54e-3 (9023)** / 4.8e-2 | 1.80e-3 (8909) / 6.7e-3 | 2.31e-3 (7651) / 1.4e-2 |
| r9_a10, normal: lTE max | 4.8e-4 | 3.4e-4 | 1.20e-3 | 3.4e-4 | 5.5e-4 |
| r9_a10, check-bounds: uTE max / Σr² | 1.61e-3 / 4.7e-3 | 3.11e-3 / 2.2e-2 | 2.59e-3 / 8.1e-3 | 1.68e-3 / 4.5e-3 | 3.63e-3 / 2.7e-2 |
| r9_a11, normal: uTE max (ell) / Σr² | 1.97e-3 (7677) / 9.9e-3 | **3.84e-3 (9009)** / 3.5e-2 | 1.69e-3 (8224) / 4.5e-3 | 2.19e-3 (5128) / 6.4e-3 | 1.72e-3 (4594) / 6.2e-3 |
| r9_a11, check-bounds: uTE max / Σr² | 2.61e-3 / 1.6e-2 | 3.41e-3 / 2.4e-2 | 1.78e-3 / 4.0e-3 | 2.03e-3 / 1.1e-2 | 1.21e-3 / 2.1e-3 |

- The worst TE residuals are in **unlensed TE in the damping tail** (ell ~4500–9000); lensed TE stays ≤ 1.2e-3 σ for both finalists.
- The raw TE fraction reaches 0.20 (r9_a10) and 0.24 (check-bounds) at TE zero crossings. Those spikes are **not** physical σ significance.

**Interpretation limits:**
- The per-ell maximum Knox residual (≲ 0.0045 σ), the diagonal TE-only sum (Σr² ≤ 0.048), the joint lensed TT/TE/EE 3×3 Gaussian distance (ΔChi² 0.08–0.14, i.e. ~0.28–0.37 CVL units) and an actual ACT/SO likelihood are **different quantities**.
- The Gaussian forms ignore non-Gaussian lensing covariance and ell–ell correlations. The CVL limit with N = 0 is the most demanding case; a real survey's noise would increase sigma_TE.

## Gradients (existing 4-parameter route, frozen switch, 9 Xe + 9 Tm; vs tight)

| profile | max elasticity difference | hscale column | nbscale column |
|---|---|---|---|
| r9_a10 | **4.9e-6** | 5.2e-6 | 7.3e-7 |
| r9_a11 | **5.2e-5** | 4.4e-5 | 9.7e-6 |

- Both are within the 1e-4 budget.
- F/A2s1s columns are ~1e-18/1e-21 in absolute value, so their relative differences are not meaningful.
- Not a certification of complete-initialisation or true-reverse (B2/B3) derivatives.

## Timing

Predeclared single pair: tight vs r9_a10 via `chunk21/guarded_tolerance_bench.jl` unchanged (FULL public call, 5 samples, warm then in-process gate,
pinned CPU 8, gate cap 10 min, unchanged matching rule), run with no own CMB/gradient job active (`A/chunk23/step7_timing_*`, 15:38–16:00).

**Result: ENVIRONMENT BLOCKED, both processes.** Neither reached a QUIET or STABLE window within its 10-minute gate (60 samples each), so **nothing was timed**
and no ratio exists.
- No foreign process > 50 % CPU was listed.
- The package throttle counters advanced in every 10 s window (24–419 per window).
- The package temperature fluctuated between 57 and 93 °C.
- Measured interval foreign load was 0.47–0.85 busy cores.
- Cold start to first result: 51.7 s / 58.8 s (single observations).

No other jobs, governor or fan settings were touched; no further pairs were run. **Established:** ~15.3× (r9_a10) and ~11.6× (r9_a11) fewer RHS evaluations
(hardware independent). Wall-time gain is **not** certified; the earlier unmatched observation for r11_a12 (~4.4×) is not transferable.

## Answer (research finding; nothing adopted)

- **Loosest profile meeting the 1e-4 relative target (positive spectra pointwise, TE zero-safe) on all five cosmologies through ell 10000, under both
  compiler flags:**
  - `reltol = 1e-9, a1 = 1e-10, aex = 1e-14`;
  - ~15.3× fewer RHS evaluations than tight;
  - worst 8.3e-5 (normal) / 8.7e-5 (check-bounds) at oc, i.e. **at the boundary**, failing the 8e-5 guard band.
- **With a guard band (≤ 8e-5):** `reltol = 1e-9, a1 = 1e-11, aex = 1e-14`, ~11.6× fewer RHS evaluations, worst 6.7e-5 / 5.6e-5.
- **With TE judged by its full Knox variance** (user preference; CVL, f_sky 1, no noise), the same `1e-9 / 1e-10` profile is the fastest tested
  candidate for TE budgets of 0.005 σ and 0.01 σ (worst 4.5e-3 σ); `1e-9 / 1e-11` has 3.8e-3 σ. A 0.001 σ TE budget is **not established** by the
  five-cosmology runs (only slower, fiducial-only profiles reach it).
- Faster settings tested fail on at least one cosmology.
- Because the error is non-monotone in the tolerances and varies between cosmologies, the guard-banded `a1 = 1e-11` profile is the safer choice; the
  `a1 = 1e-10` profile has essentially no margin at oc.

**Limits:**
- five cosmologies near the fiducial, not the full parameter space;
- one CAMB configuration;
- the literal all-spectra 1e-4 claim holds for none of the screened profiles; raw TE is ill-conditioned near its zero crossings;
- the reference is constant (tight 1e-12/1e-18/1e-14, the same input parameters per cosmology);
- the recommendation is conditional on the PRACTICAL (zero-safe TE) interpretation, not on the raw-TE literal goal;
- relative to this code's tight solution, not CosmoRec's physical accuracy;
- RHS-evaluation ratios are hardware-independent work, not wall time.

Adoption as a default or preset requires Marco's explicit decision.
