# CMB-driven ODE tolerance frontier (2026-10-04) — research only, nothing adopted

**Question (Marco, 01:53):** how far can the caller Rodas5P tolerances be pushed without the CMB becoming significantly different?

- Builds on `docs/ODE_TOLERANCE_WORK_PRECISION.md` (main's chunk21; not re-run).
- **No** source, default, helper, CI gate, clamp, grid, algorithm, dependency or API change. No tolerance was adopted. Nothing committed.
- Data, scripts and the protocol are in `A/chunk22/` (`A = cmbcheb_test/local_analysis/cosmorec_differentiability_20260930`):
  - `PROTOCOL.md` holds the tiers, which were fixed by main's 01:59 assignment before any data, plus dated addenda written before each data set they govern. Its header time was corrected.
  - `frontier_summary.{tsv,txt}` is the aggregated table.

## Setup

- **Varied:** `reltol`, `a1` (rho/ground-state abstol), `aex` (excited-state abstol). Fixed: everything else of the validated optimized callback.
- **Reference** "tight" = 1e-12 / 1e-18 / 1e-14.
- **Models:** the fiducial and four native-provenance cosmologies (ob, oc, h70, yhe), exact 10000-node grids, normal flags.
- **CMB:**
  - isolated CAMB 2.0.4 (chunk14), CosmoRec history injection, Mead2020, the chunk14 accuracy settings;
  - ell 2:3000 with the unchanged `camb_forward.py`;
  - ell 2:9000 with a local copy (`LMAX = 9000`, lmax_calc 11000; the only change, diff recorded).
- **CAMB repeatability:** two injections of the same tight history are **bitwise identical**, so the incremental metrics have a zero noise floor.

## Metrics and predeclared tiers (candidate vs tight)

- **All eight stored columns individually:** unlensed TT/EE/TE, lensed TT/EE/BB/TE, phiphi.
  - frac = max |dC|/|C| for positive spectra; TE: |dTE|/sqrt(TT EE).
  - cv = max |dC|/sigma_CV.
- **Primary:** ΔChi² of the **lensed TT/EE/TE triplet** with the per-ell 3×3 Gaussian cosmic-variance covariance (N = 2l+1):
  - Var TT 2TT²/N, Var EE 2EE²/N, Var TE (TE²+TT EE)/N, Cov TT,EE 2TE²/N, Cov TT,TE 2TT TE/N, Cov EE,TE 2EE TE/N;
  - summed over ell.
  - This is a **full-sky, noise-free, Gaussian approximation, not an experiment likelihood**: lensing non-Gaussian covariance and ell–ell correlations are ignored.
  - BB/phiphi/unlensed columns are diagnostics, not added as duplicate likelihood terms.
- **Verification:** the analyzer reproduces chunk21's `cl_compare` numbers exactly and the diagonal limit of the 3×3 form to machine precision. Main independently checked the h70/r11_a12 value with a 2×2 field-Fisher trace (identical: 2.6449569e-4).
- **STRICT (recommended target):** frac ≤ 1e-4, TE corr ≤ 1e-4, max cv ≤ 0.01 (all 8 columns), ΔChi² ≤ 0.01 (~0.1σ).
- **Also shown:** the previous chunk21 tier (frac ≤ 1e-5, cv ≤ 1e-3) and, for user choice only, ΔChi² ≤ 0.1 / ≤ 1.
- **Acceptance:** a profile counts only if STRICT holds on **all five cosmologies** and on the high-ell runs.

## Results

30 profiles were CMB-evaluated: 15 new solves plus 14 reused chunk21 fiducial histories and the reference. Work = total RHS evaluations (nf) of the three-pass
public call, hardware independent; tight averages 693,422.

| profile (reltol / a1 / aex) | mean nf (tight ÷) | worst ΔChi², ell ≤ 3000 | STRICT, 5 cosmologies, ell ≤ 3000 | worst ΔChi², ell ≤ 9000 | STRICT, ell ≤ 9000 |
|---|---|---|---|---|---|
| r6 (1e-6 / 1e-18) | — | **rejected**: solve failed, DomainError on a negative population (−5.9e-23) | — | — | — |
| r7 (1e-7 / 1e-18) | 43,608 (15.9) | 3.1e-2 (fid) | **fail** (fid) | — | — |
| r8 (1e-8 / 1e-18) | 70,706 (9.8) | 1.47e-2 (yhe) | **fail** | — | — |
| r9_a12 (1e-9 / 1e-12) | 77,072 (9.0) | 1.16e-2 (oc) | **fail** | — | — |
| r10_a12 (1e-10 / 1e-12) | 96,690 (7.2) | 3.1e-3 (yhe) | pass | **2.40e-2 (yhe)** | **fail** |
| r9_a14 (1e-9 / 1e-14) | 102,797 (6.7) | 4.8e-3 (fid) | pass | **4.10e-2 (ob)** | **fail** |
| **r11_a12 (1e-11 / 1e-12 / 1e-14)** | **114,429 (6.06)** | 2.5e-5 (h70) | **pass** | **2.64e-4 (h70)** | **pass** |
| a12 (1e-12 / 1e-12) | 128,474 (5.4) | 2.9e-6 | pass | 2.5e-5 (h70) | pass |

Other fiducial-only points (all STRICT at ell ≤ 3000; full list in `frontier_summary.txt`):
- r9 (5.7×, 9.6e-5);
- r9_a16 (6.0×, 8.9e-4);
- r9_a17 (5.8×, **1.43e-2, fail**: the error is not monotone);
- a13, a14 and the aex variants a12_x12 / a12_x10.

Relaxing aex to 1e-12 or 1e-10 changes the work by < 0.5 %, so it is not a lever.

**r11_a12 in detail** (worst over the five cosmologies, through ell 9000):
- max positive-spectrum fraction 5.0e-6;
- max cv 3.9e-4;
- ΔChi² 2.6e-4, i.e. about 40× inside the strict tier.

**The a1 axis is the efficient lever.** At reltol 1e-12, a1 = 1e-12 alone removes 81 % of the RHS work with ΔChi² ≤ 3e-6. Loosening reltol beyond 1e-11 is what produces
CMB-visible error, and that error is non-monotone and cosmology-dependent: r8 and r9_a12 pass on the fiducial but fail on yhe/oc; r10_a12 and r9_a14 pass to ell 3000 but fail at high ell.

### History-level and internal criteria (reported independently; not CMB acceptance criteria, no CI change)

- For r11_a12, the worst Xe is 1.4e-5 vs tight and 1.2e-5 vs native (h70), i.e. the **old Xe ≤ 1e-5 history budget is exceeded** at h70; Tm vs tight ≤ 9.7e-8. All faster profiles exceed it further.
- a12 stays within it (Xe ≤ 2.8e-6 vs native).
- The existing native `ob` Tm miss persists.
- The helium switch node is unchanged in every solve.
- All solves finite. Minimum state components stay positive (~6e-29) except the rejected r6.
- The original clamp is unchanged.

### Gradients (`step5_gradients_*`; existing 4-parameter route, frozen tight switch, 9 Xe + 9 Tm, direct ForwardDiff)

| profile | max elasticity difference vs tight | hscale column | nbscale column |
|---|---|---|---|
| r8 | 3.5e-6 | 6.7e-6 | 2.8e-7 |
| r9_a12 | 1.2e-5 | 6.3e-6 | 2.7e-6 |
| r11_a12 | **2.2e-6** | 1.9e-6 | 4.4e-7 |
| a12 | 2.6e-6 | 4.2e-6 | 1.7e-8 |

- All are within the 1e-4 study budget.
- The F/A2s1s columns are ~1e-18/1e-21 in absolute value, so their O(1) relative differences are not meaningful.
- This is not certification of derivatives through full cosmological initialisation or of true reverse mode through Rodas5P.

## Timing

Predeclared pair: tight vs r11_a12, via `chunk21/guarded_tolerance_bench.jl` **unchanged** (numeric-literal substitution in the validated guarded
driver; original/configured hashes recorded). FULL public call, 5 samples, evals 1, warm then in-process gate, pinned CPU 8, gate cap 10 min, unchanged
matching rule (`chunk18/match_v3.py`).

**Attempt 1** (`A/chunk22/step8_timing_*`): base **BLOCKED**, no QUIET/STABLE window within 10 min.
- Foreign > 50 %: another user's Julia (PID 3124511) and a transient python process (3158401).
- My own queued accuracy jobs had run immediately before this attempt and some were queued after it; per main, timing was deferred until all own accuracy jobs finished.
- The after process was stopped by me in its untimed cool gate (after warm-up, no timed batch). Logs are preserved; not blamed solely on outside users.

**Retry** (`A/chunk22/step10_timing_retry_*`, 02:33–02:49, no own accuracy job running; both gates QUIET). Rule verdict **REJECTED / UNMATCHED**:

| criterion | tight | r11_a12 | limit | met |
|---|---|---|---|---|
| foreign > 50 % set (ps), start / end | none / none | none / none | identical | yes |
| interval foreign busy cores (/proc/stat, reported) | 1.28 | 1.21 | — | — |
| package T start → end | 61 → 100 °C | 64 → 86 °C | ≤ 5 °C at both ends | **no** (end differs by 14 °C) |
| pinned CPU 8 frequency start → end | 4717 → 4500 MHz | 4699 → 4210 MHz | ≤ 10 % | yes |
| throttle rate (package / core) | 71.5 / 142.9 per s | 22.3 / 44.6 per s | ratio ≤ 1.5 | **no** (3.2) |

Raw observations (**not a certified speedup**):
- samples: tight 1875.00, 1907.95, 1921.09, 1911.43, 1947.13 ms; r11_a12 432.45, 450.77, 454.20, 433.33, 423.36 ms;
- median 1911.4 → 433.3 ms (observed ×0.227, ~4.4×);
- GC median 39.6 → 23.3 ms;
- memory 500.7 → 305.4 MB; allocation events 2.80 M → 0.77 M;
- cold (process start → first result) 55.7 s and 63.0 s, single observations.

The shorter candidate run heats the package less, so the thermal criteria are not met; under the predeclared rule this pair is UNMATCHED and its
ratio is not used. **What is established is hardware independent:** 6.06× fewer RHS evaluations (mean over five cosmologies) and the memory and allocation reduction. No other jobs, governor or fans were touched; no further pairs were run (bounded study).

## Bottom line (research finding, nothing adopted)

- **The credible fastest CMB-accurate profile found is reltol = 1e-11, a1 = 1e-12, aex = 1e-14.**
  - It passes STRICT (ΔChi² ≤ 0.01 and all per-column criteria) on all five tested cosmologies through ell 9000, with worst ΔChi² 2.6e-4.
  - Its gradients agree within 2.2e-6.
  - It needs 6.06× fewer RHS evaluations than tight. That is hardware-independent work, **not** a 6× wall-time claim; one observed but UNMATCHED timing pair gave ~4.4× (not certified).
- The more conservative a12 (5.4× fewer, ΔChi² ≤ 2.5e-5, and within the old history budget) is the safer alternative.
- Faster profiles fail STRICT on at least one cosmology or at high ell.

**Scope limits:**
- five cosmologies near the fiducial (not the full parameter domain);
- one CAMB configuration, ell ≤ 9000;
- a Gaussian CVL approximation, not an experiment likelihood;
- strict accuracy relative to this code's tight solution (not CosmoRec's physical accuracy, quoted by its authors at ~0.1 %);
- the old history criteria are exceeded for r11_a12 at h70.

**Recommendation:** an **opt-in** caller setting `reltol = 1e-11, a1 = 1e-12, aex = 1e-14` (the old Xe ≤ 1e-5 history budget is exceeded at h70, ~1.2e-5),
or the more conservative `reltol = 1e-12, a1 = 1e-12, aex = 1e-14`. Both pass STRICT on all five cosmologies through ell 9000 and the gradient budget.
r6 is invalid (solve failure). No default is adopted; that requires Marco's explicit decision.

Further labels: the Gaussian full-sky covariance approximation is not an exact ACT/Planck/SO likelihood or significance; five sample points are not the full
parameter domain; the gradients use the frozen-background 4-parameter route, not B2/B3 true reverse mode. The source tree is unchanged by this study and UNCOMMITTED.
