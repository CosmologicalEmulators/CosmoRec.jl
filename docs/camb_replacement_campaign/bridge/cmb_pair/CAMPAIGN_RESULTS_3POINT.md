# 3-point CMB certification: PASSED (2026-10-11)

First real end-to-end run of the frozen bridge implementation (`../IMPLEMENTATION_FREEZE.md`): original native CosmoRec vs
CosmoRec.jl through the bridge CAMB, identical worker-preset production settings (CAMB 2.0.4, lmax 9500, Mead2020 NonLinear_both,
lens_potential_accuracy 8, lens_output_margin 2050, kmax 10, k_per_logint 130, AccuracyBoost 1, lSampleBoost 1, lAccuracyBoost 1.2,
min_l_logl_sampling 6000, DoLateRadTruncation False, PPF, TCMB 2.7255, OMP 1).

**Gate (Marco 2026-10-11): pointwise |C_Julia/C_native - 1| < 5e-4 for lensed TT and EE at EVERY ell 2..9500, per cosmology.**

- fid (H0 67.36, LCDM, mnu 0.06): lensed TT max 7.674e-07 at ell 3429; lensed EE max 2.154e-06 at ell 9500. PASS.
- cornerA (H0 50, w0 -3, wa 2, ns 1.05, YHe domain minimum): lensed TT max 2.491e-06 at ell 2; lensed EE max 2.066e-06 at ell 9500. PASS.
- cornerB (H0 90, mnu 0.5, ns 0.85, YHe domain maximum): lensed TT max 9.890e-07 at ell 8739; lensed EE max 7.838e-06 at ell 9500. PASS.

Worst number anywhere: cornerB lensed EE 7.8e-6 at ell 9500 — **64x under the 5e-4 budget**, at the extreme corners of the
generation domains. Unlensed TT/EE also pass the same threshold (max 2.5e-6). Diagnostics (no threshold): lensed TE
max|dTE|/sqrt(TT*EE) 4.4e-7, lensed BB max frac 7.7e-6, phiphi max frac 1.4e-8.

## Timing (laptop, serial, OMP 1)

- Native CAMB wall: 15.5 / 15.7 / ~15 s per cosmology (fid / cornerA / cornerB).
- Julia route: first call in a process 39.4 s (one-time SciML compilation), **warm calls 0.91 / 1.16 s**; CAMB+recombination
  total elapsed 54.0 s (first) then 14.8 / 20.6 s. Per-worker steady-state overhead vs native is ~1 s of recombination.

## Evidence

- `gate.jsonl` (machine records), `gate.log` (per-spectrum lines), `compare.log` (cl_compare diagnostics incl. cosmic-variance),
  `*_native/julia_cls.txt` (9501-row spectra), `*_meta.json` (params, camblib sha, derived YHe/N_eff), `*_bridge.txt` (call counts,
  per-call seconds, inputs as received), `run_header.txt` (hashes, versions), `campaign.log`, per-step `.exit` files (all 0).
- Focused feasibility gate `../bridge_test_run4.log`: 16/16 (roundtrip, counters, native-unlinked, unregistered-abort,
  repeated A-B-A calls, bitwise determinism).
