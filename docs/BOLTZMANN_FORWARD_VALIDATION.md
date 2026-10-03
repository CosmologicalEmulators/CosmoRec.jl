# Boltzmann-facing forward validation (CAMB C_ℓ), fiducial + 4 nearby cosmologies — 2026-10-02/03

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al.

**Scope and paths.**
- Paths are relative to `A = cmbcheb_test/local_analysis/cosmorec_differentiability_20260930`. All outputs are text in `A/chunk14/`, outside git.
- No package source, test, gate, dependency or native file was changed. The original CosmoRec tree/library/data, the active CAMB install, the conda env and the shared depot were used read-only. No Git operation was performed.
- Package code is unchanged in this chunk, so the last full suite (50812 passed / 6 failed, `A/chunk12/pkgtest_a1_run1.log`) still describes it, and it was not rerun.

## 1. What was compared, and how

**Inputs.** For each cosmology, native CosmoRec (runmode 0, batch default, the CAMB-adapter route `cosmorec_calc_h_cpp_`) and the Julia port (public `recombination_history_diffusion`) receive identical inputs:
- the same CAMB H(z) table on CAMB's 10000-node grid (zmax = 1e4, zmin = 0), with identical cosmology constants;
- the same native-provenance preliminary history and initialization.

The Julia run uses the caller-supplied Rodas5P callback reltol 1e-12, a1 1e-18, aex 1e-14. This is the test/reference configuration, **not** a library default; the library takes a required solve callback.

**CAMB injection.** Both complete X_e(z), T_m(z) arrays are injected into the **same** CAMB 2.0.4 build through a minimal diagnostic hook, so everything downstream is identical:
- background, reionization (tanh with `set_tau`), and CAMB's own cubic spline of the 10000 nodes;
- Boltzmann hierarchy, accuracy settings and lensing.

**Outputs:** ℓ = 2–3000 for unlensed TT/EE/TE, lensed TT/EE/BB/TE and φφ.

**Error definitions** (`A/chunk14/cl_compare.py`):
- **max frac:** max\|ΔC\|/\|C\| for TT, EE, BB and φφ. For TE it is max\|ΔC_TE\|/√(C_TT C_EE) (correlation-normalized, **not** ΔTE/TE, since TE crosses zero).
- **max abs:** max\|ΔC\|, in μK².
- **Cosmic variance:** max\|ΔC\|/σ_CV, for full-sky cosmic variance only (no noise, beam = 1, f_sky = 1, every ℓ). σ = √(2/(2ℓ+1))\|C\| for TT, EE, BB, φφ, and σ_TE = √((C_TE² + C_TT C_EE)/(2ℓ+1)), which is finite at TE zeros.
- **Diagonal sum:** Σ_ℓ (ΔC/σ)². This is a **per-spectrum diagonal diagnostic, not a joint TT/TE/EE or lensing likelihood χ²**: covariances between spectra and the non-Gaussian lensed-BB covariance are ignored.

## 2. Infrastructure (reused, diagnostic additions recorded)

**Existing:**
- CAMB 2.0.4 at `fa3f097343fbbe427cc04b4f5f0041c22c6ec764`, built with `RECOMBINATION_FILES="recfast cosmorec"`, under `tools/CAMB-cosmorec`; env `tools/camb-cosmorec-env`.
- `tools/cosmorec_validation/compare_recfast_cosmorec.py` compares *different models* (Recfast vs CosmoRec), so it is **not** a native-vs-Julia test. Only its high-accuracy settings and its TE variance formula were reused.

**CAMB scratch copy (`A/chunk14/camb_scratch`).**
- It is an rsync of the original tree without `.git`; all 394 files are hash-identical (`hashes_camb_original_tree.txt`, `hashes_camb_scratch_before_patch.txt`).
- The diagnostic patch (`camb_diag.diff`, 37 lines, `fortran/cosmorec.f90`) acts after the native CosmoRec call:
  - `CAMB_COSMOREC_DUMP=<file>` writes z, H, X_e, T_m;
  - `CAMB_COSMOREC_INJECT=<file>` replaces X_e/T_m from a file whose z column must match CAMB's grid to \|Δz\| ≤ 1e-9 (else stop), and prints max\|Δz\|.
  - Both are inactive unless the variables are set.
- **Build:** `make python RECOMBINATION_FILES="recfast cosmorec" COSMOREC_PATH=A/chunk14/cosmorec_lib/ GSL_LINK="-L$ENV/lib -Wl,-rpath,$ENV/lib -lgsl -lcblas"`, with gfortran 13.3.0, `-O3 -march=native -fopenmp`. Only `cosmorec.o` was recompiled; the other objects are the original build's.
  - `cosmorec_lib` holds the bitwise-validated chunk11 scratch `libCosmoRec.a` behind a no-op `lib` Makefile, so CAMB's sub-make cannot touch any CosmoRec tree.
  - The NEEDED libraries and RUNPATH are identical to the original `.so`.
- Hashes are in `hashes_build.txt`: diff 2684aa74…, scratch `camblib.so` 0f16667b…, `libCosmoRec.a` c3c0ff57…, original `camblib.so` e70a6976….

**Native parameterized harnesses (`A/chunk14/harness/`).**
- `harness4b_param.cpp` and `harness10a_param.cpp` are copies of the chunk4/chunk10 captures. The cosmology comes from `COSMO_H0`, `COSMO_OMBH2`, `COSMO_OMCH2`, `COSMO_YHE`, `COSMO_TCMB`, `COSMO_NNU`, with Ω = ω/(H0/100)² as in CAMB `fortran/cosmorec.f90`. The diffs are 15 and 9 lines.
- They are built against the original `libCosmoRec.a`, read-only.
- **Validation at the fiducial:** the probe, the preliminary history and the full 10a probe (final X_e/T_b, every pass and PDE stage) are **bitwise identical** to the existing captures (`harness/val/`).

**Other drivers:**
- `gen_H.py` uses the same code path as the fixture's H generator (original CAMB `get_background`); at the fiducial its data lines are bitwise identical to the fixture.
- `julia_history_cosmo.jl` builds the Cosmos accessors from a cosmology's native 4b probe, history and H table, with θ = [F, A2s1s, Y_p, Ω_b, h100, T0] taken from the native CFG read-back. At the fiducial it reproduces the existing Julia history **bitwise**. It asserts the native configuration (nShells 3, nS_eff 500, nShellsHeI 2, nS_eff_HeI 30, spin-forbidden on, HeI feedback off, diffusion on) before comparing.
- `run_cosmo.sh` (used) and `run_cosmo_v2.sh` (aborting version for future runs) drive one cosmology.
- `accept_check.sh` labels a run RUN_COMPLETE_AND_VALID only if every step has rc 0, all vectors are complete and the native replica is bitwise. This is **operational data validity, not a numerical-criteria verdict**; all five points satisfy it.

**CAMB settings** (`camb_forward.py`, from the existing validation script):
- kmax 10, k_per_logint 130, lens_potential_accuracy 8, lAccuracyBoost 1.2, min_l_logl_sampling 6000, DoLateRadTruncation False, NonLinear_both with halofit mead2020;
- output lmax 3000, computed lmax 5000; verified in effect (max_l 5200, max_eta_k 1.44e5).
- CosmoRec settings are asserted to be the defaults (runmode 0, accuracy 0, every override −1), matching the fixture runpars.
- Parameters are exactly those of the fixture's H generator, with Y_He set **explicitly**: the old validation script's 0.2454 does not match the fixture.

**CAMB interpolation boundaries** (`fortran/cosmorec.f90:105-138`):
- X_e is held constant beyond the grid ends (z ≥ 1e4 uses the z = 1e4 value; z ≤ 0 uses the z = 0 value);
- T_m beyond z = 1e4 scales as (1+z); inside, a regular cubic spline is used;
- CAMB adds reionization on top and solves z_re for τ = 0.0568.

**Runtime.** CAMB runs take about 3 s each (OpenMP on 20 cores); the Julia history takes about 70 s, including compilation.

## 3. Fiducial cosmology (`A/chunk14/fid/`, `cl_compare.log`)

Parameters: H0 67.36, ω_b 0.02237, ω_c 0.12, m_ν 0.06 (one massive), Ω_k 0, T_CMB 2.7255, N_eff 3.046, Y_He 0.24568275240335613, ln(10¹⁰A_s) 3.044, n_s 0.9649, τ 0.0568, ppf w = −1. Derived: z* = 1089.947, 100θ* = 1.04110.

| comparison | max frac uTT / uEE / uTE* | lTT / lEE / lBB / lTE* / φφ | max \|ΔC\|/σ_CV | CAMB x_e max rel (500–1600) |
|---|---|---|---|---|
| repeat run (scratch) | **bitwise identical** | bitwise | 0 | 0 |
| hook identity (CAMB's own native arrays injected) | **bitwise identical** | bitwise | 0 | 0 |
| scratch build vs ORIGINAL CAMB (both internal native) | 1.15e-7 / 2.87e-7 / 1.68e-7 | 9.0e-8 / 2.6e-7 / 1.7e-7 / 1.5e-7 / 1.9e-9 | 1.40e-5 | 1.16e-6 |
| **native floor:** fixture native vs CAMB-internal native | 1.40e-7 / 2.04e-7 / 1.20e-7 | 8.2e-8 / 1.2e-7 / 2.0e-8 / 7.2e-8 / 1.7e-9 | 1.12e-5 | 9.1e-7 |
| **port error:** Julia vs native (identical inputs, identical CAMB path) | **4.75e-7 / 6.15e-7 / 3.12e-7** | 4.75e-7 / 5.37e-7 / 1.79e-7 / 2.24e-7 / 3.37e-9 | **3.34e-5** (uEE, ℓ = 2952) | 1.45e-6 |

- Port-error maximum \|ΔC\|: 4.3e-4 μK² in TT (at low ℓ, where C is large), 9.5e-11 μK² in EE. The largest diagonal sum is 5.6e-7 (uEE).
- **Why a native floor exists.** CAMB's own grid comes out of the `-march=native` build ulp-shifted: 3468 of 10000 nodes differ by ≤ 1.8e-12, and the last node is z = 8.3e-13 instead of 0. CAMB's H also differs from the `get_background().h_of_z` table by ≤ 3.9e-13.
  - The native responds at about 1e-6 in X_e for z ≥ 50, as expected from its known input sensitivity.
  - At z ≈ 0 it responds up to 2.3e-3 (2e-5 at z = 1, 2e-6 at z = 2): the native low-z H-table quirk. CAMB's reionization dominates x_e there, and the effect on C_ℓ is in the floor row above.
- That is why the port is judged with **both** histories computed on the identical (Python-generated) grid/H and injected into one CAMB. The floor row is the native's own reproducibility through CAMB.

## 4. Nearby cosmologies (one parameter changed each; same pipeline; all RUN_COMPLETE_AND_VALID)

History criteria (approved: X_e ≤ 1e-5, T_m < 1e-7, max pointwise over all 10000 nodes), from `parity.log`. C_ℓ from `cl_compare.log`; summary in `A/chunk14/summary_table.md`.

| point | change | X_e max pw | T_m max pw | port: max frac uTT / uEE / uTE* | port: lTT / lEE / lBB / lTE* / φφ | port max \|ΔC\|/σ_CV | port max diag sum | native floor uTT / uEE / uTE* (max σ_CV) |
|---|---|---|---|---|---|---|---|---|
| fid | — | 1.449e-6 PASS | 3.557e-8 PASS | 4.75e-7 / 6.15e-7 / 3.12e-7 | 4.75e-7 / 5.37e-7 / 1.79e-7 / 2.24e-7 / 3.4e-9 | 3.34e-5 | 5.6e-7 | 1.40e-7 / 2.04e-7 / 1.20e-7 (1.1e-5) |
| ob | ω_b 0.02337 | 2.881e-6 PASS | **1.358e-7 FAIL** (criterion < 1e-7) | 5.15e-7 / 4.86e-7 / 2.45e-7 | 5.15e-7 / 4.66e-7 / 3.23e-7 / 2.33e-7 / 5.0e-9 | 1.83e-5 | 2.6e-7 | 9.53e-7 / 1.02e-6 / 5.75e-7 (5.5e-5) |
| oc | ω_c 0.125 | 1.773e-6 PASS | 9.790e-8 PASS | 9.79e-7 / 1.27e-6 / 5.66e-7 | 8.87e-7 / 1.08e-6 / 4.59e-7 / 4.20e-7 / 8.2e-9 | 6.71e-5 | 2.4e-6 | 8.54e-7 / 1.00e-6 / 5.53e-7 (5.4e-5) |
| h70 | H0 70.0 | 1.814e-6 PASS | 5.219e-8 PASS | 7.10e-7 / 8.57e-7 / 4.35e-7 | 6.65e-7 / 7.02e-7 / 2.65e-7 / 3.57e-7 / 2.4e-9 | 4.63e-5 | 7.1e-7 | 5.23e-7 / 6.01e-7 / 2.88e-7 (3.2e-5) |
| yhe | Y_He 0.2400 | 2.129e-6 PASS | 8.649e-08 PASS | 5.06e-7 / 7.11e-7 / 3.11e-7 | 3.50e-7 / 4.30e-7 / 1.56e-7 / 2.49e-7 / 2.8e-9 | 3.89e-5 | 3.3e-7 | 6.12e-7 / 7.47e-7 / 2.88e-7 (4.1e-5) |

At every point the Julia and native helium-switch nodes agree (fid 1342, ob 1338, oc 1343, h70 1342, yhe 1340), and the native replica is bitwise equal.

**Other internal parity** (reported separately; 10a metrics; unchanged gates; not acceptance):

| point | pass-1 X1s / 2s / 2p units (gates < 1, < 10, < 10) | PDE feedback max (gate ≤ 1e-3) |
|---|---|---|
| fid | 5.21 / 16.0 / 16.2 | 5.1e-5 |
| ob | 3.85 / 16.9 / 18.4 | 2.2e-5 |
| oc | 8.05 / 298.9 / 2065.7 | 5.7e-4 |
| h70 | 3.39 / 14.1 / 15.4 | 4.8e-5 |
| yhe | 4.70 / 7.6 / 18.0 | 2.3e-5 |

- `oc`'s 2s/2p maxima are a single node at z = 500.517, the last node above the feedback cut at z = 500, where the diffusion correction switches off discontinuously. Away from \|z − 500.5\| ≤ 2 they are 9.6 / 9.8 units (`cosmo_oc/rows/pass1_outlier.log`). It was not investigated further, by instruction (no new research).
- Pass-1 X1s exceeds its gate at every point, as at the fiducial. The original fails this gate against itself (`docs/ORIGINAL_VS_PORT_ACCURACY_INVESTIGATION.md` §5.9).

## 5. What is demonstrated, and what is not

**Numerical criteria (unchanged):**
- X_e ≤ 1e-5: **PASS at all 5 points** (1.45e-6 to 2.88e-6).
- T_m < 1e-7: **PASS at 4 of 5; FAIL at `ob`** (1.358e-7). The threshold was not changed, and no published T_m accuracy claim exists (`docs/COSMOREC_PUBLISHED_ACCURACY_REVIEW.md`).

**Boltzmann-facing practical result.** Over these 5 points, the Julia-vs-native difference in CAMB 2.0.4 C_ℓ (unlensed TT/EE/TE, lensed TT/EE/BB/TE, φφ; ℓ = 2–3000) is:
- max fractional ≤ 1.27e-6 (oc uEE), with TE correlation-normalized ≤ 5.7e-7 and φφ ≤ 8.2e-9;
- max ≤ 6.7e-5 σ_CV at any single ℓ, with every diagonal sum ≤ 2.4e-6.

This is the same order as the native's own reproducibility through CAMB under ulp-level input differences (0.1–1e-6), and about 800× below the original authors' 0.1% scale for C_ℓ agreement between CosmoRec settings (Shaw & Chluba 2011). At this practical C_ℓ level, the `ob` T_m miss produces no visible effect: its C_ℓ port error is ≤ 5.2e-7. That is a practical statement, **not** a pass of the T_m criterion.

**Demonstrated:** for the default CosmoRec configuration, in CAMB 2.0.4 at the stated accuracy settings, the port reproduces the original's Boltzmann-facing output at 5 cosmologies within ~1e-6 fractional C_ℓ.

**NOT demonstrated:**
- Whole cosmological initialization in Julia (B3): every point used native-provenance preliminary history and constants, captured by the parameterized native harness. The Julia package still cannot initialize an arbitrary cosmology by itself.
- Coverage beyond one-parameter offsets around Planck-like values.
- Other Boltzmann codes (CLASS), other CAMB accuracy settings, other CosmoRec settings ('full'/high-accuracy are not ported).
- Joint likelihood or parameter-level bias: not computed, though the differences are ≪ cosmic variance per ℓ.
- Gradients and true reverse mode (B2), and the 6 remaining test failures. These are separate.
- The physics accuracy of CosmoRec itself (model claim ~0.1%).

## 6. Commands and terminal states

All commands exited with rc 0; per-run logs, `status.txt` and meta JSON are saved.

| step | command / script | output |
|---|---|---|
| CAMB scratch build | as in §2 | `build_camb_scratch*.log`; first attempt with `-lopenblas` FAILED (preserved): fixed with `-lcblas` |
| parameterized harnesses | `harness/build{4b,10a}_param.sh` | validated bitwise at the fiducial |
| fiducial | `camb_forward.py` (orig/scratch; native/dump/inject), `julia_history_fid.jl`, `cl_compare.py` | `fid/` (`runs_status.txt`, `*_cls.txt`, `*_xe.txt`, `*_meta.json`, `inj_*.txt`, `native_dump.txt`) |
| per-cosmology validation | `julia_history_cosmo.jl` at the fiducial | `fidval/` (bitwise vs `fid/`) |
| cosmologies | `run_cosmo.sh ob 67.36 0.02337 0.12 0.24568275240335613`; `oc … 0.125 …`; `h70 70.0 …`; `yhe … 0.24` | `cosmo_<name>/` (`status.txt` all rc 0, `parity.log`, `cl_compare.log`, `compare.jsonl`) |
| summary | `summarize.py` | `summary_table.md` |
