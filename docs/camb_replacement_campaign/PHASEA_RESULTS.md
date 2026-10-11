# Phase A: 36-point prior-domain certification (2026-10-11, commit 38bdf8a)

First campaign beyond the 3 certification points: 36 cosmologies spanning the three Capse.jl worker
domains exactly as `set_cosmology` receives them (PPF throughout, per the emulator-zoo workers):

- `yhe_*` (12): explicit-YHe / massive-neutrino variants (nnu 3.046, 1 massive, explicit YHe).
- `deg_*` (12): degenerate-hierarchy variants (BBN null YHe, nnu 3.044 / snn 3.044, 3 massive).
- `alen_*` (12): lensing-amplitude variants (mnu 0.06 fixed, degenerate hierarchy, A_L in [0.7, 1.4]).

Each family: 1 baseline + 6 interior prior draws (rejection-sampled, `w0+wa < -0.5` enforced) + 5 boundary
points (incl. the `w0+wa = -0.52` probe). Point definitions: `bridge/cmb_pair/POINTS_PHASEA.md`
(generator `gen_phaseA_points.py`, deterministic, seed 20261011).

## Result

**36/36 PASS** (lensed TT&EE, all ell 2..9500, budget 5e-4 max fractional difference vs native CAMB+CosmoRec).

Worst per family (max over the 12 points):

- yhe: lensed TT 2.563e-6 (`yhe_int5`), lensed EE 1.245e-5 at ell 9500 (`yhe_int3`).
- deg: lensed TT 2.125e-6 (`deg_bnd1`), lensed EE 7.255e-6 (`deg_bnd2`).
- alen: lensed TT 1.386e-6 (`alen_bnd3`), lensed EE 3.571e-6 (`alen_int6`).

The YHe domain is the most sensitive, as expected (YHe feeds the helium recombination history directly);
its worst EE is still 40x under budget. Low-ell TT maxima (~2e-6 at ell 2-3) are the absolute-scale floor:
tiny absolute differences against small C_l, not a recombination effect.

## Evidence

`bridge/cmb_pair/`: `gate.log` (per-spectrum lines for all 36), `gate.jsonl` (machine records),
`{point}_native/julia_cls.txt` (9501-row spectra), `{point}_*_meta.json` (params, camblib sha, derived YHe/N_eff),
`{point}_julia_bridge.txt` (call counts, per-call seconds), per-step `.exit` files (all 0), `ALL_DONE`.
Campaign wall: ~22 min serial on the laptop (36 native + 36 Julia routes + gates), doubles as a 36-call
bridge soak with no call-count or determinism anomalies.
