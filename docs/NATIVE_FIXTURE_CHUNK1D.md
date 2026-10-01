# Chunk 1d: native output fixture (provenance)

**Scope: anchors FINAL OUTPUTS ONLY** (x_e, T_b, H after CAMB's own resampling). It does not replace separate
original-reference fixtures for effective-rate tables, ODE RHS, PDE operator fields, correction integrals or low-z
matching. None of those internal arrays are exposed by the public CAMB API and none are captured here.

## Files
- `scripts/generate_native_fixture.py` (`generate OUT` / `compare A B`); run inside the existing activated env.
- `test/fixtures/native_camb_cosmorec_thermo.txt` (plain text, 245 rows, 29147 bytes, full double precision).
- `test/chunk1d_native_fixture.jl` (included from `test/runtests.jl`; Base + Test + SHA stdlib; no Python).

## Provenance (asserted by the generator before writing)
- CosmoRec `086769055f61ae0c244a53dd381ee65b624d0ac3`, CAMB-cosmorec `fa3f097343fbbe427cc04b4f5f0041c22c6ec764`, both clean.
- `camb.__file__` resolves inside the checked-out CAMB-cosmorec; `type(pars.Recomb).__name__ == "CosmoRec"`.
- camb version, camblib.so sha256, python/numpy versions, camb_worker.py sha256 are in the fixture header.

## Settings
- Parameters from production `camb_worker._build_params` (fiducial: ln10As 3.044, ns 0.9649, tau 0.0568, H0 67.36,
  omega_b 0.02237, omega_c 0.12, Mnu 0.06, w0 -1, wa 0), YHe from BBN (0.24568275240335613).
- Implicit adapter defaults (not overridden, asserted equal to defaults): runmode 0, fdm 0, accuracy 0, and all -1 overrides
  (diffusion iterations, shells, HeI shells, HI absorption, 2-photon, Raman -> CosmoRec preset defaults).
  Explicitly tested: only the fields above and the CosmoRec model class.
- The adapter runs CosmoRec on 10000 uniform z (1e4 -> 0) with CAMB's H(z); that grid is distinct from the fixture query grid.
  Native endpoint behaviour is recorded as-is; not patched, no impact claim.

## Query grid
z ascending, 0..10000, 245 rows: 216 `grid` nodes (dense low z 0.25/0.5/2.5 steps, 10 steps to 200, 25 steps over H recombination
550-1500, 50/100 steps over He regions to 7500, 250 to 10000; includes 50, 200, 3000, 9999, 10000) and 29 `holdout` rows at
off-node z. All rows come from ONE `get_background_redshift_evolution` call because CAMB's thermo table depends on the first
query's min time (z=0 present here). Hence values differ ~1e-9 relative from the supervisor smoke values (queried at 6 points).

## Columns / units
`z role x_e T_b H x_e_noreion T_b_noreion`. T_b in K, H in km/s/Mpc (checked against c*sqrt(8piG a^4 rho/3)/a^2 from
`get_background_densities`, relative agreement in header; H(0)=H0 to 1e-12). `x_e`,`T_b` include reionization (tau 0.0568);
`*_noreion` use `Reion.Reionization=False`. x_e includes helium (about 1.164 when fully ionized).

## Determinism
Two fresh-process runs: 0 differing fields, files byte-identical, same sha256 (below). No tolerance applied.
- fixture file sha256: `394bf61ff8c80c567d8575cfef57d72df0825ae8d6686f4696da0945a6111c96`
- rows sha256: `908d338239dcf8a21cc7d8837919b925b03ebd846c4e5111d908ad824b241939`
- run artifacts and logs: `local_analysis/cosmorec_differentiability_20260930/chunk1d/`.

## Interpolation check (tolerances calibrated from the holdouts, about 2-3x margin)
4-point Lagrange (log(1+z) for H, T_b; z for x_e no-reion) from `grid` rows vs `holdout` native queries:
max errors H 3.1e-4, T_b 1.5e-3, x_e 8.0e-4 (relative, x_e floored at 1e-3). Thresholds 1e-3 / 5e-3 / 2.5e-3.
Worst rows: H, T_b at z=0.05 (near the endpoint); x_e at z=450.

## Test results
- Focused `include("test/chunk1d_native_fixture.jl")`: 297/297.
- Full `julia --project=<repo> -e 'using Pkg; Pkg.test()'`: **900/900** = 61 + 434 + 108 + 297, exit 0, 20m35.7s (log `pkgtest.log`).
- Three earlier full runs launched without a pseudo-terminal (stdin-less background shell) were SIGTERM'd about 2 min in, inside
  chunk 1a, before reaching the new test; they are kept as `pkgtest_killed_attempt{2,3}.log` and `pkgtest` first attempt. A focused chunk 1a
  run passed 61/61, and the full run under `script` (pty) passed. Root cause of the kill was not identified; memory was ample.
- Test dependency added: stdlib `SHA` in `[extras]`/`[targets]` of `Project.toml` (hash check only).

## Captured vs unavailable
Captured: final x_e, T_b, H (with and without reionization). Unavailable via public API (need separate reference fixtures):
effective-rate tables, level populations, ODE RHS, PDE operator fields, correction integrals, low-z matching, the raw 10000-point
CosmoRec output before CAMB resampling, and the native endpoint H zero-fill effect.
