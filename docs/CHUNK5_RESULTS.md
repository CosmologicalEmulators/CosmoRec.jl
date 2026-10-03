# Chunk 5 results (phase 5 complete): the first-pass recombination ODE (native iteration 0, runmode 1)

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al.; cite Chluba & Thomas 2010, Chluba & Sunyaev 2006, Rubino-Martin et al. 2008. Native source, libraries and tables untouched; tables are read from the explicit `COSMOREC_NATIVE_DATA_DIR` (never bundled; the tests fail loudly when it is unset). Nothing is committed. Design: `CHUNK5_DESIGN.md` (D1 resolved: explicit local path; D2: runmode 1 single pass as the oracle; D3: native-tolerance units).

| stage | content | fixtures | focused | full root `Pkg.test()` |
|---|---|---|---|---|
| 5a | 12-/7-state RHS dy/dz incl. the new `flag_He = 0` branch, full-table loader | `native_ode_rhs.txt` | 1966 + 35 | 49335/49335 |
| 5b | one pass: 12-state ODE, helium switch, 7-state ODE (SciML Rodas5P) | `native_ode_pass_{nodes,summary}.txt` | 42 | 49377/49377 |
| 5c | Recfast tail (`Xe_frac_rescaled`) | `native_ode_pass_{output,summary}.txt` | 17 | 49394/49394 |
| 5d | output assembly / returned history on the CAMB grid | `native_ode_pass_grid.txt` | 13 | 49407/49407 |
| 5e | ForwardDiff + prepared Mooncake gradients of the whole pass | (5b-5d) | 37 | 49444/49444 |

Key native comparisons: dy/dz 9.5e-15; helium switch at node 1342, z = 1680.9103034346122 bitwise equal to the native; nodewise within native-tolerance units (X1s 0.11, rho 0.008, 2p 3.2); Xe 1.7e-7; tail Xe 1.03 / Te 166 native units (native Recfast error); assembly exact (2.2e-16 from native rows), Julia history inside the stored range 1.6e-7 (Xe) / 3.0e-8 (Te). Gradients: FD elasticity 1.2e-5, Mooncake VJP vs ForwardDiff 1.3e-9..1.3e-7.
Failed/rejected approaches (logs in `chunk5/`): excited-state absolute tolerances below ~1e-14 make Rodas5P noise-limited (Unstable); column-scaled FD metric; continuous adjoints (unstable). Bug fixes: `RateResult` element types; `helium_switch_reset` aliasing (4d).
Limitations carried forward: native iteration 0 only (no diffusion: phase 6+); discrete switch/grid; frozen preliminary history in gradients; machine-specific data path; one cosmology.
