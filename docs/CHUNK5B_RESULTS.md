# Chunk 5b results: one recombination pass (12-state ODE, sampled-HeI switch, 7-state ODE)

NOTICE: port of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al.; cite Chluba & Thomas 2010, Chluba & Sunyaev 2006, Rubino-Martin et al. 2008. No native source, library or table was changed or bundled; the tests require `COSMOREC_NATIVE_DATA_DIR` (decision D1, see `CHUNK5A_RESULTS.md`).

## Scope
`recombination_pass` (`src/RecombinationODE.jl`): the native `Xe_frac_effective_rates` (`main.CosmoRec.cpp:84-262`) without diffusion (iteration 0, native `runmode 1`): node grid `init_xarr_linear(3000, 50, 3000)` (accumulated like the C++ loop, last node 50.000000000188152), initial packed state from the 4a Saha initializer fed by the 4b accessors, the 4d decision `(z < 200 || fHe - He1s <= 1e-7)` tested at every node on primal values, the 4d reset at the first node that satisfies it, and a restart with the 7-state system; per node observables `[Xe, X1s, 2s, 2p, 3s, 3p, 3d, rho]` and `Te = rho TCMB`. The stiff solve is a caller-supplied callback (the package has no ODE dependency): the tests use `Rodas5P` with an explicit AD Jacobian and time gradient and error control on primal values (`test/chunk5_helpers.jl`), `reltol = 1e-12`, `abstol` 1e-16 (rho, X1s, He 1s) and 1e-14 (excited levels, effectively unmonitored: the H/He excited states are quasi-steady slaves with Jacobian entries up to 3.5e17; tighter absolute tolerances on them make the error estimator noise-limited, `retcode = Unstable`).
Never written into a fixture: no Julia solution is used as a native reference.

## Native oracle
ONE original pass: `cosmorec_calc_h_cpp_` with `runmode = 1` and the CAMB-adapter arguments (thermo-fixture cosmology, the CAMB H(z) table of 4b), run by the external harness `chunk5/native_capture/harness5b.cpp` (SHA-256 `dca906d7a9520e31ef4685a031c5770b53493b68124729b35269cecbaacd5ca8`); afterwards the native globals are read (`pass_on_the_Solution_CosmoRec` 3000 x 11, `output_CosmoRec` 3199 x 3, `HeI_was_switched_off_at_z`, the returned `Xe_arr`, `Te_arr` on the 10000-node CAMB grid, and the final `Level_I.X`). The native pass takes 0.18 s; two fresh processes give byte-identical files. Fixtures (verbatim, hashes in the headers and tested): `native_ode_pass_nodes.txt` (SHA-256 `e3f06ae0...a62a`), `native_ode_pass_output.txt` (`2a96e12d...5`), `native_ode_pass_grid.txt` (`d9bf4d2a...a`), `native_ode_pass_summary.txt` (`533a10ca...2`). Note: the native pass stores no helium populations for `runmode 1` (`Diffusion_correction_HeI` off); the helium ground state is recovered through `Xe - Xp` and the switch redshift.

## Results (focused 42/42)
- Node grid: bitwise equal to the native `z` (including the accumulated last node). Initial state equals the native Saha node 0 to 1e-13.
- Helium switch: node 1342 and redshift **1680.9103034346122 equal to the native `HeI_was_switched_off_at_z` bitwise**; 12-state vectors on nodes 1..1342, 7-state after.
- Nodewise vs the native internal arrays, differences in NATIVE tolerance units (native monitors only rho, X1s with rel 1e-7/abs 1e-12 and 2s, 2p with rel 1e-6/abs 1e-80): X1s max 0.114 (z = 1679.9), rho 0.008 (z = 50), H 2s 0.50, H 2p 3.22 (z = 87.4); medians 2e-5 to 3e-9 units. Xe max relative 1.72e-7 (1.7 units of 1e-7; z = 916.6); Te relative 8.0e-10. Unmonitored levels (3s, 3p, 3d): max relative 5.9e-7, 6.4e-7, 5.1e-7 (z = 65.7). The native error is thus at its own tolerance, the Julia solution converges: the excited levels at tolerance 1e-8 differ by 5.8e-3 (2p), at 1e-10 by 8.7e-5, at 1e-12 by 3.2e-6, at 1e-13 by 1.0e-6 (kept in `chunk5/diag5b.jl` logs).
- Solver self-consistency (reltol 1e-12 vs 1e-13): Xe 2.2e-9, Te 1.4e-11, 2p 3.3e-6. Frozen switch branch vs the detected run: agrees to solver accuracy (not bitwise: the last 12-state block ends at the switch node).
- After the switch `Xe = xp` (tested); the switch node itself still carries `XHeII` (pre-switch state), as native.
- Gates: X1s, rho < 1 native unit; 2s, 2p < 10; 3s/3p/3d < 1e-3 relative; Xe < 5e-7 relative; Te < 1e-8.

## Limitations
The switch node is discrete (a 1e-3 relative change of H(z) moves it by a few nodes; see 5e); the solution beyond the native tolerance cannot be validated against the native (the native itself is only accurate to ~1e-6..1e-7 on monitored components and unmonitored on others); one cosmology.

## Benchmarks (`benchmark/chunk5b_benchmarks.jl`, `chunk5/bench5b.log`)
Full pass (3000 nodes, Rodas5P reltol 1e-12): median 758 ms (902 MiB, 14.1 M allocations; 7-state blocks restart the solver); frozen branch 774 ms; looser pass (reltol 1e-8) 78 ms; node grid 3.0 us. The first sample includes compilation (cold).

## Files
`src/RecombinationODE.jl`, `src/CosmoRec.jl`, `test/chunk5_helpers.jl`, `test/chunk5b_ode_pass_native.jl`, `test/runtests.jl`, the four `test/fixtures/native_ode_pass_*.txt`, `benchmark/chunk5b_benchmarks.jl`, this file, `docs/ACCEPTANCE_CHUNK5B.md`. No Git operation.
