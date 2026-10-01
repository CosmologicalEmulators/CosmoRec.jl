# Chunk 3a acceptance: fixture-backed smooth hydrogen RHS components

The supervisor independently ran the two focused component suites on the persistent Chunk2 environment:

- `test/chunk3a_rhs_native.jl`: **272/272** assertions, native helper outputs bit-identical except one table path at2.4e-16.
- `test/chunk3a_rhs_ad.jl`: **323/323**, exit0; individual and combined ForwardDiff Jacobians, directional checks and prepared Mooncake VJPs match the chosen finite-difference/ForwardDiff references. Worst reported VJP errors are below1e-12.

The native text fixture `test/fixtures/native_hrhs_components.txt` comes from a harness outside all Git repositories calling the **original compiled CosmoRec C++ functions** with nonzero input-buffer sentinels. Two native processes were byte-identical. The rate fixture is the already accepted explicit AtomicRateTable input. Original native source/static libraries are unmodified. Fixture, sources/line map, citations and all errors/domain cases: `docs/CHUNK3A_RESULTS.md`. Source production code starts at `src/RHS.jl`; unit tests are unconditional.

The full root suite was also reported by the worker as **6564/6564**, `docs/CHUNK3A_RESULTS.md`, `chunk3a/final_full_pkgtest.log`. It used that task's owned pause/resume watchdog on its own Pkg.test process; no unrelated process was touched. The subsequent rate additions/full suite is 5969 tests plus 595 new Chunk3a assertions = 6564 total.

## Narrow scope: not the full hydrogen RHS

The fixture and implementation cover:

- Compton/adiabatic matter-temperature-ratio term.
- Base 2s-1s two-photon term.
- Lyman n=2,3 Sobolev channels.
- Effective hydrogen continuum Rci/Ric and resolved interlevel Rij.
- Their test harness's combined H helper, off diffusion/quadrupole.
- Ground-state electron/proton/He abundance algebra as a separate utility.

**Do not claim this is the complete default `fcn_effective` or complete default ODE.** The fixture harness uses an H-only state (`flag_He=0`) and does not initialize/fixture the helium-level populations, HI-absorption Pesc interpolation tables, or the coupled helium corrections. Default CAMB mode remaps `HI_absorption=2` to active `1`; its helium-photon-to-hydrogen feedback is a distinct active dependency below `zcrit_HI=3400`. Although the native exported H helper is directly evaluated, this harness does not establish the fixture/provenance of that He-data-dependent absorber for the configured full batch model. Treat it as explicitly excluded until captured and separately ported/tested. Likewise excluded: helium RHS, diffusion/DI1/DI2, quadrupole channels, rate rescaling/`f_t`, Saha/Recfast initialization, helium sampled switch, low-z matching, backgrounds, the ODE solve, all PDE corrections.

Known native table overflow below Tg≈55.6K yields NaNs and Mooncake gradient contamination; the probe documents this rather than hiding it. The native callback's reached temperature domain is not established here. The fixture's test-defined cosmology map is not a production background model. No ODE solve or gradient through solve was tested in this chunk.

**Accepted scope:** these source-backed hydrogen/thermal component functions and their tests. **Not accepted:** a recombination trajectory, complete default ODE, or physics solver product. Next must source-capture the missing H/He coupling/integer switch/initialization data before the connected population ODE can be claimed complete.
