# Native CosmoRec fixture milestone: accepted, limited to public outputs

## Reference producer

The plain-text fixture in the working tree, `test/fixtures/native_camb_cosmorec_thermo.txt`, was captured with the EXISTING checked-out, CosmoRec-linked CAMB 2.0.4 Python environment and public thermodynamics API, not a RECFAST substitute, copied atomic tables, locally rebuilt C++ fork, or altered source. Runtime asserts check the CAMB source/build, exact original CAMB/CosmoRec revisions, the adapter is using `CosmoRec`, CAMB batch fields are native defaults/unset, and H(0)=H0. Its production `camb_worker._build_params` is imported from the existing Capse pipeline. YHe is explicitly labeled BBN-consistent; it is not hardcoded to 0.2454. H, x_e, and T_b units/conventions and Hubble spot/independent-density checks are recorded in the fixture header.

There are **245 ascending redshift rows** across grid nodes and independent holdout queries, from z=0 to z=10,000. The file records CAMB's input-grid contract (10,000 uniform redshifts supplied to CosmoRec) separately from this output/query grid. The compiled native Hubble endpoint behavior is left unchanged and labeled as unpatched. Both reionized CAMB outputs and a no-reionization recombination-only history are captured. It contains final public outputs only, not CosmoRec level populations or PDE intermediates.

Two separate generator processes produced byte-identical fixture files, identical metadata and matching SHA256; logs and the second-run artifact are retained under `cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/native_fixture/`. Fixture raw rows are 17-digit-style round-trip text. The Julia test validates original SHAs, model/provenance metadata, row hash, shape, units, finiteness, redshift coverage, physical spot rows and interpolation against native holdouts.

## Independent checks and thermal follow-up

Two separate original-CAMB generator processes produced byte-identical complete fixtures with identical metadata and SHA256 `394bf61ff8c80c567d8575cfef57d72df0825ae8d6686f4696da0945a6111c96` (0 differing fields). Focused fixture validation passed **297/297**.

Full package validation also completed under a PTY: `julia --project=/home/marcobonici/Desktop/work/CosmologicalEmulators/CosmoRec.jl -e 'using Pkg; Pkg.test()'` → **900/900**, exit0, 20m35.7s. Log `local_analysis/cosmorec_differentiability_20260930/chunk1d/pkgtest.log` records the exact command, result and terminal exit status. This includes 603 existing tests and 297 fixture tests.

Later supervisor test attempts were deliberately stopped when the laptop reached its 100 °C CPU limit while four unrelated approved CMB-lite Julia chains were running. Only our CosmoRec test process was signaled; all four chains and native sources were left untouched. Interruption logs remain under `local_analysis/cosmorec_differentiability_20260930/native_fixture/`. Future local test suites need thermal/load headroom. The earlier completed PTY run is the valid full-suite evidence.

The supervisor separately read the fixture and generator, verified that the file exists and checked reported repeated-capture/hash evidence. This milestone establishes only a public xe/T_b/H reference target. It does NOT provide rate-table, per-level ODE, radiation operator, DPesc/correction-integral, or low-z matching fixtures. No physics component is yet represented in `src/`.

**Accepted gate:** reproducible, source-pinned, plain-text native public-thermodynamics fixture, unconditional focused/full-package tests, and provenance. **Not accepted or claimed:** validation of internal native intermediates or a Julia/native physics port.

Next, capture component-level native text references before each translation. In particular, the rate interpolant and state RHS still need their own oracle data; public CAMB history alone cannot validate them.
