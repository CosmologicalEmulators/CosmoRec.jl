# Chunk 4 acceptance summary (phase 4: initialization slice)

All four sub-chunks have native text fixtures with tested hashes, unconditional focused tests (native and AD) and a fresh full root `Pkg.test()` each:
- 4a: 765 + 280 focused; full 23701/23701 (`chunk4a/pkgtest.log`).
- 4b: 14097 + 16; full 37814/37814 (`chunk4b/pkgtest.log`).
- 4c: 8475 + 29; full 46318/46318 (`chunk4c/pkgtest.log`).
- 4d: 984 + 32; full 47334/47334, exit 0 (`chunk4d/pkgtest.log`).
Per-stage acceptance: `ACCEPTANCE_CHUNK4A.md` ... `ACCEPTANCE_CHUNK4D.md`; results and worst errors: `CHUNK4_RESULTS.md` and `CHUNK4A..D_RESULTS.md`; design: `CHUNK4_DESIGN.md`.

Accepted: explicit-input Saha initialization and packing (4a); GSL-compatible splines and Cosmos accessors from the native 6000-node history and Hubble table (4b); Recfast++ production history with node grid, closed-form Saha segments and a caller-supplied stiff solve (4c), validated against the native RHS (900 probes) and the native 6000-node history to the native solver's own tolerance; the discrete sampled-HeI switch (4d).
Not accepted / not claimed: any full CosmoRec ODE, diffusion/feedback, history or spectrum parity, differentiability across discrete branches (Saha-segment grid, z_saha/zsRe switches, helium switch), event-time sensitivities. Phase 5 has not started. Nothing is staged, committed or pushed.
