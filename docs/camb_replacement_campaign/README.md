# CAMB replacement campaign evidence (copied 2026-10-11)

Verbatim copies of the certification records. Canonical live location:
cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk36_camb_replacement/bridge/ (untracked scratch tree;
these copies exist so the evidence survives independently of it).

- IMPLEMENTATION_FREEZE.md - frozen intended implementation (callbacks, tolerances, hashes, what is NOT required)
- CAMPAIGN_RESULTS_3POINT.md - 3-point full-lens gate results (all PASS, worst 7.8e-6 = 64x under the 5e-4 budget)
- gate.jsonl - machine records of every gate evaluation
- run_header.txt - library hashes and tool versions for the campaign

Full 9501-row spectra and per-run meta live in the scratch tree (regenerable deterministically via
bridge/cmb_pair/run_cmb_pair.sh; do not commit multi-MB spectra).
