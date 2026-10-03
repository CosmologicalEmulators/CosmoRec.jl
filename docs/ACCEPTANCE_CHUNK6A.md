# Chunk 6a acceptance: HI diffusion-PDE coefficients
## Gates
- `test/chunk6a_pde_coefficients_native.jl` 825/825; `test/chunk6a_pde_coefficients_ad.jl` 6/6; full root `Pkg.test()` **50275/50275**, exit 0, `chunk6/pkgtest6a.log`.
- Native oracle with the loaded CAMB H re-armed as in production; the analytic-H capture is preserved as rejected evidence and was not used.
- Dnem judged on the source-derived residual scale (2.2e-16 from native populations; end-to-end equal to the solver-level population difference); raw and absolute errors are reported too.
## Accepted scope
Population splines, get_rates_all, pd and Rp/Rm, pd/Dnem coefficient splines for the production configuration.
## Limitations carried forward
Piecewise (spline cells, grid count, stencils); the Phase 5 gradients still freeze the Saha initial state and the preliminary history; one cosmology; required local data path (D1).
