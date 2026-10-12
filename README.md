# CosmoRec.jl

A Julia translation of CosmoRec for computing cosmological recombination histories.

## Solver configurations and examples

CosmoRec.jl requires callers to supply the ODE solve callbacks; it does not
select solver defaults. The configurations below are examples for a caller-owned
CAMB adapter, not exported package presets. The callback implementations shown
by name are integration-provided functions, not CosmoRec.jl exports. In
particular, the optimized callbacks used in the campaign are development code,
not a stable installed callback API.

```julia
using CosmoRec

# Implement these callbacks in your CAMB integration. They must satisfy the
# solve_recfast / solve_pass / solve_tail callback signatures used by the adapter.
function make_callbacks(; main_rtol, a1, aex, preliminary_rtol, tail_rtol)
    solve_recfast = (args...) -> solve_rodas_dense(
        args...; alg=Rodas5P(), reltol=preliminary_rtol,
        abstol=[1e-16, 1e-16, 1e-12],
    )
    solve_pass = (args...) -> solve_phase5(
        args...; alg=Rodas4P(), reltol=main_rtol, a1, aex, tstops=true,
    )
    solve_tail = (args...) -> solve_tail5(
        args...; alg=Rodas5P(), reltol=tail_rtol,
        abstol=[1e-16, 1e-16, 1e-12],
    )
    return (; solve_recfast, solve_pass, solve_tail)
end

# Use CAMB's own sampled background, cosmological inputs, and explicit original
# runmode-0 rate tables. Do not substitute an approximate H(z) for CAMB's table.
tables = CosmoRec.load_runmode0_tables(tables_root)
cbs = make_callbacks(; main_rtol=1e-10, a1=1e-18, aex=1e-12,
                      preliminary_rtol=1e-12, tail_rtol=1e-12)
history = CosmoRec.cosmorec_history(
    ; ombh2, omch2, omk, N_eff, H0, T_CMB, Y_He,
      z_Hz, Hz, z_out, tables, cbs...,
)
# history.Xe is n_e/n_H; history.Tm is in kelvin, on z_out.
```

The snippet shows where to bind solver settings in an existing integration; it
is not standalone until that integration supplies `solve_rodas_dense`,
`solve_phase5`, `solve_tail5`, CAMB inputs, and a valid table root. The test
helpers are not a supported callback API and may load test fixtures.

### Safe: adopted callback settings

```julia
cbs = make_callbacks(; main_rtol=1e-10, a1=1e-18, aex=1e-12,
                      preliminary_rtol=1e-12, tail_rtol=1e-12)
```

### Aggressive: fastest measured setup passing the sampled forward gate

```julia
cbs = make_callbacks(; main_rtol=1e-3, a1=1e-10, aex=1e-8,
                      preliminary_rtol=1e-8, tail_rtol=1e-8)
```

Both retain `tstops=znodes` in the main passes. Disabling node stops made loose
configurations unstable. The aggressive setup is a forward-spectrum
optimization, not a new package default and not a validated gradient/NUTS
configuration. Validate derivatives separately before gradient-based inference.

## Forward CMB validation and Genoa performance

The comparison replaced only the recombination provider: CAMB supplied the
same inputs to native CosmoRec and CosmoRec.jl. The gate checks pointwise
fractional differences in lensed TT and EE at every integer multipole
`2:9500`, requiring each maximum to be below `5e-4` (0.05%).

| Callback setup | Prior draws passing | Worst lensed TT | Worst lensed EE |
|---|---:|---:|---:|
| Safe | 100/100 | `4.178e-6` | `2.293e-6` |
| Aggressive | 100/100 | `5.202e-6` | `3.021e-6` |

Draws covered the approved domains (33 YHe, 33 degenerate, 34 A_L). This
certifies the tested runmode-0 forward CAMB setup on those draws, not every
possible cosmology, mode, or gradient.

Recombination-only timing used one bound core of a Narval AMD EPYC 9534,
identical captured fiducial CAMB inputs, and 50 warm BenchmarkTools samples per
implementation. CAMB execution, table loading, and Julia compilation were
outside the timed region.

| Configuration | Julia median | Julia p95 | Paired native C++ median |
|---|---:|---:|---:|
| Safe (`main_rtol=1e-10`) | `1.306 s` | `1.719 s` | `0.572 s` |
| Aggressive (settings above) | `0.394 s` | `0.628 s` | `0.565 s` |

The aggressive median was about 30% below its paired native median in that
benchmark. Its p95 exceeded native, so this is not a claim of uniformly lower
latency. These are hardware- and job-specific measurements, not universal
speed guarantees. Julia BenchmarkTools allocation counts do not measure native
C++ allocations.

### Earlier single-point measurements

An earlier laptop benchmark used different callbacks and one cosmology, with
spectra through `ell=10000`. Those historical `rtol=1e-5` measurements are
not the configurations or domain-wide results above, and their timings are
not directly comparable to the Genoa run.
