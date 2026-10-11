# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Chunk 35: standalone initialization (no native-provenance inputs) of the CAMB-adapter route of ORIGINAL CosmoRec v3.0b
# (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al.:
#   A1 cosmology closure       Development/Cosmology/Cosmos.cpp:207-250 (init), :286-295 (set_initials), :305-335 (calc_dependent_consts), :1141 (sigT),
#                              Cosmos.h:326-333 (calc_Orel), Development/Definitions/physical_consts.h:44-155 (CODATA-2008 cgs set; NOT the Recfast++ SI set),
#                              CosmoRec.cpp:766-782 (batch parameter assignment), Modules/setup_functions.cpp:316-321 (F, A2s1s defaults).
#   A2 preliminary history     Cosmos.cpp:390-483 (calc_coeff_X_spline: Recfast++ Xe_frac, 6000 nodes from zsRe = 2.5e4 to 0), Cosmos.cpp:104-148 (Hubble::init,
#                              node-0 quirk reproduced by `hubble_table`), Cosmos.cpp:567-578 (table inside (zmin, zmax), analytic closure outside).
#   A3/A8 entry                CosmoRec.cpp:1013-1032 (cosmorec_calc_h_cpp_), :729-986 (batch CosmoRec), CAMB fortran/cosmorec.f90:178-213 (adapter values, grid).
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# ---- native CODATA-2008 constants (physical_consts.h, cgs) used by the Cosmos closure ----
const _CR_PI = 3.1415926535897932385            # Definitions.h:17
const _CR_cl = 2.99792458e+10                   # physical_consts.h:44
const _CR_kB_J = 1.3806504e-23                  # :45
const _CR_kB = _CR_kB_J * 1.0e+7                # :46
const _CR_h = 6.62606896e-27                    # :49
const _CR_hbar = 1.054571628e-27                # :50
const _CR_G = 6.67428e-8                        # :55
const _CR_Ry_inf_erg = 2.17987197e-11           # :68 (= const_EH_inf_ergs, :101)
const _CR_me_gr = 9.10938215e-28                # :74
const _CR_me_mp = 5.4461702177e-4               # :81
const _CR_me_malp = 1.37093355570e-4            # :82
const _CR_amu = 1.660538782e-24                 # :84
const _CR_mH_gr = 1.0078250321 * _CR_amu        # :86
const _CR_R_me = 2.8179402894e-13               # :107
const _CR_sigT_0 = 8.0 * _CR_PI / 3.0 * _CR_R_me * _CR_R_me   # :111
const _CR_mHe4_mH = 3.97152594                  # :125
const _CR_fac_mHemH = _CR_mHe4_mH / 4.0         # :148
const _CR_Mpc = 3.08568025e+24                  # :155
const _CR_EionHeI_J = 3.93933e-18               # Cosmos.cpp:673, :688 (literal)
const _CR_A2s1s = 8.2206                        # physical_consts.h:17 const_HI_A2s_1s_0 (batch runpars[3] = 0 -> default; setup_functions.cpp:319-321)
const _CR_F = 1.14                              # CosmoRec.cpp:780
const _CR_zsRe = 2.5e+4                         # Cosmos.cpp:292
const _CR_nXe = 6000                            # Cosmos.cpp:291
const _CR_zsaha = 3500.0                        # Cosmos::X1s/XHeI1s/XHeII1s

# C `pow` as the native object code evaluates it (glibc libm on Linux; GCC folds pow(x, 2.0) to x*x, which is the same correctly rounded value).
# Julia's own `^` agrees to <= 1 ulp; glibc is used where available so that the closure reproduces the native doubles bitwise.
const _CR_LIBM = Sys.islinux() ? "libm.so.6" : ""
_cpow(x::Float64, y::Float64) = isempty(_CR_LIBM) ? x^y : ccall((:pow, _CR_LIBM), Float64, (Float64, Float64), x, y)

"""Native `Cosmos::calc_Orel(TCMB0, Nnu, h100)` (Cosmos.h:326-333)."""
function native_calc_Orel(T0::Float64, Nnu::Float64, h100::Float64)
    H0 = 100.0 * h100 * (1.0e+3 * 1.0e+2) / _CR_Mpc
    a = _CR_PI * _CR_PI / 15.0 * _CR_kB * T0 * _cpow(_CR_kB * T0 / (_CR_hbar * _CR_cl), 3.0) / (_CR_cl * _CR_cl)
    b = 3.0 * (H0 * H0) / (8.0 * _CR_PI * _CR_G)
    return a / b * (1.0 + Nnu * (7.0 / 8.0) * _cpow(4.0 / 11.0, 4.0 / 3.0))
end

"""
    cosmos_constants(; Omega_c, Omega_b, Omega_k, N_eff, h100, T_CMB, Y_p) -> CosmosConstants

Native Cosmos closure for the CosmoRec batch entry (`CosmoRec.cpp:766-775`: `Omega_m = Omega_c + Omega_b`, `Omega_L = 0` -> closure) evaluated as
`Cosmos::init` / `set_initials` / `calc_dependent_consts` / `sigT` do (Cosmos.cpp:207-335, 1141), with the native CODATA-2008 cgs constants:
`O_L = 1 - O_m - O_k - calc_Orel(T0, N_eff, h100)`, `H0 = 100 h 1e5 / Mpc` [1/s], `rho_c = 3 H0^2/(8 pi G)`, `rho_g = pi^2/15 kB T0 (kB T0/(hbar c))^3 / c^2`,
`O_rel = rho_g/rho_c g_rel`, `g_rel = 1 + N_eff 7/8 (4/11)^(4/3)`, `Nb0 = O_b rho_c / m_H`, `sigT0 = 8 pi/3 r_e^2` (aS = mS = 1, p = 0).
Inputs are the density PARAMETERS (not omega h^2), as passed through the C ABI; `Y_p` is the helium MASS fraction (CAMB `YHe`).
"""
function cosmos_constants(; Omega_c::Real, Omega_b::Real, Omega_k::Real, N_eff::Real, h100::Real, T_CMB::Real, Y_p::Real)
    Oc, Ob, Ok, Nnu, h, T0, Yp = Float64.((Omega_c, Omega_b, Omega_k, N_eff, h100, T_CMB, Y_p))
    O_m = Oc + Ob
    O_L = 1.0 - O_m - Ok - native_calc_Orel(T0, Nnu, h)
    g_rel = 1.0 + Nnu * (7.0 / 8.0) * _cpow(4.0 / 11.0, 4.0 / 3.0)
    H0 = 100.0 * h * (1.0e+3 * 1.0e+2) / _CR_Mpc
    rho_c = 3.0 * (H0 * H0) / (8.0 * _CR_PI * _CR_G)
    rho_g = _CR_PI * _CR_PI / 15.0 * _CR_kB * T0 * _cpow(_CR_kB * T0 / (_CR_hbar * _CR_cl), 3.0) / (_CR_cl * _CR_cl)
    O_rel = rho_g / rho_c * g_rel
    Nb0 = Ob * rho_c / _CR_mH_gr
    return CosmosConstants(Yp, T0, _CR_fac_mHemH, _CR_zsRe, Nb0, H0, Ok, O_L, O_m, O_rel, rho_g / _CR_me_gr, _CR_sigT_0, _CR_cl,
                           _CR_h * 1.0e-7, _CR_kB_J, _CR_me_gr * 1.0e-3, _CR_Ry_inf_erg, _CR_me_mp, _CR_me_malp, _CR_EionHeI_J, Float64(_CR_PI), _CR_zsaha)
end

"""The Recfast++ parameter vector `θ = [F, A2s1s, Y_p, Omega_b, h100, T0]` of the batch entry (F = 1.14, A2s1s = 8.2206 for runpars[3] = 0)."""
recfast_theta(; Omega_b::Real, h100::Real, T_CMB::Real, Y_p::Real) = Float64[_CR_F, _CR_A2s1s, Y_p, Omega_b, h100, T_CMB]

"""CAMB adapter conversion (fortran/cosmorec.f90:165-169): `h2 = (H0/100)^2`, `Omega_b = ombh2/h2`, `Omega_c = omch2/h2`."""
function camb_density_parameters(ombh2::Real, omch2::Real, H0::Real)
    h2 = (Float64(H0) / 100)^2
    return (Omega_b = Float64(ombh2) / h2, Omega_c = Float64(omch2) / h2)
end

"""
    preliminary_history(c, θ, hubble, solve_ode; npts = 6000) -> (accessors, history, grid)

In-pipeline replacement of the loaded native preliminary history: native `Cosmos::calc_coeff_X_spline` (Cosmos.cpp:390-483), i.e. the Recfast++ `Xe_frac` run
from `c.zsRe` (2.5e4) to 0 on 6000 nodes, using the Cosmos `H(z)` (the loaded table for `zmin < z < zmax`, the analytic closure of `c` outside it, e.g. for
1e4 <= z <= 2.5e4), then the six native splines and the full accessors. `solve_ode` is the caller's stiff solve callback of [`recfast_history`](@ref)
(it may return `(states, derivatives)`; the derivative of its dense output reproduces the native `Sz.dy`).
"""
function preliminary_history(c::CosmosConstants, θ, hubble, solve_ode; npts::Int = _CR_nXe, rc::RecfastConstants = NATIVE_RECFAST_CONSTANTS)
    hacc = CosmosAccessors(c, nothing, hubble)
    Hfun = z -> cosmos_H(hacc, z)
    grid = recfast_grid(θ, rc, Hfun; npts = npts, zstart = c.zsRe)
    h = recfast_history(θ, rc, Hfun, grid, solve_ode)
    sp = recfast_splines(c, h.z, h.Xe_H, h.Xe_He, h.Xe, h.dXe, h.dXe_H, h.TM)
    return (accessors = CosmosAccessors(c, sp, hubble), history = h, grid = grid)
end

# ---- A3: configuration of the CAMB adapter (fortran/cosmorec.f90:178-196) ----
"""`runpars(1:21)` exactly as the CAMB 2.0.4 adapter sets them for the defaults (fdm = 0, accuracy = 0, every override = -1)."""
const CAMB_COSMOREC_RUNPARS = (0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 2.0, 1.0, 1.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0)
# overrides runpars[16..21] (Diff_iteration_max, nShells, nShellsHeI, HI_absorption, nS_2gamma, nS_Raman): -1 or the batch default itself (CosmoRec.cpp:833-852) are neutral
const _RUNMODE0_OVERRIDE_DEFAULTS = (2.0, 3.0, 2.0, 2.0, 3.0, 2.0)
const _RUNPARS_NAMES = ("fdm", "accuracy", "output", "A2s1s", "B0", "nB", "PMF mode", "alpha scale", "me scale", "power-law p", "VFC mode", "PCA EM1", "PCA EM2",
                        "PCA EM3", "PCA EM4", "Diff_iteration_max", "nShells", "nShellsHeI", "HI_absorption", "nS_2gamma", "nS_Raman")

"""Raise `ArgumentError` unless `(runmode, runpars)` is the ported configuration (runmode 0, accuracy 0, fdm 0, all other entries at the CAMB adapter values, overrides neutral)."""
function check_runmode0_config(runmode::Integer, runpars::AbstractVector{<:Real})
    runmode == 0 || throw(ArgumentError("CosmoRec.jl implements runmode 0 only (got runmode = $runmode)"))
    length(runpars) == 21 || throw(ArgumentError("runpars must have 21 entries (CosmoRec.vX batch layout + CAMB overrides), got $(length(runpars))"))
    for k in 1:15
        Float64(runpars[k]) == CAMB_COSMOREC_RUNPARS[k] ||
            throw(ArgumentError("unsupported runpars[$k] ($(_RUNPARS_NAMES[k])) = $(runpars[k]); only the CAMB default $(CAMB_COSMOREC_RUNPARS[k]) is ported"))
    end
    for k in 16:21
        v = Float64(runpars[k])
        (v == -1.0 || v == _RUNMODE0_OVERRIDE_DEFAULTS[k - 15]) ||
            throw(ArgumentError("unsupported override runpars[$k] ($(_RUNPARS_NAMES[k])) = $v; only -1 or the batch default $(_RUNMODE0_OVERRIDE_DEFAULTS[k - 15]) is ported"))
    end
    return nothing
end

"""
    check_camb_grid(z; name = "z", n = 10000, zmax = 1e4, atol = 1e-9)

A8: the CAMB adapter grid (fortran/cosmorec.f90:24-26, 199-200): `n` nodes, strictly descending, `z[1] == zmax`, every node within `atol` of
`zmax - (i-1) zmax/(n-1)` (compiled grids differ from the exact formula by ulps), last node in `[0, atol]` (CAMB's -march=native build gives 8.3e-13;
the native Hubble loader treats an exact 0 as 1e-10 and any other value as itself, Cosmos.cpp:125-126, reproduced by `hubble_table`).
"""
function check_camb_grid(z::AbstractVector{<:Real}; name::AbstractString = "z", n::Int = 10000, zmax::Float64 = 1.0e4, atol::Float64 = 1.0e-9)
    length(z) == n || throw(ArgumentError("$name: expected the $n-node CAMB grid, got $(length(z)) nodes"))
    all(isfinite, z) || throw(ArgumentError("$name: non-finite redshift"))
    z[1] == zmax || throw(ArgumentError("$name: first node must be $zmax (descending CAMB grid), got $(z[1])"))
    all(i -> z[i + 1] < z[i], 1:(n - 1)) || throw(ArgumentError("$name: grid must be strictly descending"))
    dz = zmax / (n - 1)
    for i in 1:n
        abs(z[i] - (zmax - (i - 1) * dz)) <= atol || throw(ArgumentError("$name: node $i = $(z[i]) is not on the uniform CAMB grid (expected $(zmax - (i - 1) * dz))"))
    end
    0.0 <= z[n] <= atol || throw(ArgumentError("$name: last node must lie in [0, $atol], got $(z[n])"))
    return nothing
end

"""Explicit production tables of the runmode-0 calculation: the effective-rate model (`load_native_effective_model`) and the HI-diffusion inputs."""
struct Runmode0Tables{E,D}
    eff::E
    diffusion::D
end

"""
    load_runmode0_tables(native_root) -> Runmode0Tables

Loads the production tables from an ORIGINAL CosmoRec source tree (read-only; never bundled): `Rec_database` (HI nS_3 / nS_eff 500, HeI res_2 / nS_eff 30,
DP table), `Development/Recombination/Data.fcorr/f.corr.dat`, `Development/Line_profiles/two-photon-data` (A2s1s = 8.2206) and the HI ln Bitot table; the
PDE range is the batch one (CosmoRec.cpp:398-399 with Diff_corr_zmax = 2000, Diff_corr_zmin = 500, f_t = 1, zstart 3000, zend 50:
zs = min(2000*1.25, 3000/1.002) = 2500, ze = max(500, 50*1.002) = 500), Diff_iteration_max = 2.
"""
function load_runmode0_tables(native_root::AbstractString)
    db = joinpath(native_root, "Rec_database")
    fc = joinpath(native_root, "Development", "Recombination", "Data.fcorr", "f.corr.dat")
    tp = joinpath(native_root, "Development", "Line_profiles", "two-photon-data")
    for p in (db, fc, tp)
        ispath(p) || throw(ArgumentError("native data path does not exist: $p"))
    end
    eff = load_native_effective_model(db, fc)
    lnB = load_hi_bitot_table(joinpath(db, "Effective_Rates.HI", "Effective_Rate_Tables.nS_3"), 500, NATIVE_HYDROGEN_RESOLVED)
    setup = hi_pde_setup(load_hi_profile_data(tp; A2s1s = _CR_A2s1s))
    return Runmode0Tables(eff, HIDiffusionInputs(setup, eff.htable, lnB))
end

"""
    cosmorec_calc_h(runmode, runpars, Omega_c, Omega_b, Omega_k, N_eff, H0, T_CMB, Y_He, z_Hz, Hz, z_out, tables;
                    solve_recfast, solve_pass, solve_tail, kwargs...) -> NamedTuple

Standalone counterpart of the native `cosmorec_calc_h_cpp_` (CosmoRec.cpp:1013-1032) for the CAMB adapter route, with ZERO native-provenance inputs:
argument meaning and units as in the C ABI (`Omega_*` density parameters, `H0` in km/s/Mpc, `T_CMB` in K, `Y_He` the helium mass fraction, `z_Hz`/`Hz` the
descending CAMB table with `Hz` in 1/s, `z_out` the output grid). The configuration is checked (`check_runmode0_config`), both grids are checked
(`check_camb_grid`) and `Hz` must be finite and positive. Steps: closure constants (`cosmos_constants`), loaded Hubble table (`hubble_table`, node-0 quirk
reproduced), preliminary Recfast++ history (`preliminary_history` with `solve_recfast`), runmode-0 history (`recombination_history_diffusion` with
`solve_pass` / `solve_tail`; extra keywords are forwarded). The solve callbacks are REQUIRED (no library solver defaults).
Returns `Xe` (= n_e/n_H) and `Tm` [K] on `z_out`, plus the constants, θ, preliminary history and the full history for diagnostics.
"""
function cosmorec_calc_h(runmode::Integer, runpars::AbstractVector{<:Real}, Omega_c::Real, Omega_b::Real, Omega_k::Real, N_eff::Real, H0::Real, T_CMB::Real,
                         Y_He::Real, z_Hz::AbstractVector{<:Real}, Hz::AbstractVector{<:Real}, z_out::AbstractVector{<:Real}, tables::Runmode0Tables;
                         solve_recfast, solve_pass, solve_tail, kwargs...)
    check_runmode0_config(runmode, runpars)
    check_camb_grid(z_Hz; name = "z_Hz")
    check_camb_grid(z_out; name = "z_out")
    length(Hz) == length(z_Hz) || throw(DimensionMismatch("Hz and z_Hz differ in length"))
    all(h -> isfinite(h) && h > 0, Hz) || throw(ArgumentError("Hz must be finite and positive (1/s)"))
    for (nm, v) in (("Omega_c", Omega_c), ("Omega_b", Omega_b), ("Omega_k", Omega_k), ("N_eff", N_eff), ("H0", H0), ("T_CMB", T_CMB), ("Y_He", Y_He))
        isfinite(v) || throw(ArgumentError("$nm must be finite, got $v"))
    end
    (Omega_b > 0 && Omega_c >= 0 && H0 > 0 && T_CMB > 0 && N_eff >= 0 && 0 < Y_He < 1) ||
        throw(ArgumentError("unphysical cosmology: Omega_b = $Omega_b, Omega_c = $Omega_c, H0 = $H0, T_CMB = $T_CMB, N_eff = $N_eff, Y_He = $Y_He"))
    h100 = Float64(H0) / 100.0                                   # CosmoRec.cpp:1027 (*h0 / 100.0)
    c = cosmos_constants(; Omega_c = Omega_c, Omega_b = Omega_b, Omega_k = Omega_k, N_eff = N_eff, h100 = h100, T_CMB = T_CMB, Y_p = Y_He)
    θ = recfast_theta(; Omega_b = Omega_b, h100 = h100, T_CMB = T_CMB, Y_p = Y_He)
    pre = preliminary_history(c, θ, hubble_table(z_Hz, Hz), solve_recfast)
    rm = RecombinationModel(tables.eff, pre.accessors)
    h = recombination_history_diffusion(rm, θ, tables.diffusion, solve_pass, solve_tail, z_out; kwargs...)
    return (Xe = h.Xe, Tm = h.Te, z = z_out, k_switch = [p.k_switch for p in h.passes], constants = c, theta = θ, preliminary = pre, history = h)
end

"""
    cosmorec_history(; ombh2, omch2, omk = 0, N_eff, H0, T_CMB, Y_He, z_Hz, Hz, z_out = z_Hz, tables, solve_recfast, solve_pass, solve_tail, kwargs...)

Convenience front end with CAMB's physical inputs (`CP%ombh2`, `CP%omch2`, `CP%omk`, `CP%N_eff()`, `CP%H0`, `CP%tcmb`, `CP%Yhe`): converts like the adapter
(`camb_density_parameters`) and calls [`cosmorec_calc_h`](@ref) with the CAMB default `runpars`.
"""
function cosmorec_history(; ombh2, omch2, omk = 0.0, N_eff, H0, T_CMB, Y_He, z_Hz, Hz, z_out = z_Hz, tables, solve_recfast, solve_pass, solve_tail, kwargs...)
    d = camb_density_parameters(ombh2, omch2, H0)
    return cosmorec_calc_h(0, collect(CAMB_COSMOREC_RUNPARS), d.Omega_c, d.Omega_b, omk, N_eff, H0, T_CMB, Y_He, z_Hz, Hz, z_out, tables;
                           solve_recfast = solve_recfast, solve_pass = solve_pass, solve_tail = solve_tail, kwargs...)
end
