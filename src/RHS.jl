# Pure-Julia port of the smooth hydrogen effective-rate / Compton / population RHS components of ORIGINAL CosmoRec v3.0b
# (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3):
#   Development/Recombination/ODE_effective.cpp: evaluate_TM (:85-98), evaluate_2s_two_photon_decay (:152-162, Lambda_ind=1 at :131-133),
#     evaluate_Ly_n_channel (:~176-190, tau_S/p_ij in namespace ODE_HI_effective_funcs :34-62), evaluate_effective_Rci_Ric_terms, evaluate_effective_Rij_terms (:390-410)
#   Modules/ODEdef_CosmoRec.cpp: fcn_HI_effective (:140-210; diffusion/quadrupole branches OFF = production default of this chunk)
# CosmoRec is (c) J. Chluba et al.; use must be acknowledged and Chluba & Thomas 2010 (MNRAS 412, 748) and
# Chluba & Sunyaev 2006 (A&A 446, 39) cited (see docs/CHUNK3A_RESULTS.md NOTICE).
#
# Time variable of every right-hand side here is cosmic time (d/dt, 1/s) exactly as native; the native solver divides by dz/dt = -H(1+z) outside these terms.
# State X (length 6, per H nucleus): 1s, 2s, 2p, 3s, 3p, 3d (native Get_Level_index = l + n(n-1)/2). Resolved states m = 1..5 map to X[m+1].
# The component functions ACCUMULATE into dX like native (`+=`); only `matter_temperature_rate` assigns (native drho_dt = ...).
# Not included here: helium, H absorption, diffusion correction, quadrupole lines, rescale_rates, initialization/Saha, Recfast handoff, switches, f_t scaling of Tg.

"""Physical constants used by the RHS components (native physical_consts.h values; `sigT` is the native global `const_sigT` as read by the harness)."""
struct RHSConstants
    pi::Float64
    sigT::Float64
    c::Float64
    kB::Float64
    hbar::Float64
    me_gr::Float64
    h_kb::Float64
    A2s1s::Float64
end

const NATIVE_RHS_CONSTANTS = RHSConstants(
    3.1415926535897931, 6.6524585583988891e-25, 29979245800.0, 1.3806504000000002e-16,
    1.0545716280000001e-27, 9.1093821500000007e-28, 4.7992373449498863e-11, 8.2205999999999992,
)

"""Fixed atomic data of the production 3-shell hydrogen atom, as read from the native objects (Float64, non-differentiated).
`nu_ul[i, j]` = native `Level(i).nu_ul(n_i, n_j)` (signed, Hz) for resolved states i<j; `ly_*` are the n = 2, 3 Lyman channels."""
struct HydrogenAtom
    gw::Vector{Float64}
    nu_ul::Matrix{Float64}
    nu21_2s::Float64
    ly_n::Vector{Int}
    ly_A21::Vector{Float64}
    ly_lambda21::Vector{Float64}
    ly_nu21::Vector{Float64}
    ly_index::Vector{Int}
end

const NATIVE_HYDROGEN_ATOM = let nu = zeros(5, 5)
    for i in 1:2, j in 3:5
        nu[i, j] = -456673782179412.38
    end
    HydrogenAtom([2.0, 6.0, 2.0, 6.0, 10.0], nu, 2466038423768827.0, [2, 3],
        [626490298.80265749, 167252720.15153983], [1.2156844561319915e-05, 1.025733759861368e-05],
        [2466038423768827.0, 2922712205948239.0], [3, 5])
end

"""Factor of native `rho_g_1_cm3(Tg) = rho_g_fac Tg^4` (ODE_effective.cpp:43)."""
rho_g_fac(c::RHSConstants) = c.pi^2 / 15 * c.kB * (c.kB / (c.hbar * c.c))^3 / c.c^2 / c.me_gr

"""Electron-temperature ratio rho = Te/Tg: d rho/dt = 8/3 sigT c rho_g(Tg) Xe/(1+Xe+fHe) (1-rho) - H rho  [1/s]  (native evaluate_TM; assigned, not accumulated)."""
function matter_temperature_rate(rho, Tg, Xe, fHe, Hz, c::RHSConstants = NATIVE_RHS_CONSTANTS)
    compton = 8 / 3 * c.sigT * c.c * (rho_g_fac(c) * Tg^4) * Xe / (1 + Xe + fHe) * (1 - rho)
    return compton - Hz * rho
end

"""2s-1s two-photon decay, `dX[1] += d`, `dX[2] -= d`, `d = A2s1s (X2s - X1s exp(-h_kb nu21/Tg))` (Lambda_ind = 1)."""
function two_photon!(dX, Tg, X, atom::HydrogenAtom = NATIVE_HYDROGEN_ATOM, c::RHSConstants = NATIVE_RHS_CONSTANTS)
    d = c.A2s1s * (X[2] - X[1] * exp(-c.h_kb * atom.nu21_2s / Tg))
    dX[1] += d
    dX[2] -= d
    return dX
end

"""Single Sobolev Lyman-n channel (native evaluate_Ly_n_channel, w = 3): returns nothing, updates dX[1] and dX[idx]."""
function lyman_channel!(dX, idx::Int, Tg, X1s, Xnp, NH, Hz, A21, lambda21, nu21, c::RHSConstants = NATIVE_RHS_CONSTANTS, w = 3.0)
    tauS = A21 * lambda21^3 / (8 * c.pi * Hz) * (X1s * NH * w - Xnp * NH)
    pij = (1 - exp(-tauS)) / tauS
    efac = exp(-c.h_kb * nu21 / Tg)
    d = pij * A21 / (1 - efac) * (Xnp - w * X1s * efac)
    dX[1] += d
    dX[idx] -= d
    return nothing
end

"""All resolved Lyman channels n = 2..3 of the 3-shell atom (native loop n=2..nLy_max)."""
function lyman!(dX, Tg, X, NH, Hz, atom::HydrogenAtom = NATIVE_HYDROGEN_ATOM, c::RHSConstants = NATIVE_RHS_CONSTANTS)
    for k in eachindex(atom.ly_n)
        idx = atom.ly_index[k]
        lyman_channel!(dX, idx, Tg, X[1], X[idx], NH, Hz, atom.ly_A21[k], atom.ly_lambda21[k], atom.ly_nu21[k], c)
    end
    return dX
end

"""Effective continuum: `dX[m+1] += B[m] (A[m] Xe Np - X[m+1])`, `Np = NH Xp`."""
function continuum!(dX, Xe, Np, X, A, B)
    for m in eachindex(A)
        dX[m + 1] += B[m] * (A[m] * Xe * Np - X[m + 1])
    end
    return dX
end

"""Interlevel effective rates for resolved i<j: `d = R[i,j] (X_i - gw_i/gw_j exp(-h_kb nu_ul[i,j]/Tg) X_j)`, `dX_i -= d`, `dX_j += d`."""
function interlevel!(dX, Tg, X, R, atom::HydrogenAtom = NATIVE_HYDROGEN_ATOM, c::RHSConstants = NATIVE_RHS_CONSTANTS)
    n = length(atom.gw)
    for i in 1:n, j in (i + 1):n
        xij = c.h_kb * atom.nu_ul[i, j] / Tg
        d = R[i, j] * (X[i + 1] - atom.gw[i] / atom.gw[j] * exp(-xij) * X[j + 1])
        dX[i + 1] -= d
        dX[j + 1] += d
    end
    return dX
end

"""Combined native `fcn_HI_effective` (diffusion and quadrupole off) from already-evaluated rates `A, B, R`; accumulates into `dX` (length 6)."""
function hydrogen_rhs!(dX, Tg, Xe, Xp, NH, Hz, X, A, B, R, atom::HydrogenAtom = NATIVE_HYDROGEN_ATOM, c::RHSConstants = NATIVE_RHS_CONSTANTS)
    two_photon!(dX, Tg, X, atom, c)
    lyman!(dX, Tg, X, NH, Hz, atom, c)
    continuum!(dX, Xe, NH * Xp, X, A, B)
    interlevel!(dX, Tg, X, R, atom, c)
    return dX
end

"""As [`hydrogen_rhs!`](@ref), evaluating the effective rates from the explicit `table` at `(Tg, Te = rho Tg)`; returns `d rho/dt` from [`matter_temperature_rate`](@ref)."""
function hydrogen_rhs!(dX, table::AtomicRateTable, Tg, rho, Xe, Xp, fHe, NH, Hz, X, atom::HydrogenAtom = NATIVE_HYDROGEN_ATOM, c::RHSConstants = NATIVE_RHS_CONSTANTS)
    r = get_rates(table, Tg, rho * Tg)
    hydrogen_rhs!(dX, Tg, Xe, Xp, NH, Hz, X, r.A, r.B, r.R, atom, c)
    return matter_temperature_rate(rho, Tg, Xe, fHe, Hz, c)
end

"""Native ground-state reconstruction: `xp = 1 - X_{H,1s}`, `xHeII = fHe - X_{HeI,1s}`, `xe = xp + xHeII` (flag_He = 1)."""
proton_fraction(XH1s) = 1 - XH1s
electron_fraction(XH1s, fHe, XHeI1s) = (1 - XH1s) + (fHe - XHeI1s)
