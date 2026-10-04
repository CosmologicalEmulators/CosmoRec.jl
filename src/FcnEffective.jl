# Pointwise assembly of the first-pass default H/He right-hand side of ORIGINAL CosmoRec v3.0b
# (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3): Modules/ODEdef_CosmoRec.cpp `fcn_effective(double z, Data_Level_I &LI)` (:329-) with
#   compute_fractions (:~95-120), ODE_effective::evaluate_TM, fcn_HI_effective (Chunk 3a), fcn_HeI_effective base terms (Chunk 3b),
#   and the HI-absorption block of fcn_HeI_effective (:284-321; Chunk 3c).
# CosmoRec is (c) J. Chluba et al.; use must be acknowledged and Chluba & Thomas 2010 (MNRAS 412, 748) and Chluba & Sunyaev 2006
# (A&A 446, 39) cited (see docs/CHUNK3D_RESULTS.md NOTICE).
#
# PRODUCTION SWITCHES (native names, values of the CAMB batch runmode 0 first EMLA pass; verified at runtime by the native harness, see
# docs/CHUNK3D_RESULTS.md):
#   Diffusion_correction_is_on = 0 (H Ly-a diffusion), Diffusion_correction_HeI_is_on = 0 (HeI diffusion), HeISTfeedback = 0,
#   DM_annihilation = DM_decay = magnetic_fields = 0, FC_CosmoRec.runmode = 0 (f_t = f_b = 1), flag_He = 1,
#   parameters.CR.HI_absorption = 1 with _HI_abs_appr_flag = 1, spin_forbidden = 1, Atom_activate_HI_Quadrupole_lines = 0.
# Background (Tg = TCMB(z), NH(z), H(z), fHe) is an explicit input: the Cosmology module is not ported.
#
# State layout (native Data_Level_I, 1-based Julia indices): X[1] = Xe (ignored by fcn_effective, g[1] only receives the absorber),
#   X[2:7] = HI 1s,2s,2p,3s,3p,3d, X[8:14] = HeI levels (native index 0..6), X[15] = rho = Tm/Tg.  g = dX/dt in cosmic time [1/s].

const EFF_NEQ = 15
const EFF_IHI = 2
const EFF_IHE = 8
const EFF_NHI = 6
const EFF_NHE = 7

"""Explicit background at the evaluation redshift: `Tg = TCMB(z)` [K], `NH` [cm^-3], `Hz = H(z)` [1/s], `fHe = nHe/nH`."""
struct EffectiveBackground{T1,T2,T3,T4}
    Tg::T1
    NH::T2
    Hz::T3
    fHe::T4
end

"""Explicit model/tables workspace of [`fcn_effective`](@ref): no hidden global caches, nothing is loaded by default."""
struct EffectiveModel
    htable::AtomicRateTable
    hetable::HeliumRateTable
    dp::DPTable
    bitot::BitotSeries
    fc::FcorrSpline
    hatom::HydrogenAtom
    heatom::HeliumAtom
    hconst::RHSConstants
    heconst::HeliumConstants
    hiabs::HIAbsConstants
    spin_forbidden::Bool
    hi_absorption::Bool
end

function EffectiveModel(htable, hetable, dp, bitot, fc; hatom = NATIVE_HYDROGEN_ATOM, heatom = NATIVE_HELIUM_ATOM, hconst = NATIVE_RHS_CONSTANTS,
                        heconst = NATIVE_HELIUM_CONSTANTS, hiabs = NATIVE_HIABS_CONSTANTS, spin_forbidden::Bool = true, hi_absorption::Bool = true)
    return EffectiveModel(htable, hetable, dp, bitot, fc, hatom, heatom, hconst, heconst, hiabs, spin_forbidden, hi_absorption)
end

"""Native `compute_fractions` with `flag_He = 1`: `(Xe, Xp, XHeII)`, `Xp = 1 - X[HI 1s]`, `XHeII = fHe - X[HeI 1s]`, `XHeIII = 0`, `Xe = Xp + XHeII`."""
function effective_fractions(X, fHe)
    xp = 1 - X[EFF_IHI]
    xheii = fHe - X[EFF_IHE]
    return (xp + xheii, xp, xheii)
end

"""
    fcn_effective!(g, z, X, bg, model) -> g

Fill `g` (length 15, overwritten like native: zeroed, rho equation assigned, HI and HeI terms accumulated, absorber added last) with the
derivative of the native state `X` at redshift `z`. `g` must be able to hold the promoted element type of `X`/`bg`.
`diffusion` (an `HIDiffusionFeedback`, Chunk 9a) adds the HI PDE corrections of a diffusion iteration; `nothing` reproduces the first pass.
Throws `RateTableDomainError` where native exits; `DPTableDomainError` where native runs the explicit DP integral unless `dp_fallback` (`dpesc_fallback`-style factory `triplet -> (Tg, eta, tauS, pd) -> correction`) is supplied.
"""
function fcn_effective!(g, z, X, bg::EffectiveBackground, m::EffectiveModel; dp_fallback::F = nothing, flag_He::Bool = true, diffusion = nothing) where {F}
    (length(X) == EFF_NEQ && length(g) == EFF_NEQ) || throw(DimensionMismatch("fcn_effective! expects native neq = 15"))
    T = eltype(g)
    return _fcn_effective_core!(g, z, X, bg, m, zeros(T, EFF_NHI), zeros(T, EFF_NHE), nothing, dp_fallback, flag_He, diffusion)
end

# Private core shared by the public allocating `fcn_effective!` (rates = nothing: the allocating `get_rates`/`get_helium_rates` paths, fresh dXH/dXHe)
# and the private workspace RHS (`RHSWorkspace.jl`: rates = an `_RHSWorkspace`, dXH/dXHe from it, zeroed here). Same statements and order as before.
function _fcn_effective_core!(g, z, X, bg::EffectiveBackground, m::EffectiveModel, dXH, dXHe, rates, dp_fallback::F, flag_He::Bool, diffusion) where {F}
    Tg, NH, Hz, fHe = bg.Tg, bg.NH, bg.Hz, bg.fHe
    rho = X[EFF_NEQ]
    xe, xp, _ = effective_fractions(X, fHe)
    flag_He || (xe = xp)           # native compute_fractions with flag_He = 0: XHeII = 0, Xe = Xp (Chunk 5a)
    T = eltype(g)
    for i in 1:EFF_NEQ
        g[i] = zero(T)
    end
    rates === nothing || fill!(dXH, zero(T))
    g[EFF_NEQ] = _hydrogen_rhs_rates!(dXH, rates, m.htable, Tg, rho, xe, xp, fHe, NH, Hz, view(X, EFF_IHI:(EFF_IHI + EFF_NHI - 1)), m.hatom, m.hconst)
    for i in 1:EFF_NHI
        g[EFF_IHI + i - 1] += dXH[i]
    end
    hi_diffusion_rhs!(g, diffusion, z, X)   # native: end of fcn_HI_effective when Diffusion_correction_is_on (Chunk 9a); `nothing` = first pass
    XHe = view(X, EFF_IHE:(EFF_IHE + EFF_NHE - 1))
    if flag_He        # native: `if(flag_He==1) fcn_HeI_effective(...)` (includes the HI absorber of the HeI lines)
        rates === nothing || fill!(dXHe, zero(T))
        _helium_base_rhs_rates!(dXHe, rates, m.hetable, Tg, xe, NH, Hz, XHe, fHe, m.heatom, m.heconst, m.hconst, m.spin_forbidden)
        for i in 1:EFF_NHE
            g[EFF_IHE + i - 1] += dXHe[i]
        end
    end
    if flag_He && m.hi_absorption
        hi_absorption_rhs!(g, 1, EFF_IHI, EFF_IHE, m.dp, m.bitot, m.fc, z, Tg, NH, Hz, X[EFF_IHI], XHe;
                           spin_forbidden = m.spin_forbidden, fcorr_on = true, diffusion_correction = false, c = m.hiabs, fallback = dp_fallback)
    end
    return g
end
# rates === nothing: exactly the previous public calls (allocating rate evaluation)
_hydrogen_rhs_rates!(dX, ::Nothing, table, Tg, rho, Xe, Xp, fHe, NH, Hz, X, atom, c) = hydrogen_rhs!(dX, table, Tg, rho, Xe, Xp, fHe, NH, Hz, X, atom, c)
_helium_base_rhs_rates!(dX, ::Nothing, table, Tg, Xe, NH, Hz, X, fHe, atom, c, rc, spin_forbidden::Bool) =
    helium_base_rhs!(dX, table, Tg, Xe, NH, Hz, X, fHe, atom, c, rc; spin_forbidden = spin_forbidden)

"""Allocating form of [`fcn_effective!`](@ref); element type is the promotion of `X` and the background."""
function fcn_effective(z, X, bg::EffectiveBackground, m::EffectiveModel; dp_fallback::F = nothing, flag_He::Bool = true, diffusion = nothing) where {F}
    T = promote_type(eltype(X), typeof(bg.Tg), typeof(bg.NH), typeof(bg.Hz), typeof(bg.fHe))
    g = Vector{T}(undef, EFF_NEQ)
    return fcn_effective!(g, z, X, bg, m; dp_fallback = dp_fallback, flag_He = flag_He, diffusion = diffusion)
end
