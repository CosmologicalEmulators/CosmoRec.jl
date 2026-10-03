# Chunk 4d: the sampled-HeI state switch (native `flag_He` 1 -> 0 transition) of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3):
#   Modules/main.CosmoRec.cpp:163-197 (decision, reset, state-size change, solver restart), Modules/global_variables.cpp:37 (Xi_HeI_switch = 1e-7),
#   Modules/ODEdef_CosmoRec.cpp:87-107 (copy_LI_to_ysol, both flag_He branches).
# CosmoRec is (c) J. Chluba et al.; use must be acknowledged and Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39) and
# Rubino-Martin, Chluba & Sunyaev 2008 cited (see docs/CHUNK4D_RESULTS.md NOTICE).
#
# SCOPE: a DISCRETE operator on explicit state. The decision `(zs < 200 || fHe - X_He1s <= Xi_HeI_switch) && flag_He` is evaluated on PRIMAL values and the reset is applied with the
# branch frozen: `X_He1s = fHe`, all other helium slots `1e-300` (the native floor for the log-spline setup), `Xe`, hydrogen and `rho` untouched; the packed ODE vector then drops all helium
# entries (length `2 + nres_HI`). NOT implemented: the ODE, the solver restart, any location of the event in z, or any derivative of the event location (the decision is piecewise constant and
# the packed state changes dimension, so no continuous event sensitivity exists or is claimed). The reset depends on `fHe` (the He ground slot becomes `fHe`).

"""Native threshold `Xi_HeI_switch = 1e-7` and the redshift below which helium is switched off (`zs < 200`)."""
const NATIVE_XI_HEI_SWITCH = 1.0e-7
const NATIVE_HE_SWITCH_ZOFF = 200.0
"""Native floor assigned to the non-ground helium slots after the switch (`1.0e-300`)."""
const NATIVE_HE_FLOOR = 1.0e-300

"""
    helium_switch_condition(zs, fHe, X_He1s; flag_He = true, Xi_switch = 1e-7, zoff = 200) -> Bool

Native decision `(zs < 200.0 || cosmos.fHe() - HeI_Atoms.Xi(0) <= Xi_HeI_switch) && flag_He == 1` evaluated on primal values (`X_He1s` = the HeI ground-state population `Xi(0)` = `X[index_HeI]`).
"""
function helium_switch_condition(zs, fHe, X_He1s; flag_He::Bool = true, Xi_switch::Float64 = NATIVE_XI_HEI_SWITCH, zoff::Float64 = NATIVE_HE_SWITCH_ZOFF)
    return (_primal(zs) < zoff || _primal(fHe) - _primal(X_He1s) <= Xi_switch) && flag_He
end

"""
    helium_switch_reset(X, nH, nHe, fHe) -> X'

The reset of `main.CosmoRec.cpp:171-181`: `X[He 1s] = fHe`, `X[He 1s + 1 : end of He] = 1e-300`; `Xe` (`X[1]`), the hydrogen levels and `rho` (`X[end]`) are unchanged. Layout as in 4a
(`length(X) == 1 + nH + nHe + 1`, He ground at `2 + nH`). Differentiable with respect to `fHe` (He ground) and the retained entries (identity); the helium floor entries are constants.
"""
function helium_switch_reset(X::AbstractVector, nH::Integer, nHe::Integer, fHe)
    length(X) == 1 + nH + nHe + 1 || throw(DimensionMismatch("helium_switch_reset: length(X) = $(length(X)) != 1 + nH + nHe + 1 = $(1 + nH + nHe + 1)"))
    T = promote_type(eltype(X), typeof(fHe))
    Xn = Vector{T}(X)          # always a copy (convert would alias a Vector{Float64} input)
    Xn[2 + nH] = fHe
    for i in (3 + nH):(1 + nH + nHe)
        Xn[i] = NATIVE_HE_FLOOR
    end
    return Xn
end

"""
    helium_switch(X, zs, nH, nHe, fHe; flag_He = true, kwargs...) -> (X', switched::Bool)

Decision plus reset with a FROZEN branch: if [`helium_switch_condition`](@ref) holds on the primal values the reset state is returned (`switched = true`), otherwise `X` unchanged.
The returned state is piecewise defined; no derivative of the decision (event location) is implemented.
"""
function helium_switch(X::AbstractVector, zs, nH::Integer, nHe::Integer, fHe; flag_He::Bool = true, kwargs...)
    sw = helium_switch_condition(zs, fHe, X[2 + nH]; flag_He = flag_He, kwargs...)
    return (sw ? helium_switch_reset(X, nH, nHe, fHe) : X, sw)
end

"""Length of the packed ODE vector (native `neq` of `Xe_frac_effective_rates`): `2 + nres_HI + (helium ? 1 + nres_HeI : 0)`."""
helium_switch_ysize(nres_HI::Integer, nres_HeI::Integer, helium::Bool) = 2 + nres_HI + (helium ? 1 + nres_HeI : 0)
