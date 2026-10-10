module CosmoRecMooncakeExt

# Mooncake bridge of the portable ChainRulesCore rule (reverse mode only, DefaultCtx) and the non-differentiable status of the quadrature plans.
using CosmoRec: CosmoRec, HIQuadPlan, HIQuadPlans, HIFitScratch, HIPDESetup, HIPDEAtom, _hi_spline_parts, _hi_fit_scratched, _hi_check_levels, _check_quadplan
using ChainRulesCore: ChainRulesCore
using Mooncake: Mooncake, DefaultCtx, NoTangent, @from_rrule, @zero_derivative, @mooncake_overlay

Mooncake.tangent_type(::Type{<:HIQuadPlan}) = NoTangent
Mooncake.tangent_type(::Type{<:HIQuadPlans}) = NoTangent

@from_rrule DefaultCtx Tuple{typeof(_hi_spline_parts), HIQuadPlan{K}, Vector{Float64}, NTuple{K,Int}} where {K}
# Caller-owned fit scratch (`HIFitScratch`): NO primitive rule. The scratch is MUTATED by the fit and its arrays (y, c, z, b, d) can be read, and hence differentiated, after the call, so a rule that
# returns only the scalar parts (and NoTangent for the scratch) would silently erase those gradients. Under Mooncake the 5-argument method is therefore REPLACED by (i) the fit traced as ordinary
# mutating code, which carries the cotangents of every later read of the scratch back to F, and (ii) the existing 3-argument primitive for the quadrature parts (same values: identical
# arithmetic, so the primal is bitwise the scratch result). Outside Mooncake the fast scratch method runs unchanged. No ChainRulesCore rule for the scratch method exists.
@mooncake_overlay function _hi_spline_parts(q::HIQuadPlan{K}, F::Vector{T}, forced::NTuple{K,Int}, w::HIFitScratch{T}) where {K,T<:Real}
    length(F) == length(q.knots) || throw(DimensionMismatch("_hi_spline_parts: length(F) != number of plan knots"))      # same checks, same ORDER, as the fast scratch method: before any mutation
    _hi_check_levels(forced, "_hi_spline_parts forced levels")
    _hi_fit_scratched(w, q, F)
    return _hi_spline_parts(q, F, forced)
end
@zero_derivative DefaultCtx Tuple{typeof(_check_quadplan), HIQuadPlans, HIPDESetup, HIPDEAtom}

end
