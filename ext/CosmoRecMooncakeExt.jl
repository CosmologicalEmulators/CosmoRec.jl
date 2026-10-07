module CosmoRecMooncakeExt

# Mooncake bridge of the portable ChainRulesCore rule (reverse mode only, DefaultCtx) and the non-differentiable status of the quadrature plans.
using CosmoRec: CosmoRec, HIQuadPlan, HIQuadPlans, HIPDESetup, HIPDEAtom, _hi_spline_parts, _check_quadplan
using ChainRulesCore: ChainRulesCore
using Mooncake: Mooncake, DefaultCtx, NoTangent, @from_rrule, @zero_derivative

Mooncake.tangent_type(::Type{<:HIQuadPlan}) = NoTangent
Mooncake.tangent_type(::Type{<:HIQuadPlans}) = NoTangent

@from_rrule DefaultCtx Tuple{typeof(_hi_spline_parts), HIQuadPlan{K}, Vector{Float64}, NTuple{K,Int}} where {K}
@zero_derivative DefaultCtx Tuple{typeof(_check_quadplan), HIQuadPlans, HIPDESetup, HIPDEAtom}

end
