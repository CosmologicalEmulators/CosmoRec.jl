module CosmoRecChainRulesCoreExt

# Portable reverse rule of the private fixed-geometry spline quadrature (src/HIPDEQuadPlan.jl). Only the active ordinates `F` get a cotangent;
# the plan (grid, edges, precomputed functionals) and the integer level tuple are non-differentiable by construction (opt-in private path).
using CosmoRec: CosmoRec, HIQuadPlan, _hi_spline_parts, _hi_quad_pullback!
using ChainRulesCore: ChainRulesCore, NoTangent, AbstractZero, unthunk

_cr_float(d) = (d = unthunk(d); d isa AbstractZero ? 0.0 : Float64(d))

function ChainRulesCore.rrule(::typeof(_hi_spline_parts), q::HIQuadPlan{K}, F::Vector{Float64}, forced::NTuple{K,Int}) where {K}
    out = _hi_spline_parts(q, F, forced)
    levels = out[2]; n = length(F)
    function _hi_spline_parts_pullback(Δ)
        dF = zeros(n)                       # always a dense Vector{Float64}: never ZeroTangent/NoTangent for the active F
        Δ = unthunk(Δ)
        if !(Δ isa AbstractZero)
            Δp = unthunk(Δ[1])
            if !(Δp isa AbstractZero)
                _hi_quad_pullback!(dF, q, levels, ntuple(i -> _cr_float(Δp[i]), Val(K)))
            end
        end
        return NoTangent(), NoTangent(), dF, NoTangent()
    end
    return out, _hi_spline_parts_pullback
end

end
