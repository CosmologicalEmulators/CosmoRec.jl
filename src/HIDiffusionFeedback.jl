# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Chunk 9a: feedback of the HI PDE correction integrals into the next recombination ODE pass.
# Sources (ORIGINAL CosmoRec v3.0b, read-only): Modules/Diffusion_correction.cpp (setup_DF_spline_data with DF_SPLINES defined, interpolate_DF,
# evaluate_HI_Diffusion_correction_2_gamma / _R, evaluate_HI_DI1_DI2_correction; DIFF_CORR_STOREII undefined), Modules/ODEdef_CosmoRec.cpp:198-208
# (applied at the end of fcn_HI_effective when Diffusion_correction_is_on; the DI1 term when DI1_2s_correction_on, which production runmode 0 sets:
# set_startup_data_CR with Diffusion_flag = 1, induced_flag = 2), Modules/global_variables.cpp (Diff_corr_zmin = 500, Diff_corr_zmax = 2000).
# ----------------------------------------------------------------------------------------------------------------------------------------------------

"""Native interpolation of the PDE outputs (`DPesc_splines_Data_DF`): GSL csplines on ascending z of DI1_2s (slot 0), DF_2gamma 3s, 3d, DF_Raman 2s,
zero outside `(zmin, zmax)` with `zmin = max(Diff_corr_zmin f_t, z_last)`, `zmax = min(Diff_corr_zmax f_t, z_first)`."""
struct HIDiffusionFeedback{S}
    DI1::Union{Nothing,S}
    DF_2g::Vector{S}
    DF_R::Vector{S}
    zmin::Float64
    zmax::Float64
end

"""
    hi_diffusion_feedback(z, DI1_2s, DF_2g, DF_R; zmin = 500.0, zmax = 2000.0, f_t = 1.0, DI1_on = true)

Native `setup_DF_interpol_data` from one PDE stage (`z` descending as produced, e.g. `hi_pde_corrections(...)` outputs). `DI1_on` = `DI1_2s_correction_on`.
"""
function hi_diffusion_feedback(z::AbstractVector, DI1_2s::AbstractVector, DF_2g::AbstractVector, DF_R::AbstractVector; zmin = 500.0, zmax = 2000.0, f_t = 1.0, DI1_on::Bool = true)
    za = Float64.(reverse(z))
    sp(v) = natural_cubic_spline(za, reverse(v))
    s2g = [sp(v) for v in DF_2g]; sR = [sp(v) for v in DF_R]
    d = DI1_on ? sp(DI1_2s) : nothing
    S = eltype(s2g)
    return HIDiffusionFeedback{S}(d, s2g, sR, max(zmin * f_t, Float64(z[end])), min(zmax * f_t, Float64(z[1])))
end
hi_diffusion_feedback(out::NamedTuple; kwargs...) = hi_diffusion_feedback(out.z, out.DI1_2s, out.DF_2g, out.DF_R; kwargs...)

"""Native `interpolate_DF(z, k)`: 0 at or outside the range ends, else the spline value."""
_interp_DF(f::HIDiffusionFeedback, s, z) = (z <= f.zmin || z >= f.zmax) ? zero(eltype(s.y)) : spline_eval_native(s, z)
hi_DI1_2s(f::HIDiffusionFeedback, z) = f.DI1 === nothing ? 0.0 : _interp_DF(f, f.DI1, z)
hi_DF_2g(f::HIDiffusionFeedback, i::Int, z) = _interp_DF(f, f.DF_2g[i + 1], z)
hi_DF_R(f::HIDiffusionFeedback, i::Int, z) = _interp_DF(f, f.DF_R[i + 1], z)

"""
    hi_diffusion_rhs!(g, f, z, X)

Add the production diffusion corrections to the native 15-entry `g` (`X` the native state, HI 1s at index 2): two-photon 3s, 3d and Raman 2s
(`DR = X1s DF(z)`, removed from 1s, added to the level) and, when on, the 2s-1s term (`DR = DI1_2s(z) X1s`), in the native order.
"""
function hi_diffusion_rhs!(g, f::HIDiffusionFeedback, z, X)
    X1s = X[EFF_IHI]
    for (i, lvl) in ((0, 3), (1, 5))            # evaluate_HI_Diffusion_correction_2_gamma: (3,0), (3,2)
        DR = X1s * hi_DF_2g(f, i, z)
        g[EFF_IHI] += -DR
        g[EFF_IHI + lvl] += DR
    end
    DR = X1s * hi_DF_R(f, 0, z)                 # evaluate_HI_Diffusion_correction_R: (2,0)
    g[EFF_IHI] += -DR
    g[EFF_IHI + 1] += DR
    if f.DI1 !== nothing                        # evaluate_HI_DI1_DI2_correction
        DR = hi_DI1_2s(f, z) * X1s
        g[EFF_IHI] += -DR
        g[EFF_IHI + 1] += DR
    end
    return g
end
hi_diffusion_rhs!(g, ::Nothing, z, X) = g

"""Element type contributed by a feedback (`Float64` for none)."""
feedback_eltype(::Nothing) = Float64
feedback_eltype(f::HIDiffusionFeedback) = eltype(f.DF_2g[1].y)
