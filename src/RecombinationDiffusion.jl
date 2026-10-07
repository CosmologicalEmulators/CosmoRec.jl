# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Chunk 10a: the production CosmoRec iteration (runmode 0): ODE pass, then Diff_iteration_max = 2 times [HI diffusion PDE stage from the previous pass,
# feedback, ODE pass], the Recfast tail and the output assembly of the LAST pass.
# Source (ORIGINAL CosmoRec v3.0b, read-only): CosmoRec.cpp:340-478 (CosmoRec()), :395-408 (PDE call: zs = min(Diff_corr_zmax f_t 1.25, zstart/1.002),
# ze = max(Diff_corr_zmin f_t, zend 1.002), nS_2gamma, nS_Raman, the populations of the previous pass), Modules/main.CosmoRec.cpp:40-80
# (pass_on_the_Solution_CosmoRec stored for iterations < Diff_iteration_max, the output rows only for the last iteration).
# ----------------------------------------------------------------------------------------------------------------------------------------------------

"""Cosmos accessors with `H(z)` and `NH(z)` multiplied by `hscale`, `nbscale` (the Chunk 5e differentiable background parameters), so that the PDE stage
sees the same background as the ODE passes (`ode_background`). Only the accessors used by the PDE stage are provided."""
struct ScaledBackground{A,T1,T2}
    a::A
    hscale::T1
    nbscale::T2
end
cosmos_H(s::ScaledBackground, z) = s.hscale * cosmos_H(s.a, z)
cosmos_NH(s::ScaledBackground, z) = s.nbscale * cosmos_NH(s.a, z)
cosmos_TCMB(s::ScaledBackground, z) = cosmos_TCMB(s.a, z)

"""Native `pass_on_the_Solution_CosmoRec` rows `[z, Xe, X1s, 2s, 2p, 3s, 3p, 3d, rho]` of a `recombination_pass` (trailing zero columns dropped)."""
pass_solution_rows(pass) = reduce(vcat, [[pass.z[m] pass.observables[m]'] for m in eachindex(pass.z)])

"""Explicit inputs of the HI diffusion stages: the armed PDE setup (7a), the HI rate table and `ln Bitot` table of the PDE coefficients (6a), the resolved
levels, the PDE range and the number of diffusion iterations (`Diff_iteration_max`)."""
struct HIDiffusionInputs{T,B}
    setup::HIPDESetup
    htable::T
    lnBitot::B
    levels::HIPDELevels
    zs::Float64
    ze::Float64
    iterations::Int
end
HIDiffusionInputs(setup, htable, lnBitot; levels = NATIVE_HI_PDE_LEVELS, zstart = 3000.0, zend = 50.0, f_t = 1.0, Diff_corr_zmin = 500.0, Diff_corr_zmax = 2000.0, iterations = 2) =
    HIDiffusionInputs(setup, htable, lnBitot, levels, min(Diff_corr_zmax * f_t * 1.25, zstart / 1.002), max(Diff_corr_zmin * f_t, zend * 1.002), iterations)

"""
    hi_diffusion_stage(rm, d, rows) -> (out, feedback)

One HI diffusion stage of `CosmoRec()`: population and pd/Dnem splines of the previous pass `rows`, the PDE march with its correction integrals
(`hi_pde_corrections`) and the feedback interpolation (`hi_diffusion_feedback`).
"""
function hi_diffusion_stage(rm::RecombinationModel, d::HIDiffusionInputs, rows::AbstractMatrix; hscale = 1.0, nbscale = 1.0, levels = nothing,
                            quadplan = nothing, T::Type = promote_type(eltype(rows), typeof(hscale), typeof(nbscale), Float64))
    # always the scaled background (multiplying by 1.0 is exact): a value/type branch here would drop the explicit hscale/nbscale dependence
    # of the PDE stage for AD systems that see plain Float64 parameters (Mooncake), as found in chunk10/probe_link_pde.log
    bg = ScaledBackground(rm.cosmos, hscale, nbscale)
    pops = HIPopulationSplines(rows; zs = d.zs, ze = d.ze)
    coef = hi_pde_coefficients(rows, pops, bg, d.htable, d.lnBitot, d.levels; zs = d.zs, ze = d.ze)
    out = hi_pde_corrections(HIPDEModel(d.setup, pops, coef, bg); zs = d.zs, ze = d.ze, T = T, levels = levels, quadplan = quadplan)
    return out, hi_diffusion_feedback(out)
end

"""
    recombination_history_diffusion(rm, θ, d, solve_pass, solve_tail, zgrid; kwargs...)

The production (runmode 0) history: pass 0 without corrections, then `d.iterations` times a diffusion stage from the previous pass followed by a pass with
the feedback switched on; the last pass is completed by the Recfast tail and the output assembly (`recombination_history`). Returns the final history
plus `passes` (all passes, the last one inside `final`) and `stages` (PDE outputs). `patterson_levels` (a vector of `PattersonLevels`, one per stage)
records or replays the quadrature stopping decisions (frozen-branch verification).
"""
function recombination_history_diffusion(rm::RecombinationModel, θ, d::HIDiffusionInputs, solve_pass, solve_tail, zgrid::AbstractVector; hscale = 1.0, nbscale = 1.0,
                                         patterson_levels = nothing, quadplan = nothing, kwargs...)
    passes = Any[]; stages = Any[]
    m = rm
    kw = (; hscale = hscale, nbscale = nbscale, kwargs...)
    for it in 0:(d.iterations - 1)
        pass = recombination_pass(m, solve_pass; kw...)
        push!(passes, pass)
        out, fb = hi_diffusion_stage(rm, d, pass_solution_rows(pass); hscale = hscale, nbscale = nbscale,
                                     levels = patterson_levels === nothing ? nothing : patterson_levels[it + 1], quadplan = quadplan)
        push!(stages, out)
        m = with_diffusion(rm, fb)
    end
    final = recombination_history(m, θ, solve_pass, solve_tail, zgrid; kw...)
    push!(passes, final.pass)
    return (final = final, passes = passes, stages = stages, Xe = final.Xe, Te = final.Te, z = zgrid)
end
