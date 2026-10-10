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

Keyword `reduced_profiles` (default `false`, a `Bool`): see `hi_pde_corrections`; skips unread `hp_phi` entries inside the closed PDE stage, outputs bitwise unchanged.

Keyword `fitwork` (opt-in, requires `quadplan`; default `nothing` = the allocating fit, unchanged): `true` builds one stage-owned `CosmoRec.HIQuadWork` (reusable spline-fit scratch, ~164 kB) for the stage's ordinate
type; a caller-owned `HIQuadWork` is used as is (one per task/chain: it is mutated by every call and must not be shared by concurrent calls; plans may be shared). Results are bitwise those of the plan path.
The scratch fit's spline is a borrowed, ephemeral view valid only inside one integral. A work object of another ordinate type, or a scratch of another size, is not coerced: those fits use the allocating oracle and
are counted (`nmismatch`, `nfallback`). Reverse mode: supported through Mooncake (the scratch method is overlaid so later reads of the scratch are differentiated); there is no ChainRulesCore rule for it.
"""
function hi_diffusion_stage(rm::RecombinationModel, d::HIDiffusionInputs, rows::AbstractMatrix; hscale = 1.0, nbscale = 1.0, levels = nothing,
                            quadplan = nothing, fitwork = nothing, reduced_profiles = false, T::Type = promote_type(eltype(rows), typeof(hscale), typeof(nbscale), Float64))
    # always the scaled background (multiplying by 1.0 is exact): a value/type branch here would drop the explicit hscale/nbscale dependence
    # of the PDE stage for AD systems that see plain Float64 parameters (Mooncake), as found in chunk10/probe_link_pde.log
    bg = ScaledBackground(rm.cosmos, hscale, nbscale)
    pops = HIPopulationSplines(rows; zs = d.zs, ze = d.ze)
    coef = hi_pde_coefficients(rows, pops, bg, d.htable, d.lnBitot, d.levels; zs = d.zs, ze = d.ze)
    out = hi_pde_corrections(HIPDEModel(d.setup, pops, coef, bg); zs = d.zs, ze = d.ze, T = T, levels = levels, quadplan = quadplan, fitwork = fitwork, reduced_profiles = reduced_profiles)
    return out, hi_diffusion_feedback(out)
end

"""
    recombination_history_diffusion(rm, θ, d, solve_pass, solve_tail, zgrid; kwargs...)

The production (runmode 0) history: pass 0 without corrections, then `d.iterations` times a diffusion stage from the previous pass followed by a pass with
the feedback switched on; the last pass is completed by the Recfast tail and the output assembly (`recombination_history`). Returns the final history
plus `passes` (all passes, the last one inside `final`) and `stages` (PDE outputs). `patterson_levels` (a vector of `PattersonLevels`, one per stage)
records or replays the quadrature stopping decisions (frozen-branch verification).

Keywords `prefix_reuse` (a `CosmoRec.PrefixReuse` declaration; default `nothing` = off) and `prefix_stats` (a caller-owned `CosmoRec.PrefixReuseStats`): opt-in reuse of the feedback-free source pass's complete
12-state ODE blocks lying entirely above the actual feedback zero domain by the later passes of this ONE call (cache allocated inside the call; returns exactly the same fields). Requires the explicit declarations
documented at `PrefixReuse`; anything unproven (iterations 0, `nodes`, other solver object, signature mismatch, unknown feedback type, key/type mismatch, query excursion) recomputes and is counted in the stats.

Keyword `reduced_profiles` (default `false`, a `Bool`): see `hi_pde_corrections`; skips unread `hp_phi` entries inside the closed PDE stage, outputs bitwise unchanged.

Keyword `fitwork` (opt-in, requires `quadplan`; default `nothing` = the allocating fit, unchanged): `true` builds one stage-owned `CosmoRec.HIQuadWork` (reusable spline-fit scratch, ~164 kB) for the stage's ordinate
type; a caller-owned `HIQuadWork` is used as is (one per task/chain: it is mutated by every call and must not be shared by concurrent calls; plans may be shared). Results are bitwise those of the plan path.
The scratch fit's spline is a borrowed, ephemeral view valid only inside one integral. A work object of another ordinate type, or a scratch of another size, is not coerced: those fits use the allocating oracle and
are counted (`nmismatch`, `nfallback`). Reverse mode: supported through Mooncake (the scratch method is overlaid so later reads of the scratch are differentiated); there is no ChainRulesCore rule for it.
"""
function recombination_history_diffusion(rm::RecombinationModel, θ, d::HIDiffusionInputs, solve_pass, solve_tail, zgrid::AbstractVector; hscale = 1.0, nbscale = 1.0,
                                         patterson_levels = nothing, quadplan = nothing, fitwork = nothing, reduced_profiles = false,
                                         prefix_reuse = nothing, prefix_stats = nothing, kwargs...)
    passes = Any[]; stages = Any[]
    m = rm
    kw = (; hscale = hscale, nbscale = nbscale, kwargs...)
    pplan = prefix_reuse === nothing ? (prefix_stats === nothing || _reset!(prefix_stats); nothing) : _prefix_begin(prefix_reuse, prefix_stats, solve_pass, d, rm, kwargs)
    pcache = pplan === nothing ? nothing : PrefixCache()
    for it in 0:(d.iterations - 1)
        kw_it = pplan === nothing ? kw : merge(kw, (; prefix = _prefix_pass_binding(pplan, prefix_stats, pcache, it, rm, m)))
        pass = recombination_pass(m, solve_pass; kw_it...)
        push!(passes, pass)
        out, fb = hi_diffusion_stage(rm, d, pass_solution_rows(pass); hscale = hscale, nbscale = nbscale,
                                     levels = patterson_levels === nothing ? nothing : patterson_levels[it + 1], quadplan = quadplan, fitwork = fitwork, reduced_profiles = reduced_profiles)
        push!(stages, out)
        m = with_diffusion(rm, fb)
    end
    kw_fin = pplan === nothing ? kw : merge(kw, (; prefix = _prefix_pass_binding(pplan, prefix_stats, pcache, d.iterations, rm, m)))
    final = recombination_history(m, θ, solve_pass, solve_tail, zgrid; kw_fin...)
    push!(passes, final.pass)
    return (final = final, passes = passes, stages = stages, Xe = final.Xe, Te = final.Te, z = zgrid)
end
