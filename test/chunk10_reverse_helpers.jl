# Chunk 10 reverse-mode helpers: the composed runmode-0 objective with the feedback of each PDE stage passed to the ODE solver as DYNAMIC parameters
# (packed natural-spline coefficients y, b, c, d of DI1_2s, DF_2gamma 3s/3d, DF_Raman 2s, built differentiably from the PDE outputs), a plain-function
# RHS that rebuilds the feedback from the parameter vector, and Mooncake directly through the Rodas5P steps (the Chunk 5e reverse route).
@isdefined(SOLVEP5E) || include("chunk5e_helpers.jl")
@isdefined(SETUP7) || include("chunk7_helpers.jl")

const D10R = HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6)
const NDF10R = 199
const DFZA10R = let o = hi_pde_corrections(HIPDEModel(SETUP7, HIPopulationSplines(NODES5[:, 1:9]; zs = ZS6, ze = ZE6),
                                                     hi_pde_coefficients(NODES5[:, 1:9], HIPopulationSplines(NODES5[:, 1:9]; zs = ZS6, ze = ZE6), ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6), ACC5))
    Float64.(reverse(o.z))                                 # ascending output redshifts (fixed by the PDE step control)
end
const FBZMIN10R = max(500.0, DFZA10R[1]); const FBZMAX10R = min(2000.0, DFZA10R[end])
"""Packed coefficients of the four feedback splines (y, b, c, d for DI1, 2g 3s, 2g 3d, R 2s)."""
pack_feedback(fb::HIDiffusionFeedback) = (sp = [fb.DI1, fb.DF_2g[1], fb.DF_2g[2], fb.DF_R[1]]; reduce(vcat, [vcat(s.y, s.b, s.c, s.d) for s in sp]))
const SPLEN10R = let s = natural_cubic_spline(DFZA10R, zeros(NDF10R)); (length(s.y), length(s.b), length(s.c), length(s.d)) end
function unpack_feedback(pv, off)
    ny, nb, nc, nd = SPLEN10R; L = ny + nb + nc + nd
    sp = map(0:3) do k
        o = off + L * k
        NaturalCubicSpline(DFZA10R, pv[(o + 1):(o + ny)], pv[(o + ny + 1):(o + ny + nb)], pv[(o + ny + nb + 1):(o + ny + nb + nc)], pv[(o + ny + nb + nc + 1):(o + L)])
    end
    return HIDiffusionFeedback{eltype(sp)}(sp[1], [sp[2], sp[3]], [sp[4]], FBZMIN10R, FBZMAX10R)
end
# compact alternative: the 4 x 199 DF NODE VALUES (DI1, 2g 3s, 2g 3d, R 2s, ascending z) as solver parameters; the RHS builds the natural splines
pack_feedback_nodes(fb::HIDiffusionFeedback) = reduce(vcat, [s.y for s in (fb.DI1, fb.DF_2g[1], fb.DF_2g[2], fb.DF_R[1])])
function unpack_feedback_nodes(pv, off)
    n = NDF10R
    sp = [natural_cubic_spline(DFZA10R, pv[(off + n * k + 1):(off + n * (k + 1))]) for k in 0:3]
    return HIDiffusionFeedback{eltype(sp)}(sp[1], [sp[2], sp[3]], [sp[4]], FBZMIN10R, FBZMAX10R)
end
const FEEDBACK_PACKING10R = Ref(:nodes)     # :nodes (796 parameters) or :coefficients (3178)
_unpack10r(pv) = FEEDBACK_PACKING10R[] === :nodes ? unpack_feedback_nodes(pv, 2) : unpack_feedback(pv, 2)
_pack10r(fb) = FEEDBACK_PACKING10R[] === :nodes ? pack_feedback_nodes(fb) : pack_feedback(fb)
const RM10R = RM5
function rhs12_fb!(du, u, pv, z)
    CosmoRec._primal(u[1]) > 0 || (fill!(du, 0); return du)
    rm = RecombinationModel(RM10R.eff, RM10R.cosmos, RM10R.dp_fallback, _unpack10r(pv))
    return recombination_rhs!(du, z, u, rm; flag_He = true, hscale = pv[1], nbscale = pv[2])
end
function rhs7_fb!(du, u, pv, z)
    CosmoRec._primal(u[1]) > 0 || (fill!(du, 0); return du)
    rm = RecombinationModel(RM10R.eff, RM10R.cosmos, RM10R.dp_fallback, _unpack10r(pv))
    return recombination_rhs!(du, z, u, rm; flag_He = false, hscale = pv[1], nbscale = pv[2])
end
const FUN12F = mkfun(rhs12_fb!); const FUN7F = mkfun(rhs7_fb!)
const TOLR10_REF = Ref(1.0e-10)
# EXPLICIT sensitivity algorithm of every solve in the composed objective: without it SciMLSensitivity selects ForwardDiffSensitivity for
# length(u0) + length(p) <= 100 and GaussAdjoint above (concrete_solve.jl:252-258, 367-389), which gave a zero/unstable reverse pass for the
# packed-feedback ODEs (chunk10/probe_mc_isolate.log). Set by the probe outcome (chunk10/probe_mc_sensealg.log).
const SENSEALG10R = Ref{Any}(nothing)
function solve_rev_sa(fun, u0, z0, zn, pv, ab; reltol, dtmax = 10.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(fun, u0, (z0, zn[end]), pv)
    kw = (; reltol = reltol, abstol = ab, saveat = zn, dtmax = dtmax, internalnorm = primal_norm)
    sol = SENSEALG10R[] === nothing ? solve(prob, Rodas5P(); kw...) : solve(prob, Rodas5P(); kw..., sensealg = SENSEALG10R[])
    SciMLBase.successful_retcode(sol) || error("reverse-route ODE solve failed: $(sol.retcode)")
    return reduce(hcat, sol.u)
end
ab10r(n) = abstol5(n; a1 = 1.0e-14, aex = 1.0e-12)
# solvers of one pass: no feedback (pass 0) or with packed feedback; identical routine for the forward (ForwardDiff) and reverse (Mooncake) routes
solver10r(::Nothing) = (f!, u0, z0, zn, pp) -> solve_rev_sa(pp.flag ? FUN12 : FUN7, u0, z0, zn, [pp.hscale, pp.nbscale], ab10r(length(u0)); reltol = TOLR10_REF[])
solver10r(fb) = (f!, u0, z0, zn, pp) -> solve_rev_sa(pp.flag ? FUN12F : FUN7F, u0, z0, zn, vcat([pp.hscale, pp.nbscale], _pack10r(fb)), ab10r(length(u0)); reltol = TOLR10_REF[])
tailsolver10r(hs) = (f!, u0, z0, zn, pt) -> solve_rev_sa(FUNT, u0, z0, zn, [pt.θ[1], pt.θ[2], pt.θ[4], pt.ff, hs], [1.0e-14, 1.0e-14, 1.0e-10]; reltol = TOLR10_REF[])

"""Composed objective pieces with p = [F, A2s1s, hscale, nbscale]; `iterations` diffusion stages (2 = production); returns the final Xe, Te on ZG5E
(or, for `reduced = true`, pass-1 Xe and Te at the SEL5E nodes after ONE stage)."""
function composed10r(p; iterations = 2, reduced = false)
    hs, nb = p[3], p[4]
    kw = (; k_switch = K5E, hscale = hs, nbscale = nb)
    fb = nothing; pass = nothing
    for it in 0:iterations
        rm = fb === nothing ? RM10R : with_diffusion(RM10R, fb)
        last = it == iterations
        pass = recombination_pass(rm, solver10r(fb); kw..., (reduced && last ? (; nodes = SEL5E) : (;))...)
        (reduced && last) && return vcat(pass.Xe, pass.Te)
        last && break
        _, fb = hi_diffusion_stage(RM10R, D10R, pass_solution_rows(pass); hscale = hs, nbscale = nb)
    end
    rm = with_diffusion(RM10R, fb)
    zi = pass.z[end]; y7 = pass.states[end]
    inp = recfast_tail_inputs(rm, zi, y7; hscale = hs, nbscale = nb)
    θeff = [p[1], p[2], THETA5E[3], THETA5E[4] * nb, THETA5E[5], THETA5E[6]]
    tail = recfast_tail(θeff, NATIVE_RECFAST_CONSTANTS, z -> hs * cosmos_H(ACC5, z), zi, inp.Xe_Hi, inp.Xe_Hei, inp.Xei, inp.TMi, inp.dXei, tailsolver10r(hs))
    rows = recombination_output_rows(pass, tail)
    Xe, Te = return_solution_to_grid(rm, rows, ZG5E)
    return vcat(Xe, Te)
end
