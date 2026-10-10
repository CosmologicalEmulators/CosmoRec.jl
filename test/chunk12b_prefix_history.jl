# Chunk 12b: opt-in prefix reuse on the REAL public `recombination_history_diffusion` (c0 physics, native fixtures): OFF == ON bitwise for every returned field (Xe, Te, passes, stages), changed
# cosmology constants, iterations 0 / 1 / 3, declaration refusals / fallbacks with truthful counts, and the composed ForwardDiff Jacobian (all parameters) at the fiducial and at
# hscale = 1.0001 / nbscale = 0.9999, ON == OFF bitwise (values and partials). Real reverse-mode gates run in fresh worker processes (prefix/rev_worker.jl, one configuration at a time).
# REQUIRES COSMOREC_NATIVE_DATA_DIR (see test/chunk7_helpers.jl).
using Test
using ForwardDiff
using CosmoRec
using CosmoRec: PrefixReuse, PrefixReuseStats
@isdefined(SETUP7) || include("chunk7_helpers.jl")
@isdefined(SOLVEP5E) || include("chunk5e_helpers.jl")

const CNT12B = Ref(0)
const SOLVE12B = (a...) -> (CNT12B[] += 1; SOLVEP5E(a...))                     # declared pure in its supplied arguments; the counter is instrumentation
const OTHER12B = (a...) -> SOLVEP5E(a...)
bitsv12b(x::Float64) = reinterpret(UInt64, x)
bitsv12b(x::ForwardDiff.Dual) = (bitsv12b(ForwardDiff.value(x)), map(bitsv12b, ForwardDiff.partials(x).values))
bitsv12b(x::AbstractArray) = map(bitsv12b, x)
bitsv12b(x::Tuple) = map(bitsv12b, x)
bitsv12b(x::NamedTuple) = map(bitsv12b, values(x))
bitsv12b(x::Union{Integer,Bool,Symbol}) = x
bitsv12b(x) = isstructtype(typeof(x)) ? ntuple(i -> bitsv12b(getfield(x, i)), fieldcount(typeof(x))) : x
plan12b(sigs; excursion = 0.0, solve = SOLVE12B) = PrefixReuse(solve; signatures = sigs, max_query_excursion = excursion, pure_solver = true)
function hist12b(rm, d; plan = nothing, stats = PrefixReuseStats(), hs = 1.0, nb = 1.0, solve = SOLVE12B, kwargs...)
    CNT12B[] = 0
    kw = plan === nothing ? (;) : (; prefix_reuse = plan, prefix_stats = stats)
    return recombination_history_diffusion(rm, THETA5E_vec(P0_5E), d, solve, solve_tail5, ZG5E; hscale = hs, nbscale = nb, k_switch = K5E, kw..., kwargs...)
end
function cosmo12b(rm; nb = 1.03, t0 = 1.001)           # changed cosmology constants (baryon density, CMB temperature) with the same tables: a different `cosmos` object and different y0 / rates
    c = rm.cosmos.c
    c2 = CosmoRec.CosmosConstants(c.Y_p, c.T0 * t0, c.fac_mHemH, c.zsRe, c.Nb0 * nb, c.H0, c.O_k, c.O_L, c.O_m, c.O_rel, c.rho_g_gr, c.sigT0, c.cl, c.hPlanck, c.kBoltz, c.mElect, c.EH_inf_ergs, c.me_mp, c.me_malp, c.EionHeI, c.pi, c.z_saha)
    return RecombinationModel(rm.eff, CosmoRec.CosmosAccessors(c2, rm.cosmos.sp, rm.cosmos.hubble), rm.dp_fallback, rm.diffusion)
end

@testset "Chunk 12b: prefix reuse on the real recombination_history_diffusion" begin
    D2 = HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6)
    @testset "OFF == ON bitwise, truthful counts (fiducial and changed cosmology)" for (lab, rm, hs, nb) in (("fiducial", RM5, 1.0, 1.0), ("scales", RM5, 1.0001, 0.9999), ("cosmology", cosmo12b(RM5), 1.0, 1.0))
        off = hist12b(rm, D2; hs = hs, nb = nb); noff = CNT12B[]
        st = PrefixReuseStats(); on = hist12b(rm, D2; plan = plan12b(fill(:rodas5p, 3)), stats = st, hs = hs, nb = nb); non = CNT12B[]
        @test keys(off) == keys(on) && bitsv12b(off) == bitsv12b(on)                         # every field: final history, passes (states, observables, k_switch...), stages (PDE outputs), Xe, Te
        @test st.status === :active && st.source_blocks > 0 && st.reused_blocks > 0 && st.reused_blocks + st.recomputed_blocks == st.consumer_blocks
        @test noff - non == st.reused_blocks && non > 0                                         # exactly the skipped solve calls; none fabricated
        @test st.refused_excursion > 0                                                         # crossing / later blocks are recomputed
    end
    @test bitsv12b(hist12b(RM5, D2).Xe) != bitsv12b(hist12b(cosmo12b(RM5), D2).Xe)             # the changed-cosmology case is really different physics
    @testset "iterations 0 / 1 / 3" for it in (0, 1, 3)
        d = HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6; iterations = it)
        off = hist12b(RM5, d); st = PrefixReuseStats(); on = hist12b(RM5, d; plan = plan12b(fill(:s, it + 1)), stats = st)
        @test bitsv12b(off) == bitsv12b(on) && length(on.passes) == it + 1
        it == 0 ? (@test st.status === :fallback && st.reason === :no_consumer_pass && st.reused_blocks == 0) : (@test st.status === :active && st.reused_blocks > 0)
    end
    @testset "declarations: refusals and fallbacks with no unsafe skip" begin
        off = hist12b(RM5, D2)
        @test_throws ArgumentError hist12b(RM5, D2; plan = plan12b(fill(:s, 2)))                    # wrong number of per-pass declarations
        st = PrefixReuseStats(); on = hist12b(RM5, D2; plan = plan12b([:s, :s, :t]), stats = st)     # final pass declared differently -> its blocks recomputed
        @test bitsv12b(off) == bitsv12b(on) && st.refused_pass == 1 && st.reused_blocks > 0
        st = PrefixReuseStats(); on = hist12b(RM5, D2; plan = plan12b([:s, :t, :u]), stats = st); @test bitsv12b(off) == bitsv12b(on) && st.refused_pass == 2 && st.reused_blocks == 0
        st = PrefixReuseStats(); on = hist12b(RM5, D2; plan = plan12b(fill(:s, 3); excursion = 1.0e5), stats = st); @test bitsv12b(off) == bitsv12b(on) && st.reused_blocks == 0 && st.refused_excursion > 0
        st = PrefixReuseStats(); on = hist12b(RM5, D2; plan = plan12b(fill(:s, 3); solve = OTHER12B), stats = st)         # plan bound to a different solve object
        @test bitsv12b(off) == bitsv12b(on) && st.status === :fallback && st.reason === :solver_binding && st.reused_blocks == 0
        # a `nodes` override cannot feed the PDE stage even without reuse; the declaration check itself must fall back before any skip
        stn = PrefixReuseStats(); @test CosmoRec._prefix_begin(plan12b(fill(:s, 3)), stn, SOLVE12B, D2, RM5, (; nodes = [1000, 2000])) === nothing && stn.status === :fallback && stn.reason === :nodes_override
        stf = PrefixReuseStats(); onf = hist12b(with_diffusion(RM5, nothing), D2; plan = plan12b(fill(:s, 3)), stats = stf); @test stf.status === :active
        # a source model that already has feedback is refused
        fb = hi_diffusion_stage(RM5, D2, NODES5[:, 1:9])[2]            # a real feedback object from the native node rows
        stg = PrefixReuseStats(); ong = hist12b(with_diffusion(RM5, fb), D2; plan = plan12b(fill(:s, 3)), stats = stg); @test stg.status === :fallback && stg.reason === :source_has_feedback && stg.reused_blocks == 0
        # stats are reset on every call and never shared
        st = PrefixReuseStats(); p = plan12b(fill(:s, 3)); hist12b(RM5, D2; plan = p, stats = st); r1 = st.reused_blocks; hist12b(RM5, D2; plan = p, stats = st); @test st.reused_blocks == r1
        @test PrefixReuseStats().status === :off
    end
    @testset "composed ForwardDiff Jacobian (18 outputs x 4 parameters): ON == OFF bitwise, fiducial and hs/nb != 1" begin
        function out2(p; on)
            kw = on ? (; prefix_reuse = plan12b(fill(:s, 3)), prefix_stats = PrefixReuseStats()) : (;)
            r = recombination_history_diffusion(RM5, THETA5E_vec(p), HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6), SOLVE12B, solve_tail5, ZG5E; hscale = p[3], nbscale = p[4], k_switch = K5E, kw...)
            return vcat(r.Xe, r.Te)
        end
        for p in (P0_5E, [P0_5E[1], P0_5E[2], 1.0001, 0.9999])
            f_off = out2(p; on = false); f_on = out2(p; on = true); @test bitsv12b(f_off) == bitsv12b(f_on)
            Jo = ForwardDiff.jacobian(q -> out2(q; on = false), p); Jn = ForwardDiff.jacobian(q -> out2(q; on = true), p)
            @test size(Jn) == (18, 4) && all(isfinite, Jn) && bitsv12b(Jo) == bitsv12b(Jn)
        end
    end
end
