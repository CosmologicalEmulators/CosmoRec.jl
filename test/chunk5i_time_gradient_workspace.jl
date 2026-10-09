# Chunk 5i: the caller callback's PREPARED time gradient (`prepare_wstgrad5` / `WSTGrad5`, test/chunk5_helpers.jl; the third workspace increment after the RHS (5g) and
# the state Jacobian (5h)) vs the unprepared oracle `tgrad5!`, which is unchanged. The per-solve callable differentiates ONE Dual-redshift RHS workspace with the mutating
# `ForwardDiff.derivative!` (a DerivativeConfig built from the solve's actual u0/p/z0 types; live u, p, z written into the functor at every call).
# Accuracy gates are unchanged from 5h: ForwardDiff values/nested derivatives bitwise vs the oracle; prepared Mooncake within 1e-12 of ForwardDiff at in-table states.
# Off-table z = 1500 (DP-fallback absorber), Mooncake: the legacy `tgrad5!` AND the prepared callback both ERROR in this configuration, with the same observed
# exception (bounded reproducer inc2/offtable_repro.jl: `MethodError: no method matching increment!!(::Mooncake.RData{value::Float64, partials::RData{values::Tuple{Float64}}},
# ::Float64)`, raised in the pullback of `inner_int_appr` / `patterson_integrate` / `dpesc_coh` (src/DPescCoh.jl:368, 409, 468) whose `DPescIntegrand` holds Dual values).
# The check there is therefore OUTCOME PARITY ONLY (same observed failure class); no reverse vector exists to compare and NO precision cause or accuracy claim is made
# (the cause of this exception is unresolved; it is not attributed to the Float64 absorber Hessian precision issue).
# REQUIRES COSMOREC_NATIVE_DATA_DIR (see test/chunk5_helpers.jl); fails loudly if unset.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(RM5) || include("chunk5_helpers.jl")
@isdefined(SETUP7) || include("chunk7_helpers.jl")

const MC5I = AutoMooncake(; config = nothing)
const ROWS5I = read_fixture5a(joinpath(@__DIR__, "fixtures", "native_ode_rhs.txt"))
Tof5i(u, p, z) = promote_type(eltype(u), typeof(z), typeof(p.hscale), typeof(p.nbscale))
oracle5i(u, p, z) = (T = Tof5i(u, p, z); dT = zeros(T, length(u)); tgrad5!(dT, T.(u), p, z); dT)
function ws5i(u, p, z)                                   # a fresh prepared cache from these actual types (as prob_phase5 builds per solve)
    T = Tof5i(u, p, z); uu = T.(u); g = prepare_wstgrad5(uu, p, z); dT = zeros(T, length(u)); g(dT, uu, p, z)
    @test wstgrad5_fallbacks(g) == 0 && g.ncalls[] == 1
    return dT
end
pofh5i(flag, h) = RecombinationODEParams(RM5, flag, h, one(h))
const SEL5I = [first(r for r in ROWS5I if r.flag == f && r.z == z && r.pat == 0) for (f, z) in ((true, 2500.0), (true, 1500.0), (false, 400.0))]

@testset "Chunk 5i: prepared time-gradient callback vs the unprepared oracle tgrad5!" begin
    @testset "values, nested/mixed ForwardDiff (bitwise), $(r.flag ? 12 : 7)-state z = $(r.z)" for r in SEL5I
        p = RecombinationODEParams(RM5, r.flag); z = r.z; n = length(r.y)
        @test ws5i(r.y, p, z) == oracle5i(r.y, p, z)
        @test ForwardDiff.jacobian(u -> ws5i(u, p, z), r.y) == ForwardDiff.jacobian(u -> oracle5i(u, p, z), r.y)                  # outer Dual STATE
        @test ForwardDiff.derivative(zz -> ws5i(r.y, p, zz), z) == ForwardDiff.derivative(zz -> oracle5i(r.y, p, zz), z)          # outer Dual z (falls to the Dual-z cache)
        @test ForwardDiff.derivative(h -> ws5i(r.y, pofh5i(r.flag, h), z), 1.0) == ForwardDiff.derivative(h -> oracle5i(r.y, pofh5i(r.flag, h), z), 1.0)
        @test ForwardDiff.derivative(a -> ForwardDiff.derivative(zz -> ws5i(r.y, p, zz), a), z) ==
              ForwardDiff.derivative(a -> ForwardDiff.derivative(zz -> oracle5i(r.y, p, zz), a), z)                           # nested second z derivative
        v = (1:n) ./ n
        f_all(g) = a -> g(r.y .+ a .* v, pofh5i(r.flag, 1.0 + 0.5a), z - 3a)                                                  # state, hscale and z all carry the outer Dual
        @test ForwardDiff.derivative(f_all(ws5i), 0.0) == ForwardDiff.derivative(f_all(oracle5i), 0.0)
        # two independent caches, live (u, p, z) alternating; no stale values; return to the original inputs
        y2 = r.y .* (1 .+ 1.0e-6 .* (1:n)); p2 = pofh5i(r.flag, 1.001); z2 = z - 7.0
        ga = prepare_wstgrad5(r.y, p, z); gb = prepare_wstgrad5(r.y, p, z); A = zeros(n); B = zeros(n); C = zeros(n); D = zeros(n)
        ga(A, r.y, p, z); gb(D, y2, p2, z2); ga(B, y2, p2, z2); ga(C, r.y, p, z)
        @test A == C == oracle5i(r.y, p, z) && B == D == oracle5i(y2, p2, z2) && wstgrad5_fallbacks(ga) == 0 && wstgrad5_fallbacks(gb) == 0 && ga.ncalls[] == 3
    end
    @testset "diffusion feedback model (pass-1 feedback), changed states" begin
        S = (a...) -> solve_phase5(a...; reltol = 1.0e-12, a1 = 1.0e-18, aex = 1.0e-14, alg = Rodas5P(), tstops = false)
        p0 = recombination_pass(RM5, S); (_, fb) = hi_diffusion_stage(RM5, HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6), pass_solution_rows(p0))
        RMF = with_diffusion(RM5, fb)
        for k in (1, 500, 1342, 2000, 2900)
            y = p0.states[k]; z = p0.z[k]; p = RecombinationODEParams(RMF, length(y) == 12)
            @test ws5i(y, p, z) == oracle5i(y, p, z)
        end
        y = p0.states[500]; z = p0.z[500]; p = RecombinationODEParams(RMF, true)
        @test ForwardDiff.derivative(h -> ws5i(y, RecombinationODEParams(RMF, true, h, one(h)), z), 1.0) == ForwardDiff.derivative(h -> oracle5i(y, RecombinationODEParams(RMF, true, h, one(h)), z), 1.0)
    end
    @testset "type-mismatched calls: counted oracle fallback, never a conversion" begin
        r = SEL5I[1]; p = RecombinationODEParams(RM5, true); z = r.z; n = length(r.y)
        g = prepare_wstgrad5(r.y, p, z); dT = zeros(n)
        dfb = ForwardDiff.derivative(zz -> (d = zeros(typeof(zz), n); g(d, typeof(zz).(r.y), p, zz); d), z)         # Dual z into a Float64 cache
        @test wstgrad5_fallbacks(g) == 1 && dfb == ForwardDiff.derivative(zz -> oracle5i(r.y, p, zz), z)
        ph = RecombinationODEParams(RM5, true, 1.0f0, 1.0f0); g(dT, r.y, ph, z)                                          # p of another type (Float32 scales)
        @test wstgrad5_fallbacks(g) == 2 && dT == oracle5i(r.y, ph, z)
        g(dT, BigFloat.(r.y), p, z)                                                                                    # state of another element type (BigFloat into Float64 cache)
        @test wstgrad5_fallbacks(g) == 3
        g(dT, r.y, p, z); @test wstgrad5_fallbacks(g) == 3 && dT == oracle5i(r.y, p, z)                                # the cache still serves the prepared types
    end
    @testset "non-positive primal state: the same domain behaviour as the oracle (no clamp)" begin
        r = SEL5I[1]; p = RecombinationODEParams(RM5, true); n = length(r.y); bad = copy(r.y); bad[3] = -abs(bad[3])
        out(f) = try; ("ok", f()); catch e; ("err", typeof(e), first(sprint(showerror, e), 300)); end
        @test out(() -> ws5i(bad, p, r.z)) == out(() -> oracle5i(bad, p, r.z))
    end
    @testset "family tag: default check passes, ordering pinned (internal tagcount)" begin
        outer = ForwardDiff.Tag(identity, Float64); inner = buftg5_tag(ForwardDiff.Dual{typeof(outer),Float64,1})
        @test ForwardDiff.:≺(typeof(outer), typeof(inner)) && !ForwardDiff.:≺(typeof(inner), typeof(outer))
        r = SEL5I[1]; g = prepare_wstgrad5(r.y, RecombinationODEParams(RM5, true), r.z)
        @test ForwardDiff.checktag(g.cfg, g.f, r.z)                                  # the default check derivative! performs
        @test typeof(g.cfg).parameters[1] === ForwardDiff.Tag{BufTG5,Float64}
    end
    @testset "high precision: BigFloat(128) new == legacy, directional second derivative, worst component z = 1500" begin
        r = SEL5I[2]; p = RecombinationODEParams(RM5, true); n = 12; W1 = randn(MersenneTwister(3), n); ek = [i == 10 ? 1 : 0 for i in 1:n]
        dn, do_, vn, vo = setprecision(BigFloat, 128) do
            yb = BigFloat.(r.y)
            (ForwardDiff.derivative(t -> dot(W1, ws5i(yb .+ t .* ek, p, r.z)), big(0.0)), ForwardDiff.derivative(t -> dot(W1, oracle5i(yb .+ t .* ek, p, r.z)), big(0.0)), ws5i(yb, p, r.z), oracle5i(yb, p, r.z))
        end
        @test dn == do_ && vn == vo
    end
    @testset "prepared DI-Mooncake, ONE tape with TWO time-gradient evaluations through one prebuilt cache, $(r.flag ? 12 : 7)-state z = $(r.z)" for r in SEL5I
        p = RecombinationODEParams(RM5, r.flag); z = r.z; n = length(r.y); W1 = randn(MersenneTwister(3), n); W2 = randn(MersenneTwister(4), n); cs = 1 .+ 1.0e-6 .* (1:n)
        pre = prepare_wstgrad5(r.y, p, z)
        Tw(u) = (d = zeros(eltype(u), n); pre(d, u, p, z); d); To(u) = (d = zeros(eltype(u), n); tgrad5!(d, u, p, z); d)
        fw(u) = dot(W1, Tw(u)) + dot(W2, Tw(u .* cs)); fo(u) = dot(W1, To(u)) + dot(W2, To(u .* cs))
        @test fw(r.y) == fo(r.y)
        outc(f) = try; ("ok", f()); catch e; ("err", string(typeof(e))); end
        go = outc(() -> (po = prepare_gradient(fo, MC5I, r.y); gradient(fo, po, MC5I, r.y)))
        if z == 1500.0
            # off-table: outcome parity with the legacy oracle (currently BOTH error, same exception; a vector is compared bitwise only if both ever succeed): no accuracy/precision claim
            gn = outc(() -> (pw = prepare_gradient(fw, MC5I, r.y); gradient(fw, pw, MC5I, r.y)))
            println("Chunk5i off-table Mooncake: legacy ", go[1], " new ", gn[1])
            @test gn == go
        else
            pw = prepare_gradient(fw, MC5I, r.y); g1 = gradient(fw, pw, MC5I, r.y); g2 = gradient(fw, pw, MC5I, r.y)
            @test g1 == g2 && wstgrad5_fallbacks(pre) == 0
            ref = ForwardDiff.gradient(fo, r.y)
            @test maximum(abs.(g1 .- ref)) / maximum(abs, ref) < 1e-12        # accuracy gate, unchanged from 5h
            println("Chunk5i Mooncake in-table rel err new ", maximum(abs.(g1 .- ref)) / maximum(abs, ref), " legacy ", go[1] == "ok" ? maximum(abs.(go[2] .- ref)) / maximum(abs, ref) : "legacy-err")
        end
    end
    @testset "allocations per time-gradient call: prepared far below tgrad5!" begin
        for r in SEL5I[[1, 3]]
            p = RecombinationODEParams(RM5, r.flag); n = length(r.y); dT = zeros(n); g = prepare_wstgrad5(r.y, p, r.z)
            an() = @allocated g(dT, r.y, p, r.z); ao() = @allocated tgrad5!(dT, r.y, p, r.z)
            an(); ao(); an(); ao()
            a_n = an(); a_o = ao(); println("Chunk5i allocated bytes per time-gradient call $(n)-state: prepared $a_n, tgrad5! $a_o")
            @test a_n < a_o && a_n <= 256
        end
    end
    @testset "actual solve: prob_phase5/odefun5 use the prepared callback, cache used, 0 fallbacks, solution and solver statistics bitwise vs the oracle-tgrad problem" begin
        nblocks = Ref(0); ncalls = Ref(0); nfb = Ref(0); allsame = Ref(true)
        function S(f!, u0, z0, znodes, p)
            prob = prob_phase5(u0, z0, znodes, p)
            allsame[] &= prob.f.tgrad isa WSTGrad5
            popt = ODEFunction{true, SciMLBase.FullSpecialize}(prob.f.f; jac = prob.f.jac, tgrad = tgrad5!)               # identical problem except the unprepared oracle tgrad
            prob_o = ODEProblem{true, SciMLBase.FullSpecialize}(popt, u0, prob.tspan, p)
            kw = (; reltol = 1.0e-10, abstol = abstol5(length(u0); a1 = 1.0e-18, aex = 1.0e-12), saveat = znodes, internalnorm = primal_norm, tstops = znodes)
            sn = solve(prob, Rodas4P(); kw...); so = solve(prob_o, Rodas4P(); kw...)
            allsame[] &= sn.u == so.u && sn.t == so.t && all(s -> getfield(sn.stats, s) == getfield(so.stats, s), (:nf, :nreject, :naccept, :nsolve, :nw)) && sn.stats.njacs == so.stats.njacs
            nblocks[] += 1; ncalls[] += prob.f.tgrad.ncalls[]; nfb[] += wstgrad5_fallbacks(prob.f.tgrad)
            return reduce(hcat, sn.u)
        end
        recombination_pass(RM5, S)
        @test allsame[] && nblocks[] > 10 && ncalls[] > 0 && nfb[] == 0
        println("Chunk5i real pass: blocks ", nblocks[], " prepared time-gradient calls ", ncalls[], " fallbacks ", nfb[])
    end
end
