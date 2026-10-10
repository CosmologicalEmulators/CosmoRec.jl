# Chunk 11d: caller-owned spline-fit scratch of the opt-in quadrature plan path (`HIFitScratch`, `HIQuadWork`, `fitwork` keyword; src/HIPDEQuadPlan.jl).
# The allocating `_hi_fit_cached` and the generic `natural_cubic_spline` stay the ORACLES. Gates: (1) scratch fit bitwise equal to both for Float64/32/16, Dual, nested Dual,
# BigFloat 128/256 on synthetic irregular (n = 3 included) and the native slices, identical errors, poisoned scratch fully overwritten, live-F changes, alias/size/type fallbacks (counted,
# result still the oracle's, F never clobbered); (2) `_hi_spline_parts` with scratch vs the verbatim previous body (forced 1:8, adaptive, empty/reversed intervals); (3) PUBLIC stage path
# (`hi_diffusion_stage(...; quadplan, fitwork)`) bitwise vs the plan path incl. record/replay, ForwardDiff, caller-owned reuse, concurrent independent scratches; (4) reverse: CRC rule and
# prepared DI-Mooncake with one scratch used twice on one tape, vs the allocating path and ForwardDiff. REQUIRES COSMOREC_NATIVE_DATA_DIR (see test/chunk7_helpers.jl).
using Test
using Random
using LinearAlgebra
using ForwardDiff
using ChainRulesCore
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
using CosmoRec: HIQuadPlan, HIQuadPlans, HIFitScratch, HIQuadWork, _hi_quad_plan, _hi_quad_plans, _hi_spline_parts, _hi_fit_cached, _hi_fit_scratched, NATIVE_HI_PDE_ATOM
@isdefined(SETUP7) || include("chunk7_helpers.jl")

const MC11D = AutoMooncake(; config = nothing)
bits11d(x::Float64) = reinterpret(UInt64, x)
bits11d(x::Float32) = reinterpret(UInt32, x)
bits11d(x::Float16) = reinterpret(UInt16, x)
bits11d(x::ForwardDiff.Dual) = (bits11d(ForwardDiff.value(x)), map(bits11d, ForwardDiff.partials(x).values))
bits11d(x::BigFloat) = string(x)
bits11d(x::AbstractArray) = map(bits11d, x)
bits11d(x::Tuple) = map(bits11d, x)
same11d(a, b) = typeof(a) == typeof(b) && bits11d(a.x) == bits11d(b.x) && bits11d(a.y) == bits11d(b.y) && bits11d(a.b) == bits11d(b.b) && bits11d(a.c) == bits11d(b.c) && bits11d(a.d) == bits11d(b.d)
# verbatim previous body of `_hi_spline_parts` (natural_cubic_spline refit per call): the reference
function reference_parts11d(q::HIQuadPlan{K}, F::AbstractVector, forced::NTuple{K,Int}) where {K}
    s = natural_cubic_spline(q.knots, F)
    T = eltype(s.y); f = t -> spline_eval_native(s, t); epsabs = q.epsabs0; r = zero(T)
    parts = ntuple(_ -> zero(T), Val(K)); levs = ntuple(_ -> 0, Val(K)); oks = ntuple(_ -> true, Val(K))
    for i in 1:K
        a = q.a[i]; b = q.b[i]
        if a >= b
            ri = zero(T); li = 0; oki = true
        elseif forced[i] > 0
            ri, _ = patterson_integrate(f, a, b, CosmoRec.HI_EPSREL_DF, epsabs; level = forced[i]); li = forced[i]; oki = true
        else
            ri, li = patterson_integrate(f, a, b, CosmoRec.HI_EPSREL_DF, epsabs); oki = true
            if li >= 8
                r7, _ = patterson_integrate(f, a, b, CosmoRec.HI_EPSREL_DF, epsabs; level = 7)
                oki = abs(CosmoRec._primal(ri) - CosmoRec._primal(r7)) <= max(epsabs, abs(CosmoRec._primal(ri)) * CosmoRec.HI_EPSREL_DF)
            end
        end
        parts = Base.setindex(parts, ri, i); levs = Base.setindex(levs, li, i); oks = Base.setindex(oks, oki, i)
        r += ri; epsabs = abs(CosmoRec._primal(r)) * CosmoRec.HI_EPSREL_DF
    end
    return parts, levs, oks
end
function synth_plans11d()
    rng = Random.Xoshiro(1104); out = HIQuadPlan[]
    for n in (3, 7, 40, 300)
        h = exp.(randn(rng, n - 1) .* 4.0); x = cumsum(vcat(0.3, h ./ sum(h) .* 2.0)); lo = x[1]; hi = x[end]; w = hi - lo
        push!(out, _hi_quad_plan(x, [lo, lo + 0.13w, lo + 0.13w, lo + 0.5w, lo + 0.31w, hi], 1.0e-60))      # includes an empty and a reversed sub-interval
    end
    return out
end
function with_knots11d(q::HIQuadPlan{K}, knots) where {K}      # plan whose stored factors are those of q but whose ACTUAL knots differ (stale / invalid knots)
    return HIQuadPlan{K}(knots, q.a, q.b, q.epsabs0, q.dx, q.g, q.off, q.alpha, q.t, q.ell)
end
errinfo11d(f) = try f(); nothing catch e; (typeof(e), sprint(showerror, e)) end
function poison!(w::HIFitScratch{T}, v) where {T}
    for a in (w.y, w.c, w.z, w.b, w.d); fill!(a, convert(T, v)); end
    return w
end
# ordinate families: (label, builder F(rng, n) -> Vector{T})
function families11d(F, v)
    D(F, v) = [ForwardDiff.Dual(F[i], v[i]) for i in eachindex(F)]
    Tg = ForwardDiff.Tag(:o, Float64); Ti = ForwardDiff.Tag(:i, ForwardDiff.Dual{typeof(Tg),Float64,1})
    return [("Float64", F), ("Float32", Float32.(F)), ("Float16", Float16.(F ./ maximum(abs, F))), ("Dual", D(F, v)),
            ("nestedDual", [ForwardDiff.Dual{typeof(Ti)}(ForwardDiff.Dual{typeof(Tg)}(F[i], v[i]), ForwardDiff.Dual{typeof(Tg)}(v[i], 0.5)) for i in eachindex(F)])]
end

@testset "Chunk 11d: caller-owned spline-fit scratch of the quadrature plan" begin
    rng = Random.Xoshiro(20261110)
    nat = _hi_quad_plans(SETUP7)
    plans = vcat(synth_plans11d(), [nat.di1, nat.two3, nat.raman])
    @testset "scratch fit == _hi_fit_cached == natural_cubic_spline, all types, n = $(length(q.knots))" for q in plans
        n = length(q.knots); x = q.knots
        F = (1.0 .+ 0.3 .* sin.(40.0 .* (x .- x[1]) ./ (x[end] - x[1]))) .* exp.(-3.0 .* (x .- x[1])) .* (1.0 .+ 0.05 .* randn(rng, n)); v = randn(rng, n)
        for (lab, Fx) in families11d(F, v)
            T = eltype(Fx); w = HIFitScratch{T}(n); Fbits = bits11d(Fx); Fcopy = copy(Fx)         # bits of F saved BEFORE any fit
            poison!(w, 123.0)       # stale garbage must not survive
            s = _hi_fit_scratched(w, q, Fx)
            @test same11d(s, _hi_fit_cached(q, Fx)) && same11d(s, natural_cubic_spline(q.knots, Fx))
            @test w.nfallback[] == 0
            @test s.y === w.y && s.b === w.b && s.c === w.c && s.d === w.d && s.x === q.knots        # documented ephemeral borrowed result
            # live F changes: a second, different F fully replaces the first; returning to the original reproduces it
            F2 = copy(Fx); F2 .= F2 .* convert(T, 1.5); F2[1] += convert(T, 0.25)
            s2 = _hi_fit_scratched(w, q, F2); @test same11d(s2, _hi_fit_cached(q, F2))
            s3 = _hi_fit_scratched(w, q, Fx); @test same11d(s3, _hi_fit_cached(q, Fx)) && same11d(s3, natural_cubic_spline(q.knots, Fx))
            @test bits11d(Fx) == Fbits && Fx == Fcopy                                                 # F never clobbered by any of the fits above
        end
        for prec in (128, 256)
            setprecision(BigFloat, prec) do
                Fb = BigFloat.(F) .+ BigFloat(1) / 3; w = HIFitScratch{BigFloat}(n)
                s = _hi_fit_scratched(w, q, Fb); @test same11d(s, _hi_fit_cached(q, Fb)) && same11d(s, natural_cubic_spline(q.knots, Fb))
                Fb2 = Fb .* 2; s2 = _hi_fit_scratched(w, q, Fb2); @test same11d(s2, _hi_fit_cached(q, Fb2))
                @test same11d(_hi_fit_scratched(w, q, Fb), _hi_fit_cached(q, Fb)) && bits11d(Fb) == bits11d(BigFloat.(F) .+ BigFloat(1) / 3)
            end
        end
    end
    @testset "errors identical to _hi_fit_cached; fallbacks counted and exact; no clobber of F" begin
        q = nat.two3; n = length(q.knots); K = length(q.a); F = randn(rng, n); w = HIFitScratch{Float64}(n)
        @test errinfo11d(() -> _hi_fit_scratched(w, q, F[1:(end - 1)])) == errinfo11d(() -> _hi_fit_cached(q, F[1:(end - 1)])) !== nothing
        @test errinfo11d(() -> _hi_fit_scratched(w, q, vcat(F, 1.0))) == errinfo11d(() -> _hi_fit_cached(q, vcat(F, 1.0))) !== nothing
        for (lab, mk) in (("duplicate", k -> (k = copy(k); k[5] = k[4]; k)), ("NaN", k -> (k = copy(k); k[7] = NaN; k)), ("reversed", k -> reverse(k)),
                          ("shifted(stale dx)", k -> k .+ [i == 9 ? 1.0e-9 : 0.0 for i in eachindex(k)]), ("scaled(stale dx)", k -> k .* 1.0000001))
            qb = with_knots11d(q, mk(q.knots))
            @test errinfo11d(() -> _hi_fit_scratched(w, qb, F)) == errinfo11d(() -> _hi_fit_cached(qb, F)) !== nothing
        end
        qs = nat.di1; ns = length(qs.knots); Fs = randn(rng, ns)
        # wrong-size scratch / aliasing / element-type mismatch: served by the oracle, counted, never coerced
        ws = HIFitScratch{Float64}(n); @test same11d(_hi_fit_scratched(ws, qs, Fs), _hi_fit_cached(qs, Fs)); @test ws.nfallback[] == 1
        wa = HIFitScratch{Float64}(n); copyto!(wa.y, randn(rng, n)); Fa = wa.y; Fa0 = copy(Fa)
        ba = bits11d(Fa); sa = _hi_fit_scratched(wa, q, Fa); @test same11d(sa, _hi_fit_cached(q, Fa0)) && bits11d(Fa) == ba && wa.nfallback[] == 1        # F === scratch.y
        wv = HIFitScratch{Float64}(n); Fv = view(wv.c, 1:n); copyto!(Fv, randn(rng, n)); Fv0 = copy(Fv)
        bv = bits11d(Fv); sv = _hi_fit_scratched(wv, q, Fv); @test same11d(sv, _hi_fit_cached(q, Fv0)) && bits11d(Fv) == bv && wv.nfallback[] == 1       # F overlaps scratch.c
        w32d = HIFitScratch{Float32}(n); @test same11d(_hi_fit_scratched(w32d, q, F), _hi_fit_cached(q, F)) && w32d.nfallback[] == 1      # direct call, other ordinate type: counted oracle, no MethodError, no coercion
        w32 = HIFitScratch{Float32}(n); pa = _hi_spline_parts(q, F, ntuple(_ -> 0, K), w32); @test w32.nfallback[] == 1
        @test pa == _hi_spline_parts(q, F, ntuple(_ -> 0, K))                                                              # scratch of another eltype never coerces
        wd = HIFitScratch{Float64}(n); pv = _hi_spline_parts(q, view(copy(F), 1:n), ntuple(_ -> 0, K), wd); @test wd.nfallback[] == 1 && pv == _hi_spline_parts(q, F, ntuple(_ -> 0, K))
        @test_throws DimensionMismatch _hi_spline_parts(q, F[1:(end - 1)], ntuple(_ -> 0, K), HIFitScratch{Float64}(n))
        @test_throws ArgumentError _hi_spline_parts(q, F, ntuple(_ -> 9, K), HIFitScratch{Float64}(n))
        @test_throws ArgumentError HIFitScratch{Float64}(2)
    end
    @testset "_hi_spline_parts with scratch vs previous body: forced 1:8, adaptive, empty/reversed, Float64/Dual/BigFloat; poisoned scratch" begin
        for q in plans
            n = length(q.knots); K = length(q.a); x = q.knots
            F = (1.0 .+ 0.3 .* sin.(40.0 .* (x .- x[1]))) .* exp.(-3.0 .* (x .- x[1])); w = poison!(HIFitScratch{Float64}(n), NaN)
            for forced in (ntuple(_ -> 0, K), [ntuple(_ -> L, K) for L in 1:8]...)
                a = _hi_spline_parts(q, F, forced, w); b = reference_parts11d(q, F, forced)
                @test bits11d(collect(a[1])) == bits11d(collect(b[1])) && a[2:3] == b[2:3] && a == _hi_spline_parts(q, F, forced)
            end
            @test w.nfallback[] == 0
            v = randn(rng, n); Fd = [ForwardDiff.Dual(F[i], v[i]) for i in 1:n]; wd = HIFitScratch{eltype(Fd)}(n)
            pn = _hi_spline_parts(q, Fd, ntuple(_ -> 0, K), wd); pr = reference_parts11d(q, Fd, ntuple(_ -> 0, K))
            @test bits11d(collect(pn[1])) == bits11d(collect(pr[1])) && pn[2:3] == pr[2:3] && wd.nfallback[] == 0
            setprecision(BigFloat, 128) do
                Fb = BigFloat.(F); wb = HIFitScratch{BigFloat}(n); pn = _hi_spline_parts(q, Fb, ntuple(_ -> 5, K), wb); pr = reference_parts11d(q, Fb, ntuple(_ -> 5, K))
                @test bits11d(collect(pn[1])) == bits11d(collect(pr[1])) && pn[2:3] == pr[2:3]
            end
        end
        for q in (nat.di1, nat.two3, nat.raman)               # adaptive levels / ok flags on rough native-sized data
            K = length(q.a); F = randn(rng, length(q.knots)); w = HIFitScratch{Float64}(length(q.knots))
            @test _hi_spline_parts(q, F, ntuple(_ -> 0, K), w) == reference_parts11d(q, F, ntuple(_ -> 0, K))
        end
    end

    D11 = HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6); R11 = NODES5[:, 1:9]
    stagef(o) = (o.z, o.DI1_2s, o.DF_2g, o.DF_R, o.scale, o.unconverged, o.y)
    @testset "public stage path: fitwork vs plan path (values, scales, flags, record/replay, ForwardDiff, ownership, tasks)" begin
        for p in ([1.0, 1.0], [0.999, 1.001])
            r0 = PattersonLevels(); r1 = PattersonLevels()
            o0, _ = hi_diffusion_stage(RM5, D11, R11; hscale = p[1], nbscale = p[2], levels = r0, quadplan = nat)
            o1, _ = hi_diffusion_stage(RM5, D11, R11; hscale = p[1], nbscale = p[2], levels = r1, quadplan = nat, fitwork = true)
            @test stagef(o0) == stagef(o1) && r0.levels == r1.levels
            wk = HIQuadWork(nat, Float64)                                   # caller-owned, reused for record AND replay AND a second stage
            o2, _ = hi_diffusion_stage(RM5, D11, R11; hscale = p[1], nbscale = p[2], levels = replay(r0), quadplan = nat, fitwork = wk)
            o3, _ = hi_diffusion_stage(RM5, D11, R11; hscale = p[1], nbscale = p[2], levels = replay(r0), quadplan = nat)
            @test stagef(o2) == stagef(o3)
            o4, _ = hi_diffusion_stage(RM5, D11, R11; hscale = p[1], nbscale = p[2], quadplan = nat, fitwork = wk)
            @test stagef(o4) == stagef(o0)
            @test wk.di1.nfallback[] == 0 && wk.two3.nfallback[] == 0 && wk.raman.nfallback[] == 0 && wk.nmismatch[] == 0     # the scratch really served every fit
        end
        # a work object of a different eltype is counted and never coerced; invalid use is a clear error
        wm = HIQuadWork(nat, Float32); oM, _ = hi_diffusion_stage(RM5, D11, R11; quadplan = nat, fitwork = wm)
        @test stagef(oM) == stagef(hi_diffusion_stage(RM5, D11, R11; quadplan = nat)[1]) && wm.nmismatch[] > 0
        @test_throws ArgumentError hi_diffusion_stage(RM5, D11, R11; fitwork = true)
        @test_throws ArgumentError hi_diffusion_stage(RM5, D11, R11; quadplan = nat, fitwork = 1)
        # ForwardDiff through the whole stage (Dual state, stage-owned Dual scratch): bitwise the plan path
        function out(p; kw...)
            o, _ = hi_diffusion_stage(RM5, D11, R11; hscale = p[1], nbscale = p[2], quadplan = nat, kw...); m = (o.z .< 2000.0) .& (o.z .> 500.0)
            return vcat(o.DI1_2s[m], o.DF_2g[1][m], o.DF_2g[2][m], o.DF_R[1][m])
        end
        @test ForwardDiff.jacobian(q -> out(q), [0.999, 1.001]) == ForwardDiff.jacobian(q -> out(q; fitwork = true), [0.999, 1.001])
        # two independent caller-owned work objects on concurrent tasks sharing ONE plan: each equals the sequential result bitwise
        seq = [stagef(hi_diffusion_stage(RM5, D11, R11; hscale = h, quadplan = nat)[1]) for h in (1.0, 1.001)]
        par = fetch.([Threads.@spawn stagef(hi_diffusion_stage(RM5, D11, R11; hscale = h, quadplan = nat, fitwork = HIQuadWork(nat, Float64))[1]) for h in (1.0, 1.001)])
        @test par == seq
    end

    @testset "public hi_pde_integrals / hi_pde_corrections keyword: true, caller object, invalid, no quadplan" begin
        pops = HIPopulationSplines(R11; zs = ZS6, ze = ZE6)
        m = HIPDEModel(SETUP7, pops, hi_pde_coefficients(R11, pops, ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6), ACC5)
        cap = Ref{Any}(nothing)
        hi_pde_march(m; zs = 2500.0, ze = 2470.0, on_step = (k, zin, zout, y, S) -> (zout < 2490.0 && cap[] === nothing && (cap[] = (zout, copy(y), deepcopy(S.prev)))))
        z0, y0, st0 = cap[]
        r0 = PattersonLevels(); base = hi_pde_integrals(m, z0, y0, st0; quadplan = nat, levels = r0)
        wk = HIQuadWork(nat, Float64)
        for fw in (true, wk, nothing)
            r1 = PattersonLevels(); o = hi_pde_integrals(m, z0, y0, st0; quadplan = nat, levels = r1, fitwork = fw)
            @test o == base && r1.levels == r0.levels
        end
        @test wk.di1.nfallback[] == 0 && wk.two3.nfallback[] == 0 && wk.raman.nfallback[] == 0 && wk.nmismatch[] == 0
        @test hi_pde_integrals(m, z0, y0, st0; quadplan = nat, levels = replay(r0), fitwork = wk) == hi_pde_integrals(m, z0, y0, st0; quadplan = nat, levels = replay(r0))
        for bad in (false, 1, "x")
            @test_throws ArgumentError hi_pde_integrals(m, z0, y0, st0; quadplan = nat, fitwork = bad)
            @test_throws ArgumentError hi_pde_corrections(m; zs = 2500.0, ze = 2490.0, quadplan = nat, fitwork = bad)
        end
        @test_throws ArgumentError hi_pde_integrals(m, z0, y0, st0; fitwork = true)
        @test_throws ArgumentError hi_pde_integrals(m, z0, y0, st0; fitwork = wk)
        @test_throws ArgumentError hi_pde_corrections(m; zs = 2500.0, ze = 2490.0, fitwork = wk)
        # caller object passed to the stage is not copied: it is the one that served the fits (counters untouched, buffers reused)
        pre = (wk.two3.y, wk.raman.c); hi_pde_corrections(m; zs = 2500.0, ze = 2490.0, quadplan = nat, fitwork = wk); @test (wk.two3.y, wk.raman.c) === pre
    end
    @testset "reverse: CRC rule, one scratch used twice on one tape, prepared Mooncake vs allocating path and ForwardDiff" begin
        q = nat.two3; K = length(q.a); n = length(q.knots); forced = ntuple(_ -> 6, K); w = HIFitScratch{Float64}(n)
        F0 = randn(Random.Xoshiro(51), n); F1 = randn(Random.Xoshiro(52), n)
        @test ChainRulesCore.rrule(_hi_spline_parts, q, F0, forced, w) === nothing      # deliberately NO ChainRulesCore rule for the mutating scratch method
        # primitive status of the existing and the new rules (dispatch witnesses), plus the unchanged `_check_quadplan` zero-derivative rule
        wd0 = Base.get_world_counter(); Pq = typeof(q)
        @test !Mooncake.is_primitive(Mooncake.DefaultCtx, Mooncake.ReverseMode, Tuple{typeof(_hi_spline_parts), Pq, Vector{Float64}, NTuple{K,Int}, HIFitScratch{Float64}}, wd0)     # overlaid, not a primitive
        @test Mooncake.is_primitive(Mooncake.DefaultCtx, Mooncake.ReverseMode, Tuple{typeof(_hi_spline_parts), Pq, Vector{Float64}, NTuple{K,Int}}, wd0)
        @test Mooncake.is_primitive(Mooncake.DefaultCtx, Mooncake.ReverseMode, Tuple{typeof(CosmoRec._check_quadplan), typeof(nat), CosmoRec.HIPDESetup, CosmoRec.HIPDEAtom}, wd0)
        c = randn(Random.Xoshiro(54), K)
        parts(B, ws) = ws === nothing ? _hi_spline_parts(q, B, forced) : _hi_spline_parts(q, B, forced, ws)       # identical expressions; only the fit allocation differs
        function reuse(p, ws)                   # one scratch, two uses, plus a later use of the ordinates
            B = p .* F0; s1 = dot(c, collect(parts(B, ws)[1]))
            B2 = p .* F1; s2 = dot(c, collect(parts(B2, ws)[1]))
            return s1 + 2s2 + sum(B) + sum(B2)
        end
        fresh(p, ws) = reuse(p, nothing)                                                         # the allocating fit, same structure
        p0 = 1.0 .+ 0.1 .* randn(Random.Xoshiro(55), n)
        gr = gradient(reuse, prepare_gradient(reuse, MC11D, p0, Constant(w)), MC11D, p0, Constant(w))
        gf = gradient(fresh, prepare_gradient(fresh, MC11D, p0, Constant(w)), MC11D, p0, Constant(w))
        gfd = ForwardDiff.gradient(q2 -> fresh(q2, nothing), p0)
        @test gr == gf
        @test maximum(abs, gr .- gfd) <= 1.0e-8 * maximum(abs, gfd)
        # ---- reviewer's counterexample (preserved failing before the fix: Mooncake gave zeros): the scratch is built INSIDE the function and read after the call
        xr = [0.0, 0.1, 0.13, 0.4, 0.49, 1.0, 1.1]; qr = _hi_quad_plan(xr, [xr[1], xr[end]], 1.0e-60); Fr = collect(1.0:7.0)
        fr(p) = (ws = HIFitScratch{eltype(p)}(length(p)); _hi_spline_parts(qr, p, (6,), ws); sum(ws.y))
        gr1(p) = (ws = HIFitScratch{eltype(p)}(length(p)); _hi_spline_parts(qr, p, (6,), ws); sum(ws.y) + 2sum(ws.c) + 3sum(ws.b) + 5sum(ws.d))
        @test fr(Fr) == 28 && ForwardDiff.gradient(fr, Fr) == ones(7)
        @test gradient(fr, prepare_gradient(fr, MC11D, Fr), MC11D, Fr) == ones(7)
        for h in (gr1,)
            gm = gradient(h, prepare_gradient(h, MC11D, Fr), MC11D, Fr); gfd = ForwardDiff.gradient(h, Fr)
            @test maximum(abs, gm .- gfd) <= 1.0e-9 * maximum(abs, gfd) && all(!iszero, gm)
        end
        # all coefficient buffers, repeated calls on ONE scratch, same-tape accumulation, cross-aliased products, parts AND later reads together
        n2 = length(q.knots)
        function gmix(p)
            ws = HIFitScratch{eltype(p)}(n2); B0 = p .* F0
            s1 = dot(c, collect(_hi_spline_parts(q, B0, forced, ws)[1])); r1 = sum(ws.c) + dot(ws.b, ws.d) + sum(ws.z)
            B1 = p .* F1; s2 = dot(c, collect(_hi_spline_parts(q, B1, forced, ws)[1])); r2 = sum(ws.y .* ws.c) + sum(abs2, ws.b)
            return s1 + 2s2 + r1 + r2 + sum(ws.d) + sum(B0)
        end
        gm = gradient(gmix, prepare_gradient(gmix, MC11D, p0), MC11D, p0); gfd = ForwardDiff.gradient(gmix, p0)
        @test maximum(abs, gm .- gfd) <= 1.0e-8 * maximum(abs, gfd) && all(isfinite, gm)
        # preparation reuse with changed inputs, then back to the original inputs
        prm = prepare_gradient(gmix, MC11D, p0); g0 = gradient(gmix, prm, MC11D, p0)
        p2 = p0 .* (1 .+ 1.0e-3 .* randn(Random.Xoshiro(57), n2)); @test maximum(abs, gradient(gmix, prm, MC11D, p2) .- ForwardDiff.gradient(gmix, p2)) <= 1.0e-8 * maximum(abs, ForwardDiff.gradient(gmix, p2))
        @test gradient(gmix, prm, MC11D, p0) == g0
        # poisoning the scratch between preparation and use, and after: gradient and primal unchanged (no stale values)
        pg = prepare_gradient(reuse, MC11D, p0, Constant(w)); poison!(w, NaN)
        @test gradient(reuse, pg, MC11D, p0, Constant(w)) == gr
        p1 = p0 .* (1 .+ 1.0e-9 .* randn(Random.Xoshiro(56), n)); @test gradient(reuse, pg, MC11D, p1, Constant(w)) == gradient(fresh, prepare_gradient(fresh, MC11D, p1, Constant(w)), MC11D, p1, Constant(w))
        @test gradient(reuse, pg, MC11D, p0, Constant(w)) == gr                                    # back to the original inputs: no stale scratch dependence
        # nonzero initial accumulator through a later use of the same ordinates
        acc(p, ws) = reuse(p, ws) + sum(abs2, p .* F0) + dot(c, collect(parts(p .* F0, ws)[1]))
        acc_f(p, ws) = acc(p, nothing)
        @test gradient(acc, prepare_gradient(acc, MC11D, p0, Constant(w)), MC11D, p0, Constant(w)) == gradient(acc_f, prepare_gradient(acc_f, MC11D, p0, Constant(w)), MC11D, p0, Constant(w))
        # SubArray ordinates take the generic fallback (differentiable): gradient equals ForwardDiff
        wsub = HIFitScratch{Float64}(n); sv(p) = dot(c, collect(_hi_spline_parts(q, view(p .* F0, 1:n), forced, wsub)[1]))
        @test maximum(abs, gradient(sv, prepare_gradient(sv, MC11D, p0), MC11D, p0) .- ForwardDiff.gradient(sv, p0)) <= 1.0e-8 * maximum(abs, ForwardDiff.gradient(sv, p0))
        # two independent chains, each preparing AND evaluating with its OWN single scratch instance, on concurrent tasks (nthreads = $(Threads.nthreads()))
        xs = [p0 .* (1 .+ 1.0e-9 .* randn(Random.Xoshiro(60 + k), n)) for k in 1:2]
        chain(x) = (wk = HIFitScratch{Float64}(n); pr = prepare_gradient(reuse, MC11D, x, Constant(wk)); (gradient(reuse, pr, MC11D, x, Constant(wk)), gradient(reuse, pr, MC11D, x, Constant(wk))))
        seq = [chain(x) for x in xs]
        par = fetch.([Threads.@spawn chain(xs[k]) for k in 1:2])
        @test par == seq && all(r -> r[1] == r[2], seq)
        # dispatch witness: the Mooncake rule (not trace-through of the scratch) is what runs: its tables, not the scratch, define the gradient
        poisonq(q, s) = HIQuadPlan{length(q.a)}(q.knots, q.a, q.b, q.epsabs0, q.dx, q.g, q.off, q.alpha, q.t, s .* q.ell)
        objq(p, ws, qq) = dot(c, collect(_hi_spline_parts(qq, p .* F0, forced, ws)[1]))
        g1 = gradient(objq, prepare_gradient(objq, MC11D, p0, Constant(w), Constant(q)), MC11D, p0, Constant(w), Constant(q))
        g2 = gradient(objq, prepare_gradient(objq, MC11D, p0, Constant(w), Constant(poisonq(q, 2.0))), MC11D, p0, Constant(w), Constant(poisonq(q, 2.0)))
        @test g1 != g2 && maximum(abs, g2 .- 2 .* g1) <= 1.0e-12 * maximum(abs, g2)
        # stage-level prepared Mooncake: short production march with the stage-owned scratch vs the plan path
        function f8(x; fw = nothing)
            rows = hcat(R11[:, 1], reshape(x, size(R11, 1), 8)); pops = HIPopulationSplines(rows; zs = ZS6, ze = ZE6)
            o = hi_pde_corrections(HIPDEModel(SETUP7, pops, hi_pde_coefficients(rows, pops, ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6), ACC5);
                                   zs = 2500.0, ze = 2470.0, T = eltype(x), quadplan = nat, fitwork = fw)
            return vcat(o.DI1_2s, o.DF_2g[1], o.DF_2g[2], o.DF_R[1])
        end
        X0 = vec(R11[:, 2:9]); wv = randn(Random.Xoshiro(70), length(f8(X0)))
        @test f8(X0; fw = true) == f8(X0)
        objs(x, fw) = dot(wv, f8(x; fw = fw))
        gs = gradient(objs, prepare_gradient(objs, MC11D, X0, Constant(true)), MC11D, X0, Constant(true))
        gp = gradient(objs, prepare_gradient(objs, MC11D, X0, Constant(nothing)), MC11D, X0, Constant(nothing))
        @test gs == gp
    end
end
