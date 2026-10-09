# Chunk 11c: the quadrature plan's CACHED spline fit (`CosmoRec._hi_fit_cached`, used by `_hi_spline_parts` on the opt-in plan path) vs the unchanged generic
# `natural_cubic_spline`, which stays the independent oracle. The plan's stored Float64 factors (dx, g, off, alpha, t) come from the same expressions as the
# constructor; only the ACTIVE ordinates are fitted per call, into fresh arrays. Gates: spline fields (x, y, b, c, d) BITWISE equal for Float64 (also SubArray input),
# ForwardDiff Dual, nested Dual, BigFloat 128/256 (the stored factors stay Float64, exactly as in natural_cubic_spline), on the native 782/1728/Raman slices and on minimal-n=3 /
# highly non-uniform synthetic geometry; `_hi_spline_parts` vs a verbatim copy of the previous body (forced levels 1:8, adaptive, empty and reversed intervals);
# ownership (fresh y/b/c/d, plan unchanged, task interleaving). Reverse rules, the stale-plan guard and the full-stage gates are the unchanged chunk11 tests.
# REQUIRES COSMOREC_NATIVE_DATA_DIR (see test/chunk7_helpers.jl); fails loudly if unset.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using CosmoRec
using CosmoRec: HIQuadPlan, _hi_quad_plan, _hi_quad_plans, _hi_spline_parts, _hi_fit_cached, NATIVE_HI_PDE_ATOM
@isdefined(SETUP7) || include("chunk7_helpers.jl")

bits11c(x::Float64) = reinterpret(UInt64, x)
bits11c(x::Float32) = reinterpret(UInt32, x)
bits11c(x::Float16) = reinterpret(UInt16, x)
bits11c(x::ForwardDiff.Dual) = (bits11c(ForwardDiff.value(x)), map(bits11c, ForwardDiff.partials(x).values))
bits11c(x::BigFloat) = string(x)
bits11c(x::AbstractArray) = map(bits11c, x)
same_spline11c(a, b) = typeof(a) == typeof(b) && bits11c(a.x) == bits11c(b.x) && bits11c(a.y) == bits11c(b.y) && bits11c(a.b) == bits11c(b.b) && bits11c(a.c) == bits11c(b.c) && bits11c(a.d) == bits11c(b.d)
# verbatim previous body of `_hi_spline_parts` (natural_cubic_spline refit per call): the reference
function reference_parts11c(q::HIQuadPlan{K}, F::AbstractVector, forced::NTuple{K,Int}) where {K}
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
function synth_knots11c(rng, n)
    h = exp.(randn(rng, n - 1) .* 4.0)
    return cumsum(vcat(0.3, h ./ sum(h) .* 2.0))
end
function synth_plans11c()
    rng = Random.Xoshiro(1103); out = HIQuadPlan[]
    for n in (3, 7, 40, 300)
        x = synth_knots11c(rng, n); lo = x[1]; hi = x[end]; w = hi - lo
        push!(out, _hi_quad_plan(x, [lo, lo + 0.13w, lo + 0.13w, lo + 0.5w, lo + 0.31w, hi], 1.0e-60))      # includes an empty and a reversed sub-interval
    end
    return out
end

@testset "Chunk 11c: cached spline-fit geometry of the quadrature plan" begin
    rng = Random.Xoshiro(20261108)
    nat = _hi_quad_plans(SETUP7)
    plans = vcat(synth_plans11c(), [nat.di1, nat.two3, nat.raman])
    @testset "fit fields bitwise vs natural_cubic_spline, n = $(length(q.knots))" for q in plans
        n = length(q.knots); x = q.knots
        Fs = (randn(rng, n), (1.0 .+ 0.3 .* sin.(40.0 .* (x .- x[1]) ./ (x[end] - x[1]))) .* exp.(-3.0 .* (x .- x[1])) .* (1.0 .+ 0.05 .* randn(rng, n)), fill(2.5, n), 1.0e-30 .* randn(rng, n))
        for F in Fs
            @test same_spline11c(_hi_fit_cached(q, F), natural_cubic_spline(q.knots, F))
            Fc = copy(F); @test same_spline11c(_hi_fit_cached(q, view(Fc, 1:n)), natural_cubic_spline(q.knots, view(Fc, 1:n)))      # SubArray input
        end
        F = Fs[2]; v = randn(rng, n)
        # ForwardDiff Dual state, nested Dual, same stored Float64 factors
        D(t) = ForwardDiff.Dual(t, 1.0)
        Fd = [ForwardDiff.Dual(F[i], v[i]) for i in 1:n]
        @test same_spline11c(_hi_fit_cached(q, Fd), natural_cubic_spline(q.knots, Fd))
        Tg = ForwardDiff.Tag(:o, Float64); Ti = ForwardDiff.Tag(:i, ForwardDiff.Dual{typeof(Tg),Float64,1})
        Fn = [ForwardDiff.Dual{typeof(Ti)}(ForwardDiff.Dual{typeof(Tg)}(F[i], v[i]), ForwardDiff.Dual{typeof(Tg)}(v[i], 0.5)) for i in 1:n]
        @test same_spline11c(_hi_fit_cached(q, Fn), natural_cubic_spline(q.knots, Fn))
        for prec in (128, 256)
            setprecision(BigFloat, prec) do
                Fb = BigFloat.(F) .+ BigFloat(1) / 3
                @test same_spline11c(_hi_fit_cached(q, Fb), natural_cubic_spline(q.knots, Fb))
            end
        end
    end
    @testset "_hi_spline_parts vs the previous body: forced 1:8, adaptive, empty/reversed intervals, Float64 / Dual / BigFloat" begin
        for q in plans
            n = length(q.knots); K = length(q.a); x = q.knots
            F = (1.0 .+ 0.3 .* sin.(40.0 .* (x .- x[1]))) .* exp.(-3.0 .* (x .- x[1]))
            for forced in (ntuple(_ -> 0, K), [ntuple(_ -> L, K) for L in 1:8]...)
                @test bits11c(collect(_hi_spline_parts(q, F, forced)[1])) == bits11c(collect(reference_parts11c(q, F, forced)[1])) &&
                      _hi_spline_parts(q, F, forced)[2:3] == reference_parts11c(q, F, forced)[2:3]
            end
            v = randn(rng, n)
            Fd = [ForwardDiff.Dual(F[i], v[i]) for i in 1:n]
            pn = _hi_spline_parts(q, Fd, ntuple(_ -> 0, K)); pr = reference_parts11c(q, Fd, ntuple(_ -> 0, K))
            @test bits11c(collect(pn[1])) == bits11c(collect(pr[1])) && pn[2:3] == pr[2:3]
            setprecision(BigFloat, 128) do
                Fb = BigFloat.(F); pn = _hi_spline_parts(q, Fb, ntuple(_ -> 5, K)); pr = reference_parts11c(q, Fb, ntuple(_ -> 5, K))
                @test bits11c(collect(pn[1])) == bits11c(collect(pr[1])) && pn[2:3] == pr[2:3]
            end
        end
        # corner endpoints of the production grid and the native slices (adaptive levels, ok flags)
        for q in (nat.di1, nat.two3, nat.raman)
            K = length(q.a); F = randn(rng, length(q.knots))
            @test _hi_spline_parts(q, F, ntuple(_ -> 0, K)) == reference_parts11c(q, F, ntuple(_ -> 0, K))
        end
    end
    @testset "contract: length / geometry mismatches rejected; ownership; task interleaving" begin
        q = nat.two3; n = length(q.knots); K = length(q.a); F = randn(rng, n)
        @test_throws DimensionMismatch _hi_spline_parts(q, F[1:(end - 1)], ntuple(_ -> 0, K))
        @test_throws DimensionMismatch _hi_fit_cached(q, F[1:(end - 1)])
        @test_throws DimensionMismatch _hi_fit_cached(q, vcat(F, 1.0))
        @test_throws ArgumentError _hi_spline_parts(q, F, ntuple(i -> i == 1 ? 9 : 0, K))
        # a plan whose stored geometry is inconsistent (private misuse) is rejected, not read out of bounds
        bad = HIQuadPlan{K}(q.knots, q.a, q.b, q.epsabs0, q.dx[1:(end - 1)], q.g, q.off, q.alpha, q.t, q.ell)
        @test_throws DimensionMismatch _hi_fit_cached(bad, F)
        badx = HIQuadPlan{K}(q.knots, q.a, q.b, q.epsabs0, vcat(-abs.(q.dx[1:1]), q.dx[2:end]), q.g, q.off, q.alpha, q.t, q.ell)
        @test_throws ArgumentError _hi_fit_cached(badx, F)
        # ownership: the fit holds fresh y/b/c/d (never F, never plan arrays); mutating F afterwards or fitting again does not change an earlier fit; the plan is unchanged
        snap = (copy(q.knots), copy(q.dx), copy(q.g), copy(q.off), copy(q.alpha), copy(q.t), copy(q.ell))
        F1 = copy(F); s1 = _hi_fit_cached(q, F1); s1copy = (copy(s1.y), copy(s1.b), copy(s1.c), copy(s1.d))
        F1 .= 0.0; F2 = randn(rng, n); s2 = _hi_fit_cached(q, F2)
        @test s1.y !== F1 && s1.y !== s2.y && s1.b !== s2.b && s1.c !== s2.c && s1.d !== s2.d
        @test all(a !== b for a in (s1.y, s1.b, s1.c, s1.d) for b in (q.dx, q.g, q.off, q.alpha, q.t))
        @test (s1.y, s1.b, s1.c, s1.d) == s1copy
        @test (q.knots, q.dx, q.g, q.off, q.alpha, q.t, q.ell) == snap
        # interleaved tasks sharing ONE plan: each fit equals its own sequential result
        Fa = randn(rng, n); Fb2 = randn(rng, n); ra = natural_cubic_spline(q.knots, Fa); rb = natural_cubic_spline(q.knots, Fb2)
        ta = Threads.@spawn [same_spline11c(_hi_fit_cached(q, Fa), ra) for _ in 1:20]; tb = Threads.@spawn [same_spline11c(_hi_fit_cached(q, Fb2), rb) for _ in 1:20]
        @test all(fetch(ta)) && all(fetch(tb))
        # allocations of the cached fit are below the generic constructor
        _hi_fit_cached(q, F); natural_cubic_spline(q.knots, F)
        a_c = @allocated _hi_fit_cached(q, F); a_g = @allocated natural_cubic_spline(q.knots, F)
        println("Chunk11c fit allocated bytes (n = $n): cached $a_c, generic $a_g")
        @test a_c < a_g
    end
end

@testset "Chunk 11c repair regressions: narrow element types and malformed / stale geometry" begin
    rng = Random.Xoshiro(20261109)
    nat = _hi_quad_plans(SETUP7)
    xr = [0.0, 0.1, 0.13, 0.4, 0.49, 1.0, 1.1]
    qr = _hi_quad_plan(xr, [xr[1], xr[end]], 1.0e-60)
    plans = vcat(synth_plans11c(), [qr, nat.di1, nat.two3, nat.raman])
    @testset "Float32 / Float16 / Dual{Float32}: coefficients BITWISE as natural_cubic_spline (gamma stored in T before the solve)" begin
        for q in plans
            n = length(q.knots)
            for T in (Float32, Float16)
                F = randn(Random.MersenneTwister(1), T, n)
                @test same_spline11c(_hi_fit_cached(q, F), natural_cubic_spline(q.knots, F))
            end
            F32 = randn(rng, Float32, n); v32 = randn(rng, Float32, n)
            Fd = [ForwardDiff.Dual(F32[i], v32[i]) for i in 1:n]
            @test same_spline11c(_hi_fit_cached(q, Fd), natural_cubic_spline(q.knots, Fd))
        end
    end
    @testset "_hi_spline_parts Float32 / Float16: parts bits, levels and ok flags equal the previous body (forced 1:8 and adaptive)" begin
        for q in plans
            n = length(q.knots); K = length(q.a); x = q.knots
            for T in (Float32, Float16)
                F = T.((1.0 .+ 0.3 .* sin.(40.0 .* (x .- x[1]))) .* exp.(-3.0 .* (x .- x[1])))
                for forced in (ntuple(_ -> 0, K), [ntuple(_ -> L, K) for L in 1:8]...)
                    a = _hi_spline_parts(q, F, forced); b = reference_parts11c(q, F, forced)
                    @test bits11c(collect(a[1])) == bits11c(collect(b[1])) && a[2:3] == b[2:3]
                end
            end
        end
    end
    @testset "malformed or stale geometry is rejected like the generic fit (strictly increasing knots, stored dx == actual spacing)" begin
        q = nat.two3; n = length(q.knots); K = length(q.a); F = randn(rng, n)
        mutated(f) = (k = copy(q.knots); f(k); HIQuadPlan{K}(k, q.a, q.b, q.epsabs0, q.dx, q.g, q.off, q.alpha, q.t, q.ell))
        for (name, f) in (("duplicate knot", k -> (k[3] = k[2])), ("reversed knots", k -> reverse!(k)), ("NaN knot", k -> (k[5] = NaN)), ("Inf knot", k -> (k[5] = Inf)),
                          ("unsorted pair", k -> ((k[7], k[8]) = (k[8], k[7]))),
                          ("sorted but shifted knot (stale factors)", k -> (k[4] += 0.25 * (k[5] - k[4]))))
            qm = mutated(f)
            @test_throws ArgumentError _hi_fit_cached(qm, F)
            @test_throws ArgumentError _hi_spline_parts(qm, F, ntuple(_ -> 0, K))
            name in ("sorted but shifted knot (stale factors)",) || @test_throws ArgumentError natural_cubic_spline(qm.knots, F)      # the generic constructor's own domain error
        end
        bad_dx = copy(q.dx); bad_dx[10] *= 1.0000001
        @test_throws ArgumentError _hi_fit_cached(HIQuadPlan{K}(q.knots, q.a, q.b, q.epsabs0, bad_dx, q.g, q.off, q.alpha, q.t, q.ell), F)       # corrupt stored dx
        for (nm, qb) in (("dx", HIQuadPlan{K}(q.knots, q.a, q.b, q.epsabs0, q.dx[1:(end - 1)], q.g, q.off, q.alpha, q.t, q.ell)),
                         ("g", HIQuadPlan{K}(q.knots, q.a, q.b, q.epsabs0, q.dx, q.g[1:(end - 1)], q.off, q.alpha, q.t, q.ell)),
                         ("off", HIQuadPlan{K}(q.knots, q.a, q.b, q.epsabs0, q.dx, q.g, q.off[1:(end - 1)], q.alpha, q.t, q.ell)),
                         ("alpha", HIQuadPlan{K}(q.knots, q.a, q.b, q.epsabs0, q.dx, q.g, q.off, vcat(q.alpha, 1.0), q.t, q.ell)),
                         ("t", HIQuadPlan{K}(q.knots, q.a, q.b, q.epsabs0, q.dx, q.g, q.off, q.alpha, q.t[1:(end - 1)], q.ell)))
            @test_throws DimensionMismatch _hi_fit_cached(qb, F)
        end
        # the unmodified plan still fits; the public full-grid stale-plan guard is unchanged
        @test same_spline11c(_hi_fit_cached(q, F), natural_cubic_spline(q.knots, F))
        fresh = _hi_quad_plans(SETUP7); CosmoRec._check_quadplan(fresh, SETUP7, NATIVE_HI_PDE_ATOM)
        fresh.xfull[7] += 1.0e-9
        @test_throws ArgumentError CosmoRec._check_quadplan(fresh, SETUP7, NATIVE_HI_PDE_ATOM)
    end
end
