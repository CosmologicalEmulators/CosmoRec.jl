# Chunk 12a: opt-in high-z prefix reuse (src/PrefixReuse.jl) on a small NON-TRIVIAL block-solver toy BEFORE the real history gates: three passes (source + two consumers) of a 12-state-style
# block march with a custom RK4 block solver (recycling ONE output buffer), all inputs active (background scale, baryon scale, rate coefficient, low-z feedback amplitude), a real
# `HIDiffusionFeedback` zero domain, pass-to-pass feedback coupling, per-pass changed settings, crossing block, fresh prepared evaluations. Compared against no-reuse: primal, ForwardDiff
# (incl. nested), prepared DI-Mooncake. Independent of the native data directory.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
using CosmoRec: PrefixReuse, PrefixReuseStats, PrefixCache, PrefixBinding, _prefix_solve!, hi_diffusion_feedback, HIDiffusionFeedback

const MC12 = AutoMooncake(; config = nothing)
const NST12 = 4
const BUF12 = zeros(NST12, 64)                         # recycled solver output storage: every solve overwrites and returns a view of it
const NCALL12 = Ref(0)
function toy_fb(amp, zmax)                               # real HIDiffusionFeedback with exact zero at/above zmax
    z = collect(range(2400.0, 300.0; length = 15))
    return hi_diffusion_feedback(z, amp .* (1 .+ 0.1 .* sin.(z ./ 200)), [amp .* cos.(z ./ 300), amp .* sin.(z ./ 250)], [amp .* z ./ 1e4]; zmin = 400.0, zmax = zmax)
end

@testset "Prefix query excursion conversion never understates the declaration" begin
    solve = (args...) -> nothing
    declaration(x) = PrefixReuse(solve; signatures = (:same,), max_query_excursion = x, pure_solver = true)
    setprecision(BigFloat, 256) do
        bound = BigFloat(1) + BigFloat(2)^(-54)
        @test BigFloat(declaration(bound).max_query_excursion) >= bound
        tiny = BigFloat(2)^(-1100)
        @test BigFloat(declaration(tiny).max_query_excursion) >= tiny
        @test_throws ArgumentError declaration(BigFloat(2)^2000)
    end
    @test declaration(0.0).max_query_excursion === 0.0
end
function toy_rhs!(du, u, p, z)
    r = p.hscale * (1.0 + 1.0e-3 * z) + p.nbscale * 2.0e-4 * z^1.5 / 1.0e3
    for k in 1:NST12; du[k] = -(0.3 + 0.1k) * r * u[k] * 1.0e-3 + 1.0e-4 * (u[mod1(k + 1, NST12)] - u[k]); end
    fb = p.rm.diffusion
    fb isa HIDiffusionFeedback && (du[1] += CosmoRec.hi_DI1_2s(fb, z) * u[1] * 1.0e-2 + CosmoRec.hi_DF_R(fb, 0, z) * u[2] * 1.0e-2)
    return du
end
function toy_solve(f!, u0, z0, znodes, p)                # RK4, 4 substeps per node interval, queries only within [znodes[end], z0] (excursion 0)
    NCALL12[] += 1
    n = length(znodes); T = promote_type(eltype(u0), typeof(p.hscale), typeof(p.nbscale), (p.rm.diffusion isa HIDiffusionFeedback ? CosmoRec.feedback_eltype(p.rm.diffusion) : Float64))
    y = collect(T, u0); BUF = view(BUF12, 1:NST12, 1:n)
    fill!(BUF, zero(eltype(BUF12)))
    BUFT = Matrix{T}(undef, NST12, n)                      # typed result; BUF12 is only a poisonable scratch witness
    BUFT[:, 1] .= y
    k1 = similar(y); k2 = similar(y); k3 = similar(y); k4 = similar(y)
    for m in 2:n
        h = (znodes[m] - znodes[m - 1]) / 4; z = znodes[m - 1]
        for _ in 1:4
            f!(k1, y, p, z); f!(k2, y .+ (h / 2) .* k1, p, z + h / 2); f!(k3, y .+ (h / 2) .* k2, p, z + h / 2); f!(k4, y .+ h .* k3, p, z + h)
            y = y .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4); z += h
        end
        BUFT[:, m] .= y
    end
    return BUFT
end
const EFF12 = (:eff,); const COS12 = (:cosmos,)
struct ToyRM{D}
    eff::typeof(EFF12); cosmos::typeof(COS12); dp_fallback::Nothing; diffusion::D
end
struct ToyP{R,T}
    rm::R; flag::Bool; hscale::T; nbscale::T
end
function toy_pass(rm, solve, c0, binding, hs, nbs; zn, block = 10)
    nz = length(zn); T = promote_type(typeof(hs), typeof(nbs), typeof(c0)); u0 = T[1.0, 0.5, 0.25, 0.125] .* c0
    ys = Vector{Vector{T}}(undef, nz); ys[1] = u0; i = 1
    p = ToyP(rm, true, hs, nbs)
    while i < nz
        j = min(i + block, nz)
        sol = binding === nothing ? solve(toy_rhs!, ys[i], zn[i], zn[i:j], p) : _prefix_solve!(binding, solve, toy_rhs!, ys[i], zn, i, j, p)
        for m in (i + 1):j; ys[m] = collect(sol[:, m - i + 1]); end
        i = j
    end
    return ys
end
# 3 passes: pass 0 no feedback; pass k feedback amplitude θ[4] * mean(previous pass final state)
function toy_history(θ; reuse = false, stats = PrefixReuseStats(), excursion = 0.0, sigs = (:s, :s, :s), eff_c = EFF12, zmax = 1500.0, solve = toy_solve, perpass_scale = 0.0, fb_override = nothing)
    zn = collect(range(3000.0, 100.0; length = 61)); out = Any[]
    plan = reuse ? PrefixReuse(solve; signatures = collect(sigs), max_query_excursion = excursion, pure_solver = true) : nothing
    reuse && CosmoRec._reset!(stats); reuse && (stats.status = :active)
    cache = PrefixCache(); rm0 = ToyRM(EFF12, COS12, nothing, nothing); fbprev = nothing; ys = nothing
    for it in 0:2
        hs = θ[1] * (1 + perpass_scale * it); nbs = θ[2]
        fb = it == 0 ? nothing : (fb_override === nothing ? toy_fb(θ[4] * sum(ys[end]) / NST12 + 1e-3, zmax) : fb_override)
        rm = ToyRM(it == 0 ? EFF12 : eff_c, COS12, nothing, fb)
        b = plan === nothing ? nothing : CosmoRec._prefix_pass_binding(plan, stats, cache, it, rm0, rm)
        ys = toy_pass(rm, solve, θ[3], b, hs, nbs; zn = zn)
        push!(out, ys)
    end
    return out
end
flat12(o) = reduce(vcat, [reduce(vcat, y) for y in o])
bits12(x) = x isa Float64 ? reinterpret(UInt64, x) : x isa ForwardDiff.Dual ? (bits12(ForwardDiff.value(x)), map(bits12, ForwardDiff.partials(x).values)) : x
sc12(θ; kw...) = (o = toy_history(θ; kw...); sum(abs2, flat12(o)) + 3 * sum(o[3][end]) + 2 * sum(o[2][31]))      # consumers all read the cached prefix

@testset "Chunk 12a: prefix reuse toy (block solver, three passes, feedback zero domain)" begin
    θ0 = [1.0, 1.1, 0.9, 0.7]
    @testset "primal: OFF == ON bitwise, truthful stats, owned copies" begin
        NCALL12[] = 0; off = toy_history(θ0); noff = NCALL12[]
        st = PrefixReuseStats(); NCALL12[] = 0; on = toy_history(θ0; reuse = true, stats = st); non = NCALL12[]
        @test all(bits12.(flat12(off)) .== bits12.(flat12(on)))
        @test st.status === :active && st.source_blocks == 6 && st.consumer_blocks == 12 && st.reused_blocks > 0 && st.reused_blocks + st.recomputed_blocks == 12
        @test noff == 18 && non == noff - st.reused_blocks              # skipped solve calls are not work
        @test st.reused_blocks == 6 && st.recomputed_blocks == 6 && non == 12        # per consumer the blocks ending at nodes 11, 21, 31 (z >= 1550 > zmax = 1500) are reused; crossing block 41.. recomputed
        # excursion too large: the block that would reach below zmax is refused (counted), results unchanged
        st2 = PrefixReuseStats(); on2 = toy_history(θ0; reuse = true, stats = st2, excursion = 400.0)
        @test all(bits12.(flat12(off)) .== bits12.(flat12(on2))) && st2.reused_blocks < st.reused_blocks && st2.refused_excursion > 0
        # poisoning the recycled solver buffer between calls cannot corrupt the cached prefix
        fill!(BUF12, NaN); on3 = toy_history(θ0; reuse = true); @test all(bits12.(flat12(off)) .== bits12.(flat12(on3)))
    end
    @testset "refusals: signature, unknown feedback, other model object, per-pass settings, solver binding, iterations" begin
        off = toy_history(θ0)
        st = PrefixReuseStats(); on = toy_history(θ0; reuse = true, stats = st, sigs = (:s, :s, :t)); @test all(bits12.(flat12(off)) .== bits12.(flat12(on))) && st.refused_pass == 1 && st.reused_blocks == 3
        st = PrefixReuseStats(); on = toy_history(θ0; reuse = true, stats = st, sigs = (:s, :t, :t)); @test all(bits12.(flat12(off)) .== bits12.(flat12(on))) && st.refused_pass == 2 && st.reused_blocks == 0
        st = PrefixReuseStats(); on = toy_history(θ0; reuse = true, stats = st, eff_c = (:other,)); @test all(bits12.(flat12(off)) .== bits12.(flat12(on))) && st.reused_blocks == 0 && st.refused_key > 0
        st = PrefixReuseStats(); onp = toy_history(θ0; reuse = true, stats = st, perpass_scale = 0.01); offp = toy_history(θ0; perpass_scale = 0.01)      # changed hscale per pass: u0/scales differ -> recompute
        @test all(bits12.(flat12(offp)) .== bits12.(flat12(onp))) && st.reused_blocks == 0 && st.refused_key > 0
        @test_throws ArgumentError PrefixReuse(toy_solve; signatures = [:s], max_query_excursion = -1.0, pure_solver = true)
        @test_throws ArgumentError PrefixReuse(toy_solve; signatures = [:s], max_query_excursion = Inf, pure_solver = true)
        @test_throws ArgumentError PrefixReuse(toy_solve; signatures = [:s], max_query_excursion = 0.0, pure_solver = false)
        @test_throws ArgumentError PrefixReuse(toy_solve; signatures = [], max_query_excursion = 0.0, pure_solver = true)
        # unknown feedback object -> pass refused, never an unsafe skip
        st = PrefixReuseStats(); onu = toy_history(θ0; reuse = true, stats = st, fb_override = (:custom,))
        @test st.reused_blocks == 0 && st.refused_pass > 0
    end
    @testset "ForwardDiff and nested: OFF == ON bitwise (values and partials)" begin
        jo = ForwardDiff.jacobian(θ -> flat12(toy_history(θ)), θ0); jn = ForwardDiff.jacobian(θ -> flat12(toy_history(θ; reuse = true)), θ0)
        @test all(bits12.(jo) .== bits12.(jn)) && all(isfinite, jn) && any(!iszero, jn[:, 4])           # feedback amplitude active
        go(t) = ForwardDiff.gradient(θ -> sc12(θ), t); gn(t) = ForwardDiff.gradient(θ -> sc12(θ; reuse = true), t)
        ho = ForwardDiff.jacobian(go, θ0); hn = ForwardDiff.jacobian(gn, θ0)                           # nested Dual (second derivatives)
        @test all(bits12.(ho) .== bits12.(hn))
    end
    @testset "prepared DI-Mooncake: ON vs OFF, p0 -> p1 -> p0, fresh stats, 3 consumers accumulate into the source" begin
        offf(θ) = sc12(θ); onf(θ) = sc12(θ; reuse = true)
        go = gradient(offf, prepare_gradient(offf, MC12, θ0), MC12, θ0)
        prep = prepare_gradient(onf, MC12, θ0)
        g0 = gradient(onf, prep, MC12, θ0)
        gfd = ForwardDiff.gradient(offf, θ0)
        @test all(isfinite, g0) && maximum(abs, g0 .- go) <= 1.0e-12 * maximum(abs, go) && maximum(abs, g0 .- gfd) <= 1.0e-9 * maximum(abs, gfd)
        @test all(!iszero, g0)                                                                         # every parameter (background, baryon, coefficient, feedback amplitude) keeps its gradient
        θ1 = θ0 .* [1.01, 0.98, 1.02, 1.05]
        g1 = gradient(onf, prep, MC12, θ1); g1o = gradient(offf, prepare_gradient(offf, MC12, θ1), MC12, θ1)
        @test maximum(abs, g1 .- g1o) <= 1.0e-12 * maximum(abs, g1o)
        @test gradient(onf, prep, MC12, θ0) == g0                                                       # back to p0: no stale cache (a fresh cache is built inside every call)
        # three consumers each read the prefix; the summed objective equals the sum of per-consumer objectives' gradients
        c1(θ) = sum(toy_history(θ; reuse = true)[2][end]); c2(θ) = sum(toy_history(θ; reuse = true)[3][end]); c3(θ) = sum(toy_history(θ; reuse = true)[3][31])
        csum(θ) = c1(θ) + c2(θ) + c3(θ)
        gs = gradient(csum, prepare_gradient(csum, MC12, θ0), MC12, θ0)
        gp = gradient(c1, prepare_gradient(c1, MC12, θ0), MC12, θ0) .+ gradient(c2, prepare_gradient(c2, MC12, θ0), MC12, θ0) .+ gradient(c3, prepare_gradient(c3, MC12, θ0), MC12, θ0)
        @test maximum(abs, gs .- gp) <= 1.0e-12 * maximum(abs, gp)
        # stats are passive Ints: evaluating with a caller stats object under Mooncake does not freeze the active parent
        stm = PrefixReuseStats(); fst(θ) = sc12(θ; reuse = true, stats = stm); gst = gradient(fst, prepare_gradient(fst, MC12, θ0), MC12, θ0)
        @test maximum(abs, gst .- go) <= 1.0e-12 * maximum(abs, go)
    end
end
