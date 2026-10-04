# Chunk 5f: prepared DI-Mooncake THROUGH the caller callback's per-solve prepared state Jacobian (`PreparedJac5`, test/chunk5_helpers.jl), vs the
# unprepared oracle `jac5!` and vs ForwardDiff. Scalar projection s = sum(W .* J). A prebuilt Float64 cache is reused across calls, so Mooncake
# differentiates through its live mutation of z and its reused config/buffer. Minimal: 7-state z = 400 and 12-state z = 2500, three preparations each.
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

const MC5F = AutoMooncake(; config = nothing)
const ROWS5F = read_fixture5a(joinpath(@__DIR__, "fixtures", "native_ode_rhs.txt"))
const OBS5F = Dict{String,Float64}()
rec5f(k, v) = (OBS5F[k] = max(get(OBS5F, k, 0.0), v))

# J at (u, p, z) through the oracle or through a cache; the state is promoted to the solve element type, as a solve does
function jac5f(cache, u, p, z)
    T = promote_type(eltype(u), typeof(z), typeof(p.hscale), typeof(p.nbscale)); uu = T.(u); J = zeros(T, length(u), length(u))
    cache === nothing ? jac5!(J, uu, p, z) : cache(J, uu, p, z)
    return J
end

@testset "Chunk 5f: prepared Mooncake through the prepared state-Jacobian cache" begin
    for (flag, z0) in ((false, 400.0), (true, 2500.0))
        r = first(r for r in ROWS5F if r.flag == flag && r.z == z0 && r.pat == 0); n = length(r.y); p0 = RecombinationODEParams(RM5, flag)
        W = randn(MersenneTwister(5), n, n); label = "$(n)-state z=$(Int(z0))"
        @testset "$label" begin
            su(cache) = u -> dot(W, jac5f(cache, u, p0, z0))
            sz(cache) = x -> dot(W, jac5f(cache, r.y, p0, x[1]))
            refu = ForwardDiff.gradient(su(nothing), r.y); refz = ForwardDiff.gradient(sz(nothing), [z0])
            pre = PreparedJac5(r.y, p0, z0)                     # one Float64 cache, as a primal solve builds it, reused by all gradients below
            gor = nothing
            for (tag, f, x, ref) in (("oracle d/du", su(nothing), r.y, refu), ("cache d/du", su(pre), r.y, refu), ("cache d/dz", sz(pre), [z0], refz))
                prep = prepare_gradient(f, MC5F, x)
                g1 = gradient(f, prep, MC5F, x); g2 = gradient(f, prep, MC5F, x)
                e = maximum(abs.(g1 .- ref)) / max(maximum(abs, ref), floatmin()); rec5f("$label $tag", e)
                @test e < 1e-12                                # vs ForwardDiff through the oracle (5a AD tolerance)
                @test g1 == g2                                 # repeated prepared gradient, cache reused and mutated in between
                tag == "oracle d/du" && (gor = g1)
                tag == "cache d/du" && @test maximum(abs.(g1 .- gor)) / max(maximum(abs, gor), floatmin()) < 1e-12
            end
            @test pre.nfallback[] == 0                         # every call took the prepared path
        end
    end
    @testset "info" begin
        println("Chunk5f Mooncake ", pkgversion(Mooncake), " ForwardDiff ", pkgversion(ForwardDiff))
        foreach(kv -> println("Chunk5f ", kv[1], " = ", kv[2]), sort!(collect(OBS5F)))
    end
end
