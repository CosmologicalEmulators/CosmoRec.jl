# Chunk 4d AD tests: derivatives of the FROZEN-BRANCH reset/pack operator with respect to the state and `fHe` (ForwardDiff vs 256-bit central differences, prepared Mooncake VJPs, two
# independent preparations, changed inputs). The branch is decided on primal values; no derivative of the decision/event location is implemented or claimed.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(FX4D) || include("chunk4d_helpers.jl")

const MC4D = AutoMooncake(; config = nothing)
const OBS4DAD = Dict{String,Float64}()
rec4dad(k, v) = (OBS4DAD[k] = max(get(OBS4DAD, k, 0.0), v))

# x = (X[1:15]; fHe); branch frozen from the primal x0 (decided once)
function make_f(zs, x0)
    sw = helium_switch_condition(zs, x0[16], x0[8])
    return x -> begin
        X = sw ? helium_switch_reset(x[1:15], 6, 7, x[16]) : x[1:15]
        pack_ysol(X, 6, 7; helium = !sw)
    end, sw
end

@testset "Chunk 4d: helium switch AD (frozen branch)" begin
    for (label, r) in (("switched", first(r for r in FX4D if r.cond && r.zs == 2000.0)), ("not switched", first(r for r in FX4D if !r.cond && r.zs == 2000.0)),
                       ("switched by zs < 200", first(r for r in FX4D if r.cond && r.zs == 150.0)))
        x0 = vcat(r.Xpre, r.fHe)
        f, sw = make_f(r.zs, x0)
        @testset "$label" begin
            @test sw == r.cond
            J = ForwardDiff.jacobian(f, x0)
            Jb = setprecision(BigFloat, 256) do
                cols = map(1:16) do i
                    h = BigFloat(1.0e-20) * max(abs(BigFloat(x0[i])), BigFloat(1.0e-30)); tp = BigFloat.(x0); tm = BigFloat.(x0); tp[i] += h; tm[i] -= h
                    Float64.((f(tp) .- f(tm)) ./ (2h))
                end
                reduce(hcat, cols)
            end
            e = maximum(abs.(J .- Jb)); rec4dad("FD:" * label, e); @test e < 1e-12
            @test size(J) == (length(r.ypost), 16)
            # structure: retained entries are the identity, helium entries are constants (switched) or the identity (not switched)
            @test J[1, 15] == 1.0 && J[2, 2] == 1.0
            if sw
                @test all(J[:, 16] .== 0)                                     # y (no helium entries) does not depend on fHe after the switch
                @test all(J[:, 8:14] .== 0)
            else
                @test all(J[:, 16] .== 0) && J[8, 8] == 1.0                    # no switch: He 1s passes through
            end
            obj(x, w) = dot(w, f(x))
            rng = Random.Xoshiro(20261007)
            nout = length(r.ypost)
            prepA = prepare_gradient(obj, MC4D, x0, Constant(randn(rng, nout)))
            prepB = prepare_gradient(obj, MC4D, x0, Constant(randn(rng, nout)))
            for kk in 1:2
                x = kk == 1 ? x0 : x0 .* (1 .+ 1.0e-9 .* randn(rng, 16))      # tiny change keeps the frozen branch
                w = randn(rng, nout)
                ref = ForwardDiff.jacobian(f, x)' * w
                gA = gradient(obj, prepA, MC4D, x, Constant(w)); gB = gradient(obj, prepB, MC4D, x, Constant(w))
                e = max(maximum(abs.(gA .- ref)), maximum(abs.(gB .- ref))) / max(maximum(abs, ref), floatmin()); rec4dad("VJP:" * label, e)
                @test e < 1e-12
                @test gA == gB
            end
        end
    end
    @testset "fHe dependence of the reset state (before packing) and no event sensitivity" begin
        r = first(r for r in FX4D if r.cond && r.zs == 2000.0)
        g(fHe) = helium_switch_reset(r.Xpre, 6, 7, fHe)
        d = ForwardDiff.derivative(g, r.fHe)
        @test d[8] == 1.0 && all(d[[1:7; 9:15]] .== 0)
        # the decision is a primal Bool: its "derivative" with respect to the He 1s population is identically 0 on both sides (no event-location sensitivity is carried)
        dec(Xi0) = Float64(helium_switch_condition(2000.0, r.fHe, Xi0))
        @test ForwardDiff.derivative(dec, r.fHe - 5.0e-8) == 0.0 && ForwardDiff.derivative(dec, r.fHe - 2.0e-7) == 0.0
        # the state as a function of X[He 1s] is discontinuous across the criterion: derivative 0 (switched side) vs 1 (unswitched side), not a derivative of a smooth map
        function h(Xi0)
            X = Vector{typeof(Xi0)}(r.Xpre); X[8] = Xi0
            Xn, _ = helium_switch(X, 2000.0, 6, 7, r.fHe)
            return Xn[8]
        end
        @test ForwardDiff.derivative(h, r.fHe - 5.0e-8) == 0.0 && ForwardDiff.derivative(h, r.fHe - 2.0e-7) == 1.0
    end
    @testset "info" begin
        foreach(kv -> println("Chunk4d-AD ", kv[1], " = ", kv[2]), sort!(collect(OBS4DAD)))
    end
end
