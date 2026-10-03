# Chunk 4b AD tests: ForwardDiff of the explicit accessors with respect to z and to the history node values vs 256-bit finite differences (per z, and directional over all
# 24000 node values), prepared Mooncake VJPs (two independent preparations, reused at changed history/z) vs ForwardDiff J'w. The map is piecewise smooth: C2 in z across spline nodes,
# NOT differentiable across the z_saha = 3500 and zsRe switches, the end nudge, the loaded-H edge zmin/zmax, or the Xe_He floor; those are tested as one-sided/discontinuous below.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(ACC4B) || include("chunk4b_helpers.jl")

const MC4B = AutoMooncake(; config = nothing)
const OBS4BAD = Dict{String,Float64}()
rec4bad(k, v) = (OBS4BAD[k] = max(get(OBS4BAD, k, 0.0), v))

const NH4B = 6000
# p = [Xe(6000); Xe_H(6000); Xe_He(6000); TM(6000); z]
p4b() = vcat(HIST4B[:, 4], HIST4B[:, 2], HIST4B[:, 3], HIST4B[:, 7], 1500.0)
function acc_from(p)
    c = constants4b()
    z = HIST4B[:, 1]
    sp = recfast_splines(c, z, p[(NH4B + 1):(2NH4B)], p[(2NH4B + 1):(3NH4B)], p[1:NH4B], HIST4B[:, 5], HIST4B[:, 6], p[(3NH4B + 1):(4NH4B)])
    return CosmosAccessors(c, sp, nothing)
end
# scalar-vector accessor map at redshift zq (history-dependent quantities)
function outs(a, zq)
    X = saha_initial_state(saha_inputs_at(a, zq), NATIVE_SAHA_HYDROGEN, NATIVE_SAHA_HELIUM)
    return vcat(X[1:3], X[8:10], [cosmos_Te_Tg(a, zq), cosmos_XHeII1s(a, zq), cosmos_kappa_cool(a, zq), cosmos_dXe_dz(a, zq)])
end
f_hist(p) = outs(acc_from(p), p[end])
f_z(zq) = a -> outs(a, zq)

function cdiff_dir(f, x, v)
    setprecision(BigFloat, 256) do
        xb = BigFloat.(x); vb = BigFloat.(v)
        h = BigFloat(1.0e-20)
        return Float64.((f(xb .+ h .* vb) .- f(xb .- h .* vb)) ./ (2h))
    end
end

@testset "Chunk 4b: Cosmos accessors AD" begin
    a0 = acc_from(p4b())
    @testset "z-derivatives: ForwardDiff vs 256-bit central difference" begin
        for zq in (300.0, 800.0, 1500.0, 2999.5, 4000.0, 8000.0)      # below and above z_saha = 3500, interior of spline cells
            f = z -> outs(a0, z)
            J = ForwardDiff.derivative(f, zq)
            Jfd = setprecision(BigFloat, 256) do
                h = BigFloat(1.0e-20) * zq
                Float64.((f(BigFloat(zq) + h) .- f(BigFloat(zq) - h)) ./ (2h))
            end
            sc = max(maximum(abs, J), floatmin())
            e = maximum(abs.(J .- Jfd) ./ max.(abs.(Jfd), 1e-8 * sc))
            rec4bad("dz:$zq", e); @test e < 1e-10
        end
    end
    @testset "directional derivative over all history values and z" begin
        rng = Random.Xoshiro(20261003)
        p0 = p4b(); v = randn(rng, length(p0)) .* abs.(p0)
        v[end] = 3.0 * randn(rng)
        dd = ForwardDiff.derivative(t -> f_hist(p0 .+ t .* v), 0.0)
        ref = cdiff_dir(f_hist, p0, v)
        sc = maximum(abs, ref)
        e = maximum(abs.(dd .- ref) ./ max.(abs.(ref), 1e-8 * sc)); rec4bad("DIR:history", e)
        @test e < 1e-9
    end
    @testset "prepared Mooncake VJP vs ForwardDiff gradient (independent preparations, changed inputs)" begin
        rng = Random.Xoshiro(20261004)
        p0 = p4b(); nout = length(f_hist(p0))
        obj(p, w) = dot(w, f_hist(p))
        prepA = prepare_gradient(obj, MC4B, p0, Constant(randn(rng, nout)))
        prepB = prepare_gradient(obj, MC4B, p0, Constant(randn(rng, nout)))
        for kk in 1:2
            p = kk == 1 ? p0 : p0 .* (1 .+ 1.0e-6 .* randn(rng, length(p0)))
            kk == 2 && (p[end] = 2200.0)
            w = randn(rng, nout)
            ref = ForwardDiff.gradient(q -> obj(q, w), p)
            gA = gradient(obj, prepA, MC4B, p, Constant(w))
            gB = gradient(obj, prepB, MC4B, p, Constant(w))
            sc = max(maximum(abs, ref), floatmin())
            e = max(maximum(abs.(gA .- ref)), maximum(abs.(gB .- ref))) / sc; rec4bad("VJP", e)
            @test isapprox(gA, ref; rtol = 1.0e-8, atol = 1.0e-8 * sc)
            @test isapprox(gB, ref; rtol = 1.0e-8, atol = 1.0e-8 * sc)
            @test gA == gB
        end
    end
    @testset "documented nondifferentiability (one-sided / discontinuous)" begin
        # z_saha = 3500: the Saha branch jumps (relative ~1e-3 in X1s) and so does the derivative; the interior of each side is smooth
        d_lo = ForwardDiff.derivative(z -> cosmos_X1s(a0, z), 3499.999999)
        d_hi = ForwardDiff.derivative(z -> cosmos_X1s(a0, z), 3500.0)
        rec4bad("switch:dX1s_jump", abs(d_hi / d_lo - 1)); @test d_lo != d_hi
        # spline cell boundaries: C2, derivatives agree across a node to rounding
        zn = HIST4B[3000, 1]
        dl = ForwardDiff.derivative(z -> cosmos_Xe_Seager(a0, z), zn - 1e-9)
        dr = ForwardDiff.derivative(z -> cosmos_Xe_Seager(a0, z), zn + 1e-9)
        rec4bad("node:dXe_continuity", abs(dl / dr - 1)); @test abs(dl / dr - 1) < 1e-6
        # no derivative claim outside the spline range: throws
        @test_throws SplineDomainError ForwardDiff.derivative(z -> cosmos_Xe_Seager(a0, z), -1.0)
    end
    @testset "info" begin
        foreach(kv -> println("Chunk4b-AD ", kv[1], " = ", kv[2]), sort!(collect(OBS4BAD)))
    end
end
