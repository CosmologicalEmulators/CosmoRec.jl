# Chunk 6a AD tests: derivatives of the PDE coefficients (pd, Dnem at selected z) with respect to all stored population values of the previous pass (rows: Xe, X1s ... X3d, rho;
# 3000 x 8 = 24000 inputs): ForwardDiff directional derivative vs 256-bit central difference, prepared Mooncake VJP (two independent preparations, changed inputs) vs the ForwardDiff
# gradient. Piecewise: spline cell boundaries (C2), the coefficient-grid size (a discrete count of stored rows) and the rate-table stencils are frozen/primal; no claim across them.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(FX6A) || include("chunk6_helpers.jl")

const MC6A = AutoMooncake(; config = nothing)
const OBS6AAD = Dict{String,Float64}()
rec6aad(k, v) = (OBS6AAD[k] = max(get(OBS6AAD, k, 0.0), v))
const ROWS6 = NODES5[:, 1:9]
const ZQ6 = [2300.0, 1700.0, 1300.0, 900.0, 600.0]

function coeff_out(x)
    rows = hcat(ROWS6[:, 1], reshape(x, size(ROWS6, 1), 8))     # z column promoted to the AD type; not differentiated (primal values used)
    pops = HIPopulationSplines(rows; zs = ZS6, ze = ZE6)
    cs = hi_pde_coefficients(rows, pops, ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6)
    T = eltype(x)
    out = Vector{T}(undef, 8 * length(ZQ6))
    k = 0
    for z in ZQ6
        for m in (1, 2, 4); k += 1; out[k] = log(hi_pd(cs, z, m)); end       # pd of 3s, 3d is identically 1
        for m in 1:5; k += 1; out[k] = hi_Dnem(cs, z, m); end
    end
    return out
end
const X06 = vec(ROWS6[:, 2:9])

@testset "Chunk 6a: PDE coefficients AD" begin
    @testset "directional derivative vs 256-bit central difference" begin
        rng = Random.Xoshiro(20261101)
        v = randn(rng, length(X06)) .* abs.(X06)
        dd = ForwardDiff.derivative(t -> coeff_out(X06 .+ t .* v), 0.0)
        ref = setprecision(BigFloat, 256) do
            h = BigFloat(1.0e-20); xb = BigFloat.(X06); vb = BigFloat.(v)
            Float64.((coeff_out(xb .+ h .* vb) .- coeff_out(xb .- h .* vb)) ./ (2h))
        end
        @test all(isfinite, dd)
        e = maximum(abs.(dd .- ref) ./ max.(abs.(ref), 1e-10 * maximum(abs, ref))); rec6aad("DIR", e)
        @test e < 1e-8
    end
    @testset "prepared Mooncake VJP vs ForwardDiff gradient" begin
        rng = Random.Xoshiro(20261102)
        nout = length(coeff_out(X06))
        obj(x, w) = dot(w, coeff_out(x))
        prepA = prepare_gradient(obj, MC6A, X06, Constant(randn(rng, nout)))
        prepB = prepare_gradient(obj, MC6A, X06, Constant(randn(rng, nout)))
        for kk in 1:2
            x = kk == 1 ? X06 : X06 .* (1 .+ 1.0e-7 .* randn(rng, length(X06)))
            w = randn(rng, nout)
            ref = ForwardDiff.gradient(q -> obj(q, w), x)
            gA = gradient(obj, prepA, MC6A, x, Constant(w)); gB = gradient(obj, prepB, MC6A, x, Constant(w))
            sc = maximum(abs, ref)
            e = max(maximum(abs.(gA .- ref)), maximum(abs.(gB .- ref))) / sc; rec6aad("VJP", e)
            @test e < 1e-10
            @test gA == gB
        end
    end
    @testset "info" begin
        foreach(kv -> println("Chunk6a-AD ", kv[1], " = ", kv[2]), sort!(collect(OBS6AAD)))
    end
end
