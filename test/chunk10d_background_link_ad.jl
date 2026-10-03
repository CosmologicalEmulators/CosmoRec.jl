# Chunk 10d: regression for the explicit background dependence of the public PDE stage (`hi_diffusion_stage` with fixed native population rows) on
# hscale, nbscale. A value/type fast path at hscale = nbscale = 1 preserved the primal but erased these derivatives for Mooncake, which sees plain
# Float64 parameters (chunk10/probe_link_pde.log). Checked at exact unity and nearby values: Float vs Dual primal (bitwise), a non-zero background
# derivative, and the prepared Mooncake VJP vs ForwardDiff. Weights scale each output block by its maximum, so near-zero cancelling outputs do not
# dominate. REQUIRES COSMOREC_NATIVE_DATA_DIR.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(SETUP7) || include("chunk7_helpers.jl")

const MC10D = AutoMooncake(; config = nothing)
const D10D = HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6)
const ROWS10D = NODES5[:, 1:9]
const OBS10D = Dict{String,Float64}()
rec10d(k, v) = (OBS10D[k] = max(get(OBS10D, k, 0.0), v))
function stage10d(p)
    o, _ = hi_diffusion_stage(RM5, D10D, ROWS10D; hscale = p[1], nbscale = p[2])
    m = (o.z .< 2000.0) .& (o.z .> 500.0)
    return vcat(o.DI1_2s[m], o.DF_2g[1][m], o.DF_2g[2][m], o.DF_R[1][m])
end

@testset "Chunk 10d: background (hscale, nbscale) link of the public PDE stage" begin
    f1 = stage10d([1.0, 1.0])
    nb = length(f1) ÷ 4
    blockmax = vcat([fill(maximum(abs, f1[((k - 1) * nb + 1):(k * nb)]), nb) for k in 1:4]...)
    rng = Random.Xoshiro(1004)
    w = randn(rng, length(f1)) ./ blockmax ./ length(f1)
    obj(p, w) = dot(w, stage10d(p))
    prep = prepare_gradient(obj, MC10D, [1.0, 1.0], Constant(w))
    for p in ([1.0, 1.0], [1.001, 1.0], [1.0, 0.999], [0.998, 1.002])
        f = stage10d(p)
        PR = Ref{Vector{Float64}}()
        J = ForwardDiff.jacobian(q -> (o = stage10d(q); PR[] = ForwardDiff.value.(o); o), p)
        @test PR[] == f                                         # Float vs Dual primal, bitwise
        ref = J' * w
        @test all(abs.(ref) .> 0)                               # the background derivative exists
        g = gradient(obj, prep, MC10D, p, Constant(w))
        @test all(isfinite, g) && all(abs.(g) .> 0)             # in particular at exact unity (zero on the old fast path)
        e = maximum(abs.(g .- ref)) / maximum(abs.(ref)); rec10d("VJP_vs_ForwardDiff", e)
        println("Chunk10d p = $p: g = $g  J'w = $ref  rel = $e")
        @test e < 1e-6
    end
    foreach(kv -> println("Chunk10d ", kv[1], " = ", kv[2]), sort!(collect(OBS10D)))
end
