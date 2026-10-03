# Chunk 10c: prepared Mooncake gradient (outer reverse mode; INNER ForwardDiffSensitivity in every ODE solve, i.e. NOT true reverse through the
# Rodas5P steps, which the LinearSolve 5.18.2 / Mooncake 0.5.61 cached-solve rule blocks, see docs/STATUS_PHASES7_10.md) of the COMPOSED production objective (3 ODE passes, 2 HI PDE stages with dynamic feedback passed to the solver
# as packed spline coefficients, Recfast tail, output assembly on ZG5E) w.r.t. p = [F, A2s1s, hscale, nbscale], vs the ForwardDiff Jacobian of the
# IDENTICAL map (same solver routine, tolerances and explicit sensealg). Accuracy assertions (finite values, dot test and J'w comparison).
# FROZEN, not differentiated: k_switch, the Saha initial state and the preliminary history inside the Cosmos accessors, the grids, the Patterson stopping
# decisions and the solver step sequence (the reverse pass differentiates the realised steps). REQUIRES COSMOREC_NATIVE_DATA_DIR.
using Test
using Random
using LinearAlgebra
using CosmoRec
@isdefined(composed10r) || include("chunk10_reverse_helpers.jl")
using SciMLSensitivity: ForwardDiffSensitivity
SENSEALG10R[] = ForwardDiffSensitivity(); FEEDBACK_PACKING10R[] = :nodes

const OBS10C = Dict{String,Float64}()
rec10c(k, v) = (OBS10C[k] = max(get(OBS10C, k, 0.0), v))
@testset "Chunk 10c: prepared Mooncake gradient of the composed runmode-0 objective" begin
    f(p) = composed10r(p; iterations = 2)
    f0 = f(P0_5E)
    @test all(isfinite, f0) && length(f0) == 2 * length(ZG5E)
    Jf = ForwardDiff.jacobian(f, P0_5E)
    @test all(isfinite, Jf)
    obj(p, w) = dot(w, f(p))
    rng = Random.Xoshiro(1010)
    w0 = randn(rng, length(f0)) ./ abs.(f0)
    prepA = prepare_gradient(obj, MC5E, P0_5E, Constant(w0))
    prepB = prepare_gradient(obj, MC5E, P0_5E, Constant(randn(rng, length(f0)) ./ abs.(f0)))
    for (k, p) in enumerate((P0_5E, P0_5E .* [1.0, 1.0, 1.0 + 1e-4, 1.0 - 1e-4]))
        Jp = k == 1 ? Jf : ForwardDiff.jacobian(f, p)
        for j in 1:2
            w = randn(rng, length(f0)) ./ abs.(f0)
            gA = gradient(obj, prepA, MC5E, p, Constant(w)); gB = gradient(obj, prepB, MC5E, p, Constant(w))
            @test all(isfinite, gA) && gA == gB
            ref = Jp' * w
            e = maximum(abs.(gA .- ref)) / maximum(abs.(ref)); rec10c("Jtw_rel", e)
            println("Chunk10c point $k block $j: g = $gA  J'w = $ref")
            @test e < 1e-6
        end
    end
    foreach(kv -> println("Chunk10c ", kv[1], " = ", kv[2]), sort!(collect(OBS10C)))
end
