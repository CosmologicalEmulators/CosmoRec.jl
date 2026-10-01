# Chunk 3d benchmark: assembled pointwise fcn_effective (explicit model/background, evals = 1). BenchmarkTools only; no @time/@elapsed.
#   julia --project=<test_env with CosmoRec + BenchmarkTools + ForwardDiff + Mooncake + DifferentiationInterface> benchmark/chunk3d_benchmarks.jl
using BenchmarkTools, CosmoRec, ForwardDiff, Mooncake, LinearAlgebra
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
include(joinpath(@__DIR__, "..", "test", "chunk3d_helpers.jl"))
const S = first(s for s in FX3D.states if s.label == "z2600")
const BG = background3d(S)
const X0 = copy(S.X)
const P0 = vcat(S.X, S.Tg, S.NH, S.Hz, S.fHe)
const GBUF = zeros(15)
fcn!(g, X) = fcn_effective!(g, S.z, X, BG, MODEL3D_ON)
fvec(p) = fcn_effective(S.z, view(p, 1:15), EffectiveBackground(p[16], p[17], p[18], p[19]), MODEL3D_ON)
fvec_off(p) = fcn_effective(S.z, view(p, 1:15), EffectiveBackground(p[16], p[17], p[18], p[19]), MODEL3D_OFF)
fcn!(GBUF, X0)   # warm-up
println("fcn_effective! (H + He + rho + absorber, in place)"); display(@benchmark fcn!(gg, $X0) setup = (gg = zeros(15)) samples = 2000 evals = 1)
println("\nfcn_effective! (absorber off)"); display(@benchmark fcn_effective!(gg, $(S.z), $X0, $BG, $MODEL3D_OFF) setup = (gg = zeros(15)) samples = 2000 evals = 1)
println("\nForwardDiff.jacobian (15 outputs x 19 inputs)"); display(@benchmark ForwardDiff.jacobian($fvec, $P0) samples = 500 evals = 1)
obj(x, w) = dot(w, fvec(x)); w = collect(range(0.3, 1.7; length = 15)); backend = AutoMooncake(; config = nothing)
println("\nMooncake prepare_gradient (first sample includes compilation)"); display(@benchmark prepare_gradient($obj, $backend, $P0, Constant($w)) samples = 3 evals = 1)
prep = prepare_gradient(obj, backend, P0, Constant(w)); gradient(obj, prep, backend, P0, Constant(w))   # warm
println("\nMooncake hot gradient (prepared, reused)"); display(@benchmark gradient($obj, $prep, $backend, $P0, Constant($w)) samples = 500 evals = 1)
println()
