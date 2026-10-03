# Chunk 4d benchmark: sampled-HeI switch decision/reset/packing (BenchmarkTools only; evals = 1; setup outside the timed region).
#   julia --project=<env with CosmoRec + BenchmarkTools + ForwardDiff + Mooncake + DifferentiationInterface> benchmark/chunk4d_benchmarks.jl
using BenchmarkTools, CosmoRec, ForwardDiff, Mooncake, LinearAlgebra
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
include(joinpath(@__DIR__, "..", "test", "chunk4d_helpers.jl"))
const R = first(r for r in FX4D if r.cond && r.zs == 2000.0)
const X0 = copy(R.Xpre)
println("helium_switch_condition"); display(@benchmark helium_switch_condition($(R.zs), $(R.fHe), $(X0[8])) samples = 5000 evals = 1)
println("\nhelium_switch_reset"); display(@benchmark helium_switch_reset($X0, 6, 7, $(R.fHe)) samples = 5000 evals = 1)
println("\ndecision + reset + packing (switched)"); display(@benchmark pack_ysol(first(helium_switch($X0, $(R.zs), 6, 7, $(R.fHe)); ), 6, 7; helium = false) samples = 5000 evals = 1)
x0 = vcat(X0, R.fHe)
f(x) = pack_ysol(helium_switch_reset(x[1:15], 6, 7, x[16]), 6, 7; helium = false)
println("\nForwardDiff.jacobian (7 x 16)"); display(@benchmark ForwardDiff.jacobian($f, $x0) samples = 500 evals = 1)
w = collect(range(0.5, 1.5; length = 7)); obj(x, w) = dot(w, f(x)); backend = AutoMooncake(; config = nothing)
println("\nMooncake prepare_gradient (cold, includes compilation)"); display(@benchmark prepare_gradient($obj, $backend, $x0, Constant($w)) samples = 3 evals = 1)
prep = prepare_gradient(obj, backend, x0, Constant(w)); gradient(obj, prep, backend, x0, Constant(w))
println("\nMooncake hot gradient (prepared, reused)"); display(@benchmark gradient($obj, $prep, $backend, $x0, Constant($w)) samples = 500 evals = 1)
println()
