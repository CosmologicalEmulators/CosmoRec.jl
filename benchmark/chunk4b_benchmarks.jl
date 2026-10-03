# Chunk 4b benchmark: explicit Cosmos accessors (BenchmarkTools only, evals = 1 where mutating/allocating setup is involved).
#   julia --project=<env with CosmoRec + BenchmarkTools + ForwardDiff + Mooncake + DifferentiationInterface> benchmark/chunk4b_benchmarks.jl
using BenchmarkTools, CosmoRec, ForwardDiff, Mooncake, LinearAlgebra
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
include(joinpath(@__DIR__, "..", "test", "chunk4b_helpers.jl"))
const A = accessors4b()
const C = constants4b()
println("setup: spline construction from the 6000-node history (the 6 splines)")
display(@benchmark splines4b($C, $HIST4B) samples = 200 evals = 1)
println("\nHubble table (10000 nodes)"); display(@benchmark hubble_table($(HUB4B[:, 1]), $(HUB4B[:, 2])) samples = 200 evals = 1)
println("\nspline_eval_native, Xe_Seager(1500)"); display(@benchmark cosmos_Xe_Seager($A, 1500.0) samples = 5000 evals = 1)
println("\nH(z) loaded table"); display(@benchmark cosmos_H($A, 1500.0) samples = 5000 evals = 1)
println("\nsaha_inputs_at (Xe, Xp, NH, Te, NHe*)"); display(@benchmark saha_inputs_at($A, 1500.0) samples = 5000 evals = 1)
println("\nsaha_inputs_at + saha_initial_state + pack_ysol"); display(@benchmark pack_ysol(saha_initial_state(saha_inputs_at($A, 1500.0), NATIVE_SAHA_HYDROGEN, NATIVE_SAHA_HELIUM), 6, 7) samples = 5000 evals = 1)
const NH = 6000
p0 = vcat(HIST4B[:, 4], HIST4B[:, 2], HIST4B[:, 3], HIST4B[:, 7], 1500.0)
function acc_from(p)
    sp = recfast_splines(C, HIST4B[:, 1], p[(NH + 1):(2NH)], p[(2NH + 1):(3NH)], p[1:NH], HIST4B[:, 5], HIST4B[:, 6], p[(3NH + 1):(4NH)])
    return CosmosAccessors(C, sp, nothing)
end
function f(p)
    a = acc_from(p)
    X = saha_initial_state(saha_inputs_at(a, p[end]), NATIVE_SAHA_HYDROGEN, NATIVE_SAHA_HELIUM)
    return vcat(X[1:3], X[8:10])
end
println("\nForwardDiff.gradient of w'f over 24001 inputs (cold first sample included)"); w = collect(1.0:6.0)
display(@benchmark ForwardDiff.gradient(q -> dot($w, f(q)), $p0) samples = 3 evals = 1)
obj(x, w) = dot(w, f(x)); backend = AutoMooncake(; config = nothing)
println("\nMooncake prepare_gradient (cold, includes compilation)"); display(@benchmark prepare_gradient($obj, $backend, $p0, Constant($w)) samples = 3 evals = 1)
prep = prepare_gradient(obj, backend, p0, Constant(w)); gradient(obj, prep, backend, p0, Constant(w))
println("\nMooncake hot gradient (prepared, reused)"); display(@benchmark gradient($obj, $prep, $backend, $p0, Constant($w)) samples = 50 evals = 1)
println()
