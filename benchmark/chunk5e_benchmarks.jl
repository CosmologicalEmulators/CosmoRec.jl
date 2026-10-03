# Chunk 5e benchmark: derivatives of the single pass (BenchmarkTools only; cold and hot Mooncake reported separately). REQUIRES COSMOREC_NATIVE_DATA_DIR.
using BenchmarkTools, CosmoRec, ForwardDiff, Mooncake, LinearAlgebra
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
include(joinpath(@__DIR__, "..", "test", "chunk5e_helpers.jl"))
println("primal complete pipeline (frozen switch, 3000 nodes, 9-node grid output)"); display(@benchmark full_out($P0_5E) samples = 5 evals = 1)
println("\nForwardDiff.jacobian of the complete pipeline (18 outputs x 4 parameters)"); display(@benchmark ForwardDiff.jacobian(full_out, $P0_5E) samples = 3 evals = 1)
f_fd = p -> sel_out(p, :fd); f_rv = p -> sel_out(p, :rev)
println("\nprimal selected-node pipeline (forward route)"); display(@benchmark $f_fd($P0_5E) samples = 5 evals = 1)
println("\nForwardDiff.jacobian selected-node pipeline"); display(@benchmark ForwardDiff.jacobian($f_fd, $P0_5E) samples = 3 evals = 1)
n = length(f_fd(P0_5E)); w = collect(range(0.5, 1.5; length = n)); obj(p, w) = dot(w, f_rv(p)); backend = AutoMooncake(; config = nothing)
println("\nMooncake prepare_gradient (cold, includes compilation; first sample)"); display(@benchmark prepare_gradient($obj, $backend, $P0_5E, Constant($w)) samples = 1 evals = 1)
prep = prepare_gradient(obj, backend, P0_5E, Constant(w)); gradient(obj, prep, backend, P0_5E, Constant(w))
println("\nMooncake hot gradient (prepared, reused)"); display(@benchmark gradient($obj, $prep, $backend, $P0_5E, Constant($w)) samples = 3 evals = 1)
println()
