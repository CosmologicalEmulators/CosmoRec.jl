# Chunk 5a benchmark: ODE right-hand side (BenchmarkTools only; explicit arguments; setup = table loading and model construction outside the timed region).
# REQUIRES COSMOREC_NATIVE_DATA_DIR (see test/chunk5_helpers.jl).
using BenchmarkTools, CosmoRec, ForwardDiff, Mooncake, LinearAlgebra
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
println("setup: loading the production tables (cold)"); t0 = time_ns(); include(joinpath(@__DIR__, "..", "test", "chunk5_helpers.jl")); println("  table load + model construction: ", (time_ns() - t0) / 1.0e9, " s (single cold measurement, includes compilation)")
include(joinpath(@__DIR__, "..", "test", "chunk5a_ode_rhs_native.jl"))
const R1 = first(r for r in ROWS5A if r.flag && r.z == 1500.0 && r.pat == 0)
const R0 = first(r for r in ROWS5A if !r.flag && r.z == 400.0 && r.pat == 0)
f1 = zeros(12); f0 = zeros(7)
println("\nrecombination_rhs! flag_He = 1 (z = 1500, off-table DP fallback)"); display(@benchmark recombination_rhs!($f1, $(R1.z), $(R1.y), $RM5; flag_He = true) samples = 500 evals = 1)
println("\nrecombination_rhs! flag_He = 0 (z = 400)"); display(@benchmark recombination_rhs!($f0, $(R0.z), $(R0.y), $RM5; flag_He = false) samples = 2000 evals = 1)
println("\nForwardDiff.jacobian 12 x 12 (flag_He = 1)"); display(@benchmark ForwardDiff.jacobian($(y -> recombination_rhs(R1.z, y, RM5)), $(R1.y)) samples = 100 evals = 1)
w = collect(range(0.5, 1.5; length = 12)); obj(y, w) = dot(w, recombination_rhs(R1.z, y, RM5)); backend = AutoMooncake(; config = nothing)
println("\nMooncake prepare_gradient (cold, includes compilation)"); display(@benchmark prepare_gradient($obj, $backend, $(R1.y), Constant($w)) samples = 2 evals = 1)
prep = prepare_gradient(obj, backend, R1.y, Constant(w)); gradient(obj, prep, backend, R1.y, Constant(w))
println("\nMooncake hot gradient (prepared, reused)"); display(@benchmark gradient($obj, $prep, $backend, $(R1.y), Constant($w)) samples = 100 evals = 1)
println()
