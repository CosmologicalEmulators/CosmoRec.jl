# Chunk 4a benchmark: explicit-input Saha initialization + ysol packing (z = 1500 native state), evals = 1, BenchmarkTools only; setup/cold/hot separate.
#   julia --project=<env with CosmoRec + BenchmarkTools + ForwardDiff + Mooncake + DifferentiationInterface> benchmark/chunk4a_benchmarks.jl
using BenchmarkTools, CosmoRec, ForwardDiff, Mooncake, LinearAlgebra
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
include(joinpath(@__DIR__, "..", "test", "chunk4a_helpers.jl"))
const INP = inputs4a(1500.0)
const H = NATIVE_SAHA_HYDROGEN
const HE = NATIVE_SAHA_HELIUM
const P0 = (s = FX4A["SAHAIN"][1500.0]; [FX4A["fHe"], s["Xe_clipped"], s["Xp_raw"], s["NH"], s["Te"], s["TCMB"], s["XHeII"], s["XHeI_ground"]])
fy(p) = pack_ysol(saha_initial_state(SahaInputs(p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]), H, HE), 6, 7)
println("setup: fixture parsing and inputs are outside the timed region")
println("\nsaha_initial_state (15-vector)"); display(@benchmark saha_initial_state($INP, $H, $HE) samples = 5000 evals = 1)
const X0 = saha_initial_state(INP, H, HE)
println("\npack_ysol (12-vector)"); display(@benchmark pack_ysol($X0, 6, 7) samples = 5000 evals = 1)
println("\nfull public path: initialize + pack"); display(@benchmark pack_ysol(saha_initial_state($INP, $H, $HE), 6, 7) samples = 5000 evals = 1)
println("\nForwardDiff.jacobian (12 x 8)"); display(@benchmark ForwardDiff.jacobian($fy, $P0) samples = 500 evals = 1)
obj(x, w) = dot(w, fy(x)); w = collect(range(0.3, 1.7; length = 12)); backend = AutoMooncake(; config = nothing)
println("\nMooncake prepare_gradient (cold, includes compilation)"); display(@benchmark prepare_gradient($obj, $backend, $P0, Constant($w)) samples = 3 evals = 1)
prep = prepare_gradient(obj, backend, P0, Constant(w)); gradient(obj, prep, backend, P0, Constant(w))
println("\nMooncake hot gradient (prepared, reused)"); display(@benchmark gradient($obj, $prep, $backend, $P0, Constant($w)) samples = 500 evals = 1)
println()
