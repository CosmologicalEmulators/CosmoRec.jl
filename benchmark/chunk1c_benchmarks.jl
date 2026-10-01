# Bounded BenchmarkTools measurements, accepted Chunk 1c configuration: Rodas5P(autodiff=AutoFiniteDiff()),
# GaussAdjoint(MooncakeVJP) via Mooncake/DifferentiationInterface, abstol 1e-12 / reltol 1e-11, 3 saveat times.
# primal = nl_observe (no sensealg); prepare = WARMED prepare_gradient (evals=1); hot = prepared gradient (evals=1).
# Compile time is not measured. Toy 2-state system; NOT production performance.
using BenchmarkTools, LinearAlgebra, Random, Mooncake, Printf
using SciMLBase: ODEFunction, ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using SciMLSensitivity: GaussAdjoint, MooncakeVJP
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoFiniteDiff
include(joinpath(@__DIR__, "..", "test", "chunk1c_helpers.jl"))

println("Julia $VERSION; threads=$(Threads.nthreads())")
const backend = AutoMooncake(; config = nothing)
const sa = GaussAdjoint(autojacvec = MooncakeVJP())
const alg = Rodas5P(autodiff = AutoFiniteDiff())
const w = randn(MersenneTwister(3), 6)
summ(b) = @sprintf("median %.3f ms, min %.3f ms, mean %.3f ms, memory %.1f KiB, allocs %d, samples %d",
    median(b).time / 1e6, minimum(b).time / 1e6, mean(b).time / 1e6, b.memory / 1024, b.allocs, length(b.times))
nl_observe(THETA1C, alg)
primal = @benchmark nl_observe($THETA1C, $alg) samples = 100 seconds = 10
println("primal (hot, no sensealg):             ", summ(primal))
obj = th -> dot(w, nl_observe(th, alg; sensealg = sa))
prepare_gradient(obj, backend, THETA1C)
prep_b = @benchmark prepare_gradient($obj, $backend, $THETA1C) evals = 1 samples = 15 seconds = 60
println("prepare_gradient (warmed, evals=1):    ", summ(prep_b))
prep = prepare_gradient(obj, backend, THETA1C)
gradient(obj, prep, backend, THETA1C)
hot = @benchmark gradient($obj, $prep, $backend, $THETA1C) evals = 1 samples = 60 seconds = 30
println("hot prepared reverse (evals=1):        ", summ(hot))
@printf("ratio hot reverse / primal (median):   %.1f\n", median(hot).time / median(primal).time)
