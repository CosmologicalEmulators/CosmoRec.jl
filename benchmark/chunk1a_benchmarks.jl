# Bounded BenchmarkTools measurements for the tested Chunk 1a configuration. Run:
#   julia --project=/home/marcobonici/Desktop/work/CosmologicalEmulators/CosmoRec.jl/benchmark \
#       /home/marcobonici/Desktop/work/CosmologicalEmulators/CosmoRec.jl/benchmark/chunk1a_benchmarks.jl
# Definitions (all BenchmarkTools @benchmark, interpolated args, explicit samples/seconds):
#   primal      : numeric_g at tight tol (1e-12/1e-11), no sensealg, per solver
#   prepare     : prepare_gradient(obj) construction, evals=1, WARMED (compile excluded: a first
#                 call is made before sampling). NOT cold compilation.
#   hot reverse : gradient(obj, prep, backend, theta) with the reused preparation, evals=1
# First-call (compile-inclusive) wall time is not measured here; BenchmarkTools discards it.
using BenchmarkTools, LinearAlgebra, Random, Mooncake, Printf
using SciMLBase: ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqBDF: QNDF
using SciMLSensitivity: GaussAdjoint, MooncakeVJP
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoFiniteDiff
include(joinpath(@__DIR__, "..", "test", "chunk1a_helpers.jl"))

println("Julia $VERSION; threads=$(Threads.nthreads())")
backend = AutoMooncake(; config = nothing)
sa = GaussAdjoint(autojacvec = MooncakeVJP())
w = randn(MersenneTwister(20260930), 6)
algs = ("Rodas5P(autodiff=AutoFiniteDiff)" => Rodas5P(autodiff = AutoFiniteDiff()), "QNDF" => QNDF())
summ(b) = @sprintf("median %.3f ms, min %.3f ms, mean %.3f ms, memory %.1f KiB, allocs %d, samples %d",
    median(b).time / 1e6, minimum(b).time / 1e6, mean(b).time / 1e6, b.memory / 1024, b.allocs, length(b.times))
for (name, alg) in algs
    println("\n## $name")
    primal = @benchmark numeric_g($THETA_NOMINAL, $alg) samples = 200 seconds = 10
    println("primal (hot, no sensealg):             ", summ(primal))
    obj = th -> dot(w, numeric_g(th, alg; sensealg = sa))
    prepare_gradient(obj, backend, THETA_NOMINAL)  # warm-up (compile)
    prep_b = @benchmark prepare_gradient($obj, $backend, $THETA_NOMINAL) evals = 1 samples = 20 seconds = 60
    println("prepare_gradient (warmed, evals=1):    ", summ(prep_b))
    prep = prepare_gradient(obj, backend, THETA_NOMINAL)
    gradient(obj, prep, backend, THETA_NOMINAL)
    hot = @benchmark gradient($obj, $prep, $backend, $THETA_NOMINAL) evals = 1 samples = 100 seconds = 30
    println("hot prepared reverse (evals=1):        ", summ(hot))
    @printf("ratio hot reverse / primal (median):   %.1f\n", median(hot).time / median(primal).time)
end
