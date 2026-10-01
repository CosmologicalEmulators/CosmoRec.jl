# Bounded BenchmarkTools measurements for the ACCEPTED Chunk 1b configuration:
#   Rodas5P(autodiff=AutoFiniteDiff, linsolve=SparspakFactorization) + analytic sparse jac,
#   GaussAdjoint(MooncakeVJP) via Mooncake/DifferentiationInterface, tolerances abstol 1e-10 / reltol 1e-9.
# Run: julia --project=<repo>/benchmark <repo>/benchmark/chunk1b_benchmarks.jl
# primal = observe (no sensealg); prepare = WARMED prepare_gradient (evals=1); hot = prepared gradient (evals=1).
# Toy 7-11 unknown system; these are NOT production performance numbers. Compile time is not measured.
using BenchmarkTools, LinearAlgebra, SparseArrays, Random, Mooncake, Printf
using SciMLBase: ODEFunction, ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using LinearSolve, Sparspak
using SciMLSensitivity: GaussAdjoint, MooncakeVJP
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoFiniteDiff
include(joinpath(@__DIR__, "..", "test", "chunk1b_helpers.jl"))

println("Julia $VERSION; threads=$(Threads.nthreads())")
backend = AutoMooncake(; config = nothing)
sa = GaussAdjoint(autojacvec = MooncakeVJP())
alg = Rodas5P(autodiff = AutoFiniteDiff(), linsolve = SparspakFactorization())
summ(b) = @sprintf("median %.3f ms, min %.3f ms, mean %.3f ms, memory %.1f KiB, allocs %d, samples %d",
    median(b).time / 1e6, minimum(b).time / 1e6, mean(b).time / 1e6, b.memory / 1024, b.allocs, length(b.times))
for (mname, m) in (("mesh_a (9 points, 7 unknowns)", mesh_a()), ("mesh_b (12 points, 10 unknowns)", mesh_b()))
    ops = MeshOps(m)
    println("\n## $mname")
    w = randn(MersenneTwister(20260930), 3 * (ops.N - 1))
    observe(THETA1B, ops, alg)
    primal = @benchmark observe($THETA1B, $ops, $alg) samples = 100 seconds = 10
    println("primal (hot, no sensealg):             ", summ(primal))
    obj = th -> dot(w, observe(th, ops, alg; sensealg = sa))
    prepare_gradient(obj, backend, THETA1B)
    prep_b = @benchmark prepare_gradient($obj, $backend, $THETA1B) evals = 1 samples = 15 seconds = 60
    println("prepare_gradient (warmed, evals=1):    ", summ(prep_b))
    prep = prepare_gradient(obj, backend, THETA1B)
    gradient(obj, prep, backend, THETA1B)
    hot = @benchmark gradient($obj, $prep, $backend, $THETA1B) evals = 1 samples = 60 seconds = 30
    println("hot prepared reverse (evals=1):        ", summ(hot))
    @printf("ratio hot reverse / primal (median):   %.1f\n", median(hot).time / median(primal).time)
end
