# Reduced reproducers for non-accepted candidates on the Chunk 1b problem (diagnostic log only).
using LinearAlgebra, SparseArrays, Random, ForwardDiff, Mooncake
using SciMLBase: ODEFunction, ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqBDF: QNDF
using LinearSolve, Sparspak
using SciMLSensitivity: GaussAdjoint, MooncakeVJP
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoFiniteDiff
include(joinpath(@__DIR__, "..", "test", "chunk1b_helpers.jl"))
ops = MeshOps(mesh_a())
w = randn(MersenneTwister(1), 3 * (ops.N - 1))
sa = GaussAdjoint(autojacvec = MooncakeVJP())
cases = ("QNDF+Sparspak, Mooncake/Gauss" => QNDF(linsolve = SparspakFactorization()),
         "QNDF default linsolve, Mooncake/Gauss" => QNDF(),
         "QNDF(autodiff=AutoFiniteDiff)+Sparspak, Mooncake/Gauss" => QNDF(autodiff = AutoFiniteDiff(), linsolve = SparspakFactorization()),
         "Rodas5P default autodiff + Sparspak, Mooncake/Gauss" => Rodas5P(linsolve = SparspakFactorization()))
for (label, alg) in cases
    try
        obj = th -> dot(w, observe(th, ops, alg; sensealg = sa))
        prep = prepare_gradient(obj, AutoMooncake(; config = nothing), THETA1B)
        gradient(obj, prep, AutoMooncake(; config = nothing), THETA1B)
        println(label, " => OK")
    catch e
        frames = unique(string.(getfield.(stacktrace(catch_backtrace()), :func)))
        println(label, " => EXCEPTION ", typeof(e), ": ", first(replace(sprint(showerror, e), "\n" => " "), 400))
        println("    frames: ", join(first(frames, 14), " <- "))
    end
end
println("ForwardDiff with default-autodiff Rodas5P + Sparspak:")
try
    alg = Rodas5P(linsolve = SparspakFactorization())
    J = ForwardDiff.jacobian(t -> observe(t, ops, alg), THETA1B)
    Jr = ForwardDiff.jacobian(t -> exact_observe(t, ops), THETA1B)
    println("  OK, max abs err vs analytic = ", maximum(abs.(J .- Jr)))
catch e
    println("  EXCEPTION ", typeof(e), ": ", first(replace(sprint(showerror, e), "\n" => " "), 300))
end
