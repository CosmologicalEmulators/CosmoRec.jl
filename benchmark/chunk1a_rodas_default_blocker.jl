# Minimal reproducer: Rodas5P with default solver-internal autodiff under Mooncake + MooncakeVJP adjoint.
using LinearAlgebra: dot
using Random, Mooncake
using SciMLBase: ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using SciMLSensitivity: GaussAdjoint, MooncakeVJP
using DifferentiationInterface: prepare_gradient
using ADTypes: AutoMooncake
include(joinpath(@__DIR__, "..", "test", "chunk1a_helpers.jl"))
w = randn(MersenneTwister(1), 6)
obj = th -> dot(w, numeric_g(th, Rodas5P(); sensealg = GaussAdjoint(autojacvec = MooncakeVJP())))
try
    prepare_gradient(obj, AutoMooncake(; config = nothing), THETA_NOMINAL)
    println("NO ERROR (blocker not reproduced)")
catch e
    println("reproduced: ", typeof(e), "\n", first(sprint(showerror, e), 300))
    frames = unique(string.(getfield.(stacktrace(catch_backtrace()), :func)))
    println("frames (outermost last, first 40 unique): ", join(first(frames, 40), " <- "))
end
