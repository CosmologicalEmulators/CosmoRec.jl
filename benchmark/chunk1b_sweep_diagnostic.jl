# Calibration sweeps (diagnostic, not a test): primal / ForwardDiff / Mooncake errors vs the manufactured oracle.
using Printf, LinearAlgebra, SparseArrays, Random, ForwardDiff, Mooncake
using SciMLBase: ODEFunction, ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqBDF: QNDF
using LinearSolve, Sparspak
using SciMLSensitivity: GaussAdjoint, QuadratureAdjoint, MooncakeVJP
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoFiniteDiff
include(joinpath(@__DIR__, "..", "test", "chunk1b_helpers.jl"))
const TOLS = ((1e-6, 1e-5), (1e-8, 1e-7), (1e-10, 1e-9), (1e-12, 1e-11))
algs = ("Rodas5P(FD)+Sparspak" => Rodas5P(autodiff = AutoFiniteDiff(), linsolve = SparspakFactorization()),
        "QNDF(FD)+Sparspak" => QNDF(autodiff = AutoFiniteDiff(), linsolve = SparspakFactorization()))
for (mname, mk) in (("mesh_a", mesh_a), ("mesh_b", mesh_b))
    ops = MeshOps(mk()); th = THETA1B
    ref = exact_observe(th, ops); Jref = ForwardDiff.jacobian(t -> exact_observe(t, ops), th)
    colmax = [maximum(abs.(Jref[:, j])) for j in 1:6]
    println("== $mname; ref |u| range ", extrema(abs.(ref)), "; column max |J|: ", round.(colmax; sigdigits = 3))
    w = randn(MersenneTwister(5), length(ref))
    for (an, alg) in algs, (at, rt) in TOLS
        sol = solve_rad(th, ops, alg; abstol = at, reltol = rt)
        num = vcat(sol.u...)
        pe = maximum(abs.(num .- ref) ./ (at .+ rt .* abs.(ref)))   # error in units of (atol + rtol|u|)
        J = ForwardDiff.jacobian(t -> observe(t, ops, alg; abstol = at, reltol = rt), th)
        D = abs.(J .- Jref) ./ colmax'
        idx = argmax(D)
        g = try
            obj = t -> dot(w, observe(t, ops, alg; sensealg = GaussAdjoint(autojacvec = MooncakeVJP()), abstol = at, reltol = rt))
            prep = prepare_gradient(obj, AutoMooncake(; config = nothing), th)
            gg = gradient(obj, prep, AutoMooncake(; config = nothing), th)
            r = Jref' * w
            @sprintf("%.2e", maximum(abs.(gg .- r)) / maximum(abs.(r)))
        catch e
            "EXC " * string(typeof(e))
        end
        @printf("%s %s tol=%g: retcode=%s primal err/(atol+rtol|u|)=%.2f  FD(solve) col-scaled=%.2e worst(row %d,col %d)  Mooncake-Gauss proj err=%s\n",
            mname, an, rt, sol.retcode, pe, maximum(D), idx[1], idx[2], g)
    end
end
