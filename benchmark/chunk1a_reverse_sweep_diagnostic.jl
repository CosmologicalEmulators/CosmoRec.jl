# Reverse (Mooncake + continuous adjoint) error vs analytic oracle across tolerances/solvers/adjoints.
using Printf, ForwardDiff, LinearAlgebra, Random, Mooncake
using SciMLBase: ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqBDF: QNDF
using SciMLSensitivity: GaussAdjoint, QuadratureAdjoint, MooncakeVJP
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoFiniteDiff
include(joinpath(@__DIR__, "..", "test", "chunk1a_helpers.jl"))
backend = AutoMooncake(; config = nothing)
w = randn(MersenneTwister(20260930), 6)
thetas = [THETA_NOMINAL, [1200.0, 0.8, 40.0, 1.3, 0.7], [900.0, 1.5, 60.0, 0.6, 1.4]]
for (aname, alg) in (("Rodas5P(autodiff=AutoFiniteDiff)", Rodas5P(autodiff = AutoFiniteDiff())), ("QNDF", QNDF()))
  for (sname, sa) in (("Gauss", GaussAdjoint(autojacvec = MooncakeVJP())), ("Quadrature", QuadratureAdjoint(autojacvec = MooncakeVJP())))
    for (at, rt) in ((1e-6, 1e-5), (1e-8, 1e-7), (1e-10, 1e-9), (1e-12, 1e-11))
        res = try
            obj = th -> dot(w, numeric_g(th, alg; sensealg = sa, abstol = at, reltol = rt))
            prep = prepare_gradient(obj, backend, THETA_NOMINAL)
            errs = map(thetas) do th
                g = gradient(obj, prep, backend, th)
                ref = ForwardDiff.jacobian(analytic_g, th)' * w
                maximum(abs.(g .- ref)) / maximum(abs.(ref))
            end
            "max-norm-scaled err vs analytic over 3 thetas: " * join([@sprintf("%.2e", e) for e in errs], ", ")
        catch e
            "EXCEPTION " * first(sprint(showerror, e), 200)
        end
        println("$aname | $sname | abstol=$at reltol=$rt | $res")
    end
  end
end
