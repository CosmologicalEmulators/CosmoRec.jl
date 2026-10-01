# Chunk 1a: SciML stiff-solve <-> AD compatibility probe (toy problem, analytic oracle).
# NO CosmoRec physics. Native-code text fixtures are NOT part of this chunk; the
# future fixture contract is docs/fixture_contract.md. Results: docs/CHUNK1A_RESULTS.md.
#
# Gates (all against the closed-form oracle `analytic_g`, never only against another AD path):
#   1. primal solve vs analytic, tolerance sweep;
#   2. ForwardDiff through the full public `solve` vs ForwardDiff of the closed form,
#      with a column-scaled absolute error, a tolerance sweep and seeded projections;
#   3. central finite differences: a gate only where the solve output is smooth in theta
#      (Rodas5P); for QNDF they are reported as diagnostics because the step-size/order
#      selection makes the numerical output non-smooth at FD step scales;
#   4. Mooncake (DifferentiationInterface prepare_gradient, reused) through the public
#      `solve(...; sensealg=GaussAdjoint/QuadratureAdjoint(autojacvec=MooncakeVJP()))`
#      with theta-dependent u0/p, tolerance sweep vs the analytic gradient, cache reuse.
# That the MooncakeVJP continuous-adjoint hook is actually reached is shown by the
# process-local sentinel diagnostic benchmark/chunk1a_route_diagnostic.jl (not here,
# so that no method override exists in the package tests).

using Test
using LinearAlgebra: dot
using Random
using ForwardDiff
using SciMLBase: ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqBDF: QNDF
using SciMLSensitivity: GaussAdjoint, QuadratureAdjoint, MooncakeVJP
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoFiniteDiff

include("chunk1a_helpers.jl")

# Rodas5P with default solver-internal autodiff fails inside the continuous adjoint's BACKWARD
# solve (stack: adjoint_sensitivity_backpass -> _adjoint_sensitivities -> solve -> Rosenbrock
# calc_tderivative!; reproducer benchmark/chunk1a_rodas_default_blocker.jl): the adjoint RHS
# contains the Mooncake VJP, which ForwardDiff cannot differentiate for the Rosenbrock time
# gradient. The solver's own Jacobian/time-gradient is therefore set to AutoFiniteDiff;
# Mooncake remains the outer AD and MooncakeVJP the RHS-VJP engine.
# CAVEAT: the toy RHS is linear, so a finite-difference Jacobian is exact up to rounding;
# this does not demonstrate adequacy for a nonlinear RHS.
const MOONCAKE_RODAS = Rodas5P(autodiff = AutoFiniteDiff())
const ALGS_FD = Dict("Rodas5P" => Rodas5P(), "QNDF" => QNDF())
const ALGS_REV = Dict("Rodas5P(autodiff=AutoFiniteDiff)" => MOONCAKE_RODAS, "QNDF" => QNDF())

const TOLS = ((1.0e-6, 1.0e-5), (1.0e-8, 1.0e-7), (1.0e-10, 1.0e-9), (1.0e-12, 1.0e-11))
const JAN = ForwardDiff.jacobian(analytic_g, THETA_NOMINAL)
const COLSCALE = [maximum(abs.(JAN[:, j])) for j in 1:5]

# Column-scaled absolute error: |dJ_ij| / max_i |Jref_ij|. Exact-zero / 1e-218 true entries are
# judged by absolute deviation relative to their column's largest sensitivity (no relative floor).
col_scaled_err(J, Jref = JAN) = maximum(abs.(J .- Jref) ./ max.(COLSCALE', 1.0e-300))
# Gradient of a projection: error scaled by the largest reference component.
proj_err(g, ref) = maximum(abs.(g .- ref)) / maximum(abs.(ref))

@testset "Chunk 1a: stiff solve-through-AD compatibility probe" begin
    @testset "primal vs analytic, refinement sweep" begin
        for (name, alg) in ALGS_FD, (at, rt) in TOLS
            err = maximum(abs.(numeric_g(THETA_NOMINAL, alg; abstol = at, reltol = rt) .-
                               analytic_g(THETA_NOMINAL)))
            @info "primal" name at rt err
            @test err < 20 * max(at, rt)
        end
    end

    @testset "ForwardDiff(solve) vs analytic Jacobian, tolerance sweep" begin
        # Calibrated from the observed sweep (docs/CHUNK1A_RESULTS.md): err <= factor*reltol,
        # factor 50 (Rodas5P) / 500 (QNDF), plus monotone refinement and an absolute floor
        # ceiling at the tightest pair.
        factor = Dict("Rodas5P" => 50.0, "QNDF" => 500.0)
        final_ceiling = Dict("Rodas5P" => 1.0e-9, "QNDF" => 1.0e-7)
        for (name, alg) in ALGS_FD
            errs = Float64[]
            for (at, rt) in TOLS
                J = ForwardDiff.jacobian(th -> numeric_g(th, alg; abstol = at, reltol = rt), THETA_NOMINAL)
                e = col_scaled_err(J)
                push!(errs, e)
                @test e < factor[name] * rt
            end
            @info "ForwardDiff(solve) vs analytic, col-scaled error by tolerance" name errs
            @test issorted(errs; rev = true)
            @test errs[end] < final_ceiling[name]
        end
    end

    @testset "central FD vs analytic: gate for smooth Rodas5P, diagnostic for QNDF" begin
        tight = (1.0e-12, 1.0e-11)
        for (name, alg) in ALGS_FD
            g = th -> numeric_g(th, alg; abstol = tight[1], reltol = tight[2])
            errs = [col_scaled_err(fd_jacobian(g, THETA_NOMINAL, h)) for h in (1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5)]
            @info "central FD vs analytic (h=1e-2,1e-3,1e-4,1e-5), col-scaled" name errs
            if name == "Rodas5P"
                # truncation regime: error ~ h^2 (factor ~100 per decade), then noise at h<=1e-5
                @test errs[2] < errs[1] / 20
                @test errs[3] < errs[2] / 20
                @test errs[3] < 1.0e-8
            else
                @test all(isfinite, errs)   # QNDF: not a gate; see docs/CHUNK1A_RESULTS.md
            end
        end
    end

    @testset "ForwardDiff Jacobian: seeded non-cancelling projections vs analytic" begin
        rng = MersenneTwister(20260930)
        for name in ("Rodas5P", "QNDF")
            alg = ALGS_FD[name]
            J = ForwardDiff.jacobian(th -> numeric_g(th, alg; abstol = 1.0e-12, reltol = 1.0e-11), THETA_NOMINAL)
            for _ in 1:3
                w = randn(rng, 6)
                ref = JAN' * w
                @test maximum(abs.(ref)) > 1.0e-6                 # not a degenerate projection
                @test proj_err(J' * w, ref) < (name == "Rodas5P" ? 1.0e-9 : 1.0e-7)
            end
        end
    end

    @testset "Mooncake + continuous adjoint (MooncakeVJP) vs analytic, tolerance sweep" begin
        backend = AutoMooncake(; config = nothing)
        w = randn(MersenneTwister(20260930), 6)
        thetas = [THETA_NOMINAL, [1200.0, 0.8, 40.0, 1.3, 0.7], [900.0, 1.5, 60.0, 0.6, 1.4]]
        ceilings = Dict("Rodas5P(autodiff=AutoFiniteDiff)" => 1.0e-10, "QNDF" => 1.0e-8)
        for (name, alg) in ALGS_REV
            sa = GaussAdjoint(autojacvec = MooncakeVJP())
            sweep = Float64[]
            for (at, rt) in TOLS[2:end]
                obj = th -> dot(w, numeric_g(th, alg; sensealg = sa, abstol = at, reltol = rt))
                prep = prepare_gradient(obj, backend, THETA_NOMINAL)
                # one preparation reused across distinct theta, each vs the analytic gradient
                errs = [proj_err(gradient(obj, prep, backend, th), ForwardDiff.jacobian(analytic_g, th)' * w) for th in thetas]
                @test all(<(1.0e3 * rt), errs)
                push!(sweep, maximum(errs))
            end
            @info "Mooncake+GaussAdjoint(MooncakeVJP) max projected-gradient error vs analytic over 3 thetas, tol sweep (1e-8,1e-10,1e-12)" name sweep
            @test issorted(sweep; rev = true)
            @test sweep[end] < ceilings[name]
        end
        # QuadratureAdjoint, tightest pair
        for (name, alg) in ALGS_REV
            sa = QuadratureAdjoint(autojacvec = MooncakeVJP())
            obj = th -> dot(w, numeric_g(th, alg; sensealg = sa))
            prep = prepare_gradient(obj, backend, THETA_NOMINAL)
            errs = [proj_err(gradient(obj, prep, backend, th), ForwardDiff.jacobian(analytic_g, th)' * w) for th in thetas]
            @info "Mooncake+QuadratureAdjoint(MooncakeVJP), tightest tol" name errs
            @test maximum(errs) < ceilings[name]
        end
    end

    @testset "Prepared-cache reuse: two independent preparations, interleaved, changed theta" begin
        backend = AutoMooncake(; config = nothing)
        sa = GaussAdjoint(autojacvec = MooncakeVJP())
        rng = MersenneTwister(7)
        w1 = randn(rng, 6)
        w2 = randn(rng, 6)
        obj1 = th -> dot(w1, numeric_g(th, MOONCAKE_RODAS; sensealg = sa))
        obj2 = th -> dot(w2, numeric_g(th, MOONCAKE_RODAS; sensealg = sa))
        prep1 = prepare_gradient(obj1, backend, THETA_NOMINAL)
        prep2 = prepare_gradient(obj2, backend, THETA_NOMINAL)
        thetas = [THETA_NOMINAL, [1200.0, 0.8, 40.0, 1.3, 0.7], [900.0, 1.5, 60.0, 0.6, 1.4], THETA_NOMINAL]
        first_nominal = nothing
        for (k, th) in enumerate(thetas)
            g1 = gradient(obj1, prep1, backend, th)
            g2 = gradient(obj2, prep2, backend, th)
            g1b = gradient(obj1, prep1, backend, th)
            @test g1 == g1b                                      # bitwise repeat, no stale state
            @test proj_err(g1, ForwardDiff.jacobian(analytic_g, th)' * w1) < 1.0e-10
            @test proj_err(g2, ForwardDiff.jacobian(analytic_g, th)' * w2) < 1.0e-10
            k == 1 && (first_nominal = g1)
            k == 4 && @test g1 == first_nominal                  # returning to theta1 reproduces its gradient bitwise
        end
    end
end
