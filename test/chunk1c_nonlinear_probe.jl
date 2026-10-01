# Chunk 1c: nonlinear-state stiff ODE compatibility probe (two states, test-only).
#   x' = -a x^2, y' = -c y, theta = [a, c, x0, y0]; oracle x = x0/(1+a x0 s), y = y0 exp(-c s).
# Gates are calibrated from the observed sweep (docs/CHUNK1C_RESULTS.md); the solver-internal
# autodiff policy (AutoFiniteDiff) is distinct from the OUTER AD (Mooncake / ForwardDiff).
using Test
using LinearAlgebra: dot
using Random
using ForwardDiff
using SciMLBase: ODEFunction, ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqBDF: QNDF
using SciMLSensitivity: GaussAdjoint, QuadratureAdjoint, MooncakeVJP
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoFiniteDiff

include("chunk1c_helpers.jl")

const RODAS_FD1C = Rodas5P(autodiff = AutoFiniteDiff())
const QNDF_DEF1C = QNDF()
const TOLS1C = ((1.0e-6, 1.0e-5), (1.0e-8, 1.0e-7), (1.0e-10, 1.0e-9), (1.0e-12, 1.0e-11))
const THETAS1C = [THETA1C, [800.0, 1.5, 1.3, 0.7], [1500.0, 0.6, 0.8, 1.4]]
const JREF1C = ForwardDiff.jacobian(nl_exact, THETA1C)
colscaled1c(J, Jr = JREF1C) = maximum(abs.(J .- Jr) ./ maximum(abs.(Jr), dims = 1))
proj1c(g, ref) = maximum(abs.(g .- ref)) / maximum(abs.(ref))

# scalar objective of a seeded projection of the RHS; rhs! mutates a freshly allocated du (no global caches)
rhs_proj1c(w, z) = (du = similar(z, 2); nl_rhs!(du, view(z, 1:2), view(z, 3:4), 0.0); dot(w, du))

@testset "Chunk 1c: nonlinear quadratic-loss stiff ODE probe" begin
    @testset "RHS, analytic Jacobian, derivatives" begin
        u = [0.7, 0.3]; p = [1000.0, 1.0]
        du = zeros(2); nl_rhs!(du, u, p, 0.0)
        @test du ≈ [-1000.0 * 0.49, -0.3] rtol = 1.0e-15
        @test nl_rhs!(du, u, p, 123.0) === nothing && isapprox(du, [-490.0, -0.3]; rtol = 1.0e-15)      # autonomous: s-independent
        J = zeros(2, 2); nl_jac!(J, u, p, 0.0)
        @test J == [-2 * 1000.0 * 0.7 0.0; 0.0 -1.0]
        @test J ≈ ForwardDiff.jacobian(z -> (d = similar(z); nl_rhs!(d, z, p, 0.0); d), u) atol = 1.0e-12
        dp = ForwardDiff.jacobian(q -> (d = similar(q, 2); nl_rhs!(d, u, q, 0.0); d), p)
        @test dp ≈ [-0.49 0.0; 0.0 -0.3] rtol = 1.0e-15 atol = 1.0e-300                                        # d f/d(a,c)
        # the hook mutates its argument and does not depend on s; the time gradient of the ORIGINAL rhs is exactly zero
        dT = [9.0, 9.0]; nl_tgrad!(dT, u, p, 0.5); @test dT == [0.0, 0.0]
        @test ForwardDiff.derivative(s -> (d = zeros(typeof(s), 2); nl_rhs!(d, u, p, s); d), 0.3) == [0.0, 0.0]
        # Jacobian is state dependent (differs between two states)
        J2 = zeros(2, 2); nl_jac!(J2, [0.1, 0.3], p, 0.0); @test J2[1, 1] == -200.0 && J2[1, 1] != J[1, 1]
        # oracle: closed-form values, ODE satisfied by the oracle, active sensitivities
        a, c, x0, y0 = THETA1C
        @test nl_exact(THETA1C)[1:2] ≈ [x0 / (1 + a * x0 * 0.0005), y0 * exp(-c * 0.0005)]
        for s in (0.0005, 0.02, 0.5)
            @test ForwardDiff.derivative(t -> x0 / (1 + a * x0 * t), s) ≈ -a * (x0 / (1 + a * x0 * s))^2 atol = 1.0e-14
        end
        @test all(>(1.0e-5), maximum(abs.(JREF1C), dims = 1))                    # all four parameters resolvable
        @test JREF1C[1, 3] ≈ 1 / (1 + a * x0 * 0.0005)^2 atol = 1.0e-13          # dx/dx0 at s_early, hand-derived
        @test JREF1C[5, 1] ≈ -x0^2 * 0.5 / (1 + a * x0 * 0.5)^2 atol = 1.0e-15   # dx/da at s_late
    end

    @testset "Mooncake RHS projections vs ForwardDiff (one preparation, reused)" begin
        backend = AutoMooncake(; config = nothing)
        zs = ([0.7, 0.3, 1000.0, 1.0], [0.1, 0.9, 40.0, 3.0], [1.2, 0.05, 5.0, 0.1])
        for seed in 1:4
            w = randn(MersenneTwister(seed), 2)
            obj = z -> rhs_proj1c(w, z)
            prep = prepare_gradient(obj, backend, zs[1])
            for z in zs
                g = gradient(obj, prep, backend, z)
                ref = ForwardDiff.gradient(obj, z)
                @test maximum(abs.(g .- ref)) / maximum(abs.(ref)) < 1.0e-13
            end
        end
    end

    @testset "primal vs analytic (componentwise, atol+rtol|u| units)" begin
        ref = nl_exact(THETA1C)
        for (at, rt) in TOLS1C
            for hooks in (Val(:none), Val(:jac))
                sol = nl_observe(THETA1C, RODAS_FD1C; abstol = at, reltol = rt, hooks = hooks)
                @test maximum(abs.(sol .- ref) ./ (at .+ rt .* abs.(ref))) < 0.5       # observed <= 0.11
            end
            q = nl_observe(THETA1C, QNDF_DEF1C; abstol = at, reltol = rt)
            @test maximum(abs.(q .- ref) ./ (at .+ rt .* abs.(ref))) < 50              # QNDF observed <= 30: weaker error control
        end
    end

    @testset "ForwardDiff through the stiff solve vs analytic Jacobian (+ central FD secondary)" begin
        for (name, alg, factor) in (("Rodas5P(default autodiff)", Rodas5P(), 1.0), ("Rodas5P(AutoFiniteDiff)", RODAS_FD1C, 1.0), ("QNDF", QNDF_DEF1C, 50.0))
            errs = Float64[]
            for (at, rt) in TOLS1C
                J = ForwardDiff.jacobian(t -> nl_observe(t, alg; abstol = at, reltol = rt), THETA1C)
                push!(errs, colscaled1c(J))
                @test errs[end] < factor * rt                  # observed: Rodas <= 0.45 rt, QNDF <= 36 rt (1e-11)
            end
            @info "ForwardDiff(solve) vs analytic col-scaled error, tol 1e-5,1e-7,1e-9,1e-11" name errs
            @test issorted(errs; rev = true)
        end
        # analytic jac hook (original ODE only): Dual-compatible, same answer as the hook-free solve
        Jh = ForwardDiff.jacobian(t -> nl_observe(t, RODAS_FD1C; hooks = Val(:jac)), THETA1C)
        @test colscaled1c(Jh) < 1.0e-10
        # central differences of the numerical solve: truncation O(h^2) (x100 per decade) down to h=1e-4, noise floor after
        g = t -> nl_observe(t, RODAS_FD1C)
        fd(h) = (J = zeros(6, 4); for j in 1:4
            hh = h * max(abs(THETA1C[j]), 1.0); tp = copy(THETA1C); tp[j] += hh; tm = copy(THETA1C); tm[j] -= hh
            J[:, j] = (g(tp) .- g(tm)) ./ (2hh)
        end; J)
        errs = [colscaled1c(fd(h)) for h in (1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5)]
        @info "central FD vs analytic (h=1e-2..1e-5) col-scaled (Rodas5P, tol 1e-12/1e-11)" errs
        @test errs[2] < errs[1] / 50 && errs[3] < errs[2] / 50   # observed 1.85e-5, 1.85e-7, 1.9e-9
        @test errs[3] < 1.0e-8
    end

    @testset "Mooncake + GaussAdjoint(MooncakeVJP) through nonlinear stiff solve" begin
        backend = AutoMooncake(; config = nothing)
        w = randn(MersenneTwister(3), 6)
        sa = GaussAdjoint(autojacvec = MooncakeVJP())
        for (name, alg, factor) in (("Rodas5P(AutoFiniteDiff)", RODAS_FD1C, 0.5), ("QNDF", QNDF_DEF1C, 10.0))
            sweep = Float64[]
            for (at, rt) in TOLS1C
                obj = th -> dot(w, nl_observe(th, alg; sensealg = sa, abstol = at, reltol = rt))
                prep = prepare_gradient(obj, backend, THETA1C)    # one preparation per objective and tolerance
                errs = Float64[]
                for th in THETAS1C
                    Jt = ForwardDiff.jacobian(nl_exact, th)
                    push!(errs, proj1c(gradient(obj, prep, backend, th), Jt' * w))
                    # ForwardDiff through the numerical solve at the same theta
                    Jfd = ForwardDiff.jacobian(t -> nl_observe(t, alg; abstol = at, reltol = rt), th)
                    @test proj1c(Jfd' * w, Jt' * w) < max(10 * rt, 1.0e-9)
                end
                @test maximum(errs) < factor * rt             # observed: Rodas <= 0.23 rt; QNDF <= 4.6 rt
                push!(sweep, maximum(errs))
            end
            @info "Mooncake+Gauss max projected-gradient error over 3 thetas, tol 1e-5,1e-7,1e-9,1e-11" name sweep
            @test issorted(sweep; rev = true)
        end
        # analytic jac hook on the ORIGINAL ode (it does not reach the backward adjoint solve)
        obj = th -> dot(w, nl_observe(th, RODAS_FD1C; sensealg = sa, hooks = Val(:jac)))
        prep = prepare_gradient(obj, backend, THETA1C)
        @test proj1c(gradient(obj, prep, backend, THETA1C), JREF1C' * w) < 1.0e-11
        # QuadratureAdjoint secondary: two tolerances
        for (at, rt) in ((1.0e-10, 1.0e-9), (1.0e-12, 1.0e-11))
            objq = th -> dot(w, nl_observe(th, RODAS_FD1C; sensealg = QuadratureAdjoint(autojacvec = MooncakeVJP()), abstol = at, reltol = rt))
            prepq = prepare_gradient(objq, backend, THETA1C)
            @test maximum(proj1c(gradient(objq, prepq, backend, th), ForwardDiff.jacobian(nl_exact, th)' * w) for th in THETAS1C) < 0.5 * rt
        end
        # default Rodas5P autodiff fails in the backward adjoint solve (retained blocker, as in Chunk 1a)
        msg = try
            objd = th -> dot(w, nl_observe(th, Rodas5P(); sensealg = sa))
            gradient(objd, prepare_gradient(objd, backend, THETA1C), backend, THETA1C)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("time gradient", msg)
    end

    @testset "prepared-cache reuse: two independent preparations, interleaved, changed theta" begin
        backend = AutoMooncake(; config = nothing)
        sa = GaussAdjoint(autojacvec = MooncakeVJP())
        rng = MersenneTwister(7)
        w1 = randn(rng, 6); w2 = randn(rng, 6)
        obj1 = th -> dot(w1, nl_observe(th, RODAS_FD1C; sensealg = sa))
        obj2 = th -> dot(w2, nl_observe(th, RODAS_FD1C; sensealg = sa))
        prep1 = prepare_gradient(obj1, backend, THETA1C)
        prep2 = prepare_gradient(obj2, backend, THETA1C)
        seq = [THETAS1C[1], THETAS1C[2], THETAS1C[3], THETAS1C[1]]
        first_g = nothing
        for (k, th) in enumerate(seq)
            g1 = gradient(obj1, prep1, backend, th)
            g2 = gradient(obj2, prep2, backend, th)
            g1b = gradient(obj1, prep1, backend, th)
            @test g1 == g1b
            Jt = ForwardDiff.jacobian(nl_exact, th)
            @test proj1c(g1, Jt' * w1) < 1.0e-10
            @test proj1c(g2, Jt' * w2) < 1.0e-10
            k == 1 && (first_g = g1)
            k == 4 && @test g1 == first_g
        end
    end
end
