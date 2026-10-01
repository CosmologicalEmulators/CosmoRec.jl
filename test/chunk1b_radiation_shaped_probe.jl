# Chunk 1b: method-of-lines / structured-sparse-linear-solve probe on a small radiation-SHAPED problem.
# Toy only (linear in u): not the CosmoRec PDE, no native physics, no native fixtures. The manufactured solution
# (see chunk1b_helpers.jl) is the toolchain oracle; native-original text fixtures remain mandatory before
# any physics stage. Results/definitions: docs/CHUNK1B_RESULTS.md.
using Test
using LinearAlgebra
using SparseArrays
using Random
using ForwardDiff
using SciMLBase: ODEFunction, ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqBDF: QNDF
using LinearSolve, Sparspak
using SciMLSensitivity: GaussAdjoint, QuadratureAdjoint, MooncakeVJP
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoFiniteDiff

include("chunk1b_helpers.jl")

# Solver-internal autodiff must be AutoFiniteDiff under Mooncake (default fails: FunctionWrappersWrapper error on the
# state-carrying RHS; see benchmark/chunk1b_candidate_probe.jl). Mooncake stays the outer AD; MooncakeVJP the RHS-VJP.
# CAVEAT: the problem is linear in u; nonlinear state-dependent RHS is a separate future gate.
const RODAS_SP = Rodas5P(autodiff = AutoFiniteDiff(), linsolve = SparspakFactorization())
const QNDF_SP = QNDF(autodiff = AutoFiniteDiff(), linsolve = SparspakFactorization())
const MESHES = (("mesh_a", MeshOps(mesh_a())), ("mesh_b", MeshOps(mesh_b())))
const THETAS1B = [THETA1B, [1.5, 3.5, 0.8, 1.2, 0.4, 0.5], [2.5, 2.0, 1.3, 0.7, 0.8, 0.3]]
const TOLS1B = ((1.0e-6, 1.0e-5), (1.0e-8, 1.0e-7), (1.0e-10, 1.0e-9), (1.0e-12, 1.0e-11))

# column-scaled absolute Jacobian error and its worst (row, col)
function worst_col_scaled(J, Jref)
    cs = [maximum(abs.(Jref[:, j])) for j in axes(Jref, 2)]
    D = abs.(J .- Jref) ./ cs'
    idx = argmax(D)
    return D[idx], Tuple(idx)
end
proj_err1b(g, ref) = maximum(abs.(g .- ref)) / maximum(abs.(ref))

@testset "Chunk 1b: radiation-shaped MOL structured-solve probe" begin
    @testset "A. numerical primitives" begin
        @testset "Fornberg weights vs hand-checkable uniform 5-point stencils" begin
            h = 0.1
            xs = collect(-2h:h:2h)
            W = fornberg(0.0, xs, 2)
            @test W[:, 2] ≈ [1, -8, 0, 8, -1] ./ (12h) atol = 1.0e-10
            @test W[:, 3] ≈ [-1, 16, -30, 16, -1] ./ (12h^2) atol = 1.0e-9
            @test sum(W[:, 1]) ≈ 1                                   # interpolation weights sum to 1
        end
        for (mname, ops) in MESHES
            @testset "$mname" begin
                x = ops.x; N = ops.N; xi = x[2:(end - 1)]
                @test issorted(x) && x[1] == 0 && x[end] == 1 && length(unique(diff(x))) > 3   # genuinely nonuniform
                @test eltype(ops.A1) == Float64 && eltype(ops.A2) == Float64
                # stencil topology: 5 nodes per row, centered inside, one-sided at the two ends
                for i in 1:(N - 1)
                    want = collect(stencil_nodes(i, N))
                    @test length(want) == 5
                    got = findall(!iszero, ops.D2[i + 1, :]) .- 1
                    @test got == want
                end
                @test collect(stencil_nodes(1, N)) == 0:4
                @test collect(stencil_nodes(N - 1, N)) == (N - 4):N
                @test collect(stencil_nodes(3, N)) == 1:5
                # interior-block band structure: lower/upper half-bandwidth <= 3, pattern covers both blocks and I
                r, c, _ = findnz(ops.pattern)
                @test maximum(abs.(r .- c)) == 3
                @test all(i -> ops.pattern[i, i] != 0, 1:(N - 1))
                @test all(ops.pattern[i, j] != 0 for (i, j) in zip(findnz(ops.A1)[1], findnz(ops.A1)[2]))
                @test all(ops.pattern[i, j] != 0 for (i, j) in zip(findnz(ops.A2)[1], findnz(ops.A2)[2]))
                # sign/diagonal structure: diffusion block has negative diagonal
                @test all(<(0), diag(ops.A2))
                # polynomial differentiation exact to rounding for degree <= 4, and NOT exact for degree 5
                for d in 0:4
                    d1 = d == 0 ? zeros(N - 1) : d .* xi .^ (d - 1)
                    d2 = d < 2 ? zeros(N - 1) : d * (d - 1) .* xi .^ (d - 2)
                    @test maximum(abs.(ops.D1[2:N, :] * x .^ d .- d1)) < 1.0e-11
                    @test maximum(abs.(ops.D2[2:N, :] * x .^ d .- d2)) < 1.0e-10
                end
                @test maximum(abs.(ops.D2[2:N, :] * x .^ 5 .- 20 .* xi .^ 3)) > 1.0e-4
                # boundary lifting: interior-block + lifted boundary columns reproduce full-row derivative of a cubic
                p(x) = 1 + 0.5x + 2x^2 - x^3
                u = p.(xi)
                @test ops.A1 * u .+ ops.L1L .* p(0.0) .+ ops.L1R .* p(1.0) ≈ (0.5 .+ 4 .* xi .- 3 .* xi .^ 2) atol = 1.0e-10
                @test ops.A2 * u .+ ops.L2L .* p(0.0) .+ ops.L2R .* p(1.0) ≈ (4 .- 6 .* xi) atol = 1.0e-9
                # boundary node enters only the first two / last two interior rows (centered stencils start at node 1 from row 3 on)
                for L in (ops.L1L, ops.L2L); @test findall(!iszero, L) == 1:2; end
                for L in (ops.L1R, ops.L2R); @test findall(!iszero, L) == (N - 2):(N - 1); end
            end
        end
    end

    @testset "B. coefficient/boundary/initial/source/operator maps: ForwardDiff vs analytic" begin
        for th in THETAS1B, s in (0.0, 0.3, 1.0)
            r = rate_r(th); dr = [0.3, 0.2, 0.1, 0.0, 0.0, 1.0]
            @test ForwardDiff.gradient(t -> kappa(s, t), th) ≈ [1 + 0.1 * sin(s), 0, 0, 0, 0, 0] atol = 1.0e-14
            @test ForwardDiff.gradient(t -> vel(s, t), th) ≈ [0, 1 + 0.5 * s, 0, 0, 0, 0] atol = 1.0e-14
            @test ForwardDiff.gradient(t -> gam(s, t), th) ≈ [0, 0, 1 + 0.2 * cos(s), 0, 0, 0] atol = 1.0e-14
            @test all(>(0), (kappa(s, th), vel(s, th), gam(s, th)))
            gL = ForwardDiff.gradient(t -> bc_left(s, t), th)
            @test gL ≈ th[4] * exp(-r * s) .* (-s .* dr) .+ [0, 0, 0, exp(-r * s), 0, 0] atol = 1.0e-13
            gR = ForwardDiff.gradient(t -> bc_right(s, t), th)
            @test gR ≈ 3 * th[4] * exp(-r * s) .* (-s .* dr) .+ [0, 0, 0, 3 * exp(-r * s), exp(-s), 0] atol = 1.0e-13
            for x in (0.0, 0.37, 1.0)
                @test ForwardDiff.gradient(t -> init_u(x, t), th) ≈ [0, 0, 0, 1 + x + x^2, x^3, 0] atol = 1.0e-14
                # source: ForwardDiff vs central difference with per-parameter step (truncation ~h^2 ~1e-10)
                gs = ForwardDiff.gradient(t -> src(s, x, t), th)
                fdg = [(src(s, x, replace_at(th, j, th[j] + 1.0e-5)) - src(s, x, replace_at(th, j, th[j] - 1.0e-5))) / 2.0e-5 for j in 1:6]
                @test maximum(abs.(gs .- fdg)) < 1.0e-7 * max(1.0, maximum(abs.(gs)))
                # manufactured identity: u_s = kappa u_xx + v u_x - gamma u + src (continuum, hand-written derivatives)
                @test us_exact(s, x, th) ≈ kappa(s, th) * uxx_exact(s, x, th) + vel(s, th) * ux_exact(s, x, th) -
                      gam(s, th) * u_exact(s, x, th) + src(s, x, th) atol = 1.0e-12
                @test ForwardDiff.derivative(ss -> u_exact(ss, x, th), s) ≈ us_exact(s, x, th) atol = 1.0e-13
            end
        end
        # assembled operator: value and derivatives (u, boundary values, k, v, g) vs hand-derived continuum
        for (mname, ops) in MESHES
            xi = ops.x[2:(end - 1)]
            p(x) = 1 + 0.5x + 2x^2 - x^3; dp(x) = 0.5 + 4x - 3x^2; d2p(x) = 4 - 6x
            q(x) = (1 - x)^3; dq(x) = -3 * (1 - x)^2; d2q(x) = 6 * (1 - x)
            k, v, g = 1.7, 2.3, 0.6
            val = apply_ops(ops, p.(xi), p(0.0), p(1.0), k, v, g)
            @test val ≈ k .* d2p.(xi) .+ v .* dp.(xi) .- g .* p.(xi) atol = 1.0e-9
            # total derivative along u = p + d q with boundary values following q (exercises lifting with Duals)
            dd = ForwardDiff.derivative(d -> apply_ops(ops, p.(xi) .+ d .* q.(xi), p(0.0) + d * q(0.0), p(1.0) + d * q(1.0), k, v, g), 0.0)
            @test dd ≈ k .* d2q.(xi) .+ v .* dq.(xi) .- g .* q.(xi) atol = 1.0e-9
            for (j, expect) in ((1, d2p.(xi)), (2, dp.(xi)), (3, -p.(xi)))
                kk = [k, v, g]
                dj = ForwardDiff.derivative(a -> apply_ops(ops, p.(xi), p(0.0), p(1.0), (j == 1 ? a : kk[1]), (j == 2 ? a : kk[2]), (j == 3 ? a : kk[3])), kk[j])
                @test dj ≈ expect atol = 1.0e-9
            end
            # full RHS at the manufactured state equals the continuum u_s at every interior node
            for th in THETAS1B, s in (0.0, 0.4, 1.0)
                @test rhs_interior(ops, u_exact.(s, xi, Ref(th)), s, th) ≈ us_exact.(s, xi, Ref(th)) atol = 1.0e-9
            end
            # analytic Jacobian equals the matrix of the RHS (u-derivative) for the linear problem
            th = THETA1B; s = 0.3
            Jm = copy(convert(SparseMatrixCSC{Float64, Int}, ops.pattern))
            RadJac(ops)(Jm, nothing, th, s)
            Jad = ForwardDiff.jacobian(u -> rhs_interior(ops, u, s, th), u0_interior(ops, th))
            @test Matrix(Jm) ≈ Jad atol = 1.0e-9
        end
    end

    @testset "C. primal public solve vs manufactured solution" begin
        for (mname, ops) in MESHES, (an, alg, bound) in (("Rodas5P", RODAS_SP, 1.0), ("QNDF", QNDF_SP, 5.0))
            ref = exact_observe(THETA1B, ops)
            for (at, rt) in TOLS1B
                sol = solve_rad(THETA1B, ops, alg; abstol = at, reltol = rt)
                @test string(sol.retcode) == "Success"
                @test length(sol.u) == 3 && sol.t == S1B_SAVE
                err = abs.(vcat(sol.u...) .- ref) ./ (at .+ rt .* abs.(ref))   # componentwise, in units of atol+rtol|u|
                # bounds calibrated from sweep: observed max Rodas5P 0.12, QNDF 1.9 (benchmark/chunk1b_sweep_diagnostic.jl)
                @test maximum(err) < bound
            end
        end
        # structured/default/dense-diagnostic routes agree on the primal
        ops = MESHES[1][2]; ref = exact_observe(THETA1B, ops)
        for (label, kw) in (("analytic jac + Sparspak", (; alg = RODAS_SP, jac = :analytic)),
                            ("sparse prototype only + Sparspak", (; alg = RODAS_SP, jac = :proto)),
                            ("analytic sparse jac + default linsolve", (; alg = Rodas5P(autodiff = AutoFiniteDiff()), jac = :analytic)),
                            ("DENSE diagnostic (no prototype)", (; alg = Rodas5P(autodiff = AutoFiniteDiff()), jac = :none)))
            num = observe(THETA1B, ops, kw.alg; jac = kw.jac)
            @test maximum(abs.(num .- ref) ./ (1.0e-10 .+ 1.0e-9 .* abs.(ref))) < 1.0
        end
        # early / interpolated / late are all non-trivial and distinct
        @test length(unique(round.(exact_observe(THETA1B, ops); digits = 6))) == length(ops.x) * 3 - 6
    end

    @testset "C2. structured linear-solve route evidence (matrix type and factorization)" begin
        using SciMLBase: init, step!
        ops = MESHES[1][2]
        for T in (Float64, ForwardDiff.Dual{Nothing, Float64, 2})
            th = T === Float64 ? THETA1B : ForwardDiff.Dual{Nothing}.(THETA1B, Ref(ForwardDiff.Partials((1.0, 0.0))))
            f = ODEFunction(RadRHS(ops); jac = RadJac(ops), jac_prototype = convert(SparseMatrixCSC{eltype(th), Int}, ops.pattern))
            prob = ODEProblem(f, u0_interior(ops, th), (0.0, 1.0), th)
            integ = init(prob, RODAS_SP; abstol = 1.0e-10, reltol = 1.0e-9)
            step!(integ)
            ls = integ.cache.linsolve
            @test f.jac_prototype isa SparseMatrixCSC{eltype(th)}
            @test integ.cache.W isa SparseMatrixCSC{eltype(th)}
            @test ls.alg isa SparspakFactorization
            @test ls.A isa SparseMatrixCSC{eltype(th)}
            @test nnz(ls.A) < length(ls.A) / 1.4                       # genuinely sparse (31 of 49 on mesh_a)
            @test occursin("SparseSolver", string(typeof(ls.cacheval)))   # Sparspak sparse-LU object, not a dense LU
        end
    end

    @testset "D. ForwardDiff Jacobian of the complete sampled solve vs analytic" begin
        for (mname, ops) in MESHES
            Jref = ForwardDiff.jacobian(t -> exact_observe(t, ops), THETA1B)
            cs = [maximum(abs.(Jref[:, j])) for j in 1:6]
            @test all(>(1.0e-2), cs)                                 # every active parameter has a non-trivial sensitivity
            for (an, alg, factor) in (("Rodas5P", RODAS_SP, 10.0), ("QNDF", QNDF_SP, 100.0))
                errs = Float64[]; worst = Tuple[]
                for (at, rt) in TOLS1B
                    J = ForwardDiff.jacobian(t -> observe(t, ops, alg; abstol = at, reltol = rt), THETA1B)
                    e, w = worst_col_scaled(J, Jref)
                    push!(errs, e); push!(worst, w)
                    @test e < factor * rt                            # calibrated: observed <= 1.0*rt (Rodas5P), 20*rt (QNDF)
                end
                @info "ForwardDiff(solve) vs analytic col-scaled error by tol (worst row,col)" mname an errs worst
                @test issorted(errs; rev = true)
            end
            # directional derivative along a seeded direction vs analytic
            d = randn(MersenneTwister(11), 6)
            dd = ForwardDiff.derivative(a -> observe(THETA1B .+ a .* d, ops, RODAS_SP; abstol = 1.0e-12, reltol = 1.0e-11), 0.0)
            @test maximum(abs.(dd .- Jref * d)) / maximum(abs.(Jref * d)) < 1.0e-9
        end
    end

    @testset "E. Mooncake on reusable maps vs ForwardDiff-transpose products" begin
        backend = AutoMooncake(; config = nothing)
        sgrid = (0.0, 0.3, 1.0)
        for (mname, ops) in MESHES
            xi = ops.x[2:(end - 1)]
            coeff_map(th) = vcat([[kappa(s, th), vel(s, th), gam(s, th), bc_left(s, th), bc_right(s, th)] for s in sgrid]...,
                u0_interior(ops, th), [src(0.3, x, th) for x in xi])
            op_map(th) = rhs_interior(ops, u0_interior(ops, th), 0.3, th)
            for (label, fmap) in (("coefficient/boundary/initial/source map", coeff_map), ("assembled operator RHS map", op_map))
                n = length(fmap(THETA1B))
                for seed in (1, 2, 3)
                    w = randn(MersenneTwister(seed), n)
                    obj = th -> dot(w, fmap(th))
                    prep = prepare_gradient(obj, backend, THETA1B)
                    for th in THETAS1B     # one preparation reused across changed same-shape theta
                        g = gradient(obj, prep, backend, th)
                        ref = ForwardDiff.jacobian(fmap, th)' * w
                        @test proj_err1b(g, ref) < 1.0e-11
                    end
                end
            end
        end
    end

    @testset "F. Mooncake + continuous adjoint (MooncakeVJP) through the structured solve vs analytic" begin
        backend = AutoMooncake(; config = nothing)
        for (mname, ops) in MESHES
            nout = 3 * (ops.N - 1)
            w = randn(MersenneTwister(20260930), nout)
            sa = GaussAdjoint(autojacvec = MooncakeVJP())
            sweep = Float64[]
            for (at, rt) in TOLS1B[2:end]
                obj = th -> dot(w, observe(th, ops, RODAS_SP; sensealg = sa, abstol = at, reltol = rt))
                prep = prepare_gradient(obj, backend, THETA1B)
                errs = [proj_err1b(gradient(obj, prep, backend, th), ForwardDiff.jacobian(t -> exact_observe(t, ops), th)' * w) for th in THETAS1B]
                @test maximum(errs) < 0.5 * rt                         # calibrated: observed <= 0.05*rt
                push!(sweep, maximum(errs))
            end
            @info "Mooncake+GaussAdjoint(MooncakeVJP), Rodas5P(AutoFiniteDiff)+Sparspak: max proj error over 3 thetas, tol 1e-8,1e-10,1e-12" mname sweep
            @test issorted(sweep; rev = true)
            @test sweep[end] < 1.0e-11
        end
        # QuadratureAdjoint, tight tolerance, both meshes
        for (mname, ops) in MESHES
            w = randn(MersenneTwister(20260930), 3 * (ops.N - 1))
            sa = QuadratureAdjoint(autojacvec = MooncakeVJP())
            obj = th -> dot(w, observe(th, ops, RODAS_SP; sensealg = sa))
            prep = prepare_gradient(obj, backend, THETA1B)
            errs = [proj_err1b(gradient(obj, prep, backend, th), ForwardDiff.jacobian(t -> exact_observe(t, ops), th)' * w) for th in THETAS1B]
            @info "Mooncake+QuadratureAdjoint, default tolerances (abstol 1e-10, reltol 1e-9)" mname errs
            @test maximum(errs) < 1.0e-10    # calibrated: observed 2.6e-11 / 2.8e-11
        end
        # QNDF reverse: LIMITED candidate. Error decreases under refinement but stays orders above Rodas5P; no accuracy claim.
        let ops = MESHES[1][2]
            w = randn(MersenneTwister(20260930), 3 * (ops.N - 1))
            sa = GaussAdjoint(autojacvec = MooncakeVJP())
            sweep = Float64[]
            for (at, rt) in TOLS1B[2:end]
                obj = th -> dot(w, observe(th, ops, QNDF_SP; sensealg = sa, abstol = at, reltol = rt))
                prep = prepare_gradient(obj, backend, THETA1B)
                push!(sweep, proj_err1b(gradient(obj, prep, backend, THETA1B), ForwardDiff.jacobian(t -> exact_observe(t, ops), THETA1B)' * w))
            end
            @info "QNDF(AutoFiniteDiff)+Sparspak reverse (LIMITED; not an accepted accuracy path)" sweep
            @test issorted(sweep; rev = true)
        end
    end

    @testset "G. prepared-cache reuse: independent preparations, interleaved, changed theta" begin
        backend = AutoMooncake(; config = nothing)
        ops = MESHES[1][2]
        sa = GaussAdjoint(autojacvec = MooncakeVJP())
        rng = MersenneTwister(7)
        w1 = randn(rng, 3 * (ops.N - 1)); w2 = randn(rng, 3 * (ops.N - 1))
        obj1 = th -> dot(w1, observe(th, ops, RODAS_SP; sensealg = sa))
        obj2 = th -> dot(w2, observe(th, ops, RODAS_SP; sensealg = sa))
        prep1 = prepare_gradient(obj1, backend, THETA1B)
        prep2 = prepare_gradient(obj2, backend, THETA1B)
        seq = [THETAS1B[1], THETAS1B[2], THETAS1B[3], THETAS1B[1]]
        first_g = nothing
        for (k, th) in enumerate(seq)
            g1 = gradient(obj1, prep1, backend, th)
            g2 = gradient(obj2, prep2, backend, th)
            g1b = gradient(obj1, prep1, backend, th)
            @test g1 == g1b
            Jr = ForwardDiff.jacobian(t -> exact_observe(t, ops), th)
            @test proj_err1b(g1, Jr' * w1) < 1.0e-10
            @test proj_err1b(g2, Jr' * w2) < 1.0e-10
            k == 1 && (first_g = g1)
            k == 4 && @test g1 == first_g
        end
    end
end
