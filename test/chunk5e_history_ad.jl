# Chunk 5e: derivatives of the complete single pass (12-state ODE, sampled-HeI switch, 7-state ODE, Recfast tail, output assembly) with respect to the documented parameter set
# p = [F, A2s1s, hscale, nbscale]: Recfast fudge factor and H 2s-1s rate (tail), and multiplicative scales of H(z) and of the baryon density NH(z) (Omega_b) applied consistently in the
# pass and the tail. FROZEN (primal) discrete/preliminary parts: the node grid, the helium-switch node (no event-location derivative), the initial state and the 4b/4c preliminary
# history (Saha initialization, accessors, Yp, h100, T0, rate tables: not differentiated). REQUIRES COSMOREC_NATIVE_DATA_DIR.
# Routes: (i) ForwardDiff through `solve` (Rodas5P, analytic Jacobian/time gradient, primal error control) vs finite differences with step convergence, for the full pipeline incl. the
# spline assembly; (ii) prepared Mooncake VJP directly through the Rodas5P solve (plain-function RHS with global constant model, explicit Jacobian/time gradient) vs the ForwardDiff projection, on selected nodes.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
include("chunk5e_helpers.jl")

@testset "Chunk 5e: gradients of the complete single pass" begin
    @testset "frozen branch and node selection" begin
        @test K5E == 1342 && NZ5E == 3000
        @test issorted(SEL5E) && SEL5E[3] < K5E < SEL5E[4]            # 12-state nodes before, 7-state nodes after the switch
        # the subset-node pass reproduces the full-grid pass at the selected nodes (same solver; frozen switch)
        ps = recombination_pass(RM5, SOLVEP5E; k_switch = K5E, nodes = SEL5E)
        for (m, i) in enumerate(ps.node_index)
            i in SEL5E || continue
            @test isapprox(ps.Xe[m], PASS0_5E.Xe[i]; rtol = 5e-8)
            @test isapprox(ps.Te[m], PASS0_5E.Te[i]; rtol = 5e-8)
        end
        @test_throws ArgumentError recombination_pass(RM5, SOLVEP5E; nodes = SEL5E)
    end

    @testset "ForwardDiff through the complete pipeline vs finite differences" begin
        f0 = full_out(P0_5E)
        @test all(isfinite, f0) && length(f0) == 2 * length(ZG5E)
        J = ForwardDiff.jacobian(full_out, P0_5E)
        @test all(isfinite, J)
        errs = Float64[]
        for hrel in (1.0e-2, 1.0e-3, 1.0e-4)
            Jfd = fd_jac(full_out, P0_5E, hrel)
            push!(errs, elast_err(J, Jfd, f0, P0_5E))
        end
        rec5e("FD_errors_max", maximum(errs)); rec5e("FD_best", minimum(errs)); @info "ForwardDiff vs central FD, elasticity error, h = 1e-2, 1e-3, 1e-4" errs
        @test minimum(errs) < 1.0e-4                                    # finite differences are limited by the solver tolerance (observed best 1.2e-5, worst step 1.2e-4)
        # structure: F and A2s1s act only through the Recfast tail (Xe, Te at z >= 50 do not depend on them)
        nz = length(ZG5E); hi = [i for i in 1:nz if ZG5E[i] > 60]
        @test all(abs.(J[hi, 1]) .< 1e-9 * maximum(abs, J[:, 3])) && all(abs.(J[hi, 2]) .< 1e-9 * maximum(abs, J[:, 3]))
        # JVP along a random direction vs symmetric difference
        v = randn(Random.Xoshiro(20261020), 4) .* [0.02, 0.05, 0.01, 0.01]
        dd = ForwardDiff.derivative(t -> full_out(P0_5E .+ t .* v), 0.0)
        es = Float64[]
        for h in (1e-1, 1e-2, 1e-3)
            ref = (full_out(P0_5E .+ h .* v) .- full_out(P0_5E .- h .* v)) ./ (2h)
            push!(es, maximum(abs.(dd .- ref) ./ abs.(f0)))
        end
        rec5e("JVP_dir", minimum(es)); @info "JVP vs symmetric difference (relative output change per unit direction), h = 1e-1, 1e-2, 1e-3" es; @test minimum(es) < 1e-5
    end

    @testset "prepared Mooncake VJP (through the Rodas5P steps) vs ForwardDiff projection on selected nodes" begin
        f_fd = p -> sel_out(p, :fd)
        f_rv = p -> sel_out(p, :rev)
        nout = length(f_fd(P0_5E))
        rng = Random.Xoshiro(20261021)
        obj(p, w) = dot(w, f_rv(p))
        prepA = prepare_gradient(obj, MC5E, P0_5E, Constant(randn(rng, nout)))
        prepB = prepare_gradient(obj, MC5E, P0_5E, Constant(randn(rng, nout)))
        for kk in 1:2
            p = kk == 1 ? P0_5E : P0_5E .* (1 .+ 1.0e-3 .* randn(rng, 4))
            Jp = ForwardDiff.jacobian(f_fd, p)
            nsel = length(SEL5E)
            blocks = (("all", collect(1:nout)), ("Xe_pass", 1:(nsel + 1)), ("Te_pass", (nsel + 2):(2nsel + 2)), ("tail", (2nsel + 3):nout))
            for (bn, rows) in blocks
                w = zeros(nout); bn == "all" ? (w .= randn(rng, nout)) : (w[rows] .= 1.0)
                ref = Jp' * w
                gA = gradient(obj, prepA, MC5E, p, Constant(w)); gB = gradient(obj, prepB, MC5E, p, Constant(w))
                e = max(maximum(abs.(gA .- ref)), maximum(abs.(gB .- ref))) / maximum(abs, ref); rec5e("VJP:" * bn, e)
                @test e < 1.0e-6
                @test gA == gB
            end
        end
    end

    @testset "limitations are real: the helium-switch node is discrete" begin
        # perturbing the physics by a small amount can move the first node satisfying the criterion: the frozen k_switch is the primal one
        p2 = recombination_pass(RM5, SOLVEP5E; hscale = 1.0 + 1e-3)
        rec5e("k_switch_shift_hscale_1e-3", Float64(p2.k_switch - K5E)); @test abs(p2.k_switch - K5E) <= 5
    end

    @testset "info" begin
        foreach(kv -> println("Chunk5e-AD ", kv[1], " = ", kv[2]), sort!(collect(OBS5E)))
    end
end
