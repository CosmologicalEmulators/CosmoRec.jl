# Chunk 5g: the PRIVATE primal RHS workspace (CosmoRec._rhs_workspace / CosmoRec._recombination_rhs_ws!, docs/ODE_CORE_WORKSPACE_DESIGN.md increment 1)
# vs the unchanged public allocating `recombination_rhs!` (the oracle): bitwise values (all 5a native rows, diffusion feedback, mixed/nested Dual,
# BigFloat), element types of the rate buffers vs the public rate functions, counted fallback, ForwardDiff through the private path, prepared
# DI-Mooncake with ONE workspace reused by two RHS evaluations inside one taped function, independent workspaces, allocations.
# REQUIRES COSMOREC_NATIVE_DATA_DIR (see test/chunk5_helpers.jl); fails loudly if unset.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(RM5) || include("chunk5_helpers.jl")
@isdefined(SETUP7) || include("chunk7_helpers.jl")

const MC5G = AutoMooncake(; config = nothing)
const ROWS5G = read_fixture5a(joinpath(@__DIR__, "fixtures", "native_ode_rhs.txt"))
# private RHS with a workspace built from the actual argument types (as the per-solve callback does); returns (f, used)
function rhs5g(z, y, rm; flag_He, hscale = 1.0, nbscale = 1.0, ws = CosmoRec._rhs_workspace(y, z, rm; hscale = hscale, nbscale = nbscale))
    T = promote_type(eltype(y), typeof(z), typeof(hscale), typeof(nbscale), CosmoRec.feedback_eltype(rm.diffusion), Float64)
    f = zeros(T, length(y)); used = CosmoRec._recombination_rhs_ws!(f, z, y, rm, ws; flag_He = flag_He, hscale = hscale, nbscale = nbscale)
    return f, used
end
pub5g(z, y, rm; flag_He, hscale = 1.0, nbscale = 1.0) = recombination_rhs(z, y, rm; flag_He = flag_He, hscale = hscale, nbscale = nbscale)

@testset "Chunk 5g: private primal RHS workspace vs the public allocating RHS" begin
    @testset "bitwise at every 5a native row (12/7-state, all patterns incl. the off-table DP fallback)" begin
        for r in ROWS5G
            f, used = rhs5g(r.z, r.y, RM5; flag_He = r.flag)
            @test used && f == pub5g(r.z, r.y, RM5; flag_He = r.flag)
        end
    end
    @testset "diffusion feedback model (pass-1 feedback)" begin
        S = (a...) -> solve_phase5(a...; reltol = 1.0e-12, a1 = 1.0e-18, aex = 1.0e-14, alg = Rodas5P(), tstops = false)
        p0 = recombination_pass(RM5, S); (_, fb) = hi_diffusion_stage(RM5, HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6), pass_solution_rows(p0))
        RMF = with_diffusion(RM5, fb)
        for k in (1, 500, 1342, 2000, 2900)
            y = p0.states[k]; z = p0.z[k]; flag = length(y) == 12
            f, used = rhs5g(z, y, RMF; flag_He = flag)
            @test used && f == pub5g(z, y, RMF; flag_He = flag)
        end
    end
    r12 = first(r for r in ROWS5G if r.flag && r.z == 2500.0 && r.pat == 0); r7 = first(r for r in ROWS5G if !r.flag && r.z == 400.0 && r.pat == 0)
    r12b = first(r for r in ROWS5G if r.flag && r.z == 1500.0 && r.pat == 0)                   # below the detailed-balance switch (z < 2000)
    @testset "mixed / nested element types (bitwise vs the public path) and rate-buffer element types" begin
        for r in (r12, r12b, r7)
            D(v, s) = ForwardDiff.Dual{typeof(ForwardDiff.Tag(rhs5g, Float64))}(v, s)
            yd = [D(v, i == 1 ? 1.0 : 0.0) for (i, v) in enumerate(r.y)]; hd = D(1.0, 1.0)
            for (y, z, hs, nb) in ((yd, r.z, hd, one(hd)), (r.y, r.z, D(1.0, 1.0), D(1.0, 0.0)), (r.y, D(r.z, 1.0), 1.0, 1.0), (yd, r.z, 1.0, 1.0))
                f, used = rhs5g(z, y, RM5; flag_He = r.flag, hscale = hs, nbscale = nb)
                @test used && f == pub5g(z, y, RM5; flag_He = r.flag, hscale = hs, nbscale = nb)
            end
            fb_, ub = setprecision(() -> rhs5g(r.z, BigFloat.(r.y), RM5; flag_He = r.flag), BigFloat, 256)
            @test ub && fb_ == setprecision(() -> pub5g(r.z, BigFloat.(r.y), RM5; flag_He = r.flag), BigFloat, 256)
            # buffer element types equal the public rate functions' output types for this (Dual state) workspace, in the branch this state takes
            ws = CosmoRec._rhs_workspace(yd, r.z, RM5); bg = CosmoRec.ode_background(RM5, r.z); Tg = bg.Tg
            rr = CosmoRec.get_rates(RM5.eff.htable, Tg, yd[1] * Tg)
            db = Tg / CosmoRec.T0_detailed_balance - 1.0 > CosmoRec.z_detailed_balance
            @test eltype(rr.A) === (db ? eltype(ws.HAdb) : eltype(ws.HA)) && eltype(rr.B) === eltype(ws.HB) && eltype(rr.R) === eltype(ws.HR)
            r.flag && @test eltype(CosmoRec.get_helium_rates(RM5.eff.hetable, Tg).A) === eltype(ws.HeA)
        end
    end
    @testset "type-mismatched call uses the public path (counted by the caller), never a conversion" begin
        ws = CosmoRec._rhs_workspace(r12.y, r12.z, RM5)
        yd = ForwardDiff.Dual{:t5g}.(r12.y, 1.0); f = zeros(eltype(yd), 12)
        @test !CosmoRec._recombination_rhs_ws!(f, r12.z, yd, RM5, ws; flag_He = true) && f == pub5g(r12.z, yd, RM5; flag_He = true)
    end
    @testset "ForwardDiff through the private path (workspace built from the Dual types) equals the public path" begin
        for r in (r12, r12b, r7)
            fw(y) = rhs5g(r.z, y, RM5; flag_He = r.flag)[1]; fp(y) = pub5g(r.z, y, RM5; flag_He = r.flag)
            @test ForwardDiff.jacobian(fw, r.y) == ForwardDiff.jacobian(fp, r.y)
            @test ForwardDiff.derivative(z -> rhs5g(z, r.y, RM5; flag_He = r.flag)[1], r.z) == ForwardDiff.derivative(z -> pub5g(z, r.y, RM5; flag_He = r.flag), r.z)
            @test ForwardDiff.derivative(h -> rhs5g(r.z, r.y, RM5; flag_He = r.flag, hscale = h, nbscale = one(h))[1], 1.0) ==
                  ForwardDiff.derivative(h -> pub5g(r.z, r.y, RM5; flag_He = r.flag, hscale = h, nbscale = one(h)), 1.0)
        end
    end
    @testset "prepared DI-Mooncake: one workspace reused by two RHS evaluations in one taped function" begin
        for r in (r7, r12)
            n = length(r.y); rng = MersenneTwister(9); w1 = randn(rng, n); w2 = randn(rng, n); y2s = 1 .+ 1.0e-6 .* (1:n)
            ws = CosmoRec._rhs_workspace(r.y, r.z, RM5)                 # Float64 workspace: the second call overwrites the first call's buffers
            fws(y) = dot(w1, rhs5g(r.z, y, RM5; flag_He = r.flag, ws = ws)[1]) + dot(w2, rhs5g(r.z, y .* y2s, RM5; flag_He = r.flag, ws = ws)[1])
            fpb(y) = dot(w1, pub5g(r.z, y, RM5; flag_He = r.flag)) + dot(w2, pub5g(r.z, y .* y2s, RM5; flag_He = r.flag))
            @test fws(r.y) == fpb(r.y)
            ref = ForwardDiff.gradient(fpb, r.y)
            prep = prepare_gradient(fws, MC5G, r.y); g1 = gradient(fws, prep, MC5G, r.y); g2 = gradient(fws, prep, MC5G, r.y)
            @test maximum(abs.(g1 .- ref)) / max(maximum(abs, ref), floatmin()) < 1e-12 && g1 == g2
        end
    end
    @testset "independent workspaces (alternating states)" begin
        wa = CosmoRec._rhs_workspace(r12.y, r12.z, RM5); wb = CosmoRec._rhs_workspace(r12.y, r12.z, RM5); y2 = r12.y .* 1.001
        A = rhs5g(r12.z, r12.y, RM5; flag_He = true, ws = wa)[1]; B = rhs5g(r12.z - 5.0, y2, RM5; flag_He = true, ws = wb)[1]
        C = rhs5g(r12.z - 5.0, y2, RM5; flag_He = true, ws = wa)[1]; D = rhs5g(r12.z, r12.y, RM5; flag_He = true, ws = wa)[1]
        @test A == D == pub5g(r12.z, r12.y, RM5; flag_He = true) && B == C == pub5g(r12.z - 5.0, y2, RM5; flag_He = true)
    end
    @testset "model contract: same model type + table dimensions use the workspace (no model values cached); other types fall back" begin
        e = RM5.eff
        eff2 = CosmoRec.EffectiveModel(e.htable, e.hetable, e.dp, e.bitot, e.fc, e.hatom, e.heatom, e.hconst, e.heconst, e.hiabs, !e.spin_forbidden, e.hi_absorption)
        rm2 = RecombinationModel(eff2, ACC5)
        @test rm2 isa typeof(RM5)
        ws = CosmoRec._rhs_workspace(r12.y, r12.z, RM5)
        f2 = zeros(12); used = CosmoRec._recombination_rhs_ws!(f2, r12.z, r12.y, rm2, ws; flag_He = true)
        @test used && f2 == pub5g(r12.z, r12.y, rm2; flag_He = true) && f2 != pub5g(r12.z, r12.y, RM5; flag_He = true)
        p0 = recombination_pass(RM5, (a...) -> solve_phase5(a...; reltol = 1.0e-12, a1 = 1.0e-18, aex = 1.0e-14, alg = Rodas5P(), tstops = false))
        rmd = with_diffusion(RM5, hi_diffusion_stage(RM5, HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6), pass_solution_rows(p0))[2])   # different model type
        f3 = zeros(12)
        @test !(rmd isa typeof(RM5)) && !CosmoRec._recombination_rhs_ws!(f3, r12.z, r12.y, rmd, ws; flag_He = true) && f3 == pub5g(r12.z, r12.y, rmd; flag_He = true)
    end
    @testset "prepared DI-Mooncake with a FRESH workspace constructed inside the objective" begin
        for r in (r7, r12)
            n = length(r.y); rng = MersenneTwister(11); w1 = randn(rng, n); w2 = randn(rng, n); y2s = 1 .+ 1.0e-6 .* (1:n)
            fnew(y) = dot(w1, rhs5g(r.z, y, RM5; flag_He = r.flag)[1]) + dot(w2, rhs5g(r.z, y .* y2s, RM5; flag_He = r.flag)[1])   # ws built in each call
            fpb(y) = dot(w1, pub5g(r.z, y, RM5; flag_He = r.flag)) + dot(w2, pub5g(r.z, y .* y2s, RM5; flag_He = r.flag))
            ref = ForwardDiff.gradient(fpb, r.y)
            prep = prepare_gradient(fnew, MC5G, r.y); g1 = gradient(fnew, prep, MC5G, r.y); g2 = gradient(fnew, prep, MC5G, r.y)
            @test maximum(abs.(g1 .- ref)) / max(maximum(abs, ref), floatmin()) < 1e-12 && g1 == g2
        end
    end
    @testset "allocations per call: private workspace path below the public path" begin
        for r in (r12, r7)
            f = similar(r.y); ws = CosmoRec._rhs_workspace(r.y, r.z, RM5)
            aw() = @allocated CosmoRec._recombination_rhs_ws!(f, r.z, r.y, RM5, ws; flag_He = r.flag)
            ap() = @allocated recombination_rhs!(f, r.z, r.y, RM5; flag_He = r.flag)
            aw(); ap(); aw(); ap()
            a_w = aw(); a_p = ap(); println("Chunk5g allocated bytes per call $(length(r.y))-state: workspace $a_w, public $a_p")
            @test a_w < a_p
        end
    end
end
