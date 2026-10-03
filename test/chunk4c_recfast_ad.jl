# Chunk 4c AD tests through the history solve: ForwardDiff (through `solve`, Rodas5P) vs finite differences (Float64 with step convergence; the closed-form Saha part also vs
# 256-bit central differences) and a prepared Mooncake VJP (Rodas5P(autodiff = AutoFiniteDiff) + GaussAdjoint(MooncakeVJP), the Chunk 1a route) vs the ForwardDiff projection.
# DIFFERENTIATED: θ = [F, A2s1s, Yp, Omega_b, h100, T0] -> history values at fixed nodes. The node grid and the Saha-segment breaks are a DISCRETE operator evaluated on primal values and
# frozen here (T0 enters physics but the frozen nodes do not move); no derivative of the segment-break location, of the solver steps or of the native Saha switches is claimed.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(FX4C) || include("chunk4c_helpers.jl")

const MC4C = AutoMooncake(; config = nothing)
const OBS4CAD = Dict{String,Float64}()
rec4cad(k, v) = (OBS4CAD[k] = max(get(OBS4CAD, k, 0.0), v))

const C4 = NATIVE_RECFAST_CONSTANTS
const H4 = hfun4c()
const GRID4 = recfast_grid(theta4c(), C4, H4)
const ZSEL4 = [2500.0, 1800.0, 1200.0, 1000.0, 800.0, 500.0, 200.0, 50.0, 5.0]
const GRIDSEL4 = let j = GRID4.jode
    zs = [GRID4.z[argmin(abs.(GRID4.z .- zt))] for zt in ZSEL4]
    RecfastGrid(vcat(GRID4.z[1:j], zs), GRID4.seg_end, j)
end
const KSAHA4 = 1000                       # a Saha-segment-2 node: closed-form dependence on θ

hist_out(θ, alg, sensealg, tol) = begin
    h = recfast_history(θ, C4, H4, GRIDSEL4, (rhs!, u0, zs, zn, p) -> solve_rodas(rhs!, u0, zs, zn, p; alg = alg, sensealg = sensealg, reltol = tol[1], abstol = tol[2]))
    j = GRIDSEL4.jode; ode = (j + 1):length(h.z)
    vcat(h.Xe[KSAHA4], h.TM[KSAHA4], h.Xe_He[ode], h.Xe_H[ode], h.TM[ode])
end
f_fd(θ) = hist_out(θ, Rodas5P(), nothing, (1.0e-12, [1.0e-16, 1.0e-16, 1.0e-12]))
f_ad(θ) = hist_out(θ, Rodas5P(), nothing, (1.0e-12, [1.0e-16, 1.0e-16, 1.0e-12]))
function f_rev(θ)
    h = recfast_history(θ, C4, H4, GRIDSEL4, (rhs!, u0, zs, zn, p) -> solve_rodas(ODEFUN4C, u0, zs, zn, p[1]; alg = ALG_REV4C, sensealg = GaussAdjoint(autojacvec = MooncakeVJP()),
                                                                                  reltol = 1.0e-10, abstol = [1.0e-14, 1.0e-14, 1.0e-10]))
    j = GRIDSEL4.jode; ode = (j + 1):length(h.z)
    vcat(h.Xe[KSAHA4], h.TM[KSAHA4], h.Xe_He[ode], h.Xe_H[ode], h.TM[ode])
end

function fd_jac(f, θ, hrel)
    cols = map(eachindex(θ)) do i
        h = hrel * abs(θ[i]); tp = copy(θ); tm = copy(θ); tp[i] += h; tm[i] -= h
        (f(tp) .- f(tm)) ./ (2h)
    end
    return reduce(hcat, cols)
end
colscale4c(J) = [max(maximum(abs, J[:, j]), floatmin()) for j in axes(J, 2)]

@testset "Chunk 4c: Recfast++ history AD (ForwardDiff through solve, finite differences, prepared Mooncake VJP)" begin
    θ0 = theta4c()
    @test length(f_ad(θ0)) == 2 + 3 * length(ZSEL4)
    @testset "ForwardDiff through the solve vs finite differences (step convergence)" begin
        J = ForwardDiff.jacobian(f_ad, θ0)
        @test all(isfinite, J)
        errs = Float64[]
        for hrel in (1.0e-3, 1.0e-4, 1.0e-5)
            Jfd = fd_jac(f_fd, θ0, hrel)
            push!(errs, maximum(abs.(J .- Jfd) ./ colscale4c(J)'))
        end
        rec4cad("FD_step_errors_max", maximum(errs)); rec4cad("FD_best_step_error", minimum(errs))
        @info "ForwardDiff vs central FD, col-scaled error for h = 1e-3, 1e-4, 1e-5" errs
        @test minimum(errs) < 1.0e-6                                      # best step is the coarsest h = 1e-3 (observed 1.6e-7): smaller steps are limited by the solver tolerance noise (1e-5 at h = 1e-5)
        @test errs[1] < 1.0e-6
        # directional derivative along a random direction vs a symmetric difference
        v = randn(Random.Xoshiro(20261005), 6) .* abs.(θ0)
        dd = ForwardDiff.derivative(t -> f_ad(θ0 .+ t .* v), 0.0)
        ref = (f_fd(θ0 .+ 1e-5 .* v) .- f_fd(θ0 .- 1e-5 .* v)) ./ 2e-5
        e = maximum(abs.(dd .- ref)) / maximum(abs, ref); rec4cad("JVP_dir", e); @test e < 1e-5
    end
    @testset "closed-form Saha node: ForwardDiff vs 256-bit central difference" begin
        g(θ) = (h = recfast_history(θ, C4, H4, RecfastGrid(GRIDSEL4.z[1:(KSAHA4 + 1)], (GRIDSEL4.seg_end[1], GRIDSEL4.seg_end[2], KSAHA4 + 1), KSAHA4 + 1), (a...) -> error("no ODE")); [h.Xe[KSAHA4], h.TM[KSAHA4]])
        J = ForwardDiff.jacobian(g, θ0)
        Jb, Jbad = setprecision(BigFloat, 256) do
            cols = map(1:6) do i
                hh = BigFloat(1.0e-20) * abs(BigFloat(θ0[i])); tp = BigFloat.(θ0); tm = BigFloat.(θ0); tp[i] += hh; tm[i] -= hh
                Float64.((g(tp) .- g(tm)) ./ (2hh))
            end
            (reduce(hcat, cols), Float64.(ForwardDiff.jacobian(g, BigFloat.(θ0))))
        end
        e = maximum(abs.(Jbad .- Jb) ./ colscale4c(Jb)'); rec4cad("BIGAD_vs_BIGFD_saha", e); @test e < 1e-13          # AD formulas exact (256-bit AD vs 256-bit FD)
        e = maximum(abs.(J .- Jb) ./ colscale4c(Jb)'); rec4cad("F64AD_vs_BIGFD_saha", e); @test e < 1e-8                # Float64: conditioning of d + sqrt(d^2 + ...) at |d| ~ 1e5
        @test J[1, 1] == 0 && J[1, 2] == 0           # the HeIII Saha node does not depend on F or A2s1s
    end
    @testset "prepared Mooncake VJP (GaussAdjoint, Rodas5P with analytic Jacobian/tgrad) vs ForwardDiff projection, independent preparations" begin
        rng = Random.Xoshiro(20261006)
        nout = length(f_ad(θ0))
        obj(θ, w) = dot(w, f_rev(θ))
        prepA = prepare_gradient(obj, MC4C, θ0, Constant(randn(rng, nout)))
        prepB = prepare_gradient(obj, MC4C, θ0, Constant(randn(rng, nout)))
        nsel = length(ZSEL4)
        blocks = (("all", collect(1:nout)), ("Xe_He", 3:(2 + nsel)), ("Xe_H", (3 + nsel):(2 + 2nsel)), ("TM", (3 + 2nsel):nout), ("saha_nodes", 1:2))
        for kk in 1:2
            θ = kk == 1 ? θ0 : θ0 .* (1 .+ 1.0e-4 .* randn(rng, 6))
            Jθ = ForwardDiff.jacobian(f_ad, θ)
            for (bn, rows) in blocks
                w = zeros(nout)
                bn == "all" ? (w .= randn(rng, nout)) : (w[rows] .= 1.0)
                ref = Jθ' * w
                gA = gradient(obj, prepA, MC4C, θ, Constant(w)); gB = gradient(obj, prepB, MC4C, θ, Constant(w))
                sc = maximum(abs, ref)
                e = max(maximum(abs.(gA .- ref)), maximum(abs.(gB .- ref))) / sc; rec4cad("VJP:" * bn, e)
                @test e < (bn == "all" || bn == "saha_nodes" || bn == "TM" ? 1.0e-8 : 2.0e-6)      # Xe_He/Xe_H blocks: tiny gradients, limited by the adjoint solve tolerance (observed 4.7e-7, 1.3e-7)
                @test gA == gB
            end
        end
    end
    @testset "documented limitations are real: the grid is discrete" begin
        @test abs(recfast_grid(θ0 .* (1 + 1e-12), C4, H4).jode - GRID4.jode) <= 1
    end
    @testset "info" begin
        foreach(kv -> println("Chunk4c-AD ", kv[1], " = ", kv[2]), sort!(collect(OBS4CAD)))
    end
end
