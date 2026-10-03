# Chunk 5e helpers (definitions shared by the AD test and the benchmark): parameter set p = [F, A2s1s, hscale, nbscale], full-pipeline and selected-node outputs, forward and reverse routes.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
include("chunk5_helpers.jl")

const MC5E = AutoMooncake(; config = nothing)
const OBS5E = Dict{String,Float64}()
rec5e(k, v) = (OBS5E[k] = max(get(OBS5E, k, 0.0), v))

# a1 = 1e-18 (was 1e-16, approved 2026-10-02): at a1 = 1e-16 the X1s weight a1 + reltol|X1s| is absolute-dominated for z > 1813 (X1s ~ 7e-10 at z = 3000), docs/ORIGINAL_VS_PORT_ACCURACY_INVESTIGATION.md 5.6
const TOL5E = (1.0e-12, 1.0e-18, 1.0e-14)
const SOLVEP5E = (a...) -> solve_phase5(a...; reltol = TOL5E[1], a1 = TOL5E[2], aex = TOL5E[3])
const THETA5E = theta5()
const P0_5E = [THETA5E[1], THETA5E[2], 1.0, 1.0]
const PASS0_5E = recombination_pass(RM5, SOLVEP5E)
const K5E = PASS0_5E.k_switch
const NZ5E = length(PASS0_5E.z)
const ZG5E = [2500.0, 1500.0, 1000.0, 800.0, 500.0, 100.0, 40.0, 5.0, 1.0]           # grid redshifts for the assembled history (inside the stored range)
const SEL5E = [findmin(abs.(PASS0_5E.z .- zt))[2] for zt in (2500.0, 2000.0, 1700.0, 1200.0, 600.0, 100.0)]
const TAILSEL5E = [20, 60, 120, 180, 199]

# ---- (i) ForwardDiff route: the complete pipeline, frozen switch node ----
function full_out(p)
    h = recombination_history(RM5, THETA5E_vec(p), SOLVEP5E, solve_tail5, ZG5E; hscale = p[3], nbscale = p[4], k_switch = K5E)
    return vcat(h.Xe, h.Te)
end
THETA5E_vec(p) = [p[1], p[2], THETA5E[3], THETA5E[4], THETA5E[5], THETA5E[6]]

# ---- (ii) reverse route: plain-function RHS with global constant model ----
# SciMLSensitivity probes the RHS with a placeholder (zero) state while preparing the Mooncake pullback; rho = 0 is outside the native tables (RateTableDomainError), so the test-only
# wrapper returns zeros for a non-positive rho (never reached by a physical state; a branch on the primal value).
function rhs12_rev!(du, u, pv, z)
    CosmoRec._primal(u[1]) > 0 || (fill!(du, 0); return du)
    return recombination_rhs!(du, z, u, RM5; flag_He = true, hscale = pv[1], nbscale = pv[2])
end
function rhs7_rev!(du, u, pv, z)
    CosmoRec._primal(u[1]) > 0 || (fill!(du, 0); return du)
    return recombination_rhs!(du, z, u, RM5; flag_He = false, hscale = pv[1], nbscale = pv[2])
end
function tail_rev!(du, u, pv, z)
    θ = [pv[1], pv[2], THETA5E[3], pv[3], THETA5E[5], THETA5E[6]]
    f = recfast_rhs(θ, z, u, pv[5] * cosmos_H(ACC5, z), NATIVE_RECFAST_CONSTANTS)
    du[1] = f[1]; du[2] = f[2] * pv[4]; du[3] = f[3]
    return du
end
function mkfun(rhs!)
    jac!(Jm, u, pv, z) = (Jm .= ForwardDiff.jacobian(uu -> (du = similar(uu); rhs!(du, uu, pv, z); du), u); nothing)
    tgr!(dT, u, pv, z) = (dT .= ForwardDiff.derivative(zz -> (du = zeros(promote_type(eltype(u), typeof(zz)), length(u)); rhs!(du, u, pv, zz); du), z); nothing)
    return ODEFunction(rhs!; jac = jac!, tgrad = tgr!)
end
const FUN12 = mkfun(rhs12_rev!); const FUN7 = mkfun(rhs7_rev!); const FUNT = mkfun(tail_rev!)
# Reverse mode DIFFERENTIATES THE SOLVER STEPS directly (discretize-then-differentiate): Mooncake through `solve` of Rodas5P with explicit Jacobian and time gradient and
# `FullSpecialize` (no sensealg). The continuous adjoints (GaussAdjoint/QuadratureAdjoint with MooncakeVJP, the Chunk 4c route) fail on this system: the adjoint solve is unstable
# (Jacobian entries up to 3e21, `dt` forced below epsilon, states growing to 1e59; diagnostics kept in chunk5/diag5e.log).
function solve_rev(fun, u0, z0, zn, pv, ab; reltol, dtmax = 10.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(fun, u0, (z0, zn[end]), pv)
    sol = solve(prob, Rodas5P(); reltol = reltol, abstol = ab, saveat = zn, dtmax = dtmax, internalnorm = primal_norm)
    SciMLBase.successful_retcode(sol) || error("reverse-route ODE solve failed: $(sol.retcode)")
    return reduce(hcat, sol.u)
end
const TOLREV5E = 1.0e-10

function sel_out(p, mode)
    hs, nb = p[3], p[4]
    solvep = mode == :fd ? SOLVEP5E : (f!, u0, z0, zn, pp) -> solve_rev(pp.flag ? FUN12 : FUN7, u0, z0, zn, [pp.hscale, pp.nbscale], abstol5(length(u0); a1 = 1.0e-14, aex = 1.0e-12); reltol = TOLREV5E)
    pass = recombination_pass(RM5, solvep; k_switch = K5E, nodes = SEL5E, hscale = hs, nbscale = nb)
    zi = pass.z[end]; y7 = pass.states[end]
    inp = recfast_tail_inputs(RM5, zi, y7; hscale = hs, nbscale = nb)
    θeff = [p[1], p[2], THETA5E[3], THETA5E[4] * nb, THETA5E[5], THETA5E[6]]
    Hfun = z -> hs * cosmos_H(ACC5, z)
    solvet = mode == :fd ? solve_tail5 : (f!, u0, z0, zn, pt) -> solve_rev(FUNT, u0, z0, zn, [pt.θ[1], pt.θ[2], pt.θ[4], pt.ff, hs], [1.0e-14, 1.0e-14, 1.0e-10]; reltol = TOLREV5E)
    tail = recfast_tail(θeff, NATIVE_RECFAST_CONSTANTS, Hfun, zi, inp.Xe_Hi, inp.Xe_Hei, inp.Xei, inp.TMi, inp.dXei, solvet)
    return vcat(pass.Xe, pass.Te, tail.Xe[TAILSEL5E], tail.TM[TAILSEL5E])
end

function fd_jac(f, p, hrel)
    cols = map(eachindex(p)) do i
        h = hrel * abs(p[i]); pp = copy(p); pm = copy(p); pp[i] += h; pm[i] -= h
        (f(pp) .- f(pm)) ./ (2h)
    end
    return reduce(hcat, cols)
end
colscale5e(J) = [max(maximum(abs, J[:, j]), floatmin()) for j in axes(J, 2)]
# error in the RELATIVE change of every output per RELATIVE change of every parameter (elasticity error): robust where a column is (nearly) zero (F cancels in the rescaled Recfast tail)
elast_err(J, Jfd, f0, p) = maximum(abs.(J .- Jfd) .* abs.(p)' ./ abs.(f0))

