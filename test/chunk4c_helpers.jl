# Chunk 4c helpers: native Recfast++ RHS/Saha fixture, native history, SciML solve callbacks (test-only: the package itself has no ODE dependency).
using CosmoRec
using SHA
@isdefined(HIST4B) || include("chunk4b_helpers.jl")

const FIX4C = joinpath(@__DIR__, "fixtures", "native_recfast_rhs.txt")
const FIX4C_SHA256 = "81be03e49f5270c76121bddd75846294d8739c592bd970e98cb21476de350ea6"

function read_fixture4c(path = FIX4C)
    kv(f, s = 2) = Dict(String(f[k]) => parse(Float64, f[k + 1]) for k in s:2:length(f))
    d = Dict{String,Any}("RFR" => Vector{Vector{Float64}}(), "RFS" => Vector{Vector{Float64}}(), "RFH" => Vector{Vector{Float64}}(), "RFC" => Dict{String,Float64}(), "RFV" => Dict{String,Float64}(), "CFG" => Dict{String,Float64}())
    for l in eachline(path)
        (isempty(l) || startswith(l, "#")) && continue
        f = split(l, '\t')
        t = f[1]
        if t == "RFR"          # i pert z y1..y4 | Hz | NHrf | TCMB | f1..f4
            push!(d[t], parse.(Float64, f[[2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 17, 18, 19]]))
        elseif t == "RFS"      # z Xe TeQS SahaHeIII SahaHeII nH TR
            push!(d[t], parse.(Float64, f[[2, 3, 5, 7, 9, 11, 13]]))
        elseif t == "RFH"      # z Hz Hloc t_cos_rad
            push!(d[t], parse.(Float64, f[[2, 4, 6, 8]]))
        elseif t in ("RFC", "RFV")
            merge!(d[t], kv(f))
        elseif t == "CFG"
            d[t][String(f[2])] = parse(Float64, f[3])
        end
    end
    return d
end

const FX4C = read_fixture4c()

# native production parameters θ = [F, A2s1s, Yp, Omega_b, h100, T0]
theta4c() = [FX4C["RFC"]["F_rf"], FX4C["RFC"]["A2s1s_rf"], FX4C["CFG"]["Yp"], FX4C["CFG"]["Omega_b"], FX4C["CFG"]["h100"], FX4C["CFG"]["T0"]]
hfun4c() = (a = CosmosAccessors(constants4b(), nothing, hubble_table(HUB4B[:, 1], HUB4B[:, 2])); z -> cosmos_H(a, z))

using OrdinaryDiffEqRosenbrock: Rodas5P
using SciMLBase: ODEProblem, solve, Val
import SciMLBase
using SciMLSensitivity: GaussAdjoint, QuadratureAdjoint, MooncakeVJP

using ForwardDiff
# error control on PRIMAL values only: the default norm differentiates sqrt(sum(abs2)) at an exactly zero local error estimate (Inf*0 = NaN partials, `retcode = Unstable`)
_pv(x) = x isa ForwardDiff.Dual ? _pv(ForwardDiff.value(x)) : x
primal_norm(u::AbstractArray, t) = sqrt(sum(abs2 ∘ _pv, u) / max(length(u), 1))
primal_norm(u::Number, t) = abs(_pv(u))

"""Stiff SciML callback for `recfast_history` (descending z): `Rodas5P`, states saved at `znodes`."""
function solve_rodas(rhs!, u0, zs, znodes, p; alg = Rodas5P(), reltol = 1.0e-10, abstol = [1.0e-14, 1.0e-14, 1.0e-10], sensealg = nothing)
    prob = ODEProblem(rhs!, u0, (zs, znodes[end]), p)
    sol = sensealg === nothing ? solve(prob, alg; reltol = reltol, abstol = abstol, saveat = znodes, internalnorm = primal_norm) :
          solve(prob, alg; reltol = reltol, abstol = abstol, saveat = znodes, sensealg = sensealg, internalnorm = primal_norm)
    SciMLBase.successful_retcode(sol) || error("ODE solve failed: retcode = $(sol.retcode) after $(length(sol.u)) saved nodes")
    return reduce(hcat, sol.u)
end

"""Rodas5P with dense output (no saveat): returns `(states, derivatives)` at `znodes` (derivative of the dense interpolant), the analogue of the native solver's `Sz.dy`."""
function solve_rodas_dense(rhs!, u0, zs, znodes, p; alg = Rodas5P(), reltol = 1.0e-10, abstol = [1.0e-14, 1.0e-14, 1.0e-10])
    prob = ODEProblem(rhs!, u0, (zs, znodes[end]), p)
    sol = solve(prob, alg; reltol = reltol, abstol = abstol, dense = true)
    u = reduce(hcat, [sol(z) for z in znodes])
    du = reduce(hcat, [sol(z, Val{1}) for z in znodes])
    return u, du
end

# ---- reverse-mode (Mooncake via GaussAdjoint) route: plain-function RHS with explicit Jacobian and time-gradient (no captured state, no solver-internal AD of the adjoint RHS)
using SciMLBase: ODEFunction
using ADTypes: AutoFiniteDiff
const RFCONST4C = NATIVE_RECFAST_CONSTANTS
const RFHFUN4C = hfun4c()
rhs_ctx4c!(du, u, θ, z) = recfast_rhs3!(du, u, (θ, RFCONST4C, RFHFUN4C), z)
function jac_ctx4c!(Jm, u, θ, z)
    Jm .= ForwardDiff.jacobian(uu -> (du = similar(uu); rhs_ctx4c!(du, uu, θ, z); du), u); return nothing
end
function tgrad_ctx4c!(dT, u, θ, z)
    dT .= ForwardDiff.derivative(zz -> (du = zeros(promote_type(eltype(u), typeof(zz), eltype(θ)), 3); rhs_ctx4c!(du, u, θ, zz); du), z); return nothing
end
const ODEFUN4C = ODEFunction(rhs_ctx4c!; jac = jac_ctx4c!, tgrad = tgrad_ctx4c!)
# Rodas5P with finite-difference solver-internal autodiff is only used where the analytic Jacobian/time gradient above are supplied (the finite-difference Jacobian of the raw RHS is too noisy: Unstable)
const ALG_REV4C = Rodas5P(autodiff = AutoFiniteDiff())
