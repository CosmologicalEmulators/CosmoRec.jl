# Chunks 7-10 AD tests on dynamic inputs: the 24000 stored population values of the previous pass (rows Xe, X1s ... 3d, rho) feeding
#   (7b) the PDE coefficients A, B, C, D and pd_eff / Dnem_eff at z = 1300,
#   (7c) one theta-scheme step from a native spectrum,
#   (8a) a short production march (zs = 2500 -> 2470, 3 steps) with the correction integrals at its outputs (z = 2480, 2470),
# and (9a/10a) the PDE outputs (DF node values) feeding the next ODE right-hand side through the feedback splines.
# ForwardDiff directional derivatives vs 256-bit central differences; prepared Mooncake VJPs vs ForwardDiff (dot-product test <g, v> = <w, J v>,
# two independent preparations, changed inputs). Piecewise: spline cells, the Patterson stopping decisions (primal), the coefficient-grid size, the
# polint_JC stencil choice and the feedback range ends are frozen; no claim across them. REQUIRES COSMOREC_NATIVE_DATA_DIR.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(SETUP7) || include("chunk7_helpers.jl")

const MC710 = AutoMooncake(; config = nothing)
const OBS710 = Dict{String,Float64}()
rec710(k, v) = (OBS710[k] = max(get(OBS710, k, 0.0), v))
const ROWS710 = NODES5[:, 1:9]
const X0710 = vec(ROWS710[:, 2:9])
const IDX710 = collect(1:59:2353)
const FX7C710 = read_fixture7(FIX7C)
const YIN710 = vec7(FX7C710["YIN_50"][1])
const STEP710 = only(r for r in FX7C710["STEP"] if r[1] == 50.0)

function model710(x)
    rows = hcat(ROWS710[:, 1], reshape(x, size(ROWS710, 1), 8))
    pops = HIPopulationSplines(rows; zs = ZS6, ze = ZE6)
    cs = hi_pde_coefficients(rows, pops, ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6)
    return HIPDEModel(SETUP7, pops, cs, ACC5)
end
function f7b(x)
    st = hi_pde_rhs_coefficients(model710(x), 1300.0)
    T = eltype(x); out = T[]
    for v in (st.A, st.B, st.C, st.D), i in IDX710
        push!(out, v[i])
    end
    for v in (st.pd_eff, st.Dnem_eff), j in 1:2
        push!(out, v[j])
    end
    return out
end
function f7c(x)
    m = model710(x); T = eltype(x)
    S = HIPDEStepper{T}(SETUP7.x)
    y = T.(YIN710)
    hi_pde_step!(y, S, m, STEP710[4], STEP710[2], STEP710[3], STEP710[8], zero(T))
    return y[IDX710]
end
function f8a(x)
    o = hi_pde_corrections(model710(x); zs = 2500.0, ze = 2470.0, T = eltype(x))
    return vcat(o.DI1_2s, o.DF_2g[1], o.DF_2g[2], o.DF_R[1], o.y[IDX710])
end
# 9a/10a: DF node values -> feedback splines -> ODE right-hand side at explicit states (flag_He = 0 and 1)
const FX8A710 = read_fixture7(FIX8A)
const DFZ710 = vec7(FX8A710["DF_z"][1])
const DF0710 = vcat(vec7(FX8A710["DI1_2s"][1]), vec7(FX8A710["DF_2g"][1]), vec7(FX8A710["DF_2g"][2]), vec7(FX8A710["DF_R"][1]))
const STATES710 = let r5 = read_fixture7(joinpath(@__DIR__, "fixtures", "native_ode_rhs.txt"))["F5A"]
    [(z = r[3], flag = r[2] == 1.0, y = r[9:(8 + ode_nstate(r[2] == 1.0))]) for r in r5 if 600.0 <= r[3] <= 1900.0]
end
function f9(df)
    n = length(DFZ710)
    fb = hi_diffusion_feedback(DFZ710, df[1:n], [df[(n + 1):(2n)], df[(2n + 1):(3n)]], [df[(3n + 1):(4n)]])
    rm = with_diffusion(RM5, fb); T = eltype(df); out = T[]
    for s in STATES710
        append!(out, recombination_rhs(s.z, s.y, rm; flag_He = s.flag))
    end
    return out
end

function dir_check(name, f, x0, v; h = BigFloat(1.0e-30), floor_rel = 1e-10)
    dd = ForwardDiff.derivative(t -> f(x0 .+ t .* v), 0.0)
    ref = setprecision(BigFloat, 256) do
        xb = BigFloat.(x0); vb = BigFloat.(v)
        Float64.((f(xb .+ h .* vb) .- f(xb .- h .* vb)) ./ (2h))
    end
    @test all(isfinite, dd)
    e = maximum(abs.(dd .- ref) ./ max.(abs.(ref), floor_rel * maximum(abs, ref)))
    rec710("DIR:" * name, e)
    return e
end
function vjp_check(name, f, x0; seed = 1)
    rng = Random.Xoshiro(seed)
    nout = length(f(x0))
    obj(x, w) = dot(w, f(x))
    prepA = prepare_gradient(obj, MC710, x0, Constant(randn(rng, nout)))
    prepB = prepare_gradient(obj, MC710, x0, Constant(randn(rng, nout)))
    for kk in 1:2
        x = kk == 1 ? x0 : x0 .* (1 .+ 1.0e-9 .* randn(rng, length(x0)))
        w = randn(rng, nout); v = randn(rng, length(x0)) .* abs.(x0)
        gA = gradient(obj, prepA, MC710, x, Constant(w)); gB = gradient(obj, prepB, MC710, x, Constant(w))
        @test gA == gB
        jv = ForwardDiff.derivative(t -> f(x .+ t .* v), 0.0)
        lhs = dot(gA, v); rhs = dot(w, jv)
        e = abs(lhs - rhs) / max(sum(abs.(gA .* v)), sum(abs.(w .* jv)))      # relative to the magnitude of the summands
        rec710("VJP:" * name, e)
    end
end

@testset "Chunks 7-10: PDE stage and feedback AD (dynamic inputs)" begin
    rng = Random.Xoshiro(20261002)
    v = randn(rng, length(X0710)) .* abs.(X0710)
    @testset "7b coefficients: directional vs 256-bit central difference" begin
        @test dir_check("7b", f7b, X0710, v) < 1e-8
    end
    @testset "7c one step: directional vs 256-bit central difference" begin
        @test dir_check("7c", f7c, X0710, v) < 1e-8
    end
    @testset "8a short march + integrals: directional vs 256-bit central difference" begin
        @test dir_check("8a", f8a, X0710, v) < 1e-6
    end
    @testset "9a feedback -> ODE RHS: directional vs 256-bit central difference" begin
        vd = randn(rng, length(DF0710)) .* abs.(DF0710)
        @test dir_check("9a", f9, DF0710, vd) < 1e-8
    end
    @testset "prepared Mooncake VJPs vs ForwardDiff (dot-product test)" begin
        vjp_check("7b", f7b, X0710; seed = 11)
        vjp_check("7c", f7c, X0710; seed = 12)
        vjp_check("8a", f8a, X0710; seed = 13)
        vjp_check("9a", f9, DF0710; seed = 14)
        for k in ("7b", "7c", "8a", "9a")
            @test OBS710["VJP:" * k] < 1e-10
        end
    end
    @testset "info" begin
        foreach(kv -> println("Chunk7-10-AD ", kv[1], " = ", kv[2]), sort!(collect(OBS710)))
    end
end
