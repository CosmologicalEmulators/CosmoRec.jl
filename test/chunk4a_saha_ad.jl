# Chunk 4a AD tests: ForwardDiff of the Saha initializer and ysol packing w.r.t. the explicit scalar inputs vs 256-bit central differences and directional
# derivatives; prepared Mooncake VJPs (two independent preparations, reused at changed inputs) vs ForwardDiff J'w. Interior points only: `min` clips (Xe at 1+fHe, Xp at 1)
# and the ground-state overrides are piecewise definitions (a clipped input has zero derivative); no derivative is claimed AT a clip edge.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(FX4A) || include("chunk4a_helpers.jl")

const MC4A = AutoMooncake(; config = nothing)
const OBS4AAD = Dict{String,Float64}()
rec4aad(k, v) = (OBS4AAD[k] = max(get(OBS4AAD, k, 0.0), v))

# p = (fHe, Xe_Seager, Xp_raw, NH, Te, Tg, XHeII, XHeI1s)
p_of(z) = (s = FX4A["SAHAIN"][z]; [FX4A["fHe"], s["Xe_clipped"], s["Xp_raw"], s["NH"], s["Te"], s["TCMB"], s["XHeII"], s["XHeI_ground"]])
inp4a(p) = SahaInputs(p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8])
f_X(p) = saha_initial_state(inp4a(p), NATIVE_SAHA_HYDROGEN, NATIVE_SAHA_HELIUM)
f_y(p) = (X = f_X(p); pack_ysol(X, 6, 7))

function cdiff4a(f, x)
    setprecision(BigFloat, 256) do
        xb = BigFloat.(x)
        cols = map(eachindex(x)) do j
            h = BigFloat(1.0e-20) * abs(xb[j]); xp = copy(xb); xm = copy(xb); xp[j] += h; xm[j] -= h
            Float64.((f(xp) .- f(xm)) ./ (2h))
        end
        return reduce(hcat, cols)
    end
end
colscale4a(J) = [max(maximum(abs, J[:, j]), floatmin()) for j in axes(J, 2)]

function check4a(name, f, x0; fdtol = 1.0e-12, nseeds = 3, varied = 2, pert = 1.0e-5, rng = Random.Xoshiro(20261002))
    nout = length(f(x0))
    J = ForwardDiff.jacobian(f, x0)
    @test all(isfinite, J)
    Jfd = cdiff4a(f, x0)
    Jbig = setprecision(BigFloat, 256) do
        Float64.(ForwardDiff.jacobian(f, BigFloat.(x0)))
    end
    ebig = maximum(abs.(Jbig .- Jfd) ./ colscale4a(Jbig)')
    rec4aad("BIGAD:" * name, ebig); @test ebig < 1.0e-13
    efd = maximum(abs.(J .- Jfd) ./ colscale4a(J)')
    rec4aad("FD:" * name, efd); @test efd < fdtol
    v = randn(rng, length(x0)) .* abs.(x0)
    dd = ForwardDiff.derivative(t -> f(x0 .+ t .* v), 0.0)
    rec4aad("DIR:" * name, maximum(abs.(dd .- J * v)) / max(maximum(abs, J * v), floatmin()))
    @test isapprox(dd, J * v; rtol = 1.0e-10, atol = 1.0e-10 * maximum(abs, J * v))
    obj(x, w) = dot(w, f(x))
    prepA = prepare_gradient(obj, MC4A, x0, Constant(randn(rng, nout)))
    prepB = prepare_gradient(obj, MC4A, x0, Constant(randn(rng, nout)))
    for kk in 1:(varied + 1)
        x = kk == 1 ? x0 : x0 .* (1 .+ pert .* randn(rng, length(x0)))
        Jx = ForwardDiff.jacobian(f, x)
        for s in 1:nseeds
            w = randn(rng, nout)
            ref = Jx' * w
            gA = gradient(obj, prepA, MC4A, x, Constant(w))
            gB = gradient(obj, prepB, MC4A, x, Constant(w))
            sc = max(maximum(abs, ref), floatmin())
            rec4aad("VJP:" * name, max(maximum(abs.(gA .- ref)) / sc, maximum(abs.(gB .- ref)) / sc))
            @test isapprox(gA, ref; rtol = 1.0e-9, atol = 1.0e-9 * sc)
            @test isapprox(gB, ref; rtol = 1.0e-9, atol = 1.0e-9 * sc)
            @test gA == gB
        end
    end
    return J
end

@testset "Chunk 4a: Saha initialization AD (ForwardDiff, prepared Mooncake VJP)" begin
    for z in (3000.0, 2000.0, 1500.0, 1200.0, 800.0)
        # z = 3000: Xp_raw = 1 - 7e-10 sits 7e-10 below the clip, so perturbed points would cross it: base point only
        kw = z == 3000.0 ? (; varied = 0) : (;)
        @testset "z = $z state vector" begin check4a("X:$z", f_X, p_of(z); kw...) end
        @testset "z = $z packed ysol" begin check4a("Y:$z", f_y, p_of(z); kw...) end
    end
    @testset "analytic structure: overrides and constant rho" begin
        p0 = p_of(1500.0)
        J = ForwardDiff.jacobian(f_X, p0)
        @test all(J[15, :] .== 0)                                   # rho == 1: no dependence
        @test J[2, 3] == -1.0 && count(!iszero, J[2, :]) == 1          # H 1s = 1 - Xp
        @test J[8, 8] == 1.0 && count(!iszero, J[8, :]) == 1           # He 1s = NHeI/NH input
        @test J[1, 2] == 1.0                                         # Xe slot follows Xe_Seager when unclipped
        # clipped Xe: X[1] independent of Xe_Seager (derivative 0 through min) while the helium levels keep their Xe dependence
        pc = copy(p0); pc[2] = 2.0
        Jc = ForwardDiff.jacobian(f_X, pc)
        @test Jc[1, 2] == 0.0 && Jc[1, 1] == 1.0 && Jc[9, 2] != 0.0 && Jc[3, 2] == 0.0
        # clipped Xp: H levels lose the Xp dependence
        pp = copy(p0); pp[3] = 1.5
        Jp = ForwardDiff.jacobian(f_X, pp)
        @test Jp[3, 3] == 0.0 && Jp[2, 3] == 0.0
    end
    @testset "info" begin
        foreach(kv -> println("Chunk4a-AD ", kv[1], " = ", kv[2]), sort!(collect(OBS4AAD)))
    end
end
