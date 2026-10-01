# Chunk 3e AD tests of the explicit DPesc_coh fallback: ForwardDiff vs 256-bit central differences and directional derivatives, prepared Mooncake
# VJPs (DifferentiationInterface, two independent preparations reused across changed native fixture parameters), singlet / triplet / combined.
# DIFFERENTIATED: the smooth map (tauS, eta_c, T, pd) -> Pesc - PS with the Patterson rule of each of the seven sub-intervals PINNED to the order the
# native stopping rule selects at the primal point. NOT claimed: differentiability across Patterson order switches (the stopping test is a discrete
# decision on primal values), across the wing switch |x| = 30, the Doppler-core/edge cuts of sigma(nu), or the DP-table seams.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(FX3E) || include("chunk3e_helpers.jl")

const MC3E = AutoMooncake(; config = nothing)
const OBS3EAD = Dict{String,Float64}()
rec3ead(k, v) = (OBS3EAD[k] = max(get(OBS3EAD, k, 0.0), v))

pinned_orders(q) = begin
    line = q.triplet ? NATIVE_DPESC.triplet : NATIVE_DPESC.singlet
    Pd = sobolev_p(q.pd * q.tauS)
    dpesc_appr_I_sym(q.T, q.pd * q.tauS, q.eta, line, NATIVE_DPESC.hlyc, Pd * 1e-5)[2]
end

# x = (tauS, eta_c, T, pd)
f_chan(q, ord) = x -> [dpesc_coh(NATIVE_DPESC, q.triplet, x[1], x[2], x[3], x[3], x[4]; orders = ord)[1]]
x_of(q) = [q.tauS, q.eta, q.T, q.pd]

function cdiff3e(f, x)
    setprecision(BigFloat, 256) do
        xb = BigFloat.(x)
        cols = map(eachindex(x)) do j
            h = BigFloat(1.0e-20) * abs(xb[j]); xp = copy(xb); xm = copy(xb); xp[j] += h; xm[j] -= h
            Float64.((f(xp) .- f(xm)) ./ (2h))
        end
        return reduce(hcat, cols)
    end
end
colscale3e(J) = [max(maximum(abs, J[:, j]), floatmin()) for j in axes(J, 2)]

function check3e(name, f, x0; fdtol, nseeds = 3, varied = 2, pert = 1.0e-5, rng = Random.Xoshiro(20261001))
    nout = length(f(x0))
    J = ForwardDiff.jacobian(f, x0)
    @test all(isfinite, J)
    Jfd = cdiff3e(f, x0)
    Jbig = setprecision(BigFloat, 256) do
        Float64.(ForwardDiff.jacobian(f, BigFloat.(x0)))
    end
    ebig = maximum(abs.(Jbig .- Jfd) ./ colscale3e(Jbig)')
    rec3ead("BIGAD:" * name, ebig)
    @test ebig < 1.0e-13
    efd = maximum(abs.(J .- Jfd) ./ colscale3e(J)')
    rec3ead("FD:" * name, efd)
    @test efd < fdtol
    v = randn(rng, length(x0)) .* abs.(x0)
    dd = ForwardDiff.derivative(t -> f(x0 .+ t .* v), 0.0)
    rec3ead("DIR:" * name, maximum(abs.(dd .- J * v)) / max(maximum(abs, J * v), floatmin()))
    @test isapprox(dd, J * v; rtol = 1.0e-10, atol = 1.0e-10 * maximum(abs, J * v))
    obj(x, w) = dot(w, f(x))
    prepA = prepare_gradient(obj, MC3E, x0, Constant(randn(rng, nout)))
    prepB = prepare_gradient(obj, MC3E, x0, Constant(randn(rng, nout)))
    for kk in 1:(varied + 1)
        x = kk == 1 ? x0 : x0 .* (1 .+ pert .* randn(rng, length(x0)))
        Jx = ForwardDiff.jacobian(f, x)
        for s in 1:nseeds
            w = randn(rng, nout)
            ref = Jx' * w
            gA = gradient(obj, prepA, MC3E, x, Constant(w))
            gB = gradient(obj, prepB, MC3E, x, Constant(w))
            sc = max(maximum(abs, ref), floatmin())
            rec3ead("VJP:" * name, max(maximum(abs.(gA .- ref)) / sc, maximum(abs.(gB .- ref)) / sc))
            @test isapprox(gA, ref; rtol = 1.0e-8, atol = 1.0e-8 * sc)
            @test isapprox(gB, ref; rtol = 1.0e-8, atol = 1.0e-8 * sc)
            @test gA == gB
        end
    end
    return J
end

@testset "Chunk 3e: DPesc_coh fallback AD (ForwardDiff, prepared Mooncake VJP)" begin
    for lab in ("z1500", "z1200", "z800")
        qS = first(filter(q -> !q.triplet && q.label == lab && q.code == "000", FX3E.queries))
        qT = first(filter(q -> q.triplet && q.label == lab && q.code == "000", FX3E.queries))
        oS, oT = pinned_orders(qS), pinned_orders(qT)
        @testset "$lab singlet" begin check3e("S:$lab", f_chan(qS, oS), x_of(qS); fdtol = 1.0e-6) end
        @testset "$lab triplet" begin check3e("T:$lab", f_chan(qT, oT), x_of(qT); fdtol = 1.0e-5) end
        @testset "$lab combined" begin
            fc = x -> vcat(f_chan(qS, oS)(x[1:4]), f_chan(qT, oT)(x[5:8]))
            check3e("C:$lab", fc, vcat(x_of(qS), x_of(qT)); fdtol = 1.0e-5)
        end
    end
    @testset "adaptive (primal-selected order) and pinned-order derivatives agree" begin
        q = first(filter(q -> !q.triplet && q.label == "z1200" && q.code == "000", FX3E.queries))
        x0 = x_of(q); ord = pinned_orders(q)
        Ja = ForwardDiff.jacobian(f_chan(q, nothing), x0)
        Jp = ForwardDiff.jacobian(f_chan(q, ord), x0)
        @test Ja == Jp
    end
    @testset "one-sided behaviour across a Patterson order switch" begin
        # the primal order of sub-interval 1 changes with the tolerance; the derivative of the pinned-order map is different at different orders
        q = first(filter(q -> !q.triplet && q.label == "z1200" && q.code == "000", FX3E.queries))
        x0 = x_of(q); ord = pinned_orders(q)
        lo = ntuple(k -> k == 1 ? max(ord[1] - 1, 1) : ord[k], 7)
        Jhi = ForwardDiff.jacobian(f_chan(q, ord), x0); Jlo = ForwardDiff.jacobian(f_chan(q, lo), x0)
        rec3ead("switch:jump", maximum(abs.(Jhi .- Jlo)) / maximum(abs, Jhi))
        @test all(isfinite, Jlo)            # documented: derivative differs by the quadrature error; no global claim
    end
    @testset "info" begin
        foreach(kv -> println("Chunk3e-AD ", kv[1], " = ", kv[2]), sort!(collect(OBS3EAD)))
    end
end
