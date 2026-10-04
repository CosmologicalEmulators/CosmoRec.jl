# Chunk 5a AD tests of the ODE right-hand side (full tables): ForwardDiff Jacobians in the state (12 x 12, 7 x 7) vs 256-bit central differences, and prepared Mooncake VJPs (two independent
# preparations, changed states). Interior points only (table stencil seams, the DP table domain, the Patterson order and the absorber switch z = 3400 are piecewise).
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(RM5) || include("chunk5_helpers.jl")
@isdefined(ROWS5A) || include("chunk5a_ode_rhs_native.jl")

const MC5A = AutoMooncake(; config = nothing)
const OBS5AAD = Dict{String,Float64}()
rec5aad(k, v) = (OBS5AAD[k] = max(get(OBS5AAD, k, 0.0), v))
colscale5a(J) = [max(maximum(abs, J[:, j]), floatmin()) for j in axes(J, 2)]

function cdiff5a(f, x)
    setprecision(BigFloat, 256) do
        xb = BigFloat.(x)
        cols = map(eachindex(x)) do j
            h = BigFloat(1.0e-20) * max(abs(xb[j]), BigFloat(1.0e-30)); xp = copy(xb); xm = copy(xb); xp[j] += h; xm[j] -= h
            Float64.((f(xp) .- f(xm)) ./ (2h))
        end
        return reduce(hcat, cols)
    end
end

@testset "Chunk 5a: ODE right-hand side AD (ForwardDiff vs 256-bit FD, prepared Mooncake VJP)" begin
    for (label, flag, z) in (("flag_He=1 z=2500", true, 2500.0), ("flag_He=1 z=1500 (off-table DP fallback)", true, 1500.0), ("flag_He=1 z=1200", true, 1200.0), ("flag_He=0 z=400", false, 400.0), ("flag_He=0 z=100", false, 100.0))
        r = first(r for r in ROWS5A if r.flag == flag && r.z == z && r.pat == 0)
        f = y -> recombination_rhs(z, y, RM5; flag_He = flag)
        @testset "$label" begin
            J = ForwardDiff.jacobian(f, r.y)
            @test all(isfinite, J)
            Jfd = cdiff5a(f, r.y)
            Jbig = setprecision(BigFloat, 256) do
                Float64.(ForwardDiff.jacobian(f, BigFloat.(r.y)))
            end
            ebig = maximum(abs.(Jbig .- Jfd) ./ colscale5a(Jbig)'); rec5aad("BIGAD:" * label, ebig); @test ebig < 1e-12
            efd = maximum(abs.(J .- Jfd) ./ colscale5a(J)'); rec5aad("FD:" * label, efd); @test efd < 1e-12
            rng = Random.Xoshiro(20261010)
            nout = length(r.y)
            obj(x, w) = dot(w, f(x))
            prepA = prepare_gradient(obj, MC5A, r.y, Constant(randn(rng, nout)))
            prepB = prepare_gradient(obj, MC5A, r.y, Constant(randn(rng, nout)))
            for kk in 1:2
                y = kk == 1 ? r.y : r.y .* (1 .+ 1.0e-6 .* randn(rng, nout))
                w = randn(rng, nout)
                ref = ForwardDiff.jacobian(f, y)' * w
                gA = gradient(obj, prepA, MC5A, y, Constant(w)); gB = gradient(obj, prepB, MC5A, y, Constant(w))
                e = max(maximum(abs.(gA .- ref)), maximum(abs.(gB .- ref))) / max(maximum(abs, ref), floatmin()); rec5aad("VJP:" * label, e)
                @test e < 1e-12
                @test gA == gB
            end
        end
    end
    @testset "prepared per-solve state Jacobian (PreparedJac5) vs the unprepared oracle jac5!" begin
        # element type a solve would carry for this (u, p, z): the solver state is promoted to it
        Tof(u, p, z) = promote_type(eltype(u), typeof(z), typeof(p.hscale), typeof(p.nbscale))
        oracle(u, p, z) = (T = Tof(u, p, z); J = zeros(T, length(u), length(u)); jac5!(J, T.(u), p, z); J)
        function prepared(u, p, z)                       # a fresh cache built from these actual types, as prob_phase5 does per solve
            T = Tof(u, p, z); uu = T.(u); j = PreparedJac5(uu, p, z); J = zeros(T, length(u), length(u)); j(J, uu, p, z)
            @test j.nfallback[] == 0
            return J
        end
        pofh(flag, h) = RecombinationODEParams(RM5, flag, h, one(h))
        for (flag, z) in ((true, 2500.0), (true, 1500.0), (false, 400.0))     # 12-state in-table, 12-state off-table DP fallback, 7-state
            r = first(r for r in ROWS5A if r.flag == flag && r.z == z && r.pat == 0); p = RecombinationODEParams(RM5, flag)
            @test prepared(r.y, p, z) == oracle(r.y, p, z)                                                   # primal Jacobian, bitwise
            @test ForwardDiff.jacobian(u -> vec(prepared(u, p, z)), r.y) == ForwardDiff.jacobian(u -> vec(oracle(u, p, z)), r.y)          # d J / d y
            @test ForwardDiff.derivative(zz -> vec(prepared(r.y, p, zz)), z) == ForwardDiff.derivative(zz -> vec(oracle(r.y, p, zz)), z)  # d J / d z
            @test ForwardDiff.derivative(h -> vec(prepared(r.y, pofh(flag, h), z)), 1.0) == ForwardDiff.derivative(h -> vec(oracle(r.y, pofh(flag, h), z)), 1.0)   # d J / d hscale
            @test setprecision(() -> prepared(BigFloat.(r.y), p, z) == oracle(BigFloat.(r.y), p, z), BigFloat, 256)       # BigFloat state (Float64 z, as in the 256-bit reference above)
            # one cache, LIVE z/p alternating (no stale values), and two independent caches interleaved
            z2 = z - 7.0; y2 = r.y .* (1 .+ 1.0e-6 .* (1:length(r.y))); p2 = pofh(flag, 1.001)
            j = PreparedJac5(r.y, p, z); k = PreparedJac5(r.y, p, z); A = similar(r.y, length(r.y), length(r.y)); B = similar(A); C = similar(A); D = similar(A)
            j(A, r.y, p, z); k(D, y2, p2, z2); j(B, y2, p2, z2); j(C, r.y, p, z)
            @test A == C == oracle(r.y, p, z) && B == D == oracle(y2, p2, z2) && A != B
            @test j.nfallback[] == 0 && k.nfallback[] == 0
            # a call whose types differ from the prepared ones (Dual z into a Float64-z cache) takes the unprepared oracle, counted, no conversion
            jd = PreparedJac5(r.y, p, z)
            dfb = ForwardDiff.derivative(zz -> (Jd = zeros(typeof(zz), length(r.y), length(r.y)); jd(Jd, typeof(zz).(r.y), p, zz); vec(Jd)), z)
            @test jd.nfallback[] == 1 && dfb == ForwardDiff.derivative(zz -> vec(oracle(r.y, p, zz)), z)
        end
    end
    @testset "info" begin
        foreach(kv -> println("Chunk5a-AD ", kv[1], " = ", kv[2]), sort!(collect(OBS5AAD)))
    end
end
