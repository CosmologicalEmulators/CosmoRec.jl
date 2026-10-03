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
    @testset "info" begin
        foreach(kv -> println("Chunk5a-AD ", kv[1], " = ", kv[2]), sort!(collect(OBS5AAD)))
    end
end
