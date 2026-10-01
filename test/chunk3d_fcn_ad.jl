# Chunk 3d AD tests for the assembled pointwise RHS: ForwardDiff Jacobians (vs 256-bit central differences and directional derivatives) and prepared
# Mooncake VJPs (DifferentiationInterface, two independent preparations) vs ForwardDiff J'w, per component and for the full composition.
# Only the RIGHT-HAND-SIDE MAPPING is differentiated at smooth interior points; z is never differentiated (it only selects the absorber branch via its primal).
# The mapping is piecewise smooth: it is NOT claimed differentiable across the absorber switch z = 3400, the table stencil seams, or the DP-table domain edge.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(FX3D) || include("chunk3d_helpers.jl")

const MC3D = AutoMooncake(; config = nothing)
const OBS3DAD = Dict{String,Float64}()
rec3d(key, v) = (OBS3DAD[key] = max(get(OBS3DAD, key, 0.0), v))

# parameter vector p = (X[1:15] native layout, Tg, NH, Hz, fHe)
param3d(s) = vcat(s.X, s.Tg, s.NH, s.Hz, s.fHe)
bg3d(p) = EffectiveBackground(p[16], p[17], p[18], p[19])
fracs3d(p) = effective_fractions(p, p[19])
zeros3d(p, n) = zeros(promote_type(eltype(p), Float64), n)

f_full_on(i, p) = fcn_effective(FX3D.states[i].z, view(p, 1:15), bg3d(p), MODEL3D_ON)
f_full_off(i, p) = fcn_effective(FX3D.states[i].z, view(p, 1:15), bg3d(p), MODEL3D_OFF)
f_rho(i, p) = [matter_temperature_rate(p[15], p[16], fracs3d(p)[1], p[19], p[18])]
function f_hyd(i, p)
    xe, xp, _ = fracs3d(p); d = zeros3d(p, 6)
    dr = hydrogen_rhs!(d, FX3D.htable, p[16], p[15], xe, xp, p[19], p[17], p[18], view(p, 2:7))
    return vcat(d, dr)
end
function f_hel(i, p)
    xe, _, _ = fracs3d(p); d = zeros3d(p, 7)
    helium_base_rhs!(d, FX3D.hetable, p[16], xe, p[17], p[18], view(p, 8:14), p[19])
    return d
end
function f_abs(i, p)
    g = zeros3d(p, 15)
    hi_absorption_rhs!(g, 1, 2, 8, FX3D.dp, FX3D.bitot, FX3D.fc, FX3D.states[i].z, p[16], p[17], p[18], p[2], view(p, 8:14))
    return g
end

function cdiff3d(f, x)
    setprecision(BigFloat, 256) do
        xb = BigFloat.(x)
        cols = map(eachindex(x)) do j
            h = BigFloat(1.0e-20) * max(abs(xb[j]), BigFloat(1.0e-30)); xp = copy(xb); xm = copy(xb); xp[j] += h; xm[j] -= h
            Float64.((f(xp) .- f(xm)) ./ (2h))
        end
        return reduce(hcat, cols)
    end
end
colscale3d(J) = [max(maximum(abs, J[:, j]), floatmin()) for j in axes(J, 2)]

function check3d(name, f, x0; fdtol = 1.0e-7, nseeds = 3, varied = 2, pert = 1.0e-5, rng = Random.Xoshiro(20261001), bigfd = true)
    nout = length(f(x0))
    J = ForwardDiff.jacobian(f, x0)
    @test all(isfinite, J)
    if bigfd
        Jfd = cdiff3d(f, x0)
        # (1) AD formulas are exact: ForwardDiff run in 256-bit arithmetic equals the 256-bit central difference (error ~ h^2 ~ 1e-40)
        Jbig = setprecision(BigFloat, 256) do
            Float64.(ForwardDiff.jacobian(f, BigFloat.(x0)))
        end
        ebig = maximum(abs.(Jbig .- Jfd) ./ colscale3d(Jbig)')
        rec3d("BIGAD:" * name, ebig)
        @test ebig < 1.0e-13      # both rounded to Float64 for comparison
        # (2) Float64 ForwardDiff vs the exact derivative: bounded by the Float64 conditioning of the near-cancelling terms of the RHS (observed <= 2e-8)
        efd = maximum(abs.(J .- Jfd) ./ colscale3d(J)')
        rec3d("FD:" * name, efd)
        @test efd < fdtol
    end
    v = randn(rng, length(x0)) .* abs.(x0)
    dd = ForwardDiff.derivative(t -> f(x0 .+ t .* v), 0.0)
    rec3d("DIR:" * name, maximum(abs.(dd .- J * v)) / max(maximum(abs, J * v), floatmin()))
    @test isapprox(dd, J * v; rtol = 1.0e-11, atol = 1.0e-11 * maximum(abs, J * v))
    obj(x, w) = dot(w, f(x))
    prepA = prepare_gradient(obj, MC3D, x0, Constant(randn(rng, nout)))
    prepB = prepare_gradient(obj, MC3D, x0, Constant(randn(rng, nout)))
    for kk in 1:(varied + 1)
        x = kk == 1 ? x0 : x0 .* (1 .+ pert .* randn(rng, length(x0)))
        Jx = ForwardDiff.jacobian(f, x)
        for s in 1:nseeds
            w = randn(rng, nout)
            ref = Jx' * w
            gA = gradient(obj, prepA, MC3D, x, Constant(w))
            gB = gradient(obj, prepB, MC3D, x, Constant(w))
            sc = max(maximum(abs, ref), floatmin())
            rec3d("VJP:" * name, max(maximum(abs.(gA .- ref)) / sc, maximum(abs.(gB .- ref)) / sc))
            @test isapprox(gA, ref; rtol = 1.0e-9, atol = 1.0e-9 * sc)
            @test isapprox(gB, ref; rtol = 1.0e-9, atol = 1.0e-9 * sc)
            @test gA == gB
        end
    end
    return J
end


@testset "Chunk 3d: assembled fcn_effective AD (ForwardDiff, prepared Mooncake VJP)" begin
    @testset "components and full composition at native states" begin
        for (i, s) in enumerate(FX3D.states)
            p0 = param3d(s)
            @testset "$(s.label)" begin
                check3d("$(s.label):rho", (p -> f_rho(i, p)), p0)
                check3d("$(s.label):hydrogen", (p -> f_hyd(i, p)), p0)
                check3d("$(s.label):helium", (p -> f_hel(i, p)), p0)
                check3d("$(s.label):full_off", (p -> f_full_off(i, p)), p0)
                if s.dp_in
                    J = check3d("$(s.label):absorber", (p -> f_abs(i, p)), p0)
                    s.z <= 3400 && @test any(!iszero, J)
                    check3d("$(s.label):full_on", (p -> f_full_on(i, p)), p0)
                else
                    @test_throws DPTableDomainError f_abs(i, p0)
                    @test_throws DPTableDomainError f_full_on(i, p0)
                end
            end
        end
    end

    @testset "ForwardDiff Jacobian equals the composition of component Jacobians (absorber on = off + absorber)" begin
        for (i, s) in enumerate(FX3D.states)
            s.dp_in || continue
            p0 = param3d(s)
            Jon = ForwardDiff.jacobian(p -> f_full_on(i, p), p0)
            Joff = ForwardDiff.jacobian(p -> f_full_off(i, p), p0)
            Jab = ForwardDiff.jacobian(p -> f_abs(i, p), p0)
            @test isapprox(Jon, Joff .+ Jab; rtol = 1e-11, atol = 1e-11 * maximum(abs, Jon))
        end
    end

    @testset "X[Xe slot] has zero sensitivity; absorber inactive above zcrit has zero Jacobian" begin
        s = FX3D.states[1]; p0 = param3d(s)
        @test all(iszero, ForwardDiff.jacobian(p -> f_abs(1, p), p0))
        J = ForwardDiff.jacobian(p -> f_full_on(1, p), p0)
        @test all(iszero, J[:, 1])
    end

    @testset "observed maxima" begin
        for k in sort(collect(keys(OBS3DAD)))
            startswith(k, "FD:") && @test OBS3DAD[k] < 1e-7
            startswith(k, "BIGAD:") && @test OBS3DAD[k] < 1e-13
        end
        mx(prefix) = maximum(v for (k, v) in OBS3DAD if startswith(k, prefix))
        @info "Chunk 3d AD maxima: Float64-AD-vs-BigFloat-FD $(mx("FD:")) at $(argmax(Dict(k => v for (k, v) in OBS3DAD if startswith(k, "FD:")))), directional $(mx("DIR:")), BigFloat-AD vs BigFloat-FD $(mx("BIGAD:")), Mooncake VJP $(mx("VJP:"))"
    end
end
