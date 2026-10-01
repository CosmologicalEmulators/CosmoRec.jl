# Chunk 3a AD tests for the hydrogen RHS components: ForwardDiff Jacobians (vs central finite differences) and prepared Mooncake VJPs
# (DifferentiationInterface; preparation reused at varied points, plus independent preparations) vs ForwardDiff J'w, per component and for the combined RHS.
# These differentiate the RIGHT-HAND-SIDE MAPPING only at smooth interior points; no ODE solve is differentiated.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec

@isdefined(FIXDIR2) || include("chunk2_helpers.jl")
@isdefined(FIX3A) || include("chunk3a_helpers.jl")

const TBL3 = first(read_table_window2())
const G3, CASES3 = read_fixture3a()
const C0 = first(c for c in CASES3 if c.label == "tbl_int_a")      # Tg = 3000, Te/Tg = 0.9: interior of the accepted table window
const MC3 = AutoMooncake(; config = nothing)

zeros_like(x, n) = zeros(promote_type(eltype(x), Float64), n)
upper(R) = [R[i, j] for i in 1:5 for j in (i + 1):5]
function unupper(v, T)
    R = zeros(T, 5, 5); k = 0
    for i in 1:5, j in (i + 1):5
        k += 1; R[i, j] = v[k]
    end
    return R
end

# argument-vector mappings; each returns a fresh output vector (promoted element type)
f_tm(x) = [matter_temperature_rate(x[1], x[2], x[3], x[4], x[5])]                                           # rho, Tg, Xe, fHe, Hz
f_2ph(x) = (d = zeros_like(x, 6); two_photon!(d, x[1], @view x[2:7]); d)                                       # Tg, X(6)
f_ly(x) = (d = zeros_like(x, 6); lyman!(d, x[1], @view(x[4:9]), x[2], x[3]); d)                                # Tg, NH, Hz, X(6)
f_rci(x) = (d = zeros_like(x, 6); continuum!(d, x[1], x[2], @view(x[3:8]), @view(x[9:13]), @view(x[14:18])); d)  # Xe, Np, X(6), A(5), B(5)
f_rij(x) = (d = zeros_like(x, 6); interlevel!(d, x[1], @view(x[2:7]), unupper(@view(x[8:17]), eltype(x))); d)   # Tg, X(6), Rupper(10)
function f_all_rates(x)   # Tg, Xe, Xp, NH, Hz, X(6), A(5), B(5), Rupper(10)
    d = zeros_like(x, 6)
    hydrogen_rhs!(d, x[1], x[2], x[3], x[4], x[5], @view(x[6:11]), @view(x[12:16]), @view(x[17:21]), unupper(@view(x[22:31]), eltype(x)))
    return d
end
function f_table(x)   # Tg, rho, Xp, fHe, NH, Hz, X(6)   (Xe = Xp)
    d = zeros_like(x, 6)
    dr = hydrogen_rhs!(d, TBL3, x[1], x[2], x[3], x[3], x[4], x[5], x[6], @view x[7:12])
    return vcat(d, dr)
end

# test-defined cosmological background (the background module is not part of Chunk 3a): theta = (Tcmb, omega_b, Omega_m, h); z fixed
const Z3 = C0.Tg / 2.7255 - 1
function bg(theta, z)
    Tcmb, wb, Om, h = theta
    Orad = 2.4728e-5 * (Tcmb / 2.7255)^4 / h^2 * (1 + 0.2271 * 3.046)
    Hz = 100 * h * 1.0e5 / 3.0856775814913673e24 * sqrt(Om * (1 + z)^3 + Orad * (1 + z)^4 + 1 - Om - Orad)
    NH = 1.886e-7 * (wb / 0.0224) * (1 + z)^3 * (1 - 0.2454) / 0.76
    return Tcmb * (1 + z), NH, Hz
end
# theta(4), X(6), rho, Xp, fHe
function f_cosmo(x)
    d = zeros_like(x, 6)
    Tg, NH, Hz = bg(@view(x[1:4]), Z3)
    dr = hydrogen_rhs!(d, TBL3, Tg, x[11], x[12], x[12], x[13], NH, Hz, @view x[5:10])
    return vcat(d, dr)
end

# 256-bit central differences (h = 1e-20 |x_j|): removes the Float64 cancellation that makes plain finite differences useless here
function cdiff(f, x)
    setprecision(BigFloat, 256) do
        xb = BigFloat.(x)
        cols = map(eachindex(x)) do j
            h = BigFloat(1.0e-20) * max(abs(xb[j]), BigFloat(1.0e-30)); xp = copy(xb); xm = copy(xb); xp[j] += h; xm[j] -= h
            Float64.((f(xp) .- f(xm)) ./ (2h))
        end
        return reduce(hcat, cols)
    end
end

colscale(J) = [max(maximum(abs, J[:, j]), floatmin()) for j in axes(J, 2)]
const OBS = Dict{String,Float64}()
function record(key, v)
    OBS[key] = max(get(OBS, key, 0.0), v)
end

function check_component(name, f, x0; fdtol = 1.0e-12, nseeds = 3, varied = 3, rng = Random.Xoshiro(20260930))
    nout = length(f(x0))
    J = ForwardDiff.jacobian(f, x0)
    Jfd = cdiff(f, x0)
    # BigFloat central-difference check, each column scaled by its largest entry
    efd = maximum(abs.(J .- Jfd) ./ colscale(J)')
    record("FD:" * name, efd)
    @test efd < fdtol
    @test all(isfinite, J)
    # directional derivative: ForwardDiff with a seeded direction equals J * v
    v = randn(rng, length(x0)) .* abs.(x0)
    dd = ForwardDiff.derivative(t -> f(x0 .+ t .* v), 0.0)
    record("DIR:" * name, maximum(abs.(dd .- J * v)) / max(maximum(abs, J * v), floatmin()))
    @test isapprox(dd, J * v; rtol = 1.0e-12, atol = 1.0e-12 * maximum(abs, J * v))
    # prepared Mooncake VJP vs J'w, preparation reused at varied points; and two independent preparations
    obj(x, w) = dot(w, f(x))
    w1 = randn(rng, nout)
    prepA = prepare_gradient(obj, MC3, x0, Constant(w1))
    prepB = prepare_gradient(obj, MC3, x0, Constant(randn(rng, nout)))
    for k in 1:(varied + 1)
        x = k == 1 ? x0 : x0 .* (1 .+ 1.0e-3 .* randn(rng, length(x0)))
        Jx = ForwardDiff.jacobian(f, x)
        for s in 1:nseeds
            w = randn(rng, nout)
            ref = Jx' * w
            gA = gradient(obj, prepA, MC3, x, Constant(w))
            gB = gradient(obj, prepB, MC3, x, Constant(w))
            sc = maximum(abs, ref)
            record("VJP:" * name, max(maximum(abs.(gA .- ref)) / sc, maximum(abs.(gB .- ref)) / sc))
            @test isapprox(gA, ref; rtol = 1.0e-9, atol = 1.0e-9 * sc)
            @test isapprox(gB, ref; rtol = 1.0e-9, atol = 1.0e-9 * sc)
            @test gA == gB   # independent preparations of the same objective give identical VJPs
        end
    end
    return J
end

@testset "Chunk 3a: RHS AD (ForwardDiff Jacobians, prepared Mooncake VJPs)" begin
    c = C0
    atX = c.X .* [1.0, 1.1, 0.9, 1.2, 0.8, 1.05]   # off detailed-balance equilibrium so no output cancels exactly
    @testset "Compton / matter-temperature equation" begin
        check_component("TM", f_tm, [c.rho, c.Tg, c.Xe, c.fHe, c.Hz])
    end
    @testset "two-photon 2s-1s" begin
        J = check_component("2PH", f_2ph, vcat(c.Tg, atX))
        @test J[1, :] == -J[2, :]        # opposite transfer
    end
    @testset "Lyman Sobolev channels" begin
        check_component("LY", f_ly, vcat(c.Tg, c.NH, c.Hz, atX))
    end
    @testset "effective continuum Rci/Ric" begin
        check_component("RCI", f_rci, vcat(c.Xe, c.NH * c.Xp, atX, c.A, c.B))
    end
    @testset "interlevel Rij" begin
        J = check_component("RIJ", f_rij, vcat(c.Tg, atX, upper(c.R)))
        @test all(abs.(sum(J[2:6, :]; dims = 1)) .<= 1.0e-10 * sum(abs, J[2:6, :]; dims = 1) .+ 1e-300)   # conservation holds for the Jacobian too
    end
    @testset "combined RHS from fixed rates (all scalar, state and rate inputs)" begin
        J = check_component("ALL_RATES", f_all_rates, vcat(c.Tg, c.Xe, c.Xp, c.NH, c.Hz, atX, c.A, c.B, upper(c.R)))
    end
    @testset "combined RHS through the explicit rate table (rate-query coordinates Tg, Te = rho Tg)" begin
        check_component("TABLE", f_table, vcat(c.Tg, c.rho, c.Xp, c.fHe, c.NH, c.Hz, atX); fdtol = 1.0e-10)
    end
    @testset "combined RHS w.r.t. cosmological scalars (Tcmb, omega_b, Omega_m, h) via a test-defined background" begin
        x0 = vcat(2.7255 * c.Tg / c.Tg, 0.0224, 0.315, 0.674, atX, c.rho, c.Xp, c.fHe)
        x0[1] = c.Tg / (1 + Z3)
        J = check_component("COSMO", f_cosmo, x0; fdtol = 1.0e-10)
        @test all(J[:, 1] .!= 0) || any(J[:, 1] .!= 0)
        @test any(J[1:6, 2] .!= 0) && any(J[1:6, 3] .!= 0) && any(J[1:6, 4] .!= 0)
    end
    @testset "documented rate switches probed one-sided only (detailed-balance switch)" begin
        # rates are piecewise at Tg/2.725 - 1 = 2000 (Tg = 5452.725...): derivatives are compared on each side separately, never across
        cb = first(cc for cc in CASES3 if cc.label == "tbl_db_last_below"); ca = first(cc for cc in CASES3 if cc.label == "tbl_db_above")
        for cc in (cb, ca)
            x = vcat(cc.Tg, cc.rho, cc.Xp, cc.fHe, cc.NH, cc.Hz, cc.X)
            J = ForwardDiff.jacobian(f_table, x)
            @test all(isfinite, J)
            w = randn(Random.Xoshiro(7), 7)
            obj(x, w) = dot(w, f_table(x))
            prep = prepare_gradient(obj, MC3, x, Constant(w))
            g = gradient(obj, prep, MC3, x, Constant(w))
            @test isapprox(g, J' * w; rtol = 1.0e-9, atol = 1.0e-9 * maximum(abs, J' * w))
            record("VJP:DB_one_sided", maximum(abs.(g .- J' * w)) / maximum(abs, J' * w))
        end
    end
    @testset "known limitation: RHS and gradients below Tg ~ 55.6 K (native A overflow)" begin
        x = vcat(31.0, 0.5, 0.5, 0.0817, 1.0e-3, 1.0e-13, 0.9, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12)
        y = f_table(x)
        @test isnan(y[2]) && isnan(y[3]) && all(isfinite, y[[1, 4, 5, 6, 7]])   # only the n = 2 rates A are NaN (native A overflow); others finite
        @test !all(isfinite, ForwardDiff.jacobian(f_table, x))                    # ForwardDiff: NaN/Inf entries
        obj(x, w) = dot(w, f_table(x))
        e1 = [1.0; zeros(6)]                                                        # weight only on the finite 1s output
        prep = prepare_gradient(obj, MC3, x, Constant(e1))
        g = gradient(obj, prep, MC3, x, Constant(e1))
        @test any(isnan, g)                                                        # Mooncake: NaN leaks into rate-coordinate/background gradients even for a finite output
    end
    println("  chunk3a AD observed worst errors:")
    for k in sort(collect(keys(OBS)))
        println("    ", rpad(k, 22), OBS[k])
    end
end
