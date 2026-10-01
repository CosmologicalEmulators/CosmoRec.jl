# Chunk 3b AD tests for the helium base RHS components: ForwardDiff Jacobians (vs 256-bit central differences) and prepared Mooncake VJPs
# (DifferentiationInterface; preparation reused at varied points, plus independent preparations) vs ForwardDiff J'w, per component and for the combined RHS.
# These differentiate the RIGHT-HAND-SIDE MAPPING only at smooth interior points (and one-sided at stencil seams); no ODE solve is differentiated.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec

@isdefined(FIX3B) || include("chunk3b_helpers.jl")

const G3B_AD, LG3B_AD, CASES3B_AD, W3B_AD = read_fixture3b()
const TABH = window_table(LG3B_AD, W3B_AD)
const CH0 = first(c for c in CASES3B_AD if c.label == "phys_z2300")
const MCH = AutoMooncake(; config = nothing)
const ATH = NATIVE_HELIUM_ATOM

zeros_like(x, n) = zeros(promote_type(eltype(x), Float64), n)
upperH(R) = [R[i, j] for i in 1:4 for j in (i + 1):4]
function unupperH(v, T)
    R = zeros(T, 4, 4); k = 0
    for i in 1:4, j in (i + 1):4
        k += 1; R[i, j] = v[k]
    end
    return R
end

f_rci(x) = (d = zeros_like(x, 7); helium_continuum!(d, x[1], x[2], @view(x[3:9]), @view(x[10:13]), @view(x[14:17])); d)   # Xe, Nc, X(7), A(4), B(4)
f_rij(x) = (d = zeros_like(x, 7); helium_interlevel!(d, x[1], @view(x[2:8]), unupperH(@view(x[9:14]), eltype(x))); d)         # Tg, X(7), Rupper(6)
f_2ph(x) = (d = zeros_like(x, 7); helium_two_photon!(d, x[1], @view x[2:8]); d)                                              # Tg, X(7)
f_ly(x) = (d = zeros_like(x, 7); helium_lyman!(d, x[1], @view(x[4:10]), x[2], x[3]); d)                                     # Tg, NH, Hz, X(7)
function f_all_rates(x)   # Tg, Xe, NH, Hz, fHe, X(7), A(4), B(4), Rupper(6)
    d = zeros_like(x, 7)
    helium_base_rhs!(d, x[1], x[2], x[3], x[4], @view(x[6:12]), x[5], @view(x[13:16]), @view(x[17:20]), unupperH(@view(x[21:26]), eltype(x)))
    return d
end
function f_table(x)       # Tg, Xe, NH, Hz, fHe, X(7)
    d = zeros_like(x, 7)
    helium_base_rhs!(d, TABH, x[1], x[2], x[3], x[4], @view(x[6:12]), x[5])
    return d
end
const ZH = CH0.Tg / 2.7255 - 1
function bgH(theta, z)
    Tcmb, wb, Om, h = theta
    Orad = 2.4728e-5 * (Tcmb / 2.7255)^4 / h^2 * (1 + 0.2271 * 3.046)
    Hz = 100 * h * 1.0e5 / 3.0856775814913673e24 * sqrt(Om * (1 + z)^3 + Orad * (1 + z)^4 + 1 - Om - Orad)
    NH = 1.886e-7 * (wb / 0.0224) * (1 + z)^3 * (1 - 0.2454) / 0.76
    return Tcmb * (1 + z), NH, Hz
end
function f_cosmo(x)       # theta(4), Xe, fHe, X(7)
    d = zeros_like(x, 7)
    Tg, NH, Hz = bgH(@view(x[1:4]), ZH)
    helium_base_rhs!(d, TABH, Tg, x[5], NH, Hz, @view(x[7:13]), x[6])
    return d
end

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
const OBSH = Dict{String,Float64}()
recordH(key, v) = (OBSH[key] = max(get(OBSH, key, 0.0), v))

function check_componentH(name, f, x0; fdtol = 1.0e-12, nseeds = 3, varied = 3, pert = 1.0e-3, rng = Random.Xoshiro(20260930))
    nout = length(f(x0))
    J = ForwardDiff.jacobian(f, x0)
    Jfd = cdiff(f, x0)
    efd = maximum(abs.(J .- Jfd) ./ colscale(J)')
    recordH("FD:" * name, efd)
    @test efd < fdtol
    @test all(isfinite, J)
    v = randn(rng, length(x0)) .* abs.(x0)
    dd = ForwardDiff.derivative(t -> f(x0 .+ t .* v), 0.0)
    recordH("DIR:" * name, maximum(abs.(dd .- J * v)) / max(maximum(abs, J * v), floatmin()))
    @test isapprox(dd, J * v; rtol = 1.0e-12, atol = 1.0e-12 * maximum(abs, J * v))
    obj(x, w) = dot(w, f(x))
    prepA = prepare_gradient(obj, MCH, x0, Constant(randn(rng, nout)))
    prepB = prepare_gradient(obj, MCH, x0, Constant(randn(rng, nout)))
    for k in 1:(varied + 1)
        x = k == 1 ? x0 : x0 .* (1 .+ pert .* randn(rng, length(x0)))
        Jx = ForwardDiff.jacobian(f, x)
        for s in 1:nseeds
            w = randn(rng, nout)
            ref = Jx' * w
            gA = gradient(obj, prepA, MCH, x, Constant(w))
            gB = gradient(obj, prepB, MCH, x, Constant(w))
            sc = maximum(abs, ref)
            recordH("VJP:" * name, max(maximum(abs.(gA .- ref)) / sc, maximum(abs.(gB .- ref)) / sc))
            @test isapprox(gA, ref; rtol = 1.0e-9, atol = 1.0e-9 * sc)
            @test isapprox(gB, ref; rtol = 1.0e-9, atol = 1.0e-9 * sc)
            @test gA == gB
        end
    end
    return J
end

@testset "Chunk 3b: helium RHS AD (ForwardDiff Jacobians, prepared Mooncake VJPs)" begin
    c = CH0
    atX = c.X .* [1.0, 1.1, 0.9, 1.2, 0.8, 1.05, 0.95]
    Nc = c.NH * c.XHeII
    @testset "two-photon 2^1S-1^1S" begin
        J = check_componentH("2PH", f_2ph, vcat(c.Tg, atX))
        @test J[1, :] == -J[2, :]
        # analytic: d/dX2 = A2s1s, d/dX1 = -A2s1s exp(-h_kb Dnu/Tg)
        @test J[1, 3] == 51.3 && isapprox(J[1, 2], -51.3 * exp(-4.7992373449498863e-11 * 4984872582556676.0 / c.Tg); rtol = 1e-14)
    end
    @testset "Lyman / intercombination Sobolev channels" begin
        # native p_ij = (1 - exp(-tau))/tau is evaluated without expm1; for the 2^3S line tau ~ 1.3e-7, so Float64 (native and Julia alike)
        # carries a relative roundoff ~ eps/tau ~ 2e-9 in p and its derivatives: the 256-bit reference therefore differs at 1e-11, not at 1e-12.
        check_componentH("LY", f_ly, vcat(c.Tg, c.NH, c.Hz, atX); fdtol = 1.0e-10)
    end
    @testset "effective continuum Rci/Ric" begin
        J = check_componentH("RCI", f_rci, vcat(c.Xe, Nc, atX, c.A, c.B))
        # analytic: d dX_i / dX_i = -B_m, d dX_i / dA_m = B_m Xe Nc, d dX_i/dXe = B A Nc
        for m in 1:4
            i = ATH.res_index[m]
            @test J[i, 2 + i] == -c.B[m]
            @test isapprox(J[i, 9 + m], c.B[m] * c.Xe * Nc; rtol = 1e-14)
            @test isapprox(J[i, 13 + m], c.A[m] * c.Xe * Nc - atX[i]; rtol = 1e-12)
        end
        @test all(J[[1, 5, 7], :] .== 0)
    end
    @testset "interlevel Rij" begin
        J = check_componentH("RIJ", f_rij, vcat(c.Tg, atX, upperH(c.R)))
        @test all(abs.(sum(J; dims = 1)) .<= 1.0e-10 * sum(abs, J; dims = 1) .+ 1e-300)    # transfer: Jacobian columns sum to zero
    end
    @testset "combined RHS from fixed rates (all scalar, state and rate inputs)" begin
        check_componentH("ALL_RATES", f_all_rates, vcat(c.Tg, c.Xe, c.NH, c.Hz, c.fHe, atX, c.A, c.B, upperH(c.R)))
    end
    @testset "combined RHS through the explicit helium table (rate coordinate Tg)" begin
        # Tg perturbations are kept inside one native stencil interval (the sparse test table has finite rows only at captured windows)
        x0 = vcat(c.Tg, c.Xe, c.NH, c.Hz, c.fHe, atX)
        @test CosmoRec.helium_stencil_start(TABH, log(c.Tg * (1 + 1e-3))) == CosmoRec.helium_stencil_start(TABH, log(c.Tg * (1 - 1e-3)))
        J = check_componentH("TABLE", f_table, x0; fdtol = 1.0e-10, pert = 1.0e-4)
        @test any(J[:, 1] .!= 0)
    end
    @testset "combined RHS w.r.t. cosmological scalars (Tcmb, omega_b, Omega_m, h) via a test-defined background" begin
        x0 = vcat(c.Tg / (1 + ZH), 0.0224, 0.315, 0.674, c.Xe, c.fHe, atX)
        J = check_componentH("COSMO", f_cosmo, x0; fdtol = 1.0e-10, pert = 1.0e-4)
        @test any(J[:, 1] .!= 0) && any(J[:, 2] .!= 0) && any(J[:, 3] .!= 0) && any(J[:, 4] .!= 0)
    end
    @testset "stencil seams probed one-sided only (derivative jumps across a node)" begin
        wfix = randn(Random.Xoshiro(7), 7)
        # Points are 1e-6 (relative in Tg) from the node: closer than that, the Float64 rounding of log(Tg) (~1e-15) is a 1e-9 relative error in (x - x_node).
        seams = [("node100", 100, 1 - 1e-6), ("node100", 100, 1 + 1e-6), ("node250", 250, 1 - 1e-6), ("node250", 250, 1 + 1e-6), ("edge_low", 0, 1 + 1e-6), ("top_lx496", 497, 1 - 1e-6)]
        for (lab, k, f) in seams
            cc = first(q for q in CASES3B_AD if q.label == (lab == "edge_low" ? "edge_low_plus" : lab == "top_lx496" ? "top_lx496_minus" : lab * "_exact"))
            @test (CosmoRec.helium_stencil_start(TABH, log(exp(LG3B_AD[k + 1]) * f)) - 1) in (k - 1, k, k + 1, 496)
            x = vcat(exp(LG3B_AD[k + 1]) * f, cc.Xe, cc.NH, cc.Hz, cc.fHe, cc.X .* [1.0, 1.1, 0.9, 1.2, 0.8, 1.05, 0.95])
            J = ForwardDiff.jacobian(f_table, x)
            @test all(isfinite, J)
            # one-sided 256-bit difference in Tg along the side of the stencil (1e-9 from the node: far larger than h = 1e-20)
            Jfd = cdiff(f_table, x)
            # d(log rate)/d logTg = sum_k a_k' f_k with |f_k| up to ~70 and a_k' ~ 1/dlogTg ~ 175: Float64 roundoff ~ 70*175*eps ~ 1e-12 of the rate,
            # which is a 1e-7 relative error where the rate derivative is small (low-Tg end); the 256-bit reference uses the same Float64 table values.
            ef = maximum(abs.(J[:, 1] .- Jfd[:, 1])) / max(maximum(abs, J[:, 1]), floatmin())
            recordH("FD:SEAM_Tg", ef)
            @test ef < 1.0e-6
            obj(x, w) = dot(w, f_table(x))
            prep = prepare_gradient(obj, MCH, x, Constant(wfix))
            g = gradient(obj, prep, MCH, x, Constant(wfix))
            recordH("VJP:SEAM", maximum(abs.(g .- J' * wfix)) / maximum(abs, J' * wfix))
            @test isapprox(g, J' * wfix; rtol = 1.0e-9, atol = 1.0e-9 * maximum(abs, J' * wfix))
        end
        # exactly at a node the primal selects the upper stencil (largest j with node <= x): the Dual Jacobian equals the right-sided one
        k = 250
        cc = first(q for q in CASES3B_AD if q.label == "node$(k)_exact")
        x = vcat(exp(LG3B_AD[k + 1]), cc.Xe, cc.NH, cc.Hz, cc.fHe, cc.X .* [1.0, 1.1, 0.9, 1.2, 0.8, 1.05, 0.95])
        @test CosmoRec.helium_stencil_start(TABH, log(x[1])) in (k, k + 1)
        Jn = ForwardDiff.jacobian(f_table, x)
        xp = copy(x); xp[1] *= 1 + 1e-9
        Jp = ForwardDiff.jacobian(f_table, xp)
        @test isapprox(Jn[:, 1], Jp[:, 1]; rtol = 1e-4) || isapprox(Jn[:, 1], ForwardDiff.jacobian(f_table, x .* [1 - 1e-9; ones(11)])[:, 1]; rtol = 1e-4)
    end
    @testset "known limitation: native lookup domain" begin
        @test_throws RateTableDomainError f_table(vcat(exp(LG3B_AD[1]) * (1 - 1e-9), c.Xe, c.NH, c.Hz, c.fHe, atX))
        @test_throws RateTableDomainError f_table(vcat(exp(LG3B_AD[end]), c.Xe, c.NH, c.Hz, c.fHe, atX))
    end
    println("  chunk3b AD observed worst errors:")
    for k in sort(collect(keys(OBSH)))
        println("    ", rpad(k, 22), OBSH[k])
    end
end
