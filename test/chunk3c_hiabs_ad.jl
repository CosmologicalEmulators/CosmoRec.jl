# Chunk 3c AD tests for the H-I absorption of HeI photons: ForwardDiff Jacobians (vs 256-bit central differences), independent invariants,
# and prepared Mooncake VJPs (DifferentiationInterface; preparation reused at varied points, plus independent preparations) vs ForwardDiff J'w,
# per channel and for the combined update. Only the right-hand-side mapping is differentiated (no ODE solve); the table values and axes are constants.
# Branches (z vs zcrit, spin switch, table cell/T-stencil seams, p_ij tau branches, native clamps) are non-smooth and are probed one-sided only.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec

@isdefined(FIX3C) || include("chunk3c_helpers.jl")

const G3C_AD, LNT3C_AD, LG3C_AD, F3C_AD, CASES3C_AD = read_fixture3c()
const FC3C_AD = fcorr_spline3c(F3C_AD)
const MC3C = AutoMooncake(; config = nothing)
case3c(label) = first(c for c in CASES3C_AD if c.label == label)

struct Tab3C
    dp::DPTable
    bt::BitotSeries
    z::Float64
    NH0::Float64
    Hz0::Float64
end
tab3c(c::Case3c) = Tab3C(window_dp_table(LNT3C_AD, c), window_bitot(LG3C_AD, c), c.z, c.NH, c.Hz)
function ctx3c(label)
    c = case3c(label)
    @assert c.inS && c.inT
    return c, tab3c(c)
end

# x = (Tg, NH, Hz, XH1s, X[7])
f_single(k, x) = (d = hi_abs_singlet(k.dp, k.bt, FC3C_AD, x[1], x[5], x[7], x[2], x[3], x[4]); [d, d * x[1]])
f_triple(k, x) = (d = hi_abs_triplet(k.dp, x[1], x[5], x[10], x[2], x[3], x[4]); [d, d * x[1]])
function f_comb(k, x; spin = true)
    g = zeros(eltype(x), 9)
    hi_absorption_rhs!(g, 1, 2, 3, k.dp, k.bt, FC3C_AD, k.z, x[1], x[2], x[3], x[4], @view(x[5:11]); spin_forbidden = spin)
    return g
end
f_comb_nospin(k, x) = f_comb(k, x; spin = false)
# test-defined background: relative scalings of the captured (NH, Hz) and Tg = Tcmb (1+z); keeps the sparse windows valid
function f_cosmo(k, x)   # Tcmb, omega_b, h, XH1s, X[7]
    Tg = x[1] * (1 + k.z); NH = k.NH0 * x[2] / 0.0224; Hz = k.Hz0 * x[3] / 0.674
    return f_comb(k, vcat(Tg, NH, Hz, x[4:end]))
end
state_x(c) = vcat(c.Tg, c.NH, c.Hz, c.XH1s, c.X)

function cdiff3c(f, x)
    setprecision(BigFloat, 256) do
        xb = BigFloat.(x)
        cols = map(eachindex(x)) do j
            h = BigFloat(1.0e-20) * max(abs(xb[j]), BigFloat(1.0e-30)); xp = copy(xb); xm = copy(xb); xp[j] += h; xm[j] -= h
            Float64.((f(xp) .- f(xm)) ./ (2h))
        end
        return reduce(hcat, cols)
    end
end
colscale3c(J) = [max(maximum(abs, J[:, j]), floatmin()) for j in axes(J, 2)]
const OBS3C = Dict{String,Float64}()
rec3c(key, v) = (OBS3C[key] = max(get(OBS3C, key, 0.0), v))

function check3c(name, f, x0; fdtol = 1.0e-8, nseeds = 3, varied = 3, pert = 1.0e-4, rng = Random.Xoshiro(20261001))
    nout = length(f(x0))
    J = ForwardDiff.jacobian(f, x0)
    Jfd = cdiff3c(f, x0)
    # columns that are identically zero (inactive inputs) must be exactly zero in both; others compared column-wise
    efd = maximum(abs.(J .- Jfd) ./ colscale3c(J)')
    rec3c("FD:" * name, efd)
    @test efd < fdtol
    @test all(isfinite, J)
    v = randn(rng, length(x0)) .* abs.(x0)
    dd = ForwardDiff.derivative(t -> f(x0 .+ t .* v), 0.0)
    rec3c("DIR:" * name, maximum(abs.(dd .- J * v)) / max(maximum(abs, J * v), floatmin()))
    @test isapprox(dd, J * v; rtol = 1.0e-12, atol = 1.0e-12 * maximum(abs, J * v))
    obj(x, w) = dot(w, f(x))
    prepA = prepare_gradient(obj, MC3C, x0, Constant(randn(rng, nout)))
    prepB = prepare_gradient(obj, MC3C, x0, Constant(randn(rng, nout)))
    for kk in 1:(varied + 1)
        x = kk == 1 ? x0 : x0 .* (1 .+ pert .* randn(rng, length(x0)))
        Jx = ForwardDiff.jacobian(f, x)
        for s in 1:nseeds
            w = randn(rng, nout)
            ref = Jx' * w
            gA = gradient(obj, prepA, MC3C, x, Constant(w))
            gB = gradient(obj, prepB, MC3C, x, Constant(w))
            sc = maximum(abs, ref)
            rec3c("VJP:" * name, max(maximum(abs.(gA .- ref)) / sc, maximum(abs.(gB .- ref)) / sc))
            @test isapprox(gA, ref; rtol = 1.0e-9, atol = 1.0e-9 * sc)
            @test isapprox(gB, ref; rtol = 1.0e-9, atol = 1.0e-9 * sc)
            @test gA == gB
        end
    end
    return J
end

@testset "Chunk 3c: H-I absorption AD (ForwardDiff Jacobians, prepared Mooncake VJPs)" begin
    cA, kA = ctx3c("Tmid20_fe0.5_ft0.5")
    cB, kB = ctx3c("Tmid148_fe0.5_ft0.5")
    cC, kC = ctx3c("Tmid100_fe0.5_ft0.5")
    xA = state_x(cA)

    @testset "per channel and combined, three temperatures" begin
        for (lab, c, k) in (("A", cA, kA), ("B", cB, kB), ("C", cC, kC))
            x = state_x(c)
            check3c("SINGLET_" * lab, y -> f_single(k, y), x)
            check3c("TRIPLET_" * lab, y -> f_triple(k, y), x)
            check3c("COMBINED_" * lab, y -> f_comb(k, y), x)
        end
        check3c("COMBINED_nospin", y -> f_comb_nospin(kA, y), xA)
    end

    @testset "structure of the combined Jacobian and conservation" begin
        J = ForwardDiff.jacobian(y -> f_comb(kA, y), xA)
        @test all(iszero, J[[4, 6, 7, 9], :])
        @test J[1, :] == J[3, :] && J[2, :] == -J[1, :]
        Js = ForwardDiff.jacobian(y -> f_comb_nospin(kA, y), xA)
        @test isapprox(J[5, :], -Js[1, :]; rtol = 1e-14)                                                 # singlet row is the spin-independent part
        @test maximum(abs.(sum(J[3:9, :]; dims = 1))) <= 1e-12 * maximum(abs, J)                         # HeI transfer: Jacobian columns sum to zero
        @test isapprox(J[1, :], -J[5, :] - J[8, :]; rtol = 1e-14)                                        # electron production = HeI 2P losses
        @test any(J[3, :] .!= 0)
        @test all(iszero, J[:, [6, 8, 9, 11]])                        # X[2], X[4], X[5], X[7] do not enter either line
        @test all(any(J[:, j] .!= 0) for j in (5, 7, 10))
        # parameters (Tg, NH, Hz, XH1s) all influence the update
        for j in 1:4
            @test any(J[:, j] .!= 0)
        end
    end

    @testset "analytic: DP == 0 table gives closed-form dependence on tau; pd = 1 gives no correction" begin
        zT = DPTable([2.0, 1.0, 0.0, -1.0, -2.0], [DPSheet(collect(-3.0:1.0:3.0), collect(-8.0:1.0:6.0), collect(-8.0:1.0:6.0), zeros(7, 15), zeros(7, 15)) for _ in 1:5])
        pf(u) = (1 - exp(-u)) / u
        dpf(u) = (u * exp(-u) - (1 - exp(-u))) / u^2
        for (pd, tau) in ((1.0, 0.5), (0.01, 0.5), (0.3, 3.0), (0.9, 20.0))
            f(t) = dp_correction(zT, false, exp(0.2), t, 1.0; pd = pd, fc = 1.0)
            Pe(t) = pd * pf(pd * t) / (1 - (1 - pd) * pf(pd * t))
            @test f(tau) ≈ Pe(tau) - pf(tau) rtol = 1e-12 atol = 1e-15
            ana = pd * (pd * dpf(pd * tau)) / (1 - (1 - pd) * pf(pd * tau))^2 - dpf(tau)
            @test ForwardDiff.derivative(f, tau) ≈ ana rtol = 1e-10 atol = 1e-14
        end
        @test dp_correction(zT, false, exp(0.2), 1.0, 1.0; pd = 1.0, fc = 1.0) == 0.0
        # derivative wrt eta through a DP linear in ln eta
        sl = DPTable([2.0, 1.0, 0.0, -1.0, -2.0], [DPSheet(collect(-3.0:1.0:3.0), collect(-8.0:1.0:6.0), collect(-8.0:1.0:6.0), [0.1e-3 * e for e in -3.0:1.0:3.0, _ in 1:15], zeros(7, 15)) for _ in 1:5])
        h(eta) = dp_correction(sl, false, exp(0.2), 1.0, eta; pd = 1.0, fc = 1.0)
        @test ForwardDiff.derivative(h, 0.7) ≈ 0.1e-3 / 0.7 rtol = 1e-10
    end

    @testset "combined update w.r.t. cosmological scalars (Tcmb, omega_b, h) via a test-defined scaling" begin
        c = cA
        x0 = vcat(c.Tg / (1 + c.z), 0.0224, 0.674, c.XH1s, c.X)
        J = check3c("COSMO", y -> f_cosmo(kA, y), x0; fdtol = 1.0e-7)
        @test any(J[:, 1] .!= 0) && any(J[:, 2] .!= 0) && any(J[:, 3] .!= 0) && any(J[:, 4] .!= 0)
    end

    @testset "z = zcrit branch, probed one-sided: active side differentiable, inactive side exactly zero" begin
        ca, ka = ctx3c("zc_below_1e-6"); cb, kb = ctx3c("zc_exact")
        Ja = ForwardDiff.jacobian(y -> f_comb(ka, y), state_x(ca))
        Jb = ForwardDiff.jacobian(y -> f_comb(kb, y), state_x(cb))
        @test any(Ja .!= 0) && any(Jb .!= 0)
        cu = case3c("zc_above_1e-6")
        ku = Tab3C(ka.dp, ka.bt, cu.z, cu.NH, cu.Hz)       # same captured windows, z just above zcrit
        Ju = ForwardDiff.jacobian(y -> f_comb(ku, y), state_x(cu))
        @test all(iszero, Ju) && all(iszero, f_comb(ku, state_x(cu)))
        # the jump across zcrit is a primal-valued branch, not a derivative; its size is recorded
        rec3c("JUMP:zcrit_value", maximum(abs, f_comb(ka, state_x(ca))))
        objz(y, w) = dot(w, f_comb(ku, y))
        prz = prepare_gradient(objz, MC3C, state_x(cu), Constant(ones(9)))
        @test all(iszero, gradient(objz, prz, MC3C, state_x(cu), Constant(ones(9))))
    end

    @testset "T-stencil seam probed one-sided (derivative jumps across a sheet node)" begin
        wfix = randn(Random.Xoshiro(11), 9)
        Jside = Dict{String,Vector{Float64}}()
        for lab in ("Tnode60p_fe0.5_ft0.5", "Tnode60m_fe0.5_ft0.5")
            c, k = ctx3c(lab); x = state_x(c)
            J = ForwardDiff.jacobian(y -> f_comb(k, y), x)
            Jfd = cdiff3c(y -> f_comb(k, y), x)
            ef = maximum(abs.(J[:, 1] .- Jfd[:, 1])) / max(maximum(abs, J[:, 1]), floatmin())
            rec3c("FD:SEAM_Tg", ef)
            @test ef < 1.0e-6
            obj(y, w) = dot(w, f_comb(k, y))
            prep = prepare_gradient(obj, MC3C, x, Constant(wfix))
            g = gradient(obj, prep, MC3C, x, Constant(wfix))
            rec3c("VJP:SEAM", maximum(abs.(g .- J' * wfix)) / maximum(abs, J' * wfix))
            @test isapprox(g, J' * wfix; rtol = 1.0e-9, atol = 1.0e-9 * maximum(abs, J' * wfix))
            Jside[lab] = J[1, :]
        end
        rec3c("JUMP:T_seam_rel", maximum(abs.(Jside["Tnode60p_fe0.5_ft0.5"] .- Jside["Tnode60m_fe0.5_ft0.5"])) / maximum(abs, Jside["Tnode60m_fe0.5_ft0.5"]))
    end

    @testset "known limitation: lookup domain (native would run the explicit integral)" begin
        x = state_x(cA)
        @test_throws DPTableDomainError f_comb(kA, vcat(x[1], x[2] * 1e6, x[3:end]))       # eta (via NH, tau) outside the table
        @test_throws DPTableDomainError f_comb(kA, vcat(x[1] * 0.1, x[2:end]))             # T below the sheets
        @test_throws DPTableDomainError f_comb(kA, vcat(x[1:3], x[4] * 1e9, x[5:end]))     # eta high
    end

    println("  chunk3c AD observed worst errors:")
    for k in sort(collect(keys(OBS3C)))
        println("    ", rpad(k, 22), OBS3C[k])
    end
end
