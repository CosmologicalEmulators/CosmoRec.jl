using Test
using CosmoRec

include(joinpath(@__DIR__, "chunk3c_helpers.jl"))

@testset "Chunk 3c: H-I absorption of HeI photons vs original CosmoRec fixture" begin
    G, lnT, lgTg, F, cases = read_fixture3c()
    fc = fcorr_spline3c(F)
    relerr(a, b) = a == b ? 0.0 : abs(a - b) / max(abs(a), abs(b), floatmin())
    gval(key) = parse(Float64, grow3c(G, key)[2])

    @testset "fixture structure and constants" begin
        @test length(cases) == 63
        @test length(lnT) == 150 && issorted(lnT; rev = true)
        @test length(lgTg) == 500 && length(F) == 40
        @test gval("const_cl") == CosmoRec.NATIVE_HIABS_CONSTANTS.cl
        @test gval("const_h_kb") == CosmoRec.NATIVE_HIABS_CONSTANTS.h_kb
        @test gval("sig_c") == CosmoRec.NATIVE_HIABS_CONSTANTS.sig_c
        @test gval("PI") == CosmoRec.NATIVE_HIABS_CONSTANTS.pi
        for (key, line) in (("S_trans", NATIVE_HIABS_SINGLET), ("T_trans", NATIVE_HIABS_TRIPLET))
            r = grow3c(G, key)
            @test parse(Float64, r[3]) == line.A21
            @test parse(Float64, r[5]) == line.lambda21
            @test parse(Float64, r[7]) == line.Dnu
            @test line.gw == 3.0 && line.gwp == parse(Float64, r[9])
        end
        @test parse(Float64, grow3c(G, "nP_S_profile_A21")[2]) == NATIVE_HIABS_PROFILE_A21
        @test count(c -> c.inS && c.inT, cases) == 50
        @test count(c -> !c.inS && !c.inT, cases) == 9
        @test count(c -> c.inS && !c.inT, cases) == 4
        # production default: flag 1 (HI_absorption = 2 remapped), spin_forbidden 1
        @test all(haskey(c.V, 1) && haskey(c.V, 0) for c in cases)
    end

    @testset "scalar intermediates (fcorr, Bitot, pd, eta, efac, tauS)" begin
        for c in cases
            X = c.X
            @test relerr(fcorr(fc, c.Tg), c.fcorr) < 1e-12
            @test relerr(hi_abs_eta(c.NH, c.XH1s, c.Hz), c.eta) < 1e-14
            @test relerr(exp_nu(NATIVE_HIABS_SINGLET.Dnu, c.Tg), c.efacS) < 1e-14
            @test relerr(exp_nu(NATIVE_HIABS_TRIPLET.Dnu, c.Tg), c.efacT) < 1e-14
            @test relerr(tau_S_gw(NATIVE_HIABS_SINGLET.A21, NATIVE_HIABS_SINGLET.lambda21, 3.0, X[3] * c.NH, 1.0, X[1] * c.NH, c.Hz), c.tauS_S) < 1e-13
            @test relerr(tau_S_gw(NATIVE_HIABS_TRIPLET.A21, NATIVE_HIABS_TRIPLET.lambda21, 3.0, X[6] * c.NH, 1.0, X[1] * c.NH, c.Hz), c.tauS_T) < 1e-13
            if c.Bwin !== nothing
                bt = window_bitot(lgTg, c)
                @test relerr(helium_Bitot(bt, c.Tg), c.Bitot) < 1e-12
                @test relerr(pd_singlet(bt, c.Tg), c.pdS) < 1e-12
            end
            @test c.pdT == 1.0
        end
    end

    @testset "per-channel DP_interpol and dXe, both approximation flags" begin
        for c in cases, flag in (1, 0)
            v = c.V[flag]
            X = c.X
            bt = window_bitot(lgTg, c)
            if c.inS
                dp = window_dp_table(lnT, c)
                pd = pd_singlet(bt, c.Tg)
                fcv = flag == 1 ? fcorr(fc, c.Tg) : 1.0
                @test relerr(dp_correction(dp, false, c.Tg, c.tauS_S, c.eta; pd = pd, fc = fcv), v[1]) < 1e-9
                dS = hi_abs_singlet(dp, bt, fc, c.Tg, X[1], X[3], c.NH, c.Hz, c.XH1s; fcorr_on = flag == 1)
                @test relerr(dS, v[2]) < 1e-9
            else
                @test all(isnan, v)
            end
            if c.inT
                dp = window_dp_table(lnT, c)
                @test relerr(dp_correction(dp, true, c.Tg, c.tauS_T, c.eta; pd = 1.0, fc = 1.0), v[3]) < 1e-9
                dT = hi_abs_triplet(dp, c.Tg, X[1], X[6], c.NH, c.Hz, c.XH1s)
                @test relerr(dT, v[4]) < 1e-9
            else
                @test isnan(v[3]) && isnan(v[4])
            end
        end
    end

    @testset "combined derivative vector (all flag/spin/start variants)" begin
        for c in cases, ((flag, spin, kind), (aS, aT, gnat)) in c.K
            g = copy(gnat)           # placeholder to size
            g .= kind == 0 ? 0.0 : [(-1.0)^i * 1.25e-3 / 2.0^i for i in 0:8]
            g0 = copy(g)
            bt = window_bitot(lgTg, c)
            dp = window_dp_table(lnT, c)
            res = nothing
            if aS && (c.inS || dp === nothing)
                if dp === nothing
                    @test_throws Exception hi_absorption_rhs!(g, 1, 2, 3, window_dp_table(lnT, c), bt, fc, c.z, c.Tg, c.NH, c.Hz, c.XH1s, c.X)
                    continue
                end
            end
            if dp === nothing
                # nothing captured natively (outside the table or inactive): the native vector is untouched
                if aS
                    @test_throws DPTableDomainError hi_absorption_rhs!(g, 1, 2, 3, DPTable(lnT, DPSheet[]), bt, fc, c.z, c.Tg, c.NH, c.Hz, c.XH1s, c.X; spin_forbidden = spin == 1, fcorr_on = flag == 1)
                else
                    res = hi_absorption_rhs!(g, 1, 2, 3, DPTable(lnT, DPSheet[]), bt, fc, c.z, c.Tg, c.NH, c.Hz, c.XH1s, c.X; spin_forbidden = spin == 1, fcorr_on = flag == 1)
                    @test res == (0.0, 0.0)
                end
                @test g == gnat
                continue
            end
            call() = hi_absorption_rhs!(g, 1, 2, 3, dp, bt, fc, c.z, c.Tg, c.NH, c.Hz, c.XH1s, c.X; spin_forbidden = spin == 1, fcorr_on = flag == 1)
            if aT && !c.inT
                @test_throws DPTableDomainError call()      # native keeps the singlet update only; Julia accumulates it before throwing
            else
                res = call()
            end
            scale = kind == 0 ? 0.0 : 1.25e-3
            for i in 1:9
                @test abs(g[i] - gnat[i]) <= 1e-9 * abs(gnat[i] - g0[i]) + 4 * eps() * scale
            end
            if res !== nothing
                @test (res[1] != 0.0) == aS || !aS
                @test res[2] == 0.0 || aT
                # signed conservation
                @test abs(sum(g[3:9]) - sum(g0[3:9])) <= 8 * eps() * max(scale, 1e-300) + 1e-12 * (abs(res[1]) + abs(res[2]))
                # output slots: HeI 1s (+), 2^1P (-S), 2^3P1 (-T); untouched elsewhere
                tol = 8 * eps() * scale + 1e-12 * (abs(res[1]) + abs(res[2]))
                @test abs((g[3] - g0[3]) - (res[1] + res[2])) <= tol
                @test abs((g[5] - g0[5]) + res[1]) <= tol
                @test abs((g[8] - g0[8]) + res[2]) <= tol
                @test abs((g[1] - g0[1]) - (res[1] + res[2])) <= tol
                @test abs((g[2] - g0[2]) + (res[1] + res[2])) <= tol
                @test all(g[i] == g0[i] for i in (4, 6, 7, 9))
            end
        end
    end

    @testset "zcrit, diffusion and spin switches" begin
        byid(i) = cases[i + 1]
        flags(c) = (k = c.K[(1, 1, 0)]; (k[1], k[2]))
        @test flags(byid(14)) == (true, true)      # z = 3400 - 1e-6
        @test flags(byid(15)) == (true, true)      # z = 3400 exactly: native `z <= zcrit`
        @test flags(byid(16)) == (false, false)    # z = 3400 + 1e-6
        @test flags(byid(17)) == (false, false) && flags(byid(18)) == (false, false)
        @test byid(14).K[(1, 0, 0)][2] == false    # spin_forbidden = 0 disables only the intercombination channel
        c = byid(15)
        dp = window_dp_table(lnT, c); bt = window_bitot(lgTg, c)
        g = zeros(9)
        @test hi_absorption_rhs!(g, 1, 2, 3, dp, bt, fc, c.z, c.Tg, c.NH, c.Hz, c.XH1s, c.X; diffusion_correction = true) == (0.0, 0.0)
        @test all(iszero, g)
        @test hi_absorption_rhs!(g, 1, 2, 3, dp, bt, fc, 3400.0 * (1 + 4eps()), c.Tg, c.NH, c.Hz, c.XH1s, c.X) == (0.0, 0.0)
        @test all(iszero, g)
        res = hi_absorption_rhs!(g, 1, 2, 3, dp, bt, fc, c.z, c.Tg, c.NH, c.Hz, c.XH1s, c.X; spin_forbidden = false)
        @test res[1] != 0.0 && res[2] == 0.0 && g[8] == 0.0
        # accumulation (native `+=`): a second call doubles the contribution
        g2 = copy(g)
        hi_absorption_rhs!(g2, 1, 2, 3, dp, bt, fc, c.z, c.Tg, c.NH, c.Hz, c.XH1s, c.X; spin_forbidden = false)
        @test g2 ≈ 2 .* g rtol = 1e-14
        # relocated indices (iXe, iHI1s, iHeI) land in the matching slots
        g3 = zeros(14)
        hi_absorption_rhs!(g3, 4, 12, 5, dp, bt, fc, c.z, c.Tg, c.NH, c.Hz, c.XH1s, c.X; spin_forbidden = false)
        @test g3[4] == g[1] && g3[12] == g[2] && g3[5] == g[3] && g3[7] == g[5]
        @test count(!iszero, g3) == 4 && g3[12] == -g3[4] && g3[7] == -g3[5]
    end

    @testset "sobolev, tau, efac closed forms" begin
        @test sobolev_p(1e-11) ≈ 1 - 0.5e-11 atol = 1e-25
        @test sobolev_p(1e-10) == 1.0 - 0.5e-10 + 1e-20 / 6
        @test sobolev_p(2.0) == (1 - exp(-2.0)) / 2.0
        @test sobolev_p(500.0) == 1 / 500.0
        @test sobolev_p(600.0) == 1 / 600.0
        @test sobolev_p(499.999999) == (1 - exp(-499.999999)) / 499.999999
        @test exp_nu(1e30, 1.0) == exp(700.0)
        @test exp_nu(2.0e11, 2.0) == exp(CosmoRec.NATIVE_HIABS_CONSTANTS.h_kb * 1.0e11)
        # |Ni (Nj/Ni gwi/gwj - 1)| is symmetric in sign of the inversion
        a = tau_S_gw(1.0, 1.0, 3.0, 2.0, 1.0, 1.0, 1.0)
        b = tau_S_gw(1.0, 1.0, 3.0, 2.0, 1.0, 3.0, 1.0)
        @test a == (1 / (8 * pi)) * abs(2.0 * (1.0 / 2.0 * 3.0 - 1.0)) && b == (1 / (8 * pi)) * abs(2.0 * (3.0 / 2.0 * 3.0 - 1.0))
    end

    @testset "synthetic cubic table: node and off-node reproduction, edges, domain errors" begin
        nT = 12; ne = 9; nt = 10
        lnTs = collect(range(9.0, 8.0; length = nT))        # descending
        eta = collect(range(-3.0, 1.0; length = ne))
        tS = collect(range(-2.0, 4.0; length = nt))
        tT = collect(range(-1.0, 5.0; length = nt))
        fS(T, e, t) = 0.3 + 0.7T - 0.2T^2 + 0.05T^3 + e * (1.1 - 0.4t) + 0.3e^2 * t - 0.02e^3 + 0.5t^2 - 0.01t^3 + 0.1T * e * t
        fT(T, e, t) = -0.1 + 0.2T^3 - 0.3e^2 + 0.7t + 0.04e * t^2 + 0.1T^2 * e
        sheets = [DPSheet(eta, tS, tT, [fS(T - 8.5, e, t) for e in eta, t in tS], [fT(T - 8.5, e, t) for e in eta, t in tT]) for T in lnTs]
        tab = DPTable(lnTs, sheets)
        pts = [(8.0, eta[1], tS[1]), (9.0, eta[end], tS[end]), (lnTs[5], eta[3], tS[4]), (8.4321, -0.37, 1.234), (8.9999, 0.99, 3.99),
               (8.0001, -2.99, -1.99), (8.51, -1.5, 2.5), (8.77, 0.5, 0.0)]
        for (T, e, t) in pts
            @test dp_lookup(tab, false, T, e, t) ≈ fS(T - 8.5, e, t) rtol = 1e-11 atol = 1e-12
            if tT[1] <= t <= tT[end]
                @test dp_lookup(tab, true, T, e, t) ≈ fT(T - 8.5, e, t) rtol = 1e-11 atol = 1e-12
            end
        end
        @test dp_sheet_start(tab, lnTs[1]) == 1
        @test dp_sheet_start(tab, lnTs[end]) == nT - 3
        for T in range(8.0, 9.0; length = 37)
            @test 1 <= dp_sheet_start(tab, T) <= nT - 3
        end
        @test_throws DPTableDomainError dp_lookup(tab, false, 9.0 + 1e-9, 0.0, 0.0)
        @test_throws DPTableDomainError dp_lookup(tab, false, 8.0 - 1e-9, 0.0, 0.0)
        @test_throws DPTableDomainError dp_lookup(tab, false, 8.5, eta[end] + 1e-9, 0.0)
        @test_throws DPTableDomainError dp_lookup(tab, false, 8.5, eta[1] - 1e-9, 0.0)
        @test_throws DPTableDomainError dp_lookup(tab, false, 8.5, 0.0, tS[end] + 1e-9)
        @test_throws DPTableDomainError dp_lookup(tab, false, 8.5, 0.0, tS[1] - 1e-9)
        @test_throws DPTableDomainError dp_lookup(tab, true, 8.5, 0.0, tT[end] + 1e-9)
        @test_throws DPTableDomainError dp_lookup(tab, false, NaN, 0.0, 0.0)
        # Bitot series: ln B cubic in ln T is reproduced; out-of-range errors
        lg = collect(range(6.0, 9.5; length = 40)); lb = [1.0 + 0.5x - 0.02x^2 + 0.003x^3 for x in lg]
        bs = BitotSeries(lg, lb)
        for T in (exp(6.0), exp(7.123), exp(9.0), exp(lg[end - 4]))
            x = log(T)
            @test helium_Bitot(bs, T) ≈ exp(1.0 + 0.5x - 0.02x^2 + 0.003x^3) rtol = 1e-11
        end
        @test_throws DPTableDomainError helium_Bitot(bs, exp(5.999))
        @test_throws DPTableDomainError helium_Bitot(bs, exp(9.6))
        @test_throws DPTableDomainError helium_Bitot(bs, exp(lg[end] - 1e-12))      # stencil would exceed the table
        # pd limits
        @test pd_singlet(BitotSeries(lg, fill(log(1e9), 40)), exp(7.0); A21 = 1e9) ≈ 0.5 rtol = 1e-14
        @test pd_singlet(BitotSeries(lg, fill(log(1e9), 40)), exp(7.0); A21 = 1e9, f_b = 2.0) ≈ 2 / 3 rtol = 1e-14
        # fcorr: natural spline reproduces linear data exactly and clamps outside the knots
        z = collect(range(0.0, 10.0; length = 8)); sp = FcorrSpline(z, 2 .* z .+ 1)
        @test fcorr(sp, 2.725 * (1 + 3.3)) ≈ 2 * 3.3 + 1 atol = 1e-12
        @test fcorr(sp, 2.725 * 1000) == fcorr(sp, 2.725 * (1 + 10.0))
        @test fcorr(sp, 0.0) == fcorr(sp, 2.725 * (1 + 0.0 - 5.0))
    end
end
