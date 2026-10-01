# Chunk 3a: pure-Julia hydrogen RHS components vs ORIGINAL CosmoRec function outputs (text fixture), plus independent algebraic invariants.
# Fixture: test/fixtures/native_hrhs_components.txt (provenance, NOTICE and hashes in its header). Native sentinel "before" buffers are used
# so accumulation (+=) versus assignment (TM) is verified.
using Test
using CosmoRec

@isdefined(FIXDIR2) || include("chunk2_helpers.jl")
@isdefined(FIX3A) || include("chunk3a_helpers.jl")

const G3A, CASES3A = read_fixture3a()
const TOL3A = 1e-13   # observed errors are recorded in docs/CHUNK3A_RESULTS.md

relerr3a(x, y) = maximum(abs.(x .- y) ./ max.(abs.(y), floatmin(Float64)))

@testset "Chunk 3a: native hydrogen RHS components" begin
    @testset "fixture provenance and atomic data" begin
        h = [l for l in eachline(FIX3A) if startswith(l, "#")]
        @test any(startswith("# cosmorec_sha: 086769055f61ae0c244a53dd381ee65b624d0ac3"), h)
        @test any(contains("Chluba & Thomas 2010"), h) && any(contains("Chluba & Sunyaev 2006"), h)
        @test any(contains("static GSL 2.8"), h)
        @test length(CASES3A) == 16
        @test all(length(c.X) == 6 && haskey(c.after, "FCN_HI") && haskey(c.after, "TM") for c in CASES3A)
        K = NATIVE_RHS_CONSTANTS
        val(key) = parse(Float64, grow(G3A, key)[2])
        @test val("PI") == K.pi && val("const_sigT") == K.sigT && val("const_cl") == K.c && val("const_kB") == K.kB
        @test val("const_hbar") == K.hbar && val("const_me_gr") == K.me_gr && val("const_h_kb") == K.h_kb && val("const_HI_A2s_1s") == K.A2s1s
        at = NATIVE_HYDROGEN_ATOM
        @test parse(Float64, grow(G3A, "Dnu_1s_level1")[2]) == at.nu21_2s
        res = [r for r in G3A if r[1] == "resolved"]
        @test [parse(Float64, r[end]) for r in res] == at.gw
        for (k, r) in enumerate(r for r in G3A if r[1] == "Ly")
            @test parse(Int, r[3]) == at.ly_n[k] && parse(Int, r[5]) + 1 == at.ly_index[k]
            @test parse(Float64, r[7]) == at.ly_A21[k] && parse(Float64, r[9]) == at.ly_lambda21[k] && parse(Float64, r[11]) == at.ly_nu21[k]
        end
        for r in G3A
            r[1] == "nu_ul" || continue
            i, j = parse(Int, r[2]) + 1, parse(Int, r[3]) + 1
            @test parse(Float64, r[4]) == at.nu_ul[i, j]
        end
        @test [(r[6], r[8]) for r in res] == [("2", "0"), ("2", "1"), ("3", "0"), ("3", "1"), ("3", "2")]
    end

    errs = Dict{String,Float64}()
    upd(k, e) = (errs[k] = max(get(errs, k, 0.0), e))

    @testset "components at every captured state" begin
        for c in CASES3A
            K = NATIVE_RHS_CONSTANTS
            tm = matter_temperature_rate(c.rho, c.Tg, c.Xe, c.fHe, c.Hz)
            @test isapprox(tm, c.after["TM"][1]; rtol = TOL3A)
            upd("TM", relerr3a([tm], c.after["TM"]))

            for (name, f!) in (
                ("2PH", d -> two_photon!(d, c.Tg, c.X)),
                ("LY2", d -> lyman_channel!(d, 3, c.Tg, c.X[1], c.X[3], c.NH, c.Hz, NATIVE_HYDROGEN_ATOM.ly_A21[1], NATIVE_HYDROGEN_ATOM.ly_lambda21[1], NATIVE_HYDROGEN_ATOM.ly_nu21[1])),
                ("LY3", d -> lyman_channel!(d, 5, c.Tg, c.X[1], c.X[5], c.NH, c.Hz, NATIVE_HYDROGEN_ATOM.ly_A21[2], NATIVE_HYDROGEN_ATOM.ly_lambda21[2], NATIVE_HYDROGEN_ATOM.ly_nu21[2])),
                ("RCI", d -> continuum!(d, c.Xe, c.NH * c.Xp, c.X, c.A, c.B)),
                ("RIJ", d -> interlevel!(d, c.Tg, c.X, c.R)),
                ("FCN_HI", d -> hydrogen_rhs!(d, c.Tg, c.Xe, c.Xp, c.NH, c.Hz, c.X, c.A, c.B, c.R)),
            )
                d = copy(c.before[name])
                f!(d)
                @test all(isapprox.(d, c.after[name]; rtol = TOL3A, atol = 0))
                upd(name, relerr3a(d, c.after[name]))
            end
            # native combined function equals the sum of the native pieces, applied sequentially to the same buffer
            d = copy(c.before["FCN_HI"])
            two_photon!(d, c.Tg, c.X); lyman!(d, c.Tg, c.X, c.NH, c.Hz); continuum!(d, c.Xe, c.NH * c.Xp, c.X, c.A, c.B); interlevel!(d, c.Tg, c.X, c.R)
            @test all(isapprox.(d, c.after["FCN_HI"]; rtol = TOL3A, atol = 0))
        end
        for k in sort(collect(keys(errs)))
            println("  chunk3a native max rel err ", rpad(k, 7), errs[k])
        end
    end

    @testset "table-driven RHS reproduces the native combined call at accepted-window queries" begin
        tbl, _, _ = read_table_window2()
        for c in filter(c -> startswith(c.label, "tbl_"), CASES3A)
            r = get_rates(tbl, c.Tg, c.rho * c.Tg)
            @test isapprox(r.A, c.A; rtol = 1e-13) && isapprox(r.B, c.B; rtol = 1e-13) && isapprox(r.R, c.R; rtol = 1e-13)
            d = copy(c.before["FCN_HI"])
            drho = hydrogen_rhs!(d, tbl, c.Tg, c.rho, c.Xe, c.Xp, c.fHe, c.NH, c.Hz, c.X)
            @test all(isapprox.(d, c.after["FCN_HI"]; rtol = 1e-13, atol = 0))
            @test isapprox(drho, c.after["TM"][1]; rtol = 1e-13)
            println("  chunk3a table-driven ", rpad(c.label, 18), "max rel err ", relerr3a(d, c.after["FCN_HI"]))
        end
    end
end

@testset "Chunk 3a: algebraic invariants" begin
    c0 = CASES3A[1]
    K = NATIVE_RHS_CONSTANTS
    @testset "transfers are opposite and conserve population" begin
        for c in CASES3A
            d = zeros(6); two_photon!(d, c.Tg, c.X); @test d[1] == -d[2] && d[3:6] == zeros(4)
            d = zeros(6); lyman!(d, c.Tg, c.X, c.NH, c.Hz); @test d[1] + d[3] + d[5] == 0 || abs(d[1] + d[3] + d[5]) <= 4eps() * (abs(d[1]) + abs(d[3]) + abs(d[5]))
            @test d[2] == 0 && d[4] == 0 && d[6] == 0
            d = zeros(6); interlevel!(d, c.Tg, c.X, c.R); @test abs(sum(d)) <= 8eps() * sum(abs, d) && d[1] == 0
            d = zeros(6); hydrogen_rhs!(d, c.Tg, c.Xe, c.Xp, c.NH, c.Hz, c.X, c.A, c.B, c.R)
            dc = zeros(6); continuum!(dc, c.Xe, c.NH * c.Xp, c.X, c.A, c.B)
            @test abs(sum(d) - sum(dc)) <= 1e-12 * sum(abs, d)   # bound-bound channels conserve; only the continuum changes the total
        end
    end
    @testset "state order and ground-state electron fraction" begin
        @test NATIVE_HYDROGEN_ATOM.ly_index == [3, 5]
        d = zeros(6); continuum!(d, 0.5, 2.0, ones(6), ones(5), 2 .* ones(5)); @test d[1] == 0 && d[2:6] == fill(2 * (1.0 * 0.5 * 2.0 - 1.0), 5)
        @test proton_fraction(0.3) == 0.7
        @test electron_fraction(0.3, 0.08, 0.05) == (1 - 0.3) + (0.08 - 0.05)
        # NOT 1 - sum of all populations: excited populations are not subtracted
        X = [0.3, 0.1, 0.1, 0.1, 0.1, 0.1]
        @test proton_fraction(X[1]) != 1 - sum(X)
    end
    @testset "thermal equation sign and units" begin
        Tg = 3000.0; Hz = 1.0e-13
        @test matter_temperature_rate(1.0, Tg, 0.5, 0.08, Hz) == -Hz
        @test matter_temperature_rate(0.0, Tg, 0.5, 0.08, Hz) > 0      # Compton heats a cold electron gas
        @test matter_temperature_rate(2.0, Tg, 0.5, 0.08, 0.0) < 0     # Compton cools a hot one
        @test matter_temperature_rate(1.0, Tg, 0.0, 0.08, Hz) == -Hz
        # independent form: Compton rate coefficient = (8/3) sigT c a_rad Tg^4 /(me c^2) * Xe/(1+Xe+fHe); a_rad = pi^2 kB^4/(15 hbar^3 c^3)
        arad = K.pi^2 * K.kB^4 / (15 * K.hbar^3 * K.c^3)
        coef = 8 / 3 * K.sigT * arad * Tg^4 / (K.me_gr * K.c) * 0.5 / (1 + 0.5 + 0.08)
        @test isapprox(matter_temperature_rate(0.0, Tg, 0.5, 0.08, 0.0), coef; rtol = 1e-12)
        @test isapprox(rho_g_fac(K), arad / (K.c^2 * K.me_gr); rtol = 1e-14)
    end
    @testset "detailed balance equilibrium (independent)" begin
        # 2s-1s: zero net flow when X2s = X1s exp(-h_kb nu / Tg); Lyman: zero when Xnp = w X1s exp(-h_kb nu/Tg)
        Tg = 4000.0; at = NATIVE_HYDROGEN_ATOM
        X = [0.6, 0.0, 0.0, 0.0, 0.0, 0.0]; X[2] = X[1] * exp(-K.h_kb * at.nu21_2s / Tg)
        d = zeros(6); two_photon!(d, Tg, X); @test abs(d[1]) <= 1e-12 * K.A2s1s * X[1] * exp(-K.h_kb * at.nu21_2s / Tg)
        X[3] = 3 * X[1] * exp(-K.h_kb * at.ly_nu21[1] / Tg)
        d = zeros(6); lyman_channel!(d, 3, Tg, X[1], X[3], 1.0e3, 1.0e-13, at.ly_A21[1], at.ly_lambda21[1], at.ly_nu21[1]); @test abs(d[1]) <= 1e-12 * abs(X[3]) * at.ly_A21[1]
    end
end
