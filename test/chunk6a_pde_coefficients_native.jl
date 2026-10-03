# Chunk 6a: HI diffusion-PDE coefficients (population splines, get_rates_all, pd and Rp/Rm, pd/Dnem coefficient splines) vs the ORIGINAL routines after one native pass.
# REQUIRES COSMOREC_NATIVE_DATA_DIR (fails loudly if unset). Scope: docs/ACCEPTANCE_CHUNK6A.md.
using Test
using CosmoRec
using SHA
include("chunk6_helpers.jl")

const OBS6A = Dict{String,Float64}()
rec6a(k, v) = (OBS6A[k] = max(get(OBS6A, k, 0.0), v))
rel6(a, b) = a == b ? 0.0 : abs(a - b) / max(abs(b), floatmin())

@testset "Chunk 6a: HI PDE coefficients vs native" begin
    rowsN = NODES5[:, 1:9]                           # native populations (pass_on_the_Solution_CosmoRec)
    pops = HIPopulationSplines(rowsN; zs = ZS6, ze = ZE6)
    @testset "fixture provenance and configuration" begin
        @test bytes2hex(open(sha256, FIX6A)) == FIX6A_SHA256
        c = FX6A["CFG6"]
        @test c["zs"] == ZS6 && c["ze"] == ZE6 && c["nres"] == 5 && c["nrows"] == 3000 && c["A2s1s"] == NATIVE_HI_PDE_LEVELS.A2s1s
        for (k, r) in enumerate(FX6A["LVL6"])
            lv = NATIVE_HI_PDE_LEVELS
            @test r[2] == lv.index[k] && r[3] == lv.n[k] && r[4] == lv.l[k] && r[5] == lv.A21[k] && r[6] == lv.lambda21[k] && r[7] == lv.Dnu_1s[k]
        end
        @test size(LNBITOT6) == (500, 5)
        @test_throws ArgumentError HIPopulationSplines(rowsN; zs = 3001.0, ze = ZE6)
        @test_throws ArgumentError HIPopulationSplines(rowsN; zs = ZS6, ze = 50.0)      # the last stored row (z = 50) is not used natively
    end

    @testset "population splines" begin
        for r in FX6A["POP"]
            z = r[1]
            e = max(rel6(hi_Xe(pops, z), r[2]), rel6(hi_rho(pops, z), r[3]), maximum(rel6(hi_Xi(pops, z, i), r[4 + i]) for i in 0:5)); rec6a("pop", e)
            @test e < 1e-13
        end
    end

    @testset "get_rates_all" begin
        for r in FX6A["GRA"]
            A, B, Bt, R = get_rates_all(HTAB6, LNBITOT6, r[1], r[1] * r[2])
            nat = r[3:end]
            vals = vcat(A, B, Bt, vec(permutedims(R)))
            e = maximum(rel6(vals[k], nat[k]) for k in eachindex(nat) if nat[k] != 0 || vals[k] != 0); rec6a("get_rates_all", e)
            @test e < 1e-13
            @test all(R[m, m] == 0 for m in 1:5)
        end
    end

    @testset "Rp/Rm and pd at PDE redshifts" begin
        for r in FX6A["RPD"]
            z, Tg, Te, NH = r[1], r[2], r[3], r[4]
            rec6a("bg:TCMB", rel6(cosmos_TCMB(ACC5, z), Tg)); rec6a("bg:Te", rel6(Tg * hi_rho(pops, z), Te)); rec6a("bg:NH", rel6(cosmos_NH(ACC5, z), NH)); rec6a("bg:H", rel6(cosmos_H(ACC5, z), r[5]))
            @test isapprox(cosmos_TCMB(ACC5, z), Tg; rtol = 1e-14)
            @test isapprox(Tg * hi_rho(pops, z), Te; rtol = 1e-14)
            @test isapprox(cosmos_NH(ACC5, z), NH; rtol = 1e-14)
            @test isapprox(cosmos_H(ACC5, z), r[5]; rtol = 1e-13)
            RpRm, pd = hi_rp_rm_pd(z, pops, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS, Tg, Te, NH * hi_Xe(pops, z), 1.0 - hi_Xi(pops, z, 0), NH)
            e = max(maximum(rel6.(RpRm, r[6:10])), maximum(rel6.(pd, r[11:15]))); rec6a("RpRm_pd", e)
            @test e < 1e-13
            @test pd[3] == 1.0 && pd[5] == 1.0          # 3s, 3d: native A21(1,0) = 0
        end
    end

    @testset "pd / Dnem coefficient splines" begin
        cs = hi_pde_coefficients(rowsN, pops, ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6)
        @test cs.z[1] == ZE6 / 1.0001
        @test isapprox(cs.z[end], ZS6 * 1.0001; rtol = 1e-12)
        @test length(cs.z) == count(z -> z >= ZE6 / 1.0001, rowsN[:, 1])
        for r in FX6A["PDN"]
            z = r[1]
            epd = maximum(rel6(hi_pd(cs, z, m), r[1 + m]) for m in 1:5)
            eDs = maximum(abs(hi_Dnem(cs, z, m) - r[6 + m]) / dnem_scale6(pops, ACC5, z, m, r[1 + m]) for m in 1:5)      # in units of the residual scale
            eDr = maximum(rel6(hi_Dnem(cs, z, m), r[6 + m]) for m in 1:5)                                                 # raw relative (reported)
            eDa = maximum(abs(hi_Dnem(cs, z, m) - r[6 + m]) for m in 1:5)                                                 # absolute (reported)
            rec6a("pd_spline", epd); rec6a("Dnem_spline_scaled", eDs); rec6a("Dnem_spline_raw", eDr); rec6a("Dnem_spline_abs", eDa)
            @test epd < 1e-12
            @test eDs < 1e-12
        end
    end

    @testset "end to end: coefficients from the Julia pass" begin
        pass = recombination_pass(RM5, (a...) -> solve_phase5(a...; reltol = 1.0e-12, a1 = 1.0e-18, aex = 1.0e-14))   # a1 1e-18: production forward config (was 1e-16)
        rowsJ = pass_rows6(pass)
        popsJ = HIPopulationSplines(rowsJ; zs = ZS6, ze = ZE6)
        csJ = hi_pde_coefficients(rowsJ, popsJ, ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6)
        for r in FX6A["PDN"]
            z = r[1]
            rec6a("e2e:pd", maximum(rel6(hi_pd(csJ, z, m), r[1 + m]) for m in 1:5))
            rec6a("e2e:Dnem_raw", maximum(rel6(hi_Dnem(csJ, z, m), r[6 + m]) for m in 1:5))
            rec6a("e2e:Dnem_scaled", maximum(abs(hi_Dnem(csJ, z, m) - r[6 + m]) / dnem_scale6(pops, ACC5, z, m, r[1 + m]) for m in 1:5))
            rec6a("e2e:pop_excited", maximum(rel6(hi_Xi(popsJ, z, i), hi_Xi(pops, z, i)) for i in 1:5))
        end
        # the Julia-pass coefficients differ from the native ones only through the population differences of the two passes (native solver tolerance level, Chunk 5b)
        @test OBS6A["e2e:pd"] < 1e-5
        @test OBS6A["e2e:Dnem_scaled"] < 10 * OBS6A["e2e:pop_excited"] + 1e-12
    end

    @testset "info" begin
        foreach(kv -> println("Chunk6a-native ", kv[1], " = ", kv[2]), sort!(collect(OBS6A)))
    end
end
