# Chunk 7a: HI radiation-PDE setup (hydrogenic A_SH, HI_Transition_Data, matrix elements, two-photon/Raman profiles from the ORIGINAL tables, the
# frequency grid and the arm_PDE_solver ratios) vs the ORIGINAL routines. REQUIRES COSMOREC_NATIVE_DATA_DIR (fails loudly if unset). Scope: docs/ACCEPTANCE_CHUNK7A.md.
using Test
using CosmoRec
using SHA
include("chunk7_helpers.jl")

const OBS7A = Dict{String,Float64}()
rec7a(k, v) = (OBS7A[k] = max(get(OBS7A, k, 0.0), v))
# NaN on both sides counts as agreement (the native y = 0 and y = 1 probes are 0/0 forms in the ORIGINAL code too; NaN patterns are tested separately)
rel7(a, b) = (a == b || (isnan(a) && isnan(b))) ? 0.0 : abs(a - b) / max(abs(b), floatmin())
# vector error relative to the largest native magnitude of the array (ratios cross zero; pointwise relative errors are reported separately where defined)
vrel7(a, b) = maximum(abs.(a .- b)) / max(maximum(abs.(b)), floatmin())

@testset "Chunk 7a: HI PDE setup vs native" begin
    @testset "fixture provenance and configuration" begin
        @test bytes2hex(open(sha256, FIX7A)) == FIX7A_SHA256
        @test isfile(joinpath(TPD7, "DGamma.Mnr_ns.5000.dat")) && isfile(joinpath(TPD7, "Raman.Mnr_ns.5000.dat"))
        cfg = FX7A["CFG7"][1]
        @test cfg[1:9] == [3, 3, 2353, 3, 2, 248, 25, 50, 400]
        @test cfg[10] == 0.55 && cfg[11] == 1.0e-4
        @test SETUP7.nresmax == 3 && SETUP7.n2g == 3 && SETUP7.nR == 2
        @test SETUP7.index_2 == cfg[6]
        @test SETUP7.index_emission == Int.(cfg[12:end]) == [782, 1728]
        con = read_kv7(FIX7A, "CON7")
        @test parse(Float64, con["A2s1s"]) == PROF7.A2s1s == 8.2206
        for (k, v) in (("alpha", CosmoRec.PDE_ALPHA), ("cl", CosmoRec.PDE_CL), ("a0", CosmoRec.PDE_A0), ("me_mp", CosmoRec.PDE_ME_MP), ("EH_inf_Hz", CosmoRec.PDE_EH_INF_HZ),
                       ("Ry_inf_icm", CosmoRec.PDE_RY_INF_ICM), ("PI", CosmoRec.PDE_PI), ("FOURPI", CosmoRec.PDE_FOURPI))
            @test parse(Float64, con[k]) == v
        end
        @test hi_pde_setup(PROF7; nS_2gamma = 4).n2g == 3       # capped by nShells = 3 exactly as native switch_on_two_g_corrections
        @test_throws ArgumentError hi_pde_setup(PROF7; nShells = 4, nS_2gamma = 4, nS_Raman = 3)
        @test_throws ArgumentError hi_Rksnp(4, 6)
        @test_throws ArgumentError load_hi_profile_data(joinpath(TPD7, "no_such_dir"))
    end

    @testset "hydrogenic rates and transition data" begin
        @test pde_log10factorial(0) == 0.0 && pde_log10factorial(1) == 0.0
        @test isapprox(pde_log10factorial(10), log10(3628800.0); rtol = 1e-15)
        for r in FX7A["ASH"]
            n, l, np, lp = Int.(r[1:4])
            rec7a("A_SH", rel7(hydrogen_A_SH(n, l, np, lp), r[5]))
        end
        @test hydrogen_A_SH(2, 1, 1, 1) == 0.0 && hydrogen_A_SH(2, 0, 1, 0) == 0.0 && hydrogen_A_SH(1, 0, 2, 1) == 0.0
        td = PROF7.td
        for r in FX7A["GAM"]
            rec7a("Gamma_np", rel7(gamma_np(td, Int(r[1])), r[2]))
        end
        for r in FX7A["ANPKS"], n in 0:10
            rec7a("A_npks", rel7(A_npks(td, Int(r[1]), n), r[n + 2]))
        end
        for r in FX7A["ANPKD"], n in 0:10
            rec7a("A_npkd", rel7(A_npkd(td, Int(r[1]), n), r[n + 2]))
        end
        for r in FX7A["RME"]
            n = Int(r[1])
            rec7a("R1snp", rel7(hi_R1snp(n), r[2])); rec7a("R2snp", rel7(hi_Rksnp(2, n), r[3]))
            rec7a("R3snp", rel7(hi_Rksnp(3, n), r[4])); rec7a("R3dnp", rel7(hi_Rkdnp(3, n), r[5])); rec7a("Rnsnp", rel7(hi_Rksnp(n, n), r[6]))
        end
        for k in ("A_SH", "Gamma_np", "A_npks", "A_npkd", "R1snp", "R2snp", "R3snp", "R3dnp", "Rnsnp")
            @test OBS7A[k] <= 1e-13
        end
        # the HI-atom Ly-n data of the native PDE (Atom.cpp voigt_init) agree with the Chunk 6a level data
        lyn = Dict(Int(r[1]) => r for r in FX7A["LYN"])
        @test lyn[2][4] == NATIVE_HI_PDE_LEVELS.A21[2] && lyn[3][4] == NATIVE_HI_PDE_LEVELS.A21[4]
        @test lyn[2][3] == NATIVE_HI_PDE_LEVELS.lambda21[2] == NATIVE_HI_PDE_ATOM.lyn_lambda21[1] && lyn[3][3] == NATIVE_HI_PDE_LEVELS.lambda21[4] == NATIVE_HI_PDE_ATOM.lyn_lambda21[2]
    end

    @testset "two-photon and Raman profiles (pointwise)" begin
        for r in FX7A["PRF"]
            y = r[1]
            rec7a("sigma_2s1s_2g", rel7(sigma_2s1s_2gamma(PROF7, y), r[2]))
            rec7a("sigma_3s1s_2g", rel7(sigma_ns1s_2gamma(PROF7, 3, y), r[3]))
            rec7a("sigma_3d1s_2g", rel7(sigma_nd1s_2gamma(PROF7, 3, y), r[4]))
            y > 0 && rec7a("ratio_3s_2g", rel7(sigma_ns1s_2gamma_ratio(PROF7, 3, y), r[5]))
            y > 0 && rec7a("ratio_3d_2g", rel7(sigma_nd1s_2gamma_ratio(PROF7, 3, y), r[6]))
            0 < y < 0.3 && rec7a("ratio_2s_Raman", rel7(sigma_2s1s_Raman_ratio(PROF7, y), r[7]))
            @test isnan(sigma_2s1s_2gamma(PROF7, y)) == isnan(r[2]) && isnan(sigma_ns1s_2gamma(PROF7, 3, y)) == isnan(r[3])
        end
        for k in ("sigma_2s1s_2g", "sigma_3s1s_2g", "sigma_3d1s_2g", "ratio_3s_2g", "ratio_3d_2g", "ratio_2s_Raman")
            @test OBS7A[k] <= 1e-12
        end
        @test sigma_2s1s_2gamma(PROF7, -0.1) == 0.0 && sigma_2s1s_2gamma(PROF7, 1.1) == 0.0
        @test sigma_2s1s_2gamma(PROF7, 0.3) == sigma_2s1s_2gamma(PROF7, 0.7) || rel7(sigma_2s1s_2gamma(PROF7, 0.3), sigma_2s1s_2gamma(PROF7, 0.7)) < 1e-12
    end

    @testset "frequency grid and arm_PDE_solver ratios" begin
        xg = vec7(FX7A["XGRID"][1])
        @test length(SETUP7.x) == length(xg) == 2353
        rec7a("grid", maximum(rel7.(SETUP7.x, xg)))
        @test OBS7A["grid"] <= 1e-14
        @test issorted(SETUP7.x) && SETUP7.x[1] == 1.0e-4
        @test SETUP7.resonances == vec7(FX7A["RES7"][1])
        for (tag, v) in (("RAT_Rns", SETUP7.ratio_R_ns), ("RAT_Rnd", SETUP7.ratio_R_nd), ("RAT_2gns", SETUP7.ratio_2g_ns), ("RAT_2gnd", SETUP7.ratio_2g_nd)), m in 0:1
            ref = vec7(FX7A["$(tag)_$m"][1])
            @test length(ref) == 2353
            rec7a(tag, vrel7(v[m + 1], ref))
            nz = findall(!iszero, ref)
            isempty(nz) || rec7a(tag * "_pointwise", maximum(rel7.(v[m + 1][nz], ref[nz])))
            @test all(iszero, v[m + 1][findall(iszero, ref)])
        end
        for k in ("RAT_Rns", "RAT_Rnd", "RAT_2gns", "RAT_2gnd")
            @test OBS7A[k] <= 1e-12
        end
        for m in 0:1
            @test SETUP7.A_npns[m + 1] == vec7(FX7A["ANPNS_$m"][1]) || maximum(rel7.(SETUP7.A_npns[m + 1], vec7(FX7A["ANPNS_$m"][1]))) <= 1e-13
            @test SETUP7.A_npnd[m + 1] == vec7(FX7A["ANPND_$m"][1]) || maximum(rel7.(SETUP7.A_npnd[m + 1], vec7(FX7A["ANPND_$m"][1]))) <= 1e-13
        end
    end
    for k in sort(collect(keys(OBS7A)))
        println("Chunk7a ", rpad(k, 22), " ", OBS7A[k])
    end
end
