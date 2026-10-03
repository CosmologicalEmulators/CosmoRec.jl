# Chunk 8a: correction integrals over the HI PDE spectrum (DI1_2s, DF_2gamma 3s/3d, DF_Raman 2s) vs the ORIGINAL run, single outputs from native spectra
# and the full Julia PDE stage. REQUIRES COSMOREC_NATIVE_DATA_DIR (fails loudly if unset). Scope: docs/ACCEPTANCE_CHUNK8A.md.
using Test
using CosmoRec
using SHA
using Random
@isdefined(SETUP7) || include("chunk7_helpers.jl")
@isdefined(MODEL7N) || include("chunk7b_pde_define_native.jl")
@isdefined(FX7C) || include("chunk7c_pde_march_native.jl")

const OBS8A = Dict{String,Float64}()
rec8a(k, v) = (OBS8A[k] = max(get(OBS8A, k, 0.0), v))
const FX8A = read_fixture7(FIX8A)
const NAMES8A = ("DI1_2s", "DF_2g_3s", "DF_2g_3d", "DF_R_2s")
nat8a() = (vec7(FX8A["DI1_2s"][1]), vec7(FX8A["DF_2g"][1]), vec7(FX8A["DF_2g"][2]), vec7(FX8A["DF_R"][1]))
# FIXED acceptance bound (independent of any Julia/native comparison): the ORIGINAL integrals are only converged to epsrel_HI = 1e-5 (Patterson
# stopping rule) in each of at most 9 sub-integrals of compute_integral_over_resonances, each relative to the running sum (<= residual scale);
# the native outputs therefore carry up to ~9e-5 of their residual scale. Errors up to 1e-4 of the residual scale are accepted, nothing larger.
const BOUND8A = 1.0e-4
const FIX8A_NODES = joinpath(@__DIR__, "fixtures", "native_hi_pde_coefficient_nodes.txt")
const FIX8A_NODES_SHA256 = "a4d66264f76f194966da7093be26cea723b8d6588db03aa4fd6c64095ad1dd24"

@testset "Chunk 8a: HI PDE correction integrals vs native" begin
    @testset "fixture provenance and Sobolev escape" begin
        @test bytes2hex(open(sha256, FIX8A)) == FIX8A_SHA256
        @test length(FX8A["DF_2g"]) == 2 && length(FX8A["DF_R"]) == 1 && !haskey(FX8A, "DF_Ly") && !haskey(FX8A, "DF_nD")
        zs = vec7(FX8A["DF_z"][1])
        @test length(zs) == 199 && zs[1] == 2480.0 && zs[end] == 500.0
        @test sobolev_p_ij(1e-11) == 1.0 - 0.5 * 1e-11 + 1e-22 / 6.0 && sobolev_p_ij(600.0) == 1.0 / 600.0 && sobolev_p_ij(2.0) == (1.0 - exp(-2.0)) / 2.0
    end
    nat = nat8a()
    @testset "single outputs from the native spectra" begin
        for r in FX7C["STEP"]
            k = Int(r[1]); zout = r[3]
            (k >= 1 && haskey(FX7C, "YOUT_$k")) || continue
            y = vec7(FX7C["YOUT_$k"][1])
            st = hi_pde_rhs_coefficients(MODEL7N, zout)
            res = hi_pde_integrals(MODEL7N, zout, y, st)
            @test vec7(FX8A["DF_z"][1])[k] == zout
            for q in 1:4
                ref = nat[q][k]
                rec8a("single:" * NAMES8A[q] * ":scaled", abs(res.values[q] - ref) / res.scale[q])
                rec8a("single:" * NAMES8A[q] * ":raw", abs(res.values[q] - ref) / abs(ref))
                rec8a("single:" * NAMES8A[q] * ":abs", abs(res.values[q] - ref))
                rec8a("single:" * NAMES8A[q] * ":cancellation", res.scale[q] / abs(ref))
            end
        end
        for q in NAMES8A
            @test OBS8A["single:$q:scaled"] <= BOUND8A
        end
    end
    @testset "full Julia PDE stage (march + integrals) from y = 0" begin
        out = hi_pde_corrections(MODEL7N)
        @test out.z == vec7(FX8A["DF_z"][1])
        jl = (out.DI1_2s, out.DF_2g[1], out.DF_2g[2], out.DF_R[1])
        for q in 1:4
            rec8a("full:" * NAMES8A[q] * ":raw", maximum(abs.(jl[q] .- nat[q]) ./ abs.(nat[q])))
            rec8a("full:" * NAMES8A[q] * ":abs", maximum(abs.(jl[q] .- nat[q])))
            rec8a("full:" * NAMES8A[q] * ":rel_to_max", maximum(abs.(jl[q] .- nat[q])) / maximum(abs.(nat[q])))
            rec8a("full:" * NAMES8A[q] * ":scaled", maximum(abs.(jl[q] .- nat[q]) ./ out.scale[q]))     # Julia-side residual scale per output
        end
        rec8a("full:y_final", maximum(abs.(out.y .- fx7c_vec("NATFINAL"))) / maximum(abs.(fx7c_vec("NATFINAL"))))
        @test OBS8A["full:y_final"] <= 1e-9
        for q in NAMES8A
            @test OBS8A["full:$q:scaled"] <= BOUND8A
        end
    end
    @testset "DIAGNOSTIC (no acceptance role): native pd/Dnem knots injected, PDE + integrals isolated" begin
        @test bytes2hex(open(sha256, FIX8A_NODES)) == FIX8A_NODES_SHA256
        nd = read_fixture7(FIX8A_NODES)
        z = vec7(nd["NODE_z_1"][1])
        @test z == COEF7N.z
        CN = HIPDECoefficientSplines(z, [natural_cubic_spline(z, vec7(nd["NODE_lnpd_$m"][1])) for m in 1:5], [natural_cubic_spline(z, vec7(nd["NODE_Dnem_$m"][1])) for m in 1:5])
        @test all(r -> all(hi_Dnem(CN, r[1], m) == r[6 + m] for m in 1:5), FX6A["PDN"])      # the injected knots reproduce the native Dnem exactly
        out = hi_pde_corrections(HIPDEModel(SETUP7, POPS7N, CN, ACC5))
        jl = (out.DI1_2s, out.DF_2g[1], out.DF_2g[2], out.DF_R[1])
        for q in 1:4
            rec8a("diag:native_knots:" * NAMES8A[q] * ":scaled", maximum(abs.(jl[q] .- nat[q]) ./ out.scale[q]))
        end
        rec8a("diag:native_knots:y_final", maximum(abs.(out.y .- fx7c_vec("NATFINAL"))) / maximum(abs.(fx7c_vec("NATFINAL"))))
    end
    @testset "DIAGNOSTIC (no acceptance role): ulp-level conditioning of the outputs (Julia only, native populations x (1 + eps randn))" begin
        base = hi_pde_corrections(MODEL7N)
        b = (base.DI1_2s, base.DF_2g[1], base.DF_2g[2], base.DF_R[1])
        for seed in 1:2
            rng = Random.Xoshiro(seed); rows = copy(ROWS7N); rows[:, 2:9] .*= (1 .+ eps() .* randn(rng, size(rows, 1), 8))
            pops = HIPopulationSplines(rows; zs = ZS6, ze = ZE6); cs = hi_pde_coefficients(rows, pops, ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6)
            o = hi_pde_corrections(HIPDEModel(SETUP7, pops, cs, ACC5))
            oo = (o.DI1_2s, o.DF_2g[1], o.DF_2g[2], o.DF_R[1])
            for q in 1:4
                rec8a("diag:ulp_response:" * NAMES8A[q] * ":scaled", maximum(abs.(oo[q] .- b[q]) ./ base.scale[q]))
            end
        end
        @test true
    end
    for k in sort(collect(keys(OBS8A)))
        println("Chunk8a ", rpad(k, 34), " ", OBS8A[k])
    end
end
