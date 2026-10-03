# Chunk 4a: explicit-input Saha initialization + copy_LI_to_ysol packing vs the native fixture (12 production-cosmology states). Scope: docs/ACCEPTANCE_CHUNK4A.md.
using Test
using CosmoRec
using SHA
@isdefined(FX4A) || include("chunk4a_helpers.jl")

const OBS4A = Dict{String,Float64}()
rec4a(k, v) = (OBS4A[k] = max(get(OBS4A, k, 0.0), v))
rel4a(a, b) = a == b ? 0.0 : abs(a - b) / abs(b)

@testset "Chunk 4a: native Saha initialization and ysol packing" begin
    zs = sort(collect(keys(FX4A["XLI"])))
    @testset "fixture provenance and layout" begin
        @test bytes2hex(open(sha256, FIX4A)) == FIX4A_SHA256
        @test length(zs) == 12 && Set(zs) == Set(keys(FX4A["SAHAIN"])) && Set(zs) == Set(keys(FX4A["YSOL"])) && Set(zs) == Set(keys(FX4A["LTEH"]))
        @test FX4A["LI"] == [15, 6, 7, 1, 7, 5, 4]
        @test FX4A["RESHI"] == collect(NATIVE_RESHI) && FX4A["RESHE"] == collect(NATIVE_RESHE)
        @test FX4A["fHe"] == 0.082009323117224586
        # constants and level data embedded in the source equal the native-printed ones
        C = FX4A["CONST"]
        @test NATIVE_SAHA_CONSTANTS.lambdac == C["lambdac"] && NATIVE_SAHA_CONSTANTS.kb_mec2 == C["kb_mec2"] && NATIVE_SAHA_CONSTANTS.kB == C["kB"] && NATIVE_SAHA_CONSTANTS.pi == C["PI"]
        for (i, r) in enumerate(FX4A["LVLH"])
            @test NATIVE_SAHA_HYDROGEN.g[i] == 2 * (2 * r[2] + 1) && NATIVE_SAHA_HYDROGEN.Eion[i] == r[3] && NATIVE_SAHA_HYDROGEN.mu_red[i] == r[4] && NATIVE_SAHA_HYDROGEN.ME_scale[i] == r[5]
        end
        for (i, r) in enumerate(FX4A["LVLHE"])
            @test NATIVE_SAHA_HELIUM.g[i] == r[1] && NATIVE_SAHA_HELIUM.Eion[i] == r[2] && NATIVE_SAHA_HELIUM.mu_red[i] == r[3] && NATIVE_SAHA_HELIUM.ME_scale[i] == r[4]
        end
        @test length(FX4A["LVLH"]) == 6 && length(FX4A["LVLHE"]) == 7
        # the clip is inactive at every native state (so the clipped value printed natively is the unclipped Xe_Seager used by helium); raw Xe_Seager agrees where printed
        for z in zs
            @test FX4A["SAHAIN"][z]["Xe_clipped"] < 1 + FX4A["fHe"]
            haskey(FX4A["SP"], z) && @test FX4A["SP"][z]["Xe_Seager"] == FX4A["SAHAIN"][z]["Xe_clipped"]
        end
    end

    @testset "LTE factors per level" begin
        for z in zs
            s = FX4A["SAHAIN"][z]
            for i in 1:6
                e = rel4a(saha_lte_hydrogen(NATIVE_SAHA_HYDROGEN, i, s["Te"]), FX4A["LTEH"][z][i]); rec4a("lteH", e)
                @test e < 1e-14
            end
            for i in 1:7
                e = rel4a(saha_lte_helium(NATIVE_SAHA_HELIUM, i, s["TCMB"]), FX4A["LTEHE"][z][i]); rec4a("lteHe", e)
                @test e < 1e-14
            end
        end
    end

    @testset "complete X vector at every native state" begin
        for z in zs
            X = saha_initial_state(inputs4a(z), NATIVE_SAHA_HYDROGEN, NATIVE_SAHA_HELIUM)
            ref = FX4A["XLI"][z]
            @test length(X) == 15
            @test X[1] == ref[1]                 # Xe slot
            @test X[15] == 1.0 && ref[15] == 1.0  # rho exactly 1
            for i in 1:15
                e = rel4a(X[i], ref[i]); rec4a("X", e)
                @test e < 1e-14
            end
            @test X[2] == ref[2]                  # H ground override 1 - Xp: bitwise
            @test X[8] == ref[8]                  # He ground override: bitwise
        end
    end

    @testset "native-LTE products (bitwise) and packing (exact) vs YSOL" begin
        for z in zs
            s = FX4A["SAHAIN"][z]; ref = FX4A["XLI"][z]
            for i in 1:5      # H levels 2s..3d: X[1 + i + 1]; LTEH index i + 1
                @test s["Xe_clipped"] * s["Xp_clipped"] * s["NH"] * FX4A["LTEH"][z][i + 1] == ref[2 + i]
            end
            for i in 1:6      # He levels 2..7
                @test s["Xe_clipped"] * s["XHeII"] * s["NH"] * FX4A["LTEHE"][z][i + 1] == ref[8 + i]
            end
            @test pack_ysol(ref, 6, 7) == FX4A["YSOL"][z]            # native X -> native y: exact
            X = saha_initial_state(inputs4a(z), NATIVE_SAHA_HYDROGEN, NATIVE_SAHA_HELIUM)
            y = pack_ysol(X, 6, 7)
            @test length(y) == 12 && length(FX4A["YSOL"][z]) == 12
            for i in 1:12
                e = rel4a(y[i], FX4A["YSOL"][z][i]); rec4a("y", e)
                @test e < 1e-14
            end
            @test y[1] == 1.0
        end
    end

    @testset "branch edges and overrides (synthetic explicit inputs)" begin
        H, He = NATIVE_SAHA_HYDROGEN, NATIVE_SAHA_HELIUM
        base = inputs4a(3000.0)
        # Xe clip: H uses min(Xe, 1+fHe); He uses the UNCLIPPED Xe
        big = SahaInputs(base.fHe, 2.0, base.Xp_raw, base.NH, base.Te, base.Tg, base.XHeII, base.XHeI1s)
        Xb = saha_initial_state(big, H, He)
        @test Xb[1] == 1 + base.fHe
        @test Xb[3] == (1 + base.fHe) * min(base.Xp_raw, 1.0) * base.NH * saha_lte_hydrogen(H, 2, base.Te)
        @test Xb[9] == 2.0 * base.XHeII * base.NH * saha_lte_helium(He, 2, base.Tg)
        # Xp clip and the ground override 1 - Xp
        xp = SahaInputs(base.fHe, base.Xe_Seager, 1.5, base.NH, base.Te, base.Tg, base.XHeII, base.XHeI1s)
        Xx = saha_initial_state(xp, H, He)
        @test Xx[2] == 0.0 && Xx[3] == Xx[1] * 1.0 * base.NH * saha_lte_hydrogen(H, 2, base.Te)
        # ground states are overrides: He ground = XHeI1s regardless of the Saha formula, H ground ignores the Saha product
        Xg = saha_initial_state(SahaInputs(base.fHe, base.Xe_Seager, base.Xp_raw, base.NH, base.Te, base.Tg, base.XHeII, 0.123), H, He)
        @test Xg[8] == 0.123 && Xg[2] == 1.0 - base.Xp_raw
        # H uses Te, He uses Tg
        Xt = saha_initial_state(SahaInputs(base.fHe, base.Xe_Seager, base.Xp_raw, base.NH, base.Te * 1.01, base.Tg, base.XHeII, base.XHeI1s), H, He)
        X0 = saha_initial_state(base, H, He)
        @test Xt[9:14] == X0[9:14] && Xt[3] != X0[3]
        @test X0[15] == 1.0
    end

    @testset "shape and invalid dimensions" begin
        @test_throws DimensionMismatch SahaLevels([1.0, 2.0], [1.0], [1.0, 1.0], [1.0, 1.0])
        @test_throws DimensionMismatch SahaLevels(Float64[], Float64[], Float64[], Float64[])
        @test_throws DimensionMismatch pack_ysol(zeros(14), 6, 7)
        @test_throws BoundsError pack_ysol(zeros(15), 6, 7, (1, 2, 3, 4, 6))
        @test_throws BoundsError pack_ysol(zeros(15), 6, 7, NATIVE_RESHI, (1, 2, 3, 7))
        @test_throws BoundsError pack_ysol(zeros(15), 6, 7, (-1,), (1,))
        small = SahaLevels([2.0, 2.0], [2.0e-11, 5.0e-12], [0.999, 0.999], [1.0, 1.0])
        Xs = saha_initial_state(inputs4a(3000.0), small, small)
        @test length(Xs) == 1 + 2 + 2 + 1
        @test pack_ysol(Xs, 2, 2, (1,), (1,)) == [Xs[end], Xs[2], Xs[3], Xs[4], Xs[5]]
        @test eltype(saha_initial_state(inputs4a(3000.0), NATIVE_SAHA_HYDROGEN, NATIVE_SAHA_HELIUM)) == Float64
        @inferred saha_initial_state(inputs4a(3000.0), NATIVE_SAHA_HYDROGEN, NATIVE_SAHA_HELIUM)
        @inferred pack_ysol(zeros(15), 6, 7)
    end

    @testset "info" begin
        foreach(kv -> println("Chunk4a-native ", kv[1], " = ", kv[2]), sort!(collect(OBS4A)))
    end
end
