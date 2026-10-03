# Chunk 4b: explicit Cosmos accessors / GSL splines vs the native fixtures (history nodes, Hubble input table, accessor oracle). Scope: docs/ACCEPTANCE_CHUNK4B.md.
using Test
using CosmoRec
using SHA
@isdefined(ACC4B) || include("chunk4b_helpers.jl")

const TOL4B = Dict("H" => 1e-14, "TCMB" => 1e-14, "Nb" => 1e-14, "NH" => 1e-14, "Xe_Seager" => 1e-14, "dXe_dz" => 1e-12, "Xe_b" => 1e-14, "X1s" => 1e-12, "Xp" => 1e-14,
                   "XHeI1s" => 1e-13, "XHeII1s" => 1e-13, "NHeI" => 1e-13, "NHeII" => 1e-13, "NHeIII" => 1e-5, "Te_Tg" => 1e-14, "Te" => 1e-14, "kappa_cool" => 1e-13,
                   "Ne" => 1e-14, "Ntot" => 1e-14, "sigT" => 0.0, "rho_g_1_cm3" => 1e-14)   # observed maxima in docs/CHUNK4B_RESULTS.md; NHeIII is a cancellation (fHe - X - X) at z >= 3500
const OBS4B = Dict{String,Float64}()
rec4b(k, v) = (OBS4B[k] = max(get(OBS4B, k, 0.0), v))
relerr4b(a, b) = a == b ? 0.0 : abs(a - b) / max(abs(b), floatmin())

@testset "Chunk 4b: Cosmos accessors vs native fixtures" begin
    A = accessors4b()
    @testset "fixture provenance and layout" begin
        for (p, h) in FIX4B_SHA
            @test bytes2hex(open(sha256, p)) == h
        end
        @test size(HIST4B) == (6000, 7) && size(HUB4B) == (10000, 2)
        @test issorted(HIST4B[:, 1]; rev = true) && HIST4B[1, 1] == 25000.0 && HIST4B[end, 1] == 0.0
        @test HUB4B[1, 1] == 1.0e4 && HUB4B[end, 1] == 0.0 && issorted(HUB4B[:, 1]; rev = true)
        @test length(ACC4B.acc) >= 600 && length(ACC4B.sbt) == 6
        @test ACC4B.cos2["n_Xe"] == 6000 && ACC4B.cos2["zsRe"] == 25000.0
        # native column meaning: Xe_H, Xe_He are IONIZED populations (start at 1 and fHe), not neutral fractions
        @test HIST4B[1, 2] == 1.0 && HIST4B[1, 3] == cosmos_fHe(A)
    end

    @testset "history nodes: splines reproduce the inputs, grid independently verified" begin
        zg = reverse(HIST4B[:, 1])
        @test issorted(zg) && all(diff(zg) .> 0)
        sp = A.sp
        for i in 1:97:6000
            z = HIST4B[i, 1]
            i == 1 && continue            # z = 25000 is zsRe (closed form branch); nodes are tested below zsRe
            @test isapprox(exp(spline_eval(sp.lnXe, z)), HIST4B[i, 4]; rtol = 1e-12)
            @test isapprox(spline_eval(sp.rho, z), HIST4B[i, 7] / (2.7254999999999998 * (1 + z)); rtol = 1e-12)
        end
        # natural boundary: zero curvature at both ends (c[1] = c[end] = 0)
        @test sp.lnXe.c[1] == 0.0 && sp.lnXe.c[end] == 0.0
    end

    @testset "every native accessor value (658 redshifts)" begin
        fns = (("H", cosmos_H), ("TCMB", cosmos_TCMB), ("Nb", cosmos_Nb), ("NH", cosmos_NH), ("Xe_Seager", cosmos_Xe_Seager), ("dXe_dz", cosmos_dXe_dz),
               ("Xe_b", cosmos_Xe_b), ("X1s", cosmos_X1s), ("Xp", cosmos_Xp), ("XHeI1s", cosmos_XHeI1s), ("XHeII1s", cosmos_XHeII1s), ("NHeI", cosmos_NHeI),
               ("NHeII", cosmos_NHeII), ("NHeIII", cosmos_NHeIII), ("Te_Tg", cosmos_Te_Tg), ("Te", cosmos_Te), ("kappa_cool", cosmos_kappa_cool), ("Ne", cosmos_Ne),
               ("Ntot", cosmos_Ntot), ("sigT", cosmos_sigT), ("rho_g_1_cm3", cosmos_rho_g_1_cm3))
        for r in ACC4B.acc
            z = r["z"]
            for (name, f) in fns
                v = f(A, z)
                ref = r[name]
                e = ref == 0 ? abs(v) : abs(v - ref) / abs(ref)
                rec4b(name, e)
                @test e <= TOL4B[name]
            end
        end
    end

    @testset "Saha-Boltzmann direct functions" begin
        for r in ACC4B.sbt
            T = r["T"]
            for (name, f) in (("HI1s", cosmos_SahaBoltz_HI1s), ("HeI1s", cosmos_SahaBoltz_HeI1s), ("HeII1s", cosmos_SahaBoltz_HeII1s), ("HeII", cosmos_SahaBoltz_HeII), ("HeIII", cosmos_SahaBoltz_HeIII))
                e = relerr4b(f(A, T), r[name]); rec4b("SBT:" * name, e)
                @test e < 1e-12
            end
        end
    end

    @testset "native low-z Hubble bug reproduced (lHz[0] = 0)" begin
        # H(1e-3) is 10.7x the table value in the native object; the port reproduces it from the same spline input
        r3 = first(r for r in ACC4B.acc if r["z"] == 1.0e-3)
        @test r3["H"] > 1.0e-17
        @test isapprox(cosmos_H(A, 1.0e-3), r3["H"]; rtol = 1e-9)
        @test cosmos_H(A, 1.0e-3) / HUB4B[end - 1, 2] > 5            # far above the true table value (~H0)
        z = copy(HUB4B[:, 1]); Hz = copy(HUB4B[:, 2])
        t = hubble_table(z, Hz)
        @test t.spline.y[1] == 0.0
        @test isapprox(exp(t.spline.x[1]), 1.0e-10; rtol = 1e-14)
        # above the first ~ few table nodes the effect vanishes (z >= 100: native/table agreement 1.4e-8 or better)
        lz = log.(reverse(HUB4B[:, 1])[2:end]); lH = log.(reverse(HUB4B[:, 2])[2:end])
        for zq in (100.0, 1000.0, 3000.0, 9000.0)
            r = ACC4B.acc[argmin(abs.([x["z"] for x in ACC4B.acc] .- zq))]
            j = searchsortedlast(lz, log(r["z"]))
            Hlin = exp(lH[j] + (lH[j + 1] - lH[j]) * (log(r["z"]) - lz[j]) / (lz[j + 1] - lz[j]))
            @test abs(r["H"] / Hlin - 1) < 1e-5            # away from the bug region the native spline follows the table
        end
    end

    @testset "domains, branches, endpoint nudge" begin
        c = A.c
        # z >= zsRe closed forms
        @test cosmos_Xe_Seager(A, 25000.0) == (1.0 - c.Y_p / 2.0 * (2.0 - 1.0 / c.fac_mHemH)) / (1.0 - c.Y_p)
        @test cosmos_Xe_b(A, 30000.0) == 1.0 && cosmos_dXe_dz(A, 25000.0) == 0.0
        @test cosmos_Te_Tg(A, 25000.0) == 1.0 / (1.0 + 1.0 / cosmos_kappa_cool(A, 25000.0))
        # z_saha = 3500 switch: both sides use different formulas and agree to ~1e-3 relative
        lo, hi = cosmos_X1s(A, 3499.999999), cosmos_X1s(A, 3500.0)
        @test lo != hi && abs(lo / hi - 1) < 1e-3
        @test cosmos_NHeIII(A, 3499.0) == 1.0e-20 * cosmos_NH(A, 3499.0)
        @test cosmos_NHeIII(A, 3500.0) > 1.0e-20 * cosmos_NH(A, 3500.0)
        # end nudge: z within 1e-14 of the spline ends evaluates inside; z = 0 is its own end (xmin = 0)
        @test isfinite(cosmos_Xe_Seager(A, 0.0)) && isfinite(cosmos_Xe_Seager(A, 25000.0 * (1 - 5.0e-15) - 1e-6))
        s = A.sp.lnXe
        @test spline_eval_native(s, 25000.0 * (1 - 5.0e-15)) == spline_eval(s, 25000.0 * (1 - 1.0e-14))
        # no extrapolation: out of the node range throws (native: GSL abort); H outside the loaded table is the analytic form
        @test_throws SplineDomainError spline_eval_native(s, -1.0e-3)
        @test_throws SplineDomainError spline_eval_native(s, 25000.001)
        @test_throws SplineDomainError cosmos_Te_Tg(A, -0.5)
        An = accessors4b(hubble = false)
        @test cosmos_H(An, 10001.0) == cosmos_H(A, 10001.0) && cosmos_H(A, 10001.0) == c.H0 * sqrt(c.O_L + 10002.0^2 * (c.O_k + 10002.0 * (c.O_m + c.O_rel * 10002.0)))
        @test cosmos_H(An, 100.0) != cosmos_H(A, 100.0) && abs(cosmos_H(An, 100.0) / cosmos_H(A, 100.0) - 1) < 1e-2
        # invalid constructors
        @test_throws DimensionMismatch natural_cubic_spline([0.0, 1.0, 2.0], [1.0, 2.0])
        @test_throws ArgumentError natural_cubic_spline([0.0, 1.0, 1.0, 2.0], [1.0, 2.0, 3.0, 4.0])
        @test_throws DimensionMismatch recfast_splines(c, HIST4B[:, 1], HIST4B[1:5999, 2], HIST4B[:, 3], HIST4B[:, 4], HIST4B[:, 5], HIST4B[:, 6], HIST4B[:, 7])
        # Xe_He floor (ln max(Xe_He, 1e-20)) acts on the node values
        Hz = copy(HIST4B); Hz[2500, 3] = 0.0
        @test isfinite(spline_eval(splines4b(c, Hz).lnXeHe, Hz[2500, 1]))
    end

    @testset "4a fed by 4b at native states" begin
        # accessor-derived SahaInputs at z = 3000, 1500, 800 reproduce the native Saha vector of the Chunk 4a fixture to the accessor tolerance
        @isdefined(FX4A) || include("chunk4a_helpers.jl")
        for z in (3000.0, 2500.0, 2000.0, 1500.0, 1200.0, 800.0)
            X = saha_initial_state(saha_inputs_at(A, z), NATIVE_SAHA_HYDROGEN, NATIVE_SAHA_HELIUM)
            ref = FX4A["XLI"][z]
            for i in 1:15
                e = relerr4b(X[i], ref[i]); rec4b("4a-fed:X", e)
                @test e < 1e-13
            end
        end
    end

    @testset "info" begin
        foreach(kv -> println("Chunk4b-native ", kv[1], " = ", kv[2]), sort!(collect(OBS4B)))
    end
end
