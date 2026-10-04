# Chunk 5a: the packed-state ODE right-hand side dy/dz (12-state flag_He = 1 and 7-state flag_He = 0) vs the ORIGINAL exported fcn_effective(int*, ...), with the full production tables.
# REQUIRES COSMOREC_NATIVE_DATA_DIR (native Rec_database, read-only, never bundled): fails loudly if unset (no skip). Scope: docs/ACCEPTANCE_CHUNK5A.md.
using Test
using CosmoRec
using SHA
include("chunk5_helpers.jl")

const FIX5A = joinpath(@__DIR__, "fixtures", "native_ode_rhs.txt")
const FIX5A_SHA256 = "f9b9e03dd14acc91267d97ece7599fc8c49b223508e5e221557535737ae09906"
const ROWS5A = read_fixture5a(FIX5A)
const OBS5A = Dict{String,Float64}()
rec5a(k, v) = (OBS5A[k] = max(get(OBS5A, k, 0.0), v))

@testset "Chunk 5a: ODE right-hand side vs native fcn_effective(int*, ...)" begin
    @testset "data location and fixture provenance" begin
        @test isdir(DATADIR5) && isfile(fcorr_path5(DATADIR5))
        @test bytes2hex(open(sha256, FIX5A)) == FIX5A_SHA256
        @test length(ROWS5A) == 132 && count(r -> r.flag, ROWS5A) == 60 && count(r -> !r.flag, ROWS5A) == 72
        @test all(r -> length(r.y) == ode_nstate(r.flag) && length(r.f) == ode_nstate(r.flag), ROWS5A)
        @test ode_nstate(true) == 12 && ode_nstate(false) == 7
    end

    @testset "background from the 4b accessors reproduces the native background" begin
        for r in ROWS5A
            @test isapprox(cosmos_TCMB(ACC5, r.z), r.Tg; rtol = 1e-14)
            @test isapprox(cosmos_NH(ACC5, r.z), r.NH; rtol = 1e-14)
            @test isapprox(cosmos_H(ACC5, r.z), r.Hz; rtol = 1e-13)
            @test isapprox(cosmos_Xp(ACC5, r.z), r.Xp; rtol = 1e-12, atol = 1e-15)
        end
    end

    @testset "every native state: dy/dz component-wise" begin
        for r in ROWS5A
            f = recombination_rhs(r.z, r.y, RM5; flag_He = r.flag)
            @test length(f) == length(r.f)
            sc = maximum(abs, r.f)
            for i in eachindex(f)
                e = abs(f[i] - r.f[i]) / max(abs(r.f[i]), 1e-12 * sc)
                rec5a(r.flag ? "flag1" : "flag0", e)
                @test e < 1e-12
            end
            if !r.flag      # helium off: hydrogen only, xe = xp
                @test length(f) == 7
            end
        end
    end

    @testset "flag_He = 0 branch" begin
        r = first(r for r in ROWS5A if !r.flag && r.z == 400.0 && r.pat == 0)
        f0 = recombination_rhs(r.z, r.y, RM5; flag_He = false)
        # the He state does not enter: the same hydrogen state through the 12-state RHS with flag_He = 0 semantics gives the same hydrogen derivatives
        X = ode_unpack(r.y, ode_background(RM5, r.z).fHe, false)
        @test X[1] == 1 - r.y[2] && X[15] == r.y[1] && X[2:7] == r.y[2:7]
        @test all(==(1.0e-300), X[8:14])
        @test_throws DimensionMismatch ode_unpack(r.y, 0.08, true)
        @test_throws DimensionMismatch recombination_rhs(r.z, zeros(12), RM5; flag_He = false)
        # unresolved helium triplet slots are placeholders: changing them cannot change the 12-state derivative
        r1 = first(r for r in ROWS5A if r.flag && r.z == 1500.0 && r.pat == 0)
        fa = recombination_rhs(r1.z, r1.y, RM5)
        Xa = ode_unpack(r1.y, ode_background(RM5, r1.z).fHe, true); Xb = copy(Xa); Xb[12] = 1.0e-10; Xb[14] = 1.0e-9
        ga = zeros(15); gb = zeros(15); bg = ode_background(RM5, r1.z)
        fcn_effective!(ga, r1.z, Xa, bg, EFF5; dp_fallback = RM5.dp_fallback); fcn_effective!(gb, r1.z, Xb, bg, EFF5; dp_fallback = RM5.dp_fallback)
        @test ga == gb
    end

    @testset "dp_fallback callable is specialized (no extra allocation from the pass-through)" begin
        # in-table 12-state point: the factory fallback is never called, so the RHS with it must allocate exactly as the no-fallback model does;
        # a despecialized (::Function) pass-through boxes the callable and the keyword tuples on every call and fails this
        r = first(r for r in ROWS5A if r.flag && r.z == 2500.0 && r.pat == 0)
        rmnf = RecombinationModel(RM5.eff, RM5.cosmos, nothing)
        du = similar(r.y); dunf = similar(r.y)
        alloc5a!(du, r, m) = @allocated recombination_rhs!(du, r.z, r.y, m; flag_He = true)
        alloc5a!(du, r, RM5); alloc5a!(dunf, r, rmnf)           # warm (compile) both
        alloc5a!(du, r, RM5); alloc5a!(dunf, r, rmnf)
        @test du == dunf
        @test alloc5a!(du, r, RM5) == alloc5a!(dunf, r, rmnf)
    end

    @testset "info" begin
        foreach(kv -> println("Chunk5a-native ", kv[1], " = ", kv[2]), sort!(collect(OBS5A)))
    end
end
