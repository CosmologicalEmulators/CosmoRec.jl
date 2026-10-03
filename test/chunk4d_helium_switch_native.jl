# Chunk 4d: the sampled-HeI state switch (discrete decision, reset, state-size change, packing) vs the native fixture (80 explicit switch points, 40 switched / 40 not).
# NOT tested/claimed: any ODE, the solver restart, the event location in z, or a continuous event sensitivity (none exists: the decision is piecewise constant and the state dimension changes).
using Test
using CosmoRec
using SHA
@isdefined(FX4D) || include("chunk4d_helpers.jl")

@testset "Chunk 4d: sampled-HeI state switch vs native" begin
    @testset "fixture provenance and configuration" begin
        @test bytes2hex(open(sha256, FIX4D)) == FIX4D_SHA256
        @test length(FX4D) == 80 && count(r -> r.cond, FX4D) == 40
        @test FX4D_CFG["Xi_HeI_switch"] == NATIVE_XI_HEI_SWITCH == 1.0e-7
        @test FX4D_CFG["neq"] == 15 && FX4D_CFG["index_HeI"] == 7 && FX4D_CFG["nHeIeq"] == 7 && FX4D_CFG["nres_HI"] == 5 && FX4D_CFG["nres_HeI"] == 4
        @test all(r -> r.fHe == FX4D_CFG["fHe"] && r.flag_pre == 1, FX4D)
        # switch points include both triggers (Xi criterion at zs >= 200: 24, zs < 200: 16) and the rounding boundary of the criterion
        @test count(r -> r.cond && r.zs >= 200, FX4D) == 24 && count(r -> r.cond && r.zs < 200, FX4D) == 16
        @test any(r -> r.zs >= 200 && abs(r.diff - 1.0e-7) < 2e-14, FX4D)
    end

    @testset "decision, reset, packing: every native state, bitwise" begin
        nH, nHe = 6, 7
        for r in FX4D
            @test r.Xi0 == r.Xpre[8]
            @test r.diff == r.fHe - r.Xi0
            @test helium_switch_condition(r.zs, r.fHe, r.Xpre[8]) == r.cond
            Xn, sw = helium_switch(r.Xpre, r.zs, nH, nHe, r.fHe)
            @test sw == r.cond
            @test Xn == r.Xpost                                          # bitwise: He 1s = fHe, others 1e-300, rest untouched
            @test Xn[9:14] == r.Xi[2:7] && Xn[8] == r.Xi[1]              # HeI_Atoms.Xi(i) after the reset
            @test pack_ysol(r.Xpre, nH, nHe) == r.ypre
            y = pack_ysol(Xn, nH, nHe; helium = !sw)
            @test y == r.ypost
            @test length(y) == r.neq_y == helium_switch_ysize(5, 4, !sw)
            @test r.flag_post == (r.cond ? 0 : 1)
            if r.cond
                @test all(==(1.0e-300), Xn[9:14]) && Xn[8] == r.fHe && Xn[1] == r.Xpre[1] && Xn[2:7] == r.Xpre[2:7] && Xn[15] == r.Xpre[15]
                @test length(y) == 7
            else
                @test Xn === r.Xpre || Xn == r.Xpre
                @test length(y) == 12
            end
        end
    end

    @testset "decision logic and flag_He" begin
        fHe = FX4D_CFG["fHe"]
        @test helium_switch_condition(150.0, fHe, 0.0)                                 # zs < 200 triggers regardless of the abundance
        @test !helium_switch_condition(200.0, fHe, 0.0)                                # zs = 200 does not (strict <)
        @test helium_switch_condition(1500.0, fHe, fHe - 9.0e-8)                       # criterion fHe - X_He1s <= 1e-7 (at the threshold itself the floating-point difference decides: see the native rows)
        @test !helium_switch_condition(1500.0, fHe, fHe - 1.1e-7)
        @test !helium_switch_condition(150.0, fHe, 0.0; flag_He = false)               # already switched off: no second switch
        X = zeros(15); X[8] = fHe - 1e-6
        Xn, sw = helium_switch(X, 150.0, 6, 7, fHe; flag_He = false)
        @test !sw && Xn === X
    end

    @testset "the reset depends on fHe; the operator is discrete (no event sensitivity)" begin
        r = first(r for r in FX4D if r.cond && r.zs == 2000.0)
        nH, nHe = 6, 7
        X1 = helium_switch_reset(r.Xpre, nH, nHe, r.fHe * 1.01)
        @test X1[8] == r.fHe * 1.01 && X1[9:14] == r.Xpost[9:14]           # He ground follows fHe, floors do not
        # jump across the criterion: just below / above the threshold the post-state differs by the full He 1s reset and the dimension changes
        below = copy(r.Xpre); below[8] = r.fHe - 2.0e-7                     # not switched
        above = copy(r.Xpre); above[8] = r.fHe - 5.0e-8                     # switched
        Xb, sb = helium_switch(below, 2000.0, nH, nHe, r.fHe); Xa, sa = helium_switch(above, 2000.0, nH, nHe, r.fHe)
        @test !sb && sa
        @test abs(Xa[8] - Xb[8]) ≈ 2.0e-7 rtol = 1e-6                     # He 1s: fHe vs fHe - 2e-7 (discontinuous in the input)
        @test Xb[9] > 1.0e-30 && Xa[9] == 1.0e-300                          # helium excited levels vanish
        @test length(pack_ysol(Xb, nH, nHe)) == 12 && length(pack_ysol(Xa, nH, nHe; helium = false)) == 7    # state dimension changes
    end

    @testset "shape and invalid arguments" begin
        @test_throws DimensionMismatch helium_switch_reset(zeros(14), 6, 7, 0.08)
        @test_throws DimensionMismatch pack_ysol(zeros(14), 6, 7; helium = false)
        @test_throws BoundsError pack_ysol(zeros(15), 6, 7, (1, 2, 3, 4, 6); helium = false)
        @test helium_switch_ysize(5, 4, true) == 12 && helium_switch_ysize(5, 4, false) == 7
        @test eltype(helium_switch_reset(zeros(15), 6, 7, 0.08)) == Float64
        Xin = copy(FX4D[1].Xpre); Xin0 = copy(Xin); helium_switch_reset(Xin, 6, 7, 0.08); @test Xin == Xin0      # the reset must not mutate its input
        @inferred helium_switch_reset(zeros(15), 6, 7, 0.08)
        @inferred pack_ysol(zeros(15), 6, 7; helium = false)
    end
end
