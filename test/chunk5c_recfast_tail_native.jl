# Chunk 5c: the Recfast tail (z = 50 -> 0.001, 200 log nodes) of the native pass: native `compute_Recfast_part` inputs, `Xe_frac_rescaled` initial condition and rescaling, and the native
# `output_CosmoRec` tail rows. REQUIRES COSMOREC_NATIVE_DATA_DIR; Recfast uses the loaded Cosmos H(z) (including the reproduced native low-z bug, relevant for z < 1).
# The native Recfast solver (rel 1e-8, abs (1e-10, 1e-8, 1e-10)) is the reference error scale (as in Chunk 4c).
using Test
using CosmoRec
using SHA
include("chunk5_helpers.jl")

const OBS5C = Dict{String,Float64}()
rec5c(k, v) = (OBS5C[k] = max(get(OBS5C, k, 0.0), v))
# a1 = 1e-18 (was 1e-16, approved 2026-10-02): at a1 = 1e-16 the X1s weight a1 + reltol|X1s| is absolute-dominated for z > 1813 (X1s ~ 7e-10 at z = 3000), docs/ORIGINAL_VS_PORT_ACCURACY_INVESTIGATION.md 5.6
const SOLVE5C = (a...) -> solve_phase5(a...; reltol = 1.0e-12, a1 = 1.0e-18, aex = 1.0e-14)
const PASS5C = recombination_pass(RM5, SOLVE5C)
const TAIL_NATIVE = OUT5[3001:3199, :]              # z, Xe, Te for tail nodes i = 1..199
const NATIVE_ATOL5C = (1.0e-10, 1.0e-8, 1.0e-10)

@testset "Chunk 5c: Recfast tail vs native" begin
    θ = theta5(); Hf = hfun5()
    ri = SUM5["RECFAST_INPUT"]
    @testset "native compute_Recfast_part inputs from the native final state" begin
        Xn = SUM5["FINAL_X"]                                                 # native Level_I.X at the end of the pass (flag_He restored to 1 by the native code; helium is off here)
        y7 = pack_ysol(Xn, 6, 7; helium = false)
        inp = recfast_tail_inputs(RM5, ri["ze"], y7)
        @test isapprox(inp.Xe_Hi, ri["Xe_Hi"]; rtol = 1e-14) && isapprox(inp.Xei, ri["Xei"]; rtol = 1e-14)
        @test inp.TMi == cosmos_TCMB(ACC5, ri["ze"]) * ri["rhoi"] || isapprox(inp.TMi, cosmos_TCMB(ACC5, ri["ze"]) * ri["rhoi"]; rtol = 1e-15)
        e = abs(inp.dXei - ri["dXe_dz"]) / abs(ri["dXe_dz"]); rec5c("input:dXe_dz", e); @test e < 1e-3       # dX1s/dt at z = 50 is a difference of large rates (1.8e-20 net): conditioning-limited (observed 2.1e-4); its downstream effect is tested below
        @test inp.Xe_Hei == 1.0e-30
        @test ri["ze"] == PASS5C.z[end]
    end

    @testset "tail grid" begin
        zt = init_xarr_log(ri["ze"], 0.001, 200)
        @test zt[1] == ri["ze"] && length(zt) == 200
        e = maximum(abs.(zt[2:end] .- TAIL_NATIVE[:, 1]) ./ TAIL_NATIVE[:, 1]); rec5c("grid", e)
        @test e < 1e-13
        @test isapprox(zt[end], 0.001; rtol = 1e-12)
        @test issorted(zt; rev = true)
    end

    @testset "tail from the native inputs vs the native tail rows (isolates the solver)" begin
        Xn = SUM5["FINAL_X"]; y7 = pack_ysol(Xn, 6, 7; helium = false)
        inp = recfast_tail_inputs(RM5, ri["ze"], y7)
        tl = recfast_tail(θ, RFCONST5, Hf, ri["ze"], inp.Xe_Hi, inp.Xe_Hei, inp.Xei, inp.TMi, ri["dXe_dz"], solve_tail5)
        # native-tolerance units |Delta| / (1e-8 |y| + atol) on the monitored Recfast components
        uXe = abs.(tl.Xe[2:end] .- TAIL_NATIVE[:, 2]) ./ (1e-8 .* abs.(TAIL_NATIVE[:, 2]) .+ NATIVE_ATOL5C[2]); rec5c("units:Xe", maximum(uXe))
        uTe = abs.(tl.TM[2:end] .- TAIL_NATIVE[:, 3]) ./ (1e-8 .* abs.(TAIL_NATIVE[:, 3]) .+ NATIVE_ATOL5C[3]); rec5c("units:Te", maximum(uTe))
        rec5c("rel:Xe", maximum(abs.(tl.Xe[2:end] .- TAIL_NATIVE[:, 2]) ./ TAIL_NATIVE[:, 2])); rec5c("rel:Te", maximum(abs.(tl.TM[2:end] .- TAIL_NATIVE[:, 3]) ./ TAIL_NATIVE[:, 3]))
        rec5c("worst_z:Xe", TAIL_NATIVE[argmax(uXe), 1]); rec5c("worst_z:Te", TAIL_NATIVE[argmax(uTe), 1]); rec5c("ff", tl.ff)
        @test maximum(uXe) < 5.0 && maximum(uTe) < 500.0
        @test isfinite(tl.ff) && tl.ff > 0
        # rescaling: the rescaled xp equation reproduces the native derivative at the start
        f0 = recfast_rhs(θ, ri["ze"], (1.0e-30, inp.Xe_Hi, inp.TMi, 0.0), Hf(ri["ze"]), RFCONST5)
        @test isapprox(tl.ff * f0[2], ri["dXe_dz"]; rtol = 1e-12)
    end

    @testset "end-to-end: tail from the Julia pass" begin
        zi = PASS5C.z[end]; y7 = PASS5C.states[end]
        inp = recfast_tail_inputs(RM5, zi, y7)
        tl = recfast_tail(θ, RFCONST5, Hf, zi, inp.Xe_Hi, inp.Xe_Hei, inp.Xei, inp.TMi, inp.dXei, solve_tail5)
        rXe = maximum(abs.(tl.Xe[2:end] .- TAIL_NATIVE[:, 2]) ./ TAIL_NATIVE[:, 2]); rec5c("e2e:rel:Xe", rXe); @test rXe < 1e-3
        rTe = maximum(abs.(tl.TM[2:end] .- TAIL_NATIVE[:, 3]) ./ TAIL_NATIVE[:, 3]); rec5c("e2e:rel:Te", rTe); @test rTe < 1e-3
        @test tl.Xe[1] == inp.Xei && tl.TM[1] == inp.TMi
    end

    @testset "low-z native Hubble bug acts on the tail" begin
        # below z = 1 the loaded H spline carries the native lHz[0] error: the tail uses it (native behaviour reproduced), not the analytic H
        Han = z -> cosmos_H(accessors4b(hubble = false), z)
        @test abs(Hf(0.01) / Han(0.01) - 1) > 0.1
        @test abs(Hf(10.0) / Han(10.0) - 1) < 1e-2        # CAMB table vs the Seager analytic form: 2e-3 at z = 10
    end

    @testset "info" begin
        foreach(kv -> println("Chunk5c-native ", kv[1], " = ", kv[2]), sort!(collect(OBS5C)))
    end
end
