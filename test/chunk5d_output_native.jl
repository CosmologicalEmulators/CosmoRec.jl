# Chunk 5d: output assembly and the returned history on the CAMB grid (10000 nodes, z = 1e4 -> 0) of the native pass (runmode 1): `output_CosmoRec` rows, GSL splines of ln Xe and ln Te,
# and the Recfast-history fallback outside the stored range, vs the native `Xe_arr`, `Te_arr`. REQUIRES COSMOREC_NATIVE_DATA_DIR.
using Test
using CosmoRec
using SHA
include("chunk5_helpers.jl")

const OBS5D = Dict{String,Float64}()
rec5d(k, v) = (OBS5D[k] = max(get(OBS5D, k, 0.0), v))
# a1 = 1e-18 (was 1e-16, approved 2026-10-02): at a1 = 1e-16 the X1s weight a1 + reltol|X1s| is absolute-dominated for z > 1813 (X1s ~ 7e-10 at z = 3000), docs/ORIGINAL_VS_PORT_ACCURACY_INVESTIGATION.md 5.6
const SOLVE5D = (a...) -> solve_phase5(a...; reltol = 1.0e-12, a1 = 1.0e-18, aex = 1.0e-14, alg = Rodas5P(), tstops = false)
const HIST5D = recombination_history(RM5, theta5(), SOLVE5D, solve_tail5, GRID5[:, 1])

@testset "Chunk 5d: output assembly vs native Xe_arr, Te_arr" begin
    @testset "grid and stored rows" begin
        @test GRID5[1, 1] == 1.0e4 && GRID5[end, 1] == 0.0 && issorted(GRID5[:, 1]; rev = true)
        rows = HIST5D.rows
        @test length(rows.z) == 3199 && rows.z[1:3000] == OUT5[1:3000, 1]
        e = maximum(abs.(rows.z[3001:end] .- OUT5[3001:end, 1]) ./ OUT5[3001:end, 1]); rec5d("rows:z_tail", e); @test e < 1e-13
        @test issorted(rows.z; rev = true) && length(unique(rows.z)) == 3199
    end

    @testset "returned history on the CAMB grid" begin
        Xe, Te = HIST5D.Xe, HIST5D.Te
        z = GRID5[:, 1]
        inside = (z .< 3000.0) .& (z .> HIST5D.rows.z[end])
        @test count(inside) == 2999 && !inside[1] && !inside[end]       # z >= 3000 (7000 nodes) and z = 0 use the preliminary history
        # fallbacks (outside the stored range): the native Xe_Seager, Te of the Cosmos object
        out = .!inside
        e = maximum(abs.(Xe[out] .- GRID5[out, 2]) ./ GRID5[out, 2]); rec5d("fallback:Xe", e); @test e < 1e-13
        e = maximum(abs.(Te[out] .- GRID5[out, 3]) ./ GRID5[out, 3]); rec5d("fallback:Te", e); @test e < 1e-13
        # inside: spline of the Julia rows vs the native spline of the native rows (native solver tolerance level)
        rXe = abs.(Xe[inside] .- GRID5[inside, 2]) ./ GRID5[inside, 2]; rTe = abs.(Te[inside] .- GRID5[inside, 3]) ./ GRID5[inside, 3]
        rec5d("inside:Xe", maximum(rXe)); rec5d("inside:Te", maximum(rTe)); rec5d("inside:Xe_median", sort(rXe)[length(rXe) ÷ 2]); rec5d("inside:Te_median", sort(rTe)[length(rTe) ÷ 2])
        rec5d("inside:worst_z_Xe", z[inside][argmax(rXe)]); rec5d("inside:worst_z_Te", z[inside][argmax(rTe)])
        @test maximum(rXe) < 5e-6
        @test maximum(rTe) < 5e-6
    end

    @testset "spline stage alone: native rows in, native grid out" begin
        # the native stored rows through the Julia spline/fallback code must reproduce the native returned arrays to rounding (isolates the assembly from the ODE)
        rows = (z = OUT5[:, 1], Xe = OUT5[:, 2], Te = OUT5[:, 3])
        Xe, Te = return_solution_to_grid(RM5, rows, GRID5[:, 1])
        e1 = maximum(abs.(Xe .- GRID5[:, 2]) ./ GRID5[:, 2]); e2 = maximum(abs.(Te .- GRID5[:, 3]) ./ GRID5[:, 3])
        rec5d("splines_only:Xe", e1); rec5d("splines_only:Te", e2)
        @test e1 < 1e-12 && e2 < 1e-12
    end

    @testset "limiting cases" begin
        @test_throws DimensionMismatch natural_cubic_spline([0.0, 1.0], [1.0, 2.0])       # fewer than 3 nodes cannot define a spline
        rows = HIST5D.rows
        Xe0, Te0 = return_solution_to_grid(RM5, rows, [3000.0, 0.0, 5.0e3])
        @test Xe0[1] == cosmos_Xe_Seager(ACC5, 3000.0) && Xe0[2] == cosmos_Xe_Seager(ACC5, 0.0) && Te0[3] == cosmos_Te(ACC5, 5.0e3)
        @test Xe0[1] > 1.0   # preliminary history at z = 3000 (helium still ionized)
    end

    @testset "info" begin
        foreach(kv -> println("Chunk5d-native ", kv[1], " = ", kv[2]), sort!(collect(OBS5D)))
    end
end
