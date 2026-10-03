# Chunk 4c: preliminary Recfast++ history vs the original native outputs: RHS / Saha / Te_QS probes (no solver involved), node grid, Saha segments, ODE nodes.
# The native solver (variable-order Gear, rel tol 1e-8, abs tol (1e-10, 1e-8, 1e-10) for (xHep, xp, TM)) is NOT reproduced step by step: the ODE part is compared with a
# tight-tolerance SciML (Rodas5P) solution of the same ODE; the nodewise differences measure the native solver error. Scope: docs/ACCEPTANCE_CHUNK4C.md.
using Test
using CosmoRec
using SHA
@isdefined(FX4C) || include("chunk4c_helpers.jl")

const OBS4C = Dict{String,Float64}()
rec4c(k, v) = (OBS4C[k] = max(get(OBS4C, k, 0.0), v))
rel4c(a, b) = a == b ? 0.0 : abs(a - b) / max(abs(b), floatmin())
const TOL_UNITS4C = Dict("Xe_He" => 1.0, "Xe_H" => 2.0, "TM" => 100.0, "Xe" => 2.0, "dXe" => 0.1, "dXe_H" => 0.1)   # observed 0.20, 1.04, 48.9, 1.04, 0.034, 0.034 (docs/CHUNK4C_RESULTS.md)
const NATIVE_ATOL4C = (1.0e-10, 1.0e-8, 1.0e-10)      # native abs_vector (ODE_solver.Recfast.cpp:403-405)
const NATIVE_RTOL4C = 1.0e-8                          # native tolSol

@testset "Chunk 4c: Recfast++ preliminary history vs native" begin
    θ = theta4c(); c = NATIVE_RECFAST_CONSTANTS; Hf = hfun4c()
    @testset "fixture provenance and active-term switches" begin
        @test bytes2hex(open(sha256, FIX4C)) == FIX4C_SHA256
        v = FX4C["RFV"]
        @test v["eval_ion_Tr"] == 0 && v["reion_model"] == 0 && v["f_dec"] == 0 && v["f_ann"] == 0 && v["B0"] == 0 && v["include_CF"] == 0
        @test v["aS"] == 1 && v["mS"] == 1 && v["pS"] == 0
        @test FX4C["RFC"]["A2s1s_rf"] == 8.2205999999999992          # observed native value (the adapter passes 0): not interpreted
        @test length(FX4C["RFR"]) == 900 && length(FX4C["RFS"]) == 36 && length(FX4C["RFH"]) == 7
        @test FX4C["RFC"]["fHe_rf"] == rf_fHe(c, θ[3])
    end

    @testset "native RHS at 900 states (no solver)" begin
        # f1, f2 are differences of two large nearly equal terms at equilibrium (Saha-segment) nodes: the raw relative error is amplified by that cancellation
        # (28 of 900 states exceed 1e-12 raw, worst 1.7e-6, all at native node states). The meaningful metric is the error normalised by the term magnitudes (C factors <= 1 omitted).
        for r in FX4C["RFR"]
            z = r[3]; y = r[4:7]
            @test isapprox(Hf(z), r[8]; rtol = 1e-14)
            @test isapprox(rf_NH(c, θ, z), r[9]; rtol = 1e-14)
            @test rf_TCMB(θ, z) == r[10]
            f = recfast_rhs(θ, z, y, r[8], c)       # native Hz fed back: isolates the RHS
            nH = r[9]; xe = y[1] + y[2]; fHe = rf_fHe(c, θ[3]); pref = 1 / ((1 + z) * r[8])
            sc1 = pref * (abs(rf_alphaHe(y[3]) * xe * y[1] * nH) + abs(rf_alphaHe(y[3]) / rf_sahaboltz(1.0, 2.0, 1.0, c.EionHe2s, y[3], c) * (fHe - y[1]) * rf_boltzmann(1.0, 1.0, c.E2s1sHe, y[3], c)))
            sc2 = pref * (abs(rf_alphaH(θ[1], y[3]) * xe * y[2] * nH) + abs(rf_alphaH(θ[1], y[3]) / rf_sahaboltz(2.0, 1.0, 1.0, c.EionH2s, y[3], c) * (1.0 - y[2]) * rf_boltzmann(1.0, 1.0, c.E2s1sH, y[3], c)))
            sc = (sc1, sc2)
            for k in 1:4
                rec4c("rhs$(k)_raw", rel4c(f[k], r[10 + k]))
                en = k <= 2 ? abs(f[k] - r[10 + k]) / sc[k] : rel4c(f[k], r[10 + k])
                rec4c("rhs$(k)_scaled", en)
                @test en < 1e-12
            end
        end
    end

    @testset "Saha helpers, Te_QS, H and t_cos" begin
        for r in FX4C["RFS"]
            z, Xe = r[1], r[2]
            Hz = Hf(z)
            @test isapprox(rf_Te_QS(θ, z, Xe, Hz, c), r[3]; rtol = 1e-12)
            @test isapprox(rf_SahaBoltz_HeIII(r[6], rf_fHe(c, θ[3]), 0.97 * r[7], c), r[4]; rtol = 1e-12)
            @test isapprox(rf_SahaBoltz_HeII(r[6], rf_fHe(c, θ[3]), 0.97 * r[7], c), r[5]; rtol = 1e-12)
            rec4c("TeQS", rel4c(rf_Te_QS(θ, z, Xe, Hz, c), r[3])); rec4c("SahaHeIII", rel4c(rf_SahaBoltz_HeIII(r[6], rf_fHe(c, θ[3]), 0.97 * r[7], c), r[4]))
            rec4c("SahaHeII", rel4c(rf_SahaBoltz_HeII(r[6], rf_fHe(c, θ[3]), 0.97 * r[7], c), r[5]))
        end
        for r in FX4C["RFH"]
            @test isapprox(Hf(r[1]), r[2]; rtol = 1e-13) || r[1] < 1.0     # z < 1: native loaded-H bug region reproduced by the same spline (checked in 4b)
        end
    end

    @testset "node grid and Saha segments" begin
        grid = recfast_grid(θ, c, Hf)
        zn = HIST4B[:, 1]
        @test length(grid.z) == 6000 && grid.z[1] == 25000.0 && grid.z[end] == 0.0
        e = maximum(abs.(grid.z .- zn) ./ max.(zn, 1.0)); rec4c("grid", e)
        @test e < 1e-12
        @test issorted(grid.z; rev = true)
        j = grid.jode
        @test all(HIST4B[1:j, 2] .== 1.0) && all(HIST4B[1:j, 3] .== HIST4B[1, 3]) && all(HIST4B[1:j, 5] .== 0.0) && all(HIST4B[1:j, 6] .== 0.0)
        @test HIST4B[j, 4] == 1.0 + HIST4B[1, 3] && HIST4B[j - 1, 4] != HIST4B[j, 4]
        @test grid.seg_end[1] <= grid.seg_end[2] <= grid.seg_end[3] == j
        @test HIST4B[grid.seg_end[1] - 1, 4] == 1.0 + 2 * HIST4B[1, 3]            # last pure-ionized node
        rec4c("jode", Float64(j))
        # the grid is a discrete operator: a 1e-12 relative change of θ can move the Saha-segment break by a node (observed 2033 -> 2034)
        jp = recfast_grid(θ .* (1 + 1e-12), c, Hf).jode
        rec4c("jode_shift_1e-12", Float64(jp - j)); @test abs(jp - j) <= 1
        @test recfast_grid(θ, c, Hf).z == grid.z
    end

    @testset "full history: Saha segments exact, ODE nodes vs native solver" begin
        grid = recfast_grid(θ, c, Hf)
        h = recfast_history(θ, c, Hf, grid, (a...) -> solve_rodas_dense(a...; reltol = 1e-12, abstol = [1e-16, 1e-16, 1e-12]))
        j = grid.jode
        ode = (j + 1):6000
        @test all(isfinite, h.Xe) && all(isfinite, h.TM)
        for (name, col) in (("Xe_H", 2), ("Xe_He", 3), ("Xe", 4), ("dXe", 5), ("dXe_H", 6), ("TM", 7))
            vj = getproperty(h, Symbol(name)); ref = HIST4B[:, col]
            es = maximum(abs.(vj[1:j] .- ref[1:j]) ./ max.(abs.(ref[1:j]), name in ("dXe", "dXe_H") ? 1.0 : floatmin()))
            rec4c("saha:" * name, es)
            # Xe on Saha segment 2: d + sqrt(d^2 + ...) cancels (|d| ~ 4e5 at the worst node): 1 ulp of the libm pow/exp differences is amplified to 2.5e-11; gate = a few ulp(|d|)
            name == "Xe" || @test es < 1e-14
            rec4c("ode_raw:" * name, maximum(abs.(vj[ode] .- ref[ode]) ./ max.(abs.(ref[ode]), floatmin())))
        end
        fHe = rf_fHe(c, θ[3])
        for k in (grid.seg_end[1]):(grid.seg_end[2] - 1)
            g = begin
                T = rf_TCMB(θ, HIST4B[k, 1]); c1 = 2.0 * c.pi * c.mElect * c.kBoltz * T / (c.hPlanck * c.hPlanck)
                c1^1.5 * exp(-c.EionHeII / (c.kBoltz * T)) / rf_NH(c, θ, HIST4B[k, 1])
            end
            @test abs(h.Xe[k] - HIST4B[k, 4]) <= 8 * eps() * max(abs(1.0 + fHe - g), 1.0)
        end
        for k in 1:(grid.seg_end[1] - 1); @test h.Xe[k] == HIST4B[k, 4]; end
        # native-tolerance units: |Δ| / (rtol |y| + atol) with the native solver's own tolerances
        ux(vj, ref, atol) = abs.(vj[ode] .- ref[ode]) ./ (NATIVE_RTOL4C .* abs.(ref[ode]) .+ atol)
        uHe = ux(h.Xe_He, HIST4B[:, 3], NATIVE_ATOL4C[1]); uH = ux(h.Xe_H, HIST4B[:, 2], NATIVE_ATOL4C[2]); uT = ux(h.TM, HIST4B[:, 7], NATIVE_ATOL4C[3])
        for (nm, u) in (("Xe_He", uHe), ("Xe_H", uH), ("TM", uT))
            rec4c("ode_units:" * nm, maximum(u)); rec4c("ode_units_median:" * nm, sort(u)[length(u) ÷ 2])
            rec4c("ode_worst_z:" * nm, HIST4B[j + argmax(u), 1])
            @test maximum(u) < TOL_UNITS4C[nm]
        end
        uXe = ux(h.Xe, HIST4B[:, 4], NATIVE_ATOL4C[2]); rec4c("ode_units:Xe", maximum(uXe)); @test maximum(uXe) < TOL_UNITS4C["Xe"]
        for (nm, col) in (("dXe", 5), ("dXe_H", 6))
            vj = getproperty(h, Symbol(nm)); ref = HIST4B[:, col]
            e = abs.(vj[ode] .- ref[ode]) ./ maximum(abs.(ref[ode]))       # error in units of the peak |derivative|
            rec4c("ode_d_absmax:" * nm, maximum(e)); rec4c("ode_d_absmax_median:" * nm, sort(e)[length(e) ÷ 2]); @test maximum(e) < TOL_UNITS4C[nm]
        end
        # integration consistency: tightening the solver tolerance changes the history by far less than the native-vs-Julia difference
        h2 = recfast_history(θ, c, Hf, grid, solve_rodas_dense)      # reltol 1e-10, abstol (1e-14, 1e-14, 1e-10)
        d = maximum(abs.(h2.Xe[ode] .- h.Xe[ode]) ./ h.Xe[ode]); rec4c("solver_self_consistency:Xe(1e-10 vs 1e-12)", d); @test d < 5e-8
    end

    @testset "info" begin
        foreach(kv -> println("Chunk4c-native ", kv[1], " = ", kv[2]), sort!(collect(OBS4C)))
    end
end
