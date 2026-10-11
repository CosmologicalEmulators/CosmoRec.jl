# Native parity milestone 1: faithful port of the ORIGINAL Recfast++ Gear ODE solver (src/RecfastGear.jl, ODE_solver.Recfast.cpp) driving the preliminary
# history, with the standalone closure constants (StandaloneInit.jl) and the native-order pivot redshifts (RecfastPP.jl recfast_grid, Recfast++.cpp:151).
# References: verbatim native text fixtures native_recfast_history.txt + native_hubble_input.txt (fid, 6000 nodes) and native_standalone_h50_* (h50).
# Milestone 2: native-order RHS (evaluate_recfast_system_native, evalode.Recfast.cpp:99-220) + C-library exp/log/pow in the Float64 Hubble/Saha/NH/Te_QS paths.
# Milestone 3: Hubble-table natural spline with the GSL 2.8 solve_tridiag operation order (CosmosAccessors.jl natural_cubic_spline); direct test vs native GSL
# coefficients/evaluations (native_gsl_hubble_spline_{fid,h50}.txt). Exact (bitwise) facts only.
# No solver package needed (pure Float64 port).
using Test
using CosmoRec
using Printf

const FIXP1 = joinpath(@__DIR__, "fixtures")
rowsP1(p) = [parse.(Float64, split(l)) for l in eachline(p) if !isempty(l) && !startswith(l, "#")]
function gearP1(H0, ombh2, omch2, Y, nnu, hfile)
    d = camb_density_parameters(ombh2, omch2, H0)
    c = cosmos_constants(; Omega_c = d.Omega_c, Omega_b = d.Omega_b, Omega_k = 0.0, N_eff = nnu, h100 = H0 / 100.0, T_CMB = 2.7255, Y_p = Y)
    θ = recfast_theta(; Omega_b = d.Omega_b, h100 = H0 / 100.0, T_CMB = 2.7255, Y_p = Y)
    HB = reduce(hcat, rowsP1(hfile))'
    hacc = CosmosAccessors(c, nothing, hubble_table(HB[:, 1], HB[:, 2])); Hf = z -> cosmos_H(hacc, z)
    g = recfast_grid(θ, NATIVE_RECFAST_CONSTANTS, Hf)
    return g, recfast_history(θ, NATIVE_RECFAST_CONSTANTS, Hf, g, recfast_gear_solve)
end
relP1(a, b) = maximum(i -> a[i] == b[i] ? 0.0 : abs(a[i] - b[i]) / abs(b[i]), eachindex(a))

@testset "Native parity 1: Recfast++ Gear solver port (preliminary history)" begin
    @testset "fiducial: grid, Saha segments, whole Gear preliminary history vs native (bitwise)" begin
        g, h = gearP1(67.36, 0.02237, 0.12, 0.24568275240335613, 3.046, joinpath(FIXP1, "native_hubble_input.txt"))
        M = reduce(hcat, rowsP1(joinpath(FIXP1, "native_recfast_history.txt")))'
        @test g.z == M[:, 1]                                   # node grid bitwise (native-order pivot redshifts, Recfast++.cpp:151)
        j = g.jode
        @test h.Xe_H[1:j] == M[1:j, 2] && h.Xe_He[1:j] == M[1:j, 3] && h.TM[1:j] == M[1:j, 7]
        # whole preliminary history (Gear solver + native-order RHS evaluate_recfast_system_native + C-library exp/log/pow): ALL columns bitwise
        @test h.Xe_H == M[:, 2] && h.Xe_He == M[:, 3] && h.Xe == M[:, 4] && h.TM == M[:, 7]
        @test h.dXe == M[:, 5] && h.dXe_H == M[:, 6]            # native derivative output (Sz.dy at the accepted state)
        for (nm, col, v) in (("Xe", 4, h.Xe), ("TM", 7, h.TM), ("dXe", 5, h.dXe))
            @printf("NativeParity fid %-3s bitwise %d/6000, max rel %.3e\n", nm, count(v .== M[:, col]), relP1(v, M[:, col]))
        end
    end
    @testset "h50 (prior corner): grid, Saha, whole preliminary history bitwise" begin
        g, h = gearP1(50.0, 0.025, 0.16, 0.245, 3.046, joinpath(FIXP1, "native_standalone_h50_hubble.txt"))
        R = rowsP1(joinpath(FIXP1, "native_standalone_h50_history.txt")); idx = Int.(getindex.(R, 1)); M = reduce(hcat, [r[2:end] for r in R])'
        @test g.z[idx] == M[:, 1]
        sa = idx .<= g.jode
        @test h.Xe_H[idx[sa]] == M[sa, 2] && h.Xe_He[idx[sa]] == M[sa, 3] && h.TM[idx[sa]] == M[sa, 7]
        @test h.Xe_H[idx] == M[:, 2] && h.Xe_He[idx] == M[:, 3] && h.Xe[idx] == M[:, 4] && h.TM[idx] == M[:, 7]     # all states bitwise
        # derivative output: bitwise since the GSL-order spline (milestone 3; previously 25 nodes with z < 1.6 differed through 131 non-bitwise H values)
        @test h.dXe[idx] == M[:, 5] && h.dXe_H[idx] == M[:, 6]
    end
    @testset "Hubble table spline vs native GSL cspline (coefficients and evaluation, bitwise)" begin
        for (nm, hfile) in (("fid", "native_hubble_input.txt"), ("h50", "native_standalone_h50_hubble.txt"))
            HB = reduce(hcat, rowsP1(joinpath(FIXP1, hfile)))'
            ht = hubble_table(HB[:, 1], HB[:, 2]); s = ht.spline
            L = [split(l, '\t') for l in eachline(joinpath(FIXP1, "native_gsl_hubble_spline_$nm.txt")) if !startswith(l, "#")]
            K = [parse.(Float64, l[3:5]) for l in L if l[1] == "KNOT"]; E = [parse.(Float64, l[2:3]) for l in L if l[1] == "EVAL"]
            @test length(K) == 10000 && length(E) == 5905
            @test s.x == getindex.(K, 1) && s.y == getindex.(K, 2)        # ln z, ln H knots (node-0 quirk ln H = 0 included)
            @test s.c == getindex.(K, 3)                                  # second-derivative coefficients (GSL solve_tridiag order)
            hacc = CosmosAccessors(cosmos_constants(; Omega_c = 0.26, Omega_b = 0.05, Omega_k = 0.0, N_eff = 3.046, h100 = 0.7, T_CMB = 2.7255, Y_p = 0.245), nothing, ht)
            @test all(e -> cosmos_H(hacc, e[1]) == e[2], E)               # table branch only (zmin < z < zmax); closure constants unused there
        end
    end
end
