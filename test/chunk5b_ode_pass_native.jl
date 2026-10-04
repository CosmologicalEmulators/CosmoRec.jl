# Chunk 5b: ONE native recombination pass (runmode 1, no diffusion): Julia SciML pass (12-state ODE, sampled-HeI switch, 7-state ODE) vs the native internal arrays
# (`pass_on_the_Solution_CosmoRec`, `output_CosmoRec` nodes, `HeI_was_switched_off_at_z`). REQUIRES COSMOREC_NATIVE_DATA_DIR (see test/chunk5_helpers.jl); fails loudly if unset.
# The native solver controls only rho, X1s, H 2s, 2p (rel 1e-7/1e-7/1e-6/1e-6, abs 1e-12/1e-12/1e-80/1e-80) and He 1s: differences are reported in these native-tolerance units.
using Test
using CosmoRec
using SHA
include("chunk5_helpers.jl")

const OBS5B = Dict{String,Float64}()
rec5b(k, v) = (OBS5B[k] = max(get(OBS5B, k, 0.0), v))
# a1 = 1e-18 (was 1e-16, approved 2026-10-02): at a1 = 1e-16 the X1s weight a1 + reltol|X1s| is absolute-dominated for z > 1813 (X1s ~ 7e-10 at z = 3000), docs/ORIGINAL_VS_PORT_ACCURACY_INVESTIGATION.md 5.6
const SOLVE5B = (a...) -> solve_phase5(a...; reltol = 1.0e-12, a1 = 1.0e-18, aex = 1.0e-14)
const PASS5B = recombination_pass(RM5, SOLVE5B)

@testset "Chunk 5b: one native pass vs the SciML pass" begin
    @testset "fixture provenance" begin
        for (k, p) in FIX5B
            @test bytes2hex(open(sha256, p)) == FIX5B_SHA[k]
        end
        @test size(NODES5) == (3000, 11) && size(OUT5) == (3199, 3) && size(GRID5) == (10000, 3)
        @test SUM5["NODES"] == 3000 && SUM5["OUTROWS"] == 3199
        @test SUM5["FLAG_HE_after_run"] == 1             # native resets flag_He after the run
    end

    @testset "node grid and initial state" begin
        @test PASS5B.z == NODES5[:, 1]                    # bitwise, including the accumulated last node 50.000000000188152
        @test PASS5B.z[end] == 50.000000000188152 && PASS5B.z[1] == 3000.0
        # initial state = native Saha (4a/4b), node 0
        for c in 1:8
            @test isapprox(PASS5B.observables[1][c], NODES5[1, 1 + c]; rtol = 1e-13)
        end
    end

    @testset "reference callback: FullSpecialize (no FunctionWrappersWrapper), per-solve primal RHS workspace and buffered Jacobian" begin
        # caller-side solver configuration only (the library has no solver default); AutoSpecialize would wrap recombination_ode! at runtime
        @test ODEFUN5 isa ODEFunction{true, SciMLBase.FullSpecialize}
        u0 = PASS5B.states[1]; zn = [PASS5B.z[1], PASS5B.z[2]]
        prob = prob_phase5(u0, zn[1], zn, RecombinationODEParams(RM5, true))     # the problem solve_phase5 actually builds
        @test prob isa ODEProblem && prob.f isa ODEFunction{true, SciMLBase.FullSpecialize}
        integ = SciMLBase.init(prob, Rodas5P(); reltol = 1.0e-12, abstol = abstol5(12; a1 = 1.0e-18, aex = 1.0e-14), saveat = zn, internalnorm = primal_norm)
        @test integ.f.f isa RHSWS5                                               # the per-solve primal workspace functor itself (no FunctionWrappersWrapper)
        @test integ.f isa ODEFunction{true, SciMLBase.FullSpecialize}
        # per-solve buffered state Jacobian + primal workspace: a real 12-state block solves bitwise like the allocating recombination_ode!/jac5! callback (ODEFUN5)
        @test integ.f.jac isa WSJac5 && prob_phase5(u0, zn[1], zn, RecombinationODEParams(RM5, true)).f.jac !== prob.f.jac   # buffered Jacobian, fresh per problem
        zb = PASS5B.z[1:40]; kw = (; reltol = 1.0e-12, abstol = abstol5(12; a1 = 1.0e-18, aex = 1.0e-14), saveat = zb, internalnorm = primal_norm)
        pb = prob_phase5(u0, zb[1], zb, RecombinationODEParams(RM5, true))
        sp = solve(pb, Rodas5P(); kw...)
        so = solve(ODEProblem{true, SciMLBase.FullSpecialize}(ODEFUN5, u0, (zb[1], zb[end]), RecombinationODEParams(RM5, true)), Rodas5P(); kw...)
        @test wsjac5_fallbacks(pb.f.jac) == 0 && pb.f.f.nfallback[] == 0          # every RHS and Jacobian call used its prepared workspace/cache
        @test pb.f.f !== prob_phase5(u0, zb[1], zb, RecombinationODEParams(RM5, true)).f.f   # a fresh primal workspace per problem
        @test sp.u == so.u && sp.stats.nf == so.stats.nf && sp.stats.njacs == so.stats.njacs && sp.stats.naccept == so.stats.naccept
    end

    @testset "helium switch node and redshift reproduce the native decision" begin
        @test PASS5B.zswitch == SUM5["HE_OFF"] == 1680.9103034346122
        @test PASS5B.k_switch == 1342
        @test length(PASS5B.states[1342]) == 7 || length(PASS5B.states[1342]) == 12
        @test all(m -> length(PASS5B.states[m]) == 12, 1:1342) && all(m -> length(PASS5B.states[m]) == 7, 1343:3000)
        # frozen branch (no event detection; the last 12-state block ends at the switch node instead of overshooting it) reproduces the detected run to solver accuracy
        pf = recombination_pass(RM5, SOLVE5B; k_switch = PASS5B.k_switch)
        @test pf.k_switch == PASS5B.k_switch
        d = maximum(abs.(pf.Xe .- PASS5B.Xe) ./ PASS5B.Xe); rec5b("frozen_vs_detected:Xe", d); @test d < 1e-8
        d = maximum(abs.(pf.Te .- PASS5B.Te) ./ PASS5B.Te); rec5b("frozen_vs_detected:Te", d); @test d < 1e-8
    end

    @testset "nodewise comparison with the native pass_on_the_Solution arrays" begin
        obs = reduce(hcat, PASS5B.observables)'               # 3000 x 8: Xe, X1s, 2s, 2p, 3s, 3p, 3d, rho
        nat = NODES5[:, 2:9]
        names = ("Xe", "X1s", "2s", "2p", "3s", "3p", "3d", "rho")
        native_tol = Dict("X1s" => (1e-7, 1e-12), "rho" => (1e-7, 1e-12), "2s" => (1e-6, 1e-80), "2p" => (1e-6, 1e-80))
        for (c, nm) in enumerate(names)
            rel = abs.(obs[:, c] .- nat[:, c]) ./ abs.(nat[:, c])
            rec5b("rel:" * nm, maximum(rel)); rec5b("rel_median:" * nm, sort(rel)[length(rel) ÷ 2]); rec5b("rel_worst_z:" * nm, PASS5B.z[argmax(rel)])
            if haskey(native_tol, nm)
                u = abs.(obs[:, c] .- nat[:, c]) ./ (native_tol[nm][1] .* abs.(nat[:, c]) .+ native_tol[nm][2])
                rec5b("units:" * nm, maximum(u)); rec5b("units_median:" * nm, sort(u)[length(u) ÷ 2])
                @test maximum(u) < (nm in ("X1s", "rho") ? 1.0 : 10.0)
            end
        end
        # Xe has no native tolerance of its own (X[0] is recomputed): compared in units of the native 1e-7 relative scale
        rec5b("units1e-7:Xe", maximum(abs.(obs[:, 1] .- nat[:, 1]) ./ (1e-7 .* abs.(nat[:, 1]))))
        @test maximum(abs.(obs[:, 1] .- nat[:, 1]) ./ abs.(nat[:, 1])) < 5e-7
        # unmonitored components (3s, 3p, 3d): reported, loose gate (native itself does not control them)
        for nm in ("3s", "3p", "3d")
            c = findfirst(==(nm), names)
            @test maximum(abs.(obs[:, c] .- nat[:, c]) ./ abs.(nat[:, c])) < 1e-3
        end
        # matter temperature Te = rho TCMB vs the native output_CosmoRec ODE rows
        Tn = OUT5[1:3000, 3]
        rec5b("rel:Te", maximum(abs.(PASS5B.Te .- Tn) ./ Tn)); @test maximum(abs.(PASS5B.Te .- Tn) ./ Tn) < 1e-8
        @test OUT5[1:3000, 1] == PASS5B.z && maximum(abs.(OUT5[1:3000, 2] .- NODES5[:, 2])) == 0.0     # the two native arrays agree on Xe
        # after the switch the electron fraction is the proton fraction (flag_He = 0)
        for m in (1343, 2000, 3000); @test PASS5B.Xe[m] == 1 - PASS5B.observables[m][2]; end
        @test PASS5B.Xe[1342] > 1 - PASS5B.observables[1342][2]        # the switch node itself still carries XHeII (pre-switch state)
    end

    @testset "solver self-consistency (tolerance halving)" begin
        p2 = recombination_pass(RM5, (a...) -> solve_phase5(a...; reltol = 1.0e-13, a1 = 1.0e-17, aex = 1.0e-16))
        @test p2.k_switch == PASS5B.k_switch
        d = maximum(abs.(p2.Xe .- PASS5B.Xe) ./ PASS5B.Xe); rec5b("self:Xe", d); @test d < 1e-8
        d = maximum(abs.(p2.Te .- PASS5B.Te) ./ PASS5B.Te); rec5b("self:Te", d); @test d < 1e-8
        obs1 = reduce(hcat, PASS5B.observables)'; obs2 = reduce(hcat, p2.observables)'
        d = maximum(abs.(obs2[:, 4] .- obs1[:, 4]) ./ obs1[:, 4]); rec5b("self:2p", d); @test d < 1e-4
    end

    @testset "info" begin
        foreach(kv -> println("Chunk5b-native ", kv[1], " = ", kv[2]), sort!(collect(OBS5B)))
    end
end
