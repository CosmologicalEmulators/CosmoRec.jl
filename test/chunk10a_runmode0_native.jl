# Chunk 10a: the production CosmoRec iteration (runmode 0: 3 ODE passes, 2 HI diffusion PDE stages with feedback, Recfast tail, output assembly) vs the
# ORIGINAL cosmorec_calc_h_cpp_(runmode = 0) and the intermediates of its bitwise-verified replica. REQUIRES COSMOREC_NATIVE_DATA_DIR. Scope: docs/ACCEPTANCE_CHUNK10A.md.
using Test
using CosmoRec
using SHA
@isdefined(SETUP7) || include("chunk7_helpers.jl")

const FIX10A = joinpath(@__DIR__, "fixtures", "native_cosmorec_runmode0.txt")
const FIX10A_SHA256 = "f053f9c86f47e4ecaabe361f25eca325e3d05fb8b21b3ba036005b23723bfdca"
const OBS10A = Dict{String,Float64}()
rec10a(k, v) = (OBS10A[k] = max(get(OBS10A, k, 0.0), v))

function read_fixture10(path)
    vecs = Dict{String,Vector{Vector{Float64}}}(); mats = Dict{String,Vector{Vector{Float64}}}()
    for l in eachline(path)
        (isempty(l) || startswith(l, "#")) && continue
        f = split(l, '\t'); tag = String(f[1])
        if endswith(tag, "_ROW")
            push!(get!(mats, tag[1:(end - 4)], Vector{Vector{Float64}}()), parse.(Float64, f[2:end]))
        elseif all(x -> tryparse(Float64, x) !== nothing, f[2:end]) && length(f) > 3
            push!(get!(vecs, tag, Vector{Vector{Float64}}()), vec7(parse.(Float64, f[2:end])))
        end
    end
    return vecs, Dict(k => reduce(vcat, permutedims.(v)) for (k, v) in mats)
end
const V10, M10 = read_fixture10(FIX10A)
# a1 = 1e-18 (was 1e-16, approved 2026-10-02): at a1 = 1e-16 the X1s weight a1 + reltol|X1s| is absolute-dominated for z > 1813 (X1s ~ 7e-10 at z = 3000), docs/ORIGINAL_VS_PORT_ACCURACY_INVESTIGATION.md 5.6
const SOLVE10 = (a...) -> solve_phase5(a...; reltol = 1.0e-12, a1 = 1.0e-18, aex = 1.0e-14, alg = Rodas5P(), tstops = false)
const D10 = HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6)
const ZREC10 = V10["ZREC"][1]
const RUN10 = recombination_history_diffusion(RM5, theta5(), D10, SOLVE10, solve_tail5, ZREC10)
const NAMES10 = ("Xe", "X1s", "2s", "2p", "3s", "3p", "3d", "rho")
const NTOL10 = Dict("X1s" => (1e-7, 1e-12), "rho" => (1e-7, 1e-12), "2s" => (1e-6, 1e-80), "2p" => (1e-6, 1e-80))

function compare_pass10(tag, pass, nat)
    obs = reduce(hcat, pass.observables)'
    @test pass.z == nat[:, 1]
    for (c, nm) in enumerate(NAMES10)
        rel = abs.(obs[:, c] .- nat[:, 1 + c]) ./ abs.(nat[:, 1 + c])
        rec10a("$tag:rel:$nm", maximum(rel))
        haskey(NTOL10, nm) && rec10a("$tag:units:$nm", maximum(abs.(obs[:, c] .- nat[:, 1 + c]) ./ (NTOL10[nm][1] .* abs.(nat[:, 1 + c]) .+ NTOL10[nm][2])))
    end
end

@testset "Chunk 10a: production runmode-0 iteration vs native" begin
    @testset "fixture provenance and configuration" begin
        @test bytes2hex(open(sha256, FIX10A)) == FIX10A_SHA256
        cfg = read_kv7(FIX10A, "CFG10")
        @test (cfg["runmode"], cfg["Diffusion_correction"], cfg["DI1_2s_correction_on"], cfg["Diff_iteration_max"], cfg["nS_2gamma"], cfg["nS_Raman"], cfg["induced_flag"]) ==
              ("0", "1", "1", "2", "3", "2", "2")
        rep = read_kv7(FIX10A, "REPLICA")
        @test rep["bitwise_equal_xe"] == "1" && rep["bitwise_equal_tb"] == "1"
        @test V10["NAT_XE"][1] == V10["REP_XE"][1] && M10["NAT_OUTPUT"] == M10["REP_OUTPUT"]
        @test D10.zs == 2500.0 && D10.ze == 500.0 && D10.iterations == 2
        @test ZREC10 == GRID5[:, 1]
        # the first native pass of runmode 0 is the runmode-1 pass (Chunk 5b fixture)
        @test M10["PASS0"][:, 1:9] == NODES5[:, 1:9]
        @test M10["PASS2"] == M10["PASS1"]          # the last pass does not store pass_on_the_Solution_CosmoRec
        @test length(RUN10.passes) == 3 && length(RUN10.stages) == 2
    end
    Xn = V10["NAT_XE"][1]; Tn = V10["NAT_TB"][1]
    fbn10(k) = hi_diffusion_feedback(V10["PDE$(k)_z"][1], V10["PDE$(k)_DI1"][1], [V10["PDE$(k)_2g"][1], V10["PDE$(k)_2g"][2]], [V10["PDE$(k)_R"][1]])
    natpde10(k) = (V10["PDE$(k)_DI1"][1], V10["PDE$(k)_2g"][1], V10["PDE$(k)_2g"][2], V10["PDE$(k)_R"][1])
    jlpde10(o) = (o.DI1_2s, o.DF_2g[1], o.DF_2g[2], o.DF_R[1])
    CMP10 = ("DI1", "2g3s", "2g3d", "R2s")
    @testset "NATIVE-INJECTED DIAGNOSTIC: native PDE outputs -> Julia ODE passes, feedback, tail, assembly" begin
        p1 = recombination_pass(with_diffusion(RM5, fbn10(0)), SOLVE10)
        compare_pass10("inj_pass1", p1, M10["PASS1"])
        @test OBS10A["inj_pass1:units:X1s"] < 1.0 && OBS10A["inj_pass1:units:rho"] < 1.0
        @test OBS10A["inj_pass1:units:2s"] < 10.0 && OBS10A["inj_pass1:units:2p"] < 10.0 && OBS10A["inj_pass1:rel:Xe"] < 5e-7
        h = recombination_history(with_diffusion(RM5, fbn10(1)), theta5(), SOLVE10, solve_tail5, ZREC10)
        rec10a("inj_final:Xe", maximum(abs.(h.Xe .- Xn) ./ Xn)); rec10a("inj_final:Te", maximum(abs.(h.Te .- Tn) ./ Tn))
        @test OBS10A["inj_final:Xe"] < 5e-7 && OBS10A["inj_final:Te"] < 1e-7
    end
    @testset "NATIVE-INJECTED DIAGNOSTIC: native populations -> Julia PDE stages (fixed 8a integrator bound 1e-4 of the residual scale)" begin
        for k in 0:1
            o, _ = hi_diffusion_stage(RM5, D10, M10["PASS$k"][:, 1:9])
            @test o.z == V10["PDE$(k)_z"][1]
            nat = natpde10(k); jl = jlpde10(o)
            for q in 1:4
                rec10a("inj_pde$k:" * CMP10[q] * ":scaled", maximum(abs.(jl[q] .- nat[q]) ./ o.scale[q]))
                rec10a("inj_pde$k:" * CMP10[q] * ":rel_to_max", maximum(abs.(jl[q] .- nat[q])) / maximum(abs.(nat[q])))
                @test OBS10A["inj_pde$k:" * CMP10[q] * ":scaled"] <= 1.0e-4
            end
        end
        fx8 = read_fixture7(FIX8A)
        @test V10["PDE0_DI1"][1] == vec7(fx8["DI1_2s"][1]) && V10["PDE0_2g"][1] == vec7(fx8["DF_2g"][1]) && V10["PDE0_R"][1] == vec7(fx8["DF_R"][1])
    end
    @testset "FULL JULIA PATH: populations -> coefficients -> PDE -> feedback -> passes -> final history" begin
        compare_pass10("pass0", RUN10.passes[1], M10["PASS0"])
        compare_pass10("pass1", RUN10.passes[2], M10["PASS1"])
        @test OBS10A["pass0:units:X1s"] < 1.0 && OBS10A["pass0:units:rho"] < 1.0 && OBS10A["pass0:units:2s"] < 10.0 && OBS10A["pass0:units:2p"] < 10.0 && OBS10A["pass0:rel:Xe"] < 5e-7
        @test all(p -> p.k_switch == 1342 && p.zswitch == 1680.9103034346122, RUN10.passes)
        # PDE stages: every physical correction component, fixed bounds in the feedback range (500 < z < 2000; native interpolate_DF returns 0 outside)
        for k in 0:1
            o = RUN10.stages[k + 1]
            @test o.z == V10["PDE$(k)_z"][1]
            m = (o.z .< 2000.0) .& (o.z .> 500.0)
            nat = natpde10(k); jl = jlpde10(o)
            for q in 1:4
                rec10a("pde$k:" * CMP10[q] * ":feedback_range:rel_to_max", maximum(abs.(jl[q][m] .- nat[q][m])) / maximum(abs.(nat[q][m])))
                rec10a("pde$k:" * CMP10[q] * ":feedback_range:scaled", maximum(abs.(jl[q][m] .- nat[q][m]) ./ o.scale[q][m]))
                rec10a("pde$k:" * CMP10[q] * ":z_ge_2000:rel_to_global_max", maximum(abs.(jl[q][.!m] .- nat[q][.!m])) / maximum(abs.(nat[q])))
                @test OBS10A["pde$k:" * CMP10[q] * ":feedback_range:rel_to_max"] <= 1.0e-3
            end
            rec10a("pde$k:y_final", maximum(abs.(o.y .- V10["PDE$(k)_yfinal"][1])) / maximum(abs.(V10["PDE$(k)_yfinal"][1])))
        end
        Xe, Te = RUN10.Xe, RUN10.Te
        rec10a("final:Xe", maximum(abs.(Xe .- Xn) ./ Xn)); rec10a("final:Te", maximum(abs.(Te .- Tn) ./ Tn))
        rec10a("final:rows_Xe", maximum(abs.(RUN10.final.rows.Xe .- M10["NAT_OUTPUT"][:, 2]) ./ M10["NAT_OUTPUT"][:, 2]))
        Xr1 = GRID5[:, 2]
        rec10a("diffusion_effect_native_max", maximum(abs.((Xn .- Xr1) ./ Xr1)))
        rec10a("diffusion_effect_error_rel", maximum(abs.(Xe .- Xn) ./ Xr1) / maximum(abs.((Xn .- Xr1) ./ Xr1)))
        @test length(RUN10.final.rows.z) == size(M10["NAT_OUTPUT"], 1)
        @test OBS10A["final:Te"] < 1e-7
        # Marco approved a 1e-5 maximum pointwise relative Xe error on 2026-10-02.
        # The population and other production gates remain unchanged and blocking.
        @test OBS10A["final:Xe"] <= 1e-5
        @test OBS10A["pass1:units:X1s"] < 1.0
        @test OBS10A["pass1:units:2s"] < 10.0 && OBS10A["pass1:units:2p"] < 10.0
    end
    for k in sort(collect(keys(OBS10A)))
        println("Chunk10a ", rpad(k, 30), " ", OBS10A[k])
    end
end
