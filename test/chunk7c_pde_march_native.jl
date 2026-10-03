# Chunk 7c: the HI radiation-PDE march (Lagrange O2 stencils, Step_PDE_O2t, polint_JC lower boundary, production step control) vs the ORIGINAL run.
# REQUIRES COSMOREC_NATIVE_DATA_DIR (fails loudly if unset). Scope: docs/ACCEPTANCE_CHUNK7C.md.
using Test
using CosmoRec
using SHA
@isdefined(SETUP7) || include("chunk7_helpers.jl")
@isdefined(MODEL7N) || include("chunk7b_pde_define_native.jl")

const OBS7C = Dict{String,Float64}()
rec7c(k, v) = (OBS7C[k] = max(get(OBS7C, k, 0.0), v))
srel7c(a, b) = maximum(abs.(a .- b)) / max(maximum(abs.(b)), floatmin())     # spectrum error relative to its largest native magnitude
const FX7C = read_fixture7(FIX7C)
# Dnem cancellation: near zs the pd/Dnem spline values are small differences (6a); Julia and native round them differently (raw relative difference up to
# ~1e-8, residual-scale difference 2e-16). The PDE source is linear in Dnem, so step differences along d y / d Dnem_m are this inherited rounding. They are
# separated by projecting the step difference on the five response directions (computed with a Dnem-scaled coefficient wrapper).
struct DnemScaled7c{C}; c::C; f::Vector{Float64}; end
CosmoRec.hi_Dnem(s::DnemScaled7c, z, m::Integer) = s.f[m] * hi_Dnem(s.c, z, m)
CosmoRec.hi_pd(s::DnemScaled7c, z, m::Integer) = hi_pd(s.c, z, m)
function step7c(f, yin, r)
    m = HIPDEModel(SETUP7, POPS7N, DnemScaled7c(COEF7N, f), ACC5); S = HIPDEStepper{Float64}(SETUP7.x)
    y = copy(yin); hi_pde_step!(y, S, m, r[4], r[2], r[3], r[8], 0.0); return y
end
# FIXED acceptance bound for inherited Dnem rounding (independent of any Julia/native comparison): Dnem = (population term) - exp(-x) loses
# kappa = scale/|Dnem| digits; on the fixed native input history kappa <= 2e8 (asserted below; measured 1.88e8 at z = 2499.5), so each
# implementation's Dnem carries a raw relative rounding of <= 4 eps kappa ~ 1.8e-7. The PDE source is linear in Dnem: raw spectrum differences
# up to this size are accepted, nothing larger.
const KAPPA_MAX7C = 2.0e8
const RAW_BOUND7C = 4 * eps() * KAPPA_MAX7C
fx7c_vec(tag) = (open(FIX7C) do io; for l in eachline(io); startswith(l, tag * "\t") && return vec7([parse(Float64, v) for v in split(l, '\t')[2:end]]); end; error("no $tag"); end)

@testset "Chunk 7c: HI PDE march vs native" begin
    x = SETUP7.x; np = length(x)
    @testset "fixed-bound premise: Dnem cancellation factor of the fixed native input history" begin
        κ = 0.0
        for z in range(500.5, 2499.5; length = 4000), m in 1:5
            κ = max(κ, dnem_scale6(POPS7N, ACC5, z, m, hi_pd(COEF7N, z, m)) / abs(hi_Dnem(COEF7N, z, m)))
        end
        rec7c("kappa_Dnem", κ)
        @test κ <= KAPPA_MAX7C
    end
    @testset "fixture provenance" begin
        @test bytes2hex(open(sha256, FIX7C)) == FIX7C_SHA256
        rep = read_kv7(FIX7C, "REPLICA")
        @test rep["steps"] == "200" && rep["bitwise_equal_to_original"] == "1"
        @test length(FX7C["STEP"]) == 200 && length(fx7c_vec("NATFINAL")) == np
    end
    @testset "Lagrange O2 stencils are exact for quartics" begin
        LG = lagrange_o2(x)
        for p in 0:4, i in 2:(np - 1)
            o = i == 2 ? 0 : (i == np - 1 ? np - 5 : i - 3)
            c = i == 2 ? 2 : (i == np - 1 ? 4 : 3)
            xs = x[(o + 1):(o + 5)]; xi = x[o + c]
            h = maximum(abs.(xs .- xi))
            f = ((xs .- xi) ./ h) .^ p                          # scaled monomial: exact f' = p == 1 ? 1/h : 0, f'' = p == 2 ? 2/h^2 : 0 at xi
            d1 = sum(LG.dli[i, j] * f[j] for j in 1:5) * h; d2 = sum(LG.d2li[i, j] * f[j] for j in 1:5) * h^2
            rec7c("lagrange_d1", abs(d1 - (p == 1 ? 1.0 : 0.0))); rec7c("lagrange_d2", abs(d2 - (p == 2 ? 2.0 : 0.0)))
        end
        @test OBS7C["lagrange_d1"] <= 1e-9 && OBS7C["lagrange_d2"] <= 1e-9
    end
    @testset "lower boundary (polint_JC) and single steps from native spectra" begin
        S = HIPDEStepper{Float64}(x)
        for r in FX7C["STEP"]
            k = Int(r[1]); zin, zout, th, xeval, ylow = r[2], r[3], r[4], r[7], r[8]
            haskey(FX7C, "YIN_$k") || continue
            yin = vec7(FX7C["YIN_$k"][1]); yout = vec7(FX7C["YOUT_$k"][1])
            @test x[1] * (1.0 + zin) / (1.0 + zout) == xeval
            yl, _ = polint_jc(x, yin, xeval, 6)
            rec7c("polint_ylow", abs(yl - ylow) / max(maximum(abs.(yin)), floatmin()))
            S.initialized = false
            y = copy(yin)
            hi_pde_step!(y, S, MODEL7N, th, zin, zout, ylow, 0.0)
            rec7c("single_step_raw", srel7c(y, yout))
            @test y[1] == ylow && y[end] == 0.0
            ε = 1e-4
            Jm = reduce(hcat, [(step7c([i == q ? 1 + ε : 1.0 for i in 1:5], yin, r) .- step7c([i == q ? 1 - ε : 1.0 for i in 1:5], yin, r)) ./ (2ε) for q in 1:5])
            d = yout .- y; c = Jm \ d
            # DIAGNOSTIC ONLY (attribution, not an acceptance gate): residual after projecting the step difference on the Dnem directions
            rec7c("diag:single_step_after_Dnem_projection", maximum(abs.(d .- Jm * c)) / maximum(abs.(yout)))
            k <= 12 && rec7c("diag:Dnem_projection_coefficient_early", maximum(abs.(c)))
        end
        @test OBS7C["polint_ylow"] <= 1e-13
        @test OBS7C["single_step_raw"] <= RAW_BOUND7C
    end
    @testset "full march from y = 0" begin
        steps = FX7C["STEP"]
        outs = Dict{Int,Vector{Float64}}()
        ths = Float64[]
        y, zl = hi_pde_march(MODEL7N; on_step = (k, zin, zout, y, S) -> (push!(ths, zin > 2500.0 - 100.0 ? 0.999 : 0.55); haskey(FX7C, "YOUT_$k") && (outs[k] = copy(y))))
        @test zl == [r[3] for r in steps]
        @test ths == [r[4] for r in steps]
        for (k, v) in outs
            rec7c("march_spectrum", srel7c(v, vec7(FX7C["YOUT_$k"][1])))
        end
        rec7c("march_final", srel7c(y, fx7c_vec("NATFINAL")))
        @test OBS7C["march_spectrum"] <= RAW_BOUND7C        # fixed bound: inherited Dnem rounding (see RAW_BOUND7C)
        @test OBS7C["march_final"] <= 1e-13
    end
    for k in sort(collect(keys(OBS7C)))
        println("Chunk7c ", rpad(k, 22), " ", OBS7C[k])
    end
end
