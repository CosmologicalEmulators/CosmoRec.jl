# Chunk 9a: feedback of the HI PDE corrections into the ODE right-hand side (setup_DF_interpol_data / interpolate_DF and the corrections added at the end
# of fcn_HI_effective) vs the ORIGINAL routines. REQUIRES COSMOREC_NATIVE_DATA_DIR (fails loudly if unset). Scope: docs/ACCEPTANCE_CHUNK9A.md.
using Test
using CosmoRec
using SHA
@isdefined(SETUP7) || include("chunk7_helpers.jl")
@isdefined(MODEL7N) || include("chunk7b_pde_define_native.jl")

const FIX9A = joinpath(@__DIR__, "fixtures", "native_hi_diffusion_feedback.txt")
const FIX9A_SHA256 = "462703809df1a502b0d60cbaa9de7d1e123e5addddd73dfabaacfb50f43ab78c"
const OBS9A = Dict{String,Float64}()
rec9a(k, v) = (OBS9A[k] = max(get(OBS9A, k, 0.0), v))
const FX9A = read_fixture7(FIX9A)
const FX8A9 = read_fixture7(FIX8A)
const FB9N = hi_diffusion_feedback(vec7(FX8A9["DF_z"][1]), vec7(FX8A9["DI1_2s"][1]), [vec7(FX8A9["DF_2g"][1]), vec7(FX8A9["DF_2g"][2])], [vec7(FX8A9["DF_R"][1])])
const RM9N = with_diffusion(RM5, FB9N)
# F9 record: id, on, flag, z, pattern, Y..., F... (the "Y"/"F" labels are dropped by read_fixture7)
function f9rows()
    rows = []
    for r in FX9A["F9"]
        flag = r[3] == 1.0; n = ode_nstate(flag)
        push!(rows, (id = Int(r[1]), on = r[2] == 1.0, flag = flag, z = r[4], y = r[6:(5 + n)], f = r[(6 + n):(5 + 2n)]))
    end
    return rows
end

@testset "Chunk 9a: HI diffusion feedback vs native" begin
    @testset "fixture provenance and configuration" begin
        @test bytes2hex(open(sha256, FIX9A)) == FIX9A_SHA256
        cfg = FX9A["CFG9"][1]
        @test cfg == [500.0, 2000.0, 1.0]
        @test FB9N.zmin == 500.0 && FB9N.zmax == 2000.0 && FB9N.DI1 !== nothing
        @test length(f9rows()) == 180
        @test hi_diffusion_rhs!(zeros(15), nothing, 1000.0, ones(15)) == zeros(15)
        @test RecombinationModel(EFF5, ACC5).diffusion === nothing
    end
    @testset "interpolation of the native PDE outputs" begin
        for r in FX9A["DFI"]
            z = r[1]
            jl = (hi_DI1_2s(FB9N, z), hi_DF_2g(FB9N, 0, z), hi_DF_2g(FB9N, 1, z), hi_DF_R(FB9N, 0, z))
            for q in 1:4
                rec9a("DFI", r[q + 1] == 0.0 ? abs(jl[q]) : abs(jl[q] - r[q + 1]) / abs(r[q + 1]))
                r[q + 1] == 0.0 && @test jl[q] == 0.0
            end
        end
        @test OBS9A["DFI"] <= 1e-12
    end
    rows = f9rows()
    @testset "ODE right-hand side with the correction on (native DF outputs)" begin
        for r in rows
            f = recombination_rhs(r.z, r.y, r.on ? RM9N : RM5; flag_He = r.flag)
            sc = maximum(abs, r.f)
            rec9a(r.on ? "F_on" : "F_off", maximum(abs.(f .- r.f) ./ max.(abs.(r.f), 1e-12 * sc)))
        end
        @test OBS9A["F_on"] < 1e-12 && OBS9A["F_off"] < 1e-12
    end
    @testset "the correction term itself (native F_on - F_off)" begin
        on = [r for r in rows if r.on]
        noff = length(rows) ÷ 2
        for r in on
            ro = rows[r.id + noff]
            @test ro.z == r.z && ro.y == r.y && !ro.on
            dnat = r.f .- ro.f
            fj = recombination_rhs(r.z, r.y, RM9N; flag_He = r.flag); fj0 = recombination_rhs(r.z, r.y, RM5; flag_He = r.flag)
            # Julia correction term, exactly: the added g terms divided by dz/dt
            bg = ode_background(RM5, r.z); X = ode_unpack(r.y, bg.fHe, r.flag)
            g = hi_diffusion_rhs!(zeros(15), FB9N, r.z, X); dz_dt = -bg.Hz * (1.0 + r.z)
            djl = zeros(length(r.y)); for k in 1:6; djl[1 + k] = g[1 + k] / dz_dt; end
            # the native difference of two rounded F values resolves the correction only to ulp(F_on) + ulp(F_off) (|F| exceeds the correction by up
            # to 4e12 for the explicit excited states): fixed rounding-unit bound
            units = abs.(djl .- dnat) ./ (eps.(abs.(r.f)) .+ eps.(abs.(ro.f)) .+ 1e-12 .* abs.(dnat))
            rec9a("dF_term_rounding_units", maximum(units))
            rec9a("dF_term_rel_nonzero", maximum([abs(djl[k] - dnat[k]) / abs(dnat[k]) for k in eachindex(dnat) if abs(dnat[k]) > 1e3 * eps() * abs(r.f[k])]; init = 0.0))
            (r.z >= 2000.0 || r.z <= 500.0) && @test all(iszero, djl) && fj == fj0
        end
        @test OBS9A["dF_term_rounding_units"] <= 2.0
    end
    @testset "feedback built from the Julia PDE stage" begin
        out = hi_pde_corrections(MODEL7N)
        fb = hi_diffusion_feedback(out)
        @test fb.zmin == FB9N.zmin && fb.zmax == FB9N.zmax
        for r in FX9A["DFI"]
            z = r[1]
            jl = (hi_DI1_2s(fb, z), hi_DF_2g(fb, 0, z), hi_DF_2g(fb, 1, z), hi_DF_R(fb, 0, z))
            for q in 1:4
                rec9a("DFI_julia_stage:scaled", abs(jl[q] - r[q + 1]) / maximum(out.scale[q]))
            end
        end
        # fixed bound of Chunk 8a (native quadrature tolerance, 1e-4 of the residual scale); interpolation is linear in the node values
        @test OBS9A["DFI_julia_stage:scaled"] <= 1.0e-4
    end
    for k in sort(collect(keys(OBS9A)))
        println("Chunk9a ", rpad(k, 30), " ", OBS9A[k])
    end
end
