# Chunk 3d native checks: the assembled pointwise default H/He RHS (fcn_effective) against the ORIGINAL CosmoRec exported fcn_effective,
# captured on 12 explicit physical states (production switches: diffusion/feedback/exotic off, HI absorption on with approximation 1).
using Test
using CosmoRec
@isdefined(FX3D) || include("chunk3d_helpers.jl")

const RTOL3D = 1.0e-13          # observed max over all 15 components x 12 states x {on, off}: < 1e-14
const OBS3D = Dict{String,Any}()

@testset "Chunk 3d: assembled fcn_effective vs original native" begin
    @testset "fixture integrity and native runtime configuration" begin
        @test fixture_sha3d() == FIX3D_SHA256
        c = FX3D.config
        @test c["nShells"] == "3"
        @test length(FX3D.states) == 12
        for key in ("HI_absorption_after_remap", "_HI_abs_appr_flag", "spin_forbidden", "HeI_Feedback", "HeISTfeedback", "Diffusion_correction_is_on", "Diffusion_correction_HeI_is_on", "DM_annihilation", "DM_decay", "magnetic_fields", "flag_He", "f_t")
            @test haskey(c, key)
        end
        @test c["HI_absorption_after_remap"] == "1" && c["_HI_abs_appr_flag"] == "1" && c["spin_forbidden"] == "1" && c["flag_He"] == "1" && c["f_t"] == "1"
        @test c["DM_annihilation"] == "0" && c["DM_decay"] == "0" && c["magnetic_fields"] == "0"
        @test c["Diffusion_correction_is_on"] == "0" && c["Diffusion_correction_HeI_is_on"] == "0" && c["HeI_Feedback"] == "0" && c["HeISTfeedback"] == "0"
        for s in FX3D.states
            @test all(isfinite, s.X) && all(isfinite, s.g["on"]) && all(isfinite, s.g["off"])
            @test s.g["on_cached"] == s.g["on"]            # native static (z, rho) cache must not change the result
        end
    end

    @testset "every component, every state, absorber on and off" begin
        worst = (0.0, "", "", 0)
        for s in FX3D.states, mode in ("on", "off")
            ref = s.g[mode]
            if mode == "on" && !s.dp_in
                @test_throws DPTableDomainError fcn_effective(s.z, s.X, background3d(s), MODEL3D_ON)   # native used the unported explicit DPesc integral
                continue
            end
            g = fcn_effective(s.z, s.X, background3d(s), mode == "on" ? MODEL3D_ON : MODEL3D_OFF)
            for i in 1:15
                @test isapprox(g[i], ref[i]; rtol = RTOL3D, atol = 0.0)
                e = ref[i] == 0 ? abs(g[i]) : abs(g[i] - ref[i]) / abs(ref[i])
                e > worst[1] && (worst = (e, s.label, mode, i - 1))
            end
        end
        OBS3D["worst"] = worst
        @info "Chunk 3d native max relative error $(worst[1]) at state $(worst[2]), mode $(worst[3]), component g[$(worst[4])]"
        @test worst[1] < RTOL3D
    end

    @testset "ground-state electron fraction follows native compute_fractions" begin
        for s in FX3D.states
            xe, xp, xheii = effective_fractions(s.X, s.fHe)
            @test xp ≈ s.Xp rtol = 1e-14
            @test xheii ≈ s.XHeII rtol = 1e-13
            @test xe ≈ s.Xe rtol = 1e-14
            @test xe == xp + xheii
            # X[0] (Xe slot) is ignored by native fcn_effective: scrambling it cannot change the result
            X2 = copy(s.X); X2[1] = 7.5
            @test fcn_effective(s.z, X2, background3d(s), MODEL3D_OFF) == fcn_effective(s.z, s.X, background3d(s), MODEL3D_OFF)
        end
    end

    @testset "absorber: threshold z <= 3400, slots and transfer invariants" begin
        s_above, s_at, s_below = FX3D.states[1], FX3D.states[2], FX3D.states[3]
        @test (s_above.z, s_at.z, s_below.z) == (3400.5, 3400.0, 3399.0)
        for s in FX3D.states
            d_nat = s.g["on"] .- s.g["off"]
            if s.z > 3400
                @test all(iszero, d_nat)
            else
                @test any(!iszero, d_nat[ABS_SLOTS3D])
            end
            others = setdiff(1:15, ABS_SLOTS3D)
            @test all(iszero, d_nat[others])                       # native: absorber touches only slots 0, 1, 7, 9, 12
            s.dp_in || continue
            bg = background3d(s)
            gon = fcn_effective(s.z, s.X, bg, MODEL3D_ON); goff = fcn_effective(s.z, s.X, bg, MODEL3D_OFF)
            d = zeros(15)   # absorber contribution evaluated directly (no cancellation against the large H/He terms)
            hi_absorption_rhs!(d, 1, 2, 8, FX3D.dp, FX3D.bitot, FX3D.fc, s.z, s.Tg, s.NH, s.Hz, s.X[2], s.X[8:14])
            @test isapprox(gon, goff .+ d; rtol = 1e-13)
            @test all(iszero, d[others])
            # reaction transfer: electron created = H 1s destroyed; HeI population conserved by the absorber
            @test d[1] ≈ -d[2] rtol = 1e-12 atol = 1e-12 * abs(d[1])
            @test abs(d[8] + d[10] + d[13]) <= 1e-12 * (abs(d[8]) + abs(d[10]) + abs(d[13]))
            @test d[1] ≈ d[8] rtol = 1e-12 atol = 1e-12 * abs(d[1])
            @test d[1] ≈ -(d[10] + d[13]) rtol = 1e-12 atol = 1e-12 * abs(d[1])
            for i in ABS_SLOTS3D
                @test abs(d[i] - d_nat[i]) <= 3e-13 * max(abs(gon[i]), abs(goff[i])) + 1e-9 * abs(d_nat[i])
            end
        end
        # exact-switch behaviour: the absorber is a branch in z (one-sided, NOT differentiable across zcrit); z is never differentiated
        X = s_at.X; bg = background3d(s_at)
        for (z, active) in ((3400.0, true), (prevfloat(3400.0), true), (nextfloat(3400.0), false), (3400.5, false), (3399.0, true))
            d = fcn_effective(z, X, bg, MODEL3D_ON) .- fcn_effective(z, X, bg, MODEL3D_OFF)
            @test (d[1] != 0) == active
        end
    end

    @testset "transfer invariants of the native-form components" begin
        for s in FX3D.states
            Tg = s.Tg; rho = s.X[15]
            XH = s.X[2:7]; XHe = s.X[8:14]
            rh = get_rates(FX3D.htable, Tg, rho * Tg)
            sums = Float64[]
            for fn in (d -> two_photon!(d, Tg, XH), d -> lyman!(d, Tg, XH, s.NH, s.Hz), d -> interlevel!(d, Tg, XH, rh.R))
                d = zeros(6); fn(d)
                push!(sums, abs(sum(d)) / max(sum(abs, d), floatmin()))
            end
            re = get_helium_rates(FX3D.hetable, Tg)
            for fn in (d -> helium_two_photon!(d, Tg, XHe), d -> helium_lyman!(d, Tg, XHe, s.NH, s.Hz), d -> helium_interlevel!(d, Tg, XHe, re.R))
                d = zeros(7); fn(d)
                push!(sums, abs(sum(d)) / max(sum(abs, d), floatmin()))
            end
            @test maximum(sums) < 1e-12
            # energy/temperature equation: at rho = 1 there is no Compton coupling, d rho/dt = -H(z) exactly
            if s.X[15] == 1
                g = fcn_effective(s.z, s.X, background3d(s), MODEL3D_OFF)
                @test g[15] ≈ -s.Hz rtol = 1e-13
            end
        end
    end
end
