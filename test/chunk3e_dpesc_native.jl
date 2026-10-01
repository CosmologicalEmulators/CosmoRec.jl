# Chunk 3e: native DPesc_coh (explicit coherent-scattering fallback of the H-I absorber) vs the plain-text fixture produced by the ORIGINAL library
# (call_DP_Singlet/Triplet, DPesc_appr_I_sym, DP_interpol_S/T, Voigt and H 1s cross-section objects). Tolerances are justified in docs/CHUNK3E_RESULTS.md.
using Test
using CosmoRec
using SHA
@isdefined(FX3E) || include("chunk3e_helpers.jl")
@isdefined(FX3D) || include("chunk3d_helpers.jl")

@isdefined(RTOL3D) || (const RTOL3D = 1.0e-13)
const OBS3E = Dict{String,Float64}()
rec3e(k, v) = (OBS3E[k] = max(get(OBS3E, k, 0.0), v))
relerr3e(a, b) = abs(a - b) / abs(b)
fallback3e(tr) = dpesc_fallback(NATIVE_DPESC, tr)

function patterson_nodes(lv)
    pos = [0.0]
    for l in 2:lv
        new = PATTERSON_X[l]
        pos = vcat([[pos[i], new[i]] for i in eachindex(new)]...)
        length(pos) == 2 * length(new) || error("bad interleave")
    end
    return pos
end

@testset "Chunk 3e: native DPesc_coh fixture parity" begin
    @testset "fixture provenance" begin
        @test bytes2hex(open(sha256, FIX3E)) == FIX3E_SHA256
        @test length(FX3E.queries) == 198
        @test count(q -> q.triplet, FX3E.queries) == 99
        fbq = filter(q -> q.corr == q.tab, FX3E.queries)       # native DP_interpol_S/T returned call_DP_* bit-for-bit
        @test length(fbq) == 64
        @test Set(q.label for q in fbq) == Set(["z1500", "z1200", "z800"])
        @test all(q -> q.label in ("z1500", "z1200", "z800") || q.corr != q.tab, FX3E.queries)
    end

    @testset "Patterson rules generated from first principles" begin
        @test PATTERSON_NPOINTS == (1, 3, 7, 15, 31, 63, 127, 255)
        @test PATTERSON_X[2] ≈ [0.7745966692414834] rtol = 1e-15
        @test PATTERSON_X[3] ≈ [0.4342437493468026, 0.9604912687080203] rtol = 1e-15
        for lv in 2:8
            pos = patterson_nodes(lv)
            w = PATTERSON_W[lv]
            @test length(pos) == (PATTERSON_NPOINTS[lv] + 1) ÷ 2 == length(w)
            @test issorted(pos) && all(0 .< pos[2:end] .< 1)
            @test w[1] + 2 * sum(w[2:end]) ≈ 2 rtol = 1e-15
            deg = (3 * PATTERSON_NPOINTS[lv] + 1) ÷ 2         # Patterson rule exactness degree
            for k in 0:2:min(deg, 40)
                val = w[1] * (k == 0 ? 1.0 : 0.0) + 2 * sum(w[i] * pos[i]^k for i in 2:length(pos))
                @test val ≈ 2 / (k + 1) atol = 2e-15
            end
        end
    end

    @testset "Voigt profile, line data and H 1s cross-section vs the native objects" begin
        for t in FX3E.vt
            line = t.ch == "S" ? NATIVE_DPESC.singlet : NATIVE_DPESC.triplet
            @test voigt_DnuT(line, t.Tm) ≈ t.DnuT rtol = 1e-15
            @test voigt_a(line, t.Tm) ≈ t.a rtol = 1e-15
            @test voigt_xi_int(1e4, voigt_a(line, t.Tm)) ≈ t.xi1e4 rtol = 1e-15
        end
        for t in FX3E.vx
            line = t.ch == "S" ? NATIVE_DPESC.singlet : NATIVE_DPESC.triplet
            a = voigt_a(line, t.Tm)
            @test voigt_phi(t.x, a) ≈ t.phi rtol = 1e-14
            @test isapprox(voigt_xi_int(t.x, a), t.xi; rtol = 1e-13, atol = 1e-16)
            rec3e("phi", relerr3e(voigt_phi(t.x, a), t.phi)); rec3e("xi", abs(voigt_xi_int(t.x, a) - t.xi) / max(abs(t.xi), 1e-16))
        end
        sig_nuc = 6.3111866125271474e-18
        for (r, s) in FX3E.sx
            @test lyc_ratio(NATIVE_HLYC, r * NATIVE_HLYC.nu_ion) * sig_nuc ≈ s rtol = 1e-14 atol = 0.0
        end
    end

    @testset "intermediate integral, Pesc and the correction, every captured query" begin
        for q in FX3E.queries
            line = q.triplet ? NATIVE_DPESC.triplet : NATIVE_DPESC.singlet
            Pd = sobolev_p(q.pd * q.tauS)
            @test relerr3e(Pd, q.Pd) < 1e-14
            @test relerr3e(sobolev_p(q.tauS), q.PS) < 1e-14
            Dp, ord = dpesc_appr_I_sym(q.T, q.pd * q.tauS, q.eta, line, NATIVE_DPESC.hlyc, Pd * 1e-5)
            c, _ = dpesc_coh(NATIVE_DPESC, q.triplet, q.tauS, q.eta, q.T, q.T, q.pd)
            P = Pd + Dp
            Pesc = q.pd * P / (1 - (1 - q.pd) * P)
            # singlet: agreement at rounding level; triplet: the integral (~1e-7) is a cancellation of O(1) terms, native itself is ~4e-11 from the 256-bit value
            tolD = q.triplet ? 1e-9 : 1e-12
            tolc = q.triplet ? 5e-9 : 1e-12
            @test relerr3e(Dp, q.Dpij) < tolD
            @test relerr3e(Pesc, q.Pesc) < tolc
            @test relerr3e(c, q.corr) < tolc
            rec3e((q.triplet ? "T" : "S") * ":Dpij", relerr3e(Dp, q.Dpij)); rec3e((q.triplet ? "T" : "S") * ":corr", relerr3e(c, q.corr))
        end
    end

    @testset "final channel rates with the fallback correction" begin
        for q in filter(q -> q.corr == q.tab, FX3E.queries)
            @test relerr3e(dpesc_coh(NATIVE_DPESC, q.triplet, q.tauS, q.eta, q.T, q.T, q.pd)[1], q.tab) < (q.triplet ? 5e-9 : 1e-12)
        end
    end

    @testset "composed fcn_effective at the off-table states z = 1500, 1200, 800 (absorber on and on-minus-off vs native)" begin
        worst = 0.0
        for s in FX3D.states
            s.dp_in && continue
            @test_throws DPTableDomainError fcn_effective(s.z, s.X, background3d(s), MODEL3D_ON)     # no silent fallback
            g = fcn_effective(s.z, s.X, background3d(s), MODEL3D_ON; dp_fallback = fallback3e)
            for i in 1:15
                @test isapprox(g[i], s.g["on"][i]; rtol = RTOL3D, atol = 0.0)
                worst = max(worst, s.g["on"][i] == 0 ? abs(g[i]) : abs(g[i] - s.g["on"][i]) / abs(s.g["on"][i]))
            end
            goff = fcn_effective(s.z, s.X, background3d(s), MODEL3D_OFF)
            dn = s.g["on"] .- s.g["off"]; dj = g .- goff
            for i in 1:15
                @test isapprox(dj[i], dn[i]; rtol = 1e-9, atol = 1e-12 * maximum(abs, dn))
            end
            @test dn[8] != 0 && dn[1] != 0         # the absorber contribution is not zero
        end
        rec3e("fcn_on", worst)
        @test worst < RTOL3D
        # states inside the table do not change with a fallback supplied
        for s in FX3D.states
            s.dp_in || continue
            @test fcn_effective(s.z, s.X, background3d(s), MODEL3D_ON; dp_fallback = fallback3e) == fcn_effective(s.z, s.X, background3d(s), MODEL3D_ON)
        end
    end

    @testset "integration refinement against the native integral" begin
        # forcing a higher Patterson order on every sub-interval changes the primal by a bounded amount only
        q = first(filter(q -> !q.triplet && q.label == "z1500" && q.code == "000", FX3E.queries))
        line = NATIVE_DPESC.singlet
        Pd = sobolev_p(q.pd * q.tauS)
        Dp, ord = dpesc_appr_I_sym(q.T, q.pd * q.tauS, q.eta, line, NATIVE_DPESC.hlyc, Pd * 1e-5)
        @test relerr3e(Dp, q.Dpij) < 1e-12
        for lv in 5:8
            D2, _ = dpesc_appr_I_sym(q.T, q.pd * q.tauS, q.eta, line, NATIVE_DPESC.hlyc, 0.0; orders = ntuple(_ -> lv, 7))
            rec3e("refine$lv", relerr3e(D2, q.Dpij))
        end
        D8, _ = dpesc_appr_I_sym(q.T, q.pd * q.tauS, q.eta, line, NATIVE_DPESC.hlyc, 0.0; orders = ntuple(_ -> 8, 7))
        @test relerr3e(D8, q.Dpij) < 1e-5            # the native integral is the rel 1e-6-stopped one; the 255-point everywhere value must be within the stop criterion
    end

    @testset "info" begin
        foreach(kv -> println("Chunk3e-native ", kv[1], " = ", kv[2]), sort!(collect(OBS3E)))
    end
end
