# Chunk 7b: coefficients A, B, C, D of the HI radiation PDE (def_PDE_Lyn_and_2s1s) and its side products vs the ORIGINAL routine after the production setup.
# REQUIRES COSMOREC_NATIVE_DATA_DIR (fails loudly if unset). Scope: docs/ACCEPTANCE_CHUNK7B.md.
using Test
using CosmoRec
using SHA
@isdefined(SETUP7) || include("chunk7_helpers.jl")

const OBS7B = Dict{String,Float64}()
rec7b(k, v) = (OBS7B[k] = max(get(OBS7B, k, 0.0), v))
rel7b(a, b) = a == b ? 0.0 : abs(a - b) / max(abs(b), floatmin())
vrel7b(a, b) = maximum(abs.(a .- b)) / max(maximum(abs.(b)), floatmin())
prel7b(a, b) = (nz = findall(!iszero, b); isempty(nz) ? 0.0 : maximum(rel7b.(a[nz], b[nz])))

const FX7B = read_fixture7(FIX7B)
const ROWS7N = NODES5[:, 1:9]
const POPS7N = HIPopulationSplines(ROWS7N; zs = ZS6, ze = ZE6)
const COEF7N = hi_pde_coefficients(ROWS7N, POPS7N, ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6)
const MODEL7N = HIPDEModel(SETUP7, POPS7N, COEF7N, ACC5)

@testset "Chunk 7b: HI PDE coefficients vs native" begin
    @testset "fixture provenance and atom/constants" begin
        @test bytes2hex(open(sha256, FIX7B)) == FIX7B_SHA256
        con = read_kv7(FIX7A, "CON7"); k = NATIVE_HI_PDE_CONSTANTS
        @test (parse(Float64, con["h_kb"]), parse(Float64, con["kB"]), parse(Float64, con["me_gr"]), parse(Float64, con["cl"]), parse(Float64, con["mH_gr"]), parse(Float64, con["sigT"])) ==
              (k.h_kb, k.kB, k.me_gr, k.cl, k.mH_gr, k.sigT)
        lyn = Dict(Int(r[1]) => r for r in FX7A["LYN"]); a = NATIVE_HI_PDE_ATOM
        for n in 2:3
            @test (lyn[n][2], lyn[n][4], lyn[n][5], lyn[n][6]) == (a.lyn[n - 1].nu21, a.lyn_A21[n - 1], a.lyn[n - 1].Gamma, a.lyn[n - 1].AM)
        end
        lev = Dict((Int(r[1]), Int(r[2])) => r for r in FX7A["LEV"])
        @test lev[(2, 1)][3] == a.nu21 && lev[(2, 0)][3] == a.Dnu_1s[2] && lev[(3, 0)][3] == lev[(3, 1)][3] == lev[(3, 2)][3] == a.Dnu_1s[3]
        @test (lev[(3, 1)][7], lev[(3, 1)][8]) == (a.A_3p2s, a.nu_3p2s)
        @test (lev[(3, 0)][5], lev[(3, 0)][6]) == (a.A_3s2p, a.nu_3s2p)
        @test (lev[(3, 2)][5], lev[(3, 2)][6]) == (a.A_3d2p, a.nu_3d2p)
        @test lev[(2, 1)][5] == a.lyn_A21[1] && lev[(3, 1)][5] == a.lyn_A21[2]
        @test length(FX7B["DEFZ"]) == 9 && all(length(vec7(r)) == 2353 for r in FX7B["DEF_A"])
    end

    @testset "voigt_dphi_dx" begin
        for x in (-35.0, -12.3, -2.0, -0.3, 0.0, 0.7, 3.1, 29.9, 30.0, 80.0), av in (1e-4, 3e-3)
            h = 1e-5 * max(1.0, abs(x))
            fd = (voigt_phi(x + h, av) - voigt_phi(x - h, av)) / (2h)
            abs(x) < 29.5 && rec7b("dphi_vs_fd", abs(voigt_dphi_dx(x, av) - fd) / max(abs(fd), 1e-12))
        end
        @test OBS7B["dphi_vs_fd"] <= 1e-6
    end

    st = HIPDEState{Float64}(length(SETUP7.x))
    @testset "def_PDE_Lyn_and_2s1s at 9 redshifts" begin
        for (j, r) in enumerate(FX7B["DEFZ"])
            z = r[1]
            rec7b("bg:H", rel7b(cosmos_H(ACC5, z), r[2])); rec7b("bg:NH", rel7b(cosmos_NH(ACC5, z), r[3])); rec7b("bg:TCMB", rel7b(cosmos_TCMB(ACC5, z), r[4]))
            rec7b("bg:Xe", rel7b(hi_Xe(POPS7N, z), r[5])); rec7b("bg:X1s", rel7b(hi_Xi(POPS7N, z, 0), r[6])); rec7b("bg:rho", rel7b(hi_rho(POPS7N, z), r[7]))
            hi_pde_rhs_coefficients!(st, MODEL7N, z)
            for (tag, v) in (("aV", st.aV), ("Dnu", st.Dnu), ("pd", st.pd), ("pdeff", st.pd_eff), ("Dnemeff", st.Dnem_eff))
                ref = vec7(FX7B["DEF_$tag"][j])
                @test ref[3] == 0.0
                rec7b(tag, maximum(rel7b.(v, ref[1:2])))
            end
            for (tag, v) in (("A", st.A), ("B", st.B), ("C", st.C), ("D", st.D))
                ref = vec7(FX7B["DEF_$tag"][j])
                rec7b(tag, vrel7b(v, ref)); rec7b(tag * "_pointwise", prel7b(v, ref))
                @test all(iszero, v[findall(iszero, ref)])
            end
        end
        hi_pde_rhs_coefficients!(st, MODEL7N, 1300.0)
        rec7b("z1300:exp_x", prel7b(st.exp_x, vec7(FX7B["DEF_expx"][1])))
        for n in 1:2
            rec7b("z1300:phi", prel7b(st.phi[n], vec7(FX7B["DEF_phi"][n]))); rec7b("z1300:hp_phi", prel7b(st.hp_phi[n], vec7(FX7B["DEF_HPphi"][n])))
            rec7b("z1300:P2g_ns", vrel7b(st.P2g_ns[n], vec7(FX7B["DEF_P2gns"][n]))); rec7b("z1300:P2g_nd", vrel7b(st.P2g_nd[n], vec7(FX7B["DEF_P2gnd"][n])))
            rec7b("z1300:PR_ns", vrel7b(st.PR_ns[n], vec7(FX7B["DEF_PRns"][n]))); rec7b("z1300:PR_nd", vrel7b(st.PR_nd[n], vec7(FX7B["DEF_PRnd"][n])))
        end
        @test all(iszero, vec7(FX7B["DEF_phi"][3])) && all(iszero, vec7(FX7B["DEF_PRnd"][1])) && all(iszero, vec7(FX7B["DEF_P2gnd"][1]))
        for k in ("bg:H", "bg:NH", "bg:TCMB", "bg:Xe", "bg:X1s", "bg:rho", "aV", "Dnu", "pd", "pdeff", "Dnemeff", "A", "B", "C", "D",
                  "z1300:exp_x", "z1300:phi", "z1300:hp_phi", "z1300:P2g_ns", "z1300:P2g_nd", "z1300:PR_ns", "z1300:PR_nd")
            @test OBS7B[k] <= 1e-12
        end
    end
    for k in sort(collect(keys(OBS7B)))
        println("Chunk7b ", rpad(k, 22), " ", OBS7B[k])
    end
end
