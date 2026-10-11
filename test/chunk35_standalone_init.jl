# Chunk 35: standalone initialization (A1 Cosmos closure, A2 in-pipeline preliminary Recfast++ history, A3/A8 checked ABI-shaped entry) vs the ORIGINAL.
# Native references: harness4b_param / harness10a_param readbacks (text fixtures native_standalone_*.txt; the fiducial uses the existing 4b/10a fixtures).
# Callbacks: preliminary history = Rodas5P dense, reltol 1e-12, abstol (1e-16, 1e-16, 1e-12) (the chunk4c validation callback); main ODE = the production
# default of solve_phase5 (Rodas4P, reltol 1e-10, a1 1e-18, aex 1e-12, tstops = znodes); tail = solve_tail5. REQUIRES COSMOREC_NATIVE_DATA_DIR (no skip).
# Scope and full per-cosmology results (10 points, both main-ODE callbacks, CAMB C_ell): local_analysis/.../chunk35_standalone_init/REPORT.md.
using Test
using CosmoRec
using SHA
@isdefined(SETUP7) || include("chunk7_helpers.jl")

const FIX35 = Dict("closure" => "native_standalone_closure.txt", "hubble" => "native_standalone_h50_hubble.txt",
                   "history" => "native_standalone_h50_history.txt", "output" => "native_standalone_h50_output.txt")
const FIX35_SHA = Dict("closure" => "aaa2392d5737b4485b0a1ef161ea8ad1f18589e6da4402a8f9a7952556ae18be", "hubble" => "2086d4f5f5539f11b5216bd7fbb3160d32552eb32f907f167324307483bc7a3f",
                       "history" => "fb0b5d4277a234f0fc808bf25c65eee3de47f173e6c4492cc01bed9b41a7aca9", "output" => "4e990ce08dd72907914d154e3a1d74e2ac6adaca6352df829db92194f690ff73")
fix35(k) = joinpath(@__DIR__, "fixtures", FIX35[k])
const OBS35 = Dict{String,Float64}()
rec35(k, v) = (OBS35[k] = max(get(OBS35, k, 0.0), v))
ulp35(a, b) = a == b ? 0.0 : abs(a - b) / eps(abs(b))

"""Closure fixture: point name => (inputs, Dict(tag => Dict(key => value)))."""
function read_closure35(path = fix35("closure"))
    d = Dict{String,Any}()
    for l in eachline(path)
        (isempty(l) || startswith(l, "#")) && continue
        f = split(l, '\t')
        if f[1] == "PT"
            d[f[2]] = (inp = parse.(Float64, f[3:8]), rec = Dict{String,Dict{String,Float64}}())
        else
            r = get!(d[f[1]].rec, String(f[2]), Dict{String,Float64}())
            if length(f) == 4
                v = tryparse(Float64, f[4]); v === nothing || (r[String(f[3])] = v)
            else
                for k in 3:2:(length(f) - 1); v = tryparse(Float64, f[k + 1]); v === nothing || (r[String(f[k])] = v); end
            end
        end
    end
    return d
end
read_rows35(path) = [parse.(Float64, split(l)) for l in eachline(path) if !isempty(l) && !startswith(l, "#") && !startswith(l, "SWITCH")]

const CL35 = read_closure35()
const SRF35 = (a...) -> solve_rodas_dense(a...; reltol = 1.0e-12, abstol = [1.0e-16, 1.0e-16, 1.0e-12])
const TABLES35 = load_runmode0_tables(dirname(DATADIR5))
const ZPASS35 = init_xarr_linear(3000.0, 50.0, 3000)
const FIX10A35 = joinpath(@__DIR__, "fixtures", "native_cosmorec_runmode0.txt")      # Chunk 10a fixture (fiducial ORIGINAL runmode-0 call)
"""ZREC, NAT_XE, NAT_TB vectors and the PASSINFO switch redshifts of the Chunk 10a fixture."""
function read_nat10_35(path = FIX10A35)
    v = Dict{String,Vector{Float64}}(); sw = Float64[]
    for l in eachline(path)
        f = split(l, '\t')
        (f[1] in ("ZREC", "NAT_XE", "NAT_TB")) && !haskey(v, f[1]) && (v[String(f[1])] = parse.(Float64, f[3:end]))
        f[1] == "PASSINFO" && push!(sw, parse(Float64, f[7]))
    end
    return v, sw
end
const NAT10_35, NSW10_35 = read_nat10_35()
const DFID35 = camb_density_parameters(0.02237, 0.12, 67.36)

"""Native-tolerance units of the native Recfast++ solver (rtol 1e-8, atol (1e-10, 1e-8, 1e-10) for (Xe_He, Xe_H, TM)) on the ODE nodes."""
units35(v, ref, atol) = maximum(abs.(v .- ref) ./ (1.0e-8 .* abs.(ref) .+ atol))

function check_prelim35(tag, pre, rows)       # rows: (node, z, Xe_H, Xe_He, Xe, dXe, dXe_H, TM)
    h = pre.history; idx = Int.(getindex.(rows, 1)); R = reduce(hcat, rows)'
    @test length(h.z) == 6000 && pre.grid.z[1] == 25000.0 && pre.grid.z[end] == 0.0
    @test maximum(abs.(h.z[idx] .- R[:, 2]) ./ max.(R[:, 2], 1.0)) < 1.0e-12
    j = pre.grid.jode; sa = idx .< j; od = idx .> j
    @test h.Xe_H[idx[sa]] == R[sa, 3] && h.Xe_He[idx[sa]] == R[sa, 4] && h.dXe[idx[sa]] == R[sa, 6]
    rec35("$tag:saha_Xe_rel", maximum(abs.(h.Xe[idx[sa]] .- R[sa, 5]) ./ R[sa, 5])); @test maximum(abs.(h.Xe[idx[sa]] .- R[sa, 5]) ./ R[sa, 5]) < 1.0e-10
    uHe = units35(h.Xe_He[idx[od]], R[od, 4], 1.0e-10); uH = units35(h.Xe_H[idx[od]], R[od, 3], 1.0e-8); uT = units35(h.TM[idx[od]], R[od, 8], 1.0e-10)
    rec35("$tag:units_Xe_He", uHe); rec35("$tag:units_Xe_H", uH); rec35("$tag:units_TM", uT)
    @test uHe < 1.0 && uH < 2.0 && uT < 100.0               # the chunk4c gates (native solver error, observed 0.2 / <= 1.6 / <= 53)
    rx = abs.(h.Xe[idx] .- R[:, 5]) ./ R[:, 5]; rt = abs.(h.TM[idx] .- R[:, 8]) ./ R[:, 8]; z1 = R[:, 2] .>= 1.0
    rec35("$tag:Xe_rel_all", maximum(rx)); rec35("$tag:Xe_rel_z>=1", maximum(rx[z1])); rec35("$tag:TM_rel_all", maximum(rt))
    @test maximum(rx[z1]) <= 1.0e-5                           # working gate met for z >= 1
    # NOT met for z < 1 (observed 3.9e-5 .. 9.2e-5 at the 10 chunk35 points, worst at z = 0): the native Gear solver's absolute tolerance 1e-8 on x_p
    # (~1 native unit) at Xe ~ 2e-4; TM <= 1e-7 is NOT met either (observed ~5e-7 = ~50 native units). Reported in REPORT.md, kept as envelopes here:
    @test maximum(rx) <= 2.0e-4 && maximum(rt) <= 1.0e-6
end

function check_e2e35(tag, r, Z, natXe, natTb, idx, nswitch)
    rx = abs.(r.Xe[idx] .- natXe) ./ natXe; rt = abs.(r.Tm[idx] .- natTb) ./ natTb; m = Z[idx] .> 0.001
    rec35("$tag:Xe_rel_z>0.001", maximum(rx[m])); rec35("$tag:Tm_rel_z>0.001", maximum(rt[m])); rec35("$tag:Xe_rel_all", maximum(rx)); rec35("$tag:Tm_rel_all", maximum(rt))
    @test maximum(rx[m]) <= 1.0e-5                            # approved final-history gate, all nodes served by the CosmoRec rows
    @test maximum(rt[m]) < 1.0e-7
    # the only node at z <= 0.001 (z = 0) is served by the preliminary history (CosmoRec.cpp:584-588) and inherits its native-solver offset (see check_prelim35)
    @test r.Xe[end] == cosmos_Xe_Seager(r.preliminary.accessors, Z[end]) && r.Tm[end] == cosmos_Te(r.preliminary.accessors, Z[end])
    @test maximum(rx) <= 2.0e-4 && maximum(rt) <= 1.0e-6
    @test [ZPASS35[k] for k in r.k_switch] == nswitch         # helium-switch node of every pass equals the native one
end

@testset "Chunk 35: standalone initialization vs native" begin
    @testset "fixture provenance" begin
        for k in keys(FIX35); @test bytes2hex(open(sha256, fix35(k))) == FIX35_SHA[k]; end
        @test length(CL35) == 10
    end

    @testset "A1: Cosmos closure vs native readback (10 cosmologies)" begin
        for (name, p) in CL35
            H0, ombh2, omch2, YHe, nnu, T0 = p.inp
            d = camb_density_parameters(ombh2, omch2, H0)
            cfg = p.rec["CFG"]; c2 = p.rec["COS2"]; cs = p.rec["CONST"]
            @test d.Omega_b == cfg["Omega_b"] && d.Omega_c + d.Omega_b == cfg["Omega_m"] && H0 / 100.0 == cfg["h100"]
            c = cosmos_constants(; Omega_c = d.Omega_c, Omega_b = d.Omega_b, Omega_k = 0.0, N_eff = nnu, h100 = H0 / 100.0, T_CMB = T0, Y_p = YHe)
            θ = recfast_theta(; Omega_b = d.Omega_b, h100 = H0 / 100.0, T_CMB = T0, Y_p = YHe)
            for (k, a, b) in (("Y_p", c.Y_p, c2["Y_p"]), ("T0", c.T0, c2["T_CMB0"]), ("fac_mHemH", c.fac_mHemH, c2["fac_mHemH"]), ("zsRe", c.zsRe, c2["zsRe"]),
                              ("Nb0", c.Nb0, c2["Nb0"]), ("H0", c.H0, c2["H0"]), ("O_k", c.O_k, c2["O_k"]), ("O_L", c.O_L, c2["O_L"]), ("O_m", c.O_m, c2["O_m"]),
                              ("O_rel", c.O_rel, c2["O_rel"]), ("rho_g_gr", c.rho_g_gr, c2["rho_g_gr"]), ("sigT0", c.sigT0, c2["sigT0"]), ("hPlanck", c.hPlanck, c2["hPlanck"]),
                              ("kBoltz", c.kBoltz, c2["kBoltz"]), ("mElect", c.mElect, c2["mElect"]), ("EH_inf_ergs", c.EH_inf_ergs, cs["EH_inf_ergs"]), ("me_mp", c.me_mp, cs["me_mp"]),
                              ("me_malp", c.me_malp, cs["me_malp"]), ("pi", c.pi, cs["PI"]), ("fHe", cosmos_fHe(CosmosAccessors(c, nothing, nothing)), p.rec["COS"]["fHe"]),
                              ("theta.F", θ[1], cfg["Recfast_F"]), ("theta.A2s1s", θ[2], cfg["A2s1s"]), ("theta.Yp", θ[3], cfg["Yp"]), ("theta.Omega_b", θ[4], cfg["Omega_b"]),
                              ("theta.h100", θ[5], cfg["h100"]), ("theta.T0", θ[6], cfg["T0"]))
                u = ulp35(a, b); rec35("A1:ulp:$k", u)
                @test u <= 2.0                                    # observed: 0 ulp for every field at every point (glibc pow on Linux)
            end
        end
        # the fiducial closure equals the fixture-built constants of the 4b route field by field
        @test CL35["fid"].inp[1:4] == [67.36, 0.02237, 0.12, 0.24568275240335613]
        c = cosmos_constants(; Omega_c = DFID35.Omega_c, Omega_b = DFID35.Omega_b, Omega_k = 0.0, N_eff = 3.046, h100 = 0.6736, T_CMB = 2.7255, Y_p = 0.24568275240335613)
        c4 = constants4b()
        for f in fieldnames(CosmosConstants); @test ulp35(getfield(c, f), getfield(c4, f)) <= 2.0; end
        @test recfast_theta(; Omega_b = DFID35.Omega_b, h100 = 0.6736, T_CMB = 2.7255, Y_p = 0.24568275240335613) == theta4c()
    end

    @testset "A2/A8: Hubble loader last-node treatment and analytic fallback" begin
        z = [1.0e4 - (i - 1) * (1.0e4 / 9999) for i in 1:10000]; Hz = [2.0e-18 * (1.0 + zz)^1.5 for zz in z]
        z[end] = 0.0; @test hubble_table(z, Hz).zmin == exp(log(1.0e-10))
        z[end] = 8.3e-13; @test hubble_table(z, Hz).zmin == exp(log(8.3e-13))       # CAMB -march=native grid: last node not exactly 0 (Cosmos.cpp:125-126)
        c = constants4b(); a = CosmosAccessors(c, nothing, hubble_table(HUB4B[:, 1], HUB4B[:, 2]))
        for zq in (1.0e4, 1.2e4, 2.5e4)                      # outside (zmin, zmax): the analytic closure of the constants
            zp1 = 1.0 + zq; @test cosmos_H(a, zq) == c.H0 * sqrt(c.O_L + zp1 * zp1 * (c.O_k + zp1 * (c.O_m + c.O_rel * zp1)))
        end
    end

    @testset "A2: preliminary history vs native (fiducial, h50)" begin
        c = constants4b(); θ = theta4c()
        pre = preliminary_history(c, θ, hubble_table(HUB4B[:, 1], HUB4B[:, 2]), SRF35)
        check_prelim35("A2:fid", pre, [vcat(Float64(i), HIST4B[i, :]) for i in 1:6000])
        p = CL35["h50"]; H0, ombh2, omch2, YHe, nnu, T0 = p.inp; d = camb_density_parameters(ombh2, omch2, H0)
        c5 = cosmos_constants(; Omega_c = d.Omega_c, Omega_b = d.Omega_b, Omega_k = 0.0, N_eff = nnu, h100 = H0 / 100.0, T_CMB = T0, Y_p = YHe)
        HB = reduce(hcat, read_rows35(fix35("hubble")))'
        pre5 = preliminary_history(c5, recfast_theta(; Omega_b = d.Omega_b, h100 = H0 / 100.0, T_CMB = T0, Y_p = YHe), hubble_table(HB[:, 1], HB[:, 2]), SRF35)
        check_prelim35("A2:h50", pre5, read_rows35(fix35("history")))
    end

    @testset "A3/A8: configuration and input checks" begin
        rp = collect(CAMB_COSMOREC_RUNPARS)
        @test check_runmode0_config(0, rp) === nothing
        @test check_runmode0_config(0, [rp[1:15]; 2.0; 3.0; 2.0; 2.0; 3.0; 2.0]) === nothing   # overrides equal to the batch defaults are neutral
        @test_throws ArgumentError check_runmode0_config(1, rp)
        for (k, v) in ((1, 1.0e-25), (2, 1.0), (2, -1.0), (3, 1.0), (4, 8.0), (5, 1.0), (8, 1.01), (11, 0.0), (12, 0.1))
            bad = copy(rp); bad[k] = v; @test_throws ArgumentError check_runmode0_config(0, bad)
        end
        for k in 16:21; bad = copy(rp); bad[k] = 7.0; @test_throws ArgumentError check_runmode0_config(0, bad); end
        @test_throws ArgumentError check_runmode0_config(0, rp[1:20])
        z = HUB4B[:, 1]
        @test check_camb_grid(z) === nothing && check_camb_grid(NAT10_35["ZREC"]) === nothing
        @test_throws ArgumentError check_camb_grid(z[1:9999])
        @test_throws ArgumentError check_camb_grid(reverse(z))
        zb = copy(z); zb[5000] += 1.0e-6; @test_throws ArgumentError check_camb_grid(zb)
        zb = copy(z); zb[end] = -1.0e-12; @test_throws ArgumentError check_camb_grid(zb)
        zb = copy(z); zb[end] = 1.0e-6; @test_throws ArgumentError check_camb_grid(zb)
        zb = copy(z); zb[end] = 8.3e-13; @test check_camb_grid(zb) === nothing
        # the entry rejects bad inputs BEFORE any work (the tables are never touched)
        dummy = Runmode0Tables(nothing, nothing)
        call(; runmode = 0, runpars = rp, Hz = HUB4B[:, 2], zH = z, zo = z, Ob = 0.05, Yhe = 0.245) =
            cosmorec_calc_h(runmode, runpars, 0.26, Ob, 0.0, 3.046, 67.36, 2.7255, Yhe, zH, Hz, zo, dummy; solve_recfast = nothing, solve_pass = nothing, solve_tail = nothing)
        @test_throws ArgumentError call(runmode = 2)
        @test_throws ArgumentError call(Hz = [HUB4B[1:(end - 1), 2]; 0.0])
        @test_throws ArgumentError call(Hz = [NaN; HUB4B[2:end, 2]])
        @test_throws DimensionMismatch call(Hz = HUB4B[1:9999, 2], zH = z)
        @test_throws ArgumentError call(zo = z[1:5000])
        @test_throws ArgumentError call(Ob = -0.01)
        @test_throws ArgumentError call(Yhe = 1.2)
    end

    @testset "end-to-end standalone entry vs native (production main-ODE callback)" begin
        # fiducial: the ABI-shaped entry with Omega's (identical to what the adapter passes) on the native ZREC grid
        Zf = NAT10_35["ZREC"]
        r = cosmorec_calc_h(0, collect(CAMB_COSMOREC_RUNPARS), DFID35.Omega_c, DFID35.Omega_b, 0.0, 3.046, 67.36, 2.7255, 0.24568275240335613,
                            HUB4B[:, 1], HUB4B[:, 2], Zf, TABLES35; solve_recfast = SRF35, solve_pass = solve_phase5, solve_tail = solve_tail5)
        @test r.theta == theta4c()
        check_e2e35("E2E:fid", r, Zf, NAT10_35["NAT_XE"], NAT10_35["NAT_TB"], 1:10000, NSW10_35)
        # h50 (H0 50, ombh2 0.025, omch2 0.16, YHe 0.245): the physical-parameter front end
        p = CL35["h50"]; H0, ombh2, omch2, YHe, nnu, T0 = p.inp
        HB = reduce(hcat, read_rows35(fix35("hubble")))'
        out = read_rows35(fix35("output")); O = reduce(hcat, out)'
        nsw5 = parse.(Float64, split(first(l for l in eachline(fix35("output")) if startswith(l, "SWITCH")), '\t')[2:end])
        Z = [1.0e4 - i * ((1.0e4 - 0.0) / 9999) for i in 0:9999]          # the native harness output grid
        @test Z[Int.(O[:, 1])] == O[:, 2]
        r5 = cosmorec_history(; ombh2 = ombh2, omch2 = omch2, omk = 0.0, N_eff = nnu, H0 = H0, T_CMB = T0, Y_He = YHe, z_Hz = HB[:, 1], Hz = HB[:, 2], z_out = Z,
                              tables = TABLES35, solve_recfast = SRF35, solve_pass = solve_phase5, solve_tail = solve_tail5)
        check_e2e35("E2E:h50", r5, Z, O[:, 3], O[:, 4], Int.(O[:, 1]), nsw5)
    end

    @testset "info" begin
        foreach(kv -> println("Chunk35 ", kv[1], " = ", kv[2]), sort!(collect(OBS35)))
    end
end
