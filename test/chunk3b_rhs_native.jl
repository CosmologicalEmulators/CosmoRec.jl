using Test
using CosmoRec
@isdefined(FIX3B) || include("chunk3b_helpers.jl")

const G3B, LG3B, CASES3B, W3B = read_fixture3b()
const TAB3B = window_table(LG3B, W3B)
const ATOM3B = NATIVE_HELIUM_ATOM

relerr(a, b) = abs(a - b) / max(abs(b), floatmin(Float64))
# Component increment tolerance is relative to the largest term magnitude of the vector (cancellation-safe).
vecerr(a, b) = maximum(abs.(a .- b)) / max(maximum(abs.(b)), floatmin(Float64))
# Sentinel-subtracted native increments carry an absolute rounding floor of about eps * 1.25e-3.
const SENT_FLOOR = 2 * eps(1.25e-3)

@testset "Chunk3b native helium RHS components" begin
    @testset "atom data / constants vs native header" begin
        @test parse(Int, grow3b(G3B, "N_resolved")[2]) == 4 && parse(Int, grow3b(G3B, "N_resolved")[4]) == 7 && parse(Int, grow3b(G3B, "N_resolved")[12]) == 5
        gl = [r for r in G3B if r[1] == "level"]
        @test length(gl) == 7
        for r in gl
            i = parse(Int, r[2]) + 1
            @test ATOM3B.gw[i] == parse(Float64, r[12])
            @test ATOM3B.nu_ion[i] == parse(Float64, r[14])
        end
        gres = [r for r in G3B if r[1] == "resolved"]
        @test [parse(Int, r[4]) + 1 for r in gres] == ATOM3B.res_index
        for r in gres
            @test parse(Float64, r[10]) == ATOM3B.mu_red
            @test parse(Int, r[12]) == 500 && parse(Int, r[14]) == 4
        end
        @test parse(Float64, grow3b(G3B, "Dnu_1s2_2S")[2]) == ATOM3B.Dnu_2s
        @test parse(Float64, grow3b(G3B, "const_HeI_A2s_1s")[2]) == ATOM3B.A2s1s
        @test parse(Float64, grow3b(G3B, "mu_red")[2]) == ATOM3B.mu_red
        for (nm, f) in (("PI", :twopi),)
            @test 2 * parse(Float64, grow3b(G3B, nm)[2]) == getfield(NATIVE_HELIUM_CONSTANTS, f)
        end
        @test parse(Float64, grow3b(G3B, "const_me_gr")[2]) == NATIVE_HELIUM_CONSTANTS.me_gr
        @test parse(Float64, grow3b(G3B, "const_kB")[2]) == NATIVE_HELIUM_CONSTANTS.kB
        @test parse(Float64, grow3b(G3B, "const_h")[2]) == NATIVE_HELIUM_CONSTANTS.h
        @test parse(Float64, grow3b(G3B, "const_h_kb")[2]) == NATIVE_HELIUM_CONSTANTS.h_kb
        @test parse(Float64, grow3b(G3B, "const_cl")[2]) == NATIVE_RHS_CONSTANTS.c
        lys = [r for r in G3B if r[1] == "Ly"]
        @test length(lys) == 3
        for (k, r) in enumerate(lys)
            @test parse(Int, r[5]) + 1 == ATOM3B.ly_index[k]
            @test parse(Float64, r[9]) == ATOM3B.ly_A21[k]
            @test parse(Float64, r[11]) == ATOM3B.ly_lambda21[k]
            @test parse(Float64, r[13]) == ATOM3B.ly_nu21[k]
            @test parse(Float64, r[17]) == ATOM3B.ly_w[k]
        end
        @test ATOM3B.intercombination_index == parse(Int, grow3b(G3B, "N_resolved")[12]) + 1
        @test parse(Int, grow3b(G3B, "N_resolved")[6]) == 0       # Diffusion_correction_HeI_is_on
        @test parse(Int, grow3b(G3B, "N_resolved")[8]) == 0       # HeISTfeedback
        @test length(LG3B) == 500
    end

    @testset "native sentinel and stencil/windows" begin
        for c in CASES3B
            @test CosmoRec.helium_stencil_start(TAB3B, log(c.Tg)) == c.lx + 1
        end
    end

    @testset "rate lookup A, B, R vs native get_rates_HeI" begin
        for c in CASES3B
            r = get_helium_rates(TAB3B, c.Tg)
            @test all(relerr.(r.A, c.A) .< 2e-13)
            @test all(relerr.(r.B, c.B) .< 2e-13)
            for m in 1:4, i in 1:4
                if i <= m
                    @test r.R[m, i] == 0.0 && c.R[m, i] == 0.0
                else
                    @test relerr(r.R[m, i], c.R[m, i]) < 2e-13
                end
            end
        end
    end

    sent = sentinel3b()
    function inc(c, key)
        a = c.after[key]
        return a .- sent
    end

    @testset "components vs native, accumulation (+=) from sentinel" begin
        for c in CASES3B
            r = get_helium_rates(TAB3B, c.Tg)
            X = c.X
            # continuum (Rci/Ric)
            dX = copy(sent); helium_continuum!(dX, c.Xe, c.NH * c.XHeII, X, c.A, c.B)
            @test vecerr(dX, c.after["RCI"]) < 1e-14
            dX = copy(sent); helium_continuum!(dX, c.Xe, c.NH * c.XHeII, X, r.A, r.B)
            @test vecerr(dX, c.after["RCI"]) < 2e-13
            # interlevel
            dX = copy(sent); helium_interlevel!(dX, c.Tg, X, c.R)
            @test vecerr(dX, c.after["RIJ"]) < 1e-14
            # two-photon
            dX = copy(sent); helium_two_photon!(dX, c.Tg, X)
            @test vecerr(dX, c.after["2PH"]) < 1e-14
            # Lyman channels k = 1..3 (indices 3, 4, 6)
            for k in 1:3
                idx = ATOM3B.ly_index[k]
                dX = copy(sent)
                lyman_channel!(dX, idx, c.Tg, X[1], X[idx], c.NH, c.Hz, ATOM3B.ly_A21[k], ATOM3B.ly_lambda21[k], ATOM3B.ly_nu21[k], NATIVE_RHS_CONSTANTS, ATOM3B.ly_w[k])
                @test vecerr(dX, c.after["LY$k"]) < 1e-14
            end
            # combined base RHS
            dX = copy(sent); helium_base_rhs!(dX, c.Tg, c.Xe, c.NH, c.Hz, X, c.fHe, c.A, c.B, c.R)
            @test vecerr(dX, c.after["BASE"]) < 1e-14
            dX = copy(sent); helium_base_rhs!(dX, TAB3B, c.Tg, c.Xe, c.NH, c.Hz, X, c.fHe)
            @test vecerr(dX, c.after["BASE"]) < 3e-13
            # untouched slots: unresolved triplet P0, P2 never change
            @test dX[5] == sent[5] && dX[7] == sent[7]
        end
    end

    @testset "increments (zero-start Julia) vs native after - sentinel, table-driven" begin
        worst = 0.0
        for c in CASES3B
            r = get_helium_rates(TAB3B, c.Tg)
            comps = Dict("RCI" => d -> helium_continuum!(d, c.Xe, c.NH * c.XHeII, c.X, r.A, r.B),
                "RIJ" => d -> helium_interlevel!(d, c.Tg, c.X, r.R), "2PH" => d -> helium_two_photon!(d, c.Tg, c.X),
                "BASE" => d -> helium_base_rhs!(d, TAB3B, c.Tg, c.Xe, c.NH, c.Hz, c.X, c.fHe))
            for k in 1:3
                idx = ATOM3B.ly_index[k]
                comps["LY$k"] = d -> lyman_channel!(d, idx, c.Tg, c.X[1], c.X[idx], c.NH, c.Hz, ATOM3B.ly_A21[k], ATOM3B.ly_lambda21[k], ATOM3B.ly_nu21[k], NATIVE_RHS_CONSTANTS, ATOM3B.ly_w[k])
            end
            for (k, f) in comps
                d = zeros(NLHE); f(d)
                nat = inc(c, k)
                tol = 4 * SENT_FLOOR + 2e-12 * maximum(abs.(nat))
                err = maximum(abs.(d .- nat))
                @test err <= tol
                worst = max(worst, err / max(maximum(abs.(nat)), 1e-10))
            end
        end
        @info "chunk3b worst increment error / max(|native increment|, 1e-10)" worst
    end

    @testset "sum rules: component increments sum to base and conserve helium" begin
        for c in CASES3B
            parts = [inc(c, k) for k in ("RCI", "RIJ", "2PH", "LY1", "LY2", "LY3")]
            @test maximum(abs.(sum(parts) .- inc(c, "BASE"))) <= 8 * SENT_FLOOR + 1e-12 * maximum(abs.(inc(c, "BASE")))
            # RIJ, 2PH, LY conserve sum of populations (transfer only); only RCI changes total (continuum exchange)
            for k in ("RIJ", "2PH", "LY1", "LY2", "LY3")
                v = inc(c, k)
                @test abs(sum(v)) <= 8 * SENT_FLOOR + 1e-12 * maximum(abs.(v))
            end
            v = inc(c, "RCI")
            @test v[1] == 0.0 && v[5] == 0.0 && v[7] == 0.0       # exact: resolved indices only
        end
    end

    @testset "spin_forbidden = false drops only the intercombination channel" begin
        for c in CASES3B[1:6]
            dXa = copy(sent); helium_base_rhs!(dXa, c.Tg, c.Xe, c.NH, c.Hz, c.X, c.fHe, c.A, c.B, c.R; spin_forbidden = false)
            dXb = copy(sent); helium_base_rhs!(dXb, c.Tg, c.Xe, c.NH, c.Hz, c.X, c.fHe, c.A, c.B, c.R; spin_forbidden = true)
            @test maximum(abs.((dXb .- dXa) .- inc(c, "LY3"))) <= 4 * SENT_FLOOR + 1e-12 * maximum(abs.(inc(c, "LY3")))
        end
    end

    @testset "domain errors" begin
        @test_throws RateTableDomainError get_helium_rates(TAB3B, exp(LG3B[1]) * (1 - 1e-9))
        @test_throws RateTableDomainError get_helium_rates(TAB3B, exp(LG3B[end]) * (1 + 1e-9))
        @test_throws RateTableDomainError get_helium_rates(TAB3B, exp(LG3B[end - 2]) * (1 + 1e-9))
    end
end

@testset "Chunk3b hand-derived invariants (independent of the helper implementations)" begin
    kB, h, me, twopi = 1.3806504000000002e-16, 6.62606896e-27, 9.1093821500000007e-28, 6.2831853071795864769
    hkb = 4.7992373449498863e-11
    mu = ATOM3B.mu_red
    qe(Tg) = (twopi * me * mu * kB * Tg / h / h)^1.5
    saha(i, Tg) = (ATOM3B.gw[i] / 4) * exp(hkb * ATOM3B.nu_ion[i] / Tg) / qe(Tg)
    pts = [c for c in CASES3B if (startswith(c.label, "phys_") && c.Tg >= 4000) || c.label == "mid_interval"]
    @testset "Saha/Boltzmann equilibrium annihilates every base term" begin
        for c in pts
            Xe, Nc = c.Xe, c.NH * c.XHeII
            X = [saha(i, c.Tg) * Xe * Nc for i in 1:7]
            r = get_helium_rates(TAB3B, c.Tg)
            scale = maximum(r.B[m] * X[ATOM3B.res_index[m]] for m in 1:4) + maximum(abs.(r.R)) * maximum(X)
            for (nm, f) in (("cont", d -> helium_continuum!(d, Xe, Nc, X, r.A, r.B)), ("rij", d -> helium_interlevel!(d, c.Tg, X, r.R)),
                    ("2ph", d -> helium_two_photon!(d, c.Tg, X)), ("ly", d -> helium_lyman!(d, c.Tg, X, c.NH, c.Hz)),
                    ("base", d -> helium_base_rhs!(d, TAB3B, c.Tg, Xe, c.NH, c.Hz, X, c.XHeII + X[1])))
                d = zeros(7); f(d)
                @test maximum(abs.(d)) <= 1e-11 * maximum(abs.(X)) * maximum(vcat(r.B, [51.3, 1.8e9, maximum(abs.(r.R))])) + 1e-300
            end
        end
    end
    @testset "transfer terms are antisymmetric; only the continuum changes the helium total" begin
        for c in pts
            r = get_helium_rates(TAB3B, c.Tg)
            X = c.X .* [1.0, 1.3, 0.7, 1.1, 0.9, 1.2, 0.8]
            dtr = zeros(7); helium_interlevel!(dtr, c.Tg, X, r.R); helium_two_photon!(dtr, c.Tg, X); helium_lyman!(dtr, c.Tg, X, c.NH, c.Hz)
            @test abs(sum(dtr)) <= 1e-12 * sum(abs, dtr) + 1e-300
            dc = zeros(7); helium_continuum!(dc, c.Xe, c.NH * c.XHeII, X, r.A, r.B)
            dall = zeros(7); helium_base_rhs!(dall, TAB3B, c.Tg, c.Xe, c.NH, c.Hz, X, c.fHe)
            @test abs(sum(dall) - sum(dc)) <= 1e-12 * sum(abs, dall) + 1e-300
            @test sum(dc) == sum(r.B[m] * (r.A[m] * c.Xe * c.NH * c.XHeII - X[ATOM3B.res_index[m]]) for m in 1:4) || isapprox(sum(dc), sum(r.B[m] * (r.A[m] * c.Xe * c.NH * c.XHeII - X[ATOM3B.res_index[m]]) for m in 1:4); rtol = 1e-13)
            @test dall[5] == 0.0 && dall[7] == 0.0                       # unresolved 2^3P0, 2^3P2 are decoupled from these terms
        end
    end
    @testset "two-photon and a Lyman channel hand-evaluated" begin
        Tg = 6300.0; X = [0.07, 2e-9, 3e-9, 1e-9, 0.0, 4e-9, 0.0]; NH = 1.0e3; Hz = 1.0e-13
        d = zeros(7); helium_two_photon!(d, Tg, X)
        @test d[1] == 51.3 * (X[2] - X[1] * exp(-hkb * 4984872582556676.0 / Tg)) && d[2] == -d[1]
        A21, lam, nu = 1798900000.0, 5.8433431808919833e-06, 5130495483823985.0
        tau = A21 * lam^3 / (8 * pi * Hz) * (X[1] * NH * 3 - X[3] * NH); p = (1 - exp(-tau)) / tau; ef = exp(-hkb * nu / Tg)
        dl = zeros(7); lyman_channel!(dl, 3, Tg, X[1], X[3], NH, Hz, A21, lam, nu, NATIVE_RHS_CONSTANTS, 3.0)
        @test isapprox(dl[1], p * A21 / (1 - ef) * (X[3] - 3 * X[1] * ef); rtol = 1e-14)
    end
    @testset "table branches: node values, one-sided seam values, low/high stencil ends" begin
        t = TAB3B
        for k in (100, 250, 400)
            for m in 1:4
                Tk = exp(LG3B[k + 1])
                r = get_helium_rates(t, Tk)
                @test isapprox(r.B[m], exp(t.B[k + 1, m]); rtol = 1e-11)
                for i in (m + 1):4
                    @test isapprox(r.R[m, i], exp(t.R[k + 1, m, i]); rtol = 1e-11)
                end
            end
            lo = get_helium_rates(t, exp(LG3B[k + 1]) * (1 - 1e-9)); hi = get_helium_rates(t, exp(LG3B[k + 1]) * (1 + 1e-9))
            @test all(isapprox.(lo.B, hi.B; rtol = 1e-6))                        # near-continuous across the stencil switch (derivative may jump)
            @test CosmoRec.helium_stencil_start(t, LG3B[k + 1] - 1e-9) == k && CosmoRec.helium_stencil_start(t, LG3B[k + 1] + 1e-9) == k + 1
        end
        @test CosmoRec.helium_stencil_start(t, LG3B[1]) == 1                           # lowest allowed Tg
        @test CosmoRec.helium_stencil_start(t, LG3B[497] - 1e-9) == 497 - 1            # highest allowed: j + 3 = N - 1... (see docs)
        @test isfinite(helium_A(t.gw[1], t.nu_ion[1], t.mu_red, exp(LG3B[1]))) && isfinite(helium_A(t.gw[1], t.nu_ion[1], t.mu_red, exp(LG3B[end])))
    end
end
