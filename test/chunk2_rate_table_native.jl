# Chunk 2: pure-Julia HI effective-rate lookup vs ORIGINAL CosmoRec get_rates() text fixtures (no AD here; see
# chunk2_rate_table_ad.jl). Fixtures: test/fixtures/native_hrates_*.txt (provenance/attribution in their headers).
# Independent oracles used below: Newton divided-difference interpolation, closed-form detailed-balance expression,
# exact polynomial reproduction on a synthetic table. Julia-vs-native tolerance is rtol 1e-14 (observed <= 1.6e-16).
using Test
using SHA
using Random
using CosmoRec

@isdefined(FIXDIR2) || include("chunk2_helpers.jl")

# independent oracle: Newton divided differences through 4 points
function newton4(xs, ys, x)
    c = collect(float.(ys))
    for j in 2:4, i in 4:-1:j
        c[i] = (c[i] - c[i - 1]) / (xs[i] - xs[i - j + 1])
    end
    return c[1] + (x - xs[1]) * (c[2] + (x - xs[2]) * (c[3] + (x - xs[3]) * c[4]))
end

# independent closed form of exp(log_qnl_qe): (gw/2) exp(h_kb nu/T) / (2 pi me mu kB T / h^2)^(3/2), written without log/exp pairing
qratio_db(cs, gw, nu, mu, T) = (gw / 2) * exp(cs["const_h_kb"] * nu / T) / (cs["TWOPI"] * cs["const_me_gr"] * mu * cs["const_kB"] * T / cs["const_h"]^2)^1.5

const TABLE2, ROWS2, HDR2 = read_table_window2()
const QUERIES2 = read_queries2()
const CS2 = fixture_constants2(HDR2)

@testset "Chunk 2: native HI rate lookup (get_rates) fixtures" begin
    @testset "fixture integrity and provenance" begin
        tw = joinpath(FIXDIR2, "native_hrates_table_window.txt")
        q = joinpath(FIXDIR2, "native_hrates_get_rates.txt")
        for p in (tw, q)
            h = parse_header2(p)
            @test any(startswith("# cosmorec_sha: 086769055f61ae0c244a53dd381ee65b624d0ac3"), h)
            @test any(contains("Chluba & Thomas 2010"), h) && any(contains("Chluba & Sunyaev 2006"), h)
            @test any(contains("Gas_of_Atoms(nS=3,Z=1,Np=1.0,Qlines_on=false,Rec_flag=0,mflag=-2)"), h)
            @test any(contains("static GSL 2.8"), h)
            body = join(l * "\n" for l in eachline(p) if !startswith(l, "#"))
            want = first(m[1] for m in match.(r"^# records_sha256: (\w+)", h) if m !== nothing)
            @test bytes2hex(sha256(body)) == want
        end
        @test length(QUERIES2) == 47
        @test all(length(r.A) == 5 for r in QUERIES2)
        @test TABLE2.states.n == [2, 2, 3, 3, 3] && TABLE2.states.l == [0, 1, 0, 1, 2]
        @test TABLE2.states.gw == [2.0, 6.0, 2.0, 6.0, 10.0]
        @test TABLE2.states.gw == 2 .* (2 .* TABLE2.states.l .+ 1)
        @test size(TABLE2.A) == (72, 41, 5) && size(TABLE2.R) == (72, 5, 5) && size(TABLE2.B) == (72, 5)
        @test issorted(ROWS2) && ROWS2[1] == 0 && ROWS2[end] == 499
        @test issorted(TABLE2.lgTg; lt = <=) && issorted(TABLE2.lgrho; lt = <=)
        @test length(TABLE2.lgrho) == 41 && TABLE2.lgrho[1] == -2.3025851 && TABLE2.lgrho[end] == 0.09531018
        @test TABLE2.lgTg[1] == 3.4003637 && TABLE2.lgTg[end] == 9.1632723
        @test NATIVE_CONSTANTS.twopi == CS2["TWOPI"] && NATIVE_CONSTANTS.me_gr == CS2["const_me_gr"] && NATIVE_CONSTANTS.kB == CS2["const_kB"] &&
              NATIVE_CONSTANTS.h == CS2["const_h"] && NATIVE_CONSTANTS.h_kb == CS2["const_h_kb"]
        # deterministic loading: re-reading gives identical arrays
        T2, rows2, _ = read_table_window2()
        @test T2.A == TABLE2.A && T2.B == TABLE2.B && T2.R == TABLE2.R && T2.logq_knot == TABLE2.logq_knot && rows2 == ROWS2
    end

    @testset "Julia vs ORIGINAL get_rates, every query/state/entry" begin
        worst = 0.0; nbit = 0; ntot = 0
        for q in QUERIES2
            r = get_rates(TABLE2, q.Tg, q.Te)
            @test size(r.R) == (5, 5)
            for m in 1:5
                pairs = vcat([(r.A[m], q.A[m]), (r.B[m], q.B[m])], [(r.R[m, i], q.R[m][i]) for i in 1:5])
                for (x, y) in pairs
                    ntot += 1
                    if isnan(y)
                        @test isnan(x)
                    elseif y == 0
                        @test x === 0.0
                    else
                        rel = abs(x - y) / abs(y)
                        worst = max(worst, rel); nbit += (x == y)
                        @test rel <= 1.0e-14
                    end
                end
                @test all(r.R[m, 1:m] .== 0)
            end
        end
        @info "Julia vs native get_rates" entries = ntot bit_identical = nbit worst_rel = worst
        @test ntot == 47 * 5 * 7
    end

    @testset "query order independence and in-place variant" begin
        ref = [get_rates(TABLE2, q.Tg, q.Te) for q in QUERIES2]
        for perm in (reverse(eachindex(QUERIES2)), randperm(MersenneTwister(5), length(QUERIES2)))
            for i in perm
                r = get_rates(TABLE2, QUERIES2[i].Tg, QUERIES2[i].Te)
                @test isequal(r.A, ref[i].A) && isequal(r.B, ref[i].B) && isequal(r.R, ref[i].R)
            end
        end
        A = zeros(5); B = zeros(5); R = fill(NaN, 5, 5)
        q = QUERIES2[1]
        @test get_rates!(A, B, R, TABLE2, q.Tg, q.Te) === nothing
        @test A == ref[1].A && B == ref[1].B && R == ref[1].R
        @test_throws DimensionMismatch get_rates!(zeros(4), B, R, TABLE2, q.Tg, q.Te)
    end

    @testset "interior values vs independent Newton interpolation (own oracle, native table values)" begin
        for q in QUERIES2
            startswith(q.label, "int_") || continue
            Tg, Te = q.Tg, q.Te
            x, y = log(Tg), log(Te / Tg)
            jx = searchsortedlast(TABLE2.lgTg, x); lx = max(jx - 1, 1)
            jy = searchsortedlast(TABLE2.lgrho, y); ly = max(jy - 1, 1)
            r = get_rates(TABLE2, Tg, Te)
            for m in 1:5
                @test r.B[m] ≈ exp(newton4(TABLE2.lgTg[lx:lx+3], TABLE2.B[lx:lx+3, m], x)) rtol = 1.0e-12
                for i in (m + 1):5
                    @test r.R[m, i] ≈ exp(newton4(TABLE2.lgTg[lx:lx+3], TABLE2.R[lx:lx+3, m, i], x)) rtol = 1.0e-12
                end
                isfinite(r.A[m]) || continue
                # A: bicubic of (A - B - logq_knot), logq_knot evaluated with the closed form
                if Tg / 2.725 - 1 <= 2000
                    inner = [newton4(TABLE2.lgTg[lx:lx+3],
                                     [TABLE2.A[k, jj, m] - TABLE2.B[k, m] -
                                      log(qratio_db(CS2, TABLE2.states.gw[m], TABLE2.states.nuion[m], TABLE2.states.mu_red[m], exp(TABLE2.lgTg[k])))
                                      for k in lx:lx+3], x) for jj in ly:ly+3]
                    fxy = newton4(TABLE2.lgrho[ly:ly+3], inner, y)
                    abs(expm1(fxy)) <= 1.0e-4 && (fxy = 0.0)
                    @test r.A[m] ≈ exp(fxy) * qratio_db(CS2, TABLE2.states.gw[m], TABLE2.states.nuion[m], TABLE2.states.mu_red[m], Tg) rtol = 1.0e-9
                end
            end
        end
    end

    @testset "exact knots and seams" begin
        N = length(TABLE2.lgTg)
        for k in (5, 20, 40, N - 5)
            x = TABLE2.lgTg[k]
            lx, _ = stencil_start(TABLE2, x, -0.3)
            w = CosmoRec.lagrange_weights(TABLE2.lgTg, lx, x)
            kk = k - lx + 1
            @test w[kk] == 1.0 && all(w[i] == 0.0 for i in 1:4 if i != kk)       # exact knot reproduction (weights)
        end
        # consecutive native-grid knots: weights reproduce the stored table entry for B (exact knot value)
        for rows in ((1, 2, 3, 4), (N - 3, N - 2, N - 1, N))
            for k in rows
                x = TABLE2.lgTg[k]
                lx = k <= 2 ? 1 : (k >= N - 2 ? N - 3 : k - 1)
                w = CosmoRec.lagrange_weights(TABLE2.lgTg, lx, x)
                @test sum(w[i] * TABLE2.B[lx + i - 1, 1] for i in 1:4) == TABLE2.B[k, 1]
            end
        end
        # continuity of the interpolant across a stencil seam (value; derivative jump is documented in the AD tests)
        for k in (20, 40)
            Tk = exp(TABLE2.lgTg[k])
            lo = get_rates(TABLE2, prevfloat(Tk, 4), 0.9 * Tk); hi = get_rates(TABLE2, nextfloat(Tk, 4), 0.9 * Tk)
            @test lo.B ≈ hi.B rtol = 1.0e-9
            @test lo.A ≈ hi.A rtol = 1.0e-9
        end
        # stencil clamping at the low edge: cells 1 and 2 share lx = 1; lowest cell does not extend below the table
        @test stencil_start(TABLE2, TABLE2.lgTg[1], TABLE2.lgrho[1]) == (1, 1)
        @test stencil_start(TABLE2, TABLE2.lgTg[2] - 1.0e-9, 0.0)[1] == 1
        @test stencil_start(TABLE2, TABLE2.lgTg[2], 0.0)[1] == 1
        @test stencil_start(TABLE2, TABLE2.lgTg[3], 0.0)[1] == 2
    end

    @testset "detailed-balance branch: both sides, closed form" begin
        thr = 5452.7250000000013            # smallest double with Tg/2.725 - 1 > 2000 (native fixture header)
        @test thr / 2.725 - 1.0 > 2000.0 && !(prevfloat(thr) / 2.725 - 1.0 > 2000.0)
        for q in QUERIES2
            startswith(q.label, "db_") || continue
            r = get_rates(TABLE2, q.Tg, q.Te)
            if q.Tg >= thr
                for m in 1:5
                    @test r.A[m] ≈ qratio_db(CS2, TABLE2.states.gw[m], TABLE2.states.nuion[m], TABLE2.states.mu_red[m], q.Tg) rtol = 1.0e-12
                    @test q.A[m] ≈ qratio_db(CS2, TABLE2.states.gw[m], TABLE2.states.nuion[m], TABLE2.states.mu_red[m], q.Tg) rtol = 1.0e-12   # native itself
                end
            else
                # below the switch the table value applies; at Te = Tg it coincides with detailed balance only through the eps clamp
                if q.Te == q.Tg
                    for m in 1:5
                        @test r.A[m] ≈ qratio_db(CS2, TABLE2.states.gw[m], TABLE2.states.nuion[m], TABLE2.states.mu_red[m], q.Tg) rtol = 1.0e-12
                    end
                else
                    @test any(abs(r.A[m] / qratio_db(CS2, TABLE2.states.gw[m], TABLE2.states.nuion[m], TABLE2.states.mu_red[m], q.Tg) - 1) > 1.0e-3 for m in 1:5)
                end
            end
        end
        a = get_rates(TABLE2, prevfloat(thr), 0.9 * prevfloat(thr)).A
        b = get_rates(TABLE2, thr, 0.9 * thr).A
        @test all(b .< a)                    # table value (Te != Tg) is larger than detailed balance: switch is discontinuous for Te < Tg
    end

    @testset "eps_A_effective clamp: both sides" begin
        # clamped  <=> fxy forced to 0 <=> A equals the closed-form detailed-balance value (independent of the table);
        # unclamped <=> |exp(fxy)-1| > 1e-4, so A/db - 1 = exp(fxy) - 1 exceeds 1e-4.
        cls(r, m, q) = begin
            db = qratio_db(CS2, TABLE2.states.gw[m], TABLE2.states.nuion[m], TABLE2.states.mu_red[m], q.Tg)
            dev = abs(r / db - 1)
            dev <= 1.0e-12 ? :clamped : (dev > 1.0e-4 * (1 - 1.0e-6) ? :unclamped : :ambiguous)
        end
        n_cl = 0; n_un = 0; near_un = 0
        for q in QUERIES2
            startswith(q.label, "eps_") || continue
            r = get_rates(TABLE2, q.Tg, q.Te)
            for m in 1:5
                isfinite(r.A[m]) || continue
                cj, cn = cls(r.A[m], m, q), cls(q.A[m], m, q)
                @test cj == cn && cj != :ambiguous
                n_cl += cj == :clamped; n_un += cj == :unclamped
                if cj == :unclamped
                    db = qratio_db(CS2, TABLE2.states.gw[m], TABLE2.states.nuion[m], TABLE2.states.mu_red[m], q.Tg)
                    near_un += abs(r.A[m] / db - 1) < 3.0e-4        # straddles the threshold from the unclamped side
                end
            end
        end
        @test n_cl > 0 && n_un > 0 && near_un > 0
        @info "eps clamp coverage" clamped = n_cl unclamped = n_un unclamped_within_3e-4 = near_un
    end

    @testset "native NaN/low-Tg region and positivity" begin
        for q in QUERIES2
            r = get_rates(TABLE2, q.Tg, q.Te)
            for m in 1:5
                overflow = (NATIVE_CONSTANTS.h_kb * TABLE2.states.nuion[m] / q.Tg > log(floatmax(Float64)))
                if overflow
                    @test isnan(r.A[m]) && isnan(q.A[m])
                else
                    @test isfinite(r.A[m]) && r.A[m] > 0 && !isnan(q.A[m])
                end
                @test isfinite(r.B[m]) && r.B[m] > 0
                @test all(isfinite, r.R[m, :]) && all(r.R[m, (m + 1):5] .> 0)
            end
        end
    end

    @testset "domain and failure semantics (native probes)" begin
        lines = [split(l, '\t') for l in eachline(joinpath(FIXDIR2, "native_hrates_domain_probes.txt")) if !startswith(l, "#")]
        header, probes = lines[1], lines[2:end]
        @test header[1] == "label" && length(probes) == 10
        for f in probes
            Tg, Te = parse(Float64, f[2]), parse(Float64, f[3])
            code, cls = parse(Int, f[4]), f[5]
            want = cls == "out_of_table_message_then_exit0" ? :out_of_table : :stencil_exceeds_table
            @test code == 139 || cls == "out_of_table_message_then_exit0" || startswith(f[1], "UBband_rho")
            err = try get_rates(TABLE2, Tg, Te); nothing catch e; e end
            @test err isa RateTableDomainError && err.reason == want
        end
        for bad in ((NaN, 1.0), (1000.0, NaN), (-1.0, 1.0), (0.0, 0.0), (Inf, 1.0), (1000.0, -1.0))
            err = try get_rates(TABLE2, bad...); nothing catch e; e end
            @test err isa RateTableDomainError && err.reason == :not_finite
        end
        @test sprint(showerror, RateTableDomainError(:out_of_table, 1.0, 2.0)) == "RateTableDomainError(out_of_table): (Tg, Te) = (1.0, 2.0)"
    end

    @testset "exact polynomial reproduction on a synthetic table (independent oracle)" begin
        xs = [5.7 + 0.2 * i + 0.03 * sin(i) for i in 1:14]
        ys = [-2.3 + 0.18 * j + 0.02 * cos(j) for j in 1:10]
        st = ResolvedStates([2, 3], [0, 1], [2.0, 6.0], [8.22012807922942e14, 3.6533902574352988e14], [0.99945567942448077, 0.99945567942448077])
        PB(x, m) = 0.3 - 1.2x + 0.4x^2 - 0.05x^3 + 0.1m
        PR(x, m, k) = 0.5k - 0.7x + 0.03x^2 + 0.01m * x^3
        f(x, y) = 0.5 + 0.01x + 0.05y + 0.003x * y^2 + 0.00002x^3 * y^3
        @test minimum(abs(expm1(f(x, y))) for x in xs[1]:0.05:xs[end - 2], y in ys[1]:0.05:ys[end - 2]) > 1.0e-2   # keeps clear of the 1e-4 clamp
        N, M, nres, neq = 14, 10, 2, 3
        B = [PB(xs[i], m) for i in 1:N, m in 1:nres]
        R = [PR(xs[i], m, k) for i in 1:N, m in 1:nres, k in 1:neq]
        tmp = AtomicRateTable(xs, ys, B, R, zeros(N, M, nres), st)
        A = [B[i, m] + tmp.logq_knot[i, m] + f(xs[i], ys[j]) for i in 1:N, j in 1:M, m in 1:nres]
        T = AtomicRateTable(xs, ys, B, R, A, st)
        rng = MersenneTwister(2026)
        nq = 0
        for cell in 1:(N - 3), _ in 1:6
            x = xs[cell] + rand(rng) * (xs[cell + 1] - xs[cell]); y = ys[1] + rand(rng) * (ys[M - 2] - ys[1])
            Tg = exp(x); Te = Tg * exp(y)
            r = get_rates(T, Tg, Te); nq += 1
            xl, yl = log(Tg), log(Te / Tg)
            for m in 1:nres
                @test r.B[m] ≈ exp(PB(xl, m)) rtol = 1.0e-12
                logq = log(st.gw[m] / 2) + NATIVE_CONSTANTS.h_kb * st.nuion[m] / Tg -
                       1.5 * log(NATIVE_CONSTANTS.twopi * NATIVE_CONSTANTS.me_gr * st.mu_red[m] * NATIVE_CONSTANTS.kB * Tg / NATIVE_CONSTANTS.h^2)
                @test r.A[m] ≈ exp(f(xl, yl) + logq) rtol = 1.0e-10
                for k in 1:neq
                    @test r.R[m, k] ≈ (k > m ? exp(PR(xl, m, k)) : 0.0) rtol = 1.0e-12
                end
            end
        end
        @test nq == 66
    end

    @testset "native file reader (hand-written 2x3 table)" begin
        txt = """
        2 3 2
        0.5 1.5
        -1.0 0.0 1.0
        11.0 12.0  21.0 22.0   31.0 32.0 33.0
        41.0 42.0  51.0 52.0   61.0 62.0 63.0
        """
        d = read_native_rate_file(IOBuffer(txt))
        @test d.lgTg == [0.5, 1.5] && d.lgrho == [-1.0, 0.0, 1.0]
        @test d.B == [11.0, 41.0] && d.Bitot == [12.0, 42.0]
        @test d.R == [21.0 22.0; 51.0 52.0] && d.A == [31.0 32.0 33.0; 61.0 62.0 63.0]
        @test_throws ArgumentError read_native_rate_file(IOBuffer(txt * "7.0\n"))
    end
end
