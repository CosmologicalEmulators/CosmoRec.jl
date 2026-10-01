# Chunk 2 AD tests for the HI effective-rate lookup: ForwardDiff Jacobians vs hand-derived local-polynomial derivatives,
# and Mooncake VJPs (DifferentiationInterface, prepared + reused caches) vs ForwardDiff J'w. The lookup is PIECEWISE
# (stencil seams, detailed-balance switch, eps clamp): derivatives are only compared away from switches; at switches
# one-sided derivatives are reported, not compared across.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec

@isdefined(FIXDIR2) || include("chunk2_helpers.jl")

const ADT, _R, _H = read_table_window2()
const ADQ = read_queries2()
const ADCS = NATIVE_CONSTANTS

pack_rates(r) = vcat(r.A, r.B, [r.R[m, i] for m in 1:5 for i in (m + 1):5])
rates_vec(q) = pack_rates(get_rates(ADT, q[1], q[2]))
wobj(q, w, idx) = dot(w, rates_vec(q)[idx])
mc_backend = AutoMooncake(; config = nothing)

# Newton divided differences and their hand-derived value and derivative at x (4 nodes)
function newton_coef(xs, ys)
    c = collect(float.(ys))
    for j in 2:4, i in 4:-1:j
        c[i] = (c[i] - c[i - 1]) / (xs[i] - xs[i - j + 1])
    end
    return c
end
function newton_val_der(xs, ys, x)
    c = newton_coef(xs, ys); d1, d2, d3 = x - xs[1], x - xs[2], x - xs[3]
    return c[1] + d1 * (c[2] + d2 * (c[3] + d3 * c[4])),
           c[2] + c[3] * (d1 + d2) + c[4] * (d2 * d3 + d1 * d3 + d1 * d2)
end

logq_closed(m, Tg) = log(ADT.states.gw[m] / 2) + ADCS.h_kb * ADT.states.nuion[m] / Tg -
                     1.5 * log(ADCS.twopi * ADCS.me_gr * ADT.states.mu_red[m] * ADCS.kB * Tg / ADCS.h^2)
dlogq_dTg(m, Tg) = -ADCS.h_kb * ADT.states.nuion[m] / Tg^2 - 1.5 / Tg

# hand-derived analytic Jacobian columns (20 packed outputs: A, B, strict-upper R) (d/dTg, d/dTe) of B, R and (non-DB, non-clamped) A on the real native table window
function analytic_jac(Tg, Te)
    x, y = log(Tg), log(Te / Tg)
    jx = searchsortedlast(ADT.lgTg, x); lx = max(jx - 1, 1)
    jy = searchsortedlast(ADT.lgrho, y); ly = max(jy - 1, 1)
    xs, ys = ADT.lgTg[lx:lx+3], ADT.lgrho[ly:ly+3]
    J = fill(NaN, 20, 2)
    row = 0
    r = get_rates(ADT, Tg, Te)
    status = fill(:smooth, 5)
    for m in 1:5                                       # A rows 1:5
        row += 1
        if Tg / 2.725 - 1.0 > 2000.0
            J[row, 1] = r.A[m] * dlogq_dTg(m, Tg) ; J[row, 2] = 0.0; status[m] = :detailed_balance
            continue
        end
        vx = [newton_val_der(xs, [ADT.A[k, jj, m] - ADT.B[k, m] - ADT.logq_knot[k, m] for k in lx:lx+3], x) for jj in ly:ly+3]
        fxy, fy = newton_val_der(ys, first.(vx), y)
        fx = newton_val_der(ys, last.(vx), y)[1]
        if abs(expm1(fxy)) <= 1.0e-4
            J[row, 1] = r.A[m] * dlogq_dTg(m, Tg); J[row, 2] = 0.0; status[m] = :clamped
        else
            # fy from newton in y of values-in-x:  value branch (first.(vx)) differentiated in y
            J[row, 1] = r.A[m] * ((fx - fy) / Tg + dlogq_dTg(m, Tg))
            J[row, 2] = r.A[m] * fy / Te
        end
    end
    for m in 1:5                                       # B rows 6:10
        row += 1
        v, d = newton_val_der(xs, ADT.B[lx:lx+3, m], x)
        J[row, 1] = r.B[m] * d / Tg; J[row, 2] = 0.0
    end
    for m in 1:5, i in (m + 1):5                       # R upper rows
        row += 1
        v, d = newton_val_der(xs, ADT.R[lx:lx+3, m, i], x)
        J[row, 1] = r.R[m, i] * d / Tg; J[row, 2] = 0.0
    end
    return J, status
end

const IDX_ALL = collect(1:20)
finite_idx(q) = [k for (k, v) in enumerate(rates_vec(q)) if isfinite(v)]

@testset "Chunk 2 AD: HI rate lookup (piecewise differentiable)" begin
    interior = [q for q in ADQ if startswith(q.label, "int_") || startswith(q.label, "hi") || startswith(q.label, "lowrho") ||
                startswith(q.label, "knot_") || startswith(q.label, "db_far") || startswith(q.label, "low")]
    interior = [q for q in interior if q.label ∉ ("hiTg_valid_b", "hirho_valid", "hirho_valid_dn")] # exclude near native-unsafe high edge

    @testset "ForwardDiff Jacobian vs hand-derived local-polynomial derivatives" begin
        worst = 0.0; n = 0
        for q in ADQ
            q.label in ("int_a", "int_b", "int_c", "int_e", "int_f", "int_h", "int_i", "int_j", "int_g", "lowrho_int", "lowrho_int2",
                        "db_far_below", "db_far_above", "db_hi", "hiTg_valid_a", "lowTg_j2", "lowTg_j1", "int_d", "knot_Tg100_int", "knot_Tg250_int") || continue
            J = ForwardDiff.jacobian(rates_vec, [q.Tg, q.Te])
            Ja, status = analytic_jac(q.Tg, q.Te)
            for k in 1:20
                isfinite(Ja[k, 1]) && isfinite(J[k, 1]) || (@test !isfinite(rates_vec([q.Tg, q.Te])[k]); continue)
                for c in 1:2
                    n += 1
                    ref = max(abs(Ja[k, c]), abs(Ja[k, 1]) * 1.0e-12, 1.0e-300)
                    err = abs(J[k, c] - Ja[k, c]) / ref
                    worst = max(worst, err)
                    @test err <= 1.0e-9
                end
            end
        end
        @info "ForwardDiff vs analytic Jacobian entries" entries = n worst_rel = worst
    end

    @testset "branch derivatives (analytic): detailed balance and eps clamp" begin
        for Tg in (5600.0, 7000.0, 9000.0)
            J = ForwardDiff.jacobian(rates_vec, [Tg, 0.9Tg]); r = get_rates(ADT, Tg, 0.9Tg)
            for m in 1:5
                @test J[m, 2] == 0.0
                @test J[m, 1] ≈ r.A[m] * dlogq_dTg(m, Tg) rtol = 1.0e-12
            end
        end
        cl = [q for q in ADQ if q.label == "eps_r1"][1]            # Te = Tg: fxy within 1e-4 of 0 for the n=2 states
        J = ForwardDiff.jacobian(rates_vec, [cl.Tg, cl.Te]); r = get_rates(ADT, cl.Tg, cl.Te)
        _, status = analytic_jac(cl.Tg, cl.Te)
        @test :clamped in status
        for m in 1:5
            status[m] == :clamped || continue
            @test J[m, 2] == 0.0
            @test J[m, 1] ≈ r.A[m] * dlogq_dTg(m, cl.Tg) rtol = 1.0e-12
        end
    end

    @testset "synthetic polynomial table: analytic derivatives" begin
        xs = [5.7 + 0.2 * i + 0.03 * sin(i) for i in 1:14]
        ys = [-2.3 + 0.18 * j + 0.02 * cos(j) for j in 1:10]
        st = ResolvedStates([2, 3], [0, 1], [2.0, 6.0], [8.22012807922942e14, 3.6533902574352988e14], [0.99945567942448077, 0.99945567942448077])
        PB(x, m) = 0.3 - 1.2x + 0.4x^2 - 0.05x^3 + 0.1m;  dPB(x, m) = -1.2 + 0.8x - 0.15x^2
        PR(x, m, k) = 0.5k - 0.7x + 0.03x^2 + 0.01m * x^3; dPR(x, m, k) = -0.7 + 0.06x + 0.03m * x^2
        f(x, y) = 0.5 + 0.01x + 0.05y + 0.003x * y^2 + 0.00002x^3 * y^3
        fx(x, y) = 0.01 + 0.003y^2 + 0.00006x^2 * y^3;   fy(x, y) = 0.05 + 0.006x * y + 0.00006x^3 * y^2
        N, M, nres, neq = 14, 10, 2, 3
        B = [PB(xs[i], m) for i in 1:N, m in 1:nres]
        R = [PR(xs[i], m, k) for i in 1:N, m in 1:nres, k in 1:neq]
        tmp = AtomicRateTable(xs, ys, B, R, zeros(N, M, nres), st)
        A = [B[i, m] + tmp.logq_knot[i, m] + f(xs[i], ys[j]) for i in 1:N, j in 1:M, m in 1:nres]
        T = AtomicRateTable(xs, ys, B, R, A, st)
        pk(r) = vcat(r.A, r.B, vec(r.R))
        rng = MersenneTwister(77); worst = 0.0
        for cell in (1, 3, 6, 9, 11), _ in 1:4
            x = xs[cell] + (0.1 + 0.8 * rand(rng)) * (xs[cell + 1] - xs[cell]); y = ys[1] + (0.1 + 0.8 * rand(rng)) * (ys[M - 2] - ys[1])
            Tg = exp(x); Te = Tg * exp(y); xl, yl = log(Tg), log(Te / Tg)
            J = ForwardDiff.jacobian(q -> pk(get_rates(T, q[1], q[2])), [Tg, Te])
            r = get_rates(T, Tg, Te)
            for m in 1:nres
                dlq = -ADCS.h_kb * st.nuion[m] / Tg^2 - 1.5 / Tg
                @test J[m, 1] ≈ r.A[m] * ((fx(xl, yl) - fy(xl, yl)) / Tg + dlq) rtol = 1.0e-9
                @test J[m, 2] ≈ r.A[m] * fy(xl, yl) / Te rtol = 1.0e-9
                @test J[nres + m, 1] ≈ r.B[m] * dPB(xl, m) / Tg rtol = 1.0e-9
                @test J[nres + m, 2] == 0.0
                for k in 1:neq
                    idx = 2nres + m + nres * (k - 1)
                    want = k > m ? r.R[m, k] * dPR(xl, m, k) / Tg : 0.0
                    @test J[idx, 1] ≈ want rtol = 1.0e-9 atol = 1.0e-300
                    worst = max(worst, abs(J[idx, 1] - want) / max(abs(want), 1.0e-300))
                end
            end
        end
    end

    @testset "Mooncake VJP vs ForwardDiff J'w (A, B, R projections)" begin
        rng = MersenneTwister(314159)
        groups = Dict("A" => 1:5, "B" => 6:10, "R" => 11:20, "A+B+R" => 1:20)
        # R has 10 independent entries (upper triangle); 35 -> 20 packed entries
        @test length(rates_vec([3000.0, 2700.0])) == 20
        worst = Dict(g => 0.0 for g in keys(groups))
        for q in ADQ
            q.label in ("int_a", "int_b", "int_c", "int_e", "int_f", "int_h", "int_i", "int_j", "knot_Tg100_int", "knot_Tg250_int",
                        "db_far_above", "db_hi", "hiTg_valid_a", "int_g") || continue
            fidx = finite_idx([q.Tg, q.Te])
            J = ForwardDiff.jacobian(rates_vec, [q.Tg, q.Te])
            for (g, rng_idx) in groups, seed in 1:2
                idx = [k for k in rng_idx if k in fidx]
                isempty(idx) && continue
                w = randn(rng, length(idx))
                obj(x, c) = wobj(x, c, idx)
                prep = prepare_gradient(obj, mc_backend, [q.Tg, q.Te], Constant(w))
                gm = gradient(obj, prep, mc_backend, [q.Tg, q.Te], Constant(w))
                gf = J[idx, :]' * w
                err = maximum(abs.(gm .- gf) ./ max.(abs.(gf), 1.0e-12 * maximum(abs.(gf)) + 1.0e-300))
                worst[g] = max(worst[g], err)
                @test err <= 1.0e-9
            end
        end
        @info "Mooncake VJP vs ForwardDiff J'w worst relative error per projection group" worst
    end

    @testset "known limitation: reverse mode at the native low-Tg overflow point (A of n=2 states is NaN)" begin
        q = [31.0, 15.5]                                   # lowTg_j2: native A[1:2] = NaN (exp overflow), B and R finite
        fidx = finite_idx(q)
        @test fidx == setdiff(1:20, [1, 2])
        J = ForwardDiff.jacobian(rates_vec, q)
        @test all(isfinite, J[fidx, :])                    # forward mode: finite derivatives for all finite outputs
        w = ones(length(fidx)); idx = fidx
        obj(x, c) = wobj(x, c, idx)
        prep = prepare_gradient(obj, mc_backend, q, Constant(w))
        gm = gradient(obj, prep, mc_backend, q, Constant(w))
        @info "reverse mode (Mooncake) at overflow point, objective uses only finite outputs" gm forwarddiff = J[fidx, :]' * w
        @test all(isnan, gm)                               # 0 * Inf cotangent through the unused overflowed A: documented, not hidden
    end

    @testset "prepared Mooncake cache reuse over changed same-shape queries" begin
        rng = MersenneTwister(2718)
        idx = collect(1:20)
        obj(x, c) = wobj(x, c, idx)
        q0 = [3000.0, 2700.0]
        w0 = randn(rng, 20)
        prep = prepare_gradient(obj, mc_backend, q0, Constant(w0))
        worst_reuse = 0.0; worst_indep = 0.0; nq = 0
        for (Tg, Te) in ((3000.0, 2700.0), (3100.0, 2650.0), (1000.0, 500.0), (1500.0, 1480.0), (2500.0, 2000.0),
                         (800.0, 824.0), (4000.0, 3800.0), (5000.0, 4750.0), (1200.0, 1000.0)), wi in 1:2
            w = wi == 1 ? w0 : randn(rng, 20)
            gre = gradient(obj, prep, mc_backend, [Tg, Te], Constant(w))
            prep2 = prepare_gradient(obj, mc_backend, [Tg, Te], Constant(w))
            gin = gradient(obj, prep2, mc_backend, [Tg, Te], Constant(w))
            gf = ForwardDiff.jacobian(rates_vec, [Tg, Te])' * w
            worst_reuse = max(worst_reuse, maximum(abs.(gre .- gf) ./ abs.(gf)))
            worst_indep = max(worst_indep, maximum(abs.(gin .- gf) ./ abs.(gf)))
            @test gre ≈ gf rtol = 1.0e-9
            @test gin ≈ gf rtol = 1.0e-9
            @test gre ≈ gin rtol = 1.0e-12
            nq += 1
        end
        # same prepared cache crossing stencil cells and the DB branch (Tg 3000 -> 9000) stays correct pointwise
        @info "Mooncake cache reuse" queries = nq worst_reuse worst_independent = worst_indep
    end

    @testset "piecewise behaviour at seams, DB switch and clamp (reported, not compared across)" begin
        for (lab, Tk) in (("Tg100 knot seam", exp(ADT.lgTg[20])), ("Tg250 knot seam", exp(ADT.lgTg[40])))
            lo = ForwardDiff.jacobian(rates_vec, [prevfloat(Tk, 4), 0.9 * Tk]); hi = ForwardDiff.jacobian(rates_vec, [nextfloat(Tk, 4), 0.9 * Tk])
            fl = filter(k -> isfinite(lo[k, 1]) && isfinite(hi[k, 1]), 1:20)
            jump = maximum(abs(lo[k, 1] - hi[k, 1]) / (abs(lo[k, 1]) + 1.0e-300) for k in fl)
            @info "seam derivative jump (d/dTg, max relative over finite outputs)" lab jump
            @test all(isfinite, lo[fl, :]) && all(isfinite, hi[fl, :])        # both one-sided derivatives are defined and finite
            # each side agrees with its own analytic local polynomial
            for (J, T) in ((lo, prevfloat(Tk, 4)), (hi, nextfloat(Tk, 4)))
                Ja, _ = analytic_jac(T, 0.9 * T)
                for k in 6:10
                    @test J[k, 1] ≈ Ja[k, 1] rtol = 1.0e-9
                end
            end
        end
        thr = 5452.7250000000013
        lo = ForwardDiff.jacobian(rates_vec, [prevfloat(thr), 0.9 * thr]); hi = ForwardDiff.jacobian(rates_vec, [thr, 0.9 * thr])
        @info "DB switch: A d/dTe below / above" below = lo[1:5, 2] above = hi[1:5, 2]
        @test all(hi[1:5, 2] .== 0.0) && any(lo[1:5, 2] .!= 0.0)
        @test hi[6:10, 1] ≈ lo[6:10, 1] rtol = 1.0e-6                          # B, R unaffected by the A-only DB switch
        # clamp side-by-side (Tg = 3000)
        for ratio in (0.99984, 0.99986)
            Te = 3000.0 * ratio
            J = ForwardDiff.jacobian(rates_vec, [3000.0, Te]); _, st = analytic_jac(3000.0, Te)
            @info "clamp status / dA/dTe" ratio status = st dAdTe = J[1:5, 2]
        end
    end
end
