# Chunk 11: opt-in fixed-geometry spline quadrature of the PDE correction integrals (src/HIPDEQuadPlan.jl) and its reverse rule
# (ext/CosmoRecChainRulesCoreExt.jl + reverse-only Mooncake bridge ext/CosmoRecMooncakeExt.jl).
#   K1-K6: kernel gates on synthetic (highly non-uniform) and native knot slices, bound γ_M = M u/(1 - M u), M = 3n + 2*255 + 32 (derived from the
#          operation count, not tuned); the reference for the fit adjoint is the IMPLEMENTED stored-factor map (Float64 factors for any eltype).
#   V2:    plan path vs default path: bitwise primal/levels/replay, exact ForwardDiff, full-grid stale-plan rejection.
#   V3/V4: extension loading, Mooncake dispatch of the bridged rule (is_primitive, test_rule, poisoned/zeroed-plan witness), zero cotangents,
#          SubArray fallback, preparation reuse, task concurrency, buffer reuse and accumulation, short-march dot-product test.
#   V5 (full-stage prepared gradients) lives in chunk11b_quadplan_stage_ad.jl (one large preparation at a time).
# The strict chunk27 MC-vs-ForwardDiff 1e-6 comparison is NOT part of this file (pre-existing failure, reported separately).
# REQUIRES COSMOREC_NATIVE_DATA_DIR.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using ChainRulesCore
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
using CosmoRec: HIQuadPlan, HIQuadPlans, _hi_quad_plan, _hi_quad_plans, _hi_spline_parts, _hi_fit_adjoint!, _hi_fit_factors, _hi_quad_pullback!,
    _check_quadplan, _patterson_entries, HIPDESetup, HIPDEAtom, NATIVE_HI_PDE_ATOM
@isdefined(SETUP7) || include("chunk7_helpers.jl")

const MC11 = AutoMooncake(; config = nothing)
γM11(n) = (M = 3n + 2 * 255 + 32; u = eps(Float64) / 2; M * u / (1 - M * u))
const OBS11 = Dict{String,Float64}()
rec11(k, v) = (OBS11[k] = max(get(OBS11, k, 0.0), v))

# nodes of (sub-interval i, level L) exactly as patterson_integrate places them
function nodes11(q::HIQuadPlan, i, L)
    a = q.a[i]; b = q.b[i]; xc = (a + b) / 2.0; Dx = (b - a) / 2.0
    xs = Float64[]
    for (u, _) in _patterson_entries(L)
        u == 0.0 ? push!(xs, xc) : (Del = Dx * u; push!(xs, xc + Del, xc - Del))
    end
    return xs
end
# the stored-factor fit adjoint evaluated on absolute values with all signs positive: componentwise magnitude scale of its rounding errors
function absadj11(ȳ, b̄, c̄, d̄, q)
    n = length(ȳ); N = n - 2; dx = q.dx; g = abs.(q.g); off = abs.(q.off); α = abs.(q.alpha); t = abs.(q.t)
    for j in 1:(n - 1)
        h = dx[j]; ȳ[j + 1] += b̄[j] / h; ȳ[j] += b̄[j] / h
        c̄[j + 1] += h * b̄[j] / 3.0 + d̄[j] / (3.0 * h); c̄[j] += 2.0 * h * b̄[j] / 3.0 + d̄[j] / (3.0 * h)
    end
    c̄[2] = c̄[2] / α[1]
    for i in 2:N; c̄[i + 1] = (c̄[i + 1] + off[i - 1] * c̄[i]) / α[i]; end
    for i in (N - 1):-1:1; c̄[i + 1] = c̄[i + 1] + t[i + 1] * c̄[i + 2]; end
    for i in 1:N
        G = 3.0 * c̄[i + 1]; ȳ[i + 2] += G * g[i + 1]; ȳ[i + 1] += G * (g[i + 1] + g[i]); ȳ[i] += G * g[i]
    end
    return ȳ
end
# synthetic grids: highly non-uniform (spacing ratios up to ~1e6), plus the minimal n = 3
function synth_knots11(rng, n)
    h = exp.(randn(rng, n - 1) .* 4.0)
    return cumsum(vcat(0.3, h ./ sum(h) .* 2.0))
end
function synth_plans11()
    rng = Random.Xoshiro(1101)
    out = HIQuadPlan[]
    for n in (3, 7, 40, 300)
        x = synth_knots11(rng, n)
        lo = x[1]; hi = x[end]; w = hi - lo
        edges = [lo, lo + 0.13w, lo + 0.13w, lo + 0.5w, lo + 0.31w, hi]      # includes an empty (a == b) and a reversed (a > b) sub-interval
        push!(out, _hi_quad_plan(x, edges, 1.0e-60))
    end
    return out
end

@testset "Chunk 11: fixed-geometry spline quadrature plan and reverse rule" begin
    rng = Random.Xoshiro(20261005)
    plans = synth_plans11()
    natplans = _hi_quad_plans(SETUP7)
    allplans = vcat(plans, [natplans.di1, natplans.two3, natplans.raman])

    @testset "K1: probe (lo, delta) reproduce spline_eval_native bitwise" begin
        for q in allplans, i in eachindex(q.a)
            q.a[i] >= q.b[i] && continue
            n = length(q.knots)
            s = NaturalCubicSpline{Float64}(q.knots, randn(rng, n), randn(rng, n - 1), randn(rng, n), randn(rng, n - 1))
            Slo = NaturalCubicSpline{Float64}(q.knots, Float64.(1:n), zeros(n - 1), zeros(n), zeros(n - 1))
            Sd = NaturalCubicSpline{Float64}(q.knots, zeros(n), ones(n - 1), zeros(n), zeros(n - 1))
            for L in (1, 4, 8), x in nodes11(q, i, L)
                lo = Int(spline_eval_native(Slo, x)); d = spline_eval_native(Sd, x)
                @test s.y[lo] + d * (s.b[lo] + d * (s.c[lo] + d * s.d[lo])) === spline_eval_native(s, x)
            end
        end
    end

    @testset "K6a: stored factors reproduce the constructor's c bitwise" begin
        for q in allplans
            n = length(q.knots); N = n - 2
            y = randn(rng, n); ref = natural_cubic_spline(q.knots, y)
            gam = [3.0 * ((y[i + 2] - y[i + 1]) * q.g[i + 1] - (y[i + 1] - y[i]) * q.g[i]) for i in 1:N]
            z = similar(gam); z[1] = gam[1]
            for i in 2:N; z[i] = gam[i] - q.t[i] * z[i - 1]; end
            c = zeros(n); c[N + 1] = z[N] / q.alpha[N]
            for i in (N - 1):-1:1; c[i + 1] = (z[i] - q.off[i] * c[i + 2]) / q.alpha[i]; end
            @test c == ref.c
        end
    end

    @testset "K6b/K4: fit adjoint dot test against the implemented forward fit" begin
        for q in allplans
            n = length(q.knots)
            for _ in 1:3
                ẏ = randn(rng, n); fit = natural_cubic_spline(q.knots, ẏ)          # linear in y: the JVP of the implemented map
                ȳ0 = randn(rng, n); b̄ = randn(rng, n - 1); c̄ = randn(rng, n); d̄ = randn(rng, n - 1)
                terms = vcat(ȳ0 .* ẏ, b̄ .* fit.b, c̄[2:(n - 1)] .* fit.c[2:(n - 1)], d̄ .* fit.d)
                rhs = sum(terms)
                ȳ = _hi_fit_adjoint!(copy(ȳ0), copy(b̄), copy(c̄), copy(d̄), q.dx, q.g, q.off, q.alpha, q.t)
                lhs = dot(ȳ, ẏ)
                @test abs(lhs - rhs) <= γM11(n) * (sum(abs, terms) + sum(abs, ȳ .* ẏ))
                rec11("K6b_dot_rel", abs(lhs - rhs) / (sum(abs, terms) + sum(abs, ȳ .* ẏ)))
            end
        end
    end

    @testset "K6c/K6d: adjoint vs 256-bit same-factor evaluation; symmetric-reuse gap (report)" begin
        for q in allplans
            n = length(q.knots); N = n - 2
            ȳ0 = randn(rng, n); b̄ = randn(rng, n - 1); c̄ = randn(rng, n); d̄ = randn(rng, n - 1)
            f64 = _hi_fit_adjoint!(copy(ȳ0), copy(b̄), copy(c̄), copy(d̄), q.dx, q.g, q.off, q.alpha, q.t)
            ref256 = setprecision(BigFloat, 256) do
                Float64.(_hi_fit_adjoint!(Base.big.(ȳ0), Base.big.(b̄), Base.big.(c̄), Base.big.(d̄), Base.big.(q.dx), Base.big.(q.g), Base.big.(q.off),
                                          Base.big.(q.alpha), Base.big.(q.t)))
            end
            scale = absadj11(abs.(ȳ0), abs.(b̄), abs.(c̄), abs.(d̄), q)          # same recurrences on absolute values (rounding-error magnitude)
            @test all(abs.(f64 .- ref256) .<= γM11(n) .* scale .+ floatmin())
            rec11("K6c_rel", maximum(abs.(f64 .- ref256) ./ (scale .+ floatmin())))
            # K6d (report only): forward solve re-applied to c̄ (ideal symmetric adjoint) instead of the stored-factor transpose
            w = copy(c̄); w[1] = 0.0; w[n] = 0.0
            for j in 1:(n - 1)
                h = q.dx[j]; w[j + 1] += -(h * b̄[j]) / 3.0 + d̄[j] / (3.0 * h); w[j] += -(2.0 * h * b̄[j]) / 3.0 - d̄[j] / (3.0 * h)
            end
            wv = w[2:(n - 1)]; zz = similar(wv); zz[1] = wv[1]
            for i in 2:N; zz[i] = wv[i] - q.t[i] * zz[i - 1]; end
            gs = similar(wv); gs[N] = zz[N] / q.alpha[N]
            for i in (N - 1):-1:1; gs[i] = (zz[i] - q.off[i] * gs[i + 1]) / q.alpha[i]; end
            ws = copy(c̄); ws[1] = 0.0; ws[n] = 0.0
            for j in 1:(n - 1)
                h = q.dx[j]; ws[j + 1] += -(h * b̄[j]) / 3.0 + d̄[j] / (3.0 * h); ws[j] += -(2.0 * h * b̄[j]) / 3.0 - d̄[j] / (3.0 * h)
            end
            ws2 = _hi_fit_adjoint!(zeros(n), zeros(n - 1), ws, zeros(n - 1), q.dx, q.g, q.off, q.alpha, q.t)   # transpose route on the same c̄
            ysym = zeros(n)
            for i in 1:N
                G = 3.0 * gs[i]; ysym[i + 2] += G * q.g[i + 1]; ysym[i + 1] -= G * (q.g[i + 1] + q.g[i]); ysym[i] += G * q.g[i]
            end
            rec11("K6d_symmetric_reuse_gap_rel", maximum(abs.(ysym .- ws2)) / max(maximum(abs, ws2), floatmin()))
        end
    end

    @testset "K2/K3/K5: ell vs primal parts, 256-bit Dual reference and Float64 ForwardDiff" begin
        for q in allplans
            K = length(q.a); n = length(q.knots)
            F = randn(rng, n) .* exp.(randn(rng, n))
            v = randn(rng, n)
            for L in 1:8
                forced = ntuple(_ -> L, K)
                parts, levs, oks = _hi_spline_parts(q, F, forced)
                for i in 1:K
                    if q.a[i] >= q.b[i]
                        @test levs[i] == 0 && parts[i] == 0.0 && oks[i]
                        @test all(iszero, q.ell[:, i, L])
                        continue
                    end
                    @test levs[i] == L
                    ℓ = q.ell[:, i, L]
                    @test abs(dot(ℓ, F) - parts[i]) <= γM11(n) * sum(abs, ℓ .* F)                                   # K2
                    rec11("K2_rel", abs(dot(ℓ, F) - parts[i]) / sum(abs, ℓ .* F))
                    gfd = ForwardDiff.gradient(G -> _hi_spline_parts(q, G, forced)[1][i], F)
                    @test abs(dot(ℓ .- gfd, v)) <= γM11(n) * sum(abs, gfd .* v)                                          # K5
                    gbig = setprecision(BigFloat, 256) do
                        Float64.(ForwardDiff.gradient(G -> _hi_spline_parts(q, G, forced)[1][i], big.(F)))
                    end
                    @test abs(dot(ℓ .- gbig, v)) <= γM11(n) * sum(abs, gbig .* v)                                        # K3
                    rec11("K3_rel", abs(dot(ℓ .- gbig, v)) / sum(abs, gbig .* v))
                end
            end
            # adaptive call: returned levels are the realised ones and the pullback uses them
            parts, levs, _ = _hi_spline_parts(q, F, ntuple(_ -> 0, K))
            @test parts == _hi_spline_parts(q, F, ntuple(i -> levs[i], K))[1]
            dF = zeros(n); dp = ntuple(_ -> randn(rng), K)
            _hi_quad_pullback!(dF, q, levs, dp)
            @test dF ≈ sum(levs[i] == 0 ? zeros(n) : dp[i] .* q.ell[:, i, levs[i]] for i in 1:K) rtol = 0 atol = 0
            dz = zeros(n); _hi_quad_pullback!(dz, q, levs, ntuple(_ -> 0.0, K)); @test all(iszero, dz)          # zero cotangent
            # accumulation, never overwrite: exact against the documented order (sub-interval i = 1..K, then j) starting from the given vector
            for acc0 in (ones(n), randn(rng, n))
                ref = copy(acc0)
                for i in 1:K
                    (levs[i] == 0 || iszero(dp[i])) && continue
                    for j in 1:n; ref[j] += dp[i] * q.ell[j, i, levs[i]]; end
                end
                acc = copy(acc0); _hi_quad_pullback!(acc, q, levs, dp)
                @test acc == ref
            end
        end
    end

    @testset "Contract: invalid levels and plan inputs are rejected (no unchecked table access)" begin
        good = synth_knots11(Random.Xoshiro(7), 12); lo_ = good[1]; hi_ = good[end]
        @test _hi_quad_plan(good, [lo_, hi_], 1.0e-60) isa HIQuadPlan{1}
        for (kn, ed, ep) in ((good[1:2], [lo_, good[2]], 1.0e-60),                      # too few knots
                             (vcat(good[1:5], NaN, good[7:end]), [lo_, hi_], 1.0e-60),  # non-finite knot
                             (vcat(good[1:5], good[5], good[7:end]), [lo_, hi_], 1.0e-60),  # non-increasing knots
                             (good, [lo_], 1.0e-60),                                     # no sub-interval
                             (good, [lo_, NaN], 1.0e-60),                                # non-finite edge
                             (good, [lo_ - 1.0, hi_], 1.0e-60),                          # non-empty sub-interval outside the knots
                             (good, [lo_, hi_], -1.0))                                   # invalid epsabs0
            @test_throws ArgumentError _hi_quad_plan(kn, ed, ep)
        end
        for q in (natplans.di1, natplans.two3)
            K = length(q.a); n = length(q.knots); F = randn(rng, n)
            for bad in (-1, 9)
                forced = ntuple(i -> i == 1 ? bad : 0, K)
                @test_throws ArgumentError _hi_spline_parts(q, F, forced)
                @test_throws ArgumentError ChainRulesCore.rrule(_hi_spline_parts, q, F, forced)
                dF = zeros(n)
                @test_throws ArgumentError _hi_quad_pullback!(dF, q, ntuple(i -> i == 1 ? bad : 1, K), ntuple(_ -> 1.0, K))
                @test all(iszero, dF)                                                    # nothing accumulated before rejection
            end
            @test_throws ArgumentError CosmoRec._hi_forced_levels!(PattersonLevels(:replay, [9, 9, 9, 9, 9], 0), q)
            @test_throws ArgumentError CosmoRec._hi_forced_levels!(PattersonLevels(:replay, [0, 0, 0, 0, 0], 0), q)
        end
        # malformed replay data reaching the plan path through the public stage keyword is rejected, not clamped
        m = HIPDEModel(SETUP7, HIPopulationSplines(NODES5[:, 1:9]; zs = ZS6, ze = ZE6),
                       hi_pde_coefficients(NODES5[:, 1:9], HIPopulationSplines(NODES5[:, 1:9]; zs = ZS6, ze = ZE6), ACC5, HTAB6, LNBITOT6,
                                           NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6), ACC5)
        @test_throws ArgumentError hi_pde_corrections(m; zs = 2500.0, ze = 2470.0, levels = PattersonLevels(:replay, fill(9, 64), 0), quadplan = natplans)
    end

    # ---------------------------------------------------------------- V2: native parity and stale plans
    D11 = HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6)
    R11 = NODES5[:, 1:9]
    stagef(o) = (o.z, o.DI1_2s, o.DF_2g, o.DF_R, o.scale, o.unconverged, o.y)
    function out11(p; plan = nothing, zs = nothing, ze = nothing)
        o, _ = hi_diffusion_stage(RM5, D11, R11; hscale = p[1], nbscale = p[2], quadplan = plan)
        mask = (o.z .< 2000.0) .& (o.z .> 500.0)
        vcat(o.DI1_2s[mask], o.DF_2g[1][mask], o.DF_2g[2][mask], o.DF_R[1][mask])
    end
    @testset "V2: plan path is bitwise the default path (stage, levels, replay, ForwardDiff)" begin
        for p in ([1.0, 1.0], [0.999, 1.001])
            r0 = PattersonLevels(); r1 = PattersonLevels()
            o0, f0 = hi_diffusion_stage(RM5, D11, R11; hscale = p[1], nbscale = p[2], levels = r0)
            o1, f1 = hi_diffusion_stage(RM5, D11, R11; hscale = p[1], nbscale = p[2], levels = r1, quadplan = natplans)
            @test stagef(o0) == stagef(o1)
            @test r0.levels == r1.levels
            @test all(o0.y .=== o1.y) && all(o0.DI1_2s .=== o1.DI1_2s)
            o2, _ = hi_diffusion_stage(RM5, D11, R11; hscale = p[1], nbscale = p[2], levels = replay(r0), quadplan = natplans)
            o3, _ = hi_diffusion_stage(RM5, D11, R11; hscale = p[1], nbscale = p[2], levels = replay(r0))
            # plan replay == default replay in EVERY field (incl. convergence metadata)
            @test stagef(o2) == stagef(o3)
            # property of the ORIGINAL default path (HIPDEIntegrals.jl `_hi_calc_DF`, replay branch returns `(r, true)`): replaying the recorded
            # levels reproduces every numeric output and scale of the adaptive run, while the `unconverged` metadata of replay is all zero by design
            @test (o3.z, o3.DI1_2s, o3.DF_2g, o3.DF_R, o3.scale, o3.y) == (o0.z, o0.DI1_2s, o0.DF_2g, o0.DF_R, o0.scale, o0.y)
            @test all(u -> all(iszero, u), o3.unconverged)
            @test (o2.z, o2.DI1_2s, o2.DF_2g, o2.DF_R, o2.scale, o2.y) == (o0.z, o0.DI1_2s, o0.DF_2g, o0.DF_R, o0.scale, o0.y)
            @test ForwardDiff.jacobian(q -> out11(q), p) == ForwardDiff.jacobian(q -> out11(q; plan = natplans), p)
        end
        # BigFloat (generic path, no rule) on the short production march with replayed levels
        m11(p) = HIPDEModel(SETUP7, HIPopulationSplines(R11; zs = ZS6, ze = ZE6),
                            hi_pde_coefficients(R11, HIPopulationSplines(R11; zs = ZS6, ze = ZE6), CosmoRec.ScaledBackground(ACC5, p[1], p[2]),
                                                HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6), CosmoRec.ScaledBackground(ACC5, p[1], p[2]))
        r = PattersonLevels(); hi_pde_corrections(m11([1.0, 1.0]); zs = 2500.0, ze = 2470.0, levels = r)
        pb = setprecision(BigFloat, 256) do; big.([1.0, 1.0]); end
        ob0, ob1 = setprecision(BigFloat, 256) do
            (hi_pde_corrections(m11(pb); zs = 2500.0, ze = 2470.0, T = BigFloat, levels = replay(r)),
             hi_pde_corrections(m11(pb); zs = 2500.0, ze = 2470.0, T = BigFloat, levels = replay(r), quadplan = natplans))
        end
        @test stagef(ob0) == stagef(ob1)
    end

    @testset "V2: full-grid stale-plan rejection" begin
        s = SETUP7; atom = NATIVE_HI_PDE_ATOM
        remake(x; ie = s.index_emission) = HIPDESetup(x, s.resonances, s.nresmax, s.n2g, s.nR, s.index_2, ie, s.ratio_R_ns, s.ratio_R_nd,
                                                       s.ratio_2g_ns, s.ratio_2g_nd, s.A_npns, s.A_npnd, s.prof)
        @test _check_quadplan(natplans, s, atom) === nothing
        @test _check_quadplan(natplans, remake(copy(s.x)), atom) === nothing                       # value-equal distinct copy
        ie = s.index_emission
        for j in (ie[1] ÷ 2, (ie[1] + ie[2]) ÷ 2, (ie[2] + length(s.x)) ÷ 2)                        # DI1 slice, 3s/3d-only region, Raman-only region
            x = copy(s.x); x[j] = nextfloat(x[j])
            s2 = remake(x)
            @test s2.x[1] == s.x[1] && s2.x[end] == s.x[end] && length(s2.x) == length(s.x)
            @test CosmoRec._hi_block_edges(s2, atom) == CosmoRec._hi_block_edges(s, atom)        # endpoints/edges unchanged
            @test_throws ArgumentError _check_quadplan(natplans, s2, atom)
            @test _check_quadplan(_hi_quad_plans(s2, atom), s2, atom) === nothing                  # explicit rebuild passes
        end
        x = copy(s.x); x[end] = nextfloat(x[end]); @test_throws ArgumentError _check_quadplan(natplans, remake(x), atom)
        @test_throws ArgumentError _check_quadplan(natplans, remake(copy(s.x); ie = [ie[1] + 1, ie[2]]), atom)
        # an atom change that does NOT move any edge is not stale geometry; one that moves an edge must be rejected
        mkatom(d3) = HIPDEAtom(atom.lyn, atom.lyn_A21, atom.nu21, [atom.Dnu_1s[1], atom.Dnu_1s[2], d3], atom.A_3p2s, atom.nu_3p2s,
                               atom.A_3s2p, atom.nu_3s2p, atom.A_3d2p, atom.nu_3d2p, atom.lyn_lambda21)
        atom_same = mkatom(nextfloat(atom.Dnu_1s[3]))
        if CosmoRec._hi_block_edges(s, atom_same) == CosmoRec._hi_block_edges(s, atom)
            @test _check_quadplan(natplans, s, atom_same) === nothing                              # identical geometry is accepted
        end
        atom2 = mkatom(atom.Dnu_1s[3] * (1.0 + 1.0e-6))
        @test CosmoRec._hi_block_edges(s, atom2).two3 != CosmoRec._hi_block_edges(s, atom).two3   # precondition: an actual edge changed
        @test_throws ArgumentError _check_quadplan(natplans, s, atom2)
        # the stage rejects before marching
        x = copy(s.x); x[ie[1] ÷ 2] = nextfloat(x[ie[1] ÷ 2])
        m = HIPDEModel(remake(x), HIPopulationSplines(R11; zs = ZS6, ze = ZE6),
                       hi_pde_coefficients(R11, HIPopulationSplines(R11; zs = ZS6, ze = ZE6), ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6), ACC5)
        @test_throws ArgumentError hi_pde_corrections(m; zs = 2500.0, ze = 2470.0, quadplan = natplans)
    end

    # ---------------------------------------------------------------- V3/V4: bridge
    @testset "V3: extensions, tangent contracts and Mooncake dispatch" begin
        ext = Base.get_extension(CosmoRec, :CosmoRecMooncakeExt); crx = Base.get_extension(CosmoRec, :CosmoRecChainRulesCoreExt)
        @test ext isa Module && crx isa Module
        @test Mooncake.tangent_type(typeof(natplans)) == Mooncake.NoTangent
        @test Mooncake.tangent_type(typeof(natplans.two3)) == Mooncake.NoTangent
        w = Base.get_world_counter()
        for q in (natplans.di1, natplans.two3, natplans.raman)
            K = length(q.a)
            @test Mooncake.is_primitive(Mooncake.DefaultCtx, Mooncake.ReverseMode, Tuple{typeof(_hi_spline_parts), typeof(q), Vector{Float64}, NTuple{K,Int}}, w)
            @test !Mooncake.is_primitive(Mooncake.DefaultCtx, Mooncake.ReverseMode, Tuple{typeof(_hi_spline_parts), typeof(q), SubArray{Float64,1,Vector{Float64},Tuple{UnitRange{Int}},true}, NTuple{K,Int}}, w)
        end
        @test Mooncake.is_primitive(Mooncake.DefaultCtx, Mooncake.ReverseMode, Tuple{typeof(_check_quadplan), typeof(natplans), HIPDESetup, HIPDEAtom}, w)
        # direct ChainRules pullback: zero/thunked cotangents, never ZeroTangent/NoTangent for F
        q = natplans.two3; K = length(q.a); n = length(q.knots); F = randn(rng, n)
        forced = ntuple(_ -> 6, K)
        out, pb = ChainRulesCore.rrule(_hi_spline_parts, q, F, forced)
        @test out == _hi_spline_parts(q, F, forced)
        for Δ in (ZeroTangent(), Tangent{Any}(ZeroTangent(), NoTangent(), NoTangent()), Tangent{Any}(Tangent{Any}(ntuple(_ -> 0.0, K)...), NoTangent(), NoTangent()))
            r = pb(Δ); @test r[1] isa NoTangent && r[2] isa NoTangent && r[4] isa NoTangent && r[3] isa Vector{Float64} && all(iszero, r[3])
        end
        dp = ntuple(_ -> randn(rng), K)
        r = pb(Tangent{Any}(Tangent{Any}(dp...), NoTangent(), NoTangent()))
        @test r[3] == _hi_quad_pullback!(zeros(n), q, out[2], dp)
        r2 = pb(ChainRulesCore.@thunk(Tangent{Any}(Tangent{Any}(dp...), NoTangent(), NoTangent()))); @test r2[3] == r[3]
        # Mooncake's own rule test (frozen levels so finite differences never cross a stopping branch)
        for p in (natplans.di1, natplans.two3, natplans.raman)
            Kp = length(p.a)
            Mooncake.TestUtils.test_rule(Random.Xoshiro(11), _hi_spline_parts, p, randn(Random.Xoshiro(12), length(p.knots)), ntuple(_ -> 5, Kp);
                                         is_primitive = true, mode = Mooncake.ReverseMode, print_results = false)
        end
    end

    # short production march + integrals (chunk 8a scale) with frozen dynamic inputs = the populations rows
    X011 = vec(R11[:, 2:9])
    function f8(x; plan = natplans)
        rows = hcat(R11[:, 1], reshape(x, size(R11, 1), 8)); pops = HIPopulationSplines(rows; zs = ZS6, ze = ZE6)
        o = hi_pde_corrections(HIPDEModel(SETUP7, pops, hi_pde_coefficients(rows, pops, ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6), ACC5);
                               zs = 2500.0, ze = 2470.0, T = eltype(x), quadplan = plan)
        return vcat(o.DI1_2s, o.DF_2g[1], o.DF_2g[2], o.DF_R[1])
    end
    @testset "V3: dispatch witness (poisoned / zeroed plan), dot-product test, preparation reuse, tasks" begin
        nout = length(f8(X011)); w = randn(rng, nout)
        poison(q, s) = HIQuadPlan{length(q.a)}(q.knots, q.a, q.b, q.epsabs0, q.dx, q.g, q.off, q.alpha, q.t, s .* q.ell)
        bad = HIQuadPlans(natplans.di1, poison(natplans.two3, 2.0), natplans.raman, natplans.xfull, natplans.index_emission, natplans.nresmax)
        zer = HIQuadPlans(poison(natplans.di1, 0.0), poison(natplans.two3, 0.0), poison(natplans.raman, 0.0), natplans.xfull, natplans.index_emission, natplans.nresmax)
        obj(x, w, plan) = dot(w, f8(x; plan = plan))
        pg = prepare_gradient(obj, MC11, X011, Constant(w), Constant(natplans))
        pg2 = prepare_gradient(obj, MC11, X011, Constant(w), Constant(natplans))
        pd = prepare_gradient(obj, MC11, X011, Constant(w), Constant(nothing))
        g = gradient(obj, pg, MC11, X011, Constant(w), Constant(natplans))
        gd = gradient(obj, pd, MC11, X011, Constant(w), Constant(nothing))
        gb = gradient(obj, pg, MC11, X011, Constant(w), Constant(bad))
        gz = gradient(obj, pg, MC11, X011, Constant(w), Constant(zer))
        @test f8(X011; plan = bad) == f8(X011; plan = zer) == f8(X011)                                     # ell never enters the primal
        @test ForwardDiff.derivative(t -> f8(X011 .+ t; plan = bad), 0.0) == ForwardDiff.derivative(t -> f8(X011 .+ t), 0.0)
        @test gb != g && gz != g                                                                             # the rule's tables are what Mooncake uses
        @test maximum(abs, g .- gd) <= 1e-10 * maximum(abs, gd)                                             # rule vs derived (G-reg bound)
        rec11("V3_8a_rule_vs_derived_rel", maximum(abs, g .- gd) / maximum(abs, gd))
        for k in 1:2
            x = k == 1 ? X011 : X011 .* (1 .+ 1.0e-9 .* randn(rng, length(X011)))
            gA = gradient(obj, pg, MC11, x, Constant(w), Constant(natplans)); gB = gradient(obj, pg2, MC11, x, Constant(w), Constant(natplans))
            @test gA == gB
            v = randn(rng, length(x)) .* abs.(x)
            jv = ForwardDiff.derivative(t -> f8(x .+ t .* v), 0.0)
            e = abs(dot(gA, v) - dot(w, jv)) / max(sum(abs.(gA .* v)), sum(abs.(w .* jv)))
            @test e < 1e-10
            rec11("VJP_8a_plan", e)
        end
        # independent preparations on concurrent tasks agree bitwise with the sequential results
        xs = [X011 .* (1 .+ 1.0e-9 .* randn(Random.Xoshiro(30 + k), length(X011))) for k in 1:2]
        seq = [gradient(obj, pg, MC11, x, Constant(w), Constant(natplans)) for x in xs]
        par = fetch.([Threads.@spawn gradient(obj, prepare_gradient(obj, MC11, xs[k], Constant(w), Constant(natplans)), MC11, xs[k], Constant(w), Constant(natplans)) for k in 1:2])
        @test par == seq
    end

    @testset "V4: buffer reuse, accumulation and SubArray fallback" begin
        q = natplans.two3; K = length(q.a); n = length(q.knots)
        F0 = randn(Random.Xoshiro(41), n); F1 = randn(Random.Xoshiro(42), n); forced = ntuple(_ -> 6, K)
        c = randn(Random.Xoshiro(43), K)
        function reuse(p)                     # one buffer refilled between two primitive calls, plus a second use of the buffer
            B = similar(p); B .= p .* F0
            s1 = dot(c, collect(_hi_spline_parts(q, B, forced)[1]))
            B .= p .* F1
            s2 = dot(c, collect(_hi_spline_parts(q, B, forced)[1]))
            return s1 + 2s2 + sum(B)
        end
        function fresh(p)
            B0 = p .* F0; B1 = p .* F1
            return dot(c, collect(_hi_spline_parts(q, B0, forced)[1])) + 2dot(c, collect(_hi_spline_parts(q, B1, forced)[1])) + sum(B1)
        end
        p0 = 1.0 .+ 0.1 .* randn(Random.Xoshiro(44), n)
        gr = gradient(reuse, prepare_gradient(reuse, MC11, p0), MC11, p0)
        gf = gradient(fresh, prepare_gradient(fresh, MC11, p0), MC11, p0)
        gfd = ForwardDiff.gradient(fresh, p0)
        @test gr == gf
        @test maximum(abs, gr .- gfd) <= γM11(n) * 100 * maximum(abs, gfd)
        rec11("V4_reuse_vs_FD_rel", maximum(abs, gr .- gfd) / maximum(abs, gfd))
        sv(p) = dot(c, collect(_hi_spline_parts(q, view(p .* F0, 1:n), forced)[1]))
        gs = gradient(sv, prepare_gradient(sv, MC11, p0), MC11, p0)
        @test maximum(abs, gs .- ForwardDiff.gradient(sv, p0)) <= γM11(n) * 100 * maximum(abs, gs)
    end

    @testset "info" begin
        foreach(kv -> println("Chunk11 ", kv[1], " = ", kv[2]), sort!(collect(OBS11)))
    end
end
