# Chunk 11e: opt-in `reduced_profiles` production of the closed PDE stage (skips the `hp_phi` entries no code reads) and the unchanged FULL public state.
# Gates: reduced private state == full state bitwise on every non-hp field and on the consumed hp span (Float64, Dual z, Dual scales, nested Dual, hscale/nbscale != 1, z across 2000, 1500, 2500),
# poison (NaN) exactly outside the span; public `hi_pde_march` observer sees the COMPLETE state at both endpoints; closed stage reduced == default bitwise (values, scales, flags, levels, replay,
# quadplan / fitwork combinations, ForwardDiff Jacobian, BigFloat); public stage/pipeline keyword; argument validation; prepared DI-Mooncake gradient vs default and vs ForwardDiff
# (repeated use on one tape, nonzero initial cotangent, all parameters); overlay validation order of the scratch method. REQUIRES COSMOREC_NATIVE_DATA_DIR (see test/chunk7_helpers.jl).
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
using CosmoRec
using CosmoRec: HIQuadWork, HIFitScratch, _hi_quad_plans, _hi_quad_plan, _hi_pde_rhs_coefficients!, _hi_pde_march, _hi_hp_span, _hi_poison_hp!, _hi_spline_parts
@isdefined(SETUP7) || include("chunk7_helpers.jl")

const MC11E = AutoMooncake(; config = nothing)
bits11e(x::Float64) = reinterpret(UInt64, x)
bits11e(x::ForwardDiff.Dual) = (bits11e(ForwardDiff.value(x)), map(bits11e, ForwardDiff.partials(x).values))
bits11e(x::AbstractArray) = map(bits11e, x)
bits11e(x::Tuple) = map(bits11e, x)
bits11e(x::NamedTuple) = map(bits11e, values(x))
bits11e(x::Union{Integer,Bool}) = x
bits11e(x) = isstructtype(typeof(x)) ? ntuple(i -> bits11e(getfield(x, i)), fieldcount(typeof(x))) : x
fields11e(st) = Dict("A" => st.A, "B" => st.B, "C" => st.C, "D" => st.D, "exp_x" => st.exp_x, "phi1" => st.phi[1], "phi2" => st.phi[2], "aV" => st.aV, "Dnu" => st.Dnu, "pd" => st.pd, "pd_eff" => st.pd_eff,
                     "Dnem_eff" => st.Dnem_eff, "P2g_ns1" => st.P2g_ns[1], "P2g_ns2" => st.P2g_ns[2], "P2g_nd1" => st.P2g_nd[1], "P2g_nd2" => st.P2g_nd[2], "PR_ns1" => st.PR_ns[1], "PR_ns2" => st.PR_ns[2],
                     "PR_nd1" => st.PR_nd[1], "PR_nd2" => st.PR_nd[2])
isnan11e(v) = v isa ForwardDiff.Dual ? isnan(ForwardDiff.value(v)) : isnan(v)

@testset "Chunk 11e: reduced hp_phi production of the closed PDE stage" begin
    R11 = NODES5[:, 1:9]; D11 = HIDiffusionInputs(SETUP7, HTAB6, LNBITOT6); nat = _hi_quad_plans(SETUP7)
    function model(hs, nb)
        bg = CosmoRec.ScaledBackground(RM5.cosmos, hs, nb); pops = HIPopulationSplines(R11; zs = ZS6, ze = ZE6)
        return HIPDEModel(SETUP7, pops, hi_pde_coefficients(R11, pops, bg, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6), bg)
    end
    span = _hi_hp_span(SETUP7); npts = length(SETUP7.x)
    @testset "span is the consumed set" begin
        @test span == (SETUP7.index_2 + 1):max(SETUP7.index_emission[2], SETUP7.index_2 + 1) && first(span) >= 1 && last(span) <= npts && length(span) < npts
        @test_throws ArgumentError _hi_pde_march(model(1.0, 1.0), (0:5); zs = 2500.0, ze = 2490.0)
        @test_throws ArgumentError _hi_pde_march(model(1.0, 1.0), (1:(npts + 1)); zs = 2500.0, ze = 2490.0)
    end
    function check_state(m, z)
        full = hi_pde_rhs_coefficients(m, z); red = typeof(full)(npts); _hi_poison_hp!(red, span)
        _hi_pde_rhs_coefficients!(red, m, z, span)
        ff = fields11e(full); fr = fields11e(red)
        @test all(bits11e(ff[k]) == bits11e(fr[k]) for k in keys(ff))                                  # every non-hp field bit for bit (incl. Dual partials)
        @test bits11e(full.hp_phi[1][span]) == bits11e(red.hp_phi[1][span])                                 # the consumed span
        @test all(isnan11e, red.hp_phi[2]) && all(isnan11e(red.hp_phi[1][i]) for i in eachindex(red.hp_phi[1]) if !(i in span))
        @test !any(isnan11e, full.hp_phi[1]) && !any(isnan11e, full.hp_phi[2])                              # the public full state has no poison
    end
    @testset "reduced state vs full state, Float64, hscale/nbscale != 1" begin
        for (hs, nb) in ((1.0, 1.0), (1.0001, 0.9999)), z in (1500.0, 1999.9, 2000.1, 2500.0)
            check_state(model(hs, nb), z)
        end
    end
    @testset "Dual z / Dual scales (Te, Tm product partials) / nested Dual" begin
        for z in (1500.0, 2500.0)
            check_state(model(1.0, 1.0), ForwardDiff.Dual(z, 1.0))
            check_state(model(ForwardDiff.Dual(1.0001, 1.0, 0.0), ForwardDiff.Dual(0.9999, 0.0, 1.0)), z)
            Tg = ForwardDiff.Tag(:o, Float64); Ti = ForwardDiff.Tag(:i, ForwardDiff.Dual{typeof(Tg),Float64,1})
            check_state(model(1.0, 1.0), ForwardDiff.Dual{typeof(Ti)}(ForwardDiff.Dual{typeof(Tg)}(z, 1.0), ForwardDiff.Dual{typeof(Tg)}(1.0, 0.0)))
        end
    end
    @testset "public observer sees the COMPLETE state (default march), reduced march is closed-stage only" begin
        m = model(1.0001, 0.9999); seen = Ref(0)
        y1, z1 = hi_pde_march(m; zs = 2500.0, ze = 2470.0, on_step = (k, zin, zout, y, S) -> (seen[] += 1; for stt in (S.prev, S.cur); @test !any(isnan, stt.hp_phi[1]) && !any(isnan, stt.hp_phi[2]); end))
        @test seen[] == 3
        ref = Any[]; hi_pde_march(m; zs = 2500.0, ze = 2470.0, on_step = (k, zin, zout, y, S) -> push!(ref, (copy(S.cur.hp_phi[1]), copy(S.cur.hp_phi[2]), copy(S.prev.hp_phi[1]))))
        pri = Any[]; _hi_pde_march(m, nothing; zs = 2500.0, ze = 2470.0, on_step = (k, zin, zout, y, S) -> push!(pri, (copy(S.cur.hp_phi[1]), copy(S.cur.hp_phi[2]), copy(S.prev.hp_phi[1]))))
        @test bits11e(ref) == bits11e(pri)
        @test_throws MethodError hi_pde_march(m; zs = 2500.0, ze = 2470.0, reduced_profiles = true)       # not exposed on the public march
        y2, z2 = _hi_pde_march(m, span; zs = 2500.0, ze = 2470.0); @test bits11e(y1) == bits11e(y2) && z1 == z2
    end
    stagef(o) = (o.z, o.DI1_2s, o.DF_2g, o.DF_R, o.scale, o.unconverged, o.y)
    @testset "closed stage: reduced == default bitwise (levels, replay, quadplan, fitwork)" begin
        m = model(1.0001, 0.9999)
        for kw in ((;), (; quadplan = nat), (; quadplan = nat, fitwork = true))
            r0 = PattersonLevels(); r1 = PattersonLevels()
            o0 = hi_pde_corrections(m; zs = 2500.0, ze = 2400.0, levels = r0, kw...); o1 = hi_pde_corrections(m; zs = 2500.0, ze = 2400.0, levels = r1, reduced_profiles = true, kw...)
            @test bits11e(stagef(o0)) == bits11e(stagef(o1)) && r0.levels == r1.levels
            @test bits11e(stagef(hi_pde_corrections(m; zs = 2500.0, ze = 2400.0, levels = replay(r0), reduced_profiles = true, kw...))) == bits11e(stagef(hi_pde_corrections(m; zs = 2500.0, ze = 2400.0, levels = replay(r0), kw...)))
        end
        for bad in (1, "x", nothing, 1.0)
            @test_throws ArgumentError hi_pde_corrections(m; zs = 2500.0, ze = 2490.0, reduced_profiles = bad)
            @test_throws ArgumentError hi_diffusion_stage(RM5, D11, R11; reduced_profiles = bad)
        end
        @test_throws ArgumentError hi_pde_corrections(m; zs = 2500.0, ze = 2490.0, reduced_profiles = true, on_step = (a...) -> nothing)        # no user observer can see the incomplete state
        # ForwardDiff Jacobian of the short stage w.r.t. (hscale, nbscale): bitwise reduced vs default
        out(p; kw...) = (o = hi_pde_corrections(model(p[1], p[2]); zs = 2500.0, ze = 2470.0, T = eltype(p), kw...); vcat(o.DI1_2s, o.DF_2g[1], o.DF_2g[2], o.DF_R[1], o.y))
        @test bits11e(ForwardDiff.jacobian(out, [1.0001, 0.9999])) == bits11e(ForwardDiff.jacobian(p -> out(p; reduced_profiles = true), [1.0001, 0.9999]))
        # BigFloat stage
        r = PattersonLevels(); hi_pde_corrections(m; zs = 2500.0, ze = 2470.0, levels = r)
        setprecision(BigFloat, 256) do
            mb = model(big(1.0), big(1.0))
            ob0 = hi_pde_corrections(mb; zs = 2500.0, ze = 2470.0, T = BigFloat, levels = replay(r)); ob1 = hi_pde_corrections(mb; zs = 2500.0, ze = 2470.0, T = BigFloat, levels = replay(r), reduced_profiles = true)
            @test stagef(ob0) == stagef(ob1)
        end
    end
    @testset "public stage keyword (full diffusion stage, hscale/nbscale != 1)" begin
        for kw in ((;), (; quadplan = nat, fitwork = true))
            o0, f0 = hi_diffusion_stage(RM5, D11, R11; hscale = 1.0001, nbscale = 0.9999, kw...)
            o1, f1 = hi_diffusion_stage(RM5, D11, R11; hscale = 1.0001, nbscale = 0.9999, reduced_profiles = true, kw...)
            @test bits11e(stagef(o0)) == bits11e(stagef(o1)) && bits11e(f0) == bits11e(f1)
        end
    end
    @testset "reverse: prepared DI-Mooncake of the closed stage, reduced vs default and vs ForwardDiff" begin
        X0 = vec(R11[:, 2:9]); n = length(X0)
        function stage(p, red)                      # p = [hscale, nbscale, rows...]
            rows = hcat(R11[:, 1], reshape(p[3:end], size(R11, 1), 8)); bg = CosmoRec.ScaledBackground(ACC5, p[1], p[2]); pops = HIPopulationSplines(rows; zs = ZS6, ze = ZE6)
            m = HIPDEModel(SETUP7, pops, hi_pde_coefficients(rows, pops, bg, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6), bg)
            o = hi_pde_corrections(m; zs = 2500.0, ze = 2470.0, T = eltype(p), reduced_profiles = red)
            return vcat(o.DI1_2s, o.DF_2g[1], o.DF_2g[2], o.DF_R[1])
        end
        p0 = vcat([1.0001, 0.9999], X0); w = randn(Random.Xoshiro(81), length(stage(p0, false)))
        obj(p, red) = dot(w, stage(p, red)) + sum(abs2, p[1:2]) * 3.0                      # nonzero initial cotangent on hscale/nbscale
        twice(p, red) = obj(p, red) + 2 * dot(w, stage(p .* (1 .+ 1.0e-3 .* sin.(1:length(p))), red))   # two uses of the stage on one tape
        gd = gradient(obj, prepare_gradient(obj, MC11E, p0, Constant(false)), MC11E, p0, Constant(false))
        gr = gradient(obj, prepare_gradient(obj, MC11E, p0, Constant(true)), MC11E, p0, Constant(true))
        rowsc(g) = g[3:end]; relgap(a, b) = maximum(abs, rowsc(a) .- rowsc(b)) / maximum(abs, rowsc(b))
        println("Chunk11e reverse reduced-vs-default: rows rel gap = ", relgap(gr, gd), ", hs/nb gap = ", maximum(abs, gr[1:2] .- gd[1:2]), ", #bitwise-different entries = ", count(gr .!= gd), " of ", length(gd))
        @test all(isfinite, gr) && relgap(gr, gd) <= 1.0e-10 && gr[1:2] == gd[1:2]               # all parameters (hscale, nbscale, rows): same gradient (reverse accumulation order can differ in the last bits)
        @test !iszero(gr[1]) && !iszero(gr[2])
        v = randn(Random.Xoshiro(82), length(p0)) .* abs.(p0); jv = ForwardDiff.derivative(t -> obj(p0 .+ t .* v, true), 0.0)
        @test abs(dot(gr, v) - jv) <= 1.0e-5 * max(sum(abs.(gr .* v)), abs(jv))              # Float64 stage derivative is ill-conditioned at ~1e-6..1e-5 (known, see chunk11b notes)
        g2d = gradient(twice, prepare_gradient(twice, MC11E, p0, Constant(false)), MC11E, p0, Constant(false))
        g2r = gradient(twice, prepare_gradient(twice, MC11E, p0, Constant(true)), MC11E, p0, Constant(true))
        @test all(isfinite, g2r) && relgap(g2r, g2d) <= 1.0e-10 && maximum(abs, g2r[1:2] .- g2d[1:2]) <= 1.0e-12 * maximum(abs, g2d[1:2])
    end
    @testset "scratch overlay validates BEFORE any mutation (levels, length), exactly like the fast method" begin
        x = [0.0, 0.1, 0.13, 0.4, 0.49, 1.0, 1.1]; q = _hi_quad_plan(x, [x[1], x[end]], 1.0e-60); F = collect(1.0:7.0)
        w = HIFitScratch{Float64}(7); foreach(a -> fill!(a, 42.0), (w.y, w.c, w.z, w.b, w.d))
        f(p, forced) = _hi_spline_parts(q, p, forced, w)[1][1]
        @test_throws ArgumentError f(F, (9,)); @test_throws DimensionMismatch f(F[1:6], (6,))        # fast (non-AD) method
        @test all(a -> all(==(42.0), a), (w.y, w.c, w.z, w.b, w.d))
        h(p) = f(p, (9,)); hl(p) = f(p[1:6], (6,))
        @test_throws ArgumentError gradient(h, prepare_gradient(h, MC11E, F), MC11E, F)             # overlay under Mooncake
        @test all(a -> all(==(42.0), a), (w.y, w.c, w.z, w.b, w.d))
        @test_throws DimensionMismatch gradient(hl, prepare_gradient(hl, MC11E, F), MC11E, F)
        @test all(a -> all(==(42.0), a), (w.y, w.c, w.z, w.b, w.d))
    end
end
