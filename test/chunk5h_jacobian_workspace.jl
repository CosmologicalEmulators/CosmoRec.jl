# Chunk 5h: the caller callback's BUFFERED state Jacobian (`prepare_wsjac5` / `WSJac5`, test/chunk5_helpers.jl; Step 3 increment 2,
# docs/ODE_CORE_WORKSPACE_DESIGN.md §8–§9) vs the unprepared oracle `jac5!` (and the Step-2 `PreparedJac5`).
# Accuracy gates (unchanged): ForwardDiff values/nested derivatives bitwise; prepared Mooncake within 1e-12 of ForwardDiff at in-table states.
# Off-table z = 1500 (DP-fallback absorber): Float64 second-order AD is precision-sensitive in BOTH engines (Mooncake differs from ForwardDiff by
# ~0.5 there with the legacy oracle as well; docs §8.1–§8.2), so the Mooncake check there is BITWISE PARITY WITH THE LEGACY ORACLE'S Mooncake vector:
# a no-new-regression check, NOT an accuracy claim. REQUIRES COSMOREC_NATIVE_DATA_DIR (see test/chunk5_helpers.jl); fails loudly if unset.
using Test
using Random
using LinearAlgebra
using ForwardDiff
using Mooncake
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake
using CosmoRec
@isdefined(RM5) || include("chunk5_helpers.jl")

const MC5H = AutoMooncake(; config = nothing)
const ROWS5H = read_fixture5a(joinpath(@__DIR__, "fixtures", "native_ode_rhs.txt"))
Tof5h(u, p, z) = promote_type(eltype(u), typeof(z), typeof(p.hscale), typeof(p.nbscale))
oracle5h(u, p, z) = (T = Tof5h(u, p, z); J = zeros(T, length(u), length(u)); jac5!(J, T.(u), p, z); J)
function ws5h(u, p, z)                                   # a fresh buffered cache from these actual types (as prob_phase5 builds per solve)
    T = Tof5h(u, p, z); uu = T.(u); j = prepare_wsjac5(uu, p, z); J = zeros(T, length(u), length(u)); j(J, uu, p, z)
    @test wsjac5_fallbacks(j) == 0
    return J
end
pofh5h(flag, h) = RecombinationODEParams(RM5, flag, h, one(h))

@testset "Chunk 5h: buffered state-Jacobian workspace (family tag) vs the unprepared oracle" begin
    @testset "values and nested ForwardDiff (bitwise), $(r.flag ? 12 : 7)-state z = $(r.z)" for r in
            [first(r for r in ROWS5H if r.flag == f && r.z == z && r.pat == 0) for (f, z) in ((true, 2500.0), (true, 1500.0), (false, 400.0))]
        p = RecombinationODEParams(RM5, r.flag); z = r.z
        @test ws5h(r.y, p, z) == oracle5h(r.y, p, z)
        @test ForwardDiff.jacobian(u -> vec(ws5h(u, p, z)), r.y) == ForwardDiff.jacobian(u -> vec(oracle5h(u, p, z)), r.y)
        @test ForwardDiff.derivative(zz -> vec(ws5h(r.y, p, zz)), z) == ForwardDiff.derivative(zz -> vec(oracle5h(r.y, p, zz)), z)
        @test ForwardDiff.derivative(h -> vec(ws5h(r.y, pofh5h(r.flag, h), z)), 1.0) == ForwardDiff.derivative(h -> vec(oracle5h(r.y, pofh5h(r.flag, h), z)), 1.0)
        @test ForwardDiff.derivative(a -> ForwardDiff.derivative(zz -> vec(ws5h(r.y, p, zz)), a), z) ==
              ForwardDiff.derivative(a -> ForwardDiff.derivative(zz -> vec(oracle5h(r.y, p, zz)), a), z)
        # two independent caches, live (y, p, z) alternating; no stale parameters
        n = length(r.y); y2 = r.y .* (1 .+ 1.0e-6 .* (1:n)); p2 = pofh5h(r.flag, 1.001); z2 = z - 7.0
        ja = prepare_wsjac5(r.y, p, z); jb = prepare_wsjac5(r.y, p, z); A = zeros(n, n); B = zeros(n, n); C = zeros(n, n); D = zeros(n, n)
        ja(A, r.y, p, z); jb(D, y2, p2, z2); ja(B, y2, p2, z2); ja(C, r.y, p, z)
        @test A == C == oracle5h(r.y, p, z) && B == D == oracle5h(y2, p2, z2) && wsjac5_fallbacks(ja) == 0 && wsjac5_fallbacks(jb) == 0
        # type-mismatched call (Dual z into a Float64-z cache): counted oracle fallback, derivative equal to the oracle's
        jd = prepare_wsjac5(r.y, p, z)
        dfb = ForwardDiff.derivative(zz -> (Jd = zeros(typeof(zz), n, n); jd(Jd, typeof(zz).(r.y), p, zz); vec(Jd)), z)
        @test jd.nfallback[] == 1 && dfb == ForwardDiff.derivative(zz -> vec(oracle5h(r.y, p, zz)), z)
    end
    @testset "family tag: default check passes, ordering pinned (internal tagcount)" begin
        outer = ForwardDiff.Tag(identity, Float64); inner = bufrhs5_tag(ForwardDiff.Dual{typeof(outer),Float64,1})
        @test ForwardDiff.:≺(typeof(outer), typeof(inner)) && !ForwardDiff.:≺(typeof(inner), typeof(outer))
        r = first(r for r in ROWS5H if r.flag && r.z == 2500.0 && r.pat == 0); j = prepare_wsjac5(r.y, RecombinationODEParams(RM5, true), r.z)
        @test ForwardDiff.checktag(j.cfg, j.f, r.y)                                  # the default check jacobian! performs
        @test typeof(j.cfg).parameters[1] === ForwardDiff.Tag{BufRHS5,Float64}
    end
    @testset "prepared DI-Mooncake, ONE tape with TWO Jacobian evaluations through one prebuilt cache, $(r.flag ? 12 : 7)-state z = $(r.z)" for r in
            [first(r for r in ROWS5H if r.flag == f && r.z == z && r.pat == 0) for (f, z) in ((true, 2500.0), (false, 400.0), (true, 1500.0))]
        p = RecombinationODEParams(RM5, r.flag); z = r.z; n = length(r.y)
        W1 = randn(MersenneTwister(3), n, n); W2 = randn(MersenneTwister(4), n, n); cs = 1 .+ 1.0e-6 .* (1:n)
        pre = prepare_wsjac5(r.y, p, z)
        Jw(u) = (J = zeros(eltype(u), n, n); pre(J, u, p, z); J); Jo(u) = (J = zeros(eltype(u), n, n); jac5!(J, u, p, z); J)
        fw(u) = dot(W1, Jw(u)) + dot(W2, Jw(u .* cs)); fo(u) = dot(W1, Jo(u)) + dot(W2, Jo(u .* cs))
        @test fw(r.y) == fo(r.y)
        pw = prepare_gradient(fw, MC5H, r.y); g1 = gradient(fw, pw, MC5H, r.y); g2 = gradient(fw, pw, MC5H, r.y)
        @test g1 == g2 && wsjac5_fallbacks(pre) == 0
        if z == 1500.0
            # off-table: no-new-regression parity with the legacy oracle's Mooncake vector (Float64 accuracy here is not established for either engine)
            po = prepare_gradient(fo, MC5H, r.y); @test g1 == gradient(fo, po, MC5H, r.y)
        else
            ref = ForwardDiff.gradient(fo, r.y)
            @test maximum(abs.(g1 .- ref)) / maximum(abs, ref) < 1e-12        # accuracy gate, unchanged
        end
    end
    @testset "high precision: BigFloat(128) directional second derivative, worst component k = 10 at z = 1500, new == legacy" begin
        r = first(r for r in ROWS5H if r.flag && r.z == 1500.0 && r.pat == 0); p = RecombinationODEParams(RM5, true); n = 12
        W1 = randn(MersenneTwister(3), n, n); ek = [i == 10 ? 1 : 0 for i in 1:n]
        dn, do_ = setprecision(BigFloat, 128) do
            yb = BigFloat.(r.y)
            (ForwardDiff.derivative(t -> dot(W1, ws5h(yb .+ t .* ek, p, r.z)), big(0.0)), ForwardDiff.derivative(t -> dot(W1, oracle5h(yb .+ t .* ek, p, r.z)), big(0.0)))
        end
        @test dn == do_
    end
    @testset "prepared DI-Mooncake d/dz and d/dhscale through ONE prebuilt buffered cache (live p/z mutation), $(r.flag ? 12 : 7)-state z = $(r.z)" for r in
            [first(r for r in ROWS5H if r.flag == f && r.z == z && r.pat == 0) for (f, z) in ((true, 2500.0), (false, 400.0))]     # in-table: accuracy gate unchanged
        n = length(r.y); W = randn(MersenneTwister(6), n, n); p0 = RecombinationODEParams(RM5, r.flag)
        pre = prepare_wsjac5(r.y, p0, r.z)                      # Float64 cache; every call below writes its own (p, z) into it
        Jp(pp, zz) = (J = zeros(promote_type(typeof(zz), typeof(pp.hscale)), n, n); pre(J, promote_type(typeof(zz), typeof(pp.hscale)).(r.y), pp, zz); J)
        Jo(pp, zz) = oracle5h(r.y, pp, zz)
        for (lab, fw, fo, x0, x1) in (("z", x -> dot(W, Jp(p0, x[1])), x -> dot(W, Jo(p0, x[1])), [r.z], [r.z - 3.0]),
                                       ("hscale", x -> dot(W, Jp(pofh5h(r.flag, x[1]), r.z)), x -> dot(W, Jo(pofh5h(r.flag, x[1]), r.z)), [1.0], [1.002]))
            prep = prepare_gradient(fw, MC5H, x0)
            for x in (x0, x1, x0)                               # the same preparation re-used at changed z / hscale: no stale parameter in the cache
                g = gradient(fw, prep, MC5H, x); ref = ForwardDiff.gradient(fo, x)
                @test abs(g[1] - ref[1]) / max(abs(ref[1]), floatmin()) < 1e-12
            end
            @test gradient(fw, prep, MC5H, x0) == gradient(fw, prep, MC5H, x0)
        end
        @test wsjac5_fallbacks(pre) == 0                        # Mooncake runs the primal types: every call took the buffered path
    end
    @testset "allocations per Jacobian call: buffered below the Step-2 PreparedJac5" begin
        for (f, z) in ((true, 2500.0), (false, 400.0))
            r = first(r for r in ROWS5H if r.flag == f && r.z == z && r.pat == 0); p = RecombinationODEParams(RM5, f); n = length(r.y); J = zeros(n, n)
            jw = prepare_wsjac5(r.y, p, z); jp = PreparedJac5(r.y, p, z)
            aw() = @allocated jw(J, r.y, p, z); ap() = @allocated jp(J, r.y, p, z)
            aw(); ap(); aw(); ap()
            a_w = aw(); a_p = ap(); println("Chunk5h allocated bytes per Jacobian call $(n)-state: buffered $a_w, PreparedJac5 $a_p")
            @test a_w < a_p
        end
    end
end
