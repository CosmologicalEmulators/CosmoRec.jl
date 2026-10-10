# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Chunk 11 (opt-in, private): fixed-geometry spline quadrature of the PDE correction integrals.
# For a fixed knot slice, fixed sub-interval edges and the Patterson level the primal chose, the map from the spectrum ordinates F to each
# sub-integral is linear. `_hi_spline_parts` evaluates it with the UNCHANGED primal code (natural_cubic_spline + patterson_integrate); the
# precomputed `ell[:, i, L]` is its gradient w.r.t. F (adjoint of the IMPLEMENTED stored-factor fit, see `_hi_fit_adjoint!`), used only by the
# reverse rules of the optional ChainRulesCore/Mooncake extensions. Grid, edges and atomic constants are constants of a plan (they are Float64
# in `HIPDESetup`/`HIPDEAtom`); `_check_quadplan` rejects a plan whose copied grid/edges differ from the model's.
# Design: cmbcheb_test/local_analysis/cosmorec_differentiability_20260930/chunk28/RRULE_MOONCAKE_DESIGN.md (revision 2).
# ----------------------------------------------------------------------------------------------------------------------------------------------------

"""Fixed-geometry quadrature plan of one correction-integral block: private knot copy, `K` sub-intervals `[a[i], b[i]]` (empty if `a[i] >= b[i]`),
the first `epsabs`, the stored Float64 factors of the GSL fit on these knots and `ell[:, i, L] = ∂(sub-integral i at Patterson level L)/∂F`."""
struct HIQuadPlan{K}
    knots::Vector{Float64}
    a::NTuple{K,Float64}
    b::NTuple{K,Float64}
    epsabs0::Float64
    dx::Vector{Float64}
    g::Vector{Float64}
    off::Vector{Float64}
    alpha::Vector{Float64}
    t::Vector{Float64}
    ell::Array{Float64,3}
end

"""Quadrature plans of the four production blocks (3s and 3d share `two3`) plus a private copy of the full grid for `_check_quadplan`."""
struct HIQuadPlans{K2,KR}
    di1::HIQuadPlan{1}
    two3::HIQuadPlan{K2}
    raman::HIQuadPlan{KR}
    xfull::Vector{Float64}
    index_emission::NTuple{2,Int}
    nresmax::Int
end

# block limits and edges of the default path (`_hi_limits_*`, `_hi_resonance_edges` in HIPDEIntegrals.jl, shared by both paths)
function _hi_block_edges(s::HIPDESetup, atom::HIPDEAtom)
    x = s.x; np = length(x); ie = s.index_emission; k0 = ie[1]
    l1 = _hi_limits_di1(x, ie[1])
    l2 = _hi_limits_two3(x, ie[2], atom)
    lr = _hi_limits_raman(x, k0, np - k0, atom)
    return (di1 = [l1[1], l1[2]], two3 = _hi_resonance_edges(atom, l2[1], l2[2], 2, 2), raman = _hi_resonance_edges(atom, lr[1], lr[2], 3, s.nresmax))
end

# stored factors of `natural_cubic_spline` on these knots (same expressions, CosmosAccessors.jl:29-44)
function _hi_fit_factors(xs::Vector{Float64})
    n = length(xs); N = n - 2
    dx = Vector{Float64}(undef, n - 1); g = Vector{Float64}(undef, n - 1)
    for i in 1:(n - 1)
        h = xs[i + 1] - xs[i]
        dx[i] = h; g[i] = h != 0.0 ? 1.0 / h : 0.0
    end
    off = Vector{Float64}(undef, N); diag = Vector{Float64}(undef, N)
    for i in 1:N
        off[i] = dx[i + 1]; diag[i] = 2.0 * (dx[i + 1] + dx[i])
    end
    alpha = Vector{Float64}(undef, N); t = zeros(N)
    alpha[1] = diag[1]
    for i in 2:N
        t[i] = off[i - 1] / alpha[i - 1]
        alpha[i] = diag[i] - t[i] * off[i - 1]
    end
    return dx, g, off, alpha, t
end

"""
    _hi_fit_adjoint!(ybar, bbar, cbar, dbar, dx, g, off, alpha, t) -> ybar

Adjoint of the IMPLEMENTED GSL fit y ↦ (b, c, d) with its stored Float64 factors: b/d stencils, natural ends (c̄₁, c̄ₙ dropped), the exact
TRANSPOSE of the stored-factor solve u = Û⁻¹ L̂⁻¹ γ (Ûᵀ w = c̄, then L̂ᵀ γ̄ = w — not a re-application of the forward solve, whose rounded
factors are not exactly symmetric), then the γ stencil. Accumulates into `ybar` (which may hold direct ordinate contributions); overwrites `cbar`.
"""
function _hi_fit_adjoint!(ybar, bbar, cbar, dbar, dx, g, off, alpha, t)
    n = length(ybar); N = n - 2
    @inbounds for j in 1:(n - 1)
        h = dx[j]
        ybar[j + 1] += bbar[j] / h
        ybar[j] -= bbar[j] / h
        cbar[j + 1] += -(h * bbar[j]) / 3.0 + dbar[j] / (3.0 * h)
        cbar[j] += -(2.0 * h * bbar[j]) / 3.0 - dbar[j] / (3.0 * h)
    end
    # w̄ = cbar[2:n-1] (u_i = c_{i+1}); Ûᵀ w = w̄ (forward substitution), stored in cbar[2:n-1]
    @inbounds begin
        cbar[2] = cbar[2] / alpha[1]
        for i in 2:N
            cbar[i + 1] = (cbar[i + 1] - off[i - 1] * cbar[i]) / alpha[i]
        end
        # L̂ᵀ γ̄ = w (back substitution)
        for i in (N - 1):-1:1
            cbar[i + 1] = cbar[i + 1] - t[i + 1] * cbar[i + 2]
        end
        for i in 1:N
            G = 3.0 * cbar[i + 1]
            ybar[i + 2] += G * g[i + 1]
            ybar[i + 1] -= G * (g[i + 1] + g[i])
            ybar[i] += G * g[i]
        end
    end
    return ybar
end

"""Patterson entries of level `L` (1-based into `PATTERSON_NPOINTS`): `(offset, weight)` with offset 0.0 for the centre, otherwise the
positive node offset in units of `Dx`; reproduces the interleaving of `patterson_integrate` (DPescCoh.jl:367-393)."""
function _patterson_entries(L::Int)
    L == 1 && return [(0.0, 2.0)]                 # r = 2 f(xc) Dx
    ent = [0.0]
    for it in 2:L
        xnew = PATTERSON_X[it]
        loc = Vector{Float64}(undef, 2 * length(xnew))
        for k in eachindex(xnew)
            loc[2k - 1] = ent[k]; loc[2k] = xnew[k]
        end
        ent = loc
    end
    return [(ent[j], PATTERSON_W[L][j]) for j in eachindex(ent)]
end

function _hi_quad_plan(knots::Vector{Float64}, edges::Vector{Float64}, epsabs0::Float64)
    n = length(knots); K = length(edges) - 1
    n >= 3 || throw(ArgumentError("_hi_quad_plan: need at least 3 knots, got $n"))
    all(isfinite, knots) || throw(ArgumentError("_hi_quad_plan: knots must be finite"))
    all(i -> knots[i + 1] > knots[i], 1:(n - 1)) || throw(ArgumentError("_hi_quad_plan: knots must be strictly increasing"))
    K >= 1 || throw(ArgumentError("_hi_quad_plan: need at least 2 edges (one sub-interval), got $(length(edges))"))
    all(isfinite, edges) || throw(ArgumentError("_hi_quad_plan: edges must be finite"))
    for i in 1:K
        a = edges[i]; b = edges[i + 1]
        a >= b && continue                                   # empty sub-interval (native convention), contributes zero
        (knots[1] <= a && b <= knots[end]) || throw(ArgumentError("_hi_quad_plan: non-empty sub-interval $i outside the knot range"))
    end
    (isfinite(epsabs0) && epsabs0 >= 0.0) || throw(ArgumentError("_hi_quad_plan: epsabs0 must be finite and non-negative"))
    xs = copy(knots)
    dx, g, off, alpha, t = _hi_fit_factors(xs)
    Slo = NaturalCubicSpline{Float64}(xs, Float64.(1:n), zeros(n - 1), zeros(n), zeros(n - 1))   # spline_eval_native(Slo, x) == lo exactly
    Sdel = NaturalCubicSpline{Float64}(xs, zeros(n), ones(n - 1), zeros(n), zeros(n - 1))        # spline_eval_native(Sdel, x) == δ exactly
    ell = zeros(n, K, 8)
    ybar = zeros(n); bbar = zeros(n - 1); cbar = zeros(n); dbar = zeros(n - 1)
    for i in 1:K
        a = edges[i]; b = edges[i + 1]
        a >= b && continue
        xc = (a + b) / 2.0; Dx = (b - a) / 2.0
        for L in 1:8
            fill!(ybar, 0.0); fill!(bbar, 0.0); fill!(cbar, 0.0); fill!(dbar, 0.0)
            for (u, w) in _patterson_entries(L)
                om = w * Dx
                nodes = u == 0.0 ? (xc,) : (Del = Dx * u; (xc + Del, xc - Del))
                for xq in nodes
                    lo = Int(spline_eval_native(Slo, xq)); d = spline_eval_native(Sdel, xq)
                    ybar[lo] += om; bbar[lo] += om * d; cbar[lo] += om * d * d; dbar[lo] += om * d * d * d
                end
            end
            _hi_fit_adjoint!(ybar, bbar, cbar, dbar, dx, g, off, alpha, t)
            ell[:, i, L] .= ybar
        end
    end
    return HIQuadPlan{K}(xs, Tuple(edges[1:K]), Tuple(edges[2:(K + 1)]), epsabs0, dx, g, off, alpha, t, ell)
end

"""
    _hi_quad_plans(setup, atom = NATIVE_HI_PDE_ATOM) -> HIQuadPlans

Build (outside any AD tape) the fixed-geometry plans of the production correction integrals. Pass the result as the opt-in `quadplan` keyword of
`hi_pde_corrections` / `hi_diffusion_stage` / `recombination_history_diffusion`; it is rejected (`_check_quadplan`) if grid or edges differ.
"""
function _hi_quad_plans(s::HIPDESetup, atom::HIPDEAtom = NATIVE_HI_PDE_ATOM)
    x = s.x; np = length(x); ie = s.index_emission; k0 = ie[1]
    e = _hi_block_edges(s, atom)
    di1 = _hi_quad_plan(x[1:ie[1]], e.di1, 1.0e-60)
    two3 = _hi_quad_plan(x[1:ie[2]], e.two3, 1.0e-60)
    raman = _hi_quad_plan(x[(k0 + 1):np], e.raman, 1.0e-60)
    return HIQuadPlans(di1, two3, raman, copy(x), (ie[1], ie[2]), s.nresmax)
end

"""Passive full-grid validation of a plan against the model's setup/atom (every knot, `index_emission`, block knot slices and edge tuples);
throws `ArgumentError` on any difference. Declared zero-derivative under Mooncake (no tangent flow, no tracing)."""
function _check_quadplan(q::HIQuadPlans, s::HIPDESetup, atom::HIPDEAtom)
    x = s.x; np = length(x); ie = s.index_emission
    bad(msg) = throw(ArgumentError("stale HI quadrature plan: " * msg * " (rebuild with CosmoRec._hi_quad_plans(setup, atom))"))
    length(q.xfull) == np || bad("grid length")
    for j in 1:np
        isequal(q.xfull[j], x[j]) || bad("grid node $j differs")
    end
    (length(ie) == 2 && q.index_emission == (ie[1], ie[2])) || bad("index_emission")
    q.nresmax == s.nresmax || bad("nresmax")
    k0 = ie[1]
    for (p, r) in ((q.di1, 1:ie[1]), (q.two3, 1:ie[2]), (q.raman, (k0 + 1):np))
        length(p.knots) == length(r) || bad("block knot count")
        for (j, k) in enumerate(r)
            isequal(p.knots[j], x[k]) || bad("block knot $k differs")
        end
    end
    e = _hi_block_edges(s, atom); e2 = e.two3; er = e.raman
    (isequal(q.di1.a[1], e.di1[1]) && isequal(q.di1.b[1], e.di1[2])) || bad("DI1 edges")
    (length(e2) == length(q.two3.a) + 1 && all(isequal(q.two3.a[i], e2[i]) && isequal(q.two3.b[i], e2[i + 1]) for i in eachindex(q.two3.a))) || bad("3s/3d edges")
    (length(er) == length(q.raman.a) + 1 && all(isequal(q.raman.a[i], er[i]) && isequal(q.raman.b[i], er[i + 1]) for i in eachindex(q.raman.a))) || bad("Raman edges")
    return nothing
end

"""
    _hi_fit_cached(q, F) -> NaturalCubicSpline

Private: the natural cubic spline of the ACTIVE ordinates `F` on the plan's knots, fitted with the plan's stored Float64 geometry (`dx`, `g`, `off`, `alpha`, `t`, computed once by
`_hi_fit_factors` with the same expressions as `natural_cubic_spline`) instead of recomputing it for every one of the integrals of every output. Operation order is that of
`natural_cubic_spline` (`gamma = 3 (dy g - dy' g')`, the stored LDL factors, `b = dy/dx - dx (c + 2c')/3` and `d = Δc/(3 dx)` with divisions, not products with `g`), so all fields are
BITWISE equal to it for Float64, Dual and BigFloat ordinates alike (the stored factors stay Float64, as in the constructor). `y`, `b`, `c`, `d` are fresh arrays owned by the result
(`y` is a copy of `F`); only the read-only `knots` of the plan are shared (a spline never mutates them), and the plan itself is never written, so one plan may serve concurrent calls.
Contract: the plan geometry is validated once by `_hi_quad_plan` (`_check_quadplan` guards the stage path); the plan's vectors are mutable, so each call re-checks in O(n) the lengths of the stored
factors against `F`, that the ACTUAL knots are strictly increasing (ArgumentError, as the generic constructor) and that every stored `dx[i]` equals the actual `xs[i+1] - xs[i]` (rejects stale factors
of sorted-but-shifted knots). It does not re-derive `g`, `off`, `alpha`, `t` (no refactorisation, no cache). Narrow element types: `gamma` is converted to `T` at the point where the constructor stores it.
"""
function _hi_fit_cached(q::HIQuadPlan, F::AbstractVector{T}) where {T<:Real}
    xs = q.knots; n = length(xs); N = n - 2
    dx = q.dx; g = q.g; off = q.off; alpha = q.alpha; tt = q.t
    (n >= 3 && length(F) == n && length(dx) == n - 1 && length(g) == n - 1 && length(off) == N && length(alpha) == N && length(tt) == N) ||
        throw(DimensionMismatch("_hi_fit_cached: need length(F) == length(plan knots) >= 3 and stored factors of matching length"))
    @inbounds for i in 1:(n - 1)                  # actual knots strictly increasing (the generic constructor's domain error) AND the stored spacing is theirs
        h = xs[i + 1] - xs[i]
        h > 0.0 || throw(ArgumentError("_hi_fit_cached: plan knots must be strictly increasing (natural_cubic_spline: x must be strictly increasing)"))
        h == dx[i] || throw(ArgumentError("_hi_fit_cached: stored fit factors are stale (dx[$i] differs from the knot spacing); rebuild the plan"))
    end
    y = collect(F)
    c = zeros(T, n); z = Vector{T}(undef, N)
    for i in 1:N
        gam = convert(T, 3.0 * ((y[i + 2] - y[i + 1]) * g[i + 1] - (y[i + 1] - y[i]) * g[i]))      # natural_cubic_spline stores gamma in a Vector{T} BEFORE the forward substitution
        z[i] = i == 1 ? gam : gam - tt[i] * z[i - 1]
    end
    c[N + 1] = z[N] / alpha[N]
    for i in (N - 1):-1:1
        c[i + 1] = (z[i] - off[i] * c[i + 2]) / alpha[i]
    end
    b = Vector{T}(undef, n - 1); d = Vector{T}(undef, n - 1)
    for i in 1:(n - 1)
        h = dx[i]
        b[i] = (y[i + 1] - y[i]) / h - h * (c[i + 1] + 2.0 * c[i]) / 3.0
        d[i] = (c[i + 1] - c[i]) / (3.0 * h)
    end
    return NaturalCubicSpline{T}(xs, y, b, c, d)
end

"""
    HIFitScratch{T}(n)

Caller-owned (per stage / per chain / per task) scratch of the plan-path spline fit for ONE knot slice of `n` knots and ordinate type `T`: the five work/result arrays of
`_hi_fit_cached` (`y`, `c`, `z`, `b`, `d`), allocated once. It is NOT part of the immutable plan and is never global: independent scratches may be used concurrently with one
shared plan, ONE scratch must not be used by two tasks at once. `nfallback` counts calls served by the allocating oracle `_hi_fit_cached` (element-type mismatch, wrong scratch size,
aliasing of the ordinates with a scratch array); nothing is ever coerced. Initial contents are irrelevant: every array is written completely on each call before it is read.
"""
struct HIFitScratch{T<:Real}
    y::Vector{T}
    c::Vector{T}
    z::Vector{T}
    b::Vector{T}
    d::Vector{T}
    nfallback::Base.RefValue{Int}
end
function HIFitScratch{T}(n::Integer) where {T<:Real}
    n >= 3 || throw(ArgumentError("HIFitScratch: need n >= 3 knots"))
    return HIFitScratch{T}(Vector{T}(undef, n), Vector{T}(undef, n), Vector{T}(undef, n - 2), Vector{T}(undef, n - 1), Vector{T}(undef, n - 1), Ref(0))
end
HIFitScratch(::Type{T}, n::Integer) where {T<:Real} = HIFitScratch{T}(n)

"""Scratches of the three production blocks (`di1`, `two3` shared by 3s and 3d, which run one after the other, `raman`) for ordinates of type `T`; built from a plan set."""
struct HIQuadWork{T<:Real}
    di1::HIFitScratch{T}
    two3::HIFitScratch{T}
    raman::HIFitScratch{T}
    nmismatch::Base.RefValue{Int}              # blocks served without scratch because the ordinate type differs from `T`
end
HIQuadWork(q::HIQuadPlans, ::Type{T}) where {T<:Real} =
    HIQuadWork{T}(HIFitScratch{T}(length(q.di1.knots)), HIFitScratch{T}(length(q.two3.knots)), HIFitScratch{T}(length(q.raman.knots)), Ref(0))

# stage plumbing: `fitwork = nothing` (default: allocating fit, unchanged) | `true` (a stage-owned `HIQuadWork` built here for ordinate type `T`) | a caller-owned `HIQuadWork` (used as is, never copied).
# Anything else (including `false`) is an ArgumentError, as is any `fitwork` without the opt-in `quadplan`.
_hi_stage_work(::Nothing, quadplan, ::Type) = nothing
function _hi_stage_work(fitwork, quadplan, ::Type{T}) where {T}
    quadplan === nothing && throw(ArgumentError("fitwork requires the opt-in quadplan path"))
    fitwork === true && return HIQuadWork(quadplan, T)
    fitwork isa HIQuadWork && return fitwork
    throw(ArgumentError("fitwork must be nothing, true or a CosmoRec HIQuadWork"))
end
_hi_scratch_for(::Nothing, ::Symbol, ::AbstractVector) = nothing
function _hi_scratch_for(w::HIQuadWork{T}, blk::Symbol, F::AbstractVector) where {T}
    eltype(F) === T || (w.nmismatch[] += 1; return nothing)
    return getfield(w, blk)
end

"""
    _hi_fit_scratched(w, q, F) -> NaturalCubicSpline{T}

Same spline as `_hi_fit_cached(q, F)` (same validation and error messages, same operation order, the `gamma` conversion to `T` before the forward substitution, natural zero ends),
fitted into the scratch `w`. The result is an EPHEMERAL BORROWED view of `w`'s arrays (and of the plan's knots): valid only until the next fit into `w`; `_hi_spline_parts` consumes it within
the call. Anything the scratch cannot serve identically (ordinates of another element type, scratch of another size, ordinates aliasing a scratch array) is served by the allocating `_hi_fit_cached`
and counted in `w.nfallback`; nothing is coerced. AD: under Mooncake `_hi_spline_parts(q, F, forced, w)` is overlaid (see ext/CosmoRecMooncakeExt.jl) so that every later read of the scratch is differentiated.
"""
function _hi_fit_scratched(w::HIFitScratch{T}, q::HIQuadPlan, F::AbstractVector{T}) where {T<:Real}
    xs = q.knots; n = length(xs); N = n - 2
    dx = q.dx; g = q.g; off = q.off; alpha = q.alpha; tt = q.t
    (n >= 3 && length(F) == n && length(dx) == n - 1 && length(g) == n - 1 && length(off) == N && length(alpha) == N && length(tt) == N) ||
        throw(DimensionMismatch("_hi_fit_cached: need length(F) == length(plan knots) >= 3 and stored factors of matching length"))
    @inbounds for i in 1:(n - 1)
        h = xs[i + 1] - xs[i]
        h > 0.0 || throw(ArgumentError("_hi_fit_cached: plan knots must be strictly increasing (natural_cubic_spline: x must be strictly increasing)"))
        h == dx[i] || throw(ArgumentError("_hi_fit_cached: stored fit factors are stale (dx[$i] differs from the knot spacing); rebuild the plan"))
    end
    y = w.y; c = w.c; z = w.z; b = w.b; d = w.d
    if length(y) != n || length(c) != n || length(z) != N || length(b) != n - 1 || length(d) != n - 1 ||
       Base.mightalias(F, y) || Base.mightalias(F, c) || Base.mightalias(F, z) || Base.mightalias(F, b) || Base.mightalias(F, d)
        w.nfallback[] += 1
        return _hi_fit_cached(q, F)
    end
    copyto!(y, F); fill!(c, zero(T))
    for i in 1:N
        gam = convert(T, 3.0 * ((y[i + 2] - y[i + 1]) * g[i + 1] - (y[i + 1] - y[i]) * g[i]))
        z[i] = i == 1 ? gam : gam - tt[i] * z[i - 1]
    end
    c[N + 1] = z[N] / alpha[N]
    for i in (N - 1):-1:1
        c[i + 1] = (z[i] - off[i] * c[i + 2]) / alpha[i]
    end
    for i in 1:(n - 1)
        h = dx[i]
        b[i] = (y[i + 1] - y[i]) / h - h * (c[i + 1] + 2.0 * c[i]) / 3.0
        d[i] = (c[i + 1] - c[i]) / (3.0 * h)
    end
    return NaturalCubicSpline{T}(xs, y, b, c, d)
end

"""
    _hi_spline_parts(q, F, forced) -> (parts, levels, ok)

Pure, generic: the `K` sub-integrals of the natural cubic spline of `F` on the plan's knots, computed by the unchanged primal arithmetic
(the plan's cached-geometry fit `_hi_fit_cached`, bitwise `natural_cubic_spline`; `patterson_integrate` with the native epsabs chain and the level-7 `ok` test of `_hi_calc_DF`). `forced[i] > 0` replays that
level. `levels[i] == 0` marks an empty sub-interval (contributes zero, not recorded). Reverse rules: optional extensions (Float64 `Vector` only).
"""
# private level contract: 0 = adaptive (forced) / empty sub-interval (returned levels), 1:8 = Patterson rule index (PATTERSON_NPOINTS)
function _hi_check_levels(levels::NTuple{K,Int}, what) where {K}
    for L in levels
        0 <= L <= length(PATTERSON_NPOINTS) || throw(ArgumentError("$what: Patterson level $L outside 0:$(length(PATTERSON_NPOINTS))"))
    end
    return levels
end

function _hi_spline_parts(q::HIQuadPlan{K}, F::AbstractVector, forced::NTuple{K,Int}) where {K}
    length(F) == length(q.knots) || throw(DimensionMismatch("_hi_spline_parts: length(F) != number of plan knots"))
    _hi_check_levels(forced, "_hi_spline_parts forced levels")
    s = _hi_fit_cached(q, F)                      # bitwise natural_cubic_spline(q.knots, F), with the plan's stored geometry
    return _hi_parts_from_spline(q, s, forced)
end

_hi_fit_scratched(w::HIFitScratch, q::HIQuadPlan, F::AbstractVector) = (w.nfallback[] += 1; _hi_fit_cached(q, F))      # other ordinate type / SubArray: counted allocating oracle

"""
    _hi_spline_parts(q, F, forced, w::HIFitScratch{T}) -> (parts, levels, ok)

Same value, levels and flags as `_hi_spline_parts(q, F, forced)` (it shares `_hi_parts_from_spline`), with the spline fitted into the caller-owned scratch `w` (`_hi_fit_scratched`). `F` must be a
`Vector{T}`; any other combination (SubArray ordinates, element type different from the scratch's) takes the ordinary allocating method and is counted in `w.nfallback`. Reverse mode: the optional
Mooncake extension overlays this method (traced fit + the 3-argument primitive); there is deliberately no ChainRulesCore rule for it.
"""
function _hi_spline_parts(q::HIQuadPlan{K}, F::Vector{T}, forced::NTuple{K,Int}, w::HIFitScratch{T}) where {K,T<:Real}
    length(F) == length(q.knots) || throw(DimensionMismatch("_hi_spline_parts: length(F) != number of plan knots"))
    _hi_check_levels(forced, "_hi_spline_parts forced levels")
    return _hi_parts_from_spline(q, _hi_fit_scratched(w, q, F), forced)
end
function _hi_spline_parts(q::HIQuadPlan{K}, F::AbstractVector, forced::NTuple{K,Int}, w::HIFitScratch) where {K}
    w.nfallback[] += 1
    return _hi_spline_parts(q, F, forced)
end

# the quadrature of the fitted spline (unchanged primal arithmetic; moved out of `_hi_spline_parts` so the scratch and allocating fits share it)
function _hi_parts_from_spline(q::HIQuadPlan{K}, s::NaturalCubicSpline, forced::NTuple{K,Int}) where {K}
    T = eltype(s.y)
    f = t -> spline_eval_native(s, t)
    epsabs = q.epsabs0
    r = zero(T)
    parts = ntuple(_ -> zero(T), Val(K)); levs = ntuple(_ -> 0, Val(K)); oks = ntuple(_ -> true, Val(K))
    for i in 1:K
        a = q.a[i]; b = q.b[i]
        if a >= b
            ri = zero(T); li = 0; oki = true
        elseif forced[i] > 0
            ri, _ = patterson_integrate(f, a, b, HI_EPSREL_DF, epsabs; level = forced[i])
            li = forced[i]; oki = true
        else
            ri, li = patterson_integrate(f, a, b, HI_EPSREL_DF, epsabs)
            oki = true
            if li >= 8
                r7, _ = patterson_integrate(f, a, b, HI_EPSREL_DF, epsabs; level = 7)
                oki = abs(_primal(ri) - _primal(r7)) <= max(epsabs, abs(_primal(ri)) * HI_EPSREL_DF)
            end
        end
        parts = Base.setindex(parts, ri, i); levs = Base.setindex(levs, li, i); oks = Base.setindex(oks, oki, i)
        r += ri
        epsabs = abs(_primal(r)) * HI_EPSREL_DF
    end
    return parts, levs, oks
end

"""`dF[j] += Σᵢ dparts[i] ⋅ ell[j, i, levels[i]]` over non-empty sub-intervals with non-zero cotangent (in-place accumulation; never reads F)."""
function _hi_quad_pullback!(dF::AbstractVector, q::HIQuadPlan{K}, levels::NTuple{K,Int}, dparts) where {K}
    n = length(q.knots)
    length(dF) == n || throw(DimensionMismatch("_hi_quad_pullback!: length(dF) != number of plan knots"))
    size(q.ell) == (n, K, length(PATTERSON_NPOINTS)) || throw(DimensionMismatch("_hi_quad_pullback!: plan table size"))
    _hi_check_levels(levels, "_hi_quad_pullback! levels")              # validated before the unchecked table access below
    for i in 1:K
        L = levels[i]; w = dparts[i]
        (L == 0 || iszero(w)) && continue
        @inbounds for j in 1:n
            dF[j] += w * q.ell[j, i, L]
        end
    end
    return dF
end

# record/replay bookkeeping (kept in the caller so AD sees the same PattersonLevels mutation as on the default path)
_hi_forced_levels!(::Nothing, ::HIQuadPlan{K}) where {K} = ntuple(_ -> 0, Val(K))
function _hi_forced_levels!(rec::PattersonLevels, q::HIQuadPlan{K}) where {K}
    rec.mode === :replay || return ntuple(_ -> 0, Val(K))
    return ntuple(Val(K)) do i
        q.a[i] >= q.b[i] && return 0
        rec.pos += 1
        L = rec.levels[rec.pos]
        1 <= L <= length(PATTERSON_NPOINTS) || throw(ArgumentError("replayed Patterson level $L at position $(rec.pos) outside 1:$(length(PATTERSON_NPOINTS))"))
        L
    end
end
_hi_record_levels!(::Nothing, ::HIQuadPlan, levels) = nothing
function _hi_record_levels!(rec::PattersonLevels, q::HIQuadPlan{K}, levels::NTuple{K,Int}) where {K}
    rec.mode === :replay && return nothing
    for i in 1:K
        levels[i] > 0 && push!(rec.levels, levels[i])
    end
    return nothing
end

# plan path of the single-interval DI1 integral (`_hi_calc_DF`): (value, converged)
_hi_parts(q, F, forced, ::Nothing) = _hi_spline_parts(q, F, forced)
_hi_parts(q, F, forced, w::HIFitScratch) = _hi_spline_parts(q, F, forced, w)
function _hi_plan_single(q::HIQuadPlan{1}, F::AbstractVector, rec, w = nothing)
    parts, levs, oks = _hi_parts(q, F, _hi_forced_levels!(rec, q), w)
    _hi_record_levels!(rec, q, levs)
    return parts[1], oks[1]
end

# plan path of `_hi_integral_over_resonances` (same accumulation order; parts/ok from `_hi_spline_parts`)
function _hi_plan_resonances(q::HIQuadPlan{K}, F::AbstractVector, rec, w = nothing) where {K}
    parts, levs, oks = _hi_parts(q, F, _hi_forced_levels!(rec, q), w)
    _hi_record_levels!(rec, q, levs)
    r = zero(eltype(F)); sabs = 0.0; nunc = 0
    for i in 1:K
        r += parts[i]; sabs += abs(_primal(parts[i])); nunc += !oks[i]
    end
    return r, sabs, nunc
end
