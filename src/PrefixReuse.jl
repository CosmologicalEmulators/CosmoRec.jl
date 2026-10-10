# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Increment 3 (opt-in, default OFF): reuse of the high-z 12-state ODE blocks of the feedback-free source pass by the later passes of one `recombination_history_diffusion` call.
# Above the zero domain of the HI diffusion feedback (`z >= HIDiffusionFeedback.zmax`, exact zero in `_interp_DF`) the RHS of every pass is the same function of the same inputs, so a block whose
# whole solver query interval lies above that boundary and whose complete inputs (nodes, initial state with Dual partials, scales, model objects) are identical returns the identical solution.
# Contract: chunk33_medium_implementation/PREFIX_CONTRACT.md and prefix/DESIGN_CONTRACT.md. Nothing here is exported; no result field is added; no solver work counter is fabricated.
# ----------------------------------------------------------------------------------------------------------------------------------------------------

"""
    PrefixReuse(solve; signatures, max_query_excursion, pure_solver)

DECLARATION (holds no solutions, may be shared) that the caller solve object `solve` is pure in its arguments `(f!, u0, z0, znodes, p)` (instrumentation side effects allowed; no hidden state, no
captured feedback outside `p.rm`), with `signatures` one comparable declaration per pass (`iterations + 1` of them, the final pass included) that fully describe algorithm, tolerances, stops and
configuration of that pass, and `max_query_excursion` the REQUIRED finite non-negative largest distance (redshift units) below a block's last node at which the solve may evaluate the RHS, its
Jacobian or time gradient. `pure_solver` must be `true`. Pass it as `recombination_history_diffusion(...; prefix_reuse = plan, prefix_stats = PrefixReuseStats())`; reuse is OFF otherwise.
Exact object binding (`solve === plan.solve`) only guards against an undeclared solve; it does not prove equivalence, which is the caller's declaration.
"""
struct PrefixReuse{S}
    solve::S
    signatures::Vector{Any}
    max_query_excursion::Float64
    function PrefixReuse(solve::S; signatures, max_query_excursion, pure_solver) where {S}
        pure_solver === true || throw(ArgumentError("PrefixReuse: `pure_solver = true` must be declared explicitly (an impure or undeclared solve can never be skipped)"))
        (max_query_excursion isa Real && isfinite(max_query_excursion) && max_query_excursion >= 0) ||
            throw(ArgumentError("PrefixReuse: `max_query_excursion` must be a finite non-negative redshift distance, got $max_query_excursion"))
        (signatures isa Union{Tuple,AbstractVector} && !isempty(signatures)) || throw(ArgumentError("PrefixReuse: `signatures` must be a non-empty tuple/vector, one declaration per pass"))
        # A rounded-down declaration would understate the solver's query range.
        excursion = Float64(max_query_excursion)
        isfinite(excursion) || throw(ArgumentError("PrefixReuse: query excursion must be representable as a finite Float64"))
        excursion < max_query_excursion && (excursion = nextfloat(excursion))
        isfinite(excursion) || throw(ArgumentError("PrefixReuse: conservatively rounded query excursion must remain finite"))
        return new{S}(solve, collect(Any, signatures), excursion)
    end
end

"""Caller-owned passive diagnostics of one history call (reset at its start): block counts and refusal reasons. Counts of SKIPPED solve calls are not solver work."""
mutable struct PrefixReuseStats
    status::Symbol                 # :off | :active | :fallback
    reason::Symbol                 # fallback / first refusal reason (:none if reuse ran)
    source_blocks::Int
    consumer_blocks::Int
    reused_blocks::Int
    recomputed_blocks::Int
    refused_excursion::Int
    refused_key::Int
    refused_type::Int
    refused_pass::Int              # whole consumer passes refused (signature / unknown feedback)
end
PrefixReuseStats() = PrefixReuseStats(:off, :none, 0, 0, 0, 0, 0, 0, 0, 0)
function _reset!(s::PrefixReuseStats)
    s.status = :off; s.reason = :none
    s.source_blocks = 0; s.consumer_blocks = 0; s.reused_blocks = 0; s.recomputed_blocks = 0
    s.refused_excursion = 0; s.refused_key = 0; s.refused_type = 0; s.refused_pass = 0
    return s
end

# one cached 12-state block of the source pass: complete inputs (for exact comparison) and the OWNED, differentiably copied solution matrix
struct PrefixBlock{M,U,H,N}
    i::Int
    j::Int
    zn::Vector{Float64}
    u0::U
    hscale::H
    nbscale::N
    sol::M
end
struct PrefixCache
    blocks::Vector{PrefixBlock}
end
PrefixCache() = PrefixCache(PrefixBlock[])

# ephemeral per-pass binding consulted by `recombination_pass` (role :source records, :consumer reuses); `rm_src` is the source model whose eff/cosmos/dp_fallback objects must be identical
struct PrefixBinding{R,P}
    role::Symbol
    cache::PrefixCache
    plan::P
    stats::PrefixReuseStats
    rm_src::R
    zmax::Float64
end

_same_bits(a::AbstractArray, b::AbstractArray) = typeof(a) === typeof(b) && size(a) == size(b) && all(i -> a[i] === b[i], eachindex(a))

function _prefix_solve!(b::PrefixBinding, solve_ode, f!, u0, zn, i::Int, j::Int, p)
    st = b.stats
    znodes = zn[i:j]
    if b.role === :source
        sol = solve_ode(f!, u0, zn[i], znodes, p)
        M = reduce(hcat, [collect(sol[:, m]) for m in 1:(j - i + 1)])          # owned differentiable copy: the caller's solve storage is never shared
        push!(b.cache.blocks, PrefixBlock(i, j, copy(znodes), copy(u0), p.hscale, p.nbscale, M))
        st.source_blocks += 1
        return sol
    end
    st.consumer_blocks += 1
    k = 0
    for q in eachindex(b.cache.blocks)
        b.cache.blocks[q].i == i && (k = q; break)
    end
    if k != 0
        blk = b.cache.blocks[k]
        if !(znodes[end] - b.plan.max_query_excursion >= b.zmax)                    # the block's solver query interval is not proven inside the zero domain
            st.refused_excursion += 1; st.reason === :none && (st.reason = :excursion)
        elseif !(blk.j == j && blk.zn == znodes && _same_bits(u0, blk.u0) && typeof(p.hscale) === typeof(blk.hscale) && p.hscale === blk.hscale &&
                 typeof(p.nbscale) === typeof(blk.nbscale) && p.nbscale === blk.nbscale && p.flag === true &&
                 p.rm.eff === b.rm_src.eff && p.rm.cosmos === b.rm_src.cosmos && p.rm.dp_fallback === b.rm_src.dp_fallback)
            st.refused_key += 1; st.reason === :none && (st.reason = :key)
        elseif !(eltype(blk.sol) === promote_type(eltype(blk.sol), feedback_eltype(p.rm.diffusion)))
            st.refused_type += 1; st.reason === :none && (st.reason = :type)
        else
            st.reused_blocks += 1
            return blk.sol                                                            # read-only: the caller copies columns (`collect(sol[:, m])`)
        end
    else
        st.refused_key += 1; st.reason === :none && (st.reason = :key)
    end
    st.recomputed_blocks += 1
    return solve_ode(f!, u0, zn[i], znodes, p)
end

# history-level control: returns `nothing` (reuse off for this call, `stats` says why) or the plan
function _prefix_begin(plan, stats, solve_pass, d, rm, kwargs)
    stats === nothing || _reset!(stats)
    plan === nothing && return nothing
    plan isa PrefixReuse || throw(ArgumentError("prefix_reuse must be nothing or a CosmoRec.PrefixReuse"))
    length(plan.signatures) == d.iterations + 1 ||
        throw(ArgumentError("prefix_reuse: $(length(plan.signatures)) signatures declared for $(d.iterations + 1) passes (iterations + final pass)"))
    fallback(r) = (stats === nothing || (stats.status = :fallback; stats.reason = r); nothing)
    d.iterations >= 1 || return fallback(:no_consumer_pass)
    haskey(kwargs, :nodes) && kwargs[:nodes] !== nothing && return fallback(:nodes_override)
    plan.solve === solve_pass || return fallback(:solver_binding)
    rm.diffusion === nothing || return fallback(:source_has_feedback)
    stats === nothing || (stats.status = :active)
    return plan
end
_prefix_pass_binding(::Nothing, args...) = nothing
function _prefix_pass_binding(plan::PrefixReuse, stats, cache::PrefixCache, it::Int, rm_src, m)
    it == 0 && return PrefixBinding(:source, cache, plan, stats === nothing ? PrefixReuseStats() : stats, rm_src, Inf)
    st = stats === nothing ? PrefixReuseStats() : stats
    fb = m.diffusion
    if !(fb isa HIDiffusionFeedback)
        st.refused_pass += 1; st.reason === :none && (st.reason = :unknown_feedback); return nothing
    end
    if !isequal(plan.signatures[1], plan.signatures[it + 1])
        st.refused_pass += 1; st.reason === :none && (st.reason = :signature); return nothing
    end
    return PrefixBinding(:consumer, cache, plan, st, rm_src, fb.zmax)
end
