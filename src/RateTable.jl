# Pure-Julia port of the ORIGINAL CosmoRec effective HI rate lookup get_rates().
#
# Original: CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3),
#   Rec_database/Effective_Rates.HI/get_effective_rates.HI.cpp (get_rates: lines 254-347,
#   compute_coefficients_pol: 233-241, log_qnl_qe_func: 243-249, loader: 112-217).
# CosmoRec is (c) J. Chluba et al.; use must be acknowledged and Chluba & Thomas 2010 (MNRAS 412, 748) and
# Chluba & Sunyaev 2006 (A&A 446, 39) cited (see docs/CHUNK2_RESULTS.md NOTICE). Table data are NOT bundled.
#
# The table is an explicit argument (no global state, no implicit file loading). Only get_rates (not
# get_rates_all, which is a different native API) is ported. Interpolation is piecewise: derivatives are
# discontinuous at stencil seams, at the detailed-balance switch and at the eps_A_effective clamp.

const eps_A_effective = 1.0e-4          # cpp:71
const z_detailed_balance = 2000.0       # cpp:72
const T0_detailed_balance = 2.725       # cpp:300 (literal)

"""Physical constants used by `log_qnl_qe` (native physical_consts.h / Definitions.h values)."""
struct RateConstants
    twopi::Float64
    me_gr::Float64
    kB::Float64
    h::Float64
    h_kb::Float64
end

const NATIVE_CONSTANTS = let kB = 1.3806504e-23 * 1.0e+7, h = 6.62606896e-27
    RateConstants(6.2831853071795864769, 9.10938215e-28, kB, h, h / kB)
end

"""Resolved-level data in native order (res_state_list.dat): quantum numbers, spin weight, ionisation frequency, reduced mass."""
struct ResolvedStates
    n::Vector{Int}
    l::Vector{Int}
    gw::Vector{Float64}
    nuion::Vector{Float64}
    mu_red::Vector{Float64}
end

"""Reasons: `:not_finite`, `:out_of_table` (native: message + exit(0)), `:stencil_exceeds_table` (native: exit(0) or undefined out-of-bounds read / segfault)."""
struct RateTableDomainError <: Exception
    reason::Symbol
    Tg::Float64
    Te::Float64
end

function Base.showerror(io::IO, e::RateTableDomainError)
    print(io, "RateTableDomainError(", e.reason, "): (Tg, Te) = (", e.Tg, ", ", e.Te, ")")
end

"""
Explicit table. `lgTg[i] = ln(Tg/K)`, `lgrho[j] = ln(Te/Tg)`; `B[i, m]`, `R[i, m, k]`, `A[i, j, m]` are the native log-space
table entries (native `BiVec[i]`, `RijMatrix[i][k]`, `AiMatrix[i][j]` of state m). `logq_knot[i, m] = log_qnl_qe(m, exp(lgTg[i]))`
is precomputed with the same function and inputs the native code re-evaluates on every call.
"""
struct AtomicRateTable
    lgTg::Vector{Float64}
    lgrho::Vector{Float64}
    B::Matrix{Float64}
    R::Array{Float64,3}
    A::Array{Float64,3}
    states::ResolvedStates
    consts::RateConstants
    logq_knot::Matrix{Float64}
end

function log_qnl_qe(consts::RateConstants, gw::Real, nuion::Real, mu_red::Real, Tg)   # cpp:243-249
    qe = (consts.twopi * consts.me_gr * mu_red * consts.kB * Tg / consts.h / consts.h)^1.5
    gw_gc = gw / 2.0
    qln = gw_gc * exp(consts.h_kb * nuion / Tg)
    return log(qln / qe)
end

log_qnl_qe(t::AtomicRateTable, m::Integer, Tg) =
    log_qnl_qe(t.consts, t.states.gw[m], t.states.nuion[m], t.states.mu_red[m], Tg)

function AtomicRateTable(lgTg, lgrho, B, R, A, states::ResolvedStates, consts::RateConstants = NATIVE_CONSTANTS)
    N, M, nres = length(lgTg), length(lgrho), length(states.n)
    N >= 4 && M >= 4 || throw(ArgumentError("need at least 4 knots per axis"))
    issorted(lgTg; lt = <=) && issorted(lgrho; lt = <=) || throw(ArgumentError("axes must be strictly ascending"))
    size(B) == (N, nres) || throw(DimensionMismatch("B must be (N, nres)"))
    size(R, 1) == N && size(R, 2) == nres || throw(DimensionMismatch("R must be (N, nres, neqres)"))
    size(A) == (N, M, nres) || throw(DimensionMismatch("A must be (N, M, nres)"))
    for f in (states.l, states.gw, states.nuion, states.mu_red)
        length(f) == nres || throw(DimensionMismatch("state vectors must have equal length"))
    end
    logq = [log_qnl_qe(consts, states.gw[m], states.nuion[m], states.mu_red[m], exp(lgTg[i])) for i in 1:N, m in 1:nres]
    return AtomicRateTable(Vector{Float64}(lgTg), Vector{Float64}(lgrho), Matrix{Float64}(B), Array{Float64,3}(R),
        Array{Float64,3}(A), states, consts, logq)
end

n_resolved(t::AtomicRateTable) = length(t.states.n)
n_eq(t::AtomicRateTable) = size(t.R, 3)

# routines.cpp:708-736 (locate_JC), ascending branch, 1-based: largest j with xx[j] <= x, clamped to 1 and n.
function locate(xx::Vector{Float64}, xin)
    x = _primal(xin)
    n = length(xx)
    x <= xx[1] && return 1
    x >= xx[n] && return n
    jl = 0
    ju = n
    while ju - jl > 1
        jm = (ju + jl) >> 1
        if x >= xx[jm + 1]
            jl = jm
        else
            ju = jm
        end
    end
    return jl + 1
end

# cpp:233-241
@inline function lagrange4(x, x1, x2, x3, x4)
    a1 = (x - x2) * (x - x3) * (x - x4) / ((x1 - x2) * (x1 - x3) * (x1 - x4))
    a2 = (x - x1) * (x - x3) * (x - x4) / ((x2 - x1) * (x2 - x3) * (x2 - x4))
    a3 = (x - x1) * (x - x2) * (x - x4) / ((x3 - x1) * (x3 - x2) * (x3 - x4))
    a4 = (x - x1) * (x - x2) * (x - x3) / ((x4 - x1) * (x4 - x2) * (x4 - x3))
    return (a1, a2, a3, a4)
end

@inline lagrange_weights(xx::Vector{Float64}, l::Int, x) = lagrange4(x, xx[l], xx[l + 1], xx[l + 2], xx[l + 3])

"""
    stencil_start(table, logTg, logrho) -> (lx, ly)

1-based first stencil index on each axis (native lines 263-285; native `lx` is this minus one). Throws
`RateTableDomainError` where native exits (out of table) or reads out of bounds (stencil would need row/column > N).
"""
function stencil_start(t::AtomicRateTable, logTg, logrho)
    N, M = length(t.lgTg), length(t.lgrho)
    lT, lr = _primal(logTg), _primal(logrho)
    if lT < t.lgTg[1] || lT > t.lgTg[N] || lr < t.lgrho[1] || lr > t.lgrho[M]
        throw(RateTableDomainError(:out_of_table, exp(Float64(_val(logTg))), exp(Float64(_val(logTg))) * exp(Float64(_val(logrho)))))
    end
    jx = locate(t.lgTg, logTg)
    jy = locate(t.lgrho, logrho)
    lx = jx >= 2 ? jx - 1 : jx
    ly = jy >= 2 ? jy - 1 : jy
    if lx + 3 > N || ly + 3 > M
        throw(RateTableDomainError(:stencil_exceeds_table, exp(Float64(_val(logTg))), exp(Float64(_val(logTg))) * exp(Float64(_val(logrho)))))
    end
    return lx, ly
end

# Primal value of a possibly-dual number without depending on ForwardDiff. All branch decisions use it because
# ForwardDiff.Dual comparisons break value ties by the partials (Dual(2000.0, 1) > 2000.0 is true), which would make the
# AD branch differ from the primal branch exactly at a switch point.
_primal(x::AbstractFloat) = x
_primal(x::Real) = _primal(getfield(x, :value))
_val(x) = _primal(x)

@inline function _A(t::AtomicRateTable, m::Int, lx::Int, ly::Int, a, b, Tg, db::Bool)
    if db                                                               # cpp:300
        return exp(log_qnl_qe(t, m, Tg))
    end
    fxy = zero(a[1] * b[1])
    for i in 1:4
        fx = zero(fxy)
        for k in 1:4
            fx += a[k] * (t.A[k + lx - 1, i + ly - 1, m] - t.B[k + lx - 1, m] - t.logq_knot[k + lx - 1, m])
        end
        fxy += b[i] * fx
    end
    if abs(_primal(exp(fxy)) - 1.0) <= eps_A_effective                  # cpp:320
        fxy = zero(fxy)
    end
    return exp(fxy + log_qnl_qe(t, m, Tg))
end

@inline function _B(t::AtomicRateTable, m::Int, lx::Int, a)
    fx = zero(a[1])
    for k in 1:4
        fx += a[k] * t.B[k + lx - 1, m]
    end
    return exp(fx)
end

@inline function _R(t::AtomicRateTable, m::Int, i::Int, lx::Int, a)
    i <= m && return zero(a[1])                                         # cpp:335 (native i<=m, 0-based == 1-based i<=m)
    fx = zero(a[1])
    for k in 1:4
        fx += a[k] * t.R[k + lx - 1, m, i]
    end
    return exp(fx)
end

@inline function _setup(t::AtomicRateTable, Tg, Te)
    (isfinite(Tg) && isfinite(Te) && Tg > 0 && Te > 0) || throw(RateTableDomainError(:not_finite, Float64(_val(Tg)), Float64(_val(Te))))
    logTg = log(Tg)
    logrho = log(Te / Tg)
    lx, ly = stencil_start(t, logTg, logrho)
    a = lagrange_weights(t.lgTg, lx, logTg)
    b = lagrange_weights(t.lgrho, ly, logrho)
    db = _primal(Tg) / T0_detailed_balance - 1.0 > z_detailed_balance
    return lx, ly, a, b, db
end

"""Result of [`get_rates`](@ref): `A[m]`, `B[m]`, `R[m, i]` (row = resolved state m, column i; `R[m, i] = 0` for `i <= m`, as native)."""
struct RateResult{T<:Real}
    A::Vector{T}
    B::Vector{T}
    R::Matrix{T}
end

"""
    get_rates(table, Tg, Te) -> RateResult

Non-mutating port of native `get_rates(Tg, Te, Ai, Bi, RijVec)` (Tg, Te in K). Differentiable with ForwardDiff / Mooncake
(piecewise: see file header).
"""
function get_rates(t::AtomicRateTable, Tg::Real, Te::Real)
    lx, ly, a, b, db = _setup(t, Tg, Te)
    n, neq = n_resolved(t), n_eq(t)
    A = [_A(t, m, lx, ly, a, b, Tg, db) for m in 1:n]
    B = [_B(t, m, lx, a) for m in 1:n]
    R = [_R(t, m, i, lx, a) for m in 1:n, i in 1:neq]
    return RateResult(A, B, R)
end

"""In-place variant mirroring the native signature; returns `nothing`."""
function get_rates!(A::AbstractVector, B::AbstractVector, R::AbstractMatrix, t::AtomicRateTable, Tg::Real, Te::Real)
    lx, ly, a, b, db = _setup(t, Tg, Te)
    n, neq = n_resolved(t), n_eq(t)
    (length(A) == n && length(B) == n && size(R) == (n, neq)) || throw(DimensionMismatch("A, B must have length n_resolved and R size (n_resolved, neqres)"))
    for m in 1:n
        A[m] = _A(t, m, lx, ly, a, b, Tg, db)
        B[m] = _B(t, m, lx, a)
        for i in 1:neq
            R[m, i] = _R(t, m, i, lx, a)
        end
    end
    return nothing
end

# Native file format (cpp:112-184): line1 "NTg Nrho neqres"; then NTg values of lgTg, Nrho values of lgrho;
# then per Tg row: Bi, Bitot, neqres Rij, Nrho Ai (whitespace separated).
function read_native_rate_file(io::IO)
    tok = split(read(io, String))
    pos = 1
    nextf() = (v = parse(Float64, tok[pos]); pos += 1; v)
    N = Int(nextf()); M = Int(nextf()); neq = Int(nextf())
    lgTg = [nextf() for _ in 1:N]
    lgrho = [nextf() for _ in 1:M]
    B = zeros(N); Bitot = zeros(N); R = zeros(N, neq); A = zeros(N, M)
    for i in 1:N
        B[i] = nextf(); Bitot[i] = nextf()
        for k in 1:neq; R[i, k] = nextf(); end
        for j in 1:M; A[i, j] = nextf(); end
    end
    pos == length(tok) + 1 || throw(ArgumentError("unexpected trailing data in rate file"))
    return (; lgTg, lgrho, B, Bitot, R, A)
end

"""
    load_rate_table(dir, nS_effective, states) -> AtomicRateTable

Explicit loader for a native `Rates_n{n}_l{l}.nS_{nS_effective}.dat` set. `states` fixes the resolved levels (native: res_state_list.dat);
spin weights / ionisation frequencies / reduced masses are supplied by the caller (native derives them from its atom model).
"""
function load_rate_table(dir::AbstractString, nS_effective::Integer, states::ResolvedStates, consts::RateConstants = NATIVE_CONSTANTS)
    parts = [open(read_native_rate_file, joinpath(dir, "Rates_n$(states.n[m])_l$(states.l[m]).nS_$(nS_effective).dat")) for m in eachindex(states.n)]
    N, M, neq = length(parts[1].lgTg), length(parts[1].lgrho), size(parts[1].R, 2)
    for p in parts
        (p.lgTg == parts[1].lgTg && p.lgrho == parts[1].lgrho && size(p.R, 2) == neq) || throw(ArgumentError("rate files disagree on axes"))
    end
    nres = length(parts)
    B = [parts[m].B[i] for i in 1:N, m in 1:nres]
    R = [parts[m].R[i, k] for i in 1:N, m in 1:nres, k in 1:neq]
    A = [parts[m].A[i, j] for i in 1:N, j in 1:M, m in 1:nres]
    return AtomicRateTable(parts[1].lgTg, parts[1].lgrho, B, R, A, states, consts)
end
