# PRIVATE (unexported) workspace RHS (docs/ODE_CORE_WORKSPACE_DESIGN.md, increment 1: primal only). The public allocating `recombination_rhs!`
# is unchanged and remains the oracle. A workspace is built for ONE solve from its actual argument types (state, z, hscale, nbscale, feedback);
# `_recombination_rhs_ws!` uses it only when the call's argument types are exactly those (otherwise it runs the public path and returns false),
# so no value is ever converted to a narrower/different type. Buffer element types are those the allocating path produces for these argument
# types: X and g as in `recombination_rhs!`; H-rate A in two buffers (detailed-balance branch: Tg-typed; otherwise promote(Tg, Te)-typed, exactly
# the element type of `get_rates`' comprehension in each branch), H-rate B/R Tg-typed; He rates of `get_helium_rates`' own type T (which converts
# A/B/R into T as before). Same math: the existing primitives are called with these buffers.
struct _RHSWorkspace{TM,TY,TZ,TH,TN,TF,TX,TG,TA,TAdb,TB,TR,THe}
    X::Vector{TX}
    g::Vector{TG}
    dXH::Vector{TG}
    dXHe::Vector{TG}
    HA::Vector{TA}
    HAdb::Vector{TAdb}
    HB::Vector{TB}
    HR::Matrix{TR}
    HeA::Vector{THe}
    HeB::Vector{THe}
    HeR::Matrix{THe}
end

"""Private: workspace for the solve whose state, start redshift and parameters are `y0`, `z0`, `hscale`, `nbscale` (types taken from them)."""
function _rhs_workspace(y0::AbstractVector, z0, rm::RecombinationModel; hscale = 1.0, nbscale = 1.0)
    bg = ode_background(rm, z0; hscale = hscale, nbscale = nbscale)
    TX = promote_type(eltype(y0), typeof(bg.fHe))
    TF = feedback_eltype(rm.diffusion)
    TG = promote_type(TX, typeof(bg.Tg), typeof(bg.NH), typeof(bg.Hz), typeof(z0), TF)
    m = rm.eff; ht = m.htable; Tg = bg.Tg
    # element types only, from the same return expressions at a ZERO interpolant (no table value is interpolated or exponentiated at a fictitious
    # stencil): _A non-detailed-balance = exp(fxy + log_qnl_qe(t, m, Tg)) with fxy::typeof(zero(a[1] * b[1])); _A detailed balance =
    # exp(log_qnl_qe(t, m, Tg)) (the physical value at this Tg); _B/_R = exp(fx) with fx::typeof(zero(a[1])); He T = typeof(a[1] * logTg + Tg).
    # The Lagrange weights are pure polynomial arithmetic (finite for finite input); stencil start 1 avoids any table-domain lookup, so e.g. a
    # 7-state solve below the He table never touches He rates, exactly as the public path.
    rho = convert(TX, y0[1]); logTg = log(Tg)
    a = lagrange_weights(ht.lgTg, 1, logTg); b = lagrange_weights(ht.lgrho, 1, log(rho * Tg / Tg))
    lq = log_qnl_qe(ht, 1, Tg)
    TA = typeof(exp(zero(a[1] * b[1]) + lq)); TAdb = typeof(exp(lq))
    TB = typeof(exp(zero(a[1]))); TR = TB
    he = m.hetable; ah = lagrange_weights(he.lgTg, 1, logTg)
    THe = typeof(ah[1] * logTg + Tg)
    n, neq, nh = n_resolved(ht), n_eq(ht), nres(he)
    return _RHSWorkspace{typeof(rm),eltype(y0),typeof(z0),typeof(hscale),typeof(nbscale),TF,TX,TG,TA,TAdb,TB,TR,THe}(
        Vector{TX}(undef, EFF_NEQ), Vector{TG}(undef, EFF_NEQ), Vector{TG}(undef, EFF_NHI), Vector{TG}(undef, EFF_NHE),
        Vector{TA}(undef, n), Vector{TAdb}(undef, n), Vector{TB}(undef, n), Matrix{TR}(undef, n, neq),
        zeros(THe, nh), zeros(THe, nh), zeros(THe, nh, nh))           # He R: entries i <= m stay zero, as in get_helium_rates
end

# Contract: a workspace holds element types and buffer SIZES only, no model values. It is valid for any model of the same type whose tables have the
# buffer dimensions (checked here), with the same argument types it was built for; anything else runs the public path (the caller counts it).
_ws_compatible(ws::_RHSWorkspace{TM,TY,TZ,TH,TN,TF}, y, z, hscale, nbscale, rm) where {TM,TY,TZ,TH,TN,TF} =
    rm isa TM && eltype(y) === TY && z isa TZ && hscale isa TH && nbscale isa TN && feedback_eltype(rm.diffusion) === TF &&
    length(ws.HB) == n_resolved(rm.eff.htable) && size(ws.HR, 2) == n_eq(rm.eff.htable) && length(ws.HeA) == nres(rm.eff.hetable)

"""
Private: `recombination_rhs!` through the workspace `ws` (same result as the public path). Returns `true` if the workspace was used, `false` if the
argument types differ from those `ws` was built for, in which case the public allocating path filled `f`.
"""
function _recombination_rhs_ws!(f::AbstractVector, z, y::AbstractVector, rm::RecombinationModel, ws::_RHSWorkspace; flag_He::Bool = true, hscale = 1.0, nbscale = 1.0)
    if !_ws_compatible(ws, y, z, hscale, nbscale, rm)
        recombination_rhs!(f, z, y, rm; flag_He = flag_He, hscale = hscale, nbscale = nbscale)
        return false
    end
    length(y) == ode_nstate(flag_He) || throw(DimensionMismatch("state length $(length(y)) != $(ode_nstate(flag_He)) for flag_He = $flag_He"))
    bg = ode_background(rm, z; hscale = hscale, nbscale = nbscale)
    X = ws.X
    fill!(X, eltype(X)(ODE_UNPACK_FLOOR))                              # ode_unpack's default floor
    _ode_unpack!(X, y, bg.fHe, flag_He)
    _fcn_effective_core!(ws.g, z, X, bg, rm.eff, ws.dXH, ws.dXHe, ws, rm.dp_fallback, flag_He, rm.diffusion)
    _pack_rhs!(f, ws.g, bg, z, flag_He)
    return true
end

# rates from the workspace (the allocating counterparts: `hydrogen_rhs!(dX, table, ...)` -> `get_rates`, `helium_base_rhs!(dX, table, ...)` -> `get_helium_rates`)
function _hydrogen_rhs_rates!(dX, ws::_RHSWorkspace, table::AtomicRateTable, Tg, rho, Xe, Xp, fHe, NH, Hz, X, atom, c)
    lx, ly, a, b, db = _setup(table, Tg, rho * Tg)
    if db
        _fill_hrates!(ws.HAdb, ws.HB, ws.HR, table, lx, ly, a, b, Tg, true)
        hydrogen_rhs!(dX, Tg, Xe, Xp, NH, Hz, X, ws.HAdb, ws.HB, ws.HR, atom, c)
    else
        _fill_hrates!(ws.HA, ws.HB, ws.HR, table, lx, ly, a, b, Tg, false)
        hydrogen_rhs!(dX, Tg, Xe, Xp, NH, Hz, X, ws.HA, ws.HB, ws.HR, atom, c)
    end
    return matter_temperature_rate(rho, Tg, Xe, fHe, Hz, c)
end
function _helium_base_rhs_rates!(dX, ws::_RHSWorkspace, table::HeliumRateTable, Tg, Xe, NH, Hz, X, fHe, atom, c, rc, spin_forbidden::Bool)
    logTg = log(Tg); j = helium_stencil_start(table, logTg)            # as get_helium_rates
    _fill_herates!(ws.HeA, ws.HeB, ws.HeR, table, Tg, c, j, lagrange_weights(table.lgTg, j, logTg))
    return helium_base_rhs!(dX, Tg, Xe, NH, Hz, X, fHe, ws.HeA, ws.HeB, ws.HeR, atom, c, rc; spin_forbidden = spin_forbidden)
end
